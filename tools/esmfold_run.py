"""ESMFold inference helper.

Single-sequence (no MSA) structure prediction using facebook/esmfold_v1 from
HuggingFace transformers. Outputs a PDB (with pLDDT in B-factor column) plus
a metrics JSON.

Usage:
    python -m tools.esmfold_run --fasta path.fa --out-dir path/to/out --tag mytag

The script:
- loads ESMFold v1 (~3 GB safetensors download on first run)
- uses fp16 on CUDA to fit 12 GB VRAM
- predicts each FASTA record into <out-dir>/<tag>_<record_id>.pdb
- writes <out-dir>/<tag>_<record_id>_metrics.json with mean and per-residue pLDDT
- writes <out-dir>/<tag>_summary.json with all records' mean pLDDT and timing
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37


def read_fasta(path: Path):
    records = []
    rid, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    records.append((rid, "".join(seq)))
                rid = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if rid is not None:
            records.append((rid, "".join(seq)))
    return records


def convert_outputs_to_pdb(outputs):
    """Convert ESMFold raw outputs into a list of PDB strings (one per chain in the batch)."""
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
        aa = outputs["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs["residue_index"][i] + 1
        pred = OFProtein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs["plddt"][i],
            chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--model-id", default="facebook/esmfold_v1")
    parser.add_argument("--max-tokens-per-batch", type=int, default=1024,
                        help="Skip records longer than this to avoid OOM on 12 GB VRAM.")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if not torch.cuda.is_available():
        print("ERROR: torch.cuda.is_available() is False — refuse to run on CPU.", file=sys.stderr)
        return 2

    records = read_fasta(args.fasta)
    print(f"Loaded {len(records)} record(s) from {args.fasta}")

    t_load_start = time.time()
    tokenizer = AutoTokenizer.from_pretrained(args.model_id)
    # Two-stage precision setup to live within WSL2 7.4 GiB host RAM AND avoid
    # fp16 NaN in pTM/pAE on this GPU:
    #  1) Load everything in fp16 on CPU/disk → ~6.8 GB peak host RAM (fits).
    #  2) Move to GPU, then upcast the folding trunk + heads to fp32 in place.
    # ESM-2 stem (the bulk of weights) stays fp16; trunk/heads are small enough
    # (~few hundred MB) that fp32 still fits in 12 GB VRAM for short sequences.
    model = EsmForProteinFolding.from_pretrained(args.model_id, dtype=torch.float16)
    model = model.eval().cuda()
    # Upcast everything except the ESM-2 stem to fp32. Adapter modules
    # (esm_s_mlp, embedding, etc.) and direct Parameters on the model
    # (esm_s_combine, af2_to_esm) must also be promoted; otherwise
    # silent fp16 propagation chains back through the trunk.
    for child_name, child_mod in model.named_children():
        if child_name == "esm":
            continue
        child_mod.float()
    for p_name, p in model.named_parameters(recurse=False):
        if p.is_floating_point():
            p.data = p.data.float()
    for b_name, b in model.named_buffers(recurse=False):
        if b.is_floating_point():
            b.data = b.data.float()
    model.trunk.set_chunk_size(64)
    t_load = time.time() - t_load_start
    print(f"Loaded model in {t_load:.1f}s ({args.model_id})")

    summary = {
        "model_id": args.model_id,
        "fasta": str(args.fasta),
        "load_seconds": round(t_load, 2),
        "records": [],
    }

    torch.manual_seed(42)

    for rid, seq in records:
        L = len(seq)
        if L > args.max_tokens_per_batch:
            print(f"  skip {rid}: length {L} > max_tokens_per_batch {args.max_tokens_per_batch}")
            continue
        t0 = time.time()
        tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False).to("cuda")
        with torch.no_grad():
            outputs = model(**tokenized)
        elapsed = time.time() - t0
        # HF transformers ESMFold returns pLDDT on a 0-1 scale (categorical mixture
        # over [0, 1]); convert to the canonical 0-100 scale to match AF2/ColabFold
        # conventions and the ≥85 hard gate. We expose three views:
        #   - per-residue (Cα slot = atom37 index 1) → matches AF2's per-residue convention
        #   - per-residue (mean over masked heavy atoms)
        #   - aggregate means + a "structured core" mean that drops the flexible tail
        atom_mask = outputs["atom37_atom_exists"][0]  # [L, 37]
        plddt_per_atom = outputs["plddt"][0]  # [L, 37]
        masked_sum = (plddt_per_atom * atom_mask).sum(dim=-1)
        atom_counts = atom_mask.sum(dim=-1).clamp(min=1)
        plddt_atom_mean = (masked_sum / atom_counts).cpu().numpy() * 100.0  # [L]
        plddt_ca = plddt_per_atom[:, 1].cpu().numpy() * 100.0  # Cα = atom37[1]
        mean_plddt = float(plddt_atom_mean.mean())
        mean_plddt_ca = float(plddt_ca.mean())
        per_res_plddt = [float(v) for v in plddt_atom_mean.tolist()]
        per_res_plddt_ca = [float(v) for v in plddt_ca.tolist()]

        # Promote the plddt tensor used by to_pdb to the 0-100 scale so the PDB
        # B-factor column follows AF2/ColabFold convention.
        outputs["plddt"] = outputs["plddt"] * 100.0

        # Write PDB
        pdbs = convert_outputs_to_pdb(outputs)
        pdb_path = args.out_dir / f"{args.tag}_{rid}.pdb"
        pdb_path.write_text(pdbs[0])

        # Write per-record metrics
        metrics_path = args.out_dir / f"{args.tag}_{rid}_metrics.json"
        metrics_path.write_text(json.dumps({
            "record_id": rid,
            "length": L,
            "mean_plddt": mean_plddt,
            "mean_plddt_ca": mean_plddt_ca,
            "per_residue_plddt": per_res_plddt,
            "per_residue_plddt_ca": per_res_plddt_ca,
            "inference_seconds": round(elapsed, 2),
            "plddt_scale": "0-100 (HF native 0-1 scaled x100)",
            "plddt_convention": "per_residue_plddt = mean over masked heavy atoms; per_residue_plddt_ca = atom37[1]=Cα slot",
        }, indent=2))

        print(f"  {rid} L={L} mean_pLDDT={mean_plddt:.2f} in {elapsed:.1f}s -> {pdb_path.name}")
        summary["records"].append({
            "record_id": rid,
            "length": L,
            "mean_plddt": mean_plddt,
            "inference_seconds": round(elapsed, 2),
        })

    (args.out_dir / f"{args.tag}_summary.json").write_text(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
