# KRAS G12D / HLA-A\*11:01 新抗原 mRNA 疫苗可行性研究 — 双栈结构预测与表位筛选

> 一个完全在本地完成的、面向胰腺导管腺癌(PDAC)的 KRAS G12D 新抗原疫苗可行性管线:
> ColabFold(AF2)+ ESMFold v1 双栈结构验证 + HLA-A\*11:01 P2/P9 锚点表位筛选。

## ⚠️ 重要声明 / Important Notice

**本仓库中的所有代码、工作流脚本、日志、图表与文档,均由 Anthropic 的 Claude(全程使用 Opus 4.7 模型)在 Claude Code 命令行环境下,根据作者设计的 prompt（由claude辅助设计） 辅助生成、运行与整理。作者负责科学问题的设定、prompt 的设计、对中间结果的人工核查,以及最终结论的取舍。**

请在评价本工作时,把上述"AI 辅助"事实视为方法学的一部分,而非隐含的人工劳动。

---

## 1. 背景与动机

胰腺导管腺癌(PDAC)是五年生存率最低的实体瘤之一,约 **90% 的肿瘤携带 KRAS 突变**,其中 G12D 是最常见的亚型。2023 年 Rojas 等人在 *Nature* 报告的 autogene cevumeran(BNT122)Phase 1 试验首次证明,**个体化 mRNA 新抗原疫苗在切除后 PDAC 患者中可诱导持久的新抗原特异性 T 细胞应答**[1]。这一结果为"通用型" KRAS G12D 疫苗(不再针对个体突变定制,而是面对一个高频共享突变)提供了机会窗口。

本研究的目标是,在一台普通的笔记本(Windows 11 + WSL2 + 12 GB GPU)上,从零搭建一个最小可行的疫苗候选筛选管线,回答两个问题:

1. **结构层面**:在双栈(AF2 + ESMFold)预测下,KRAS G12D 与野生型在全长(残基 1–169)和 Switch II 区域(残基 60–76)是否仍保持相同的全局折叠?
2. **表位层面**:KRAS G12D 的 9 残基新抗原中,有哪一段可以同时满足 HLA-A\*11:01 的 P2 疏水/小极性 + P9 K 锚点要求,并且穿过 G12D 的突变位点?

## 2. 流程概览

```
KRAS 1–169 (WT + G12D 序列)
        │
        ├── ColabFold (AF2, MMseqs2 在线 MSA)  ──┐
        │                                       ├── Cα RMSD (Kabsch)
        └── ESMFold v1 (单序列,本地 GPU 推理) ──┘     ↓
                                              图表 + 结构 PASS/FAIL 判定
                                                       ↓
KRAS G12D 周围 9-mer 滑窗 ── HLA-A*11:01 P2/P9 锚点打分 ── 表位 PASS 候选
```

所有结构预测都在本地完成,只有 ColabFold 的 MSA 步骤通过 MMseqs2 在线服务排队。

## 3. 仓库结构

```
GithubPublish/
├── README.md                  ← 本文件
├── LICENSE                    ← MIT
├── RESEARCH_README.md         ← 原始研究说明(Step 6 自审计前的版本)
├── .gitignore
├── prompts/
│   └── RESEARCH_PROMPT.md     ← 研究使用的原始 prompt 占位文件
├── tools/                     ← 自写工具(mhc.py、rmsd.py、figures 等)
│   └── tests/                 ← pytest 单元测试
├── data/
│   ├── sequences/             ← KRAS WT / G12D FASTA(残基 1–169)
│   └── structures/            ← 仅包含每个栈的 rank_001 PDB
├── results/
│   ├── kras_wt_vs_g12d_structural.json
│   └── neoantigen_9mers.json
├── figures/                   ← 3 张 PNG + 对应说明文本
├── refs/
│   ├── citations.bib          ← 已逐条核验的 14 条参考文献
│   └── notes.md               ← 引用校对日志(包含 3 处更正)
└── logs/                      ← 关键 Markdown 日志白名单(见下)
```

`logs/` 内含:`00_environment.md`、`01_install.md`、`01_footprint.md`、`02_helper_toolkit.md`、`03_core_workflow.md`、`04_figures.md`、`06_self_audit.md`、`downloads.md`、`FINAL_RESEARCH_SUMMARY.md`、`PRIVACY_SCAN.md`。其余中间过程日志已剔除,以避免发布无关噪声。

## 4. 环境与安装

**最低硬件**(实际跑通本管线的笔记本配置):

| 项目 | 配置 |
|---|---|
| GPU | RTX 5070 Ti Laptop GPU(Blackwell sm\_120,12 GB GDDR7 VRAM) |
| 系统内存 | 16 GB |
| 系统 | Windows 11 + **WSL2** (Ubuntu 24.04) |
| Python | 3.10 – 3.12(本研究用 3.12) |
| NVIDIA 驱动 | ≥ 570(CUDA 12.x)或 ≥ 580(CUDA 13) |

**重要**:本管线只在 Linux/WSL 环境下验证通过。Windows 原生 Python 下 ColabFold 的依赖(JAX/CUDA)安装路径与本仓库的脚本不一致,**强烈建议在 WSL2 中运行**。

完整的依赖列表与逐条 `pip` 命令见 `logs/01_install.md`,核心依赖为:

- `colabfold` + `jax[cuda12_pip]`(AF2 后端)
- `fair-esm` 或 `esm`(ESMFold v1)
- `biopython`(SVDSuperimposer)、`matplotlib`、`numpy`、`pandas`

## 5. 复现步骤

```bash
# 1. 在 WSL2 中
git clone <this-repo>
cd <this-repo>
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt    # 见 logs/01_install.md

# 2. 结构预测(各自约 30 – 60 分钟,GPU)
colabfold_batch data/sequences/kras_wt_1-169.fasta   colabfold/kras_wt_1-169/
colabfold_batch data/sequences/kras_g12d_1-169.fasta colabfold/kras_g12d_1-169/
python tools/run_esmfold.py --fasta data/sequences/kras_wt_1-169.fasta   --out esmfold/kras_wt_1-169/
python tools/run_esmfold.py --fasta data/sequences/kras_g12d_1-169.fasta --out esmfold/kras_g12d_1-169/

# 3. 结构比对 + 表位筛选
python tools/rmsd.py        # 写 results/kras_wt_vs_g12d_structural.json
python tools/mhc.py         # 写 results/neoantigen_9mers.json
python tools/make_figures.py
```

完整命令、参数、随机种子、recycle 次数等已逐条记录在 `logs/03_core_workflow.md`。

## 6. 关键结果

**6.1 结构层面 — KRAS G12D 不改变全局折叠**

| 度量 | AF2 (ColabFold) | ESMFold v1 |
|---|---|---|
| rank-1 平均 pLDDT(WT) | **93.66** | 86.06 |
| rank-1 平均 pLDDT(G12D) | **93.28** | 86.08 |
| 全长 Cα RMSD,WT vs G12D(同栈) | **0.102 Å** | 0.081 Å |
| Switch II (60–76) Cα RMSD,WT vs G12D(同栈) | **0.175 Å** | 0.134 Å |

同栈对比下,WT 与 G12D 的全局 Cα RMSD < 0.2 Å,与 Krengel 1990 *Cell* 报告的 H-Ras G12 突变体晶体结果一致[9]。

**6.2 双栈交叉验证 — Switch II 仍稳定**

跨栈(AF2 vs ESMFold)的同序列对比:

- WT, Switch II Cα RMSD ≈ **1.81 Å**
- G12D, Switch II Cα RMSD ≈ **1.87 Å**

跨栈 Switch II 差异约为 1.8 Å,反映两个模型对柔性环区的预测分歧,但仍处于 H-Ras Switch II 实验晶体结构间已知的构象波动范围内[8,9]。**没有任何一种预测在 G12D 引入了大于 0.2 Å 的额外构象偏离**——即从结构层面,G12D 不是"全局折叠破坏型"突变,而是"局部催化破坏型"突变[8]。

**6.3 表位层面 — VVGADGVGK 是唯一 PASS**

对 KRAS G12D 围绕第 12 位的所有 9-mer 滑窗,应用 HLA-A\*11:01 P2/P9 锚点规则(P2 = T/V/M/L/I/S/A,P9 = K)[5,6,7]:

| 9-mer | 残基范围 | P2 | P9 | 覆盖 G12D? | 判定 |
|---|---|---|---|---|---|
| **VVGADGVGK** | 8–16 | V ✓ | K ✓ | ✓ | **PASS** |
| 其他 9-mer | — | — | — | — | FAIL |

`VVGADGVGK` 对应 Wang 2016 *Cancer Immunol Res* 在 HLA-A\*11:01 转基因小鼠中已实验验证的、可被 G12D-特异性 TCR 识别的 9 残基表位[2]。请注意,Tran 2016 *NEJM*[3] 与 Leidner 2022 *NEJM*[4] 的 KRAS G12D-定向 TCR 治疗病例使用的是 **HLA-C\*08:02 限制性表位 GADGVGKSA / GADGVGKSAL**,而非 HLA-A\*11:01——这两个工作只作为"G12D-TCR 临床概念验证"被引用,而不是 A\*11:01/9-mer 的来源。

## 7. 注意事项与局限

1. **9-mer vs 10-mer 之争。** Zhu 等人 2026 年在 *Communications Biology* 上的 NMR + X-ray + MD 联合研究指出,VVGADGVGK 9-mer 虽然在锚点规则上可通过,但在 A\*11:01 凹槽内的稳定构象**显著弱于 10-mer 延伸 VVVGADGVGK**(对应 PDB 7OW4)[14]。本仓库的 `tools/mhc.py` 仅做锚点合规筛选,未做结构层级评分;**若要把这个候选推向湿实验,应当同时合成 9-mer 与 10-mer 并以 10-mer 为主要假设**。
2. **MSA 在线依赖。** ColabFold 通过 MMseqs2 在线服务排队进行 MSA;`esmfold` 完全本地。这意味着 AF2 端的可复现性依赖于 MMseqs2 服务的版本与序列数据库快照。
3. **没有 AFDB baseline。** 本研究**没有**与 AlphaFold Database (AFDB) 中的 KRAS 公开模型做对比,所有 AF2 结构均为本机重跑。
4. **没有湿实验。** 全部结果都是 in-silico 推理,**未经任何细胞或动物实验验证**。
5. **HLA 仅覆盖 A\*11:01。** 没有评估 G12D 在其他 HLA 等位基因(尤其 HLA-C\*08:02 / A\*03:01)上的呈递。
6. **`Copyright (c) 2026 [Your Name]`。** `LICENSE` 文件保留了 `[Your Name]` 占位符,首次发布前请替换为作者真实姓名或机构名。


## 8. 如何引用本仓库

```
@software{kras_g12d_neoantigen_feasibility_2026,
  title  = {KRAS G12D / HLA-A*11:01 neoantigen mRNA vaccine feasibility:
            dual-stack structural prediction and epitope screening},
  author = {[Your Name]},
  year   = {2026},
  url    = {https://github.com/<owner>/<repo>},
  note   = {AI-assisted by Claude Opus 4.7 in Claude Code.}
}
```

## 9. 许可证

本仓库以 **MIT License** 发布(见 `LICENSE`)。

---

# KRAS G12D / HLA-A\*11:01 Neoantigen mRNA Vaccine Feasibility — Dual-Stack Structural Prediction and Epitope Screening

> An end-to-end, laptop-scale feasibility pipeline for a KRAS-G12D-directed mRNA neoantigen vaccine in pancreatic ductal adenocarcinoma (PDAC):
> ColabFold (AF2) + ESMFold v1 cross-validation + HLA-A\*11:01 P2/P9 anchor-based 9-mer screening.

## ⚠️ Important Notice

**Every line of code, workflow script, log, figure, and documentation file in this repository was generated, executed, and curated by Anthropic's Claude (using the Opus 4.7 model throughout) inside the Claude Code command-line environment, under prompts authored by the human author. The author is responsible for the scientific framing, the prompt design（assistant by Claude）, the manual verification of intermediate results, and the final scientific conclusions.**

When evaluating this work, please treat the "AI-assisted" fact above as part of the methodology, not as implicit human labour.

## 1. Background and motivation

PDAC has one of the lowest 5-year survival rates among solid tumours, and roughly **90% of PDAC tumours carry a KRAS mutation**, of which G12D is the most common subtype. The 2023 *Nature* report by Rojas et al. on autogene cevumeran (BNT122) was the first Phase 1 trial to show that **individualised mRNA neoantigen vaccines can elicit durable neoantigen-specific T-cell responses in resected PDAC patients**[1]. That result opens a window for a "shared-mutation" KRAS G12D vaccine — one that does not need to be re-personalised per tumour, but targets a single high-prevalence driver.

The goal of this study is to build a minimum-viable candidate-screening pipeline from scratch on a consumer laptop (Windows 11 + WSL2 + 12 GB GPU) that answers two questions:

1. **Structurally**: under a dual-stack prediction (AF2 + ESMFold), does KRAS G12D retain the same global fold as wild-type across the full 1–169 G-domain and the Switch II region (residues 60–76)?
2. **Epitopically**: among the KRAS-G12D-spanning 9-mers, is there a peptide that satisfies the HLA-A\*11:01 anchor preferences at both P2 (hydrophobic / small polar) and P9 (K)?

## 2. Pipeline overview

```
KRAS 1–169 (WT and G12D sequences)
        │
        ├── ColabFold (AF2, MMseqs2 online MSA) ──┐
        │                                        ├── Cα RMSD (Kabsch)
        └── ESMFold v1 (single-sequence, local GPU) ─┘  ↓
                                              figures + structural PASS/FAIL
                                                        ↓
KRAS G12D 9-mer sliding window ── HLA-A*11:01 P2/P9 anchor scoring ── epitope PASS
```

All structure inference runs locally; only the ColabFold MSA step queues through the public MMseqs2 service.

## 3. Repository layout

```
GithubPublish/
├── README.md                  ← this file
├── LICENSE                    ← MIT
├── RESEARCH_README.md         ← original research README (pre-publish snapshot)
├── .gitignore
├── prompts/
│   └── RESEARCH_PROMPT.md     ← placeholder for the original research prompt
├── tools/                     ← in-house Python (mhc.py, rmsd.py, figures, …)
│   └── tests/                 ← pytest unit tests
├── data/
│   ├── sequences/             ← KRAS WT / G12D FASTA (residues 1–169)
│   └── structures/            ← rank_001 PDBs only, per stack
├── results/
│   ├── kras_wt_vs_g12d_structural.json
│   └── neoantigen_9mers.json
├── figures/                   ← 3 PNGs + per-figure captions
├── refs/
│   ├── citations.bib          ← 14 hand-verified BibTeX entries
│   └── notes.md               ← citation-verification log (3 corrections)
└── logs/                      ← whitelisted Markdown logs (see below)
```

`logs/` contains: `00_environment.md`, `01_install.md`, `01_footprint.md`, `02_helper_toolkit.md`, `03_core_workflow.md`, `04_figures.md`, `06_self_audit.md`, `downloads.md`, `FINAL_RESEARCH_SUMMARY.md`, `PRIVACY_SCAN.md`. Other intermediate logs were dropped from the bundle to keep the public release focused.

## 4. Environment and install

**Minimum tested hardware** (the laptop this pipeline actually ran on):

| Item | Value |
|---|---|
| GPU | RTX 5070 Ti Laptop GPU (Blackwell sm\_120, 12 GB GDDR7 VRAM) |
| System RAM | 16 GB |
| OS | Windows 11 + **WSL2** (Ubuntu 24.04) |
| Python | 3.10 – 3.12 (3.12 used here) |
| NVIDIA driver | ≥ 570 (CUDA 12.x) or ≥ 580 (CUDA 13) |

**Important**: this pipeline has only been validated under Linux / WSL. Native Windows Python paths for ColabFold's JAX/CUDA dependencies diverge from what the scripts in this repository assume — **strongly recommended to run inside WSL2**.

Full per-package install commands live in `logs/01_install.md`. Core dependencies:

- `colabfold` + `jax[cuda12_pip]` (AF2 backend)
- `fair-esm` or `esm` (ESMFold v1)
- `biopython` (SVDSuperimposer), `matplotlib`, `numpy`, `pandas`

## 5. Reproducing the run

```bash
# 1. Inside WSL2
git clone <this-repo>
cd <this-repo>
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt   # see logs/01_install.md

# 2. Structure prediction (each ~30–60 min on the laptop GPU)
colabfold_batch data/sequences/kras_wt_1-169.fasta   colabfold/kras_wt_1-169/
colabfold_batch data/sequences/kras_g12d_1-169.fasta colabfold/kras_g12d_1-169/
python tools/run_esmfold.py --fasta data/sequences/kras_wt_1-169.fasta   --out esmfold/kras_wt_1-169/
python tools/run_esmfold.py --fasta data/sequences/kras_g12d_1-169.fasta --out esmfold/kras_g12d_1-169/

# 3. Structural comparison + epitope screening
python tools/rmsd.py        # writes results/kras_wt_vs_g12d_structural.json
python tools/mhc.py         # writes results/neoantigen_9mers.json
python tools/make_figures.py
```

All commands, parameters, random seeds, and recycle counts are recorded line-by-line in `logs/03_core_workflow.md`.

## 6. Key results

**6.1 Structural — KRAS G12D does not perturb the global fold.**

| Metric | AF2 (ColabFold) | ESMFold v1 |
|---|---|---|
| rank-1 mean pLDDT, WT | **93.66** | 86.06 |
| rank-1 mean pLDDT, G12D | **93.28** | 86.08 |
| Full Cα RMSD, WT vs G12D (within stack) | **0.102 Å** | 0.081 Å |
| Switch II (60–76) Cα RMSD, WT vs G12D (within stack) | **0.175 Å** | 0.134 Å |

Within each stack, the global Cα RMSD between WT and G12D stays below 0.2 Å — consistent with the crystallographic record of H-Ras G12 mutants in Krengel et al. 1990 *Cell*[9].

**6.2 Cross-stack validation — Switch II is stable.**

For the same sequence across the two stacks:

- WT, Switch II Cα RMSD ≈ **1.81 Å**
- G12D, Switch II Cα RMSD ≈ **1.87 Å**

The ~1.8 Å cross-stack spread reflects model-level disagreement on a flexible loop, but stays within the conformational range known for H-Ras Switch II among experimental structures[8,9]. **In neither stack does G12D introduce an additional > 0.2 Å deviation** — i.e., G12D is a *catalysis-impairing* mutation, not a *global-fold-disrupting* one[8].

**6.3 Epitope — VVGADGVGK is the single PASS.**

Sliding-window 9-mers around residue 12 of KRAS G12D were filtered by HLA-A\*11:01 P2/P9 anchor rules (P2 ∈ {T,V,M,L,I,S,A}, P9 = K)[5,6,7]:

| 9-mer | Residues | P2 | P9 | Spans G12D? | Verdict |
|---|---|---|---|---|---|
| **VVGADGVGK** | 8–16 | V ✓ | K ✓ | ✓ | **PASS** |
| Others | — | — | — | — | FAIL |

`VVGADGVGK` matches the experimentally validated HLA-A\*11:01-restricted G12D 9-mer in Wang et al. 2016 *Cancer Immunol Res*[2]. Note that the KRAS-G12D-directed adoptive cell therapy reports from Tran 2016 *NEJM*[3] and Leidner 2022 *NEJM*[4] target **HLA-C\*08:02-restricted GADGVGKSA / GADGVGKSAL**, not HLA-A\*11:01 — they are cited here only as clinical proof-of-concept for KRAS-G12D-TCR therapy, not as the A\*11:01 9-mer source.

## 7. Caveats and limitations

1. **9-mer vs 10-mer.** Zhu et al. 2026 *Communications Biology* (NMR + X-ray + MD + TCR binding) reports that while VVGADGVGK passes anchor rules, **the 10-mer extension VVVGADGVGK adopts a substantially more stable, immunogenically productive conformation in the A\*11:01 groove** (cf. PDB 7OW4)[14]. `tools/mhc.py` does anchor-rule screening only — it does not score structural compatibility. **If you intend to take this candidate into wet-lab work, synthesise both the 9-mer and the 10-mer, and treat the 10-mer as the primary hypothesis.**
2. **MSA online dependency.** ColabFold queues MMseqs2 MSAs against the public service; ESMFold is fully local. AF2-side reproducibility therefore depends on the MMseqs2 service version and sequence-database snapshot at run time.
3. **No AFDB baseline.** This study did **not** compare against any public KRAS model from the AlphaFold Database; all AF2 structures were re-folded locally.
4. **No wet-lab validation.** All findings are *in silico*. Nothing here has been tested in cells or animals.
5. **HLA scope is A\*11:01 only.** No assessment of G12D presentation on other HLA alleles (notably HLA-C\*08:02 / A\*03:01).
6. **`Copyright (c) 2026 [Your Name]`.** The `LICENSE` file still carries the `[Your Name]` placeholder — replace it with the actual author / institution name before the first public release.


## 8. How to cite

```
@software{kras_g12d_neoantigen_feasibility_2026,
  title  = {KRAS G12D / HLA-A*11:01 neoantigen mRNA vaccine feasibility:
            dual-stack structural prediction and epitope screening},
  author = {[Your Name]},
  year   = {2026},
  url    = {https://github.com/<owner>/<repo>},
  note   = {AI-assisted by Claude Opus 4.7 in Claude Code.}
}
```

## 9. License

Released under the **MIT License** (see `LICENSE`).

---

## 10. References / 参考文献

References are rendered in Nature numbered style. DOIs are clickable links. All entries below have been verified against the publisher / PubMed / PMC record (see `refs/notes.md`).

1. Rojas, L. A. *et al.* Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer. *Nature* **618**, 144–150 (2023). [doi:10.1038/s41586-023-06063-y](https://doi.org/10.1038/s41586-023-06063-y)
2. Wang, Q. J. *et al.* Identification of T-cell receptors targeting KRAS-mutated human tumors. *Cancer Immunol. Res.* **4**, 204–214 (2016). [doi:10.1158/2326-6066.CIR-15-0188](https://doi.org/10.1158/2326-6066.CIR-15-0188)
3. Tran, E. *et al.* T-cell transfer therapy targeting mutant KRAS in cancer. *N. Engl. J. Med.* **375**, 2255–2262 (2016). [doi:10.1056/NEJMoa1609279](https://doi.org/10.1056/NEJMoa1609279)
4. Leidner, R. *et al.* Neoantigen T-cell receptor gene therapy in pancreatic cancer. *N. Engl. J. Med.* **386**, 2112–2119 (2022). [doi:10.1056/NEJMoa2119662](https://doi.org/10.1056/NEJMoa2119662)
5. Sidney, J., Peters, B., Frahm, N., Brander, C. & Sette, A. HLA class I supertypes: a revised and updated classification. *BMC Immunol.* **9**, 1 (2008). [doi:10.1186/1471-2172-9-1](https://doi.org/10.1186/1471-2172-9-1)
6. Falk, K. *et al.* Peptide motifs of HLA-A1, -A11, -A31, and -A33 molecules. *Immunogenetics* **40**, 238–241 (1994). [doi:10.1007/BF00167086](https://doi.org/10.1007/BF00167086)
7. Zhang, Q. J., Gavioli, R., Klein, G. & Masucci, M. G. An HLA-A11-specific motif in nonamer peptides derived from viral and cellular proteins. *Proc. Natl. Acad. Sci. USA* **90**, 2217–2221 (1993). [doi:10.1073/pnas.90.6.2217](https://doi.org/10.1073/pnas.90.6.2217)
8. Pai, E. F., Krengel, U., Petsko, G. A., Goody, R. S., Kabsch, W. & Wittinghofer, A. Refined crystal structure of the triphosphate conformation of H-ras p21 at 1.35 Å resolution: implications for the mechanism of GTP hydrolysis. *EMBO J.* **9**, 2351–2359 (1990). [doi:10.1002/j.1460-2075.1990.tb07409.x](https://doi.org/10.1002/j.1460-2075.1990.tb07409.x)
9. Krengel, U. *et al.* Three-dimensional structures of H-ras p21 mutants: molecular basis for their inability to function as signal switch molecules. *Cell* **62**, 539–548 (1990). [doi:10.1016/0092-8674(90)90018-A](https://doi.org/10.1016/0092-8674(90)90018-A)
10. Jumper, J. *et al.* Highly accurate protein structure prediction with AlphaFold. *Nature* **596**, 583–589 (2021). [doi:10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
11. Mirdita, M. *et al.* ColabFold: making protein folding accessible to all. *Nat. Methods* **19**, 679–682 (2022). [doi:10.1038/s41592-022-01488-1](https://doi.org/10.1038/s41592-022-01488-1)
12. Lin, Z. *et al.* Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* **379**, 1123–1130 (2023). [doi:10.1126/science.ade2574](https://doi.org/10.1126/science.ade2574)
13. Kabsch, W. A solution for the best rotation to relate two sets of vectors. *Acta Crystallogr. A* **32**, 922–923 (1976). [doi:10.1107/S0567739476001873](https://doi.org/10.1107/S0567739476001873)
14. Zhu, J. *et al.* Structure-guided analysis of KRAS G12 mutants in HLA-A\*11:01 reveals a length-encoded immunogenic advantage in G12D. *Commun. Biol.* **9**, 26 (2026). [doi:10.1038/s42003-025-09285-0](https://doi.org/10.1038/s42003-025-09285-0)
