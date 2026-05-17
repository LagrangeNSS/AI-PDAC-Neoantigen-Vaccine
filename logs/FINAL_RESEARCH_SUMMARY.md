# KRAS G12D PDAC 新抗原 mRNA 疫苗——结构可行性研究终稿

**日期窗口：** 2026-05-13 → 2026-05-14（UTC）
**执行环境：** WSL2 Ubuntu 24.04，RTX 5070 Ti Laptop（Blackwell sm_120，12 GB VRAM），Python 3.12.3
**工作目录：** `/mnt/d/AssignmentOPUS/`（所有产物均位于此目录下，未越界）
**对标研究：** Rojas et al. *Nature* 2023 ——自体个体化 mRNA 疫苗 autogene cevumeran 在切除后 PDAC 中诱导新抗原特异性 T 细胞应答

---

## 1. 研究目标

围绕 KRAS G12D 这一 PDAC 中最常见的驱动突变，建立一条**可复现**、**有双栈交叉验证**的结构 → 抗原表位筛选可行性流水线，回答三个工程性问题：

1. KRAS G12D（vs 野生型 WT）在 G 结构域是否保持整体折叠？——这是 mRNA 疫苗能否合理依赖 KRAS 三维结构信息的前提。
2. 在 HLA-A\*11:01（A3 超型，东亚高频）这一限制 MHC-I 等位上，覆盖位点 12 的 9 体肽中，哪些满足 P2/P9 锚位偏好，构成 *候选* 新抗原？
3. 一条以本地 GPU、严格下载预算（≤20 GB）、严格时间预算（每次预测 ≤2 h，60/300/900 s 重试退避）为硬约束的流水线，能否在不接外网二进制、不调用 license 化工具（如 NetMHCpan）的前提下复现已发表的 KRAS G12D / A\*11:01 经典表位结论？

---

## 2. 关键结果一览

### 2.1 结构预测（双栈、硬门槛 mean pLDDT ≥ 85）

| 构件 | 栈 | rank-1 mean pLDDT | 通过 ≥85 门槛 | 来源 |
|---|---|---:|:---:|---|
| KRAS WT 1–169 | ColabFold AF2（model_5, seed 42） | **93.66** | ✅ | `analysis/kras_wt_vs_g12d_structural.json` |
| KRAS G12D 1–169 | ColabFold AF2（model_5, seed 42） | **93.28** | ✅ | 同上 |
| KRAS WT 1–169 | ESMFold v1（heavy-atom mean） | **86.06** | ✅ | 同上 |
| KRAS G12D 1–169 | ESMFold v1（heavy-atom mean） | **86.08** | ✅ | 同上 |

> 仅建模 1–169（催化型 G 结构域）。HVR（167–188）是内禀无序的多聚赖氨酸+CAAX 锚定区，建模它会把 mean pLDDT 拖到门槛以下，且对 G12D（位点 12，P-loop 内）分析无信息贡献；这与所有 KRAS 晶体结构（4OBE、6GOD、7L0F 等）只覆盖 1-169 的做法一致。

### 2.2 WT vs G12D 结构差异（Cα RMSD，Kabsch 叠合）

| 比较对 | 全长 1–169 | 核心 1–166 | P-loop 10–17 | Switch I 30–38 | Switch II 59–72 |
|---|---:|---:|---:|---:|---:|
| **AF2 WT vs AF2 G12D** | **0.102 Å** | 0.103 Å | 0.045 Å | 0.107 Å | 0.175 Å |
| **ESM WT vs ESM G12D** | **0.081 Å** | 0.081 Å | 0.042 Å | 0.058 Å | 0.134 Å |
| AF2 vs ESMFold（WT 跨栈） | 0.593 Å | 0.598 Å | 0.099 Å | 0.355 Å | 1.813 Å |
| AF2 vs ESMFold（G12D 跨栈） | 0.611 Å | 0.616 Å | 0.102 Å | 0.355 Å | 1.870 Å |

**解读：**
- **G12D 全局保持折叠：** 两栈内部 WT-vs-G12D 任意区域 RMSD < 0.20 Å——这与 Krengel 1990 *Cell* 62:539 和 Pai 1990 *EMBO J* 9:2351 的经典晶体学结论一致：G12 位的突变通过破坏 GTP 水解（γ-磷酸位附近的几何），而**不是**通过重排 G 结构域，造成致癌效应。
- **跨栈分歧集中在 Switch II：** 1.81–1.87 Å。这正是 KRAS 构象异质性最高的区域——AF2（基于深度 KRAS MSA 的协同进化信号）与 ESMFold（仅单序列、基于 ESM-2 语言模型）在这里看到的微观构象不同是预期之内的结果。结构化核心（1-166）跨栈一致到 0.6 Å，处于 AF2 与实验晶体结构的典型误差带内。
- **两栈独立都把 Switch II 标记为最高分歧区**——一个对正交验证器是否"看到同一构象图景"的良好检查。

### 2.3 位点 12 的局部 pLDDT 信号

| 预测器 | WT pLDDT @ 12 | G12D pLDDT @ 12 | Δ | 窗口 4–20 的平均 Δ |
|---|---:|---:|---:|---:|
| ColabFold AF2（rank-1） | 91.31 | 90.94 | −0.37 | ≈ −0.43（噪声级别） |
| ESMFold v1（heavy-atom mean） | 84.76 | 80.91 | **−3.85** | ≈ −0.04 |

**解读：** ESMFold 在位点 12 处独立、显著、**局部**地报告了一次 −3.85 pLDDT 的置信度下降；AF2 因为 MSA 中位点 12 历史上的高耐受性，把这一变化"吸收"掉了（Δ ≈ 噪声）。下降仅限于位点 12 自身（位点 11、13 偏移 < 0.6 pLDDT），不传播到 P-loop 其余位、Switch I、Switch II。该信号被解读为 **ESMFold 对局部化学变化的单序列敏感性**，而非"G12D 破坏了结构"——后者已由 Cα RMSD（< 0.11 Å）排除。

### 2.4 HLA-A\*11:01 锚位筛选（9 体肽，覆盖位点 12）

`tools/mhc.py` 实现 A\*11:01 的 P2 / P9 锚位 motif 评分（A3 超型；保守的 motif 评分器，**不是** IC50 预测器；参考 Sidney 2008 *BMC Immunol* 9:1 / Falk 1994 *Immunogenetics* 40:238 / Zhang 1993 *PNAS* 90:2217）。9 个候选 9 体肽中通过的唯一一条：

> **VVGADGVGK**（KRAS 残基 8–16；P2 = V 偏好，P9 = K 偏好；与 WT 9 体不同）

这正是 **Wang QJ et al. 2016 *Cancer Immunology Research* 4:204–214** 首次确证的 HLA-A\*11:01 / KRAS G12D 经典 9 体新抗原——盲扫的 motif 评分器独立复现了该已发表表位，是流水线一致性的一项强 sanity check。

> **重要免责声明（详见 §4 caveat A）：** Zhu et al. 2026 *Commun Biol* 9:26 用 NMR + 晶体 + MD + TCR 结合实验系统刻画了 A\*11:01/KRAS-G12 复合物，结论是：A\*11:01 同时能呈递 9 体（VVGADGVGK）和 10 体（VVVGADGVGK），但**只有 10 体**形成稳定、TCR 可识别的构象（Asp12 在 p6 凸出，相关 PDB 7OW4）；9 体在凹槽内构象被"压扁"。本流水线的 motif 评分器构造上只看 9 体，因此可以正确**定位**到位点 12 上下的肽窗口，但**不能比较** 9 体与 10 体的相对优劣。

---

## 3. 流水线工程要点

| 维度 | 落实方式 |
|---|---|
| 双栈正交预测 | ColabFold AF2（MSA 模式 `mmseqs2_uniref_env`，5 模型，3 recycles，seed 42，model type `alphafold2_ptm`，rank by plddt）+ ESMFold v1（单序列、HF `facebook/esmfold_v1`，fp16 ESM-2 stem + fp32 折叠 trunk 混合精度，适配 WSL2 7.4 GiB 主机 RAM 上限） |
| brief-锁定参数 | 全部硬编码于 `tools/msa.LOCKED_FLAGS`，任何驱动脚本都不可修改；硬件调优参数 `--disable-unified-memory` 与模型本身正交 |
| 时间/重试预算 | `tools/retry.py`：60/300/900 s 指数退避，2 h SIGALRM 单次超时；预测产物若被中断，部分输出保留在 `logs/timeouts/<jobid>/` |
| 硬质量门槛 | 任一 rank-1 预测必须 mean pLDDT ≥ 85 才视为成功；4/4 通过 |
| 下载预算 | 累计 19.62 GB / 20 GB（见 `logs/downloads.md` 与 `logs/01_footprint.md`）；其中 ESMFold 权重 8.44 GB、AF2 参数 3.47 GB、pip wheels 7.70 GB、MSA 0.003 GB |
| 缓存隔离 | `env/activate_project.sh` 把 `PIP_CACHE_DIR / HF_HOME / TORCH_HOME / JAX_COMPILATION_CACHE_DIR / XDG_CACHE_HOME / TMPDIR / MPLCONFIGDIR` 全部重定向到 `./cache/`；未污染 `~/.cache`、`/tmp` |
| 复现性测试 | `tools/tests/` 43 个单测全绿（~3.4 s）；覆盖 FASTA I/O、HGVS 点突变应用、retry/SIGALRM、pLDDT 加载（AF2 + ESMFold）、Cα Kabsch RMSD、A\*11:01 motif 评分（含 VVGADGVGK 偏好双锚验证） |

修补的上游库（保留在 `./env` 内，附行号在 `logs/01_install.md` §5）：
- `alphafold/model/modules.py:1934` 与 `modules_multimer.py:521,543`：JAX 0.10 取消 `a_min/a_max` 别名，改写为 `min=/max=`。
- `transformers/models/esm/modeling_esmfold.py:2173`：防御性 `ptm_logits.float()` 升精，避免 fp16 softmax 进 NaN 而被 `compute_tm` 触发 IndexError。

---

## 4. 流水线**不**做的事（reader 必须看到）

以下九条 caveat 在 `logs/06_self_audit.md` §D 中逐条对应到论据/反证；这里给读者一份摘要：

**A. 9 体 vs 10 体（Zhu 2026 *Commun Biol* 9:26）：** 本流水线只筛 9 体；A\*11:01 上免疫生产性的 KRAS G12D 表位是 **10 体 VVVGADGVGK**。我们正确定位到该窗口，但不评分 10 体延伸。

**B. motif 评分器 ≠ IC50 预测器：** `tools/mhc.py` 返回锚位等级（preferred / tolerated / poor），不是 nM 级亲和力。brief 不允许装 NetMHCpan / MHCflurry。

**C. Tran 2016 *NEJM* 不是 A\*11:01 的论据：** 在 Step-5 文献核对中纠正了一处早期草稿错误。Tran 2016 用的是 **HLA-C\*08:02** 限制性 TIL 在转移性结直肠癌（非 PDAC）患者上的过继回输；其识别肽是 GADGVGKSA(L)，**不是** VVGADGVGK。A\*11:01 / VVGADGVGK 的首要论据是 Wang QJ 2016 *CIR* 4:204。

**D. Hunt 1992 *Science* 不是 A\*11:01 motif 的论据：** 在 Step-5 文献核对中纠正了 `tools/mhc.py` 文档字符串中的一处错误引用——Hunt 1992 报告的是 HLA-A2.1 的 motif，不是 A\*11:01。正确替换为 Zhang QJ 1993 *PNAS* 90:2217。

**E. heavy-atom-mean vs Cα-only pLDDT 约定差异：** ESMFold 用 heavy-atom mean per residue（HF 默认输出 37 原子 × L），AF2 的 per-residue pLDDT 在概念上更接近 Cα-only。两个约定在同一模型上差 3–6 pLDDT，因此 ≈86 vs ≈93 的绝对差高估了真实分歧。两者都以两种约定通过 ≥85 门槛。

**F. JAX RNG 在 Blackwell 上非确定性：** 即便 `--random-seed 42` 锁定，消费级 Blackwell（sm_120）的 JAX RNG 不是位级确定的（上游 JAX/XLA release notes），逐残基 pLDDT 可能在重运行间偏移 ~0.1。均值门槛的余量远大于此。

**G. WSL2 主机 RAM 上限 ≈ 7.4 GiB：** ESMFold 强制走 fp16 ESM-2 stem + fp32 trunk/head 的混合精度配方才能容纳。具体补丁见 §3。

**H. 未做 OpenMM/AMBER 松弛：** 会再下载 > 2 GB 突破 20 GB 预算；全程使用 ColabFold 的 `unrelaxed_rank_001_*.pdb`。由于本研究所有结论都建立在主链 Cα RMSD + per-residue pLDDT 上（而不是侧链 packing 能量），这一限制不影响结论。

**I. Rojas 2023 是个体化疫苗，不是固定 G12D/A\*11:01 panel：** autogene cevumeran (BNT122) 是按患者肿瘤特异性新抗原图谱定制的。本研究把 Rojas 2023 作为"PDAC mRNA 疫苗临床可行性"的**整体动因**引用，不据此声称 `VVGADGVGK` 被该 16 例患者中的某一具体患者使用过（那需要 trial 附件，不在主文中）。

---

## 5. 复现指南（快速）

```bash
source env/activate_project.sh        # 重定向所有缓存到 ./cache/
pytest                                # 43 个单测 ≈ 3 s
python -m tools.analyze_kras_step3    # Step 3 数值结果（基于既有预测）
python -m tools.make_step4_figures    # Step 4 四张论文级图
```

重跑结构预测（每条 5–7 min，GPU）：

```bash
python -m tools.run_kras_colabfold \
    --fasta data/sequences/kras_g12d_human_1-169.fasta \
    --out-dir colabfold/kras_g12d_1-169 \
    --log-path logs/colabfold/kras_g12d_1-169.log

python -m tools.esmfold_run \
    --fasta data/sequences/kras_g12d_human_1-169.fasta \
    --out-dir esmfold/kras_g12d_1-169
```

详见 `README.md` 与 `logs/00..06_*.md` 时序记录。

---

## 6. 阅读地图

| 我想知道… | 看 |
|---|---|
| 主机审计、CUDA/JAX/torch 版本、平台 caveat | `logs/00_environment.md` |
| 安装、smoke 测试、踩到的 5 个上游 bug 与修复 | `logs/01_install.md` |
| 20 GB 下载预算明细 | `logs/01_footprint.md`、`logs/downloads.md` |
| 工具模块清单、43 个单测覆盖矩阵 | `logs/02_helper_toolkit.md`、`tools/`、`tools/tests/` |
| 双栈预测、RMSD、pLDDT、9 体筛选——**核心结果** | `logs/03_core_workflow.md`、`analysis/*.json` |
| 4 张论文级图的设计、数值与来源溯源 | `logs/04_figures.md`、`figures/*.txt` 边注 |
| 文献核验日志（含三处更正） | `refs/notes.md`、`refs/citations.bib` |
| 反幻觉自审计（数值溯源 + 9 条 caveat） | `logs/06_self_audit.md` |
| 流水线总述 + 工程门槛兑现 | 本文件 |
