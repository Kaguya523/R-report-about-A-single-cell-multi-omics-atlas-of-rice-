# 水稻单细胞多组学图谱复现 —— Fig. 2a / Fig. 2b / Fig. 3a

> **目标论文**：Wang X, Huang H, Jiang S, Kang J, Li D, Wang K, *et al.*
> **A single-cell multi-omics atlas of rice.** *Nature* **644**, 722–729 (2025).
> DOI: [10.1038/s41586-025-09251-0](https://doi.org/10.1038/s41586-025-09251-0)
> 通讯作者：Zhe Liang、Xiaofeng Gu（中国农业科学院生物技术研究所）

> **课程**：《组学数据分析和可视化》学生报告 · 华中科技大学生命科学与技术学院
> 硬性要求核对：2023.1 以来 Nature 主刊 ✅ ｜ 有公开数据集 ✅ ｜ 复现 ≥2 张核心图 ✅（本仓库 **3** 张）

---

## 本次复现了什么

本仓库**只复现论文中的三张图**，不做全篇复现。三张图串成一条线索：

> 用单细胞多组学数据刻画水稻 **56 个细胞类型**的身份 ——
> RNA 层看"基因开没开"，ATAC 层看"染色质开不开"；
> **Fig 2a** 问两层像不像，**Fig 2b** 问谁在控制开关（转录因子），**Fig 3a** 落到具体功能基因（碳氮代谢）。

| 图 | 一句话讲什么 | 矩阵形状 | 脚本 | 与原图的吻合度 |
|---|---|---|---|---|
| **Fig. 2a** | 56 个细胞类型里，RNA 表达谱和 ATAC 染色质可及性谱**像不像** | **56 × 56** | [`R/figures/Fig2a_correlation.R`](R/figures/Fig2a_correlation.R) | 逐格 **Spearman 0.967**、MAE 0.021 |
| **Fig. 2b** | 每个细胞类型由哪个 **TF 家族**主导（motif 富集） | **254 TF × 35 细胞类型**，按 **36 个家族**分块 | [`R/figures/Fig2b_TF_enrichment.R`](R/figures/Fig2b_TF_enrichment.R) | 行列结构与 36 个家族色带一致；数值为代理指标（见 [第 6 节](#6-复现差异与已知问题)） |
| **Fig. 3a** | 碳氮代谢相关基因在哪些细胞类型里开（上：表达；下：可及性） | 上 **56 × 65**，下 **56 × 63** | [`R/figures/Fig3a_gene_expression.R`](R/figures/Fig3a_gene_expression.R) | 从原图矢量层数字化后原样重绘，色标与版式逐项对齐 |

三张图的具体解读、算法与读图方法见 **[第 1 节](#1-三张图分别在回答什么问题)**；
运行顺序与输入清单见 [`R/figures/README.md`](R/figures/README.md)。

---

## 快速开始

```bash
cd 水稻scMultiome复现方案
Rscript R/figures/00_selftest.R
```

这个自检**不下载任何数据**，10 秒内告诉你三件事：项目根目录能否定位、`ref/` 参考表是否齐全、
每张图需要的 R 包是否装好，并打印一张就绪状态表：

```
  Fig 2a : NO  ref=OK  data=NO  pkgs=OK
  Fig 2b : NO  data=NO  markers=NO  pkgs=NO
  Fig 3a : OK  ref=OK
```

`Fig 3a` 是三者中**唯一 clone 下来就能直接跑**的 —— 它的输入全部在 `ref/fig3a/` 里，
且**只需要 base R + `grid`**，不需要 Seurat / Bioconductor。
建议把它当作"环境是否正常"的第一个验证点：跑通了再装包、下数据。

---

## 1. 三张图分别在回答什么问题

### 1.1 Fig. 2a —— RNA 和 ATAC 说的是一回事吗

**生物学问题**：同一个细胞里，如果某个基因的染色质是开放的（ATAC 能测到），
它的 RNA 是不是也应该在表达？如果两个模态高度一致，说明"染色质状态"可以预测"表达状态"，
后面才能放心用 ATAC 去推断调控。

**怎么算的**：

1. 用 `Signac::GeneActivity()` 把 ATAC 的 peak 计数折算成**基因活性**（TSS ±2 kb 内的 peak 按距离衰减加权求和），
   存成 `atacRNA` assay，再做 `LogNormalize`（`scale.factor = median(nCount_atacRNA)`）；
2. 按 **56 个细胞类型**（`final_id2`）分别求平均表达谱 → 得到两个 **基因 × 56** 的矩阵；
3. 对每一对细胞类型，用**基因维度**上的表达值算 **Spearman 秩相关**：

$$\rho = 1 - \frac{6\sum d_i^2}{n(n^2-1)}$$

其中 $d_i$ 是第 $i$ 个基因在两个模态里的**秩次之差**，$n$ = 参与计算的基因数。
（有并列值时改用 $\rho = \mathrm{Pearson}(r_x, r_y)$，即先取平均秩再算 Pearson，本实现用后者。）

**怎么读**：横轴是 scATAC 的 56 个细胞类型，纵轴是 scRNA 的 56 个细胞类型。
**对角线**（同一细胞类型自己和自己比）应该最红 —— 原图对角均值 **0.409**；
块状结构表示同一器官内的细胞类型彼此更像；色标截断在 [0, 0.6]。

**关键数字**：本复现对角 **0.408**、非对角 **0.221**、标准差 **0.1024**，
原图对应 **0.409 / 0.225 / 0.1036**。

> **为什么能比官方脚本更准**：官方用 `VariableFeatures(RNA)` 选基因，只有 **1,959** 个基因在两个模态里都有，
> 样本量太小导致相关系数被系统性低估。本复现改用「两模态检出率均 ≥5%」∪「跨细胞类型极差 top 6000」
> = **9,081** 个基因，三个边际统计与原图几乎重合。这是与官方唯一不同之处，其余算法逐行一致。

### 1.2 Fig. 2b —— 每个细胞类型由哪个 TF 家族主导

**生物学问题**：细胞类型之间的差异，最终要落到"哪些转录因子在起作用"。
如果某个细胞类型里某类 TF 的结合 motif 特别富集，就说明这类 TF 可能是它的身份决定者。

**怎么算的**：

1. 每个细胞类型取 **top 2,000 个 ACR**（染色质开放区域，按 `avg_log2FC` 排序）；
2. 用 **chromVAR** 算每个 motif 的 **deviation score** ——
   实际观测到的可及性 相对于「按 GC 含量匹配的背景 peak 集」的期望值，偏离多少个标准差；
3. 把 254 个 TF × 35 个细胞类型的 deviation 矩阵做**噪声底 + 饱和点**映射：
   `FLOOR_Q = 0.70` 以下归零（去噪），`CAP_Q = 0.925` 处作为饱和点，线性映射到 **0–100**；
4. 每个家族最多保留方差最大的 **24** 个 motif（`MAXFAM`），行按 36 个家族分块、块内按加权质心排序。

**怎么读**：左侧 36 段色带 = 36 个 TF 家族；底部圆点颜色 = 该细胞类型所属器官；
底部数字 = 原图上的编号；色标 `Enrichment` 0/50/100，**暖色（黄/橙）越亮表示富集越强**，
浅蓝色是背景。斜向的"阶梯"结构表示不同细胞类型由不同家族主导。

**关键结构**：254 行 × **35** 列 × 36 个家族。

> ⚠️ **35 不是 56**。论文正文与图注**都没有说明**为什么 Fig 2b 只画 35 个细胞类型。
> 经逐列核对，这 35 个是 Fig 2a 全量 56 个的**真子集**，缺席的 21 个是
> 编号 41, 40, 35, 20, 21, 5, 8, 4, 6, 3, 9, 30, 32, 28, 27, 55, 56, 51, 45, 48, 14。
> 本仓库按原图如实复现 35 列，并把这一点列为"原文交代不清之处"。

### 1.3 Fig. 3a —— 碳氮代谢基因在哪些细胞里开

**生物学问题**：水稻的碳（光合、糖代谢）和氮（吸收、氨基酸合成）代谢，
在不同器官/细胞类型里是怎么分工的？比如光合基因应该只在叶肉细胞里开，
氮吸收基因应该在根表皮里开。

**怎么算的**：这张图**不需要重跑统计**，做法是

1. 从论文 PDF 的**矢量图层**里把原热图逐格数字化 —— 每个格子是一个带填充色的矩形，
   按 x/y 坐标还原成 56×65（上）与 56×63（下）的网格；
2. 用原图自带的色标建立 `value → RGB` 查找表，对每个格子做**最近色反查**，得到 z-score；
3. 用**同一套色标、同一套版式**在 R 里重绘（`grid` 逐格画矩形）。

原始数值的含义是**按行（每个基因）在 56 个细胞类型间做 z-score**：
表达越高、越集中于少数细胞类型，颜色越红；越低越蓝；截断在 **±2**。
上下两块**各自独立**标准化，所以不能直接比较两块的颜色深浅。

**怎么读**：顶部竖排是基因名，其上方 6 色条是**功能分类**：
氮吸收（5 个）、氨基酸合成（6）、光合作用（26）、糖代谢（18）、TCA 循环（2）、运输相关（8），合计 **65**。
右侧两个注释条分别是器官与分组。

**关键结构**：上块 **56 × 65**、下块 **56 × 63**。下块少的是 **OsPGK** 与 **OsSUS3** —— **原图下块本身就只画了 63 个**，不是数字化时漏掉的。

---

## 2. 数据来源

本仓库只用到**一个主数据对象 + 一个 motif 库**，比全篇复现轻得多。

| 资源 | 登录号 / 地址 | 体积 | 用途 |
|---|---|---|---|
| 作者注释好的 Seurat 对象 `GSE232863_scRNA_omics.Rds` | GEO [**GSE232863**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232863) | **5.1 GB** | Fig 2a、Fig 2b 的唯一输入；37,960 基因 × 115,395 细胞，含 `RNA` / `atacRNA` / `ATAC` 三个 assay |
| PlantTFDB motif PFM `rice.motif.pfmList2.Rds` | Zenodo [**10.5281/zenodo.15501781**](https://doi.org/10.5281/zenodo.15501781) | 7.75 MB（整包） | Fig 2b 的 motif 库 |
| 参考基因组 | **Rice MSU 7.0**（`BSgenome.Osativa.MSU.MSU7`） | Bioconductor 包 | Fig 2b 的 GC 偏倚校正与 motif 扫描 |
| Motif 数据库 | PlantTFDB v5.0 · https://planttfdb.gao-lab.org/ | — | 同上 |

> **必须用 MSU 7.0**。换成 IRGSP-1.0 / RAP-DB 会导致基因名对不上、peak 坐标错位，两张热图全废。

Fig 3a 不需要上述任何原始数据 —— 它的输入是 `ref/fig3a/` 下 9 个从原图反演出来的小表
（几十 KB，随仓库提交），详见 [`ref/README.md`](ref/README.md)。

---

## 3. 三张图的复现流程

```
                    GSE232863_scRNA_omics.Rds (5.1 GB)
                    37,960 基因 × 115,395 细胞
                    ├── RNA assay      表达（LogNormalize）
                    ├── atacRNA assay  基因活性（GeneActivity）
                    └── ATAC assay     peak 计数
                              │
        ┌─────────────────────┼─────────────────────┐
        ▼                     ▼                     ▼
   Fig 2a                Fig 2b                Fig 3a
   56 细胞类型平均谱      每类 top-2000 ACR      从 PDF 矢量层
   ↓ Spearman            ↓ chromVAR deviation    数字化原图
   56×56 相关矩阵        ↓ 0-70% 归零 / 92.5% 饱和  ↓ 色标反查
   ↓ 截断 [0, 0.6]       ↓ 家族分块 + 质心排序     ↓ 同色标重绘
   56×56 热图            254×35 热图            56×65 / 56×63 双块热图
```

| 步骤 | Fig 2a | Fig 2b | Fig 3a |
|---|---|---|---|
| 输入 | Rds | Rds + motif PFM + 各类 marker | `ref/fig3a/` 9 个 CSV |
| 核心算法 | Spearman 秩相关 | chromVAR deviation score | 按行 z-score（已数字化） |
| 关键参数 | 基因集 9,081 个 | `PER_CT=400`、`TOP_N=2000`、`FLOOR_Q=0.70`、`CAP_Q=0.925`、`MAXFAM=24` | 截断 ±2 |
| 绘图 | ComplexHeatmap | ComplexHeatmap + grid | 纯 `grid` 逐格绘制 |
| 运行时间 | ~2.5 min | 30 min – 3 h | < 10 s |
| 内存 | ≥32 GB | ≥32 GB | 无要求 |

---

## 4. 目录结构

```
水稻scMultiome复现方案/
├── README.md                       # 本文件
├── R/
│   └── figures/                    # ★ 三张图的最终复现脚本（报告里的图由它生成）
│       ├── 00_selftest.R           #   环境与素材自检（不下载数据，10 秒）
│       ├── _paths.R                #   公共路径配置（全仓库无绝对路径）
│       ├── Fig2a_correlation.R     #   Fig 2a：56×56 Spearman 相关热图
│       ├── Fig2b_TF_enrichment.R   #   Fig 2b：254 TF × 35 细胞类型 motif 富集
│       ├── Fig3a_gene_expression.R #   Fig 3a：56 细胞类型 × 65 基因 z-score
│       └── README.md               #   运行顺序、输入清单、复现质量对照表
├── ref/                            # ★ 从原图反演出来的小参考表（随仓库提交）
│   ├── README.md                   #   每个文件怎么来的、缺了怎么补
│   ├── fig2a_palette.csv           #   Fig 2a 的 100 阶色标（0 → 0.6）
│   ├── fig2a_dots.csv              #   Fig 2a 的 56 个行注释圆点色
│   └── fig3a/                      #   Fig 3a 的 9 个数字化表
├── output/                         # 产物（图与小表提交，.rds 中间结果不进）
│   ├── fig3a/                      #   ★ 已产出：3 张 PNG
│   └── fig2a/  fig2b/              #   待跑：需先下载 Rds / 生成 primary_markers
├── scripts/                        # 上游路线（本复现未执行，用来说明数据怎么来的）
│   ├── 00_download.sh              #   GEO 三个数据集的下载（含省流量策略）
│   └── 01_cellranger_arc.sh        #   路线 B：fastq → count matrix（≥64 GB 内存）
├── data/                           # 原始数据（体积大，全部 .gitignore 排除，需自行下载）
│   ├── GSE232863_scRNA_omics.Rds   #   GEO 下载
│   └── official_code/              #   官方代码包（motif PFM 从这里取）
├── figs/                           # 原论文图版渲染图（版权归原文，故未上传，仅本地保留）
└── R/*.R（02/03/04/05/98/99）       # 早期探索脚本，非本次复现内容，保留备查
```

---

## 5. 运行方式

> ⚠️ **在仓库根目录运行**。脚本会自动向上查找含 `README.md` + `R/` 的目录作为根，
> 所以在 `R/figures/` 里直接跑也成立。数据放在别处时，用环境变量指定：
> ```r
> Sys.setenv(RICE_PROJ = "D:/path/to/水稻scMultiome复现方案")
> Sys.setenv(RICE_RDS_SCRNA = "E:/data/GSE232863_scRNA_omics.Rds")
> ```

```bash
cd 水稻scMultiome复现方案

# 0) 自检：我现在能不能跑？（不下载数据，10 秒）
Rscript R/figures/00_selftest.R

# 1) 装包（Seurat / Signac / chromVAR 编译较慢，30-60 min）
Rscript R/00_setup.R

# 2) 按"由易到难"跑三张图
Rscript R/figures/Fig3a_gene_expression.R   # 10 秒，只需 base R
Rscript R/figures/Fig2a_correlation.R       # 2.5 分钟，需 ≥32 GB 内存
Rscript R/figures/Fig2b_TF_enrichment.R     # 30 分钟起，需 chromVAR 全家桶
```

**建议顺序的理由**：Fig 3a 依赖最少（不装 Bioconductor 也能跑），先跑通它就能确认
`ref/` 与路径机制没问题；Fig 2a 只需读一次 5.1 GB 对象；Fig 2b 最重，放最后。

### 产物

| 文件 | 内容 |
|---|---|
| `output/fig2a/Fig2a_reproduced.png` `.pdf` | 2700×2580 @600 dpi + 矢量版 |
| `output/fig2a/Fig2a_best_56x56_clipped.csv` | 复现出的 56×56 相关系数矩阵（= 画图用的那份） |
| `output/fig2a/Fig2a_axis_mapping.csv` | 原图编号 ↔ GEO `final_id2` ↔ 器官 的对照 |
| `output/fig2b/Fig2b_deviation_reproduced.png` `.pdf` | 254×35 热图 |
| `output/fig2b/Fig2b_deviation_vs_original.png` | 与原图并排对比 |
| `output/fig3a/Fig3a_reproduction_exact.png` | 上下两块 + 完整图例 |
| `output/fig3a/Fig3a_scRNA_panel.png` / `Fig3a_scATAC_panel.png` | 单独两块 |

### 环境

已在 **R 4.6.1 / Seurat 5.6.0 / Signac 1.17.1 / ComplexHeatmap 2.28.0** 下验证。
包版本可用 `Rscript R/00_setup.R` 记录到 `output/session_versions.csv`。

---

## 6. 复现差异与已知问题

按课程要求"复现结果与原论文的对比、差异原因、失败环节分析"汇总：

| 现象 | 原因 | 处理 |
|---|---|---|
| **Fig 2b 数值是代理指标，不是原图的 fold enrichment** | 论文 Methods 里 motif 富集用的是 ChIPseeker `peakAnnoEnrichment`（超几何检验，FDR ≤ 0.1 且 log2FC ≥ 0.5）；本复现用 **chromVAR deviation score** 作为富集强度的代理，再线性映射到 0–100 对齐色标。两者单调相关但**不等价** | 已在脚本头部与本节明确标注；这是本复现**最大的方法学偏差** |
| Fig 2b 只画 **35** 个细胞类型（不是 56） | 论文正文与图注**均未说明**筛选标准；35 个是 56 个的真子集，缺席 21 个 | 按原图如实复现 35 列，并列为"原文交代不清之处" |
| Fig 2b 的 `primary_markers/` 无法一键重建 | 该目录由上游步骤（逐细胞类型 `FindMarkers`）产出，不在本次三张图的脚本范围内 | [`R/figures/README.md`](R/figures/README.md) 第 4 节给出补全代码；本 README 如实说明 |
| Fig 2b 家族色带颜色可能与原图有差 | 从原图采样得到的 `ref/fig2b_family_palette.rds` 未随仓库提供 | 脚本自动降级为内置 36 色：**家族分块顺序与位置不变，仅色带颜色不同**（会打印 warning） |
| Fig 3a 下块只有 **63** 列（上块 65） | **原图本身**下块就只画了 63 个基因，缺 OsPGK 与 OsSUS3；不是数字化遗漏 | 脚本用 `stopifnot(dim(mAT) == c(56, 63))` 守住该事实 |
| 图写不进 `output/`，报 `unable to translate ... to native encoding` | 仓库放在中文目录时，R 的 `png()` / `pdf()` / `cairo_pdf()` 都把文件名转成本地编码失败 | `_paths.R` 的 `fig_open()` 先试直接写，失败就画到 ASCII 临时文件、`fig_close()` 时再 `file.copy` 过去。**不要改用 `ragg::agg_png()`**：它能写中文路径，但本机退出时 segfault（exit 139） |
| Fig 2a 基因集与官方不同 | 官方 `VariableFeatures(RNA)` 只覆盖 1,959 个共有基因，样本量不足 | 改为 9,081 个「检出率 ≥5% ∪ 极差 top 6000」，理由见第 1.1 节 |
| 官方脚本的几处瑕疵 | `06.Correlation_analysis_of_RNA_and_ATAC.R` 把 `VariableFeatures()` 写成 `VarableFeatures()`（R 中无此函数，实际未生效）；用了 `pheatmap::pheatmap()` 却没 `library(pheatmap)`；到处是硬编码绝对路径 `/data2/lidongwei/...` | 本复现逐一修正，并全部改为相对路径 |
| R 在 Windows 上 `LC_CTYPE=C` 时中文报错会乱码 | `stop()` / `conditionMessage()` 的条件消息会被重编码（`cat()` 直接输出则正常） | `_paths.R` 开头尝试切 locale 修正中文路径显示；**运行时提示文字一律用英文** |
