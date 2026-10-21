# 水稻单细胞多组学图谱复现 —— Fig. 1b / 1c / 2a / 2b

> **目标论文**：Wang X, Huang H, Jiang S, Kang J, Li D, Wang K, *et al.*
> **A single-cell multi-omics atlas of rice.** *Nature* **644**, 722–729 (2025).
> DOI: [10.1038/s41586-025-09251-0](https://doi.org/10.1038/s41586-025-09251-0)
> 通讯作者：Zhe Liang、Xiaofeng Gu（中国农业科学院生物技术研究所）

> **课程**：《组学数据分析和可视化》学生报告 · 华中科技大学生命科学与技术学院
> 硬性要求核对：2023.1 以来 Nature 主刊 ✅ ｜ 有公开数据集 ✅ ｜ 复现 ≥2 张核心图 ✅（本仓库 **4** 张）

---

## 本次复现了什么

本仓库**只复现论文中的四张图**，不做全篇复现。四张图串成一条线索：

> **Fig 1b** 先把 115,395 个细胞核铺在 UMAP 上，看 8 个器官怎么分布 ——
> 结果发现"细胞类型相似性"比"器官身份"更能决定位置；
> **Fig 1c** 顺着这个线索，用多组学当标准答案，量化"分开测再拼"到底有多不准（**37%–82%**）；
> **Fig 2a** 问 RNA 与 ATAC 在细胞类型层面像不像；
> **Fig 2b** 再问一句：谁在控制这些开关（转录因子 motif）。

| 图 | 一句话讲什么 | 矩阵形状 | 脚本 | 与原图的吻合度 |
|---|---|---|---|---|
| **Fig. 1b** | 8 个器官的细胞在 UMAP 上怎么分布，各有多少 | 115,395 点 × 2 维 + 8 根条 | [`R/figures/Fig1b_reproduce.R`](R/figures/Fig1b_reproduce.R) | 版式与配色一致；细胞数差 **1.0%**（口径差异，见 [第 6 节](#6-复现差异与已知问题)） |
| **Fig. 1c** | 把 RNA 与 ATAC 拆成"两批分开测的数据"再拼，能拼多准 | 8 器官 × 2 根横条 | [`R/figures/Fig1c_reproduce.R`](R/figures/Fig1c_reproduce.R) | **定性一致**：冠根最差、剑叶最好；绝对值受口径影响，不可逐项比 |
| **Fig. 2a** | 56 个细胞类型里，RNA 表达谱和 ATAC 染色质可及性谱**像不像** | **56 × 56** | [`R/figures/Fig2a_correlation.R`](R/figures/Fig2a_correlation.R) | 逐格 **Spearman 0.967**、MAE 0.021 |
| **Fig. 2b** | 每个细胞类型由哪个 **TF 家族**主导（motif 富集） | **254 TF × 35 细胞类型**，按 **36 个家族**分块 | [`R/figures/Fig2b_TF_enrichment.R`](R/figures/Fig2b_TF_enrichment.R) | 行列结构与 36 个家族色带一致；数值为代理指标（见 [第 6 节](#6-复现差异与已知问题)） |

四张图的具体解读、算法与读图方法见 **[第 1 节](#1-四张图分别在回答什么问题)**；
运行顺序与输入清单见 [`R/figures/README.md`](R/figures/README.md)。

> **交付文档**：`reports/Fig1b复现报告.docx`、`reports/Fig1c复现报告_与原论文对比_差异归因与失败环节分析.docx`

---

## 快速开始

```bash
cd 水稻scMultiome复现方案
Rscript R/figures/00_selftest.R
```

这个自检**不下载任何数据**，10 秒内告诉你四件事：项目根目录能否定位、`ref/` 参考表是否齐全、
主数据对象是否就位、每张图需要的 R 包是否装好，并打印一张就绪状态表：

```
  Fig 1b : OK  data=OK  pkgs=OK
  Fig 1c : OK  data=OK  pkgs=OK
  Fig 2a : NO  ref=OK  data=NO  pkgs=OK
  Fig 2b : NO  data=NO  markers=NO  pkgs=NO
```

四张图都依赖同一个主数据对象 `GSE232863_scRNA_omics.Rds`（5.1 GB），
所以第一步永远是把它放进 `data/`，详见 [第 2 节](#2-数据来源)。

---

## 1. 四张图分别在回答什么问题

### 1.1 Fig. 1b —— 8 个器官在 UMAP 上怎么分布

**生物学问题**：不同器官的细胞，在转录组 + 染色质联合空间里是各自成岛，还是混在一起？

**怎么算的**：

1. 取作者发布的 **WNN 降维坐标**（`wnn.umap`）—— 由 RNA 的 PCA 前 50 维与 ATAC 的 LSI 第 2–50 维
   做 `FindMultiModalNeighbors`，再用 `RunUMAP(n.neighbors=30, n.epochs=200, min.dist=0.5, a=0.03, b=1.5)` 得到；
2. 按 `tissue` 列（8 个器官）着色；
3. 按器官计数，画成右栏的横向条形图。

**关键点：器官标签来自实验设计，不是聚类结果。** 8 个器官分别取材、分别建库，
细胞身份由建库时的样本来源直接确定。真正由算法产生的是 **56 个细胞类型簇**（`FindClusters`，resolution 0.1），
注释后得到 54 个类型 —— 那是 Fig 1d 的内容。

**怎么读**：图上**不是 8 个分离的岛**。原文 p60-63 指出：不同器官的细胞**倾向于**分开，
但生理相关的器官彼此靠近甚至部分重叠；跨器官共有的细胞类型（如叶肉，簇 19 与 25）也出现重叠，
说明 **"细胞类型相似性"是比"器官身份"更强的因素**。这直接解释了 Fig 1c 里为什么会判错。

**复现边界**：Fig 1b 属于**可视化层面复现** —— 坐标与器官注释均为作者给定，
复现目标是版式、配色与细胞数比例。脚本另提供 `WNN_SOURCE = "recompute"` 用于一致性检验。

### 1.2 Fig. 1c —— 分开测再拼，能拼多准

**生物学问题**：scRNA-seq 和 scATAC-seq 通常是**分开测两批细胞**，再靠算法配对。
这套做法成立的前提是"两个模态一致"。这个前提站得住吗？

**怎么算的**（原文 Methods "Comparison of multi-omics and paired data"）：

1. **标准答案**：多组学同一个核同时测了 RNA 与 ATAC，注释标签（`tissue_cluster`）是真值；
2. **伪造两批数据**：把同一批细胞的 `RNA` 矩阵当"只测了 RNA 的那批"，
   把 `atacRNA` 矩阵（ATAC peak 按启动子 + 基因体汇总出的**基因活性**）当"只测了 ATAC 的那批"；
3. 两边各自 `NormalizeData → FindVariableFeatures → ScaleData → RunPCA`；
4. `FindTransferAnchors` + `TransferData` 把 RNA 的标签迁移到 ATAC 侧；
5. 与真值**逐细胞比对**，`预测 == 真值` 的比例就是 accuracy。

**怎么读**：8 个器官各两根横条 —— **实心 = 计算配对准确率**，**斜纹 = 多组学（恒 100%，是天花板不是方法）**，
虚线 = 均值。原图数值：冠根 **37.1%** 最低，剑叶 **82.4%** 最高。

**为什么低**：三个机制叠加 ——
① **时间差**：染色质先开、RNA 后起（过渡态细胞）；
② **转录后调控**：RNA 水平还受 mRNA 周转影响，且这种影响有很强的细胞类型特异性；
③ **测序深度差**：原文指出冠根深度最低、剑叶最高，准确率与深度正相关。
另外，"类型之间越像"越容易串（Fig 2a 的非对角高值就是证据）。

**关键结构**：**斜纹条不是复现出来的方法**，它是"自己跟自己比"，数学上必然 100%，只作参照。

### 1.3 Fig. 2a —— RNA 和 ATAC 说的是一回事吗

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

**非对角才是信息量所在**：非对角高 = A 类型的 RNA 谱与 B 类型的 ATAC 谱也像，
说明这两类在某一模态上难以区分 —— 这正是 Fig 1c 判错的结构性原因。
注意色标截断到 0–0.6 的代价：>0.6 全部饱和、负相关被压成 0，这张图读不出"反向"。

**关键数字**：本复现对角 **0.408**、非对角 **0.221**、标准差 **0.1024**，
原图对应 **0.409 / 0.225 / 0.1036**。

> **为什么能比官方脚本更准**：官方用 `VariableFeatures(RNA)` 选基因，只有 **1,959** 个基因在两个模态里都有，
> 样本量太小导致相关系数被系统性低估。本复现改用「两模态检出率均 ≥5%」∪「跨细胞类型极差 top 6000」
> = **9,081** 个基因，三个边际统计与原图几乎重合。这是与官方唯一不同之处，其余算法逐行一致。

### 1.4 Fig. 2b —— 每个细胞类型由哪个 TF 家族主导

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

> **这张图不回答"哪个 TF 调控哪个基因"**。它的输出是 TF × 细胞类型矩阵，没有基因维度；
> "TF → 靶基因"是后续 pySCENIC regulon（251 个）与 CellOracle 的工作，不在本仓库复现范围内。

---

## 2. 数据来源

本仓库只用到**一个主数据对象 + 一个 motif 库**，比全篇复现轻得多。

| 资源 | 登录号 / 地址 | 体积 | 用途 |
|---|---|---|---|
| 作者注释好的 Seurat 对象 `GSE232863_scRNA_omics.Rds` | GEO [**GSE232863**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232863) | **5.1 GB** | **四张图共同**的唯一主输入；37,960 基因 × 115,395 细胞，含 `RNA` / `atacRNA` / `ATAC` 三个 assay |
| PlantTFDB motif PFM `rice.motif.pfmList2.Rds` | Zenodo [**10.5281/zenodo.15501781**](https://doi.org/10.5281/zenodo.15501781) | 7.75 MB（整包） | Fig 2b 的 motif 库 |
| 参考基因组 | **Rice MSU 7.0**（`BSgenome.Osativa.MSU.MSU7`） | Bioconductor 包 | Fig 2b 的 GC 偏倚校正与 motif 扫描 |
| Motif 数据库 | PlantTFDB v5.0 · https://planttfdb.gao-lab.org/ | — | 同上 |

> **必须用 MSU 7.0**。换成 IRGSP-1.0 / RAP-DB 会导致基因名对不上、peak 坐标错位，两张热图全废。

**四张图各自用到对象里的哪部分**：

| | 用到的 assay | 用到的 metadata |
|---|---|---|
| Fig 1b | **不用任何 assay** | `wnn.umap` 坐标 + `tissue` |
| Fig 1c | `RNA`（参考集）+ `atacRNA`（查询集） | `tissue_cluster`（真值）+ `tissue`（分批） |
| Fig 2a | `RNA` + `atacRNA` | `final_id2`（56 个细胞类型） |
| Fig 2b | `ATAC`（peak 矩阵）+ motif PFM | `tissue_cluster` |

> `atacRNA` 是作者已算好的**基因活性矩阵**（`GeneActivity`，upstream 2 kb / downstream 0），
> 不需要重跑。Fig 1c 与 Fig 2a 都直接读它。

---

## 3. 四张图的复现流程

```
                    GSE232863_scRNA_omics.Rds (5.1 GB)
                    37,960 基因 × 115,395 细胞
                    ├── RNA assay      表达（LogNormalize）
                    ├── atacRNA assay  基因活性（GeneActivity）
                    ├── ATAC assay     peak 计数（82,219 峰）
                    └── wnn.umap       WNN 联合降维坐标
                              │
      ┌───────────┬───────────┼───────────┬───────────┐
      ▼           ▼           ▼           ▼           ▼
   Fig 1b      Fig 1c      Fig 2a      Fig 2b
   读坐标+标签   RNA→atacRNA  56 类平均谱  每类 top-2000 ACR
   ↓ 按器官着色  ↓ 锚点+迁移  ↓ Spearman   ↓ chromVAR deviation
   ↓ 按器官计数  ↓ 逐细胞比对 ↓ 截断[0,0.6] ↓ 0-70% 归零 / 92.5% 饱和
   UMAP+条形图   8×2 横条     56×56 热图   254×35 热图（36 家族分块）
```

| 步骤 | Fig 1b | Fig 1c | Fig 2a | Fig 2b |
|---|---|---|---|---|
| 输入 | Rds（坐标 + 标签） | Rds（RNA + atacRNA） | Rds | Rds + motif PFM + 各类 marker |
| 核心算法 | 直接读取 WNN-UMAP | `FindTransferAnchors` + `TransferData` | Spearman 秩相关 | chromVAR deviation score |
| 关键参数 | `n.neighbors=30, n.epochs=200, min.dist=0.5` | `nfeatures=2000, npcs=30, reduction=cca, k.anchor=5` | 基因集 9,081 个 | `PER_CT=400`、`TOP_N=2000`、`FLOOR_Q=0.70`、`CAP_Q=0.925`、`MAXFAM=24` |
| 绘图 | ggplot2 + patchwork | ggplot2 + ggpattern | ComplexHeatmap | ComplexHeatmap + grid |
| 运行时间 | ~6 min | 30 min – 2 h | ~2.5 min | 30 min – 3 h |
| 内存 | ≥16 GB | ≥16 GB | ≥32 GB | ≥32 GB |

---

## 4. 目录结构

```
水稻scMultiome复现方案/
├── README.md                       # 本文件
├── R/
│   └── figures/                    # ★ 四张图的最终复现脚本（报告里的图由它生成）
│       ├── 00_selftest.R           #   环境与素材自检（不下载数据，10 秒）
│       ├── 00b_make_primary_markers.R  # Fig 2b 前置：逐类 ATAC marker 表
│       ├── _paths.R                #   公共路径配置（Fig2a/2b 用；全仓库无绝对路径）
│       ├── Fig1b_reproduce.R       #   Fig 1b：WNN-UMAP + 各器官细胞数
│       ├── Fig1c_reproduce.R       #   Fig 1c：RNA→ATAC 标签迁移的逐器官准确率
│       ├── Fig2a_correlation.R     #   Fig 2a：56×56 Spearman 相关热图
│       ├── Fig2b_TF_enrichment.R   #   Fig 2b：254 TF × 35 细胞类型 motif 富集
│       └── README.md               #   运行顺序、输入清单、复现质量对照表
├── reports/                        # ★ 复现报告（Word）
│   ├── Fig1b复现报告.docx
│   └── Fig1c复现报告_与原论文对比_差异归因与失败环节分析.docx
├── ref/                            # ★ 从原图反演出来的小参考表（随仓库提交）
│   ├── README.md                   #   每个文件怎么来的、缺了怎么补
│   ├── fig2a_palette.csv           #   Fig 2a 的 100 阶色标（0 → 0.6）
│   └── fig2a_dots.csv              #   Fig 2a 的 56 个行注释圆点色
├── output/                         # 产物（图与小表提交，.rds 中间结果不进）
│   ├── fig2a/                      #   ★ Fig2a_vs_original.jpg（与原图并排）
│   └── fig2b/                      #   ★ Fig2b_vs_original.png（与原图并排）
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

> ⚠️ **在仓库根目录运行**。Fig2a / Fig2b 会自动向上查找含 `README.md` + `R/` 的目录作为根；
> Fig1b / Fig1c 默认读仓库相对路径。数据放在别处时，用环境变量指定：
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

# 2) 按"由易到难"跑四张图
Rscript R/figures/Fig1b_reproduce.R         # ~6 min，只读坐标与标签
Rscript R/figures/Fig1c_reproduce.R         # 30 min 起，跑锚点
Rscript R/figures/Fig2a_correlation.R       # 2.5 分钟，需 ≥32 GB 内存
Rscript R/figures/00b_make_primary_markers.R  # Fig 2b 前置：30 min – 1 h
Rscript R/figures/Fig2b_TF_enrichment.R     # 30 分钟起，需 chromVAR 全家桶
```

**建议顺序的理由**：Fig 1b 最轻（不碰表达值），先跑通它能确认数据与路径没问题；
Fig 1c 最耗时（逐器官跑锚点）；Fig 2a 只需读一次 5.1 GB 对象；Fig 2b 最重，放最后。

### 产物

| 文件 | 内容 | 状态 |
|---|---|---|
| `reports/Fig1b复现报告.docx` | Fig 1b 的复现说明与差异分析 | ✅ 已提交 |
| `reports/Fig1c复现报告_与原论文对比_差异归因与失败环节分析.docx` | Fig 1c 的差异归因与失败环节 | ✅ 已提交 |
| `output/fig2a/Fig2a_vs_original.jpg` | 与原图并排对比（56×56，Spearman 0.9668 / Pearson 0.9706 / MAE 0.0211） | ✅ 已提交 |
| `output/fig2b/Fig2b_vs_original.png` | 与原图并排对比（254×54 vs 253×56，口径见第 6 节） | ✅ 已提交 |
| `output/fig1b/Fig1b_reproduced.png` `.pdf` | UMAP + 细胞数条形图 | 跑 Fig1b 脚本后生成 |
| `output/fig1b/Fig1b_cell_counts.csv` | 8 器官细胞数 + 与原文差值 | 同上 |
| `output/fig1c/results/Fig1c_reproduced.png` `.pdf` | 8 器官 × 2 根横条 | 跑 Fig1c 脚本后生成 |
| `output/fig1c/results/Fig1c_confusion_<organ>.csv` | 逐器官混淆矩阵 | 同上 |
| `output/fig2a/Fig2a_reproduced.png` `.pdf` | 2700×2580 @600dpi + 矢量版 | 跑 Fig2a 脚本后生成 |
| `output/fig2b/Fig2b_deviation_reproduced.png` `.pdf` | 254×35 热图（当前脚本口径） | 跑 Fig2b 脚本后生成 |

> **Fig 2a / 2b 的成品对比图已提交**；Fig 1b / 1c 的图与报告已提交（报告内含结果图），
> 脚本本身可重跑但需要 5.1 GB 主对象与相应 R 包。

### 环境

已在 **R 4.6.1 / Seurat 5.6.0 / Signac 1.17.1 / ComplexHeatmap 2.28.0** 下验证。
包版本可用 `Rscript R/00_setup.R` 记录到 `output/session_versions.csv`。

---

## 6. 复现差异与已知问题

按课程要求"复现结果与原论文的对比、差异原因、失败环节分析"汇总：

| 现象 | 原因 | 处理 |
|---|---|---|
| **Fig 1b 细胞数 115,395 vs 原文 116,564（差 1.0%）** | 正文报的是 QC 后最终纳入分析的数；提交的 Rds 是另一过滤版本。去 doublet（DoubletFinder）是随机算法，不同版本/seed 会去掉不同细胞 | 报告里以**器官间相对比例**对比（占比差 <1 个百分点），不用绝对数比较；`Fig1b_cell_counts.csv` 列出逐项差值 |
| **Fig 1c 的绝对值不可与原文逐项比较** | 三处口径差异：① 图注写"迁移到 ATAC **簇**"、Methods 写"迁移到 **细胞**"，本复现取细胞级；② 参考集范围（逐器官 subset vs 全局）；③ 标签体系（`tissue_cluster` 35 类 vs `wsnn_res.0.1` 56 簇） | 只比**器官间相对高低顺序**（剑叶 > 穗 > 分蘖芽 > 茎 > 种子 > 茎尖 > 叶 > 冠根），不比绝对值 |
| **Fig 1c 中叶（Leaf）复现值偏高约 +33.7 个百分点** | 逐器官建参考集在**结构上排除了跨器官同源类型的混淆**：原文指出簇 19/25 是"叶与剑叶共有的叶肉"、在 UMAP 上部分重叠；逐器官 subset 后，"叶的叶肉被认成剑叶叶肉"这条错误路径不存在 | **这是口径放宽导致的虚高，不是复现优于原文**。已如实写进 `reports/Fig1c复现报告…docx`；脚本可通过 `GLOBAL` 开关切换为全细胞参考集做对照 |
| **Fig 1c 中冠根 / 茎尖是真失败** | 原文准确率 37.1% / 48.8%，多数类基线约 27.3% / 30.1%，仅高出约 10–19 个百分点。原文给了两个并列解释：测序深度最低（冠根）与发育动态强 | 复现忠实重现了"该方法在这两个器官上失效"这一结论，**这本身就是复现成功的证据** |
| **Fig 2b 数值是代理指标，不是原图的 fold enrichment** | 原文对 Fig 2b 的统计口径**三处不一致**：图注写 `average deviation scores`（chromVAR）、正文写 Fisher 精确检验（FDR<0.01，平均 17 个 motif）、Methods 写 ChIPseeker `peakAnnoEnrichment`（FDR≤0.1 且 log2FC≥0.5） | 按**图注**采用 chromVAR deviation score（唯一与热图连续颜色对应的口径），再线性映射到 0–100 对齐色标。三处不一致已写进本节与报告 |
| **Fig 2b 选 peak 规则与官方脚本不一致** | 正文与图注称 `top 2,000 differential ACRs`；官方 `07.Analyzing_the_TF_and_motif_by_Signac.R` 实际用 `p<0.05 & log2FC>0.5` 阈值筛选，不足 10 个时退回 `top_n(300)`，**代码中没有 2000** | 本复现按正文采用 `TOP_N=2000`，并保留该参数便于切换；该差异影响单个 TF 的分数绝对值，不改变家族分块的整体模式 |
| Fig 2b 只画 **35** 个细胞类型（不是 56） | 论文正文与图注**均未说明**筛选标准；35 个是 56 个的真子集，缺席 21 个 | 按原图如实复现 35 列，并列为"原文交代不清之处" |
| **`output/fig2b/Fig2b_vs_original.png` 是早期版本的对比图** | 该图生成时按 56 个细胞类型（原图侧标 54）绘制、行数 253；当前脚本已收窄到原图口径 **254 TF × 35 细胞类型**。统计量与配色相同，差别只在列的取舍 | 图内标注保持原样以便追溯版本；重跑 `Fig2b_TF_enrichment.R` 即得 254×35 版，可直接替换该文件 |
| Fig 2b 的 `primary_markers/` 需要先跑一次前置步骤 | 该目录由逐细胞类型 `FindMarkers` 产出，与画图无关且耗时较长，故拆成独立脚本 | 已补 [`R/figures/00b_make_primary_markers.R`](R/figures/00b_make_primary_markers.R)，参数逐行对齐官方 07 号脚本；跑一次即缓存。未跑时 `Fig2b_TF_enrichment.R` 会明确报错并指向该脚本 |
| Fig 2b 家族色带颜色可能与原图有差 | 从原图采样得到的 `ref/fig2b_family_palette.rds` 未随仓库提供 | 脚本自动降级为内置 36 色：**家族分块顺序与位置不变，仅色带颜色不同**（会打印 warning） |
| **Fig 2b 的绘图代码在官方仓库中缺失** | `07.…R` 只做到 `RunChromVAR()` 与 `FindMotifs()`，第 89–91 行仅切换 `DefaultAssay`，**没有画 Fig 2b 的代码** | 家族分块与块间顺序依据原图图例复原；**块内排序规则原文未披露**，本复现按 TF 活性在细胞类型轴上的加权质心排序，故块内行序与原文不完全一致 |
| 图写不进 `output/`，报 `unable to translate ... to native encoding` | 仓库放在中文目录时，R 的 `png()` / `pdf()` / `cairo_pdf()` 都把文件名转成本地编码失败 | `_paths.R` 的 `fig_open()` 先试直接写，失败就画到 ASCII 临时文件、`fig_close()` 时再 `file.copy` 过去。**不要改用 `ragg::agg_png()`**：它能写中文路径，但本机退出时 segfault（exit 139） |
| Fig 2a 基因集与官方不同 | 官方 `VariableFeatures(RNA)` 只覆盖 1,959 个共有基因，样本量不足 | 改为 9,081 个「检出率 ≥5% ∪ 极差 top 6000」，理由见第 1.3 节 |
| 官方脚本的几处瑕疵 | `06.Correlation_analysis_of_RNA_and_ATAC.R` 把 `VariableFeatures()` 写成 `VarableFeatures()`（R 中无此函数，实际未生效）；用了 `pheatmap::pheatmap()` 却没 `library(pheatmap)`；到处是硬编码绝对路径 `/data2/lidongwei/...` | 本复现逐一修正，并全部改为相对路径 |
| R 在 Windows 上 `LC_CTYPE=C` 时中文报错会乱码 | `stop()` / `conditionMessage()` 的条件消息会被重编码（`cat()` 直接输出则正常） | `_paths.R` 开头尝试切 locale 修正中文路径显示；**运行时提示文字一律用英文**（源码里的中文字节会被 GBK 误读，连 `\n` 都会失效） |

---

## 7. 数据与代码来源声明

本节对应课程对"数据来源说明"与学术诚信的要求，逐项说明本仓库中**哪些是别人的、哪些是自己的**。

### 7.1 原始数据（非本项目产生）

| 数据 | 来源 | 获取方式 |
|---|---|---|
| 单细胞多组学测序数据（RNA + ATAC，8 器官） | NCBI GEO [**GSE232863**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232863) | 作者公开，可匿名下载 |
| 穗发育阶段数据（SP3–SP8） | GEO **GSE273875** | 同上 |
| 突变体 RNA-seq | GEO **GSE245410** | 同上 |

本仓库**未重新测序**，全部结果基于上述公开数据的再分析。

### 7.2 官方代码（作者提供，非本项目原创）

- Zenodo：[10.5281/zenodo.15501781](https://doi.org/10.5281/zenodo.15501781)（单细胞多组学分析代码）
- GitHub：`dongwei-2023/Single_cell_multiomics_in_rice`
- 本仓库 `data/official_code/` 下保留原始脚本副本，**仅用于对照**，不参与运行

### 7.3 本仓库的原创部分

| 内容 | 说明 |
|---|---|
| `R/figures/*.R` | 参考官方代码改写。**脚本头部逐处标注了对应的官方脚本与行号**；改动点（Fig 2a 的基因集选择、Fig 2b 的色标映射、Fig 1b/1c 的路径与内存适配）在第 1 节与第 6 节说明 |
| `ref/*.csv` | 从论文 PDF 矢量图层反演得到的色标与数值表，方法见 [`ref/README.md`](ref/README.md) 第 3 节，可复核 |
| `reports/*.docx` | 自行撰写 |

### 7.4 图版版权

原论文图版（Fig 1–5、Extended Data 各图）版权归 Nature / 原作者所有。
本仓库**不存储原图**，仅保存从图中反演出的**色标数值**（几十 KB 的 CSV）用于复现配色。
`figs/` 目录已在 `.gitignore` 中排除。

### 7.5 生成式 AI 辅助说明

脚本的语法检查、文档整理与语言润色过程中使用了 AI 辅助工具；
**分析方案的选择、参数的确定、结果的核验与结论的判断由作者完成**。
本仓库所有数值结论均可用其中的脚本重新生成并核对。

### 7.6 可复现性边界（如实说明）

| 类别 | 内容 |
|---|---|
| **可直接复现** | Fig 1b / 1c / 2a / 2b 的图，前提是已下载 5.1 GB 主数据对象并装好相应 R 包 |
| **无法完全复现** | Fig 1c 的**绝对准确率**（原文口径有三处歧义，见第 6 节）；Fig 2b 的**原始数值**（原文统计口径三处不一致，本复现用 chromVAR 代理） |
| **未复现** | 论文 Fig 3–5、Extended Data 各图、pySCENIC 的 251 个 regulon、CellOracle 扰动模拟 |

> **一句话**：这个仓库复现的是**四张图的生成过程与差异归因**，不是整篇论文的结论。
> 凡是无法完全对齐的地方，第 6 节逐条写明了原因，没有隐去。
