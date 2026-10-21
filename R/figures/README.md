# `R/figures/` —— 四张核心图的最终复现脚本

这个子目录放的是**报告里那四张图的实际生成脚本**，也就是"最终交付版"。

和 `R/` 根目录下的 `02_` / `03_` / `04_` / `05_` 的区别：

| | `R/figures/`（本目录） | `R/02_…` `R/04_…`（根目录） |
|---|---|---|
| 定位 | **最终交付版**，报告里的图由它生成 | 入门版 / 教学版，用来理解流程 |
| 复现精度 | 与原图逐格比对（Fig 2a: Spearman 0.967） | 只求结构一致，不追数值 |
| 依赖 | 需要 `ref/` 下的数字化参考表 | 只需要主数据对象 |
| 代码量 | 200–400 行，参数逐项对齐原图 | 100–150 行，重在可读 |

> 换句话说：**要交的图跑这个目录；想快速看懂"这一步在干什么"看根目录。**

---

## 1. 五个脚本

| 脚本 | 对应图 | 内容 | 运行时间 | 依赖 |
|---|---|---|---|---|
| `00_selftest.R` | — | **先跑这个**：检查路径 / 参考表 / 原始数据 / R 包，打印四张图的就绪状态 | 10 s | 无 |
| `Fig1b_reproduce.R` | **Fig. 1b** | 全细胞 WNN-UMAP（按 8 个器官着色）+ 各器官细胞数条形图 | ~6 min | Seurat 5 / Signac / ggplot2 / patchwork / dplyr，内存 ≥16 GB |
| `Fig1c_reproduce.R` | **Fig. 1c** | RNA→ATAC 标签迁移的逐器官准确率（8 器官 × 2 根横条） | 30 min – 2 h | Seurat 5 / Signac / ggpattern / dplyr，内存 ≥16 GB |
| `Fig2a_correlation.R` | **Fig. 2a** | scRNA vs scATAC 细胞类型平均谱的 Spearman 相关热图，56×56 | ~2.5 min | Seurat 5 / Signac / ComplexHeatmap，内存 ≥32 GB |
| `Fig2b_TF_enrichment.R` | **Fig. 2b** | 254 TF × 35 细胞类型的 motif deviation score 热图，按 36 个 TF 家族分块 | 30 min – 3 h | 上面这些 + chromVAR / motifmatchr / BSgenome.Osativa.MSU.MSU7 |
| `00b_make_primary_markers.R` | （Fig 2b 前置） | 为每个细胞类型在 ATAC assay 上跑 `FindMarkers`，产出 `primary_markers/*.rds` | 30 min – 1 h | Seurat / Signac |

**建议的验证顺序**：先跑 `00_selftest.R` 看就绪状态 → 再跑 `Fig1b`（只读元数据与坐标，最快出图）
→ 然后 `Fig1c`（要跑锚点，最久）→ 再 `Fig2a`（2.5 分钟）→ 最后 `Fig2b`（需要 chromVAR 全家桶）。

---

## 2. 怎么跑

```bash
# 在仓库根目录
cd 水稻scMultiome复现方案

Rscript R/figures/00_selftest.R              # 10 秒，先看能不能跑
Rscript R/figures/Fig1b_reproduce.R          # 约 6 分钟
Rscript R/figures/Fig1c_reproduce.R          # 30 分钟起
Rscript R/figures/Fig2a_correlation.R        # 2.5 分钟
Rscript R/figures/Fig2b_TF_enrichment.R      # 30 分钟起
```

也可以在 `R/figures/` 目录内直接 `Rscript Fig2a_correlation.R`，脚本会自动找到 `_paths.R`。

**数据不在默认位置时**（比如放在移动硬盘），不用改脚本，设一个环境变量即可：

```r
Sys.setenv(RICE_PROJ = "D:/水稻scMultiome复现方案")
```

或单独指定主数据对象：

```r
Sys.setenv(RICE_RDS_SCRNA = "E:/data/GSE232863_scRNA_omics.Rds")
```

---

## 3. 路径约定（重要）

四张图的脚本都**不写死绝对路径**，但定位方式分两种：

- `Fig2a` / `Fig2b`：`source("_paths.R")`，由它向上查找仓库根目录；
- `Fig1b` / `Fig1c`：默认读仓库相对路径 `data/…` 与 `output/…`，可用环境变量覆盖
  （`RICE_RDS_SCRNA` / `RICE_OUT_FIG1B` / `RICE_OUT_FIG1C`）。

约定如下：

| 变量 | 指向 | 是否进仓库 |
|---|---|---|
| `DIR_DATA` | `data/` 原始数据（Rds / h5 / tar） | ❌ `.gitignore` 排除（体积太大） |
| `DIR_REF` | `ref/` 从原图反演的小参考表 | ✅ 提交（几十 KB） |
| `DIR_OUT` | `output/` 产物 | ✅ 只提交图片，`.rds` 中间结果排除 |
| `RDS_SCRNA` | `data/GSE232863_scRNA_omics.Rds` | ❌ 需自行从 GEO 下载 |

产物统一写到 `output/<图号>/`，例如 `output/fig2a/Fig2a_reproduced.png`。
`out_dir("fig2a")` 会自动建目录。

> **为什么不用 `here::here()`**：多一个依赖，而且 `here` 依赖 `.Rproj` 或 `.git` 的位置，
> 在"整个文件夹打包发给老师"的场景下会失效。这里用「从工作目录逐级向上找 `README.md` + `R/`」
> 的土办法，零依赖，且对"在仓库根目录跑"和"在脚本目录跑"都成立。

---

## 4. 各脚本的输入清单

### Fig 1b

| 文件 | 必需 | 说明 |
|---|---|---|
| `data/GSE232863_scRNA_omics.Rds` | ✅ | GEO 下载，5.1 GB。**只读取元数据与降维坐标**：`obj@reductions` 里的 WNN-UMAP，以及 `meta.data` 的 `tissue` 列 |
| `output/fig1b/` | — | 产物目录，脚本自动创建 |

> Fig 1b **不读任何表达值**，属于可视化层面复现：坐标与器官注释均为作者发布结果，
> 复现目标是版式、配色与细胞数比例。脚本可选地用 `WNN_SOURCE = "recompute"` 重算 UMAP 做一致性检验。

### Fig 1c

| 文件 | 必需 | 说明 |
|---|---|---|
| `data/GSE232863_scRNA_omics.Rds` | ✅ | 同上。需要 `RNA`（参考集）与 `atacRNA`（查询集，作者已算好的基因活性矩阵）两个 assay；`ATAC` 峰矩阵在脚本内被丢弃以省内存 |
| `output/fig1c/results/` | — | 产物目录，脚本自动创建 |

> `atacRNA` 不需要重跑 `GeneActivity()` —— 作者已把它存在对象里。
> 脚本按器官分批做锚点，避免 115,395 × 115,395 的全局锚点把内存打爆。

### Fig 2a

| 文件 | 必需 | 说明 |
|---|---|---|
| `data/GSE232863_scRNA_omics.Rds` | ✅ | GEO 下载，5.1 GB |
| `ref/fig2a_palette.csv` | ✅ | 已提交 |
| `ref/fig2a_dots.csv` | ✅ | 已提交 |
| `ref/fig2a_matrix_values_exact.csv` | ⬜ | 可选，仅用于打印逐格比对指标 |

### Fig 2b

| 文件 | 必需 | 说明 |
|---|---|---|
| `data/GSE232863_scRNA_omics.Rds` | ✅ | 同上 |
| `data/official_code/dongwei-2023-Single_Cell_Multiomics_in_Rice-6020ed1/data/rice.motif.pfmList2.Rds` | ✅ | 官方代码包（Zenodo `10.5281/zenodo.15501781`）里的 PlantTFDB motif PFM |
| `output/fig2b/primary_markers/*.rds` | ✅ | **前置步骤产出**：对每个细胞类型在 ATAC assay 上跑 `FindMarkers`，存成含 `avg_log2FC` 列、行名为 peak 坐标的表 |
| `ref/fig2b_family_palette.rds` | ⬜ | 可选；缺失时用内置 36 色备用配色 |
| `data/Fig2b_crop.png` | ⬜ | 可选，仅 Step 4 并排对比用 |

> **`primary_markers/` 由 `00b_make_primary_markers.R` 生成**（30 min – 1 h）：
>
> ```bash
> Rscript R/figures/00b_make_primary_markers.R
> ```
>
> 该脚本逐行对齐官方 `07.Analyzing_the_TF_and_motif_by_Signac.R` 的 `FindMarkers` 参数
> （`test.use = "LR"`、`only.pos = TRUE`、`logfc.threshold = 0.05`、
> `max.cells.per.ident = 300`、`latent.vars = "nCount_ATAC"`），
> 每类存一张表到 `output/fig2b/primary_markers/<细胞类型>.rds`。
> 表里的行名就是 peak 坐标，`Fig2b_TF_enrichment.R` 用 `rownames()` 直接 match。
> 跑过一次即缓存，重画 Fig 2b 不需要再跑（`FORCE=1` 可强制重算）。

---

## 5. 产物一览

```
output/fig1b/
├── Fig1b_reproduced.png / .pdf        左：WNN-UMAP（按 8 器官着色）；右：各器官细胞数
├── Fig1b_cell_counts.csv              8 器官细胞数 + 与原文 116,564 的逐项差值
├── Fig1b_consistency_vs_authors.csv   重算 UMAP 与作者存档坐标的一致性检验
├── Fig1b_umap_coords.csv              画图用的二维坐标（便于第三方复核）
└── Fig1b_weights_check.csv            WNN 权重矩阵的来源核对

output/fig1c/results/
├── Fig1c_reproduced.png / .pdf        8 器官 × 2 根横条（实心 = 配对准确率，斜纹 = 多组学 100%）
├── Fig1c_reproduced.csv               逐器官 accuracy
└── Fig1c_confusion_<organ>.csv        逐器官混淆矩阵（失败环节分析用）

output/fig2a/
├── Fig2a_reproduced.png          2700×2580 @600dpi（报告用）
├── Fig2a_reproduced.pdf          4.5×4.3 in 矢量（投稿/放大用）
├── Fig2a_best_56x56_raw.csv      未截断的相关系数矩阵
├── Fig2a_best_56x56_clipped.csv  截断到 [0, 0.6] 后的矩阵（= 画图用的那份）
├── Fig2a_axis_mapping.csv        原图编号 ↔ GEO final_id2 ↔ 器官 的对照
├── Fig2a_best_geneset_summary.csv 基因集各步的基因数
└── Fig2a_C_per_tissue.csv        逐组织的对角/组织内/跨组织偏差

output/fig2b/
├── Fig2b_deviation_reproduced.png / .pdf
├── Fig2b_deviation_vs_original.png   与原图并排
├── cv_top2000_mat.rds                中间结果（不进仓库）
└── Fig2b_deviation_reproduced_data.rds
```

---

## 6. 复现质量（可写进报告的对比数据）

Fig 2a，基准 = 从原图反演出的 56×56 矩阵：

| 指标 | 本复现（C 版基因集） | 官方脚本直译 | 说明 |
|---|---|---|---|
| Spearman | **0.9668** | 0.7652 | 逐格相关系数 |
| Pearson | **0.9706** | 0.7791 | 同上 |
| MAE | **0.0211** | 0.0540 | 平均绝对误差 |
| 对角均值 | 0.408（原图 0.409） | — | 同一细胞类型 RNA vs ATAC 自相关 |
| 非对角均值 | 0.221（原图 0.225） | — | 不同细胞类型之间 |
| 标准差 | 0.1024（原图 0.1036） | — | 整体离散度 |

提升的来源是**基因集选择**：官方 `VariableFeatures(RNA)` 只有 1,959 个基因在两模态中都有，
样本量不足导致相关系数被系统性低估（低相关区抬高、高相关区压平）。
改用「两模态检出率均 ≥5%」∪「跨细胞类型极差 top 6000」= 9,081 个基因后，三个边际统计与原图几乎重合。
