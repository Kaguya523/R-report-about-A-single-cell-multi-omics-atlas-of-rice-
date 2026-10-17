# `R/figures/` —— 三张核心图的最终复现脚本

这个子目录放的是**报告里那三张图的实际生成脚本**，也就是"最终交付版"。

和 `R/` 根目录下的 `02_` / `03_` / `04_` / `05_` 的区别：

| | `R/figures/`（本目录） | `R/02_…` `R/04_…`（根目录） |
|---|---|---|
| 定位 | **最终交付版**，报告里的图由它生成 | 入门版 / 教学版，用来理解流程 |
| 复现精度 | 与原图逐格比对（Fig 2a: Spearman 0.967） | 只求结构一致，不追数值 |
| 依赖 | 需要 `ref/` 下的数字化参考表 | 只需要主数据对象 |
| 代码量 | 200–400 行，参数逐项对齐原图 | 100–150 行，重在可读 |

> 换句话说：**要交的图跑这个目录；想快速看懂"这一步在干什么"看根目录。**

---

## 1. 三个脚本

| 脚本 | 对应图 | 内容 | 运行时间 | 依赖 |
|---|---|---|---|---|
| `00_selftest.R` | — | **先跑这个**：检查路径 / 参考表 / 原始数据 / R 包，打印三张图的就绪状态 | 10 s | 无 |
| `Fig2a_correlation.R` | **Fig. 2a** | scRNA vs scATAC 细胞类型平均谱的 Spearman 相关热图，56×56 | ~2.5 min | Seurat 5 / Signac / ComplexHeatmap，内存 ≥32 GB |
| `Fig2b_TF_enrichment.R` | **Fig. 2b** | 254 TF × 35 细胞类型的 motif deviation score 热图，按 36 个 TF 家族分块 | 30 min – 3 h | 上面这些 + chromVAR / motifmatchr / BSgenome.Osativa.MSU.MSU7 |
| `Fig3a_gene_expression.R` | **Fig. 3a** | 56 细胞类型 × 65 基因的碳氮代谢基因 z-score 热图（上 scRNA / 下 scATAC） | < 10 s | **只要 base R + grid**，不需要任何 Bioconductor 包 |

**建议的验证顺序**：先跑 `00_selftest.R` 看就绪状态 → 再跑 `Fig3a`（最快、装包最少，10 秒出图）
→ 然后 `Fig2a`（2.5 分钟）→ 最后 `Fig2b`（最慢，且需要 chromVAR 全家桶）。

---

## 2. 怎么跑

```bash
# 在仓库根目录
cd 水稻scMultiome复现方案

Rscript R/figures/00_selftest.R              # 10 秒，先看能不能跑
Rscript R/figures/Fig3a_gene_expression.R    # 10 秒
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

三个脚本都 `source("_paths.R")`，**脚本正文里不出现任何绝对路径**。约定如下：

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

> ⚠️ **`primary_markers/` 的生成脚本目前没有随仓库提供**（它属于上游分析步骤，不在本次三张图的范围里）。
> 这一点在报告的"复现差异与已知问题"里已如实说明。
> 补法：对 `obj` 的每个 `tissue_cluster`，`FindMarkers(obj, ident.1 = ct, assay = "ATAC", only.pos = TRUE)`
> 后 `saveRDS(as.data.frame(mk), file.path(DIR_MK, paste0(gsub("[^A-Za-z0-9]", "_", ct), ".rds")))`。
> 文件名只需能被 `norm_key()` 规范化后与细胞类型名对上即可。

### Fig 3a

只需 `ref/fig3a/` 下的 9 个 CSV，清单见 [`../../ref/README.md`](../../ref/README.md)。
**不需要主数据对象** —— 这张图是从数字化矩阵重绘的。

---

## 5. 产物一览

```
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

output/fig3a/
├── Fig3a_scRNA_panel.png
├── Fig3a_scATAC_panel.png
└── Fig3a_reproduction_exact.png      上下两块 + 完整图例
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
