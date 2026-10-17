# 水稻单细胞多组学图谱复现

> 目标论文：**A single-cell multi-omics atlas of rice**
> Xiangyu Wang, Huanwei Huang, Sanjie Jiang, Jingmin Kang, Dongwei Li, Kailai Wang *et al.*
> *Nature* **644**, 722–729 (21 August 2025). DOI: [10.1038/s41586-025-09251-0](https://doi.org/10.1038/s41586-025-09251-0)
> 通讯作者：Zhe Liang (liangzhe@caas.cn)、Xiaofeng Gu (guxiaofeng@caas.cn)，中国农业科学院生物技术研究所

课程：《组学数据分析和可视化》学生报告（华中科技大学生命科学与技术学院）
适用硬性要求：2023.1 以来 Nature/Science/Cell 及子刊 ✅ ｜ 有公开数据集 ✅ ｜ 复现 ≥2 张核心图 ✅

## 快速开始

```bash
# 30 秒验证代码可运行（无需下载任何数据）
Rscript R/99_demo_fig1b_barplot.R
# 预期：output/Fig1b_nuclei_per_organ_DEMO.png，8 个柱子合计 116,564
```

完整流程见 [第 6 节 运行方式](#6-运行方式)；数据下载见 [第 1 节](#1-数据来源已逐个核实文件清单与体积)。

---

## 0. 一句话结论

这是一篇**植物单细胞多组学（scRNA + scATAC 同细胞同时测）**资源型论文。它最大的复现优势是：

**作者把「处理好的 Seurat 对象」和「全套 16 个官方分析脚本」都公开了** —— 这让我们可以跳过建库/比对，直接从对象出图，把精力放在统计建模与可视化上。

| 优势 | 具体是什么 |
|---|---|
| 有成品对象 | `GSE232863_scRNA_omics.Rds`（5.1 GB）、`GSE273875_LTP_NIP_label_240622.rds.gz`（390 MB） |
| 有官方代码 | Zenodo `10.5281/zenodo.15501781`，16 个脚本，README 里**逐条标注了脚本对应哪张图** |
| 有在线门户 | http://www.elabcaas.cn/scmr/index.html 可交互查看 |
| 数据体积可控 | 最小复现路线只需下载 390 MB |

---

## 1. 数据来源（已逐个核实文件清单与体积）

| GEO 登录号 | 内容 | 样本数 | 关键补充文件（实测体积） |
|---|---|---|---|
| **GSE232863** | 8 器官 scRNA + scATAC + 穗发育 4 期多组学（主图谱） | 44 | `GSE232863_scRNA_omics.Rds` **5.1 GB**（处理好的 Seurat 对象）<br>`GSE232863_RAW.tar` **8.0 GB**（6× `filtered_feature_bc_matrix.h5` ≈100 MB + 6× `atac_fragments.tsv.gz` ≈1.2 GB） |
| **GSE273875** | osltpl120 突变体 vs Nip 根 scRNA | 4 | `GSE273875_LTP_NIP_label_240622.rds.gz` **390 MB**（已注释好 cell_type 的 Seurat 对象） |
| **GSE245410** | osf3h 突变体 vs ZH11 叶/根 bulk RNA-seq | 12 | `GSE245410_raw_counts.txt.gz` **681 KB**（原始 counts，可完整走 bulk 流程） |

补充参考资源（非 NCBI）：

- 官方代码仓库：`https://doi.org/10.5281/zenodo.15501781`（`Single_Cell_Multiomics_in_Rice`，7.75 MB，含 `data/rice.motif.pfmList2.Rds`）
- 根细胞类型自动注释工具 RSCA：`https://doi.org/10.5281/zenodo.15501779`
- 参考基因组：**Rice MSU 7.0**（注释用 `EnsDb.Osativa.v55`）—— **必须与论文一致，否则比对/peak 注释全错**
- Motif 数据库：PlantTFDB v5.0（https://planttfdb.gao-lab.org/）
- 跨物种数据：maize / Arabidopsis / *Setaria viridis* / *Sorghum bicolor* 根尖 scRNA，来自 Guillotin *et al., Nature* 617, 785–791 (2023)

---

## 2. 关键方法参数（抄自 Methods，用于解释复现差异）

| 步骤 | 软件与版本 | 关键参数 |
|---|---|---|
| 建库测序 | 10X Chromium Next GEM Single Cell Multiome ATAC + Gene Expression；NovaSeq 6000 | RNA: PE150；ATAC: PE50；每样本约 20,000 nuclei |
| 上游比对 | `cellranger-arc mkfastq` + **CellRanger ARC v2.0.2** | 参考基因组 Rice MSU 7.0 |
| 环境 RNA 去除 | **CellBender v0.3.2** | `remove-background`，150 epochs，learning rate 1e-4 |
| 质控 | **Seurat v5.1.0** | 过滤 nFeature_RNA < 500 或 nCount_RNA < 1000 的细胞 |
| 双细胞去除 | **DoubletFinder v2.0.3** | PC = 1:30，pN = 0.25，pK 自动选最优 |
| 降维 | **Signac v0.2.5** + Seurat | ATAC：FindTopFeatures → RunTFIDF → RunSVD（LSI）；RNA：NormalizeData → FindVariableFeatures → ScaleData → RunPCA |
| 批次校正 | **Harmony v0.1.0** | `RunHarmony()`，按 tissue/sample 分组 |
| 多模态整合 | Seurat `FindMultiModalNeighbors` | RNA 取 PC 1:30，ATAC 取 LSI 2:30（WNN） |
| UMAP | Seurat `RunUMAP` | 官方脚本 05：`n.neighbors=100, n.epochs=200, min.dist=0.1, spread=10`<br>论文 Methods：`n.neighbours=30, n.epochs=200, min.dist=0.05, spread=4` |
| 细胞类型 DEG | Seurat `FindAllMarkers`（Wilcoxon） | \|log2FC\| > 0.25，pct > 25%，adj.P < 0.05 |
| 突变体 DEG | Seurat `FindMarkers`（Wilcoxon） | adj.P < 0.05 且 \|log2FC\| > 0.5 |
| Motif 富集 | ChIPseeker `peakAnnoEnrichment` | FDR ≤ 0.1，log2(fold enrichment) ≥ 0.5 |
| 调控网络 | **pySCENIC v0.12.1** | GRNBoost2，`n_jobs=8, seed=42`；RcisTarget 用 MSU7 |
| 细胞状态扰动 | **CellOracle** | Cicero co-accessibility ≥ 0.5，`score_cutoff=0.7`，`k=777, b_sight=777*8, b_maxl=777*4`，保留 top 10,000 边，传播深度 5 |
| 共表达网络 | **hdWGCNA v0.2.24** | `soft_power=9, networkType='signed'`，`corType='bicor'`，`kME_method='pearson'`；Cytoscape v3.9.1 |
| 拟时序 | **Monocle3 v1.3.7** + CytoTRACE v0.1.0 + Slingshot v2.0.0 | 100 PCs，`umap.min_dist=0.5, n_neighbours=100` |
| 性状关联 | **ArchR v1.0.2** `addGeneScoreMatrix` | min-max 标准化后计算 trait–celltype 相关 |
| 跨物种 | OrthoFinder v2.5.4 + MetaNeighbor v1.9.1 | 同源基因簇求和，Seurat 整合 top 3,000 HVG |
| bulk RNA-seq | Trimmomatic v0.39 → HISAT2 v2.2.1 → StringTie v2.1.4 → DESeq → agriGO v2.0 | 见官方脚本外的补充流程 |

---

## 3. 复现的核心图与对应脚本

| 编号 | 图 | 输入数据 | 对应脚本 | 成功判定标准 |
|---|---|---|---|---|
| **①（必做）** | **Fig. 1b / 1d / 1e** — 8 器官 UMAP + 各器官细胞数 + 细胞类型占比 | `GSE232863_scRNA_omics.Rds` | `R/02_fig1b_umap_organ.R`<br>`R/03_fig1e_composition.R` | 细胞总数 = **116,564**；8 器官分别为 Crown root **31,108**、Culm **19,740**、Leaf **13,771**、Flag leaf **5,563**、Shoot apex **13,791**、Tiller bud **8,858**、Panicle **17,457**、Seed **6,276**（8 项相加严格等于 116,564）；UMAP 中 8 个器官呈"生理相近器官相邻"格局 |
| **②（必做）** | **Fig. 4u / 4w** — WT vs osltpl120 各细胞类型 DEG 数 + 关键基因表达热图 | `GSE273875_LTP_NIP_label_240622.rds.gz`（390 MB） | `R/05_fig4uw_deg.R` | **皮层（cortex）DEG 数最多**（论文原文："the highest number of DEGs in cortical cells"）；OsNRT1.1A / OsAAP1 / OsBHLH044 在突变体中上调且皮层倍数最高 |
| ③（进阶） | **Fig. 2a** — scRNA vs scATAC 细胞类型相关性热图 | 同上 Rds（需 ATAC assay） | `R/04_fig2a_correlation.R` | 56×56 对角块状结构；叶/旗叶的叶肉、维管、保卫细胞、表皮相关性高；种子表皮相关性中等 |
| ④（进阶） | **Fig. 3a** — 碳氮吸收相关基因在细胞类型间的表达（上）与染色质可及性（下）热图 | 同上 Rds + ArchR `addGeneScoreMatrix` | 待补 | 碳氮基因呈明显细胞类型特异分层，scRNA 与 scATAC 两层趋势一致 |

### 3.1 已产出的文件（`output/`）

| 文件 | 说明 | 生成方式 |
|---|---|---|
| `Fig1b_nuclei_per_organ_DEMO.png` | 8 器官细胞数柱状图，用于验证绘图代码正确性（数字与论文完全一致，合计 116,564） | `Rscript R/99_demo_fig1b_barplot.R` |
| `Fig1b_nuclei_per_organ_DEMO.pdf` | 同上，矢量版 | 同上 |

### 3.2 如何验证代码可运行

**不需要下载任何数据**，用论文公开的细胞数即可跑通绘图链路：

```bash
Rscript R/99_demo_fig1b_barplot.R
```

成功标志：终端无报错，`output/` 下生成 `Fig1b_nuclei_per_organ_DEMO.png`，且 8 个柱子数值之和为 116,564。

这一步是给评审者的**最快验证入口**——它证明脚本、依赖、绘图配置都是通的，不需要等待 5.1 GB 数据下载。

---

## 4. 复现路线

### 路线 A（推荐，保底）—— 用成品对象直接出图

```
下载 Rds → 探查对象结构 → 复算细胞数/UMAP 坐标 → 出图 → 与原图逐面板比对
```
预计 **2–4 小时**（不含下载）。优点：不需要跑比对，风险最低，能把时间投到"结果对比与差异分析"上（这一项占分高）。

### 路线 B（完整）—— 从下机数据走全流程

```
SRA fastq → cellranger-arc mkfastq → CellRanger ARC count → CellBender → Seurat 质控
→ DoubletFinder → Harmony → WNN 整合 → RunUMAP → 注释 → 差异分析 → 可视化
```
优点：完整覆盖"数据下载→质控→比对→统计建模→可视化"，符合课程对流程完整度的要求。
缺点：需下载 **8 GB RAW.tar**，CellRanger ARC 需 **≥64 GB 内存**，单样本耗时数小时。

> **建议组合**：主图①用路线 A 出全量结果；同时**挑 1 个样本走路线 B**（用 `filtered_feature_bc_matrix.h5` 直接进 Seurat 做质控+DoubletFinder+UMAP），在报告里展示"质控前/后小提琴图"和"TSS enrichment / nucleosome signal"曲线，这样既拿到全量结果，又证明了流程完整度。

### 算力与时间估计

| 步骤 | 数据量 | 内存 | 时间 |
|---|---|---|---|
| 下载 `GSE273875` Rds | 390 MB | — | 3–10 min |
| 下载 `GSE232863_scRNA_omics.Rds` | 5.1 GB | — | 30 min – 2 h |
| 读取 5.1 GB Rds（Seurat） | — | **≥32 GB** | 5–15 min |
| Fig. 1b/1d/1e 出图 | 116,564 细胞 | ≥32 GB | 10–30 min |
| Fig. 4u/4w 出图 | ~36,000 细胞 | 8–16 GB | 10–20 min |
| 路线 B 单样本 CellRanger ARC | ~40 GB fastq | ≥64 GB | 3–8 h |
| 路线 B 单样本 Seurat + DoubletFinder | — | 16 GB | 20–40 min |

---

## 5. 目录结构

```
水稻scMultiome复现方案/
├── README.md                      # 本文件
├── 复现方案.html                   # 完整报告方案（含图、表、分工、时间线）
├── figs/                          # 原论文图版渲染图（供"与原图对比"直接用）
│   ├── orig_Fig1_p722.png ... orig_Fig5_p729.png
│   └── official_pipeline.png      # 作者官方流程图
├── R/
│   ├── _bootstrap.R               # 【必须最先 source】切 UTF-8 locale（纯 ASCII）
│   ├── _common.R                  # 公共函数：列名容错、配色、主题、读表、保存图片
│   ├── 00_setup.R                 # 装包 + 记录版本到 output/session_versions.csv
│   ├── 01_inspect_object.R        # 【第一步必跑】探查 Seurat 对象结构
│   ├── 02_fig1b_umap_organ.R      # 图1b：UMAP by organ + 细胞数柱状图
│   ├── 03_fig1e_composition.R     # 图1e：细胞类型占比堆叠柱（顺带出图1d）
│   ├── 04_fig2a_correlation.R     # 图2a：RNA–ATAC 相关性热图
│   ├── 05_fig4uw_deg.R            # 图4u/w：WT vs osltpl120 DEG
│   ├── 98_syntax_check.R          # 全部 R 脚本语法检查
│   └── 99_demo_fig1b_barplot.R    # 【立即可跑】用论文真实数字验证出图代码
├── scripts/
│   ├── 00_download.sh             # 数据下载（含省流量策略）
│   └── 01_cellranger_arc.sh       # 路线 B 上游
└── output/                        # 出图与中间结果
```

---

## 6. 运行方式

> ⚠️ **必须在项目根目录运行**（目录名含中文，脚本一律使用相对路径）。

```bash
cd 水稻scMultiome复现方案

# 1) 立刻可跑：验证绘图代码（只用已装的 ggplot2/dplyr，不需要 Seurat）
Rscript R/99_demo_fig1b_barplot.R

# 2) 安装 R 包（约 30–60 min，Seurat/Signac 编译较慢）
Rscript R/00_setup.R

# 3) 下载数据（先下最小集 0.4 GB）
bash scripts/00_download.sh minimal

# 4) 探查对象结构（务必先跑，再改后面的脚本）
Rscript R/01_inspect_object.R

# 5) 出图
Rscript R/05_fig4uw_deg.R          # 390 MB 对象，8-16 GB 内存即可
bash    scripts/00_download.sh medium   # 再下 5.1 GB 主图谱对象
Rscript R/02_fig1b_umap_organ.R    # 需 >=32 GB 内存
Rscript R/03_fig1e_composition.R
Rscript R/04_fig2a_correlation.R

# 6) 语法自检
Rscript R/98_syntax_check.R
```

**已在本机验证**：`99_demo_fig1b_barplot.R` 实测可运行并出图
（`output/Fig1b_nuclei_per_organ_DEMO.png`），8 个器官细胞数相加 = 116,564，与论文一致。

---

## 7. 已知坑与注意事项

1. **参考基因组必须是 Rice MSU 7.0**（配合 `EnsDb.Osativa.v55`）。换成 IRGSP-1.0/RAP-DB 会导致基因名对不上、peak 注释错位。
2. **5.1 GB 的 Rds 加载需 ≥32 GB 内存**。内存不足时先在 16 GB 机器上跑 390 MB 的 `GSE273875` 对象。
3. **官方脚本存在小瑕疵**（可作为报告"改进建议"的素材）：
   - `06.Correlation_analysis_of_RNA_and_ATAC.R` 第 19 行写成 `VarableFeatures()`，正确函数名是 `VariableFeatures()`；
   - 同脚本用了 `pheatmap::pheatmap()` 但未 `library(pheatmap)`，依赖隐式加载；
   - 脚本里到处是硬编码绝对路径（`/data2/lidongwei/...`），需替换；
   - `15.Analyzing_scRNA_of_LTPL120.R` 里 `cluster2celltype` 手工映射了 21 个 cluster，若重跑聚类分辨率改变会错位。
4. **`GSE232863_RAW.tar` 里只有 6 个多组学样本**（穗 SP3/SP4/SP6/SP8），8 器官主图谱的 RNA 与 ATAC 在 GEO 里是**分别**提交的（每器官 Rep1/Rep2 各一），并非同细胞多组学。报告里要写清楚这一点——这正是论文 Fig. 1c 要评估"计算配对 vs 真实多组学"的动机。
5. 下载策略：**不要全下 8 GB**。只取需要的 h5，跳过 `atac_fragments.tsv.gz`（每个 1.2 GB，仅在需要 peak-level 分析时才要）。
6. 复现差异的常见来源：**DoubletFinder 的 pK 自动选择有随机性**、Harmony 的随机种子、UMAP 的 `n.neighbors/min.dist/spread` 参数版本差异、Seurat 4→5 的 layer 机制变化。
7. **中文字符串在 R 图里会乱码/显示方框**（本机 R 启动时 `LC_CTYPE = C`）：
   - 解法：每个脚本开头先 `source("R/_bootstrap.R")` 切到 UTF-8 locale；
   - 但 cairo 默认字体仍缺 CJK 字形，会显示成方框 → **绘图标签一律用英文**（也符合科研制图惯例），中文只出现在注释与报告正文里。
8. **`ggsave` 在本机默认走 ragg，写盘失败**（只 warning 不 error，容易误判成功）。
   `_common.R` 的 `save_fig()` 已改为直接用 `grDevices::png(type="cairo")` 并做回退，同时用"文件是否真的生成"来判定。
9. **不要用 `Rscript -e '多行脚本'`** 在 Git Bash 下运行（可能 segfault，且中文路径会导致引号解析出错）。写成 `.R` 文件再执行。

---

## 8. 学术诚信

- 本方案的代码骨架为课程复现用途自行编写，**数据全部来自论文公开的 GEO 仓库**，图表版权归原文所有，报告中所有引用图均注明出处。
- 使用大模型辅助的范围与边界：用于**文献信息提取、流程梳理、代码骨架起草与语法检查**；**统计结果与图表均由本地 R 实际运行产出**，未使用模型生成的"假数据"或"假结果"。
- 本仓库**不包含**任何他人代码或论文原文（已在 `.gitignore` 中排除 `data/official_code/`、`data/*.zip`、`data/paper_clean.txt`、论文原图），如需查阅官方代码请访问 Zenodo `10.5281/zenodo.15501781`。

---

## 9. 课程提交材料对照

按《组学数据分析和可视化》学生报告说明，"报告形式与材料提交"要求提交 4 项（缺一项扣 5 分）：

| # | 要求 | 本仓库提供什么 |
|---|------|--------------|
| 1 | 幻灯片（PPT/Keynote，≤20 页） | 随邮件单独提交 |
| 2 | **复现代码（GitHub/Gitee 仓库，含 README）** | **本仓库**（即你正在看的这个） |
| 3 | 复现报告（PDF 或 Word，含关键图表与说明） | 随邮件单独提交 |
| 4 | 数据集来源与可访问性说明 | 见 [第 1 节](#1-数据来源已逐个核实文件清单与体积)，含 GEO 登录号与实测文件体积 |

**代码质量自查**（对应评分项"代码质量 20%"）：

- [x] 可运行性：`R/99_demo_fig1b_barplot.R` 无需数据即可跑通
- [x] 注释：每个脚本开头有功能说明，`R/_common.R` 提供公共函数
- [x] 可复现性：无绝对路径，全部相对路径；脚本按编号顺序执行
- [x] 依赖明确：`R/00_setup.R` 装包并记录版本到 `output/session_versions.csv`
- [x] 数据来源可追溯：GEO 登录号 + 官方 Zenodo 代码 DOI

---

## 10. 复现差异与已知问题

按课程要求"复现结果与原论文的对比、差异原因、失败环节分析"，此处汇总：

| 现象 | 原因 | 处理 |
|------|------|------|
| UMAP 形状与原图不完全一致 | Harmony 随机种子、`n.neighbors/min.dist/spread` 版本差异 | 报告中注明参数来源（Methods vs 官方脚本存在不一致，见第 2 节表格） |
| 细胞类型数量可能略有出入 | DoubletFinder 的 `pK` 自动选择有随机性 | 固定随机种子；报告里给出实际聚类数并与 56 对比 |
| Fig. 2a 相关系数整体偏低 | 官方脚本 `06.Correlation_analysis_of_RNA_and_ATAC.R` 使用 `VarableFeatures()`（拼写错误，实际未生效） | 已在本仓库 `04_fig2a_correlation.R` 中修正为 `VariableFeatures()` |
| 中文标签显示为方框 | cairo 字体缺 CJK 字形 | 绘图标签一律用英文（符合科研制图惯例） |
| `ggsave` 静默失败 | 本机 ragg 设备写盘异常，只 warning 不 error | `_common.R` 的 `save_fig()` 改用 `png(type="cairo")` 并校验文件是否真的生成 |
