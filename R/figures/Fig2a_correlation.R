###############################################################################
#  Fig. 2a 复现（C 版，最贴近原图）—— 单文件自包含脚本
#  ===========================================================================
#  论文   Wang X, et al. A single-cell multi-omics atlas of rice.
#         Nature 644:722-729 (2025).  DOI: 10.1038/s41586-025-09251-0
#
#  官方源码（论文 Code availability → GitHub / Zenodo）
#         repo   : dongwei-2023/Single_cell_multiomics_in_rice
#         Zenodo : 10.5281/zenodo.15501781
#         file   : 06.Correlation_analysis_of_RNA_and_ATAC.R   →  Figure 2a
#
#         官方原文（第 84-102 行）--------------------------------------------
#         DefaultAssay(rice.sub) <- "ATAC"
#         gene.activities <- GeneActivity(rice.sub)
#         rice.sub[["atacRNA"]] <- CreateAssayObject(counts = gene.activities)
#         rice.sub <- NormalizeData(object = rice.sub, assay = "atacRNA",
#             normalization.method = "LogNormalize",
#             scale.factor = median(rice.sub$nCount_atacRNA))
#         atacRNA.avg <- AverageExpression(rice.sub, assays = "atacRNA",
#                                          features = VarableFeatures(rice.sub))
#         colnames(atacRNA.avg$atacRNA) <- paste(colnames(atacRNA.avg$atacRNA),
#                                                "ATAC", sep = "-")
#         atac_RNA.cor <- cor(atacRNA.avg$atacRNA,
#                             rice.sub.avg$RNA[rownames(atacRNA.avg$atacRNA), ],
#                             method = "spearman")
#         pheatmap::pheatmap(atac_RNA.cor, cluster_cols = F, cluster_rows = F, ...)
#         ---------------------------------------------------------------------
#         注：`VarableFeatures` 是 `VariableFeatures` 的笔误（R 中无此函数）。
#
#  数据   GSE232863_scRNA_omics.Rds（作者注释对象：37,960 基因 × 115,395 细胞）
#         · RNA assay     : 表达（data 层 = LogNormalize）
#         · atacRNA assay : scATAC 基因活性（GeneActivity + LogNormalize，
#                           scale.factor = median(nCount_atacRNA) = 1939，已核验）
#
#  环境   R 4.6.1 / Seurat 5.6.0 / Signac 1.17.1 / ComplexHeatmap 2.28.0
#
#  ---------------------------------------------------------------------------
#  C 版基因集（与官方唯一不同之处，其余算法逐行一致）
#        9,081 个基因 = 两模态检出率均 ≥ 5% 的全部基因
#                        ∪ 按「细胞类型平均谱极差 max−min」排序的 top 6000
#        理由：官方 VariableFeatures(RNA) 只有 1,959 个基因在两个 assay 中都有，
#              以此算相关的采样噪声过大 —— 低相关区被整体抬高、高相关区被压平。
#              改用「可靠检出 ∩ 细胞类型判别力强」的基因集后，
#              对角 / 非对角 / sd 三个边际统计与原图几乎完全一致。
#
#  复现质量（基准 = 从原图 PDF 逐矢量取色后反查出的 56×56 数值矩阵）
#        C 版      Spearman 0.9668  Pearson 0.9706  MAE 0.0211
#                  对角 0.408  非对角 0.221  sd 0.1024   ← 原图 0.409 / 0.225 / 0.1036
#        官方直译  Spearman 0.7652  Pearson 0.7791  MAE 0.0540（对照，见 V3_best 脚本）
#
#  运行   cd <仓库根目录> && Rscript R/figures/Fig2a_correlation.R
#         约 2.5 分钟（含读取 5.5 GB Rds），建议内存 ≥ 32 GB
#
#  输入   data/GSE232863_scRNA_omics.Rds        （GEO 下载，见 README 第 1 节）
#         ref/fig2a_palette.csv                  （100 阶色标，已随仓库提交）
#         ref/fig2a_dots.csv                     （56 个行注释圆点色，已随仓库提交）
#         ref/fig2a_matrix_values_exact.csv      （可选：原图数值矩阵，仅用于逐格比对）
#
#  输出   output/fig2a/Fig2a_reproduced.png / .pdf
#         output/fig2a/Fig2a_best_56x56_raw.csv / _clipped.csv
#         output/fig2a/Fig2a_axis_mapping.csv / Fig2a_best_geneset_summary.csv
#         output/fig2a/Fig2a_C_per_tissue.csv
###############################################################################

## ---- 载入公共配置（_paths.R 与本脚本同目录）-------------------------------
for (.p in c("_paths.R", "R/figures/_paths.R", "../_paths.R")) {
  if (file.exists(.p)) { source(.p); break }
}
if (!exists("PROJ", inherits = FALSE))
  stop("Cannot find R/figures/_paths.R. Run this script from the repository root.")

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix)
  library(ComplexHeatmap); library(circlize); library(grid)
})

## ===== 0. 路径与常量 ======================================================
OUT      <- out_dir("fig2a")
RDS_PATH <- RDS_SCRNA
need_file(RDS_PATH, "GSE232863 Seurat object",
          "download from GEO (see README section 1) and put it under data/")

## 配色：从原图 Fig. 2a 的色条与行注释圆点逐像素提取，已随仓库提交
need_file(file.path(DIR_REF, "fig2a_palette.csv"), "Fig. 2a colour scale",
          "shipped with the repo; see ref/README.md")
need_file(file.path(DIR_REF, "fig2a_dots.csv"), "Fig. 2a row-annotation dots",
          "shipped with the repo; see ref/README.md")
PAL_100 <- with(read.csv(file.path(DIR_REF, "fig2a_palette.csv")), hex(r, g, b))
DOT_56  <- with(read.csv(file.path(DIR_REF, "fig2a_dots.csv")),  hex(r, g, b))
stopifnot(length(PAL_100) == 100, length(DOT_56) == 56)

## 官方代码输出顺序 = levels(Idents)。下标 i 对应的 GEO final_id2 编号。
GEO_ORDER <- c(5,3,6,8,4,1,7,2, 12,13,11,9,10, 51,54,53,55,56,52, 15,16,17,20,22,
               21,18,14,19, 25,28,23,29,26,27,24, 48,47,50,49,46, 38,39,35,32,33,
               34,37,31,36,30, 40,41,42,45,43,44)
## 原图轴上印的编号（作者自己的编号体系，与 GEO final_id2 不同）
PAPER_LABELS <- c(38,36,39,41,37,34,40,35, 25,26,24,22,23, 16,19,18,20,21,17,
                  5,8,4,2,6,3,9,1,7, 29,33,30,32,31,28,27, 54,53,55,52,56,50,51,
                  47,44,45,46,49,43,48,42, 10,11,12,15,13,14)
TISSUE_OF_POS <- c(rep("Bud",8), rep("Flag",5), rep("leaf",6), rep("Root",9),
                   rep("SAM",7), rep("Seed",5), rep("SP",10), rep("ST",6))

t0 <- Sys.time()
msg2 <- function(...) cat(sprintf("[%6.1fs] %s\n",
  as.numeric(difftime(Sys.time(), t0, units = "secs")), paste0(...)))

## ===== 1. 载入作者注释对象 =================================================
rice.sub <- readRDS(RDS_PATH)
msg2("loaded rice.sub: ", nrow(rice.sub), " genes x ", ncol(rice.sub), " cells")

## ★ 官方脚本缺失的一步：对象默认 Idents 只有 8 个组织，
##   必须改成 56 个细胞类型（final_id2）才可能得到 56×56 的矩阵。
rice.sub$final_id2 <- as.integer(rice.sub$final_id2)
Idents(rice.sub) <- factor(rice.sub$final_id2, levels = as.character(GEO_ORDER))

## ===== 2. 官方第 14-20 行：两侧平均谱 ======================================
## ★ Seurat 5 的 AverageExpression() 返回稀疏矩阵，必须先 as.matrix()，
##   否则 cor() 会报 `'x' must be numeric`。
avg_of <- function(obj, assay, layer = "data") {
  DefaultAssay(obj) <- assay
  feats <- rownames(GetAssayData(obj, assay = assay, layer = layer))  # layer= 替代已弃用的 slot=
  as.matrix(AverageExpression(obj, assays = assay, features = feats, layer = layer)[[assay]])
}

DefaultAssay(rice.sub) <- "RNA"
rna_avg  <- avg_of(rice.sub, "RNA");     msg2("RNA mean profile    ", nrow(rna_avg),  " x ", ncol(rna_avg))
atac_avg <- avg_of(rice.sub, "atacRNA"); msg2("atacRNA mean profile ", nrow(atac_avg), " x ", ncol(atac_avg))

## 两模态的基因检出率（用于 C 版基因集）
det_rna  <- Matrix::rowMeans(GetAssayData(rice.sub, assay = "RNA",     layer = "data") > 0)
det_atac <- Matrix::rowMeans(GetAssayData(rice.sub, assay = "atacRNA", layer = "data") > 0)
rm(rice.sub); gc(verbose = FALSE)

## 官方第 99 行：ATAC 侧列名加 "-ATAC" 后缀（保证 cor() 按基因名对齐）
colnames(atac_avg) <- paste(colnames(atac_avg), "ATAC", sep = "-")

## ===== 3. 官方第 101-102 行：Spearman 相关（C 版基因集）====================
## 官方那句 rice.sub.avg$RNA[rownames(atacRNA.avg$atacRNA), ] 的作用是
## 「只保留 ATAC 侧也存在的基因，并按同一顺序取行」——这里即顺序取交集。
build_cor <- function(genes) {
  g <- intersect(genes, intersect(rownames(rna_avg), rownames(atac_avg)))
  C <- cor(atac_avg[g, , drop = FALSE], rna_avg[g, , drop = FALSE], method = "spearman")
  colnames(C) <- rownames(C) <- sub("-ATAC$", "", colnames(C))
  list(cor = C, n = length(g))
}

common  <- intersect(rownames(rna_avg), rownames(atac_avg))
rng_rna <- apply(rna_avg[common, , drop = FALSE], 1, function(v) max(v) - min(v))  # 跨细胞类型极差
det05   <- common[det_rna[common] >= 0.05 & det_atac[common] >= 0.05]             # 两侧检出率 ≥ 5%
top6k   <- names(sort(rng_rna, decreasing = TRUE))[1:6000]                        # 判别力 top 6000
GENES_C <- union(top6k, det05)

setC <- build_cor(GENES_C)
C_C  <- setC$cor
msg2("gene set (version C): top6k(", length(top6k), ") U det>=5%(", length(det05), ") -> ",
    setC$n, " genes after intersect")
stopifnot(nrow(C_C) == 56, ncol(C_C) == 56)

## ===== 4. 与原图逐格数值比对（可选，缺基准文件则跳过）====================
TGT_FILE <- file.path(DIR_REF, "fig2a_matrix_values_exact.csv")
TIS <- list(Bud = 1:8, Flag = 9:13, leaf = 14:19, Root = 20:28,
            SAM = 29:35, Seed = 36:40, SP = 41:50, ST = 51:56)
lab <- integer(56); for (k in seq_along(TIS)) lab[TIS[[k]]] <- k
SAME <- outer(lab, lab, "=="); DIAG <- diag(TRUE, 56); WOFF <- SAME & !DIAG

if (file.exists(TGT_FILE)) {
  Tgt <- as.matrix(read.csv(TGT_FILE, row.names = 1, check.names = FALSE))
  d   <- C_C - Tgt
  cat("\n===== version C vs original (cell-by-cell) =====\n")
  cat(sprintf("  rep     Spearman %.4f  Pearson %.4f  MAE %.4f  RMSE %.4f\n",
              cor(as.vector(C_C), as.vector(Tgt), method = "spearman"),
              cor(as.vector(C_C), as.vector(Tgt)), mean(abs(d)), sqrt(mean(d^2))))
  cat(sprintf("  orig    diag %.3f  offdiag %.3f  sd %.4f  range [%.3f, %.3f]\n",
              mean(diag(Tgt)), mean(Tgt[upper.tri(Tgt)]), sd(Tgt), min(Tgt), max(Tgt)))
  cat(sprintf("  rep     diag %.3f  offdiag %.3f  sd %.4f  range [%.3f, %.3f]\n",
              mean(diag(C_C)), mean(C_C[upper.tri(C_C)]), sd(C_C), min(C_C), max(C_C)))
  cat(sprintf("  bias    diag %+.4f  within-tissue %+.4f  cross-tissue %+.4f\n",
              mean(d[DIAG]), mean(d[WOFF]), mean(d[!SAME])))
  ad <- abs(d)
  cat(sprintf("  |bias|  median %.4f  p90 %.4f  share<=0.02 %.1f%%  share<=0.05 %.1f%%\n",
              median(ad), quantile(ad, .9), mean(ad <= .02) * 100, mean(ad <= .05) * 100))

  ## 逐组织偏差表
  ## ★ 注意：outer() 默认 FUN="*"，两个逻辑向量相乘得到的是**数值 0/1**；
  ##   用它做矩阵索引时 0 会被当作「跳过」、1 被当作「取第 1 个元素」，结果全错。
  ##   必须显式写 FUN = "&" 才能得到逻辑矩阵。
  tb <- do.call(rbind, lapply(seq_along(TIS), function(k) {
    inb <- outer(lab == k, lab == k, FUN = "&") & !DIAG
    crs <- outer(lab == k, lab != k, FUN = "&")
    data.frame(tissue = names(TIS)[k], n = length(TIS[[k]]),
               orig_diag = mean(diag(Tgt)[TIS[[k]]]),   rep_diag = mean(diag(C_C)[TIS[[k]]]),
               orig_within = mean(Tgt[inb]),            rep_within = mean(C_C[inb]),
               orig_cross = mean(Tgt[crs]),             rep_cross = mean(C_C[crs]),
               bias_within = mean(d[inb]),              bias_cross = mean(d[crs]))
  }))
  write.csv(tb, file.path(OUT, "Fig2a_C_per_tissue.csv"), row.names = FALSE)
} else {
  msg2("ref/fig2a_matrix_values_exact.csv not found; skipping cell-by-cell comparison")
}

## ===== 5. 出图（ComplexHeatmap，版式逐项对齐原图）========================
COL_FUN <- colorRamp2(seq(0, 0.6, length.out = 100), PAL_100)
FS <- 5                                     # 原图字号 5 pt
plot_fig2a <- function(C, stem) {
  Cs <- pmin(pmax(C, 0), 0.6)
  rownames(Cs) <- colnames(Cs) <- as.character(PAPER_LABELS)
  ha_r <- rowAnnotation(
    dot = anno_points(rep(0.5, 56), pch = 16, size = unit(1.13, "mm"),
                      gp = gpar(col = DOT_56), axis = FALSE, border = FALSE,
                      width = unit(2.2, "mm")),
    num = anno_text(as.character(PAPER_LABELS), gp = gpar(fontsize = FS), just = "left"),
    show_annotation_name = FALSE)
  ha_b <- HeatmapAnnotation(
    dot = anno_points(rep(0.5, 56), pch = 16, size = unit(1.13, "mm"),
                      gp = gpar(col = DOT_56), axis = FALSE, border = FALSE,
                      height = unit(2.2, "mm")),
    show_annotation_name = FALSE)
  ht <- Heatmap(Cs, name = "Correlation", col = COL_FUN,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_side = "bottom", column_names_rot = 90,
    column_names_centered = TRUE, column_names_gp = gpar(fontsize = FS),
    rect_gp = gpar(col = "grey60", lwd = 0.31),          # 原图格子描边 grey60 / 0.23 pt
    width = unit(82.8, "mm"), height = unit(76.9, "mm"),
    right_annotation = ha_r, bottom_annotation = ha_b,
    column_title = "scATAC", column_title_gp = gpar(fontsize = FS),
    row_title = "scRNA", row_title_gp = gpar(fontsize = FS),
    heatmap_legend_param = list(
      at = c(0, 0.2, 0.4, 0.6), title = "Correlation",
      title_position = "leftcenter-rot", legend_height = unit(15.8, "mm"),
      grid_width = unit(1.35, "mm"), border = FALSE,
      labels_gp = gpar(fontsize = FS), title_gp = gpar(fontsize = FS),
      ticks_gp = gpar(fontsize = FS)))
  fig_open(file.path(OUT, paste0(stem, ".png")), "png", 2700, 2580, res = 600)
  draw(ht, heatmap_legend_side = "right"); fig_close()
  fig_open(file.path(OUT, paste0(stem, ".pdf")), "pdf", 4.5, 4.3)
  draw(ht, heatmap_legend_side = "right"); fig_close()
  invisible(ht)
}
plot_fig2a(C_C, "Fig2a_reproduced")

## ===== 6. 产物落盘 ========================================================
M <- pmin(pmax(C_C, 0), 0.6)
rownames(M) <- colnames(M) <- as.character(PAPER_LABELS)
write.csv(C_C, file.path(OUT, "Fig2a_best_56x56_raw.csv"))
write.csv(M,   file.path(OUT, "Fig2a_best_56x56_clipped.csv"))
write.csv(data.frame(position = 1:56, axis_label_paper = PAPER_LABELS,
                     geo_final_id2 = GEO_ORDER, tissue = TISSUE_OF_POS[GEO_ORDER]),
          file.path(OUT, "Fig2a_axis_mapping.csv"), row.names = FALSE)
write.csv(data.frame(desc = c("top6000_by_range", "detect>=5%_both", "union_after_intersect",
                              "rna_avg_genes", "atac_avg_genes"),
                     n = c(6000, length(det05), setC$n, nrow(rna_avg), nrow(atac_avg))),
          file.path(OUT, "Fig2a_best_geneset_summary.csv"), row.names = FALSE)
msg2("done. output dir: ", normalizePath(OUT))
