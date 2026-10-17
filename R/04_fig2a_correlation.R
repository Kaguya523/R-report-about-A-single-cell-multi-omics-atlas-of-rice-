# =============================================================================
# 04_fig2a_correlation.R —— 复现 Fig. 2a（scRNA 与 scATAC 的细胞类型间相关性）
# -----------------------------------------------------------------------------
# 原图构成：热图，行 = scATAC 细胞类型，列 = scRNA 细胞类型，色标 = Correlation (0-1)
#   论文原文："Correlations between scRNA-seq and scATAC for each cell type."
# 思路（改写自官方脚本 06，已修正其拼写错误与隐式依赖）：
#   1) scRNA：按细胞类型求平均表达 -> 两两 Pearson 相关
#   2) scATAC：先找每个细胞类型的 top 500 差异 peak -> 求平均可及性 -> 两两 Pearson 相关
#   3) 交叉：用 GeneActivity 把 peak 活性映射回基因，再与 scRNA 平均表达求 Spearman
# =============================================================================

source("R/_bootstrap.R")   # 先设 UTF-8 locale，否则中文乱码
source("R/_common.R")
suppressPackageStartupMessages({ library(Seurat); library(Signac) })

RDS_MAIN <- "data/GSE232863_scRNA_omics.Rds"
if (!file.exists(RDS_MAIN)) stop("找不到 ", RDS_MAIN)
obj <- readRDS(RDS_MAIN)

if (!all(c("RNA", "ATAC") %in% names(obj@assays)))
  stop("该对象缺少 RNA 或 ATAC assay，无法复现 Fig.2a。实际: ",
       paste(names(obj@assays), collapse = ", "))

ct_col <- guess_meta_col(obj@meta.data, "celltype")
if (is.null(ct_col)) ct_col <- guess_meta_col(obj@meta.data, "cluster")
if (is.null(ct_col)) stop("找不到细胞类型/簇列，请先跑 01_inspect_object.R")
message("使用细胞类型列: ", ct_col)
Idents(obj) <- obj@meta.data[[ct_col]]

n_celltype <- length(unique(Idents(obj)))
message("细胞类型数: ", n_celltype)

# ---- 1. RNA 平均表达 + 相关性 -----------------------------------------------
DefaultAssay(obj) <- "RNA"
feat_rna <- VariableFeatures(obj)
if (length(feat_rna) == 0) { obj <- FindVariableFeatures(obj, nfeatures = 2000); feat_rna <- VariableFeatures(obj) }

rna_avg <- AverageExpression(obj, assays = "RNA", features = feat_rna,
                             group.by = ct_col, slot = "data")$RNA
cor_rna <- cor(rna_avg, method = "pearson")

# ---- 2. ATAC 平均可及性 + 相关性 --------------------------------------------
DefaultAssay(obj) <- "ATAC"
# 找每个细胞类型的 top 500 差异 peak（用 LR 检验，控制测序深度）
atac_markers <- lapply(levels(Idents(obj)), function(x) {
  message("  ATAC FindMarkers: ", x)
  tryCatch({
    m <- FindMarkers(obj, ident.1 = x, only.pos = TRUE, test.use = "LR",
                     max.cells.per.ident = 300, latent.vars = "nCount_ATAC")
    data.frame(m, peak = rownames(m), cluster = x)
  }, error = function(e) { warning("跳过 ", x, ": ", conditionMessage(e)); NULL })
})
atac_markers <- do.call(rbind, atac_markers)
peaks_id <- atac_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 500) %>%
  pull(peak) %>% unique()
message("用于相关性的 peak 数: ", length(peaks_id))

atac_avg <- AverageExpression(obj, assays = "ATAC", features = peaks_id,
                              group.by = ct_col, slot = "data")$ATAC
cor_atac <- cor(atac_avg, method = "pearson")

# ---- 3. 交叉相关：ATAC gene activity vs RNA --------------------------------
message("计算 GeneActivity（较慢，视 peak 数而定）...")
gene_act <- GeneActivity(obj)
obj[["atacRNA"]] <- CreateAssayObject(counts = gene_act)
obj <- NormalizeData(obj, assay = "atacRNA", normalization.method = "LogNormalize",
                     scale.factor = median(obj$nCount_atacRNA))
act_avg <- AverageExpression(obj, assays = "atacRNA", features = feat_rna,
                             group.by = ct_col, slot = "data")$atacRNA
common <- intersect(rownames(act_avg), rownames(rna_avg))
cor_cross <- cor(act_avg[common, , drop = FALSE],
                 rna_avg[common, , drop = FALSE], method = "spearman")

# ---- 4. 出图（热图）--------------------------------------------------------
plot_heat <- function(m, title, name, show_x = TRUE, show_y = TRUE) {
  df <- as.data.frame(as.table(m)); colnames(df) <- c("ATAC", "RNA", "value")
  p <- ggplot(df, aes(RNA, ATAC, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                         limits = c(-1, 1), name = "Correlation") +
    labs(title = title, x = if (show_x) "scRNA (cell type)" else NULL,
         y = if (show_y) "scATAC (cell type)" else NULL) +
    theme_pub(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          axis.text.y = element_text(size = 5),
          panel.grid = element_blank(),
          aspect.ratio = 1)
  save_fig(p, name, width = 9, height = 8.5)
}

plot_heat(cor_rna,   "Fig. 2a (part) — RNA vs RNA, cell-type correlation", "Fig2a_RNA_correlation")
plot_heat(cor_atac,  "Fig. 2a (part) — ATAC vs ATAC, cell-type correlation", "Fig2a_ATAC_correlation")
plot_heat(cor_cross, "Fig. 2a — scATAC vs scRNA correlation per cell type", "Fig2a_ATAC_RNA_correlation")

write.csv(as.data.frame(as.table(cor_cross)), file.path(DIR_FIG, "fig2a_atac_rna_cor.csv"),
          row.names = FALSE)
saveRDS(list(cor_rna = cor_rna, cor_atac = cor_atac, cor_cross = cor_cross),
        file.path(DIR_FIG, "fig2a_cor_matrices.rds"))

cat("\n--- 交叉相关性摘要（应与论文一致：叶/旗叶的叶肉、维管、保卫细胞、表皮偏高）---\n")
print(summary(as.vector(cor_cross)))
message("\n完成。")
