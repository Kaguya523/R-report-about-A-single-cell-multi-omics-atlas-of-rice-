# =============================================================================
# 05_fig4uw_deg.R —— 复现 Fig. 4u / 4w（WT vs osltpl120 突变体，各细胞类型 DEG）
# -----------------------------------------------------------------------------
# 原图构成：
#   Fig.4u: UMAP + 每个细胞类型的 DEG 数（柱状/点图），突出 cortex 最多
#   Fig.4v: 皮层 DEG 数散点，标注 OsNRT1.1A / OsAAP1 / OsBHLH044
#   Fig.4w: 热图，上述 3 个基因在 cortex 与其他细胞类型中的表达（WT vs 突变体）
# 输入：GSE273875_LTP_NIP_label_240622.rds.gz（390 MB，已带 cell_type 注释）
# 成功判据：cortex 的 DEG 数在各细胞类型中最多
# 参考：官方脚本 15.Analyzing_scRNA_of_LTPL120.R
# =============================================================================

source("R/_bootstrap.R")   # 先设 UTF-8 locale，否则中文乱码
source("R/_common.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

RDS_LTP <- "data/GSE273875_LTP_NIP_label_240622.rds.gz"
if (!file.exists(RDS_LTP)) stop("找不到 ", RDS_LTP, "，请先运行 scripts/00_download.sh")

message("读取突变体对象（390 MB）...")
obj <- readRDS(gzfile(RDS_LTP))
print(table(Idents(obj)))

# ---- 1. 定位细胞类型列与条件列 ----------------------------------------------
ct_col <- guess_meta_col(obj@meta.data, "celltype")
if (is.null(ct_col)) stop("找不到 cell_type 列，请先跑 01_inspect_object.R")
message("细胞类型列: ", ct_col, " | 类型数: ", length(unique(obj@meta.data[[ct_col]])))

# 条件列：官方对象用 condition（"WT" / "LTP"），否则从 orig.ident 推断
cond_col <- intersect(c("condition", "group", "genotype"), colnames(obj@meta.data))
if (length(cond_col) > 0) {
  obj$condition2 <- as.character(obj@meta.data[[cond_col[1]]])
} else {
  obj$condition2 <- ifelse(grepl("^NIP|WT|Nip", obj@meta.data$orig.ident), "WT", "LTP")
}
obj$condition2 <- dplyr::recode(obj$condition2, "LTP" = "osltpl120", "NIP" = "WT", "Nip" = "WT")
print(table(obj$condition2))

obj$compare <- paste(obj$condition2, obj@meta.data[[ct_col]], sep = "_")
Idents(obj) <- "compare"
DefaultAssay(obj) <- "RNA"

# ---- 2. 逐细胞类型做差异分析（突变体 vs WT）---------------------------------
ct_levels <- sort(unique(as.character(obj@meta.data[[ct_col]])))
message("逐细胞类型差异分析，共 ", length(ct_levels), " 个...")

deg_list <- lapply(ct_levels, function(ct) {
  g1 <- paste("osltpl120", ct, sep = "_")
  g2 <- paste("WT", ct, sep = "_")
  if (!all(c(g1, g2) %in% levels(Idents(obj)))) {
    warning("跳过 ", ct, "（缺少分组: ", g1, " / ", g2, "）"); return(NULL)
  }
  message("  DEG: ", ct)
  m <- FindMarkers(obj, ident.1 = g1, ident.2 = g2,
                   logfc.threshold = 0.1, min.pct = 0.1, verbose = FALSE)
  m$gene <- rownames(m); m$cell_type <- ct
  m
})
deg_list <- deg_list[!sapply(deg_list, is.null)]
names(deg_list) <- ct_levels[!sapply(lapply(ct_levels, function(ct) {
  g1 <- paste("osltpl120", ct, sep = "_"); g2 <- paste("WT", ct, sep = "_")
  all(c(g1, g2) %in% levels(Idents(obj)))
}), isFALSE)]

saveRDS(deg_list, file.path(DIR_FIG, "fig4_deg_list.rds"))

# ---- 3. 显著 DEG 统计 --------------------------------------------------------
sig_list <- lapply(deg_list, function(x) {
  subset(x, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
})
deg_stat <- data.frame(
  cell_type = names(sig_list),
  n_up      = sapply(sig_list, function(x) sum(x$avg_log2FC > 0)),
  n_down    = sapply(sig_list, function(x) sum(x$avg_log2FC < 0)),
  row.names = NULL
)
deg_stat$n_total <- deg_stat$n_up + deg_stat$n_down
deg_stat <- deg_stat[order(-deg_stat$n_total), ]
cat("\n--- 各细胞类型显著 DEG 数（论文：cortex 最多）---\n"); print(deg_stat)
write.csv(deg_stat, file.path(DIR_FIG, "fig4u_deg_counts.csv"), row.names = FALSE)

if (nrow(deg_stat) > 0 && deg_stat$cell_type[1] != "cortex") {
  warning("⚠️ DEG 最多的细胞类型是 '", deg_stat$cell_type[1],
          "'，而论文是 'cortex'。请排查：注释列是否用错、条件分组是否反了、阈值是否一致。")
}

# ---- 4. Fig.4u：DEG 数柱状图 -------------------------------------------------
p_u <- ggplot(deg_stat, aes(x = reorder(cell_type, n_total), y = n_total)) +
  geom_col(fill = "#DC143C", colour = "black", linewidth = 0.25, width = 0.7) +
  geom_text(aes(label = n_total), hjust = -0.15, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(x = NULL, y = "Number of DEGs",
       title = "Fig. 4u reproduction — DEGs between WT and osltpl120 per cell type") +
  theme_pub()
save_fig(p_u, "Fig4u_deg_counts", width = 7.5, height = 5)

# ---- 5. Fig.4w：关键基因表达热图 --------------------------------------------
key_genes <- c("OsNRT1.1A", "OsAAP1", "OsBHLH044")
# 基因名在不同注释里可能带 LOC 前缀，做一次模糊匹配
all_genes <- rownames(obj)
resolve_gene <- function(g) {
  hit <- all_genes[toupper(all_genes) == toupper(g)]
  if (length(hit) == 0) hit <- all_genes[grepl(g, all_genes, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}
gene_map <- setNames(sapply(key_genes, resolve_gene), key_genes)
cat("\n--- 关键基因名映射 ---\n"); print(gene_map)

found <- gene_map[!is.na(gene_map)]
if (length(found) > 0) {
  obj$group_lab <- paste(obj$condition2, obj@meta.data[[ct_col]], sep = "_")
  Idents(obj) <- "group_lab"
  avg <- AverageExpression(obj, assays = "RNA", features = unname(found),
                           group.by = "group_lab", slot = "data")$RNA
  rownames(avg) <- names(gene_map)[match(rownames(avg), unname(gene_map))]

  hd <- as.data.frame(as.table(avg)); colnames(hd) <- c("gene", "group", "expr")
  hd <- hd %>% mutate(
    cond = ifelse(grepl("^osltpl120", group), "osltpl120", "WT"),
    ct   = sub("^(osltpl120|WT)_", "", group)
  )
  # 让 cortex 排在最前
  hd$ct <- factor(hd$ct, levels = c("cortex", setdiff(unique(hd$ct), "cortex")))

  p_w <- ggplot(hd, aes(x = cond, y = gene, fill = expr)) +
    geom_tile(colour = "grey85", linewidth = 0.3) +
    facet_wrap(~ ct, nrow = 1) +
    scale_fill_gradientn(colours = c("#EAEFD7", "#F4A582", "#B2182B"),
                         name = "Expression") +
    labs(x = NULL, y = NULL,
         title = "Fig. 4w reproduction — key root-development genes, WT vs osltpl120") +
    theme_pub(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(), strip.text = element_text(size = 7))
  save_fig(p_w, "Fig4w_key_gene_heatmap", width = 11, height = 3.6)
} else {
  warning("三个关键基因一个都没匹配上，请检查基因命名体系（MSU vs RAP-DB）")
}

message("\n完成。请把 output/fig4u_deg_counts.csv 与论文 Fig.4u 并排比对。")
