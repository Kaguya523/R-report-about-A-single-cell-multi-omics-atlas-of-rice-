# =============================================================================
# 03_fig1e_composition.R —— 复现 Fig. 1e（各器官中细胞类型的占比）
# -----------------------------------------------------------------------------
# 原图构成：横向堆叠柱状图，每行一个器官，x 轴为占比 0-100%
#   论文原文："The proportion of cell numbers for each cell type in different organs.
#              The colour represents the cell type as shown in panel d."
# 输入：优先读 02 脚本存下的 fig1b_cell_embedding.rds；没有则读原始对象
# =============================================================================

source("R/_bootstrap.R")   # 先设 UTF-8 locale，否则中文乱码
source("R/_common.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

EMB_RDS <- "output/fig1b_cell_embedding.rds"

if (file.exists(EMB_RDS)) {
  message("读取 ", EMB_RDS, "（轻量，推荐）")
  emb <- readRDS(EMB_RDS)
} else {
  RDS_MAIN <- "data/GSE232863_scRNA_omics.Rds"
  if (!file.exists(RDS_MAIN)) stop("找不到输入，请先跑 R/02_fig1b_umap_organ.R")
  obj <- readRDS(RDS_MAIN)
  red_use <- intersect(c("wsnnumap", "umap", "wnn.umap", "harmony"), names(obj@reductions))[1]
  emb <- as.data.frame(Embeddings(obj, red_use)); colnames(emb)[1:2] <- c("UMAP1", "UMAP2")
  emb$organ <- norm_organ(obj@meta.data[[guess_meta_col(obj@meta.data, "organ")]])
}

# 需要细胞类型列 —— 若 emb 里没有，则重新从对象里取
if (!"celltype" %in% colnames(emb)) {
  RDS_MAIN <- "data/GSE232863_scRNA_omics.Rds"
  if (!file.exists(RDS_MAIN)) stop("需要从原始对象读取细胞类型列，但找不到 ", RDS_MAIN)
  message("从原始对象读取细胞类型列...")
  obj <- readRDS(RDS_MAIN)
  ct_col <- guess_meta_col(obj@meta.data, "celltype")
  if (is.null(ct_col)) stop("找不到细胞类型列，请先跑 01_inspect_object.R")
  emb$celltype <- as.character(obj@meta.data[[ct_col]])
  saveRDS(emb, EMB_RDS)
}

# ---- 统计占比 ---------------------------------------------------------------
comp <- emb %>%
  filter(!is.na(organ), !is.na(celltype)) %>%
  count(organ, celltype, name = "n") %>%
  group_by(organ) %>%
  mutate(prop = n / sum(n) * 100) %>%
  ungroup()

# 只保留在≥1个器官中占比 >0.5% 的细胞类型，否则图例太长（论文图例是两栏长列表）
keep_ct <- comp %>% group_by(celltype) %>% summarise(m = max(prop)) %>%
  filter(m >= 0.5) %>% pull(celltype)
comp_f <- comp %>% filter(celltype %in% keep_ct)

# 细胞类型过多，用调色板循环生成
n_ct <- length(unique(comp_f$celltype))
pal <- setNames(rep_len(c(ggsci::pal_igv("default")(33),
                          scales::hue_pal()(n_ct),
                          RColorBrewer::brewer.pal(12, "Paired")), n_ct),
                unique(comp_f$celltype))

comp_f$organ <- factor(comp_f$organ, levels = names(PAPER_NUCLEI))

p <- ggplot(comp_f, aes(x = organ, y = prop, fill = celltype)) +
  geom_col(width = 0.75, colour = "grey20", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(x = NULL, y = "Proportion of nuclei (%)",
       title = "Fig. 1e reproduction — cell-type composition per organ") +
  theme_pub() +
  theme(legend.position = "right",
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(0.28, "cm"),
        legend.spacing.y = unit(0.02, "cm"))

save_fig(p, "Fig1e_celltype_composition", width = 12, height = 6.5)
write.csv(comp, file.path(DIR_FIG, "fig1e_composition_table.csv"), row.names = FALSE)

# 顺手出 Fig.1d：全细胞 UMAP 按细胞类型着色
if (requireNamespace("ggrepel", quietly = TRUE)) {
  cen <- emb %>% filter(!is.na(celltype)) %>%
    group_by(celltype) %>% summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
  p_d <- ggplot(emb, aes(UMAP1, UMAP2, colour = celltype)) +
    geom_point(size = 0.06, stroke = 0, alpha = 0.6) +
    ggrepel::geom_text_repel(data = cen, aes(label = celltype), size = 2.2, max.overlaps = 60) +
    scale_color_manual(values = pal, guide = "none") +
    labs(title = "Fig. 1d reproduction — UMAP coloured by cell type") +
    theme_pub()
  save_fig(p_d, "Fig1d_umap_celltype", width = 8, height = 7)
}
message("\n完成。")
