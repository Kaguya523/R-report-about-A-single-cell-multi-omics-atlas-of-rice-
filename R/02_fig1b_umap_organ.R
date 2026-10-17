# =============================================================================
# 02_fig1b_umap_organ.R —— 复现 Fig. 1b
# -----------------------------------------------------------------------------
# 原图构成：左侧 = 全细胞 UMAP（按器官着色）；右侧 = 各器官细胞数横向柱状图
# 输入数据：GSE232863_scRNA_omics.Rds
# 成功判据：8 器官细胞数相加 = 116,564
#   Crown root 31,108 / Culm 19,740 / Leaf 13,771 / Flag leaf 5,563
#   Shoot apex 13,791 / Tiller bud 8,858 / Panicle 17,457 / Seed 6,276
# =============================================================================

source("R/_bootstrap.R")   # 先设 UTF-8 locale，否则中文乱码
source("R/_common.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

RDS_MAIN <- "data/GSE232863_scRNA_omics.Rds"
if (!file.exists(RDS_MAIN)) stop("找不到 ", RDS_MAIN, "，请先运行 scripts/00_download.sh")

message("读取对象（5.1 GB，需要 ≥32 GB 内存，请耐心等待）...")
obj <- readRDS(RDS_MAIN)

# ---- 1. 定位 UMAP 与器官列 ---------------------------------------------------
red_use <- intersect(c("wsnnumap", "umap", "wnn.umap", "harmony"), names(obj@reductions))[1]
if (is.na(red_use)) stop("对象里没有 UMAP 降维结果，可用: ",
                         paste(names(obj@reductions), collapse = ", "))
message("使用降维: ", red_use)

organ_col <- guess_meta_col(obj@meta.data, "organ")
if (is.null(organ_col)) stop("找不到器官列，请先跑 01_inspect_object.R 查看真实列名")
message("使用器官列: ", organ_col)

# ---- 2. 统一器官名（命名归一化函数在 _common.R 里）--------------------------
emb  <- as.data.frame(Embeddings(obj, red_use))
colnames(emb)[1:2] <- c("UMAP1", "UMAP2")
emb$organ <- norm_organ(obj@meta.data[[organ_col]])

# ---- 3. 细胞数统计 + 与论文核对 ----------------------------------------------
cnt <- emb %>% count(organ, name = "n")
cnt$paper <- PAPER_NUCLEI[cnt$organ]
cat("\n--- 各器官细胞数（本复现 vs 论文）---\n")
print(cnt)
check_total(nrow(emb))
if (sum(!is.na(cnt$paper)) > 0) {
  cnt$diff <- cnt$n - cnt$paper
  cat("\n差值:\n"); print(cnt[, c("organ", "n", "paper", "diff")])
}

# 保持论文顺序
lv <- names(PAPER_NUCLEI)
cnt$organ <- factor(cnt$organ, levels = lv)
emb$organ <- factor(emb$organ, levels = lv)

# ---- 4. 左：UMAP -------------------------------------------------------------
p_umap <- ggplot(emb, aes(UMAP1, UMAP2, colour = organ)) +
  geom_point(size = 0.06, stroke = 0, alpha = 0.7) +
  scale_color_manual(values = ORGAN_COLORS, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
  labs(x = "UMAP1", y = "UMAP2", title = "UMAP of 116,564 nuclei, coloured by organ") +
  theme_pub() + theme(legend.position = "right")

# ---- 5. 右：细胞数柱状图 -----------------------------------------------------
p_bar <- ggplot(cnt, aes(x = organ, y = n, fill = organ)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.12, size = 3.1) +
  scale_fill_manual(values = ORGAN_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_flip() +
  labs(x = NULL, y = "Number of nuclei",
       title = "Nuclei per organ") +
  theme_pub()

# ---- 6. 拼图（用 patchwork；若未安装则退回 cowplot）-------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_all <- patchwork::wrap_plots(p_umap, p_bar, widths = c(1.55, 1)) +
    patchwork::plot_annotation(
      title = "Fig. 1b reproduction | A single-cell multi-omics atlas of rice (Nature 644:722, 2025)",
      theme = theme(plot.title = element_text(face = "bold", size = 12)))
} else {
  message("未安装 patchwork，改用 cowplot::plot_grid")
  p_all <- cowplot::plot_grid(p_umap, p_bar, nrow = 1, rel_widths = c(1.55, 1))
}

save_fig(p_all, "Fig1b_umap_by_organ", width = 13, height = 5.6)

# 单独存 UMAP 数据，便于 Fig.1d/1e 复用，避免反复读 5 GB 对象
saveRDS(emb, file.path(DIR_FIG, "fig1b_cell_embedding.rds"))
write.csv(cnt, file.path(DIR_FIG, "fig1b_nuclei_per_organ.csv"), row.names = FALSE)
message("\n完成。中间结果已存 output/fig1b_cell_embedding.rds（后续脚本直接读它，省内存）")
