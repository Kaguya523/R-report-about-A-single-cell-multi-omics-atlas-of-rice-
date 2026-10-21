###############################################################################
#  Fig. 1b 复现  ——  A single-cell multi-omics atlas of rice
#  Wang X. et al., Nature 644:722-729 (2025), DOI 10.1038/s41586-025-09251-0
#
#  图注原文：b, UMAP visualization integrating scRNA-seq and scATAC of all nuclei
#           (left). Nuclei are coloured by organ type. The barplots show the
#           number of nuclei for each organ (right).
#
#  官方源代码（已找到，不是自行瞎猜）：
#    github.com/dongwei-2023/Single_cell_multiomics_in_rice
#    └── 04.Integration_based_on_ATAC_and_RNA_data.R   ← README 注明对应 Figure 1b,1d
#         run_wnn(): RNA Normalize→HVG→Scale→PCA(50)→UMAP
#                    ATAC TFIDF→TopFeatures(q0)→SVD→UMAP(lsi 2:50)
#                    FindMultiModalNeighbors(list("pca","lsi"), list(1:50,2:50))
#                    RunUMAP(nn.name="weighted.nn") + FindClusters(wsnn, res=0.1)
#         末尾:      RunUMAP(nn.name="weighted.nn", dims=1:50, n.neighbors=30,
#                            n.epochs=200, min.dist=0.5, a=0.03, b=1.5)
#
#  本脚本：
#    1) 读 GSE232863 的 Seurat 对象（115,395 核，RNA + ATAC 双组学）
#    2) 取得 WNN 近邻图（作者存档的 weighted.nn，或由 pca+lsi 重算）
#    3) 用官方那组 RunUMAP 参数重算二维嵌入
#    4) 按 organ 着色画 UMAP + 各器官细胞数条形图，拼成 Fig.1b
#    5) 与对象中作者存档的 wsnnumap 做一致性检验
#
#  依赖：Seurat(>=5) Signac ggplot2 patchwork dplyr
#  运行：Rscript Fig1b_reproduce.R
###############################################################################

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})

## ======================= 0. 参数 =======================
# 路径：默认按仓库相对路径；数据放在别处时用环境变量覆盖，例如
#   Sys.setenv(RICE_RDS_SCRNA = "E:/data/GSE232863_scRNA_omics.Rds")
RDS_PATH <- Sys.getenv("RICE_RDS_SCRNA", "data/GSE232863_scRNA_omics.Rds")
OUT_DIR  <- Sys.getenv("RICE_OUT_FIG1B",  "output/fig1b")

# WNN 近邻图来源：
#   "authors"   —— 直接用对象里存好的 weighted.nn（默认，作者原始图就是用它出的）
#   "recompute" —— 用对象里的 pca / lsi 按官方参数重跑 FindMultiModalNeighbors
#                  （约 6 分钟；本机 15.7 GB 内存可跑通，见 README 的差异分析）
WNN_SOURCE <- "authors"
SEED       <- 42

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf(...), "\n")

## 论文 Fig.1b 的器官顺序 / 配色 / 显示名
## tissue 字段的 8 个水平与论文器官名的对应关系，由细胞数逐一比对确认
CODES  <- c("Root", "ST", "leaf", "Flag", "SAM", "Bud", "SP", "Seed")
ORGANS <- c("Crown root", "Culm", "Leaf", "Flag leaf",
            "Shoot apex", "Tiller bud", "Panicle", "Seed")
## 配色：由论文 PDF 第 2 页 Fig.1b 柱子的矢量填充色反推（见 extract_fig1b_ref.py）
PAL <- c(Root = "#6BBC6C", ST = "#DF61C7", leaf = "#8AAB28", Flag = "#BD9423",
         SAM  = "#6BBDC1", Bud = "#DB756D", SP   = "#B37BFA", Seed = "#5FA6FB")
LAB <- setNames(ORGANS, CODES)

## ======================= 1. 读入并瘦身 =======================
t0  <- Sys.time()
obj <- readRDS(RDS_PATH)
msg("readRDS: %.1f min, %d nuclei",
    as.numeric(difftime(Sys.time(), t0, units = "min")), ncol(obj))

# 原始对象 5.1 GB，读入峰值约 10 GB。ATAC(82,219 peaks)/atacRNA/chromvar
# 对 Fig.1b 无用，先释放，避免后面建图/UMAP 时 OOM
DefaultAssay(obj) <- "RNA"
for (a in c("ATAC", "atacRNA", "chromvar")) if (a %in% Assays(obj)) obj[[a]] <- NULL
invisible(gc(full = TRUE))
msg("assays kept: %s", paste(Assays(obj), collapse = ", "))

## ======================= 2. 取 WNN 近邻图 =======================
set.seed(SEED)

auth_nn <- if ("weighted.nn" %in% names(obj@neighbors)) obj@neighbors$weighted.nn else NULL
auth_w  <- data.frame(RNA = obj$RNA.weight, ATAC = obj$ATAC.weight)

if (WNN_SOURCE == "recompute") {
  stopifnot(all(c("pca", "lsi") %in% Reductions(obj)))
  t1  <- Sys.time()
  obj <- FindMultiModalNeighbors(obj,
                                 reduction.list = list("pca", "lsi"),
                                 dims.list      = list(1:50, 2:50))
  msg("FindMultiModalNeighbors: %.1f min",
      as.numeric(difftime(Sys.time(), t1, units = "min")))
  # 与作者存档的模态权重比对（FindMultiModalNeighbors 会覆盖 RNA.weight/ATAC.weight）
  wchk <- data.frame(
    metric = c("pearson_r_RNA.weight", "pearson_r_ATAC.weight",
               "MAE_RNA.weight", "MAE_ATAC.weight"),
    value  = c(cor(obj$RNA.weight, auth_w$RNA), cor(obj$ATAC.weight, auth_w$ATAC),
               mean(abs(obj$RNA.weight - auth_w$RNA)),
               mean(abs(obj$ATAC.weight - auth_w$ATAC)))
  )
  write.csv(wchk, file.path(OUT_DIR, "Fig1b_weights_check.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  msg("weights vs authors: r_RNA=%.3f r_ATAC=%.3f", wchk$value[1], wchk$value[2])
} else {
  ## 复用作者存档的 weighted.nn，保持 obj@neighbors 不被改动
  stopifnot(!is.null(auth_nn))
  msg("reuse authors' stored weighted.nn")
}

## ======================= 3. 官方参数重算二维 UMAP =======================
## 与官方脚本最后一行的参数完全一致。
## ⚠ 兼容性：官方脚本写的是 RunUMAP(..., nn.name="weighted.nn", dims=1:50, ...)，
##   Seurat 5 起 nn.name 与 dims 不能同时给（报 "Only one parameter among
##   'dims','nn.name','graph','features' should be used"），故此处去掉 dims；
##   nn.name 模式下 dims 本来就不参与计算，结果不受影响。
t2  <- Sys.time()
obj <- RunUMAP(obj,
               nn.name        = "weighted.nn",
               reduction.name = "wnn.umap",
               reduction.key  = "wnnUMAP_",
               n.neighbors    = 30,
               n.epochs       = 200,
               min.dist       = 0.5,
               a              = 0.03,
               b              = 1.5,
               seed.use       = SEED)
msg("RunUMAP: %.1f min", as.numeric(difftime(Sys.time(), t2, units = "min")))

emb <- Embeddings(obj, "wnn.umap")[, 1:2]
colnames(emb) <- c("UMAP_1", "UMAP_2")

## ======================= 4. 每个器官的细胞数 =======================
df <- data.frame(cell   = rownames(emb),
                 UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                 tissue = factor(as.character(obj$tissue), levels = CODES),
                 stringsAsFactors = FALSE)

cnt <- df %>% count(tissue, .drop = FALSE) %>%
  mutate(organ = factor(LAB[as.character(tissue)], levels = rev(ORGANS)),
         ynum  = as.numeric(organ))

## 论文报告的各器官细胞数（Fig.1b 柱子的数值标签，PDF 文本层直接读出）
PAPER_N <- c("Crown root" = 31108, "Culm" = 19740, "Leaf" = 13771,
             "Flag leaf" = 5563, "Shoot apex" = 13791, "Tiller bud" = 8858,
             "Panicle" = 17457, "Seed" = 6276)

cmp <- data.frame(order = seq_along(ORGANS), organ = ORGANS, tissue_code = CODES,
                  n_reproduced = cnt$n[match(CODES, as.character(cnt$tissue))])
cmp$n_paper <- as.integer(PAPER_N[cmp$organ])
cmp$diff    <- cmp$n_reproduced - cmp$n_paper
write.csv(cmp, file.path(OUT_DIR, "Fig1b_cell_counts.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
msg("total nuclei = %d (paper 116,564)", sum(cnt$n))

## ======================= 5. 出图 =======================
YLIM <- c(-1.95, 8.85)
YB   <- 0.45

## --- 5.1 左：按器官着色的 UMAP ---
rx <- range(df$UMAP_1); ry <- range(df$UMAP_2)
ax <- rx[1] + 0.01 * diff(rx); ay <- ry[1] + 0.03 * diff(ry)

p_umap <- ggplot(df, aes(UMAP_1, UMAP_2, colour = tissue)) +
  geom_point(size = 0.2, stroke = 0, shape = 16) +
  scale_colour_manual(values = PAL, guide = "none") +
  coord_fixed(expand = FALSE, clip = "off") +
  annotate("segment", x = ax, xend = ax + 0.15 * diff(rx), y = ay, yend = ay,
           arrow = arrow(length = unit(0.13, "cm"), type = "closed"), linewidth = 0.4) +
  annotate("segment", x = ax, xend = ax, y = ay, yend = ay + 0.15 * diff(ry),
           arrow = arrow(length = unit(0.13, "cm"), type = "closed"), linewidth = 0.4) +
  annotate("text", x = ax - 0.05 * diff(rx), y = ay + 0.10 * diff(ry),
           label = "UMAP2", angle = 90, size = 3.1) +
  annotate("text", x = ax + 0.10 * diff(rx), y = ay - 0.045 * diff(ry),
           label = "UMAP1", size = 3.1) +
  theme_void() + theme(plot.margin = margin(6, 6, 6, 6))

## --- 5.2 中：器官图例（与条形图共用同一连续 y 轴，保证逐行对齐）---
p_leg <- ggplot(cnt, aes(x = 0, y = ynum, colour = tissue)) +
  geom_point(size = 3.0) +
  geom_text(aes(label = organ), x = 0.16, hjust = 0, size = 3.2, colour = "black") +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_y_continuous(limits = YLIM, expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.10, 1)) +
  theme_void() + theme(plot.margin = margin(6, -6, 6, 0))

## --- 5.3 右：细胞数条形图（手绘坐标轴，与原图一致只画到 30,000）---
p_bar <- ggplot(cnt, aes(x = n, y = ynum, fill = tissue)) +
  geom_col(width = 0.5, orientation = "y") +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = 0, nudge_x = 800,
            size = 3.2, colour = "grey25") +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_y_continuous(limits = YLIM, expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 35600), expand = c(0, 0)) +
  annotate("segment", x = 0, xend = 30000, y = YB, yend = YB, linewidth = 0.4) +
  annotate("segment", x = c(0, 10000, 20000, 30000),
           xend = c(0, 10000, 20000, 30000),
           y = YB, yend = YB - 0.16, linewidth = 0.3) +
  annotate("text", x = c(0, 10000, 20000, 30000), y = YB - 0.26,
           label = c("0", "10,000", "20,000", "30,000"),
           angle = 90, size = 2.75, hjust = 1, vjust = 0.5) +
  annotate("text", x = 15000, y = -1.05, label = "Cell number", size = 3.2) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_void() + theme(plot.margin = margin(6, 16, 6, 4))

fig <- (p_umap | p_leg | p_bar) + plot_layout(widths = c(2.9, 0.56, 1.06))

ggsave(file.path(OUT_DIR, "Fig1b_reproduced.png"), fig,
       width = 8.6, height = 3.9, dpi = 330, bg = "white")
ggsave(file.path(OUT_DIR, "Fig1b_reproduced.pdf"), fig,
       width = 8.6, height = 3.9, bg = "white")
msg("figure saved")

## ======================= 6. 与作者存档 wsnnumap 的一致性 =======================
if ("wsnnumap" %in% Reductions(obj)) {
  ref <- Embeddings(obj, "wsnnumap")[, 1:2]
  norm <- function(m) apply(m, 2, function(v) (v - min(v)) / max(1e-12, diff(range(v))))
  a <- norm(emb); b <- norm(ref)
  nb <- 60
  bin <- function(m) {
    i <- pmin(nb, pmax(1, floor(m[, 1] * (nb - 1)) + 1))
    j <- pmin(nb, pmax(1, floor(m[, 2] * (nb - 1)) + 1))
    i + (j - 1) * nb
  }
  ia <- bin(a); ib <- bin(b); tis <- as.character(df$tissue)
  dens <- function(idx) sapply(CODES, function(k) as.numeric(tabulate(idx[tis == k], nbins = nb * nb)))
  da <- dens(ia); db <- dens(ib)
  cors <- sapply(CODES, function(k) cor(sqrt(da[, k]), sqrt(db[, k])))
  chk <- data.frame(organ = ORGANS, density_r = round(cors, 3))
  write.csv(chk, file.path(OUT_DIR, "Fig1b_consistency_vs_authors.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  msg("consistency vs authors' wsnnumap: mean density r = %.3f", mean(cors))
  print(chk)
}

## ======================= 7. 落盘 =======================
write.csv(df, file.path(OUT_DIR, "Fig1b_umap_coords.csv"), row.names = FALSE, quote = FALSE)
msg("done: %s", format(Sys.time()))
