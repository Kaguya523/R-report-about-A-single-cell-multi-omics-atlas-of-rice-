## ============================================================================
##  Fig. 1c 复现  --  Single-cell multi-omics vs. paired scATAC-seq + scRNA-seq
##  Wang X. et al., Nature 644:722-729 (2025)  DOI: 10.1038/s41586-025-09251-0
## ----------------------------------------------------------------------------
##  论文 Methods（"Comparison of multi-omics and paired data"）原文只给了流程：
##    1) 用 Signac/Seurat 处理同时测 RNA 与 ATAC 的 multiome 数据，得到"真值"注释；
##    2) 把同一批细胞拆成两个"伪模态"：
##         scRNA 端 —— 按 Seurat 标准流程（NormalizeData→HVG→ScaleData→RunPCA）
##         scATAC 端 —— 按 gene body + promoter 统计每个基因的 ATAC 计数，得到
##                      "基因活性矩阵"（atacRNA assay），随后 NormalizeData + RunPCA。
##         【重要澄清】Methods 中的 RunTFIDF/FindTopFeatures/RunSVD 是 ATAC *peak*
##          矩阵的降维（用于 ATAC 侧聚类/UMAP），并不参与基因活性矩阵的标签转移；
##          基因活性矩阵按官方 pbmc_multiomic 教程做 LogNormalize，而非 TF-IDF/LSI。
##    3) FindTransferAnchors + TransferData 把 scRNA 的注释转移到 scATAC 细胞；
##    4) 与 multiome 注释逐细胞比对，算 accuracy。
##  原文未给出：HVG 个数、dims、anchor 的 reduction、基因活性窗口宽度。本脚本
##  采用 Seurat 默认量级（2000 HVG / 1:30 dims / CCA anchors），见下方 PARAMS。
## ----------------------------------------------------------------------------
##  输入：GSE232863_scRNA_omics.Rds（GEO 公开，含 RNA / ATAC / atacRNA 三个 assay）
##        atacRNA = 作者已算好的"基因活性"矩阵（gene body + promoter），本脚本直接使用
##  输出：results/Fig1c_reproduced.csv        逐器官 accuracy
##        results/Fig1c_reproduced.png/.pdf   复现图
##        results/Fig1c_confusion_<organ>.csv 逐器官混淆矩阵（失败环节分析用）
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix)
  library(ggplot2); library(ggpattern); library(dplyr)
})
options(future.globals.maxSize = 16000 * 1024^2)
set.seed(42)

## ------------------------------- 路径 --------------------------------------
# 默认按仓库相对路径；数据放在别处时用环境变量覆盖，例如
#   Sys.setenv(RICE_RDS_SCRNA = "E:/data/GSE232863_scRNA_omics.Rds")
RDS <- Sys.getenv("RICE_RDS_SCRNA", "data/GSE232863_scRNA_omics.Rds")
OUT <- Sys.getenv("RICE_OUT_FIG1C", "output/fig1c/results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ------------------------------ 参数 ---------------------------------------
PARAMS <- list(
  nfeatures = 2000,      # scRNA 端高变基因数
  npcs      = 30,        # PCA 维数 / 转移使用 1:30
  reduction = "cca",     # FindTransferAnchors 的 reduction
  k.anchor  = 5,         # 默认
  k.weight  = 50         # TransferData 默认邻域权重数
)

## 8 个器官：tissue 字段 -> 图上名称 -> 原图取色（从 PDF 矢量逐条取色）
ORGANS <- data.frame(
  tissue = c("Root", "ST",   "leaf", "Flag", "SAM",   "Bud",   "SP",    "Seed"),
  organ  = c("Crown root", "Culm", "Leaf", "Flag leaf", "Shoot apex",
             "Tiller bud", "Panicle", "Seed"),
  hex    = c("#6BBC6C", "#DF61C7", "#8AAB28", "#BD9423",
             "#6BBDC1", "#DB756D", "#B37BFA", "#5FA6FB"),
  stringsAsFactors = FALSE
)

## ============================ 1. 读入 ======================================
cat("[1/4] reading multiome object ...\n")
t0 <- Sys.time()
obj <- readRDS(RDS)
## 内存优化：本研究只用 RNA（参考）与 atacRNA（基因活性，查询）两套 assay，
## 读入后立即丢弃 82,219 峰的 ATAC 矩阵与 chromvar，可省下约一半内存。
DefaultAssay(obj) <- "RNA"
obj[["ATAC"]] <- NULL; obj[["chromvar"]] <- NULL
invisible(gc(verbose = FALSE))
cat("      cells:", ncol(obj), " assays:", paste(Assays(obj), collapse = "/"),
    " 用时", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

## ============================ 2. 逐器官标签转移 =============================
cat("[2/4] label transfer scRNA -> scATAC, per organ ...\n")
res <- list(); conf_list <- list()

for (i in seq_len(nrow(ORGANS))) {
  org <- ORGANS$tissue[i]; name <- ORGANS$organ[i]
  cells <- colnames(obj)[obj$tissue == org]
  if (length(cells) < 100) { cat("   skip", name, "\n"); next }
  t1 <- Sys.time()

  ## ---- 2a. scRNA 端（参考集）----
  rna <- subset(obj, cells = cells)
  DefaultAssay(rna) <- "RNA"
  rna[["ATAC"]] <- NULL; rna[["atacRNA"]] <- NULL; rna[["chromvar"]] <- NULL
  rna <- NormalizeData(rna, verbose = FALSE)                       # LogNormalize
  rna <- FindVariableFeatures(rna, nfeatures = PARAMS$nfeatures, verbose = FALSE)
  rna <- ScaleData(rna, features = VariableFeatures(rna), verbose = FALSE)
  rna <- RunPCA(rna, features = VariableFeatures(rna),
                npcs = PARAMS$npcs, verbose = FALSE)
  rna$label <- rna$tissue_cluster                                  # multiome 注释

  ## ---- 2b. scATAC 端（查询集）：基因活性矩阵 ----
  atac <- subset(obj, cells = cells)
  DefaultAssay(atac) <- "atacRNA"     # = gene body + promoter 的 ATAC 计数
  atac[["RNA"]] <- NULL; atac[["chromvar"]] <- NULL; atac[["ATAC"]] <- NULL
  atac <- NormalizeData(atac, verbose = FALSE)   # LogNormalize（官方 vignette 做法）
  ## 两侧行名同为 Os 位点号，HVG 交集约 98%（1952–1961 / 2000）；取交集只为保证
  ## FindTransferAnchors 两侧特征行完全对齐（少数 RNA 高变基因在 atacRNA 无对应行）。
  feat <- intersect(VariableFeatures(rna), rownames(atac))         # 两侧共同基因
  atac <- ScaleData(atac, features = feat, verbose = FALSE)
  atac <- RunPCA(atac, features = feat, npcs = PARAMS$npcs, verbose = FALSE)
  atac$label <- atac$tissue_cluster

  ## ---- 2c. 锚点 + 标签转移 ----
  dims <- seq_len(PARAMS$npcs)
  anchors <- FindTransferAnchors(reference = rna, query = atac, features = feat,
                                 reduction = PARAMS$reduction, dims = dims,
                                 k.anchor = PARAMS$k.anchor, verbose = FALSE)
  pred <- TransferData(anchorset = anchors, refdata = rna$label, dims = dims,
                       weight.reduction = atac[["pca"]],
                       k.weight = PARAMS$k.weight, verbose = FALSE)

  ## ---- 2d. 与真值比对 ----
  truth <- as.character(atac$label)
  call  <- as.character(pred$predicted.id)
  acc   <- mean(call == truth)
  lv    <- sort(unique(truth))
  cm    <- table(factor(truth, levels = lv), factor(call, levels = lv))
  write.csv(as.data.frame.matrix(cm),
            file.path(OUT, paste0("Fig1c_confusion_", org, ".csv")))
  conf_list[[org]] <- cm

  res[[i]] <- data.frame(order = i, organ = name, tissue = org, hex = ORGANS$hex[i],
                         n_cells = length(cells), n_labels = length(lv),
                         accuracy = round(acc * 100, 2))
  cat(sprintf("   %-11s n=%5d  labels=%2d  accuracy=%.2f%%  (%.1f min)\n",
              name, length(cells), length(lv), acc * 100,
              as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  rm(rna, atac, anchors, pred); gc(verbose = FALSE)
}

df <- do.call(rbind, res)
write.csv(df, file.path(OUT, "Fig1c_reproduced.csv"), row.names = FALSE)
rm(obj); gc(verbose = FALSE)

## ============================ 3. 出图 ======================================
cat("[3/4] plotting ...\n")
n <- nrow(df); df$y <- rev(seq_len(n))
mu <- mean(df$accuracy)
half <- 0.195          # 半条柱高（行距 = 1，对应原图 4.15pt/5.19pt ≈ 0.80）
gap  <- 0.022          # 上下两段之间的留白

solid <- data.frame(xmin = 0, xmax = df$accuracy,
                    ymin = df$y - half, ymax = df$y - gap, hex = df$hex)
hatch <- data.frame(xmin = 0, xmax = 100,
                    ymin = df$y + gap, ymax = df$y + half, hex = df$hex)

p <- ggplot() +
  ## 上段：斜线填充（Single-cell multi-omics，恒 100%）
  geom_rect_pattern(
    data = hatch,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = hex),
    colour = "black", linewidth = 0.26,
    pattern = "stripe", pattern_angle = 45, pattern_density = 0.55,
    pattern_spacing = 0.009, pattern_fill = "black", pattern_colour = "black",
    pattern_size = 0.16
  ) +
  ## 下段：纯色填充（Paired scATAC-seq + scRNA-seq 的 accuracy）
  geom_rect(data = solid,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = hex),
            colour = "black", linewidth = 0.26) +
  scale_fill_identity() +
  ## 均值虚线
  geom_segment(aes(x = mu, xend = mu, y = 0.60, yend = n + 0.20),
               linetype = "22", linewidth = 0.4, colour = "grey25") +
  ## 右侧星号
  geom_text(data = df, aes(x = 101.8, y = y, label = "*"),
            hjust = 0, size = 4.4, fontface = "bold") +
  scale_x_continuous(limits = c(0, 106), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = df$y, labels = df$organ,
                     limits = c(0.45, n + 1.85), expand = c(0, 0)) +
  labs(x = "Accuracy (%)", y = NULL, tag = "c") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    axis.line.x  = element_line(linewidth = 0.45),
    axis.line.y  = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.45),
    axis.ticks.length.x = unit(4, "pt"),
    axis.ticks.y = element_line(linewidth = 0.45),
    axis.ticks.length.y = unit(4, "pt"),
    axis.text.x  = element_text(colour = "black", size = 11, margin = margin(t = 3)),
    axis.text.y  = element_text(colour = "black", size = 11, hjust = 1, margin = margin(r = 4)),
    axis.title.x = element_text(colour = "black", size = 12, margin = margin(t = 7)),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.tag = element_text(face = "bold", size = 19),
    plot.tag.position = c(0.005, 1.02),
    plot.margin = margin(t = 42, r = 24, b = 6, l = 6)
  )

## 手工图例（与原图同位置：面板号右侧上方）
leg1 <- n + 1.50; leg2 <- n + 0.96
p <- p +
  geom_rect_pattern(data = data.frame(xmin = 12.5, xmax = 16.2,
                                      ymin = leg1 - 0.16, ymax = leg1 + 0.16),
                    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                    fill = "white", colour = "black", linewidth = 0.3,
                    pattern = "stripe", pattern_angle = 45, pattern_density = 0.55,
                    pattern_spacing = 0.009, pattern_fill = "black",
                    pattern_colour = "black", pattern_size = 0.16) +
  annotate("text", x = 19.5, y = leg1, label = "Single-cell multi-omics",
           hjust = 0, size = 3.7) +
  annotate("rect", xmin = 12.5, xmax = 16.2, ymin = leg2 - 0.16, ymax = leg2 + 0.16,
           fill = "white", colour = "black", linewidth = 0.3) +
  annotate("text", x = 19.5, y = leg2, label = "Paired scATAC-seq and scRNA-seq",
           hjust = 0, size = 3.7)

ggsave(file.path(OUT, "Fig1c_reproduced.png"), p, width = 5.0, height = 3.6,
       dpi = 400, bg = "white")
ggsave(file.path(OUT, "Fig1c_reproduced.pdf"), p, width = 5.0, height = 3.6,
       bg = "white")

## ============================ 4. 摘要 ======================================
cat("[4/4] done\n")
print(df[, c("organ", "n_cells", "n_labels", "accuracy")], row.names = FALSE)
cat(sprintf("mean accuracy = %.2f%%\n", mu))
