#!/usr/bin/env Rscript
# =============================================================================
# Fig2b_TF_enrichment.R  —— Nature 644, 724-729  Fig. 2b 完整复现（deviation score 版）
#
# 图的内容：254 个转录因子（TF）× 35 个细胞类型 的 motif 富集热图；
#           行按 36 个 TF 家族分块（左侧色带），列按主导家族分块，
#           底部圆点标注该细胞类型所属器官，底部数字为原图编号。
#
# 一站式脚本，串联原 4 个步骤：
#   Step 1  准备输入（降采样细胞 + ATAC counts + peak 坐标）        -> output/fig2b/cv_input.rds
#   Step 2  chromVAR 计算 deviation score（top-2,000 ACR 口径）      -> output/fig2b/cv_top2000_mat.rds
#   Step 3  按原图版式绘制热图                                        -> output/fig2b/Fig2b_deviation_reproduced.png/.pdf
#   Step 4  与原图并排对比（需 data/Fig2b_crop.png，缺失则跳过）      -> output/fig2b/Fig2b_deviation_vs_original.png
#
# 运行：cd <仓库根目录> && Rscript R/figures/Fig2b_TF_enrichment.R
# 可选环境变量（带默认值，均可省略）：
#   PER_CT   每细胞类型最多抽样的细胞数（默认 400）
#   TOP_N    top-N ACR 数量              （默认 2000）
#   MAXCT    >0 时只跑前 N 个细胞类型（调试用，默认 0 = 全部）
#   FLOOR_Q  噪声底：正值分位           （默认 0.70）
#   CAP_Q    饱和点：值分位             （默认 0.925）
#   MAXFAM   每个家族最多保留的 motif 数（默认 24）
#
# 依赖包：Seurat, Signac, Matrix, GenomicRanges, GenomeInfoDb, SummarizedExperiment,
#        chromVAR, motifmatchr, BSgenome.Osativa.MSU.MSU7, BiocParallel, ComplexHeatmap,
#        circlize, png
#
# 依赖文件（相对仓库根目录）：
#   data/GSE232863_scRNA_omics.Rds                        GEO 原始对象（5.5 GB）
#   data/official_code/dongwei-2023-.../rice.motif.pfmList2.Rds   PlantTFDB motif PFM
#   output/fig2b/primary_markers/*.rds                    各细胞类型 marker（前置步骤产出）
#   ref/fig2b_family_palette.rds                          36 家族色（可选，缺失用内置备用配色）
#   data/Fig2b_crop.png                                   原图裁剪（可选，仅 Step 4 用）
# =============================================================================

## ---- 载入公共配置（_paths.R 与本脚本同目录）-------------------------------
for (.p in c("_paths.R", "R/figures/_paths.R", "../_paths.R")) {
  if (file.exists(.p)) { source(.p); break }
}
if (!exists("PROJ", inherits = FALSE))
  stop("Cannot find R/figures/_paths.R. Run this script from the repository root.")

RES     <- out_dir("fig2b")                 # 中间结果与产物
DIR_MK  <- file.path(RES, "primary_markers")# 各细胞类型 marker（前置步骤产出）
PFM_FILE <- file.path(DIR_DATA, "official_code",
                      "dongwei-2023-Single_Cell_Multiomics_in_Rice-6020ed1",
                      "data", "rice.motif.pfmList2.Rds")

log <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
mem <- function() round(sum(gc()[, 2]) / 1024, 1)

## 参数
PER_CT  <- as.integer(Sys.getenv("PER_CT",  "400"))
TOP_N   <- as.integer(Sys.getenv("TOP_N",   "2000"))
MAXCT   <- as.integer(Sys.getenv("MAXCT",   "0"))
FLOOR_Q <- as.numeric(Sys.getenv("FLOOR_Q", "0.70"))
CAP_Q   <- as.numeric(Sys.getenv("CAP_Q",   "0.925"))
MAXFAM  <- as.integer(Sys.getenv("MAXFAM",  "24"))

# =============================================================================
# Step 1 —— 准备输入
# =============================================================================
log("========== Step 1: prepare input ==========")
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix); library(GenomicRanges)
})

need_file(RDS_SCRNA, "GSE232863 Seurat object",
          "download from GEO (see README section 1) and put it under data/")
obj <- readRDS(RDS_SCRNA)
log("loaded | dim =", nrow(obj), "x", ncol(obj), "| assays:",
    paste(Assays(obj), collapse = ","), "| mem(GB) =", mem())

## 只留 ATAC
DefaultAssay(obj) <- "ATAC"
for (a in setdiff(Assays(obj), "ATAC")) obj[[a]] <- NULL
gc(); log("after drop assays | mem(GB) =", mem())

ct <- as.character(obj$tissue_cluster)
log("cell types:", length(unique(ct)), "| cells:", length(ct))

## 每类最多 PER_CT 个细胞
set.seed(1)
sel <- unlist(lapply(split(colnames(obj), ct), function(x)
  if (length(x) <= PER_CT) x else sample(x, PER_CT)), use.names = FALSE)
log("down-sampled cells:", length(sel))

cnt <- GetAssayData(obj, assay = "ATAC", layer = "counts")
gr  <- granges(obj[["ATAC"]])
cnt <- as(cnt[, sel], "dgCMatrix")
meta <- data.frame(cell = sel, ct = ct[match(sel, colnames(obj))],
                   stringsAsFactors = FALSE)
rm(obj); gc()

saveRDS(list(counts = cnt, meta = meta, gr = gr), file.path(RES, "cv_input.rds"))
log("saved output/fig2b/cv_input.rds | size(MB) =",
    round(file.size(file.path(RES, "cv_input.rds")) / 1048576, 1))

# =============================================================================
# Step 2 —— chromVAR 计算 deviation score
# =============================================================================
log("========== Step 2: chromVAR deviation score ==========")
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(chromVAR); library(motifmatchr)
  library(GenomeInfoDb); library(BSgenome.Osativa.MSU.MSU7); library(BiocParallel)
})
register(SerialParam())   # 强制串行，避免并行复制 counts 致内存爆

inp <- readRDS(file.path(RES, "cv_input.rds"))
cnt <- inp$counts; meta <- inp$meta; gr <- inp$gr
rm(inp); gc()

std  <- paste0("chr", 1:12)
keep <- as.character(seqnames(gr)) %in% std
log("peaks on standard chromosomes:", sum(keep), "/", length(keep))
se <- SummarizedExperiment(assays = list(counts = cnt[keep, , drop = FALSE]),
                           rowRanges = gr[keep])
rm(cnt, gr); gc()

GENOME <- GenomeInfoDb::renameSeqlevels(BSgenome.Osativa.MSU.MSU7,
                                        setNames(paste0("chr", 1:12), paste0("Chr", 1:12)))
se <- addGCBias(se, genome = GENOME)

need_file(PFM_FILE, "PlantTFDB motif PFM list",
          "from the official code bundle (Zenodo 10.5281/zenodo.15501781), unzip into data/official_code/")
pfm <- readRDS(PFM_FILE)
log("PFM:", length(pfm), "motifs")

log("matchMotifs (all peaks) ...")
mm <- matchMotifs(pfm, se, genome = GENOME)
log("matchMotifs done | dim:",
    paste(dim(assay(mm, "motifMatches")), collapse = " x "))

## (B) 全局偏差分数（对照）
cts <- sort(unique(meta$ct))
cells_of <- split(meta$cell, meta$ct)
avg_ct <- function(d, cells) {
  m <- Matrix::rowMeans(d[, cells, drop = FALSE])
  m[!is.finite(m)] <- NA_real_
  m
}
GLOB_RDS <- file.path(RES, "cv_global_mat.rds")
if (!file.exists(GLOB_RDS)) {
  log("computeDeviations (global) ...")
  dv <- computeDeviations(object = se, annotations = mm)
  ds <- assay(dv, "z")
  matG <- sapply(cts, function(cc) avg_ct(ds, cells_of[[cc]]))
  colnames(matG) <- cts
  saveRDS(matG, GLOB_RDS)
  log("saved cv_global_mat.rds :", paste(dim(matG), collapse = " x "))
  rm(dv, ds); gc()
} else log("skip global (cv_global_mat.rds already exists)")

## (A) 每类 top 2,000 ACR 偏差分数
if (!dir.exists(DIR_MK))
  stop("\n[MISSING per-cell-type marker tables]\n  ", DIR_MK, "\n",
       "  -> produced by an upstream step (FindMarkers per cell type, saved as .rds\n",
       "     with an avg_log2FC column and peak coordinates as row names).\n",
       "     See R/figures/README.md section 4 for the code to rebuild them.\n",
       call. = FALSE)
mfs <- list.files(DIR_MK, pattern = "\\.rds$", full.names = TRUE)
names(mfs) <- sub("\\.rds$", "", basename(mfs))
norm_key <- function(x) gsub("[^A-Za-z0-9]", "_", x)
names(mfs) <- norm_key(names(mfs))
want <- if (MAXCT > 0) cts[seq_len(min(MAXCT, length(cts)))] else cts

matT <- matrix(NA_real_, nrow = length(pfm), ncol = length(cts),
               dimnames = list(names(pfm), cts))
mm_rows <- rownames(se)
log("computing top", TOP_N, "ACR deviation scores per cell type | n =", length(want))

for (k in seq_along(want)) {
  cc <- want[k]
  f <- mfs[[norm_key(cc)]]
  if (is.null(f) || !file.exists(f)) { log("  [skip]", cc, "- no marker table"); next }
  mk <- readRDS(f)
  ord <- order(mk$avg_log2FC, decreasing = TRUE)
  top <- head(rownames(mk)[ord], TOP_N)
  idx <- match(top, mm_rows); idx <- idx[!is.na(idx)]
  if (length(idx) < 200) { log("  [skip]", cc, "- too few peaks matched:", length(idx)); next }

  se_s <- se[idx, ]
  mm_s <- mm[idx, ]
  ok <- Matrix::colSums(assay(mm_s, "motifMatches")) > 0
  res <- rep(NA_real_, nrow(matT)); names(res) <- rownames(matT)

  r <- tryCatch({
    set.seed(1)
    d <- computeDeviations(object = se_s, annotations = mm_s[, ok, drop = FALSE])
    ds_s <- assay(d, "z")
    avg_ct(ds_s, cells_of[[cc]])
  }, error = function(e) { log("  [FAIL]", cc, ":", conditionMessage(e)); NULL })

  if (!is.null(r)) res[names(r)] <- r
  matT[, cc] <- res
  log(sprintf("  [%2d/%2d] %-30s peaks=%4d motifs=%3d mean=%+.3f max=%+.2f",
              k, length(want), cc, length(idx), sum(ok),
              mean(r, na.rm = TRUE), max(r, na.rm = TRUE)))
  if (k %% 10 == 0) { saveRDS(matT, file.path(RES, "cv_top2000_mat.rds")); gc() }
}
saveRDS(matT, file.path(RES, "cv_top2000_mat.rds"))
log("Step 2 done | dim:", paste(dim(matT), collapse = " x "), "| mem(GB) =", mem())

# =============================================================================
# Step 3 —— 按原图版式绘图
# =============================================================================
log("========== Step 3: draw the deviation-score heatmap ==========")
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid)
})

OUT <- "Fig2b_deviation_reproduced"

## 家族定义（顺序取自原图底部图例）
FAM <- c("AP-2","ARF","ARR-B","B3","BBR-BPC","BES1","bHLH","bZIP","C2H2","C3H",
         "CAMTA","CPP","Dof","E2F/DP","EIL","ERF","FAR1","G2-like","GATA",
         "HD-ZIP","HSF","LBD","LFY","M-type MADS","MIKC MADS","MYB","MYB related",
         "NAC","RAV","SBP","TALE","TCP","Trihelix","WOX","WRKY","YABBY")

## 36 个家族的颜色：优先用从原图采样得到的调色板；缺失则用内置备用配色
FAM_PAL <- file.path(DIR_REF, "fig2b_family_palette.rds")
if (file.exists(FAM_PAL)) {
  FAM_COL <- readRDS(FAM_PAL)
} else {
  warning("ref/fig2b_family_palette.rds not found; falling back to a built-in ",
          "36-colour palette.\n  Family blocks keep their order and position; ",
          "only the strip colours differ from the paper.", call. = FALSE)
  .h <- grDevices::hcl(seq(15, 375, length.out = 37)[1:36], c = 70, l = 62)
  FAM_COL <- t(grDevices::col2rgb(.h) / 255)
  rownames(FAM_COL) <- FAM
}
FAM_HEX <- apply(FAM_COL, 1, function(p) sprintf("#%02X%02X%02X",
                       round(p[1]*255), round(p[2]*255), round(p[3]*255)))
names(FAM_HEX) <- FAM

norm_fam <- function(x) {
  k <- gsub("[^a-z0-9]", "", tolower(FAM))
  i <- match(gsub("[^a-z0-9]", "", tolower(x)), k)
  ifelse(is.na(i), NA_character_, FAM[i])
}

M <- readRDS(file.path(RES, "cv_top2000_mat.rds"))
raw_fam <- sub(":.*", "", rownames(M))
fam <- norm_fam(raw_fam)
if (anyNA(fam)) {
  M <- M[!is.na(fam), , drop = FALSE]; fam <- fam[!is.na(fam)]
}

if (MAXFAM > 0) {
  v <- apply(M, 1, var, na.rm = TRUE)
  keep <- unlist(lapply(split(seq_along(fam), fam), function(ix)
    ix[order(-v[ix])][seq_len(min(MAXFAM, length(ix)))]), use.names = FALSE)
  M <- M[sort(keep), , drop = FALSE]; fam <- fam[sort(keep)]
}
log("matrix:", nrow(M), "motifs x", ncol(M), "cell types | families:", length(unique(fam)))

## 行排序：家族序
ro <- order(match(fam, FAM))
M <- M[ro, , drop = FALSE]; fam <- fam[ro]

## 列排序：主导家族 + 质心
Spre <- { v <- pmax(M, 0); cp <- as.numeric(quantile(v, 0.99, na.rm = TRUE)); pmin(v / cp * 100, 100) }
fams_present <- intersect(FAM, unique(fam))
SC <- sapply(seq_len(ncol(M)), function(j)
  vapply(fams_present, function(f) mean(Spre[fam == f, j], na.rm = TRUE), 0))
peak_fam <- fams_present[apply(SC, 2, which.max)]
co <- order(match(peak_fam, FAM), -apply(SC, 2, max))
M <- M[, co, drop = FALSE]; peak_fam <- peak_fam[co]

row_pos <- apply(M, 1, function(v) {
  v <- pmax(v, 0)
  if (sum(v, na.rm = TRUE) == 0) return(0.5)
  sum(v * seq_len(ncol(M)) / ncol(M), na.rm = TRUE) / sum(v, na.rm = TRUE)
})
ro2 <- order(match(fam, FAM), row_pos)
M <- M[ro2, , drop = FALSE]; fam <- fam[ro2]

## 数值映射（噪声底 + 饱和点）
V <- pmax(M, 0)
pos <- V[V > 0]
floor_v <- as.numeric(quantile(pos, FLOOR_Q, na.rm = TRUE))
V[V < floor_v] <- 0
cap <- as.numeric(quantile(V, CAP_Q, na.rm = TRUE)) / 0.45
V[is.na(V)] <- 0
V <- pmin(V / cap * 100, 100)
log(sprintf("noise floor=%.3f saturation=%.3f | empty %.1f%% mid-blue %.1f%% warm %.1f%%",
            floor_v, cap, 100 * mean(V < 5), 100 * mean(V >= 5 & V < 45),
            100 * mean(V >= 45)))

stair_idx <- function(X, w = 0.12) {
  nr <- nrow(X); nc <- ncol(X)
  band <- outer(seq_len(nr) / nr, seq_len(nc) / nc, function(r, c) abs(c - r) < w)
  mean(X[band], na.rm = TRUE) / mean(X[!band], na.rm = TRUE)
}
SI <- stair_idx(V)
log(sprintf("staircase index = %.2f | warm cells = %.1f%%", SI, 100 * mean(V > 50)))

## 注释
organ <- sub(":.*", "", colnames(M))
ORG <- c("Bud" = "#2E9B8F", "Flag" = "#7B4EA3", "leaf" = "#9A8B2E", "Root" = "#2C6FB5",
         "SAM" = "#43A047", "SP" = "#D9647E", "Seed" = "#8C6D3F", "ST" = "#C79A2E")
ORG <- ORG[names(ORG) %in% unique(organ)]
if (any(!organ %in% names(ORG)))
  ORG <- c(ORG, setNames(rep("#777777", sum(!unique(organ) %in% names(ORG))),
                         setdiff(unique(organ), names(ORG))))
canon <- sort(colnames(M))
num <- setNames(match(colnames(M), canon), colnames(M))

## 色标（逐像素采样自原图）
COL_FUN <- colorRamp2(
  c(0, 13.8, 29.5, 41.3, 48.5, 49.2, 55.1, 61, 65, 70.9, 76.8, 80.7, 87.8, 94.5, 100),
  c("#94C6D3","#85B9CD","#79A9C5","#6C90BB","#6892BD","#C0604E","#D26236","#D57038",
    "#D77938","#DC8539","#DF9239","#E19F3B","#EABA40","#F0D543","#F6E845"))

ht <- Heatmap(
  V, name = "Enrichment", col = COL_FUN,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = factor(fam, levels = FAM),
  column_split = factor(peak_fam, levels = intersect(FAM, unique(peak_fam))),
  row_title = NULL, column_title = NULL,
  show_row_names = FALSE, show_column_names = FALSE,
  row_gap = unit(0.14, "mm"), column_gap = unit(0.05, "mm"),
  width = unit(205, "mm"), height = unit(168, "mm"),
  border = FALSE,
  left_annotation = rowAnnotation(
    `TF family` = fam,
    col = list(`TF family` = FAM_HEX),
    show_annotation_name = FALSE, show_legend = FALSE,
    simple_anno_size = unit(3.2, "mm")),
  bottom_annotation = HeatmapAnnotation(
    Organ = anno_points(rep(0.5, ncol(V)), ylim = c(0, 1), pch = 16,
                        size = unit(1.7, "mm"), gp = gpar(col = ORG[organ]),
                        height = unit(3.6, "mm"), axis = FALSE),
    ` ` = anno_text(num[colnames(V)], rot = 90, just = c("center", "top"),
                    gp = gpar(fontsize = 6.4), height = unit(11, "mm")),
    show_annotation_name = FALSE, show_legend = FALSE),
  heatmap_legend_param = list(
    at = c(0, 50, 100), labels = c("0", "50", "100"),
    legend_height = unit(30, "mm"), grid_width = unit(3.2, "mm"),
    title = "Enrichment", title_position = "lefttop-rot",
    title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8),
    border = NA, ticks = TRUE))

W <- 10.2; H <- 8.90
for (ext in c("png", "pdf")) {
  fn <- file.path(RES, sprintf("%s.%s", OUT, ext))
  if (ext == "png") png(fn, width = W, height = H, units = "in", res = 400) else
                     pdf(fn, width = W, height = H)
  grid.newpage()
  pushViewport(viewport(x = 0.018, y = 0.118, width = 0.975, height = 0.845, just = c("left", "bottom")))
  draw(ht, newpage = FALSE, heatmap_legend_side = "right",
       padding = unit(c(1, 1, 0, 1), "mm"))
  popViewport()
  grid.text("b", x = unit(0.016, "npc"), y = unit(0.985, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))

  x0 <- 0.075; x1 <- 0.985; y0 <- 0.022; y1 <- 0.118
  pushViewport(viewport(x = x0, y = y0, width = x1 - x0, height = y1 - y0,
                        just = c("left", "bottom"), xscale = c(0, 36), yscale = c(0, 1)))
  for (i in seq_along(FAM)) {
    grid.rect(x = unit((i - 1), "native"), y = unit(1, "npc"),
              width = unit(1, "native"), height = unit(0.40, "npc"),
              just = c("left", "top"), gp = gpar(fill = FAM_HEX[i], col = "white", lwd = 0.5))
    grid.text(FAM[i], x = unit((i - 0.5), "native"), y = unit(0.53, "npc"),
              just = c("center", "top"), rot = 90, gp = gpar(fontsize = 6.4))
  }
  popViewport()
  grid.text("TF family", x = unit(x0 - 0.010, "npc"), y = unit((y0 + y1) / 2, "npc"),
            just = c("right", "center"), gp = gpar(fontsize = 8.5))
  dev.off()
  log("[ok]", fn)
}
saveRDS(list(mat = M, v = V, fam = fam, peak_fam = peak_fam, cap = cap,
             stair = SI, num = num), file.path(RES, sprintf("%s_data.rds", OUT)))

# =============================================================================
# Step 4 —— 与原图并排对比（需要原图裁剪；缺失则跳过）
# =============================================================================
log("========== Step 4: side-by-side comparison with the original ==========")
CROP <- file.path(DIR_DATA, "Fig2b_crop.png")
if (!file.exists(CROP)) {
  log("[skip] ", CROP, " not found - skipping the side-by-side panel (main figure unaffected)")
} else {
suppressPackageStartupMessages({ library(png) })

A <- readPNG(CROP)
hA <- dim(A)[1]; wA <- dim(A)[2]
A <- A[round(hA * 0.075):hA, round(wA * 0.09):round(wA * 0.975), , drop = FALSE]
aspA <- dim(A)[2] / dim(A)[1]

B <- readPNG(file.path(RES, sprintf("%s.png", OUT)))
aspB <- dim(B)[2] / dim(B)[1]

H <- 6.4
wA2 <- H * aspA; wB2 <- H * aspB
GAP <- 0.45; M0 <- 0.75
W2 <- M0 * 2 + wA2 + GAP + wB2
Ht <- H + 1.25

png(file.path(RES, "Fig2b_deviation_vs_original.png"), width = W2, height = Ht,
    units = "in", res = 300)
grid.newpage()
toX <- function(inch) inch / W2
toY <- function(inch) inch / Ht
yTop <- Ht - 0.55
xA <- M0; xB <- M0 + wA2 + GAP

for (k in 1:2) {
  img <- if (k == 1) A else B
  asp <- if (k == 1) aspA else aspB
  x0  <- if (k == 1) xA else xB
  h2  <- H; w2 <- h2 * asp
  pushViewport(viewport(x = unit(toX(x0), "npc"),
                        y = unit(toY(yTop - h2), "npc"),
                        width = unit(toX(w2), "npc"),
                        height = unit(toY(h2), "npc"),
                        just = c("left", "bottom")))
  grid.raster(img, interpolate = TRUE)
  grid.rect(gp = gpar(col = "grey65", fill = NA, lwd = 0.7))
  popViewport()
}
grid.text("a   Original  Fig. 2b (Nature 644, 724-729)",
          x = unit(toX(xA), "npc"), y = unit(toY(yTop + 0.06), "npc"),
          just = c("left", "bottom"), gp = gpar(fontsize = 11.5, fontface = "bold"))
grid.text("b   Our reproduction",
          x = unit(toX(xB), "npc"), y = unit(toY(yTop + 0.06), "npc"),
          just = c("left", "bottom"), gp = gpar(fontsize = 11.5, fontface = "bold"))
D <- readRDS(file.path(RES, sprintf("%s_data.rds", OUT)))
lab <- sprintf("%d TF x %d cell types;  %d TF families",
               nrow(D$mat), ncol(D$mat), length(unique(D$fam)))
grid.text(lab,
          x = unit(toX(xB), "npc"), y = unit(toY(yTop - H - 0.06), "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 8.5, col = "grey30"))
dev.off()
log("[ok] output/fig2b/Fig2b_deviation_vs_original.png")
}

log("========== all done ==========")
