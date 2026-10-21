#!/usr/bin/env Rscript
# =============================================================================
#  00_selftest.R —— 四张图的"环境 + 素材"自检（不下载任何数据，10 秒内跑完）
#
#  作用：在正式跑图之前，一次性回答"我现在能不能跑 Fig 1b / 1c / 2a / 2b"。
#        检查四件事：
#          1. 项目根目录能否定位（决定所有相对路径是否成立）
#          2. ref/ 下的参考表是否齐全、格式是否正确
#          3. data/ 下的原始数据是否就位
#          4. 每张图需要的 R 包是否装好
#
#  运行：cd <仓库根目录> && Rscript R/figures/00_selftest.R
#  输出：一张就绪状态表；退出码 0 = 至少 Fig 1b 可以跑
# =============================================================================

for (.p in c("_paths.R", "R/figures/_paths.R", "../_paths.R")) {
  if (file.exists(.p)) { source(.p); break }
}
if (!exists("PROJ", inherits = FALSE))
  stop("Cannot find R/figures/_paths.R. Run this script from the repository root.")

ok   <- function(x) if (isTRUE(x)) "OK" else "NO"
mark <- function(x) if (isTRUE(x)) "[ok]   " else "[MISS] "

cat("\n=============================================================\n")
cat("  self-test for Fig. 1b / 1c / 2a / 2b\n")
cat("=============================================================\n")
cat("project root : ", PROJ, "\n\n", sep = "")

## ---- 1. 参考表 ref/ --------------------------------------------------------
cat("--- ref/ (digitised reference tables) ---\n")

pal_f <- file.path(DIR_REF, "fig2a_palette.csv")
dot_f <- file.path(DIR_REF, "fig2a_dots.csv")
pal_ok <- dot_ok <- FALSE
if (file.exists(pal_f)) {
  p <- tryCatch(read.csv(pal_f), error = function(e) NULL)
  pal_ok <- !is.null(p) && all(c("value", "r", "g", "b") %in% names(p)) && nrow(p) == 100
  cat(mark(pal_ok), "fig2a_palette.csv  ",
      if (pal_ok) sprintf("100 steps, %s -> %s", hex(p$r[1], p$g[1], p$b[1]),
                          hex(p$r[100], p$g[100], p$b[100])) else "expected 100 rows of value,r,g,b",
      "\n", sep = "")
} else cat(mark(FALSE), "fig2a_palette.csv  not found\n")

if (file.exists(dot_f)) {
  d <- tryCatch(read.csv(dot_f), error = function(e) NULL)
  dot_ok <- !is.null(d) && all(c("r", "g", "b") %in% names(d)) && nrow(d) == 56
  cat(mark(dot_ok), "fig2a_dots.csv     ",
      if (dot_ok) "56 dots" else "expected 56 rows of r,g,b", "\n", sep = "")
} else cat(mark(FALSE), "fig2a_dots.csv     not found\n")

cat(mark(TRUE), "Fig 1b / 1c need no ref/ tables (metadata + coordinates only)\n", sep = "")

## ---- 2. 原始数据 data/ -----------------------------------------------------
cat("\n--- data/ (raw inputs, not tracked by git) ---\n")
rds_ok <- file.exists(RDS_SCRNA)
cat(mark(rds_ok), basename(RDS_SCRNA), " ",
    if (rds_ok) sprintf("%.2f GB", file.size(RDS_SCRNA) / 1024^3) else
      "download from GEO GSE232863 (see README section 2)", "\n", sep = "")

pfm_ok <- file.exists(PFM_FILE <- file.path(DIR_DATA, "official_code",
  "dongwei-2023-Single_Cell_Multiomics_in_Rice-6020ed1", "data", "rice.motif.pfmList2.Rds"))
cat(mark(pfm_ok), "rice.motif.pfmList2.Rds  ",
    if (pfm_ok) "present" else "unzip official code bundle from Zenodo 10.5281/zenodo.15501781", "\n", sep = "")

mk_ok <- dir.exists(file.path(DIR_OUT, "fig2b", "primary_markers"))
cat(mark(mk_ok), "output/fig2b/primary_markers/  ",
    if (mk_ok) "present" else "produced by an upstream FindMarkers step (see R/figures/README.md)", "\n", sep = "")

## ---- 3. R 包 ---------------------------------------------------------------
cat("\n--- R packages ---\n")
pk <- function(p) requireNamespace(p, quietly = TRUE)
P1B <- c("Seurat", "ggplot2", "patchwork", "dplyr")
P1C <- c("Seurat", "Signac", "Matrix", "ggplot2", "ggpattern", "dplyr")
P2A <- c("Seurat", "Signac", "Matrix", "ComplexHeatmap", "circlize")
P2B <- c("chromVAR", "motifmatchr", "SummarizedExperiment", "GenomeInfoDb",
         "BSgenome.Osativa.MSU.MSU7", "BiocParallel", "png")
have1b <- vapply(P1B, pk, logical(1))
have1c <- vapply(P1C, pk, logical(1))
have2a <- vapply(P2A, pk, logical(1))
have2b <- vapply(P2B, pk, logical(1))

cat(mark(all(have1b)), "Fig 1b : ", sum(have1b), "/", length(P1B), " packages",
    if (any(!have1b)) paste0("   missing: ", paste(P1B[!have1b], collapse = ", ")) else "",
    "\n", sep = "")
cat(mark(all(have1c)), "Fig 1c : ", sum(have1c), "/", length(P1C), " packages",
    if (any(!have1c)) paste0("   missing: ", paste(P1C[!have1c], collapse = ", ")) else "",
    "\n", sep = "")
cat(mark(all(have2a)), "Fig 2a : ", sum(have2a), "/", length(P2A), " packages",
    if (any(!have2a)) paste0("   missing: ", paste(P2A[!have2a], collapse = ", ")) else "",
    "\n", sep = "")
cat(mark(all(have2b)), "Fig 2b : ", sum(have2b), "/", length(P2B), " packages",
    if (any(!have2b)) paste0("   missing: ", paste(P2B[!have2b], collapse = ", ")) else "",
    "\n", sep = "")

## ---- 4. 就绪状态汇总 -------------------------------------------------------
cat("\n=============================================================\n")
cat("  readiness\n")
cat("=============================================================\n")
cat(sprintf("  Fig 1b : %s  data=%s  pkgs=%s\n",
            ok(rds_ok && all(have1b)), ok(rds_ok), ok(all(have1b))))
cat(sprintf("  Fig 1c : %s  data=%s  pkgs=%s\n",
            ok(rds_ok && all(have1c)), ok(rds_ok), ok(all(have1c))))
cat(sprintf("  Fig 2a : %s  ref=%s  data=%s  pkgs=%s\n",
            ok(pal_ok && dot_ok && rds_ok && all(have2a)),
            ok(pal_ok && dot_ok), ok(rds_ok), ok(all(have2a))))
cat(sprintf("  Fig 2b : %s  data=%s  markers=%s  pkgs=%s\n",
            ok(rds_ok && pfm_ok && mk_ok && all(have2a) && all(have2b)),
            ok(rds_ok && pfm_ok), ok(mk_ok), ok(all(have2a) && all(have2b))))
cat("-------------------------------------------------------------\n")
if (rds_ok && all(have1b)) {
  cat("  next: Rscript R/figures/Fig1b_reproduce.R   (UMAP 重算，约 6 min)\n")
} else {
  cat("  next: put GSE232863_scRNA_omics.Rds into data/ (see README section 2)\n")
}
cat("=============================================================\n\n")

invisible(if (rds_ok) 0L else 1L)
