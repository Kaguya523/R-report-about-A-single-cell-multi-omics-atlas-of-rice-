#!/usr/bin/env Rscript
# =============================================================================
#  00b_make_primary_markers.R —— Fig 2b 的前置步骤
#  ===========================================================================
#  作用：为每个细胞类型在 ATAC assay 上跑一次 FindMarkers，产出
#        output/fig2b/primary_markers/<细胞类型>.rds
#        Fig2b_TF_enrichment.R 的 Step 2 会读这些表，按 avg_log2FC 取 top N 个 ACR。
#
#  为什么单独一个脚本：
#        Fig 2b 的输入不是"原始对象"，而是"每类细胞的差异 ACR 表"。
#        这一步与画图无关、耗时又长（每类一次逻辑回归），所以拆出来单独跑、
#        单独缓存。跑过一次之后，重画 Fig 2b 不需要再跑。
#
#  官方对应代码（07.Analyzing_the_TF_and_motif_by_Signac.R 第 31-52 行）：
#        rice.sub.atac.markers <- parallel::mclapply(unique(rice.sub$tissue_cluster),
#          function(x) FindMarkers(rice.sub, ident.1 = x, only.pos = T,
#                                  test.use = "LR", logfc.threshold = 0.05,
#                                  max.cells.per.ident = 300,
#                                  latent.vars = "nCount_ATAC"), mc.cores = 10)
#        本脚本逐行对齐，只把 mclapply 换成跨平台实现（Windows 上 mclapply 只跑单核）。
#
#  运行：cd <仓库根目录> && Rscript R/figures/00b_make_primary_markers.R
#
#  可选环境变量：
#    CELLS_PER_IDENT  每类最多抽多少细胞（默认 300，同官方）
#    MAXCT            >0 时只跑前 N 个细胞类型（调试用，默认 0 = 全部）
#    FORCE            1 = 已存在的 .rds 也重跑（默认 0，跳过已完成的）
#
#  输出：output/fig2b/primary_markers/*.rds     每类一张表（行名 = peak 坐标）
#        output/fig2b/primary_markers_log.csv   每类的 peak 数、耗时、状态
#
#  依赖：Seurat, Signac, Matrix
# =============================================================================

for (.p in c("_paths.R", "R/figures/_paths.R", "../_paths.R")) {
  if (file.exists(.p)) { source(.p); break }
}
if (!exists("PROJ", inherits = FALSE))
  stop("Cannot find R/figures/_paths.R. Run this script from the repository root.")

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix)
})

RES    <- out_dir("fig2b")
DIR_MK <- file.path(RES, "primary_markers")
dir.create(DIR_MK, recursive = TRUE, showWarnings = FALSE)

## ---- 参数（与官方 07 号脚本一致）-------------------------------------------
CELLS_PER_IDENT <- as.integer(Sys.getenv("CELLS_PER_IDENT", "300"))
MAXCT <- as.integer(Sys.getenv("MAXCT", "0"))
FORCE <- Sys.getenv("FORCE", "0") == "1"
MIN_LFC <- 0.05          # 官方 logfc.threshold
TEST_USE <- "LR"         # 官方 test.use
LATENT <- "nCount_ATAC"  # 官方 latent.vars

log <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
mem <- function() round(sum(gc()[, 2]) / 1024, 1)

## 文件名规范化：必须与 Fig2b_TF_enrichment.R 里的 norm_key() 完全一致
norm_key <- function(x) gsub("[^A-Za-z0-9]", "_", x)

need_file(RDS_SCRNA, "GSE232863 Seurat object",
          "download from GEO (see README section 2) and put it under data/")

log("========== loading object ==========")
obj <- readRDS(RDS_SCRNA)
log("loaded | dim =", nrow(obj), "x", ncol(obj),
    "| assays:", paste(Assays(obj), collapse = ","), "| mem(GB) =", mem())

if (!"ATAC" %in% Assays(obj))
  stop("The object has no 'ATAC' assay; this script needs the peak matrix.")

DefaultAssay(obj) <- "ATAC"
## 只保留 ATAC：RNA / atacRNA / chromvar 与这一步无关，删掉可省一半以上内存
for (a in setdiff(Assays(obj), "ATAC")) obj[[a]] <- NULL
invisible(gc(verbose = FALSE))
log("after dropping other assays | mem(GB) =", mem())

if (!"tissue_cluster" %in% colnames(obj@meta.data))
  stop("metadata column 'tissue_cluster' not found.")
Idents(obj) <- as.character(obj$tissue_cluster)

cts <- sort(unique(as.character(obj$tissue_cluster)))
cts <- cts[!is.na(cts) & nzchar(cts)]
if (MAXCT > 0) cts <- cts[seq_len(min(MAXCT, length(cts)))]
log("cell types:", length(cts), "| cells:", ncol(obj))

## ---- 逐类跑 FindMarkers ----------------------------------------------------
## Windows 上 mclapply 只跑单核，所以这里用统一接口：能并行就并行，不能就顺序
ncore <- if (.Platform$OS.type == "windows") 1L else
  max(1L, min(8L, parallel::detectCores() - 2L))
log("cores:", ncore)

run_one <- function(ct) {
  fn <- file.path(DIR_MK, paste0(norm_key(ct), ".rds"))
  if (file.exists(fn) && !FORCE) {
    return(data.frame(cell_type = ct, n_cells = NA_integer_, n_peaks = NA_integer_,
                      sec = NA_real_, status = "cached"))
  }
  t0 <- Sys.time()
  n_cells <- sum(as.character(obj$tissue_cluster) == ct, na.rm = TRUE)

  mk <- tryCatch(
    FindMarkers(obj, ident.1 = ct, only.pos = TRUE,
                test.use = TEST_USE,
                logfc.threshold = MIN_LFC,
                max.cells.per.ident = CELLS_PER_IDENT,
                latent.vars = LATENT,
                verbose = FALSE),
    error = function(e) { log("  [error]", ct, "->", conditionMessage(e)); NULL })

  sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  if (is.null(mk) || nrow(mk) == 0) {
    log("  [skip ]", ct, "| no peak passed the threshold |", sec, "s")
    return(data.frame(cell_type = ct, n_cells = n_cells, n_peaks = 0L,
                      sec = sec, status = "empty"))
  }

  ## FindMarkers 在 ATAC assay 上的行名就是 peak 坐标（chr:start-end），
  ## 这正是 Fig2b_TF_enrichment.R 用 rownames() 去 match 的东西。
  mk <- as.data.frame(mk)
  if (!"avg_log2FC" %in% colnames(mk))
    stop("FindMarkers output has no 'avg_log2FC' column; cannot continue.")

  saveRDS(mk, fn)
  log("  [ok   ]", ct, "| peaks =", nrow(mk), "| cells =", n_cells, "|", sec, "s")
  data.frame(cell_type = ct, n_cells = n_cells, n_peaks = nrow(mk),
             sec = sec, status = "written")
}

res <- if (ncore > 1L) {
  parallel::mclapply(cts, run_one, mc.cores = ncore)
} else {
  lapply(cts, run_one)
}
res <- do.call(rbind, res)

## ---- 汇总 ------------------------------------------------------------------
logf <- file.path(DIR_MK, "..", "primary_markers_log.csv")
write.csv(res, logf, row.names = FALSE)

n_ok <- sum(res$status %in% c("written", "cached"))
log("========== summary ==========")
log("written/cached:", n_ok, "/", nrow(res),
    "| empty:", sum(res$status == "empty"),
    "| total peaks:", sum(res$n_peaks, na.rm = TRUE))
log("tables ->", DIR_MK)
log("log    ->", normalizePath(logf, mustWork = FALSE))

if (n_ok == 0)
  stop("No marker table was produced. Check the ATAC assay and tissue_cluster column.")
cat("\nNext: Rscript R/figures/Fig2b_TF_enrichment.R\n")
