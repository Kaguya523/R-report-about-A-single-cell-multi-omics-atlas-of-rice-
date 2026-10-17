# =============================================================================
# _common.R —— 全项目公共函数与常量
# -----------------------------------------------------------------------------
# ⚠️ 调用方式（顺序不能反！）：
#     source("R/_bootstrap.R")   # 先切 UTF-8 locale
#     source("R/_common.R")      # 再读本文件
#   原因：R 的 source() 会先把整个文件解析完再执行，所以 locale 必须在
#         读取本文件【之前】就设好，否则本文件里的中文注释与字符串会乱码。
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ---- 路径 -------------------------------------------------------------------
# ⚠️ 坑：项目目录名含中文时，normalizePath(".") 在 C locale 下会把路径编码搞乱
#    （报错形如 Cannot find directory 'D:/fg+ /f04g(...'）。
#    因此这里一律使用相对路径，并要求所有脚本都从项目根目录运行：
#        cd 水稻scMultiome复现方案 && Rscript R/xx.R
DIR_FIG <- "output"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

# ---- 论文 Fig.1b 的官方数字（复现成功判据，硬编码用于校验）------------------
# 来源：Nature 644, 722-729 (2025) Fig.1b，八项相加 = 116,564
PAPER_NUCLEI <- c(
  "Crown root" = 31108, "Culm" = 19740, "Leaf" = 13771, "Flag leaf" = 5563,
  "Shoot apex" = 13791, "Tiller bud" = 8858, "Panicle" = 17457, "Seed" = 6276
)
PAPER_TOTAL <- 116564L

# ---- 配色：与论文 Fig.1b 器官配色近似 ---------------------------------------
ORGAN_COLORS <- c(
  "Crown root" = "#7FC97F", "Culm" = "#F781BF", "Leaf" = "#8DB600",
  "Flag leaf"  = "#C8A415", "Shoot apex" = "#6FC7D0", "Tiller bud" = "#E8796B",
  "Panicle"    = "#B384F0", "Seed" = "#4C9BF0"
)

# ---- 工具函数 ---------------------------------------------------------------

#' 从候选列名里挑出实际存在的第一个（列名容错，避免硬编码崩溃）
#' @param df data.frame
#' @param candidates 候选列名向量
#' @param what 用途描述，用于报错信息
pick <- function(df, candidates, what = "column") {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop(sprintf("[pick] 找不到 %s。候选: %s\n实际列名: %s",
                 what, paste(candidates, collapse = ", "),
                 paste(colnames(df), collapse = ", ")))
  }
  hit[1]
}

#' 读表（自动兼容 csv / tsv / rds / xlsx）
read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    return(readRDS(path))
  }
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("需要 readxl 包")
    return(as.data.frame(readxl::read_excel(path)))
  }
  if (ext %in% c("gz")) {
    # .csv.gz / .tsv.gz
    inner <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
    if (inner == "tsv") {
      return(read.delim(gzfile(path), check.names = FALSE))
    }
    return(read.csv(gzfile(path), check.names = FALSE))
  }
  if (ext == "tsv") return(read.delim(path, check.names = FALSE))
  read.csv(path, check.names = FALSE)
}

#' 统一主题（深色/浅色通用，出图给报告用浅底更清楚）
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
      axis.text  = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      legend.key.size = unit(0.4, "cm")
    )
}

#' 保存图片（同时出 pdf + png，便于报告与PPT各取所需）
#' 坑：Windows 上 ggsave 默认可能走 ragg/agg 设备并失败（"agg could not write"），
#'     这里显式指定 device，并对 png 做多级回退。
save_fig <- function(plot, name, width = 6, height = 5, dpi = 300) {
  pdf_path <- file.path(DIR_FIG, paste0(name, ".pdf"))
  png_path <- file.path(DIR_FIG, paste0(name, ".png"))

  ok_pdf <- tryCatch({
    ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
    TRUE
  }, error = function(e) {
    tryCatch({ ggsave(pdf_path, plot, width = width, height = height); TRUE },
             error = function(e2) { warning("PDF 保存失败: ", conditionMessage(e2)); FALSE })
  })

  # PNG：本机 ggsave 默认走 ragg，ragg 在此环境写盘失败（"agg could not write"），
  #      且只 warning 不 error，因此不走 ggsave，直接用设备。
  #      优先 cairo（能正确渲染中文），失败再退回 windows 设备（中文可能乱码）。
  if (file.exists(png_path)) unlink(png_path)
  try({
    grDevices::png(png_path, width = width * dpi, height = height * dpi, res = dpi,
                   bg = "white", type = "cairo")
    print(plot); grDevices::dev.off()
  }, silent = TRUE)
  if (!file.exists(png_path) || file.size(png_path) < 1000) {
    try({
      grDevices::png(png_path, width = width * dpi, height = height * dpi, res = dpi,
                     bg = "white", type = "windows")
      print(plot); grDevices::dev.off()
    }, silent = TRUE)
  }
  ok_png <- file.exists(png_path) && file.size(png_path) > 1000

  message("[save_fig] ", if (ok_pdf) pdf_path else "(pdf 失败)",
          " | ", if (ok_png) png_path else "(png 失败)")
  invisible(plot)
}

#' 校验细胞总数是否与论文一致
check_total <- function(n, tol = 0) {
  ok <- abs(n - PAPER_TOTAL) <= tol
  msg <- sprintf("细胞总数 = %s（论文 %s）%s", format(n, big.mark = ","),
                 format(PAPER_TOTAL, big.mark = ","),
                 if (ok) "  ✅ 一致" else "  ⚠️ 不一致，需排查过滤/去双细胞参数")
  message(msg)
  invisible(ok)
}

#' 器官名归一化：官方对象里命名不统一（Root/Crown Root、Flag/Flag leaf、ST/Culm...）
norm_organ <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("_", " ", x, fixed = TRUE)
  dplyr::recode(x,
    "Root" = "Crown root", "Crown Root" = "Crown root", "Crown root" = "Crown root",
    "ST" = "Culm", "Culm" = "Culm", "clum" = "Culm",
    "leaf" = "Leaf", "Leaf" = "Leaf", "Third Leaf" = "Leaf",
    "Flag" = "Flag leaf", "Flag Leaf" = "Flag leaf", "Flag leaf" = "Flag leaf",
    "SAM" = "Shoot apex", "Shoot apex" = "Shoot apex",
    "Bud" = "Tiller bud", "Tiller Bud" = "Tiller bud", "Tiller bud" = "Tiller bud",
    "SP" = "Panicle", "Panicle" = "Panicle", "Spikelet" = "Panicle",
    "Young Spike" = "Panicle", "Seed" = "Seed",
    .default = x
  )
}

#' 从 Seurat 对象里智能定位元数据列（不同版本列名不同）
guess_meta_col <- function(meta, kind = c("organ", "celltype", "cluster", "sample")) {
  kind <- match.arg(kind)
  cand <- switch(kind,
    organ    = c("tissue", "organ", "orig.ident", "sample", "Tissue", "Organ"),
    celltype = c("cell_type", "celltype", "cluster_names", "name_id", "annotation", "CellType"),
    cluster  = c("seurat_clusters", "final_id", "cluster", "RNA_snn_res.0.5"),
    sample   = c("orig.ident", "sample", "Sample", "batch")
  )
  hit <- intersect(cand, colnames(meta))
  if (length(hit) == 0) {
    warning(sprintf("找不到 %s 列，候选: %s；实际: %s", kind,
                    paste(cand, collapse = ", "), paste(colnames(meta), collapse = ", ")))
    return(NULL)
  }
  hit[1]
}
