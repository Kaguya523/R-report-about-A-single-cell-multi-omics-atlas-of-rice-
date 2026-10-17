# =============================================================================
# 00_setup.R —— 安装依赖并记录版本（复现性要求：版本号必须写进报告）
# 用法：Rscript R/00_setup.R
# =============================================================================

source("R/_bootstrap.R")   # 先设 UTF-8 locale，否则中文乱码

# ---- 1. CRAN 包 --------------------------------------------------------------
cran_pkgs <- c(
  # 单细胞核心
  "Seurat", "Signac", "harmony", "DoubletFinder", "clustree",
  # 可视化
  "ggplot2", "ggrepel", "ggrastr", "ggpubr", "patchwork", "cowplot",
  "pheatmap", "corrplot", "ComplexHeatmap", "RColorBrewer", "scales", "ggsci",
  # 数据处理
  "dplyr", "tidyr", "stringr", "readr", "data.table", "reshape2", "Matrix",
  # 富集 / 注释
  "clusterProfiler", "ChIPseeker", "GenomicRanges", "rtracklayer",
  # 拟时序 / 网络
  "monocle3", "hdWGCNA", "WGCNA", "igraph",
  # 基因组注释（MSU7）
  "EnsDb.Osativa.v55"
)

# Bioconductor 专属包（走 BiocManager）
bioc_pkgs <- c(
  "EnsDb.Osativa.v55", "ChIPseeker", "clusterProfiler",
  "GenomicRanges", "rtracklayer", "ComplexHeatmap",
  "monocle3", "hdWGCNA", "WGCNA"
)

install_if_missing <- function(pkgs, repos = "https://cloud.r-project.org") {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing) == 0) {
    message("全部已安装 ✅")
    return(invisible(NULL))
  }
  message("待安装: ", paste(missing, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = repos)
  }
  for (p in missing) {
    message(">>> installing ", p)
    tryCatch(
      if (p %in% bioc_pkgs) BiocManager::install(p, update = FALSE, ask = FALSE)
      else install.packages(p, repos = repos),
      error = function(e) warning("安装失败: ", p, " —— ", conditionMessage(e))
    )
  }
}

message("\n==== 开始检查/安装 CRAN 包 ====")
install_if_missing(cran_pkgs)
message("\n==== 开始检查/安装 Bioconductor 包 ====")
install_if_missing(bioc_pkgs)

# ---- 2. Python 侧（pySCENIC / CellOracle，路线进阶才需要）--------------------
# 建议用独立 conda 环境，不要污染系统 Python：
#   conda create -n scenic python=3.10 -y
#   conda activate scenic
#   pip install pyscenic celloracle
# 并把 cisTarget 数据库放在 data/ 下（见官方脚本 08）

# ---- 3. 记录会话版本（复现报告必交材料）-------------------------------------
out_dir <- "output"; dir.create(out_dir, showWarnings = FALSE)
ver_file <- file.path(out_dir, "session_versions.csv")

want <- c("R", "Seurat", "Signac", "harmony", "DoubletFinder",
          "monocle3", "hdWGCNA", "ChIPseeker", "clusterProfiler",
          "ggplot2", "pheatmap", "ComplexHeatmap", "EnsDb.Osativa.v55")
rows <- lapply(want, function(p) {
  v <- if (p == "R") paste(R.version$major, R.version$minor, sep = ".")
       else tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  data.frame(package = p, version = v, stringsAsFactors = FALSE)
})
ver_df <- do.call(rbind, rows)
write.csv(ver_df, ver_file, row.names = FALSE)
print(ver_df)

message("\n版本信息已写入: ", ver_file)
message("⚠️ 论文原环境为 R 4.2.2；若你的 R 更新（如 4.6.x），")
message("   请在报告的『复现结果与原论文的对比』一节中说明版本差异。")
