# =============================================================================
# 01_inspect_object.R —— 【第一步必跑】探查 Seurat 对象结构
# -----------------------------------------------------------------------------
# 为什么要先跑这一步？
#   官方脚本里全是硬编码列名（tissue / cluster_names / final_id ...），
#   不同 Seurat 版本、不同提交批次的元数据列名可能不一样。
#   先打印出真实结构，再改后续脚本，能省掉大量报错时间。
# 用法：Rscript R/01_inspect_object.R
# =============================================================================

source("R/_bootstrap.R")   # 先设 UTF-8 locale，否则中文乱码
source("R/_common.R")
suppressPackageStartupMessages(library(Seurat))

# ---- 配置：按实际下载位置修改 -----------------------------------------------
RDS_MAIN <- "data/GSE232863_scRNA_omics.Rds"                       # 5.1 GB，8 器官主图谱
RDS_LTP  <- "data/GSE273875_LTP_NIP_label_240622.rds.gz"           # 390 MB，osltpl120 突变体

inspect_seurat <- function(path, tag) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("### ", tag, " —— ", path, "\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")

  if (!file.exists(path)) {
    cat("⚠️ 文件不存在，跳过。请先运行 scripts/00_download.sh\n")
    return(invisible(NULL))
  }
  cat("文件大小: ", round(file.size(path) / 1024^3, 2), " GB\n")

  obj <- if (grepl("\\.gz$", path)) readRDS(gzfile(path)) else readRDS(path)

  cat("\n--- 对象类别 ---\n");        print(class(obj))
  cat("\n--- 维度 (基因/peak x 细胞) ---\n"); print(dim(obj))

  cat("\n--- Assays ---\n");          print(names(obj@assays))
  for (a in names(obj@assays)) {
    cat("  assay [", a, "] layers: ", paste(Layers(obj[[a]]), collapse = ", "), "\n", sep = "")
  }

  cat("\n--- Reductions (降维结果) ---\n"); print(names(obj@reductions))
  for (rd in names(obj@reductions)) {
    cat("  [", rd, "] dims = ", paste(dim(Embeddings(obj, rd)), collapse = " x "), "\n", sep = "")
  }

  cat("\n--- meta.data 列名与类型 ---\n")
  meta <- obj@meta.data
  print(data.frame(column = colnames(meta), class = sapply(meta, class), row.names = NULL))

  cat("\n--- 关键列的唯一取值（前 20 个）---\n")
  for (col in colnames(meta)) {
    v <- meta[[col]]
    if (is.character(v) || is.factor(v)) {
      u <- unique(as.character(v))
      cat(sprintf("  %-22s : n=%d | %s\n", col, length(u),
                  paste(head(u, 20), collapse = ", ")))
    }
  }

  cat("\n--- Idents 当前层级 ---\n")
  print(table(Idents(obj))[1:min(30, length(table(Idents(obj))))])

  # 自动识别各列
  cat("\n--- 自动识别的元数据列 ---\n")
  for (k in c("organ", "celltype", "cluster", "sample")) {
    cat(sprintf("  %-10s -> %s\n", k, guess_meta_col(meta, k)))
  }

  # 细胞总数校验（仅主图谱对象适用）
  if (grepl("GSE232863", path)) {
    cat("\n--- 细胞总数校验 ---\n")
    check_total(ncol(obj))
  }

  # 保存一份元数据摘要，便于写报告
  out <- file.path(DIR_FIG, paste0("meta_summary_", tag, ".csv"))
  write.csv(data.frame(column = colnames(meta), class = sapply(meta, class)),
            out, row.names = FALSE)
  cat("\n元数据结构已保存: ", out, "\n", sep = "")

  invisible(obj)
}

obj_main <- inspect_seurat(RDS_MAIN, "main_8organs")
obj_ltp  <- inspect_seurat(RDS_LTP,  "LTPL120_mutant")

cat("\n\n===== 下一步 =====\n")
cat("1. 把上面打印出的真实列名，填到 R/_common.R 的 guess_meta_col() 候选里；\n")
cat("2. 再运行 R/02_fig1b_umap_organ.R\n")
