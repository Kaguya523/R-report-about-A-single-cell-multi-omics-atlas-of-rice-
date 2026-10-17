# =============================================================================
#  _paths.R —— 三张核心图复现脚本的公共配置（路径 / 输入 / 小工具）
# -----------------------------------------------------------------------------
#  设计原则
#    1. 仓库里不出现任何绝对路径 —— 换电脑、换目录、换账号都不用改脚本
#    2. 大文件（Rds / h5）放 data/，不进仓库；从原图数字化的小参考表放 ref/，随仓库提交
#    3. 所有产物统一写到 output/<图号>/
#    4. 运行时的提示文字一律用英文 —— Windows 上 R 的 LC_CTYPE 常为 C，
#       stop()/warning() 里的中文会被重编码成乱码；中文说明放在 README 与注释里
#
#  用法（在每个脚本开头，二选一）
#    a) 在仓库根目录或 R/figures/ 下运行，脚本会自动找到本文件
#    b) 数据不在默认位置时，先指定根目录：
#         Sys.setenv(RICE_PROJ = "D:/path/to/水稻scMultiome复现方案")
# =============================================================================

## ---- 0. 让中文路径能正常打印 ----------------------------------------------
## Windows 上 R 默认 LC_CTYPE=C，会把 D:/文章/... 打印成 D:/<U+6587><U+7AE0>/...
## 只影响显示、不影响功能；切不过去就静默跳过。
for (.loc in c("Chinese", "C.UTF-8", "en_US.UTF-8")) {
  .ok <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", .loc),
                                   error = function(e) ""))
  if (nzchar(.ok) && !identical(.ok, "C")) break
}
rm(.loc, .ok)

## ---- 1. 项目根目录 ---------------------------------------------------------
## 优先读环境变量 RICE_PROJ；否则从当前工作目录逐级向上找「同时含 README.md 与 R/ 的目录」
.find_proj_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "README.md")) && dir.exists(file.path(d, "R"))) {
      return(d)
    }
    up <- dirname(d)
    if (identical(up, d)) {
      stop("Cannot locate the project root (a directory containing both ",
           "README.md and R/).\n",
           "  * Run this script from inside the repository, or\n",
           "  * set the environment variable first:\n",
           "      Sys.setenv(RICE_PROJ = \"<absolute path to repo>\")",
           call. = FALSE)
    }
    d <- up
  }
}

PROJ <- Sys.getenv("RICE_PROJ", unset = "")
PROJ <- if (nzchar(PROJ)) {
  normalizePath(PROJ, winslash = "/", mustWork = TRUE)
} else {
  .find_proj_root()
}

## ---- 2. 三个固定目录 -------------------------------------------------------
DIR_DATA <- file.path(PROJ, "data")    # 原始数据（体积大，.gitignore 已排除）
DIR_REF  <- file.path(PROJ, "ref")     # 从原图数字化的小参考表（随仓库提交）
DIR_OUT  <- file.path(PROJ, "output")  # 图与表的产物

## 取 output/<name>/，不存在则创建
out_dir <- function(name) {
  d <- file.path(DIR_OUT, name)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

## ---- 3. 输入数据（可用环境变量覆盖，便于把大文件放在移动硬盘等位置）--------
RDS_SCRNA <- Sys.getenv("RICE_RDS_SCRNA",
                        unset = file.path(DIR_DATA, "GSE232863_scRNA_omics.Rds"))

## ---- 4. 小工具 -------------------------------------------------------------
## 带时间戳的日志（三张图要跑几分钟到几十分钟，便于看进度）
msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

## 缺文件时给一句「照着做就行」的提示，而不是让 R 抛原始报错
need_file <- function(path, what = "input file", hint = NULL) {
  if (!file.exists(path)) {
    stop(sprintf("\n[MISSING %s]\n  %s\n%s",
                 what, path,
                 if (is.null(hint)) "" else paste0("  -> ", hint, "\n")),
         call. = FALSE)
  }
  invisible(TRUE)
}

## 0–255 整数 RGB → "#RRGGBB"
hex <- function(r, g, b) grDevices::rgb(r, g, b, maxColorValue = 255)

msg("project root: ", PROJ)
