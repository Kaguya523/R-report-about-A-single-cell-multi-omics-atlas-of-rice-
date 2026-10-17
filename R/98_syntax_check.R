# =============================================================================
# 98_syntax_check.R —— 对所有 R 脚本做语法检查（不执行，只解析）
# 用法：Rscript R/98_syntax_check.R
# 注意：不要用 Rscript -e '...' 在 Git Bash 下跑多行脚本，可能 segfault；
#       写成文件再执行更稳。
# =============================================================================

source("R/_bootstrap.R")   # 先设 UTF-8 locale，否则中文乱码

files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
files <- files[!grepl("98_syntax_check", files)]

ok <- TRUE
for (f in files) {
  res <- tryCatch({
    parse(f)
    "OK"
  }, error = function(e) conditionMessage(e))
  if (identical(res, "OK")) {
    cat(sprintf("  [PASS] %s\n", f))
  } else {
    ok <- FALSE
    cat(sprintf("  [FAIL] %s\n         %s\n", f, res))
  }
}

sh_files <- list.files("scripts", pattern = "\\.sh$", full.names = TRUE)
cat("\n--- shell 脚本（用 bash -n 单独检查）---\n")
for (f in sh_files) cat("  ", f, "\n")

cat("\n==== 结果: ", if (ok) "全部通过 ✅" else "存在语法错误 ❌", " ====\n", sep = "")
if (!ok) quit(status = 1)
