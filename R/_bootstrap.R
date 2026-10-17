# =============================================================================
# _bootstrap.R  —— 必须在任何含中文的脚本之前 source
# -----------------------------------------------------------------------------
# 坑：本机 R 启动时 LC_CTYPE = "C"，此时 UTF-8 编码的 .R 源码会被按单字节解析，
#     导致中文字符串乱码、em-dash 变成乱码字符（例如 "—" 显示成 "b"）。
#     解决办法：在读取任何含非 ASCII 字符的脚本之前，先把 locale 切成 UTF-8。
#
# ⚠️ 本文件必须保持【纯 ASCII】，否则它自己就会被错误解析。
# 用法（每个 R 脚本的前两行）：
#     source("R/_bootstrap.R")
#     source("R/_common.R")
# =============================================================================

.loc_candidates <- c(
  "Chinese (Simplified)_China.utf8",
  "zh_CN.UTF-8",
  "Chinese_China.utf8",
  "en_US.UTF-8",
  "English_United States.utf8"
)

.loc_ok <- FALSE
for (.loc in .loc_candidates) {
  .loc_ok <- tryCatch({
    suppressWarnings(Sys.setlocale("LC_CTYPE", .loc))
    grepl("utf-?8|65001", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)
  }, error = function(e) FALSE)
  if (isTRUE(.loc_ok)) break
}
rm(.loc, .loc_candidates)

if (!isTRUE(.loc_ok)) {
  warning("未能设置 UTF-8 locale，中文标签可能乱码；",
          "建议在脚本中用英文标签。当前: ", Sys.getlocale("LC_CTYPE"))
}
