# =============================================================================
# 99_demo_fig1b_barplot.R —— 【立即可跑】验证出图代码
# -----------------------------------------------------------------------------
# 目的：在还没下载 5.1 GB 数据之前，先用论文 Fig.1b 的真实数字把绘图代码跑通，
#       确认配色、坐标轴、主题、保存路径都没问题，避免"等下载完才发现代码报错"。
# 依赖：只用本机已安装的 ggplot2 / dplyr / scales（无需 Seurat）
# 用法：Rscript R/99_demo_fig1b_barplot.R
# 产出：output/Fig1b_nuclei_per_organ_DEMO.pdf / .png
# =============================================================================

source("R/_bootstrap.R")   # 先设 UTF-8 locale，否则中文乱码
source("R/_common.R")

# 论文 Fig.1b 的真实数字（从原图高分辨率读数得到，八项相加 = 116,564）
cnt <- data.frame(
  organ = factor(names(PAPER_NUCLEI), levels = rev(names(PAPER_NUCLEI))),
  n     = as.integer(PAPER_NUCLEI)
)

cat("=== 数据自检 ===\n")
print(cnt)
cat("合计 = ", format(sum(cnt$n), big.mark = ","), " | 论文 = ",
    format(PAPER_TOTAL, big.mark = ","), "\n", sep = "")
stopifnot(sum(cnt$n) == PAPER_TOTAL)
cat("✅ 与论文摘要 116,564 完全一致\n\n")

# 复刻 Fig.1b 右侧面板：横向柱状图 + 数值标签
p <- ggplot(cnt, aes(x = organ, y = n, fill = organ)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.12, size = 3.2) +
  scale_fill_manual(values = ORGAN_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_flip() +
  labs(x = NULL, y = "Number of nuclei",
       title = "Fig. 1b (right panel) — nuclei per organ",
       subtitle = "Reproduced from Wang et al., Nature 644:722 (2025); total = 116,564",
       caption = "DEMO: plotting-code self-check only; values are the paper's original numbers, not a reproduction result") +
  theme_pub()

save_fig(p, "Fig1b_nuclei_per_organ_DEMO", width = 7, height = 4.6)

# 同时导出一份 CSV，方便直接贴进报告
write.csv(cnt[order(-cnt$n), c("organ", "n")],
          file.path(DIR_FIG, "fig1b_nuclei_per_organ_PAPER.csv"), row.names = FALSE)

cat("完成。请打开 output/Fig1b_nuclei_per_organ_DEMO.png 检查样式。\n")
