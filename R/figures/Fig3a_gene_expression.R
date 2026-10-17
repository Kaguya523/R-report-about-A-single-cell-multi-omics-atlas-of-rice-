# ============================================================
# Fig3a_gene_expression.R
#   Fig. 3a 精确复现 (Wang et al., Nature 644:722-729, 2025)
#
#   图的内容：上块 scRNA、下块 scATAC，均为「56 个细胞类型 × 65/63 个碳氮代谢基因」
#             的 z-score 热图；顶部竖排基因名 + 6 类功能色条，
#             右侧为 pathway 图例与两条独立色标。
#
#   方法：从论文 PDF 矢量图层直接数字化原热图 (56 行细胞类型 x 65/63 列基因)，
#         经原图自带色标反演为数值矩阵，再用相同色标/布局/注释在 R 中重绘。
#
#   运行：cd <仓库根目录> && Rscript R/figures/Fig3a_gene_expression.R
#         ⚠ 只需 base R + grid，不需要 Seurat —— 装包负担最轻的一张图。
#
#   输入（全部在 ref/fig3a/ 下，均为从原图反演出来的小表，随仓库提交）：
#         lut_bwr.csv              scRNA 色标（蓝-白-黄-红，列 r,g,b,value）
#         lut_scATAC_grad.csv      scATAC 色标（浅黄-红渐变，-2~2）
#         fig3a_scRNA_values.csv   56 x 65  z-score 数值矩阵
#         fig3a_scATAC_values22.csv 56 x 63 z-score 数值矩阵
#         fig3a_genes.csv          65 个基因名（列 gene）
#         fig3a_atac_genes.csv     63 个基因名（列 gene）
#         row_anno_upper.csv       上块 56 行的 organ_/group_ 颜色
#         row_anno_lower.csv       下块 56 行的 organ_/group_ 颜色
#         col_strip.csv            65 列的 pathway 色条（列 r,g,b）
#
#   输出：output/fig3a/Fig3a_scRNA_panel.png
#         output/fig3a/Fig3a_scATAC_panel.png
#         output/fig3a/Fig3a_reproduction_exact.png
#
#   右侧图例 (与原图一致)：
#     1) 基因功能分类图例: 6 色块 + 文字 + 灰色外框
#     2) scRNA 色标: 蓝-白-黄-红渐变, 刻度 -2 ~ 2
#     3) scATAC 色标: 浅黄-红渐变, 刻度 -2 ~ 2 (独立色条, 从原图渲染像素提取)
# ============================================================

## ---- 载入公共配置（_paths.R 与本脚本同目录）-------------------------------
for (.p in c("_paths.R", "R/figures/_paths.R", "../_paths.R")) {
  if (file.exists(.p)) { source(.p); break }
}
if (!exists("PROJ", inherits = FALSE))
  stop("Cannot find R/figures/_paths.R. Run this script from the repository root.")

suppressPackageStartupMessages({library(grid)})

REF  <- file.path(DIR_REF, "fig3a")
OUTD <- out_dir("fig3a")

## 输入文件自检：缺哪个就直接说清楚缺哪个
.need <- c("lut_bwr.csv", "lut_scATAC_grad.csv",
           "fig3a_scRNA_values.csv", "fig3a_scATAC_values22.csv",
           "fig3a_genes.csv", "fig3a_atac_genes.csv",
           "row_anno_upper.csv", "row_anno_lower.csv", "col_strip.csv")
for (.f in .need)
  need_file(file.path(REF, .f), "Fig. 3a digitised reference table",
            "see the file list in ref/README.md; this folder ships with the repo")

hexc <- function(rgb) rgb(rgb[1], rgb[2], rgb[3])

# ---------- LUT (从原图图例色标提取) ----------
lutRNA <- read.csv(file.path(REF, "lut_bwr.csv"))          # scRNA: 蓝-白-黄-红
lutAT  <- read.csv(file.path(REF, "lut_scATAC_grad.csv"))  # scATAC: 浅黄-红渐变 (-2~2)
col_of <- function(v, lut) {
  L <- as.matrix(lut[, c("r","g","b")]); V <- lut$value
  apply(as.matrix(v), 1:2, function(x) {
    k <- which.min((V - x)^2); hexc(L[k, ])
  })
}

# ---------- 矩阵与注释 ----------
mRNA <- as.matrix(read.csv(file.path(REF, "fig3a_scRNA_values.csv"),   header = TRUE, check.names = FALSE))
mAT  <- as.matrix(read.csv(file.path(REF, "fig3a_scATAC_values22.csv"), header = TRUE, check.names = FALSE))
genesRNA <- read.csv(file.path(REF, "fig3a_genes.csv"))$gene
genesAT  <- read.csv(file.path(REF, "fig3a_atac_genes.csv"))$gene
stopifnot(nrow(mRNA) == 56, ncol(mRNA) == 65, dim(mAT) == c(56, 63))

annoRNA <- read.csv(file.path(REF, "row_anno_upper.csv"))
annoAT  <- read.csv(file.path(REF, "row_anno_lower.csv"))
strip <- function(df, pre) {
  v <- df[, paste0(pre, c("r","g","b"))]
  sapply(seq_len(nrow(df)), function(i) {
    z <- v[i, ]
    if (all(is.na(z))) NA_character_ else hexc(as.numeric(z))
  })
}
organRNA <- strip(annoRNA, "organ_"); groupRNA <- strip(annoRNA, "group_")
organAT  <- strip(annoAT,  "organ_"); groupAT  <- strip(annoAT,  "group_")

# pathway 注释条 (65 列逐列颜色, 由原图色条游程直接读取)
cs <- read.csv(file.path(REF, "col_strip.csv"))
stopifnot(nrow(cs) == 65)
pathRNA <- apply(cs[, c("r","g","b")], 1, hexc)
xrna_pos <- seq(48.32, 310.1, length.out = 65)
xata_pos <- seq(48.32, 310.0, length.out = 63)
idx <- apply(outer(xata_pos, xrna_pos, function(a, b) abs(a - b)), 1, which.min)
pathAT <- pathRNA[idx]

# pathway 图例 (6 类): 色块颜色按原图游程顺序
pathway_names <- c("Nitrogen absorption", "Amino acid synthesis", "Photosynthesis",
                   "Sugar metabolism", "TCA cycle", "Transport related")
pathway_cols  <- c(hexc(c(0.89, 0.51, 0.98)), hexc(c(0.90, 0.56, 0.53)),
                   hexc(c(0.57, 0.71, 0.98)), hexc(c(0.47, 0.83, 0.40)),
                   hexc(c(0.48, 0.85, 0.87)), hexc(c(0.79, 0.72, 0.17)))

# ---------- 布局常量 (PDF pt 坐标, 原图实测) ----------
LAY <- list(
  x0 = 48.32, x1r = 310.10, x1a = 310.00,
  y_rna0 = 91.80, y_rna1 = 199.20,   # 上块
  y_ata0 = 206.90, y_ata1 = 314.30,  # 下块
  path_y = 84.35, path_h = 1.90,
  grp_x = 317.66, grp_w = 1.80,
  org_x = 321.87, org_w = 1.72,
  panel_x = 40.0,
  # 图例区 (原图实测)
  pl_x = 328.8, pl_w = 5.2, pl_y0 = 100.8, pl_h = 4.77,   # pathway 图例色块
  pl_tx = 336.5,                                           # pathway 文字 x
  leg_x0 = 333.2, leg_x1 = 369.3,                          # 色标条 x 范围
  leg1_y = 141.5, leg1_h = 4.0,   # scRNA 色条
  leg2_y = 163.8, leg2_h = 4.5,   # scATAC 色条
  leg_gap_title = 4.3, leg_gap_tick = 6.2
)

yy <- function(y) 330 - y   # PDF 顶底 -> grid 底顶

# ---------- 右侧图例 ----------
draw_pathway_legend <- function() {
  # 6 个色块 (带灰外框, 与原图一致)
  for (k in seq_len(6)) {
    yk <- LAY$pl_y0 + (k - 1) * LAY$pl_h
    grid.rect(x = unit(LAY$pl_x, "pt"), y = unit(yy(yk + LAY$pl_h), "pt"),
              width = unit(LAY$pl_w, "pt"), height = unit(LAY$pl_h, "pt"),
              just = c("left", "bottom"),
              gp = gpar(fill = pathway_cols[k], col = NA))
    grid.text(pathway_names[k], x = unit(LAY$pl_tx, "pt"),
              y = unit(yy(yk + LAY$pl_h / 2 + 1.0), "pt"), just = c("left", "center"),
              gp = gpar(fontsize = 5.8, fontfamily = "Arial"))
  }
  # 整体灰色外框
  grid.rect(x = unit(LAY$pl_x - 0.5, "pt"),
            y = unit(yy(LAY$pl_y0 + 6 * LAY$pl_h + 0.5), "pt"),
            width = unit(LAY$pl_w + 1.0, "pt"),
            height = unit(6 * LAY$pl_h + 1.0, "pt"),
            just = c("left", "bottom"),
            gp = gpar(fill = NA, col = "grey55", lwd = 0.4))
}

# 单个色标图例: 标题 + 水平色条 + 刻度 -2~2
draw_scale_legend <- function(lut, title, leg_y, leg_h) {
  lx <- seq(LAY$leg_x0, LAY$leg_x1, length.out = 256)
  lv <- seq(-2, 2, length.out = 256)
  lc <- col_of(lv, lut)
  for (k in seq_along(lx))
    grid.rect(x = unit(lx[k], "pt"), y = unit(yy(leg_y + leg_h), "pt"),
              width = unit((LAY$leg_x1 - LAY$leg_x0) / 255, "pt"),
              height = unit(leg_h, "pt"),
              just = c("left", "bottom"), gp = gpar(fill = lc[k], col = NA))
  # 标题 (色条上方居中)
  grid.text(title, x = unit((LAY$leg_x0 + LAY$leg_x1) / 2, "pt"),
            y = unit(yy(leg_y - LAY$leg_gap_title), "pt"),
            gp = gpar(fontsize = 6.0, fontfamily = "Arial"))
  # 刻度 -2 -1 0 1 2 (色条下方)
  for (t in c(-2, -1, 0, 1, 2)) {
    tx <- LAY$leg_x0 + (t + 2) / 4 * (LAY$leg_x1 - LAY$leg_x0)
    grid.text(format(t, trim = TRUE), x = unit(tx, "pt"),
              y = unit(yy(leg_y + leg_h + LAY$leg_gap_tick), "pt"),
              just = "center",
              gp = gpar(fontsize = 6.0, fontfamily = "Arial"))
  }
}

# ---------- 热图面板 ----------
draw_heat_panel <- function(M, genes, x0, x1, y0, y1, lut,
                            path_cols, organ_cols, group_cols,
                            panel_label, draw_genes = TRUE,
                            path_y = LAY$path_y, label_y = 52) {
  nr <- nrow(M); nc <- ncol(M)
  cw <- (x1 - x0) / nc; ch <- (y1 - y0) / nr

  cols <- col_of(as.vector(M), lut)
  dim(cols) <- c(nr, nc)

  # 基因标签 (竖排, 位于 pathway 条上方)
  if (draw_genes) {
    gcent <- x0 + (seq_len(nc) - 0.5) * cw
    for (j in seq_len(nc)) {
      grid.text(genes[j], x = unit(gcent[j], "pt"), y = unit(yy(label_y), "pt"),
                rot = 90, just = c("right", "center"),
                gp = gpar(fontsize = 5.2, fontfamily = "Arial", col = "black"))
    }
  }
  # pathway 注释条
  for (j in seq_len(nc))
    grid.rect(x = unit(x0 + (j - 1) * cw, "pt"), y = unit(yy(path_y), "pt"),
              width = unit(cw, "pt"), height = unit(LAY$path_h, "pt"),
              just = c("left", "bottom"), gp = gpar(fill = path_cols[j], col = NA))
  # 热图单元格 (白色细边)
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    grid.rect(x = unit(x0 + (j - 1) * cw, "pt"), y = unit(yy(y0 + (i - 1) * ch), "pt"),
              width = unit(cw, "pt"), height = unit(ch, "pt"),
              just = c("left", "bottom"),
              gp = gpar(fill = cols[i, j], col = "#FFFFFF", lwd = 0.18))
  }
  # group 注释条 (有段才画)
  if (!is.null(group_cols)) for (i in seq_len(nr)) if (!is.na(group_cols[i]))
    grid.rect(x = unit(LAY$grp_x, "pt"), y = unit(yy(y0 + (i - 1) * ch), "pt"),
              width = unit(LAY$grp_w, "pt"), height = unit(ch, "pt"),
              just = c("left", "bottom"), gp = gpar(fill = group_cols[i], col = NA))
  # organ 注释条 (整列)
  for (i in seq_len(nr))
    grid.rect(x = unit(LAY$org_x, "pt"), y = unit(yy(y0 + (i - 1) * ch), "pt"),
              width = unit(LAY$org_w, "pt"), height = unit(ch, "pt"),
              just = c("left", "bottom"), gp = gpar(fill = organ_cols[i], col = NA))
  # 面板标签 (左侧竖排)
  grid.text(panel_label, x = unit(LAY$panel_x, "pt"),
            y = unit(yy((y0 + y1) / 2), "pt"), rot = 90,
            gp = gpar(fontsize = 8.5, fontfamily = "Arial"))
}

# ---------- 单面板成图 ----------
draw_one <- function(file, panel) {
  png(file.path(OUTD, file), width = 398 * 4, height = 290 * 4, res = 288, bg = "white")
  pushViewport(viewport(xscale = c(20, 415), yscale = c(0, 330)))
  draw_pathway_legend()
  if (panel == "RNA") {
    draw_heat_panel(mRNA, genesRNA, LAY$x0, LAY$x1r, LAY$y_rna0, LAY$y_rna1,
                    lutRNA, pathRNA, organRNA, groupRNA, "scRNA")
    draw_scale_legend(lutRNA, "scRNA", LAY$leg1_y, LAY$leg1_h)
  } else {
    draw_heat_panel(mAT, genesAT, LAY$x0, LAY$x1a, LAY$y_ata0, LAY$y_ata1,
                    lutAT, pathAT, organAT, groupAT, "scATAC",
                    path_y = LAY$y_ata0 - 7.5, label_y = LAY$y_ata0 - 40)
    draw_scale_legend(lutAT, "scATAC", LAY$leg2_y, LAY$leg2_h)
  }
  popViewport(); dev.off()
  cat("saved", file.path(OUTD, file), "\n")
}
draw_one("Fig3a_scRNA_panel.png",  "RNA")
draw_one("Fig3a_scATAC_panel.png", "ATAC")

# ---------- 整体复刻 (上下两块, 完整右侧图例) ----------
png(file.path(OUTD, "Fig3a_reproduction_exact.png"),
    width = 398 * 4, height = 300 * 4, res = 288, bg = "white")
pushViewport(viewport(xscale = c(20, 415), yscale = c(0, 330)))
draw_pathway_legend()
draw_scale_legend(lutRNA, "scRNA",  LAY$leg1_y, LAY$leg1_h)
draw_scale_legend(lutAT,  "scATAC", LAY$leg2_y, LAY$leg2_h)
draw_heat_panel(mRNA, genesRNA, LAY$x0, LAY$x1r, LAY$y_rna0, LAY$y_rna1,
                lutRNA, pathRNA, organRNA, groupRNA, "scRNA")
draw_heat_panel(mAT, genesAT, LAY$x0, LAY$x1a, LAY$y_ata0, LAY$y_ata1,
                lutAT, pathAT, organAT, groupAT, "scATAC", draw_genes = FALSE)
popViewport(); dev.off()
cat("saved", file.path(OUTD, "Fig3a_reproduction_exact.png"), "\n")
