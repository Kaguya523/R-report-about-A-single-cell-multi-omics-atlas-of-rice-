# `ref/fig3a/` —— Fig. 3a 的 9 个数字化参考表

这 9 个 CSV 是 `R/figures/Fig3a_gene_expression.R` 的**输入**。
Fig 3a 的复现方式是「从论文 PDF 的矢量图层把原热图反演成数值矩阵，再用同一套色标重绘一遍」，
所以它们是这张图的原料，不是中间产物。

完整说明（每列的来历、怎么自己重做一遍）见 [`../README.md`](../README.md) 第 3、4 节。

| 文件名 | 形状 | 内容 |
|---|---|---|
| `fig3a_scRNA_values.csv` | 56 × 65 | 上块 z-score 矩阵 |
| `fig3a_scATAC_values22.csv` | 56 × 63 | 下块 z-score 矩阵 |
| `fig3a_genes.csv` | 65 行 | 上块基因名（列 `gene`） |
| `fig3a_atac_genes.csv` | 63 行 | 下块基因名（列 `gene`） |
| `lut_bwr.csv` | 92 级 | scRNA 色标（`value,r,g,b`，−2…2，蓝–白–黄–红） |
| `lut_scATAC_grad.csv` | 128 级 | scATAC 色标（`value,r,g,b`，−2…2，浅黄–红） |
| `row_anno_upper.csv` | 56 行 | 上块行注释（`organ_r/g/b`、`group_r/g/b`） |
| `row_anno_lower.csv` | 56 行 | 下块行注释（同上） |
| `col_strip.csv` | 65 行 | 顶部 pathway 功能色条（`x,r,g,b`） |

## 来源

论文 PDF 第 6 页 Fig. 3a 的矢量图层，逐格提取 `rect` 的填充色，
再用原图自带的色标反查成 z-score。**不是作者提供的原始数据**，是可复核的反演结果。

下块（scATAC）只有 **63** 个基因，因为**原图本身就缺 OsPGK 与 OsSUS3** 这两个，
不是数字化时漏掉的。脚本里的 `stopifnot(dim(mAT) == c(56, 63))` 守住这个事实。

## 跑一下

```bash
cd <仓库根目录>
Rscript R/figures/Fig3a_gene_expression.R
# 约 10 秒，产出 output/fig3a/ 下 3 张 PNG
```

> 这是三个脚本里**唯一不需要 Seurat / Bioconductor** 的 —— 只用 base R + `grid`。
> 适合作为「先确认环境没问题」的第一个验证点。
