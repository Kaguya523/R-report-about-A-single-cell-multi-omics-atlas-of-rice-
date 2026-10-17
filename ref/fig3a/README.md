# `ref/fig3a/` —— Fig. 3a 的 9 个数字化参考表

这个目录需要放 9 个 CSV，是 `R/figures/Fig3a_gene_expression.R` 的**输入**。
它们由「从论文 PDF 的矢量图层把原热图反演成数值矩阵」这一步产出 ——
Fig 3a 的复现方式是**用同一套色标把数字化矩阵重绘一遍**，所以没有这些表就画不出来。

完整说明（每个文件的列名、怎么来的）见 [`../README.md`](../README.md) 第 3 节。

| 文件名 | 形状 | 内容 |
|---|---|---|
| `fig3a_scRNA_values.csv` | 56 × 65 | 上块 z-score 矩阵 |
| `fig3a_scATAC_values22.csv` | 56 × 63 | 下块 z-score 矩阵 |
| `fig3a_genes.csv` | 65 行 | 上块基因名（列 `gene`） |
| `fig3a_atac_genes.csv` | 63 行 | 下块基因名（列 `gene`） |
| `lut_bwr.csv` | — | scRNA 色标（列 `r,g,b,value`，−2…2） |
| `lut_scATAC_grad.csv` | — | scATAC 色标（列 `r,g,b,value`，−2…2） |
| `row_anno_upper.csv` | 56 行 | 上块行注释（`organ_r/g/b`、`group_r/g/b`） |
| `row_anno_lower.csv` | 56 行 | 下块行注释（同上） |
| `col_strip.csv` | 65 行 | 顶部 pathway 功能色条（列 `r,g,b`） |

补齐这 9 个文件后：

```bash
cd <仓库根目录>
Rscript R/figures/Fig3a_gene_expression.R
# 10 秒内产出 output/fig3a/Fig3a_reproduction_exact.png
```

> 这张图是三个脚本里**唯一不需要 Seurat / Bioconductor** 的 —— 只用 base R + `grid`。
> 适合作为「先确认环境没问题」的第一个验证点。
