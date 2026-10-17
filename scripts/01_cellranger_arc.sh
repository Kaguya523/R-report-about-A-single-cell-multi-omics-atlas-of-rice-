#!/usr/bin/env bash
# =============================================================================
# 01_cellranger_arc.sh —— 路线 B 上游：从下机 fastq 到 count matrix
# -----------------------------------------------------------------------------
# 硬件要求：>=64 GB 内存、>=32 核、>=500 GB 磁盘（单个多组学样本）
# 参考基因组：Rice MSU 7.0 —— 必须与论文一致
# 软件：CellRanger ARC v2.0.2（论文所用版本）、CellBender v0.3.2
# 论文参数：每样本约 20,000 nuclei；RNA PE150，ATAC PE50
# =============================================================================
set -euo pipefail

# ---- 配置 -------------------------------------------------------------------
CELLRANGER_ARC="cellranger-arc"          # 建议用绝对路径
REF="refdata-msu7"                        # cellranger-arc mkref 生成的参考基因组目录
SAMPLE="E03_1"                            # 样本名
FASTQ_DIR="rawdata/${SAMPLE}"             # 内含 RNA 与 ATAC 的 fastq
OUT_DIR="results"
LOCAL_CORES=32
LOCAL_MEM="120"                           # GB

mkdir -p "$OUT_DIR" "$OUT_DIR/$SAMPLE"

# ---- 1. 生成样本 CSV（CellRanger ARC 需要描述 RNA/ATAC 两类 fastq）---------
# 论文数据是 10X 多组学，一次上机同时产出 RNA 与 ATAC 两个 fastq 集合
CSV="$OUT_DIR/${SAMPLE}.csv"
cat > "$CSV" <<EOF
fastq_id,fastqs,feature_types
${SAMPLE}_GEX,${FASTQ_DIR}/GEX,Gene Expression
${SAMPLE}_ATAC,${FASTQ_DIR}/ATAC,Chromatin Accessibility
EOF
echo "[1/4] 样本表已写入: $CSV"
cat "$CSV"

# ---- 2. mkfastq（若已拿到 bcl）---------------------------------------------
# 若数据来自 SRA，可跳过此步，直接用 fasterq-dump 得到 fastq
# $CELLRANGER_ARC mkfastq --run=<bcl_dir> --csv=<samplesheet> --output-dir=fastq

# ---- 3. count：比对 + 定量 --------------------------------------------------
echo "[2/4] 运行 cellranger-arc count（约 3-8 小时/样本）"
$CELLRANGER_ARC count \
  --id="$SAMPLE" \
  --reference="$REF" \
  --libraries="$CSV" \
  --localcores=$LOCAL_CORES \
  --localmem=$LOCAL_MEM \
  --output-dir="$OUT_DIR"

# 产物：
#   $OUT_DIR/$SAMPLE/outs/filtered_feature_bc_matrix.h5   <- 表达+可及性矩阵（进 Seurat）
#   $OUT_DIR/$SAMPLE/outs/atac_fragments.tsv.gz           <- peak 分析用
#   $OUT_DIR/$SAMPLE/outs/web_summary.html                <- 报告里放这个截图
#   $OUT_DIR/$SAMPLE/outs/per_barcode_metrics.csv         <- 质控指标（TSS enrichment 等）
#   $OUT_DIR/$SAMPLE/outs/summary.csv                     <- 关键指标汇总

# ---- 4. CellBender 去环境 RNA（论文：150 epochs, lr=1e-4）------------------
echo "[3/4] CellBender remove-background"
H5="$OUT_DIR/$SAMPLE/outs/raw_feature_bc_matrix.h5"   # 注意用 raw 而非 filtered
if [ -f "$H5" ]; then
  cellbender remove-background \
    --input "$H5" \
    --output "$OUT_DIR/$SAMPLE/cellbender_out.h5" \
    --cuda \
    --epochs 150 \
    --learning-rate 0.0001
else
  echo "  未找到 raw_feature_bc_matrix.h5，跳过 CellBender（改用 filtered 矩阵）"
fi

# ---- 5. 质控指标速查（写报告用）--------------------------------------------
echo "[4/4] 关键质控指标"
if [ -f "$OUT_DIR/$SAMPLE/outs/per_barcode_metrics.csv" ]; then
  python - <<PY
import pandas as pd
df = pd.read_csv("$OUT_DIR/$SAMPLE/outs/per_barcode_metrics.csv")
df = df[df["is_cell"] == 1]
print("细胞数            :", len(df))
print("中位 UMI (RNA)    :", int(df["gex_umis_count"].median()))
print("中位基因数 (RNA)  :", int(df["gex_genes_count"].median()))
print("中位 ATAC fragments:", int(df["atac_raw_reads"].median()))
print("中位 TSS enrichment:", round(df["atac_TSS_fragments"].median(), 2) if "atac_TSS_fragments" in df else "n/a")
PY
else
  echo "  未找到 per_barcode_metrics.csv"
fi

echo
echo "完成。下一步：把 filtered_feature_bc_matrix.h5 导入 Seurat 做质控与双细胞去除。"
echo "论文 QC 阈值：nFeature_RNA > 500 且 nCount_RNA > 1000 且 percent.mt < 10"
