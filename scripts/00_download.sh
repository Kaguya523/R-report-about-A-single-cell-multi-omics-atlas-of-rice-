#!/usr/bin/env bash
# =============================================================================
# 00_download.sh —— 下载论文公开数据（含省流量策略）
# -----------------------------------------------------------------------------
# 三个 GEO 数据集：
#   GSE232863  8 器官主图谱（scRNA + scATAC + 穗多组学）
#   GSE273875  osltpl120 突变体根 scRNA          <- 推荐先下这个（390 MB）
#   GSE245410  osf3h 突变体 bulk RNA-seq         <- 最小（681 KB）
#
# 用法：
#   bash scripts/00_download.sh minimal   # 只下最小复现集（约 0.4 GB）★推荐先跑
#   bash scripts/00_download.sh medium    # 加下主图谱 Seurat 对象（约 5.5 GB）
#   bash scripts/00_download.sh full      # 再加 8 GB RAW.tar（需 CellRanger ARC）
# =============================================================================
set -euo pipefail

MODE="${1:-minimal}"
DATA_DIR="data"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

FTP="https://ftp.ncbi.nlm.nih.gov/geo/series"

dl () {  # dl <url> <outname>
  local url="$1" out="$2"
  if [ -f "$out" ]; then
    echo "  [skip] $out 已存在 ($(du -h "$out" | cut -f1))"
    return 0
  fi
  echo "  [get ] $out"
  # curl 在本机可能因 schannel 吊销检查失败，故加 -k 并放宽重试
  curl -k -L --retry 5 --retry-delay 3 -C - -o "$out" "$url" \
    || { echo "  curl 失败，改用 python urllib"; \
         python -c "
import urllib.request, ssl, shutil, sys
ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
with urllib.request.urlopen('$url', timeout=600, context=ctx) as r, open('$out','wb') as f:
    shutil.copyfileobj(r, f)
"; }
}

echo "=== 模式: $MODE ==="

# ---- 1. 最小集：突变体 scRNA 对象（390 MB）---------------------------------
echo "[1/3] GSE273875 —— osltpl120 突变体根 scRNA 处理对象"
dl "${FTP}/GSE273nnn/GSE273875/suppl/GSE273875_LTP_NIP_label_240622.rds.gz" \
   "GSE273875_LTP_NIP_label_240622.rds.gz"

# ---- 2. bulk RNA-seq counts（681 KB）---------------------------------------
echo "[2/3] GSE245410 —— osf3h bulk RNA-seq 原始 counts"
dl "${FTP}/GSE245nnn/GSE245410/suppl/GSE245410_raw_counts.txt.gz" \
   "GSE245410_raw_counts.txt.gz"

if [ "$MODE" = "minimal" ]; then
  echo; echo "最小集下载完成。可立即运行："
  echo "  Rscript R/01_inspect_object.R"
  echo "  Rscript R/05_fig4uw_deg.R"
  exit 0
fi

# ---- 3. 主图谱 Seurat 对象（5.1 GB）----------------------------------------
echo "[3/3] GSE232863 —— 8 器官主图谱处理对象（5.1 GB，请确保磁盘充足）"
dl "${FTP}/GSE232nnn/GSE232863/suppl/GSE232863_scRNA_omics.Rds" \
   "GSE232863_scRNA_omics.Rds"

if [ "$MODE" = "medium" ]; then
  echo; echo "中等集下载完成。可运行："
  echo "  Rscript R/02_fig1b_umap_organ.R   # 需 >=32 GB 内存"
  echo "  Rscript R/03_fig1e_composition.R"
  echo "  Rscript R/04_fig2a_correlation.R"
  exit 0
fi

# ---- 4. RAW.tar（8 GB）—— 仅路线 B 需要 ------------------------------------
echo "[4/4] GSE232863_RAW.tar（8.0 GB）"
echo "⚠️  强烈建议不要整包下载。里面只有 6 个穗多组学样本，且包含 6 个 1.2 GB 的"
echo "    atac_fragments.tsv.gz（仅在需要 peak-level 分析时才要）。"
echo "    推荐只取 filtered_feature_bc_matrix.h5（每个约 100 MB）："
echo
cat <<'EOS'
  # 例：只下 E03_1 这一个多组学样本的表达矩阵
  curl -k -L -O "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8865nnn/GSM8865405/suppl/GSM8865405_E03_1_filtered_feature_bc_matrix.h5"
  # 以及它的 ATAC fragments（1.2 GB，需要 peak 分析时再下）
  curl -k -L -O "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8865nnn/GSM8865406/suppl/GSM8865406_E03_1_atac_fragments.tsv.gz"
EOS
echo
read -r -p "仍要下载完整 RAW.tar（8 GB）？输入 yes 继续: " ans
[ "$ans" = "yes" ] || { echo "已取消。"; exit 0; }
dl "${FTP}/GSE232nnn/GSE232863/suppl/GSE232863_RAW.tar" "GSE232863_RAW.tar"

# ---- 5. 官方代码包（7.75 MB + 4.8 MB）--------------------------------------
echo "[5/5] 作者官方代码（Zenodo）"
dl "https://zenodo.org/records/15501781/files/dongwei-2023%2FSingle_Cell_Multiomics_in_Rice-v1.0.zip?download=1" \
   "official_code.zip"
dl "https://zenodo.org/records/15501779/files/dongwei-2023%2FRice-Root-Cell-Type-Prediction-Tool-v1.0.zip?download=1" \
   "rsca_tool.zip"
unzip -o -q official_code.zip && unzip -o -q rsca_tool.zip

echo
echo "=== 全部完成 ==="
echo "下一步：Rscript R/01_inspect_object.R"
