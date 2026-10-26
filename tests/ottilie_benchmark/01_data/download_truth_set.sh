#!/usr/bin/env bash
# Download the Ottilie et al. (2022) supplementary truth set (Sup Data 4 / 5 / 7).
#
#   Data 4 (MOESM6) — 1,405 SNV/INDEL calls across 355 clones  (SNV/INDEL truth)
#   Data 5 (MOESM7) — 24 CNV events across 23 clones            (CNV/aneuploidy truth)
#   Data 7 (MOESM9) — CRISPR/Cas9 validated causal alleles      (independent ground truth)
#
# Primary host is Springer Nature static-content (direct, no bot challenge). PMC's
# copy now sits behind a JS proof-of-work page and is not curl-able.
#
# Usage:  bash 01_data/download_truth_set.sh   [OUTDIR]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-$SCRIPT_DIR/../data/supplementary}"
mkdir -p "$OUTDIR"

SN="https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-022-03076-7/MediaObjects"

declare -A FILES=(
    ["sup4_snv_indel.xlsx"]="42003_2022_3076_MOESM6_ESM.xlsx"
    ["sup5_cnv.xlsx"]="42003_2022_3076_MOESM7_ESM.xlsx"
    ["sup7_crispr.xlsx"]="42003_2022_3076_MOESM9_ESM.xlsx"
)

for local in "${!FILES[@]}"; do
    dest="$OUTDIR/$local"
    if [[ -s "$dest" ]]; then echo "  have $local"; continue; fi
    echo "  fetching $local ..."
    curl -fSL -A "Mozilla/5.0" -o "$dest" "$SN/${FILES[$local]}"
    # sanity: xlsx files are ZIP archives (start with 'PK')
    if [[ "$(head -c2 "$dest")" != "PK" ]]; then
        echo "  ERROR: $local is not a valid xlsx (got $(file -b "$dest"))" >&2
        rm -f "$dest"; exit 1
    fi
done

echo "Truth set in: $OUTDIR"
ls -lh "$OUTDIR"/*.xlsx
