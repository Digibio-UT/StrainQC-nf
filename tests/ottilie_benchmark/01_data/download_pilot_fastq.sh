#!/usr/bin/env bash
# Download paired-end FASTQs for the four Ottilie pilot samples and write the samplesheet.
#
# Uses ENA's filereport API to resolve each SRR to its direct fastq FTP URLs, then wgets
# them (no sra-tools needed — avoids the sra-tools 3.4.x segfault). Full genome, ~2 GB total.
#
# Usage:  bash 01_data/download_pilot_fastq.sh   [OUTDIR]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-$SCRIPT_DIR/../data/fastq}"
SHEET="$SCRIPT_DIR/../data/samplesheet_pilot.csv"
mkdir -p "$OUTDIR" "$(dirname "$SHEET")"

# sample_id  SRR                (sample_id must match 04_comparison/compare_snvs.py PILOT)
SAMPLES=(
    "NODRUG-GM2:SRR10985539"
    "Doxorubicin-16-R2b:SRR10985527"
    "Carmaphycin-R9-2:SRR10985678"
    "CBR110-15-R3a:SRR10985585"
)

echo "sample,fastq_1,fastq_2" > "$SHEET"

for entry in "${SAMPLES[@]}"; do
    sid="${entry%%:*}"; srr="${entry##*:}"
    r1="$OUTDIR/${srr}_1.fastq.gz"; r2="$OUTDIR/${srr}_2.fastq.gz"
    echo "== $sid ($srr) =="

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        # Resolve fastq FTP paths from ENA (semicolon-separated R1;R2).
        ftp=$(curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr}&result=read_run&fields=fastq_ftp&format=tsv" \
              | awk 'NR==2{print $NF}')
        if [[ -z "$ftp" ]]; then echo "  ERROR: no fastq_ftp for $srr" >&2; exit 1; fi
        IFS=';' read -r f1 f2 <<< "$ftp"
        echo "  R1 $f1"; echo "  R2 $f2"
        wget -q -O "$r1" "https://${f1}"
        wget -q -O "$r2" "https://${f2}"
    else
        echo "  already downloaded"
    fi
    echo "${sid},${r1},${r2}" >> "$SHEET"
    ls -lh "$r1" "$r2" | awk '{print "  "$5, $9}'
done

echo; echo "Samplesheet: $SHEET"; cat "$SHEET"
