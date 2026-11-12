#!/bin/bash
# Step 1: pull lineage-barcode reads out of a Cell Ranger BAM.
#
# Keeps only reads that carry both a cell barcode (CB) and UMI (UB) tag AND contain one of the
# lineage-barcode anchors (either orientation), and writes   CB <tab> UMI <tab> SEQ   per read.
# Pre-filtering on the anchors keeps the output small (the upstream SydShafferLab script wrote
# every unmapped read).
#
# usage: step1_extract_reads.sh [options] <bam> <out.txt>
#   --all-reads        scan every read, not only unmapped ones (samtools -f 4). Useful as a QC check
#                      that no barcode reads get (partially) aligned and would be missed by -f 4.
#   -u, --upstream S   upstream anchor   (default AAAGGACGAAACACCG)
#   -d, --downstream S downstream anchor (default GTTTTAGAGCTAGAAA)
#   -t, --threads N    samtools threads  (default 4, or $THREADS)
#
# <bam> is e.g. per_sample_outs/<sample>/count/sample_alignments.bam (cellranger multi)
#            or <sample>/outs/possorted_genome_bam.bam (cellranger count).
# Requires samtools on PATH.
set -euo pipefail

UP=AAAGGACGAAACACCG
DOWN=GTTTTAGAGCTAGAAA
THREADS=${THREADS:-4}
FLAG=(-f 4); MODE="unmapped reads only"

usage() { sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-1}"; }
while [ $# -gt 0 ]; do
  case "$1" in
    --all-reads) FLAG=(); MODE="all reads"; shift ;;
    -u|--upstream) UP="$2"; shift 2 ;;
    -d|--downstream) DOWN="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage 0 ;;
    -*) echo "unknown option: $1" >&2; usage ;;
    *) break ;;
  esac
done
[ $# -eq 2 ] || usage
BAM="$1"; OUT="$2"
[ -s "$BAM" ] || { echo "BAM not found or empty: $BAM" >&2; exit 1; }
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

rc() { echo "$1" | rev | tr ACGTN TGCAN; }
mkdir -p "$(dirname "$OUT")"
echo "[$(date)] $BAM -> $OUT ($MODE; anchors $UP / $DOWN)"
samtools view -@ "$THREADS" ${FLAG[@]+"${FLAG[@]}"} "$BAM" | awk -v U="$UP" -v D="$DOWN" -v Urc="$(rc "$UP")" -v Drc="$(rc "$DOWN")" '
  {
    s = $10
    if (index(s,U)==0 && index(s,D)==0 && index(s,Urc)==0 && index(s,Drc)==0) next
    cb = ""; umi = ""
    for (i = 12; i <= NF; i++) {
      if ($i ~ /^CB:Z:/) cb = substr($i, 6)
      if ($i ~ /^UB:Z:/) umi = substr($i, 6)
    }
    if (cb != "" && umi != "") print cb "\t" umi "\t" s
  }' > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
echo "[$(date)] $(wc -l < "$OUT" | tr -d ' ') anchor-bearing reads with CB+UMI -> $OUT"
