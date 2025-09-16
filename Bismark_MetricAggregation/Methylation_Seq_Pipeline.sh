#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -i SAMPLE -1 R1.fastq.gz -g /path/to/genome_folder -o /path/to/outdir

  -i SAMPLE        Sample name (will prefix all outputs)
  -1 R1.fastq.gz   Single-end FASTQ
  -g GENOME_DIR    Path to bisulfite‐converted reference (Bismark indexed)
  -o OUTDIR        Output directory (will be created)
  -h               Show this help
EOF
  exit 1
}

## Parse args
while getopts "i:1:g:o:h" opt; do
  case $opt in
    i) SAMPLE=$OPTARG ;;
    1) R1=$OPTARG     ;;
    g) GENOME=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    *) usage ;;
  esac
done

# make sure required args are set
[ -z "${SAMPLE:-}" ] || [ -z "${R1:-}" ] || [ -z "${GENOME:-}" ] || [ -z "${OUTDIR:-}" ] && usage

mkdir -p "$OUTDIR"/{fastqc,trimmed,align,bismark_reports}

## 1) Initial QC
fastqc -o "$OUTDIR/fastqc" "$R1"

# Reduced Representation Bisulfite Seq library 
#non-directional prep (what you’ve told Trim Galore), Cutadapt will clip those 2 bp from both ends (start and finish) of all reads matching the CAA/CGA pattern, because it doesn’t know which strand orientation they came from

## 2) RRBS trimming
trim_galore \
  --rrbs \
  --non_directional \
  --output_dir "$OUTDIR/trimmed" \
  "$R1"

#TRIM1="$OUTDIR/trimmed/${SAMPLE}_trimmed.fq.gz"
TRIM1=$(ls "$OUTDIR/trimmed/"*trimmed.fq.gz)

## 3) Post‐trim QC
fastqc -o "$OUTDIR/fastqc" "$TRIM1"

## 4) Alignment with Bismark (single‐end)
#bismark \
#  --genome "$GENOME" \
#  --bowtie2 \
#  --non_directional \
#  -o "$OUTDIR/align" \
#  "$TRIM1"

#bismark --genome "$GENOME" "$TRIM1" -o "$OUTDIR/align"

bismark \
  --genome_folder "$GENOME" \
  --non_directional \
  --single_end "$TRIM1" \
  -p 10 \
  -o "$OUTDIR/align"

ALIGNED_BAM=$(ls "$OUTDIR/align"/*.bam)


## 5) Sort & index
samtools sort -o "${OUTDIR}/align/${SAMPLE}.sorted.bam" "$ALIGNED_BAM"
samtools index "${OUTDIR}/align/${SAMPLE}.sorted.bam"

## 6) Methylation extraction (single‐end)
bismark_methylation_extractor \
  --bedGraph \
  --cytosine_report \
  --genome_folder "$GENOME" \
  --report \
  --output "$OUTDIR/bismark_reports" \
  "${OUTDIR}/align/${SAMPLE}.sorted.bam"

## 7) Summarize
cd "$OUTDIR"
bismark2summary -o "bismark_reports/${SAMPLE}_summary" ./align/*_bismark_bt2.bam

echo "✅ Done: all outputs in $OUTDIR"
