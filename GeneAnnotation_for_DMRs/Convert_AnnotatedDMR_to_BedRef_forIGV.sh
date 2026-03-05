#!/bin/bash

# Usage:
# ./convert_dmr_to_bed.sh input_dmr.tsv output_dmr.bed

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_dmr.tsv output_dmr.bed"
    exit 1
fi

input_file="$1"
output_file="$2"

# Convert tab-separated DMR table into BED format (chr, start, end, name, score, strand)
# Use diff.Methy as score, NA for strand
# BED uses 0-based start, so subtract 1 from the start coordinate

awk 'BEGIN {OFS="\t"}
     NR==1 {next}  # Skip header
     {
       chrom=$1;
       start=$2 - 1;  # Convert to 0-based start for BED
       end=$3;
       name=$14 != "" ? $14 : "DMR_"NR;
       score=($8 ~ /^-?[0-9.]+$/) ? $8 : 0;
       strand=".";
       print chrom, start, end, name, score, strand
     }' "$input_file" > "$output_file"

echo "✅ BED file written to: $output_file"