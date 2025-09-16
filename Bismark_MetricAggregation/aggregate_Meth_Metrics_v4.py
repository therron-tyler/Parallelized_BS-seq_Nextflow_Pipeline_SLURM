#!/usr/bin/env python3
import os
import re
import glob
import gzip
import pandas as pd

def parse_flagstat(path):
    """
    From samtools flagstat output:
      ... mapped (XX.XX% : YY.YY)
    Returns the first percentage (XX.XX).
    """
    with open(path) as f:
        for line in f:
            if "mapped (" in line:
                m = re.search(r"mapped \((\d+\.?\d*)%", line)
                if m:
                    return float(m.group(1))
    return None

def parse_trimming_report(path):
    """
    Reads both the Cutadapt summary (for raw) and the Trim Galore 'removed because too short'
    count, and computes trimmed = raw - removed_short.
    """
    raw = None
    removed_short = None
    open_fn = gzip.open if path.endswith('.gz') else open
    with open_fn(path, 'rt') as f:
        for line in f:
            m = re.search(r"Total reads processed:\s*([\d,]+)", line)
            if m:
                raw = int(m.group(1).replace(",", ""))
            m = re.search(r"Sequences removed because they became shorter.*:\s*([\d,]+)", line)
            if m:
                removed_short = int(m.group(1).replace(",", ""))
    if raw is None:
        return None, None
    trimmed = raw - (removed_short or 0)
    return raw, trimmed

def parse_bismark_summary(path):
    """
    Parses the Bismark summary.txt table to compute percent methylation at CpG:
      Methylated CpGs / (Methylated CpGs + Unmethylated CpGs) * 100
    """
    with open(path) as f:
        header = f.readline().rstrip().split("\t")
        values = f.readline().rstrip().split("\t")
    meth_idx = next((i for i,h in enumerate(header) if h.strip() == "Methylated CpGs"), None)
    unmeth_idx = next((i for i,h in enumerate(header) if h.strip() == "Unmethylated CpGs"), None)
    if meth_idx is None or unmeth_idx is None:
        return None
    try:
        meth = int(values[meth_idx].replace(",", ""))
        unmeth = int(values[unmeth_idx].replace(",", ""))
        total = meth + unmeth
        return round(meth / total * 100, 2) if total > 0 else None
    except (ValueError, IndexError):
        return None

def parse_bismark_mapping(path):
    """
    Parse Bismark summary.txt to compute percent mapped for Bismark:
      Unique Reads (remaining) / Total Reads * 100, or fallback to Aligned Reads.
    """
    with open(path) as f:
        header = f.readline().rstrip().split("\t")
        values = f.readline().rstrip().split("\t")
    total_idx   = next((i for i,h in enumerate(header) if h.strip() == "Total Reads"), None)
    unique_idx  = next((i for i,h in enumerate(header) if "Unique Reads" in h), None)
    aligned_idx = next((i for i,h in enumerate(header) if h.strip() == "Aligned Reads"), None)
    if total_idx is None or (unique_idx is None and aligned_idx is None):
        return None
    try:
        total = int(values[total_idx].replace(",", ""))
    except (IndexError, ValueError):
        return None
    if unique_idx is not None:
        try:
            unique = int(values[unique_idx].replace(",", ""))
            return round(unique / total * 100, 2)
        except (IndexError, ValueError):
            pass
    if aligned_idx is not None:
        try:
            aligned = int(values[aligned_idx].replace(",", ""))
            return round(aligned / total * 100, 2)
        except (IndexError, ValueError):
            pass
    return None

def parse_methyldackel_summary(path):
    """
    Computes global percent methylation from a MethylDackel percentMeth.bedGraph:
      average of the 4th column across all positions.
    """
    total_pct = 0.0
    count = 0
    open_fn = gzip.open if path.endswith('.gz') else open
    with open_fn(path, 'rt') as f:
        for line in f:
            if line.startswith('track'):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    pct = float(parts[3])
                    total_pct += pct
                    count += 1
                except ValueError:
                    continue
    return round(total_pct / count, 2) if count > 0 else None

def collect_metrics(root):
    rows = []
    for bis_dir in glob.glob(os.path.join(root, "metrics_*")):
        sample = os.path.basename(bis_dir).replace("metrics_", "")
        bwa_dir = os.path.join(root, sample + "_BWAmeth")

        trim_files       = glob.glob(os.path.join(bis_dir, "trimmed", "*_trimming_report.txt"))
        flagstat_files   = glob.glob(os.path.join(bwa_dir,  "align",            "*.flagstat.txt"))
        bismark_summaries= glob.glob(os.path.join(bis_dir, "bismark_reports", sample + "_summary.txt"))
        meth_reports     = glob.glob(os.path.join(bwa_dir, "meth_reports", "*.percentMeth*.bedGraph"))

        raw = trimmed = pct_map = pct_meth_bismark = pct_map_bismark = pct_meth_methdackel = None
        if trim_files:
            raw, trimmed = parse_trimming_report(trim_files[0])
        if flagstat_files:
            pct_map = parse_flagstat(flagstat_files[0])
        if bismark_summaries:
            summary = bismark_summaries[0]
            pct_meth_bismark = parse_bismark_summary(summary)
            pct_map_bismark = parse_bismark_mapping(summary)
        if meth_reports:
            pct_meth_methdackel = parse_methyldackel_summary(meth_reports[0])

        pct_trimmed = round(trimmed / raw * 100, 2) if (raw and trimmed is not None) else None

        rows.append({
            "sample": sample,
            "raw_reads": raw,
            "trimmed_reads": trimmed,
            "percent_trimmed": pct_trimmed,
            "percent_mapped_BWA-Meth": pct_map,
            "percent_mapped_Bismark": pct_map_bismark,
            "percent_methylation_CpG_Bismark": pct_meth_bismark,
            "percent_methylation_CpG_MethDackel": pct_meth_methdackel
        })
    return pd.DataFrame(rows)

if __name__ == "__main__":
    df = collect_metrics(".")
    df = df.sort_values("sample")
    df.to_csv("aggregate_metrics.tsv", sep="\t", index=False)
    print("Wrote aggregate_metrics.tsv:")
    print(df.to_string(index=False))
