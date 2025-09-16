#!/usr/bin/env python3
import os
import glob
import gzip
import csv

"""
per_cpg_methylation_matrix.py

Scans all Bismark .cov(.gz) files in the directory tree,
extracts per-CpG methylation percentage (b score) for each sample,
and writes a combined TSV matrix:
  chrom\tpos\t<sample1>\t<sample2>\t...\nwhere each cell is the percent methylation at that cytosine.

Usage:
  python per_cpg_methylation_matrix.py

Output:
  per_cpg_methylation_matrix.tsv
"""

def open_cov(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')


def parse_cov(cov_path):
    """
    Generator yielding (chrom, pos, methylated, unmethylated) per line.
    """
    with open_cov(cov_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('track') or line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            # expected: chrom, start, end, percent, methylated, unmethylated
            if len(fields) < 6:
                continue
            chrom = fields[0]
            pos = fields[1]
            try:
                meth = int(fields[4])
                unmeth = int(fields[5])
            except ValueError:
                continue
            yield chrom, pos, meth, unmeth


def main():
    # locate all cov files
    patterns = ['**/*.cov', '**/*.cov.gz']
    cov_files = []
    for pat in patterns:
        cov_files.extend(glob.glob(pat, recursive=True))
    cov_files = sorted(cov_files)
    if not cov_files:
        print("No .cov or .cov.gz files found.")
        return

    # chromosome-position keys and sample order
    samples = []
    matrix = {}  # (chrom,pos) -> [b1, b2, ...]

    # first pass: initialize positions set
    positions = set()
    for cov in cov_files:
        sample = os.path.basename(cov).split('.')[0]
        samples.append(sample)
        for chrom, pos, meth, unmeth in parse_cov(cov):
            positions.add((chrom, pos))
    # sort positions
    sorted_pos = sorted(positions, key=lambda x: (x[0], int(x[1])))

    # initialize matrix rows with empty lists
    for key in sorted_pos:
        matrix[key] = []

    # second pass: fill in per-sample values
    for cov in cov_files:
        # build a lookup for this sample
        sample_map = {}
        for chrom, pos, meth, unmeth in parse_cov(cov):
            total = meth + unmeth
            min_coverage = 25  # or 5, or 3 depending on your data
            if total >= min_coverage:
                sample_map[(chrom, pos)] = round(meth / total * 100, 2)
#            if total > 0:
#                sample_map[(chrom, pos)] = round(meth / total * 100, 2)
            else:
                sample_map[(chrom, pos)] = None
        # append to each key in order
        for key in sorted_pos:
            matrix[key].append(sample_map.get(key, None))

    # Exclude CpGs missing in most samples
    min_present_samples = int(0.7 * len(samples))  # keep sites present in ≥70% of samples
    filtered_matrix = {
       key: vals for key, vals in matrix.items()
       if sum(v is not None for v in vals) >= min_present_samples
    }

    # write TSV
    out_file = 'per_cpg_methylation_matrix.tsv'
    with open(out_file, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        header = ['chrom', 'pos'] + samples
        writer.writerow(header)
        for (chrom, pos), vals in filtered_matrix.items():
            writer.writerow([chrom, pos] + vals)

    print(f"Wrote {out_file} with {len(sorted_pos)} positions and {len(samples)} samples.")


if __name__ == '__main__':
    main()
