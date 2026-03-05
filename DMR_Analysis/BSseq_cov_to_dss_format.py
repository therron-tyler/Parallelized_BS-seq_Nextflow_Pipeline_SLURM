#!/usr/bin/env python3
import os
import glob
import gzip
import pandas as pd
import argparse

def read_cov_for_dss(path):
    """
    Read a Bismark .cov or .cov.gz file and return
    a DataFrame with columns [chr, pos, N, X], where
      - chr: chromosome
      - pos: genomic coordinate
      - N: total coverage (methylated + unmethylated)
      - X: methylated read count
    """
    open_fn = gzip.open if path.endswith('.gz') else open
    rows = []
    with open_fn(path, 'rt') as fh:
        for line in fh:
            if not line.strip() or line.startswith('track') or line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) < 6:
                continue
            chrom = fields[0]
            pos   = int(fields[1])
            # fields[4] is methylated count, fields[5] is unmethylated
            try:
                meth   = int(fields[4])
                unmeth = int(fields[5])
            except ValueError:
                continue
            total = meth + unmeth
            rows.append((chrom, pos, total, meth))

    df = pd.DataFrame(rows, columns=['chr','pos','N','X'])
    return df

def main(cov_pattern, outdir):
    os.makedirs(outdir, exist_ok=True)

    # find all .cov and .cov.gz under cwd (or given pattern)
    cov_files = sorted(glob.glob(cov_pattern, recursive=True))
    if not cov_files:
        print(f"No .cov or .cov.gz files found with pattern: {cov_pattern}")
        return

    for cov in cov_files:
        sample = os.path.basename(cov).split('.')[0]
        df = read_cov_for_dss(cov)
        out_path = os.path.join(outdir, f"{sample}.cov_for_DSS.tsv")
        df.to_csv(out_path, sep='\t', index=False)
        print(f"Wrote {out_path}  ({len(df)} CpG positions)")

if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description="Convert Bismark .cov(.gz) files into per-sample DataFrames for DSS::makeBSseqData"
    )
    p.add_argument('--cov-pattern', '-p',
                   default='**/*.cov*',
                   help="Glob pattern (recursive) to find .cov or .cov.gz files (default '**/*.cov*')")
    p.add_argument('--outdir', '-o',
                   default='dss_input',
                   help='Directory in which to write per-sample TSVs (default: dss_input)')
    args = p.parse_args()

    main(args.cov_pattern, args.outdir)
