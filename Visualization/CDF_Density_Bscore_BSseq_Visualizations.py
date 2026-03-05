#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def plot_per_cpg_cdf(tsv_path, out_file=None):
    """
    Plot the empirical CDF of per-CpG B-scores for each sample,
    using seaborn.ecdfplot().
    """
    # 1) Load the matrix (chrom, pos, sample1, sample2, …)
    df = pd.read_csv(tsv_path, sep='\t')

    # 2) Get the sample columns
    samples = df.columns.tolist()[2:]

    # 3) Seaborn style
    sns.set_style("whitegrid")

    # 4) New figure
    plt.figure(figsize=(12,8))

    # 5) Plot ECDF for each sample
    for sample in samples:
        vals = df[sample].dropna()
        sns.ecdfplot(x=vals, label=sample)

    # 6) Labels
    plt.xlabel('B-score (CpG methylation %)')
    plt.ylabel('Cumulative proportion')

    # 7) Compact legend
    plt.legend(
        fontsize='small',
        frameon=True,
        bbox_to_anchor=(1.02,1),
        loc='upper left',
        borderaxespad=0.
    )
    plt.tight_layout()

    # 8) Save or show
    if out_file:
        plt.savefig(out_file)
        print(f"Wrote CDF plot to {out_file}")
    else:
        plt.show()


def plot_per_cpg_density(tsv_path, out_file=None, bins=50):
    """
    Overlay density-normalized step-histograms for each sample,
    using seaborn.histplot().
    """
    # 1) Load the matrix
    df = pd.read_csv(tsv_path, sep='\t')

    # 2) Sample columns
    samples = df.columns.tolist()[2:]

    # 3) Seaborn style
    sns.set_style("whitegrid")

    # 4) New figure
    plt.figure(figsize=(12,8))

    # 5) Plot a step histogram (density) for each sample
    for sample in samples:
        vals = df[sample].dropna()
        sns.histplot(
            x=vals,
            stat='density',   # normalize to density
            bins=bins,
            element='step',   # outline only
            fill=False,
            linewidth=1.5,
            label=sample
        )

    # 6) Labels
    plt.xlabel('B-score (CpG methylation %)')
    plt.ylabel('Density')

    # 7) Compact legend
    plt.legend(
        fontsize='small',
        frameon=True,
        bbox_to_anchor=(1.02,1),
        loc='upper left',
        borderaxespad=0.
    )
    plt.tight_layout()

    # 8) Save or show
    if out_file:
        plt.savefig(out_file)
        print(f"Wrote density histogram to {out_file}")
    else:
        plt.show()


if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description="Plot per-CpG methylation distributions (CDF and/or density)."
    )
    p.add_argument('tsv',
                   help='per_cpg_methylation_matrix.tsv (chrom, pos, sample1, sample2, …)')
    p.add_argument('-c', '--cdf',
                   action='store_true',
                   help='Make the ECDF plot')
    p.add_argument('-d', '--density',
                   action='store_true',
                   help='Make the density‐normalized step‐histogram')
    p.add_argument('-b', '--bins',
                   type=int,
                   default=50,
                   help='Number of bins for density hist (default:50)')
    p.add_argument('-o', '--out',
                   help='Output filename (PNG/PDF). If omitted, plots display on screen.')

    args = p.parse_args()

    if not (args.cdf or args.density):
        p.error("Please specify at least one of --cdf or --density")

    if args.cdf:
        plot_per_cpg_cdf(args.tsv, args.out)

    if args.density:
        plot_per_cpg_density(args.tsv, args.out, bins=args.bins)
