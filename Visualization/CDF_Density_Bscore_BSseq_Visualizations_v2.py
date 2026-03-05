#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def load_sample_group_map(map_path):
    """
    Load a sample->group mapping file with columns:
      Sample, Group, [optional] Color

    Returns:
      sample_to_group: dict {sample_name: group_name}
      group_to_color: dict {group_name: color_hex} or None
      groups: list of unique group names in order of appearance
    """
    if map_path is None:
        return None, None, None

    # Let pandas infer separator (comma, tab, space, etc.)
    mdf = pd.read_csv(map_path, sep=None, engine="python")

    mdf.columns = [c.lstrip('\ufeff') for c in mdf.columns]

    # Normalize column names a bit
    cols_lower = {c.lower(): c for c in mdf.columns}
    if "sample" not in cols_lower or "group" not in cols_lower:
        raise ValueError(
            f"Mapping file must contain 'Sample' and 'Group' columns, found: {list(mdf.columns)}"
        )

    sample_col = cols_lower["sample"]
    group_col  = cols_lower["group"]
    color_col  = cols_lower.get("color", None)

    # Drop rows with missing sample/group
    mdf = mdf.dropna(subset=[sample_col, group_col])

    mdf[sample_col] = mdf[sample_col].astype(str)
    mdf[group_col]  = mdf[group_col].astype(str)

    sample_to_group = dict(zip(mdf[sample_col], mdf[group_col]))
    groups = list(dict.fromkeys(mdf[group_col]))  # preserve order, dedupe

    group_to_color = None
    if color_col is not None:
        mdf[color_col] = mdf[color_col].astype(str)
        # One color per group: take the first non-NA color for each group
        cg = (
            mdf[[group_col, color_col]]
            .dropna(subset=[color_col])
            .drop_duplicates(subset=[group_col], keep="first")
        )
        group_to_color = dict(zip(cg[group_col], cg[color_col]))

    return sample_to_group, group_to_color, groups


def plot_per_cpg_cdf(tsv_path, out_file=None, map_file=None):
    """
    Plot the empirical CDF of per-CpG B-scores.

    If map_file is provided, curves are colored by experimental GROUP
    using the 'Color' column if present.

    Mapping file format (CSV/TSV, inferred):
      Sample   Group   [Color]

    If map_file is None, curves are colored & labeled by SAMPLE (original behavior).
    """
    # 1) Load the matrix (chrom, pos, sample1, sample2, …)
    df = pd.read_csv(tsv_path, sep='\t')

    # 2) Get the sample columns (assume first 2 are chrom + pos)
    samples = df.columns.tolist()[2:]

    # 3) Optional sample->group mapping
    sample_to_group, group_to_color, groups = load_sample_group_map(map_file)

    # 4) Seaborn style
    sns.set_style("whitegrid")

    # 5) New figure
    plt.figure(figsize=(12, 8))

    if sample_to_group is None:
        # --- Original behavior: colored by SAMPLE ---
        for sample in samples:
            vals = df[sample].dropna()
            sns.ecdfplot(x=vals, label=sample)

    else:
        # --- Group-based coloring ---
        unique_groups = groups if groups is not None else sorted(set(sample_to_group.values()))

        # If user provided colors, use them; otherwise build a palette
        if group_to_color is not None:
            # Make sure every group has *some* color
            palette = sns.color_palette(n_colors=len(unique_groups))
            auto_group_to_color = {g: palette[i] for i, g in enumerate(unique_groups)}
            # Fill missing from auto palette
            for g in unique_groups:
                if g not in group_to_color:
                    group_to_color[g] = auto_group_to_color[g]
        else:
            palette = sns.color_palette(n_colors=len(unique_groups))
            group_to_color = {g: palette[i] for i, g in enumerate(unique_groups)}

        # track which groups already have a legend entry
        seen_groups = set()

        for sample in samples:
            if sample not in sample_to_group:
                print(f"Warning: sample '{sample}' not found in mapping file; skipping in CDF plot.")
                continue

            group = sample_to_group[sample]
            color = group_to_color.get(group, "gray")

            vals = df[sample].dropna()
            # Only label the FIRST sample for a given group in the legend
            label = group if group not in seen_groups else None
            if label is not None:
                seen_groups.add(group)

            sns.ecdfplot(x=vals, label=label, color=color)

    # 6) Labels
    plt.xlabel('B-score (CpG methylation %)')
    plt.ylabel('Cumulative proportion')

    # 7) Legend
    plt.legend(
        fontsize='small',
        frameon=True,
        bbox_to_anchor=(1.02, 1),
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
    df = pd.read_csv(tsv_path, sep='\t')
    samples = df.columns.tolist()[2:]

    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 8))

    for sample in samples:
        vals = df[sample].dropna()
        sns.histplot(
            x=vals,
            stat='density',
            bins=bins,
            element='step',
            fill=False,
            linewidth=1.5,
            label=sample
        )

    plt.xlabel('B-score (CpG methylation %)')
    plt.ylabel('Density')

    plt.legend(
        fontsize='small',
        frameon=True,
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        borderaxespad=0.
    )
    plt.tight_layout()

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
                   help='Make the density-normalized step-histogram')
    p.add_argument('-b', '--bins',
                   type=int,
                   default=50,
                   help='Number of bins for density hist (default: 50)')
    p.add_argument('-o', '--out',
                   help='Output filename (PNG/PDF). If omitted, plots display on screen.')
    p.add_argument('-m', '--map',
                   help='Sample–Group mapping file (CSV/TSV) with columns '
                        'Sample, Group, [Color]. Used for CDF group coloring.')

    args = p.parse_args()

    if not (args.cdf or args.density):
        p.error("Please specify at least one of --cdf or --density")

    if args.cdf:
        plot_per_cpg_cdf(args.tsv, args.out, map_file=args.map)

    if args.density:
        plot_per_cpg_density(args.tsv, args.out, bins=args.bins)
