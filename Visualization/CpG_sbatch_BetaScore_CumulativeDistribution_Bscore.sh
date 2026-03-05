#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=30gb
#SBATCH --job-name=_SampleName_
#SBATCH -N 1
#SBATCH -n 10

module load python/3.8.4

cd /home/ttm3567/rootdir_scratch/20250520_BWA_Bismark_Combo_Pipe/out/

# ECDF only, display on screen
python CDF_Density_Bscore_BSseq_Visualizations.py per_cpg_methylation_matrix.tsv --cdf -o CDF_plot.pdf

# Density histogram only, saving to PDF
python CDF_Density_Bscore_BSseq_Visualizations.py per_cpg_methylation_matrix.tsv --density -o density_plot.pdf

# Both in one go (same output file)
python CDF_Density_Bscore_BSseq_Visualizations.py per_cpg_methylation_matrix.tsv --cdf --density -o CDF_Density_combined.pdf
