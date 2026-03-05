#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=20gb
#SBATCH --job-name=ma_interactive
#SBATCH -N 1
#SBATCH -n 10

module load R/4.2.0

cd /home/ttm3567/rootdir_scratch/20250520_BWA_Bismark_Combo_Pipe/out

Rscript Interactive_DMR_MAplot_v2.R DMRs_annotated_with_geneparts.csv DEGs_for_DMRs_20250610.csv /home/ttm3567/rootdir_scratch/20250520_BWA_Bismark_Combo_Pipe/out
