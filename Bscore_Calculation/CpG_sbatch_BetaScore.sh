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

python /home/ttm3567/rootdir_scratch/20250520_BWA_Bismark_Combo_Pipe/out/Per_CpG_BetaScoreCalculation_v1.py
