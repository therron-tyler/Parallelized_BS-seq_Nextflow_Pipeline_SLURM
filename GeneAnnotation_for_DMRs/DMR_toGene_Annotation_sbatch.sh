#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=20gb
#SBATCH --job-name=_SampleName_
#SBATCH -N 1
#SBATCH -n 10


module purge
module load R/4.2.0


cd /home/ttm3567/rootdir_scratch/20250520_BWA_Bismark_Combo_Pipe/out

Rscript DMR_toGene_Annotation.R
