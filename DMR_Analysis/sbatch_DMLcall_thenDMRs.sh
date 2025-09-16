#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=80gb
#SBATCH --job-name=callDML_thenDMR
#SBATCH -N 1
#SBATCH -n 12

module load R/4.2.0

cd /home/ttm3567/rootdir_scratch/20250617_Ronen_DMR_Workflow_ALL

# 1

Rscript DSS_callDMLs_thenDMRs_v1.R \
  --groups Rowen_BulkApp_Test_GroupFile.csv \
  --grp1 MF \
  --grp2 MFplus_IL10 \
  --bsobj BSobj_Ronen_DSSanalysis.rds \
  --out MF_vs_IL10

# 2

Rscript DSS_callDMLs_thenDMRs_v1.R \
  --groups Rowen_BulkApp_Test_GroupFile.csv \
  --grp1 MF \
  --grp2 MFplus_IL4 \
  --bsobj BSobj_Ronen_DSSanalysis.rds \
  --out MF_vs_IL4

# 3

Rscript DSS_callDMLs_thenDMRs_v1.R \
  --groups Rowen_BulkApp_Test_GroupFile.csv \
  --grp1 MF \
  --grp2 MFplus_LPS-IFN \
  --bsobj BSobj_Ronen_DSSanalysis.rds \
  --out MF_vs_LPS-IFN

# 5

Rscript DSS_callDMLs_thenDMRs_v1.R \
  --groups Rowen_BulkApp_Test_GroupFile.csv \
  --grp1 MF \
  --grp2 MFplus_PMNs_no_stim \
  --bsobj BSobj_Ronen_DSSanalysis.rds \
  --out MF_vs_PMNs_no_stim

# 6 

Rscript DSS_callDMLs_thenDMRs_v1.R \
  --groups Rowen_BulkApp_Test_GroupFile.csv \
  --grp1 MF \
  --grp2 MFplus_PMNs_stim_with_IL10 \
  --bsobj BSobj_Ronen_DSSanalysis.rds \
  --out MF_vs_PMNs_stim_with_IL10

# 7

Rscript DSS_callDMLs_thenDMRs_v1.R \
  --groups Rowen_BulkApp_Test_GroupFile.csv \
  --grp1 MF \
  --grp2 MFplus_PMNs_stim_with_TNF \
  --bsobj BSobj_Ronen_DSSanalysis.rds \
  --out MF_vs_PMNs_stim_with_TNF

# 8

Rscript DSS_callDMLs_thenDMRs_v1.R \
  --groups Rowen_BulkApp_Test_GroupFile.csv \
  --grp1 MF \
  --grp2 MFplus_TNF \
  --bsobj BSobj_Ronen_DSSanalysis.rds \
  --out MF_vs_TNF













