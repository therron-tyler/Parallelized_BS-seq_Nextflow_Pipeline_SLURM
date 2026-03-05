#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 26:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.therron@northwestern.edu
#SBATCH --output=%x.%j.out
#SBATCH --mem=40gb
#SBATCH --job-name=_SampleName_
#SBATCH -N 1
#SBATCH -n 10


module load fastqc/0.12.0
module load python/3.8.4            
module load cutadapt/4.2
module load TrimGalore/0.6.10
module load bowtie2/2.4.1
module load bismark/0.21.0
module load samtools/1.10.1

bash /home/ttm3567/63_tylert/Analysis_Algorithms/Methylation_Seq_Pipeline.sh \
  -i _SampleName_ \
  -1 _fastq_ \
  -g /home/ttm3567/b1063/Reference/mm10_bismark_20250501 \
  -o _output_/_SampleName__BisulfiteSeq_Analysis

module purge
module load fastqc/0.12.0
module load python/3.8.4            
module load cutadapt/4.2
module load TrimGalore/0.6.10
module load bowtie2/2.4.1
module load bismark/0.21.0
module load samtools/1.10.1
module load bwa/0.7.17 
module load htslib/1.8

bash /home/ttm3567/63_tylert/Analysis_Algorithms/bwa_MethylationSeq_Pipeline_v1.sh \
  -i _SampleName_ \
  -1 _fastq_ \
  -g /home/ttm3567/b1063/Reference/bwa-meth_ref_mm10-2020-A/fasta \
  -o _output_/_SampleName__BWAmeth
