#!/bin/bash
#SBATCH --job-name=Pop1_16216_submitBWA.txt
#SBATCH --mem=10Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=BWA_Out/clustOut/Pop1_16216.%j.out
#SBATCH --error=BWA_Out/clustOut/Pop1_16216.%j.err
module load lotterhos/2020-08-24
source activate lotterhos-py38
bwa mem BWA_genome/GCF_902167405.1_gadMor3.0_genomic.fna FastP_Out/Pop1_16216.R1.fq.gz FastP_Out/Pop1_16216.R2.fq.gz  > BWA_Out/Pop1_16216aln.sam -O 5 -B 3
