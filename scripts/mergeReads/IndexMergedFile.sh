#!/bin/bash
#SBATCH --job-name=IndexBamMerged
#SBATCH --mem=50Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=../labeled_bam_Out/clustOut/mergedIndex.%j.out
#SBATCH --error=../labeled_bam_Out/clustOut/mergedIndex.%j.err
module load samtools/1.9
samtools index ../labeled_bam_Out/mergedBam_n128_all_lot.bam > ../labeled_bam_Out/mergedBam_n128_all_lot.indexed.bam