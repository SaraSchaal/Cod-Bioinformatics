#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --job-name=picard_benchmark
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --output=picard_log.o
#SBATCH --error=picard_log.e

###cd ${SLURM_SUBMIT_DIR}
source activate lotterhos_utils_sara
picard MarkDuplicates I=../samtools_sortedBam_Out/Pop9_18093aln.sorted.bam O=../picard_Out/Pop9_18093aln.sorted.md.bam M=md_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TAGGING_POLICY=OpticalOnly &> md.pop9_18093aln.sorted.md.log