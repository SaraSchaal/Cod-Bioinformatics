#!/bin/bash
#SBATCH --job-name=picardMarkDupPop9_18078
#SBATCH --mem=10Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=../picard_Out/clustOut/picard.%j.out
#SBATCH --error=../picard_Out/clustOut/picard.%j.Pop9_18078.err
source activate lotterhos_utils_sara
picard MarkDuplicates I=../samtools_sortedBam_Out/Pop9_18078aln.sorted.bam O=../picard_Out/Pop9_18078aln.sorted.md.bam M=../picard_Out/logFiles/Pop9_18078md_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TAGGING_POLICY=OpticalOnly &> ../picard_Out/logFiles/Pop9_18078aln.sorted.md.log
