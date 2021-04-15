#!/bin/bash
#SBATCH --job-name=picardMarkDupPop8_18131
#SBATCH --mem=10Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=../picard_Out/clustOut/picard.%j.out
#SBATCH --error=../picard_Out/clustOut/picard.%j.Pop8_18131.err
source activate lotterhos_utils_sara
picard MarkDuplicates I=../samtools_sortedBam_Out/Pop8_18131aln.sorted.bam O=../picard_Out/Pop8_18131aln.sorted.md.bam M=../picard_Out/logFiles/Pop8_18131md_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TAGGING_POLICY=OpticalOnly &> ../picard_Out/logFiles/Pop8_18131aln.sorted.md.log