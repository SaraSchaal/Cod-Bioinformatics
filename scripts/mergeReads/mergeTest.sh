#!/bin/bash
#SBATCH --job-name=mergeFiles
#SBATCH --mem=2Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=2
#SBATCH --output=../picard_merge_Out/clustOut/mergeSamFiles.%j.out
#SBATCH --error=../picard_merge_Out/clustOut/mergeSamFiles.%j.err
source activate lotterhos_utils_sara
java -jar picard.jar MergeSamFiles I=Pop1_16216.f.rg.bam I=Pop2_17028.f.rg.bam O=test_merge2samps.bam