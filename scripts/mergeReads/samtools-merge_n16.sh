#!/bin/bash

#SBATCH -p short
#SBATCH --nodes 1
#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00
#SBATCH --constraint=zen2
#SBATCH --mem=10G
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --job-name="Samtools-merge_test"
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=schaal.s@northeastern.edu

begin=`date +%s`

BAMLIST=./labeled_bam_Out/bamlist5samps.txt

source ~/miniconda3/bin/activate lotterhos_utils_sara

samtools merge merged_n16_samp5.bam -b ${BAMLIST} -@16


end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
