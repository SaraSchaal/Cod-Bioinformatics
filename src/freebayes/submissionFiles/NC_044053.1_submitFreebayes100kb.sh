#!/bin/bash
#SBATCH --job-name=chr6freebayes
#SBATCH --partition=short
#SBATCH --mem=150Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=64
#SBATCH --constraint=zen2
#SBATCH --output=clustOut/NC_044053.1_Freebayes.%j.out
#SBATCH --error=clustOut/NC_044053.1_Freebayes.%j.err
POPFILE=/scratch/schaal.s/CodGenomes/10_freebayes/poplist.txt
REF=/scratch/schaal.s/CodGenomes/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna
source ~/miniconda3/bin/activate /work/rc/s.sekar/miniconda/envs/lotterhos_variantCallers
freebayes-parallel regionsFiles/NC_044053.1_100kbRegions.txt 128 -f ${REF} -b ../labeled_bam_Out/mergedBam_n128_all_lot.bam --populations ${POPFILE} -m 5 -q 5 -E 3 --min-repeat-entropy 1 -n 10 -F 0.1 >> outFiles/VarCall_freebayes-par.chrom_NC_044053.1.vcf
