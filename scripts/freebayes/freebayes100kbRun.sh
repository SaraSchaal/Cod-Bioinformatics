#!/bin/bash
#SBATCH -p long
#SBATCH --nodes 1
#SBATCH --cpus-per-task=128
#SBATCH -t 120:00:00
#SBATCH --constraint=zen2
#SBATCH --mem=0
#SBATCH -o noRegion.%N.%j.out
#SBATCH -e noRegion.%N.%j.err
#SBATCH --job-name="Freebayes_parallel_test"
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=schaal.s@northeastern.edu

begin=`date +%s`

POPFILE=/scratch/schaal.s/CodGenomes/10_freebayes/poplist.txt
REF=/scratch/schaal.s/CodGenomes/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna

source ~/miniconda3/bin/activate /work/rc/s.sekar/miniconda/envs/lotterhos_variantCallers

freebayes-parallel regions100kb.txt 128 -f ${REF} -b ../labeled_bam_Out/mergedBam_n128_all_lot.bam  --populations ${POPFILE} -m 5 -q 5 -E 3 --min-repeat-entropy 1 -n 10 -F 0.1 >> VarCall_freebayes-parallel.merged.chrom_NC_044051.1.vcf

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed