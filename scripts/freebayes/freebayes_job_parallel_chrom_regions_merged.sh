#!/bin/bash

#SBATCH -p lotterhos
#SBATCH --nodes 1
#SBATCH --cpus-per-task=32
#SBATCH -t 120:00:00
#SBATCH --mem=160G
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --job-name="Freebayes_parallel_test"
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=schaal.s@northeastern.edu

begin=`date +%s`

POPFILE=/scratch/schaal.s/CodGenomes/10_freebayes/poplist.txt
REF=/scratch/schaal.s/CodGenomes/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna

source ~/miniconda3/bin/activate /work/rc/s.sekar/miniconda/envs/lotterhos_variantCallers

freebayes-parallel <(fasta_generate_regions.py ${REF}.fai 100000) 128 -f ${REF} -b ../labeled_bam_Out/mergedBam_n128_all_lot.bam --region ${REGION} --populations ${POPFILE} -m 5 -q 5 -E 3 --min-repeat-entropy 1 -n 10 -F 0.1 --vcf VarCall_freebayes-parallel.merged.chrom_NC_044051.${REGION}.vcf

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
