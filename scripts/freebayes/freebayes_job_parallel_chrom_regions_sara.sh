#!/bin/bash

#SBATCH -p short
#SBATCH --nodes 1
#SBATCH --cpus-per-task=128
#SBATCH -t 24:00:00
#SBATCH --mem=450G
#SBATCH --constraint=zen2
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --job-name="Freebayes_parallel_test"
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.sekar@northeastern.edu

begin=`date +%s`

BAMLIST=/work/rc/s.sekar/Sara_pipeline_Aln_VarCall/Freebayes_run/bamlist.txt
POPFILE=/work/rc/s.sekar/Sara_pipeline_Aln_VarCall/Freebayes_run/poplist2.txt
REF=/work/rc/s.sekar/Sara_pipeline_Aln_VarCall/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna

source ~/miniconda3/bin/activate /work/rc/s.sekar/miniconda/envs/lotterhos_variantCallers

freebayes-parallel <(fasta_generate_regions.py ${REF}.fai 100000) 128 -f ${REF} --region ${REGION} -L ${BAMLIST} --populations ${POPFILE} -m 5 -q 5 -E 3 --min-repeat-entropy 1 -n 10 -F 0.1 --vcf VarCall_freebayes-parallel.chrom_NC_044051.${REGION}.vcf

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed
