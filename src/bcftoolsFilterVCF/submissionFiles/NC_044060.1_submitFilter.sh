#!/bin/bash
#SBATCH --job-name=bcftoolsFilter
#SBATCH --mem=100Mb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=../clustOut/bcftoolsFilter.%j.out
#SBATCH --error=../clustOut/bcftoolsFilter.%j.err
source activate lotterhos_utils_sara
bcftools view -O v -o ../VarCall_NC_044060.1.f.vcf --exclude-types indels --max-alleles 2 -i "MAF > 0.05 & QUAL > 30 & F_MISSING < 0.01 & FORMAT/DP > 5 & FORMAT/DP < 70"  ../../10_freebayes/outFiles/VarCall_freebayes-par.chrom_NC_044060.1.vcf &> ../logFiles/NC_044060.1-BCFtools.log
