#!/bin/bash
#SBATCH --job-name=subGenomesBWA			      
#SBATCH --mem=2Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1                        
#SBATCH --output=indexGenomeBWA.%j.out                
#SBATCH --error=indexGenomeBWA.%j.err         

module load lotterhos/2020-08-24
source activate lotterhos-py38

bwa index GCF_902167405.1_gadMor3.0_genomic.fna