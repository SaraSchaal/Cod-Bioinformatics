#!/bin/bash
#SBATCH --job-name=alnQualBWA			      
#SBATCH --mem=100Mb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1                        
#SBATCH --output=submitJobsBWA.%j.out                
#SBATCH --error=submitJobsBWA.%j.err       

samtools view -Sbt BWA_genome/GCF_902167405.1_gadMor3.0_genomic.fna BWA_Out/Pop1_16216alnDef.sam | samtools flagstat –