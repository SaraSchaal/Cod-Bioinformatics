#!/bin/bash
#SBATCH --job-name=subGenomesBWA			      
#SBATCH --mem=100Mb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1                        
#SBATCH --output=submitJobsBWA.%j.out                
#SBATCH --error=submitJobsBWA.%j.err                

for file in BWA_submissionFiles/*.sh 
do 
	echo $file
	sbatch $file
done

