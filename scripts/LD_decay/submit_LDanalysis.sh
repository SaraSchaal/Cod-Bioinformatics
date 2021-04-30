#!/bin/bash
#SBATCH --job-name=plotLDdecay
#SBATCH --mem=1Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=clustOut/plotLDdecay.%j.out
#SBATCH --error=clustOut/plotLDdecay.%j.err
source ~/miniconda3/bin/activate lotterhos_utils_sara
Rscript --vanilla LD_analysis.R

