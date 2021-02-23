#!/bin/bash
#SBATCH --job-name=plotBedCov
#SBATCH --mem=50Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=bedtools_coverage/clustOut/plotting.%j.out
#SBATCH --error=bedtools_coverage/clustOut/plotting.%j.err
module load lotterhos/2020-08-24
Rscript --vanilla subsetBedtoolsCov.R
