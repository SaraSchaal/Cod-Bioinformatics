#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --time=3:00:00
#SBATCH --job-name=slimInvScript
#SBATCH --partition=lotterhos
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=invTest.out
#SBATCH --error=invTest.err
module load lotterhos/2019-11-15
srun slim --vanilla Testcode.slim 