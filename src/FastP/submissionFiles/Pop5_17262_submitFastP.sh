#!/bin/bash
#SBATCH --job-name=Pop5_17262_submitFastP.txt
#SBATCH --mem=2Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=FastP_Out/jobSum/clustOut/Pop5_17262.%j.out
#SBATCH --error=FastP_Out/jobSum/clustOut/Pop5_17262.%j.err
module load lotterhos/2020-08-24
source activate lotterhos-py38
fastp --in1 Stacks_Out/Pop5_17262.1.fq.gz --in2 Stacks_Out/Pop5_17262.2.fq.gz --out1 FastP_Out/Pop5_17262.R1.fq.gz --out2 FastP_Out/Pop5_17262.R2.fq.gz -q 15 -u 50 --trim_front1 1 --cut_front --cut_tail --disable_adapter_trimming --cut_window_size 5 --cut_mean_quality 15 -j FastP_Out/jobSum/Pop5_17262.fp.json -h FastP_Out/jobSum/Pop5_17262.fp.html &> FastP_Out/jobSum/Pop5_17262.fp.trim.log
