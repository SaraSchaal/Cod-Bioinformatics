#!/bin/bash


for line in `cat regions.txt`
do
	echo $line
	sbatch --export=REGION=${line} freebayes_job_parallel_chrom_regions_merged.sh
done
