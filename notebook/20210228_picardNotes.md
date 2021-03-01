## Picard Notes

The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. Needed to use OPTIICAL_DUPLICATE_PIXEL_DISTANCE=2500 which is for data coming from a patterned flow cell. Our whole genome data comes from a NovaSeq run which uses the patterned flow cell. I think the default for that flag is only 100 which is appropriate for unpatterned flow cells. 


```
	#!/bin/bash
	#SBATCH --nodes=1
	#SBATCH --time=8:00:00
	#SBATCH --job-name=picard_benchmark
	#SBATCH --partition=short
	#SBATCH --cpus-per-task=8
	#SBATCH --output=picard_log.o
	#SBATCH --error=picard_log.e

	###cd ${SLURM_SUBMIT_DIR}
	source activate lotterhos_utils_sara
	picard MarkDuplicates I=../samtools_sortedBam_Out/Pop9_18093aln.sorted.bam O=../picard_Out/Pop9_18093aln.sorted.md.bam M=md_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TAGGING_POLICY=OpticalOnly &> md.pop9_18093aln.sorted.md.log
```

