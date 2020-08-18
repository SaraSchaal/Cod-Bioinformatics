#### Read header of a file

# read first few lines
	zcat file.gz | head -n 1

#### Example Output File

	Sample 	- 17_304_Gm 
	Adapter - ACTTGA
	i5 - 2 - GCCTCTAT (ATAGAGGC RevCom)
	i7 - 9 - GATCAG

## Read 1
[schaal.s@login-00 CodGenomes]$ zcat i5-2-i7-9_R1_001.fastq.gz | head -n 10

	@GWNJ-1012:218:GW191226406th:1:1101:1090:1000 1:N:0:GATCAG+ATAGAGGC
	ACTTGACTGTGCGTTGGCCTGCGGGCTGACTCGGTCCTGAGATGGACTGCTGTGTAGTTTGAACCATAGATTCATTATATAGAACACGGTCTCCTCTGCGCTGCTGGCCAATGGAGCCGAACGTCCGCACTGGCGGGCGGCCATCTTGCC
	+
	FF:,FF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF,:FFFFFFFFF:FFFFFFF:FFFFFFFFFFFFF:FFFFFFF
	@GWNJ-1012:218:GW191226406th:1:1101:1687:1000 1:N:0:GATCAG+ATAGAGGC
	ACTTGATCCCTCTCACTCTCCTCTCCGTCTCCTCTTTTGTCCTCGTCTCTCTCCTCTCTCCCTCTCTCCCATCTCCCTCTCTATCAAGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGC
	+
	FFFFFF,FFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FFFF:::F,F::FFFFF
	@GWNJ-1012:218:GW191226406th:1:1101:6840:1000 1:N:0:GATCAG+ATAGAGGC
	ACTTGAAAAAAATACATAGCGGCCATGGACAGGATGACCTCTATGACAATGATAGAAACAGAAAGGACGCGGAGACTCTTGAGTCATCAAGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTC
	
## Read 2
[schaal.s@login-00 CodGenomes]$ zcat i5-2-i7-9_R2_001.fastq.gz | head -n 10

	@GWNJ-1012:218:GW191226406th:1:1101:1090:1000 2:N:0:GATCAG+ATAGAGGC
	ACTTGAAGGCAAGATGGCCGCCCGCCAGTGCGGACGTTCGGCTCCATTGGCCAGCAGCGCAGAGGAGACCGTGTTCTATATAATGAATCTATGGTTCAAACTACACAGCAGTCCATCTCAGGACCGAGTCAGCCCGCAGGCCAACGCACA
	+
	FFFFFF,FFFF:FFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFF:,FFFFF:FF,FFFFFFFFFFFFFFF,F,:FFFF:F,FFFFF:FF,FFFF,FFFFFF,FFF,FF:FFFFF:F:FFFF,FFF:F:F:FFFFFFFFFFF,:FFFFFF
	@GWNJ-1012:218:GW191226406th:1:1101:1687:1000 2:N:0:GATCAG+ATAGAGGC
	ACTTGAAAGAGAGGGAGATGGGAGAGAGGGAGAGAGGAGAGAGACGAGGACAAAAGAGGAGACGGAGAGGAGAGTGAGAGGGATCAAGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGCCTCTATGTGTAGATCTCGGTGGTCGC
	+
	FFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FF:
	@GWNJ-1012:218:GW191226406th:1:1101:6840:1000 2:N:0:GATCAG+ATAGAGGC
	ACTTGAAGACTCAAGAGTCTCCGCGTCCTTTCTGTTTCTATCATTGTCATAGAGGTCATCCTGTCCATGGCCGCTATGTATTTTTATCAAGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGCCTCTATGTGTAGATCTCGGTGGT

## Stacks
- Searches sequence reads for barcodes that you supply and demultiplex the reads using your barcode file specifiying which barcode belongs to which sample. This will also trim off those barcodes. Move to quality control after this step. 

# use this flag for the barcode type because we have an inline barcode on each side of the read that matches
	--inline_inline

# Example use on one file that only has a single barcode in it (barcodefilei52_i79)
	ACTTGA ACTTGA Pop3_18304

# command line code for processing radtags on cluster
	process_radtags -P -p Practice/GenomeFiles/ -o Practice/Stacks_Out/ -b Practice/barcodeFilei52_i79.txt --inline_inline --disable_rad_check


### Dealing with the A-Overhang 

check the barcode to see if the A is the overhang or not. 

Confusing because in the example above the first read from Read 1 has no A or T after the barcode. Why would that be?
*Easy to clip using trimmomatic the HEADCROP argument allows you to clip the read for a specified number of bases at the start of the read*

### 20200515 - zcat issue on OSX personal computer

On OSX zcat doesn't work properly. You need to use zcat < i5-2-i7-9_R1_001.fastq.gz | head -n 10

zcat < i5-2-i7-9_R1_001.fastq.gz | head -n 10 > i5-2-i7-9_R1_001_head.fastq
zcat < i5-2-i7-9_R2_001.fastq.gz | head -n 10 > i5-2-i7-9_R2_001_head.fastq

### Practice with just head of files
# command line code for processing radtags on cluster (capital P for paired end reads)

process_radtags -P -p Practice/GenomeFiles/ -o Practice/Stacks_Out/ -b Practice/barcodeFilei52_i79.txt --inline_inline --disable_rad_check

### Output files 
# If you are processing paired-end reads, then you will get four files per barcode, two for the single-end read and two for the paired-end read. For example, given barcode ACTCG, you would see the following four files
	sample_ACTCG.1.fq
	sample_ACTCG.rem.1.fq
	sample_ACTCG.2.fq
	sample_ACTCG.rem.2.fq
# The process_radtags program wants to keep the reads in phase, so that the first read in the sample_XXX.1.fq file is the mate of the first read in the sample_XXX.2.fq file. Likewise for the second pair of reads being the second record in each of the two files and so on. When one read in a pair is discarded due to low quality or a missing restriction enzyme cut site, the remaining read can't simply be output to the sample_XXX.1.fq or sample_XXX.2.fq files as it would cause the remaining reads to fall out of phase. Instead, this read is considered a remainder read and is output into the sample_XXX.rem.1.fq file if the paired-end was discarded, or the sample_XXX.rem.2.fq file if the single-end was discarded.

### Code for running on the cluster

##### CODE STARTS HERE v ######
	#!/bin/bash
	#SBATCH --job-name=stacksTest                     # sets the job name
	#SBATCH -n 1                                 	  # reserves 1 machine
	#SBATCH -N 1
	#SBATCH --mem=10Gb                                # reserves 10 GB memory
	#SBATCH --partition=lotterhos                     # requests that the job is executed in partition my partition
	#SBATCH --mail-user=schaal.s@northeastern.edu
	#SBATCH --mail-type=FAIL
	#SBATCH --time=24:00:00                           # reserves machines/cores for 24 hours.
	#SBATCH --output=stacksTest.%j.out                # sets the standard output to be stored in file, where %j is the job id
	#SBATCH --error=stacksTest.%j.err                 # sets the standard error to be stored in file

	module load lotterhos/2019-11-15

	srun process_radtags -P -p GenomeFiles/ -o ../Stacks_Out/ -b Practice/barcodeFilei52_i79.txt --inline_inline --disable_rad_check

##### CODE ENDS HERE ^ ######

# To check status of submission on Discovery #
	squeue -u schaal.s
	squeue -p lotterhos
	seff <jobid> #output stats


### Results 
Time: 				each file took approx. 7 hrs  
Memory: 			approx. 80 MB RAM
Low coverage: 		11x 
High coverage: 		25x
Average: 			approx. 15x
Notes: 				One sample seemed to fail

Delete original files?

### Advice for next user
- only need to request mem=1GB 
- if you are writing output files to the same directory make sure to have unique names for the 
folders that the original files are each being read. Otherwise each summary stat output from
process radtags will over write the last one.

