## Stacks
- If your samples need to be demultiplexed, running the program Stacks is a necessary first step. Stacks process_radtags function searches sequence reads for barcodes that you supply and demultiplex the reads using your barcode file specifiying which barcode belongs to which sample. This program will also trim off those barcodes. 

My samples were multiplexed using 3 different unique sequences. Every sample had a unique combination of an i5 primer, i7 primer and adapter barcode. Data from the sequencing facility will demultiplex to the level of unique combinations of i5 and i7 primers. Therefore for my samples, I had 32 unique i5 and i7 primer combinations and so the sequencing facility sent me 64 fastq.gz files one for the forward sequence and one for the reverse for all 32 combinations. I will take one file as an example in the following code. 

#### process_radtags manual and specific walkthrough for process_radtags
https://catchenlab.life.illinois.edu/stacks/manual/
https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php

### Example Fastq File from Sequencing Facility 

#### Read 1
```
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
```	
#### Read 2
```
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
```

#### Stacks - process_radtags


##### barcode file
For each set of paired end files that need to be demultiplexed, a barcode file needsto be created. For example, for files i5-4-i7-11_R1_001.fastq.gz and i5-4-i7-11_R2_001.fastq.gz its barcode file is as follows:

```
ATCACG	ATCACG	Pop6_18001
CGATGT	CGATGT	Pop7_18173
TTAGGC	TTAGGC	Pop8_18148
TGGCCA	TGGCCA	Pop9_18064
ACAGTG	ACAGTG	Pop3_16228
GCCAAT	GCCAAT	Pop4_17210
CAGATC	CAGATC	Pop6_18035
ACTTGA	ACTTGA	Pop6_18038
GATCAG	GATCAG	Pop8_18130
TAGCTT	TAGCTT	Pop8_18113
GGCTAC	GGCTAC	Pop7_18190
CTTGCA	CTTGCA	Pop1_17327
```


This file contains the 12 barcodes that have the i5_4 and i7_11 index and which sample that combination belongs to. This file is tab delimited and because I had inline barcodes, meaning they were the same on both ends of the read, I have the same barcode sequence in the first and second column. Then the sample ID is on the right. 

Sample naming convention is important for downstream steps. The way my samples are names are first a Pop identifier for what a priori population they belonged to (see Sample Notes file for my population IDs). Then a 5 digit identifier for the specific sample from that population where the first two digits represent the year the sample was collected and the last 3 digits represent the indivdiual sample ID from that sampling year.

I wrote a simple R script to create the barcodefiles for each primer index combination. See file: Stacks_FileCreate.R


##### additional flags 
 
```--inline_inline ``` The barcode option flag has 6 different options and you have to use the proper flag for you your barcodes or indexes need to be demultiplexed. Because my samples have inline barcodes on either end of the sequence I need to use the --inline_inline flag. 

```-P``` If paired end sequencing was used you need to add the -P flag to indicate this.


 command line code for processing radtags on cluster
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
```
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
```


# To check status of submission on Discovery #
```
	squeue -u schaal.s
	squeue -p lotterhos
	seff <jobid> #output stats
```

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

### 20200819 - reran stacks for files with sample Pop2_17011 because it has two barcode combinations which need to be combined. Therefore, need to make sure the resulting files don't overwrite eachother. 
the 'bc' commmand takes the previous command and uses the bash calculator 

echo $(cat Pop2_17011.1.fq.gz|wc -l)/4|bc 
cat Pop2_17011a.1.fq.gz Pop2_17011b.1.fq.gz > Pop2_17011.1.fq.gz
cat i58_i711rerun/Pop2_17011a.rem.1.fq.gz i59_i76rerun/Pop2_17011b.rem.1.fq.gz > Pop2_17011.rem.1.fq.gz
cat Pop2_17011a.2.fq.gz Pop2_17011b.2.fq.gz > Pop2_17011.2.fq.gz
cat i58_i711rerun/Pop2_17011a.rem.2.fq.gz i59_i76rerun/Pop2_17011b.rem.2.fq.gz > Pop2_17011.rem.2.fq.gz

1a - 1481406
1b - 742691
2a - 1584435
2b - 786731


specifications for this sample  

	Sample 	- 17_304_Gm 
	Adapter - ACTTGA
	i5 - 2 - GCCTCTAT (ATAGAGGC Reverse Complement)
	i7 - 9 - GATCAG

## Read 1 Post Stacks

[schaal.s@login-00 Stacks_Out]$ zcat Pop3_17304.1.fq.gz | head -n 10

	@218_1_1101_1090_1000/1
	CTGTGCGTTGGCCTGCGGGCTGACTCGGTCCTGAGATGGACTGCTGTGTAGTTTGAACCATAGATTCATTATATAGAACACGGTCTCCTCTGCGCTGCTGGCCAATGGAGCCGAACGTCCGCACTGGCGGGCGGCCATCTTGCC
	+
	,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF,:FFFFFFFFF:FFFFFFF:FFFFFFFFFFFFF:FFFFFFF
	@218_1_1101_1687_1000/1
	TCCCTCTCACTCTCCTCTCCGTCTCCTCTTTTGTCCTCGTCTCTCTCCTCTCTCCCTCTCTCCCATCTCCCTCTCTATCAAGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGC
	+
	,FFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FFFF:::F,F::FFFFF
	@218_1_1101_6840_1000/1
	AAAAAATACATAGCGGCCATGGACAGGATGACCTCTATGACAATGATAGAAACAGAAAGGACGCGGAGACTCTTGAGTCATCAAGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTC

## Read 2 Post Stacks

[schaal.s@login-00 Stacks_Out]$ zcat Pop3_17304.2.fq.gz | head -n 10

	@218_1_1101_1090_1000/2
	AGGCAAGATGGCCGCCCGCCAGTGCGGACGTTCGGCTCCATTGGCCAGCAGCGCAGAGGAGACCGTGTTCTATATAATGAATCTATGGTTCAAACTACACAGCAGTCCATCTCAGGACCGAGTCAGCCCGCAGGCCAACGCACA
	+
	,FFFF:FFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFF:,FFFFF:FF,FFFFFFFFFFFFFFF,F,:FFFF:F,FFFFF:FF,FFFF,FFFFFF,FFF,FF:FFFFF:F:FFFF,FFF:F:F:FFFFFFFFFFF,:FFFFFF
	@218_1_1101_1687_1000/2
	AAGAGAGGGAGATGGGAGAGAGGGAGAGAGGAGAGAGACGAGGACAAAAGAGGAGACGGAGAGGAGAGTGAGAGGGATCAAGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGCCTCTATGTGTAGATCTCGGTGGTCGC
	+
	,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FF:
	@218_1_1101_6840_1000/2
	AGACTCAAGAGTCTCCGCGTCCTTTCTGTTTCTATCATTGTCATAGAGGTCATCCTGTCCATGGCCGCTATGTATTTTTATCAAGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGCCTCTATGTGTAGATCTCGGTGGT
