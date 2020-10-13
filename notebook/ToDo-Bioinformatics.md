## Bioinformatics To-do

* Mapping
	* BWA - Loop through samples, run for each one
	* From Buffalo page 352:	 			
		* Given that this metadata is important (and often required), you should add read group and sample metadata when you align reads to a reference. Luckily, most aligners allow you to specify this important metadata through your alignment command. For example, BWA allows (using made up files in this example):	$ bwa mem -R'@RG\tID:readgroup_id\tSM:sample_id' ref.fa in.fq 
	* Check if bwa outputs sorted sam/bam/cram files. If they are unsorted, make sure to sort them in next step.
	* Remove PCR duplicates
		* Picard Tools MarkDuplicatesWithCigar function and assign read group information with the ddOrReplaceReadGroups function 
	* (Identify poorly aligned regions
		* GATK RealignerTargetCreator and perform local realignments with IndelRealigner)
	* Jon P didn’t do that because it didn’t make a difference with FreeBayes but it does for other programs


* Samtools
	* Read Li et al 2008, review Buffalo chapter 11
	* Make sure bam files are sorted for subsequent steps
	* You may want to index bamfiles for specifically extracting locations of the inversions each individual simply needs to be identified through the SM tags in the @RG lines in the SAM header
	* Filtering on Low quality mapping filter on quality scores based on CIGAR strings that are strange or clipped, sometimes high levels of clipping will make it through
	* Talk to Jon P about filtering
	* Checking coverage across the genome
	* New version Samtools has handy built in tools, will generate histogram on command line 

* SNP calling
	* FreeBayes - accurate and works quickly, but not good when you run into n+1 situation and you have to run all your samples again
		* For oysters, it took 3 weeks to run for full dataset (Jon padded estimates -269 samples at 20x coverage 600MB) on local server Jon used 64 cores on his cluster. Need gnu-parallel
		* Designed to do all samples at once
	* GATK (standard in field) - Haplotyper
		* Is built for n+1 problem, but you generate more intermediate files
		* Gvcf file (puts into non-variant bases too)
		* Use this if ongoing project and will add samples over time
		* Jon thinks this will be slower
	* Produce 400-500 GB VCF file

* Filter VCF 
	* Jon has worked on parallelizing filtering steps
	* Use BCFtools - way faster than VCF
	* Decide on filtering criteria
		* MAF and missing across a lot of individuals
		* Filter SNPs with individual coverage < 5x or >70x or with a PHRED quality score <30
		* Remove indels
		* Filter for HWE?
			* Heterozygous across all samples need to be filtered
	* MultiQC

* LD decay analysis

* SNP thinning for LD - “core set” SNPs
	* Checking for relatedness
	* Population structure
	* Neutral parameterization

* Statistical analysis
	* Visualize haplotype structure
		* phased haplotypes and impute missing data with Beagle v4.1
		* Haplostrips (Marnetto and Huerta-Sánchez 2017) to visualize haplotype structure
	* Visual inspection of significant regions in IGV to ensure no mapping issues
		* See Buffalo Chpater 11















