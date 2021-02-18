### Using Bedtools to determine genome coverage

To calculate coverage you need your genomic elements file that have the chromosome information for your species. For the cod genome I had to subset the .gff file to just get "region" elements which are the chromosomes (Super_Scaffolds) and contigs. See this file describing these genomic elements from the cod genome: http://www.magic.re.kr/publicdb/Gadus%20morhua/GCF_902167405.1_gadMor3.0_assembly_report.txt and see R script scripts/subsetCodGFFfile.R for how I subset the original gff file.

You do not need to convert your bam file to bed. Bedtools current version is set up to take Bam files. 

Example run for coverage calculator: 
	module load lotterhos/2020-08-24
	source activate lotterhos-py38
	bedtools coverage -a Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic_scaff_contigs.gff -b Pop1_16216aln.sorted.bam -sorted > Pop1_16216switch.coverageCalc.txt


you need to use -abam flag if you are inputting a bam file instead of a bed/gff/vcf

The output file has 16 columns of information: 

	column 1 = chromosome
	column 2 = start
	column 3 = end
	column 4 = read name (first field of BAM file)
	column 5 = mapping quality (fifth field of BAM file)
	column 6 = strand
	column 7 = thickStart
	column 8 = thickEnd
	column 9 = itemRGB (three numbers comma-delimited numbers)
	column 10 = blockCount
	column 11 = blockSizes (has a trailing comma)
	column 12 = blockStarts (has a trailing comma)
	column 13 = The number of features in B that overlapped (by at least one base pair) the A interval
	column 14 = The number of bases in A that had non-zero coverage from features in B
	column 15 = The length of the entry in A
	column 16 = The fraction of bases in A that had non-zero coverage from features in B


# head of gff file: 
[schaal.s@login-01 CodGenomes]$ head -n 20 Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic_scaff_contigs.gff 
	##gff-version 3
	#!gff-spec-version 1.21
	#!processor NCBI annotwriter
	#!genome-build gadMor3.0
	#!genome-build-accession NCBI_Assembly:GCF_902167405.1
	#!annotation-source NCBI Gadus morhua Annotation Release 100
	##sequence-region NC_044048.1 1 30875876
	##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=8049
	NC_044048.1	RefSeq	region	1	30875876	.	+	.	ID=NC_044048.1:1..30875876;Dbxref=taxon:8049;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA
	NC_044049.1	RefSeq	region	1	28732775	.	+	.	ID=NC_044049.1:1..28732775;Dbxref=taxon:8049;Name=2;chromosome=2;gbkey=Src;genome=chromosome;mol_type=genomic DNA
	NC_044050.1	RefSeq	region	1	30954429	.	+	.	ID=NC_044050.1:1..30954429;Dbxref=taxon:8049;Name=3;chromosome=3;gbkey=Src;genome=chromosome;mol_type=genomic DNA

# head of output file for coverage calc:

	NC_044048.1	RefSeq	region	1	30875876	.	+	.	ID=NC_044048.1:1..30875876;Dbxref=taxon:8049;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA	6404686	29254940	30875876	0.9475015


![Bedtools Coverage Calc](../Figures/SampCoverage_chrom_genom.pdf)

