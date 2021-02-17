### Using Bedtools to determine genome coverage

To calculate coverage you need your genomic elements file that have the chromosome information for your species. For the cod genome I had to subset the .gff file to just get "region" elements which are the chromosomes (Super_Scaffolds) and contigs. See this file describing these genomic elements from the cod genome: http://www.magic.re.kr/publicdb/Gadus%20morhua/GCF_902167405.1_gadMor3.0_assembly_report.txt and see R script scripts/subsetCodGFFfile.R for how I subset the original gff file.

You do not need to convert your bam file to bed. Bedtools current version is set up to take Bam files. 

Example run for coverage calculator: 
	module load lotterhos/2020-08-24
	source activate lotterhos-py38
	bedtools coverage -abam Pop1_16216aln.sorted.bam -b Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic_scaff_contigs.gff > Pop1_16216.coverageCalc.txt


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
	column 13 = coverage depth
	column 14 = # bases at depth
	column 15 = size of A
	column 16 = % of A at depth


NC_044048.1	RefSeq	region	1	30875876	.	+	.	ID=NC_044048.1:1..30875876;Dbxref=taxon:8049;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA	6404686	29254940	30875876	0.9475015