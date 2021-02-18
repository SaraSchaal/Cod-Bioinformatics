## Samtools coverage calculator


samtools coverage -r Pop1_216aln.sorted.bam


# bedtools output
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