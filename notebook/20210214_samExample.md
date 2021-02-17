## below is an output from one sam file after alignment to the cod genome
[schaal.s@login-01 BWA_Out]$ head Pop1_16216aln.sam 
@SQ	SN:NC_044048.1	LN:30875876
@SQ	SN:NC_044049.1	LN:28732775
@SQ	SN:NC_044050.1	LN:30954429
@SQ	SN:NC_044051.1	LN:43798135
@SQ	SN:NC_044052.1	LN:25300426
@SQ	SN:NC_044053.1	LN:27762770
@SQ	SN:NC_044054.1	LN:34137969
@SQ	SN:NC_044055.1	LN:29710654
@SQ	SN:NC_044056.1	LN:26487948
@SQ	SN:NC_044057.1	LN:27234273

The sam format starts with information about the reference sequences, read and sample information. These lines all begin with @ symbol followed by a two letter code for the type of metadata. @SQ reference info @RG read group and sample metadata and @PG programs used to create and process the files 

example read information: 
1)	218_4_2678_15447_36886	
2)	99	
3)	NC_044053.1	
4)	6873831	
5)	60	
6)	122M	
7)	6873826	
8)	117	
9)  TCTGGGTGCCGAAGCTTCTGCTCTTGTTATTCTGATTTAACAGCCAATAGAAACACAGCACCCAGCACTTAGTAACAGGGGTAAACATTCAATGCTAGATTATGCACGTCTTCCCAGGGTGA	
10)	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	
	NM:i:5	
	MD:Z:6C66T17G3C25C0	
	MC:Z:122M	
	AS:i:105	
	XS:i:20

From page 360 of Vince Buffalo's Bioinformatics Data Skills:
1) QNAME, the query name (e.g., a sequence read’s name).

2) FLAG, the bitwise flag, which contains information about the alignment. 

3) RNAME, the reference name (e.g., which sequence the query aligned to, such as a specific chromosome name like “chr1”). The reference name must be in the SAM/BAM header as an SQ entry. If the read is unaligned, this entry may be `*`.

4) POS, the position on the reference sequence (using 1-based indexing) of the first mapping base (leftmost) in the query sequence. This may be zero if the read does not align.

5) MAPQ is the mapping quality

6) CIGAR is the CIGAR string which is a specialized format for describing the alignment. 122M means that this sequence is 122 bases long and is fully aligned to the reference.

7) RNEXT and PNEXT reference name and position of a paiired-end read's partner

8) TLEN is the template length for a paired-end read

9) SEQ stores the original read sequence

10) QUAL stores the quality scores
