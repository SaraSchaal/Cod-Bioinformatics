# Samtools filtering

Before we run the SNP caller we need to filter our samples. To do this we filter the bam file using samtools view.


## From dDocent: 
```
samtools view -@32 -h -q 10 -F 0x100 -F 0x400 $1-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 32 -b
```

- @ option allocates threads for the job
- h include the header in the output
- q skip alignments with MAPQ smaller than INT[0]
- F Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with 0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with 0' (i.e. /^0[0-7]+/) [0]
- b output in the bam format

In Jon's code he filters all mapping quality of less than 10, removes PCR duplicates with the -F 0x400 flag and removes none primary alignments with the 0x100. The mawk line of code then removes hard and soft clipped reads. Finally it is exported as a bam file. 

### Altered dDocent code
Changed the mawk code to not remove the hard clipped reads because I am interested in inversion break points. This would most likely filter those out. 
Also added in picard to change read group names to be unique per sample. 

**Pop1_16216_testSubmitFilter.sh**
```
#!/bin/bash
#SBATCH --job-name=samtoolsFilter
#SBATCH --mem=750Mb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=2
#SBATCH --output=../samtools_filter_Out/clustOut/samtoolsFilter.%j.out
#SBATCH --error=../samtools_filter_Out/clustOut/samtoolsFilter.%j.err
source activate lotterhos_utils_sara
samtools view -@2 -h -q 10 -F 0x100 -F 0x400 ../picard_Out/Pop1_16216aln.sorted.md.bam | mawk '$6 !~ /[1-9][0-9].[SH]/'| samtools view -@2 -b > ../samtools_filter_Out/Pop1_16216.f.bam
wait
picard AddOrReplaceReadGroups I=../samtools_filter_Out/Pop1_16216.f.bam O=../labeled_bam_Out/Pop1_16216.f.rg.bam RGID=Pop1_16216 RGSM=Pop1_16216 RGPL=Illumina RGLB=lib1 RGPU=unit1

```

to check that the read group replacement worked use the following lines:
```
source activate lotterhos_utils_sara
samtools view -H Pop1_16216.f.rg.bam | grep '^@RG'

> @RG	ID:Pop1_16216	LB:lib1	PL:Illumina	SM:Pop1_16216	PU:unit1
```

### Jon's email

Hi Sara,

Generally, this would only affect reads that directly span an inversion boundary. Otherwise, reads will map well to an inversion, just in a potentially unexpected orientation. However, I believe you are looking at https://github.com/SaraSchaal/Cod-Bioinformatics/blob/c903b557089c41491eb95fbb8e9728205ae52565/dDocent_ngs#L359 which doesn't apply to your data.  That's for RADseq alignments with a dDocent reference.  The raw alignments (those ending with -RGmd.bam) will contain all of the alignments.  It will get later filtered in https://github.com/SaraSchaal/Cod-Bioinformatics/blob/c903b557089c41491eb95fbb8e9728205ae52565/dDocent_ngs#L381 which is less stringent, basically requiring more than half of the read to be clipped to be trimmed.  Honestly, the first filter of -q 10 removes the vast majority of those reads anyway.   You could remove change line 381 to:

samtools view -@32 -h -q 10 -F 0x100 -F 0x400 $1-RGmd.bam | mawk '$6 !~ /[1-9][0-9].[SH]/' | samtools view -@ 32 -b

if you wanted to make it less stringent, but you could test this and my guess would be it affects less than 0.5% of reads left in the alignments.  


-- 
Jon Puritz, PhD