## SNP Filtering Criteria
BCFtools (faster than VCFtools)

### Filtering criteria
indels 
multi allelic site
MAF > 0.05 (MAF of 0.1)
found in 90% individuals (95%)
PHRED quality > 30
coverage >5x 70x<

Are there any plots we want to make or statistics that we want to output before we filter? 

### Good tutorial on vcftools filtering:
https://speciationgenomics.github.io/filtering_vcfs/

vcftools --gzvcf filename.vcf.gz --max_missing 0.99 --minQ 30 --maf 0.05 --max-alleles 2 --remove-indels --min-meanDP 5 --max-meanDP 70 --recode | gzip filename.vcf.gz

### Now try to convert to bcftools:

bcftools view -O b -o filename.vcf --exclude-types indels --max-alleles 2 -q 0.05:minor -e "QUAL < 30 || F_MISSING < 0.01 || DP > 5 || DP < 70" 


## KEL notes:
Sizes of eastern oyster genome vcf SNP files: (GB is gigabytes)
* found in 100% of individuals and maf 0.05 - 23GB
* found in 100% of individuals and maf 0.01 - 41GB
* found in 95% of individuals and maf 0.05 - 118 GB
* found in 90% of individuals and maf 0.01 - 259GB

I would opt for 99% of individuals and MAF of 0.05. That will give you plenty of SNPs for analysis, and hopefully a file size less than 50GB.

