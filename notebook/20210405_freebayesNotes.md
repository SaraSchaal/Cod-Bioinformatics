# Freebayes

okay something didn't work right because our run started at the correct chromosome position, but it just kept going through the chromosome and didn't stop at our endpoint. So good news we have a chromosome and a half of data nad a full chromosome takes about 1.5 - 2 days maybe (need to confirm that)?

code that made this happen: 

1) freebayes script, freebayes_job_parallel_chrom_regions_merged_long25kb.sh
```
	#!/bin/bash

	#SBATCH -p long
	#SBATCH --nodes 1
	#SBATCH --cpus-per-task=128
	#SBATCH -t 120:00:00
	#SBATCH --constraint=zen2
	#SBATCH --mem=0
	#SBATCH -o noMem25kb.%N.%j.out
	#SBATCH -e noMem25kb.%N.%j.err
	#SBATCH --job-name="Freebayes_parallel_test"
	#SBATCH --mail-type=END,FAIL
	#SBATCH --mail-user=schaal.s@northeastern.edu

	begin=`date +%s`

	POPFILE=/scratch/schaal.s/CodGenomes/10_freebayes/poplist.txt
	REF=/scratch/schaal.s/CodGenomes/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna

	source ~/miniconda3/bin/activate /work/rc/s.sekar/miniconda/envs/lotterhos_variantCallers

	freebayes-parallel <(fasta_generate_regions.py ${REF}.fai 25000) 128 -f ${REF} -b ../labeled_bam_Out/mergedBam_n128_all_lot.bam --region ${REGION} --populations ${POPFILE} -m 5 -q 5 -E 3 --min-repeat-entropy 1 -n 10 -F 0.1 >> VarCall_freebayes-parallel.merged.25kb.chrom_NC_044051.${REGION}.vcf

	end=`date +%s`
	elapsed=`expr $end - $begin`
	echo Time taken: $elapsed

```

2) launcher script, freebayes_launcher_regions.sh
```

#!/bin/bash

for line in `cat regions100kb.txt`
do
        echo $line
        sbatch --export=REGION=${line} freebayes_job_parallel_chrom_regions_merged_long25kb.sh
done

```

3) regions file, regions100kb.txt

```
NC_044051.1:1-100000
NC_044051.1:100001-200000
```
```
-r --region <chrom>:<start_position>..<end_position>
                   Limit analysis to the specified region, 0-base coordinates,
                   end_position not included (same as BED format).
```
this seems correct. 

Based on the size of each line in the vcf file. The NC_044051.1 chromosome file should be ~525GB. One line was 12KB so with the header and the number of bases as 43798135 that is 12 KB * 43798135 bases = ~525GB + header. 


filter output file by chromosome:
```
	#!/bin/bash
	#SBATCH --job-name=bcftoolsFilter
	#SBATCH --mem=50Gb
	#SBATCH --mail-user=schaal.s@northeastern.edu
	#SBATCH --mail-type=FAIL
	#SBATCH --partition=short
	#SBATCH --time=4:00:00
	#SBATCH --nodes=1
	#SBATCH --tasks-per-node=1
	#SBATCH --output=10_freebayes/clustOut/NC_044051.1.%j.out
	#SBATCH --error=10_freebayes/clustOut/NC_044051.1.%j.err

	bcftools filter --regions NC_044051.1 > NC_044051.1.vcf

```


THIS WORKED!!! :)
```

#!/bin/bash

#SBATCH -p long
#SBATCH --nodes 1
#SBATCH --cpus-per-task=128
#SBATCH -t 120:00:00
#SBATCH --constraint=zen2
#SBATCH --mem=0
#SBATCH -o noRegion.%N.%j.out
#SBATCH -e noRegion.%N.%j.err
#SBATCH --job-name="Freebayes_parallel_test"
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=schaal.s@northeastern.edu

begin=`date +%s`

POPFILE=/scratch/schaal.s/CodGenomes/10_freebayes/poplist.txt
REF=/scratch/schaal.s/CodGenomes/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna

source ~/miniconda3/bin/activate /work/rc/s.sekar/miniconda/envs/lotterhos_variantCallers

freebayes-parallel regions100kb.txt 128 -f ${REF} -b ../labeled_bam_Out/mergedBam_n128_all_lot.bam  --populations ${POPFILE} -m 5 -q 5 -E 3 --min-repeat-entropy 1 -n 10 -F 0.1 >> VarCall_freebayes-parallel.merged.chrom_NC_044051.1.vcf

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed



```

Final script that was run per chromosome:
```

#!/bin/bash
#SBATCH --job-name=chr1freebayes
#SBATCH --partition=short
#SBATCH --mem=150Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=64
#SBATCH --constraint=zen2
#SBATCH --output=clustOut/NC_044048.1_Freebayes.%j.out
#SBATCH --error=clustOut/NC_044048.1_Freebayes.%j.err

POPFILE=/scratch/schaal.s/CodGenomes/10_freebayes/poplist.txt
REF=/scratch/schaal.s/CodGenomes/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna

source ~/miniconda3/bin/activate /work/rc/s.sekar/miniconda/envs/lotterhos_variantCallers

freebayes-parallel regionsFiles/NC_044048.1_100kbRegions.txt 128 -f ${REF} -b ../labeled_bam_Out/mergedBam_n128_all_lot.bam --populations ${POPFILE} -m 5 -q 5 -E 3 --min-repeat-entropy 1 -n 10 -F 0.1 >> outFiles/VarCall_freebayes-par.chrom_NC_044048.1.vcf

```


I made regions files for each chromosome using the following R script (subsetCodGFFfile.R):
```

scaffolds <- read.table("src/alignment/GCF_902167405.1_gadMor3.0_assembly_report.txt", fill =TRUE)
chroms <- scaffolds[1:23,c(3,7,10)]
colnames(chroms) <- c("chromNum", "chromName", "chromLength")
options(scipen = 999)
for(i in 1:nrow(chroms)){
  chr.final.loc <- as.numeric(chroms$chromLength[i])
  n <- 100000
  final.sect <-  floor(chr.final.loc/n) * n
  chr.by100kb <- seq(from = 1, to = final.sect+1, by = n)
  chr.by100kb.ends <- c((chr.by100kb[1:(length(chr.by100kb)-1)])+(n-1), chr.final.loc)
  section.names <- paste0(chroms$chromName[i], ":", chr.by100kb, "-", chr.by100kb.ends)
  df.chr.100kb <- as.matrix(section.names)
  write.table(df.chr.100kb, paste0("src/alignment/", chroms$chromName[i], "_100kbRegions.txt"), row.names = FALSE, 
              sep = "\t", col.names = FALSE, quote = FALSE)
}

```

## output stats
#### JobIDs:
18263291 chrom 1
..
18263314 chrom 23


using 64 cores (cascadelake and zen2 nodes)
we had an cpu efficiency around 90% 
max memory used was 62GB of memory
and each run was between 6 and 12 hours. Total run time was about 36 hours depending on node availability.




