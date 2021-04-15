Notes on Freebayes troubleshooting

Although Jon Puritz's dDocent pipeline could handle running freebayes, the way the Discovery cluster is set up is to only have jobs running for a max of 5 days. This is because it is a shared resource and Research Computing wants to ensure that individuals troubleshoot their code with check pointing to ensure no job takes longer than 5 days. This avoids individuals from using too much of the resource than is necessary at any given time. The dDocent pipeline needs to run Freebayes for an estimate 3 weeks for the number of samples/coverage data that we have (from Jon). Because of this and the fact that the cluster is not set up to handle running jobs for that long, we switched to running freebayes-parallel so we could parallelize chunks of a chromosome across multiple cpus that were sent as separate jobs. 

We first had to figure out how to use freebayes-parallel in a way that would run a job on a single computational node. To do this, we needed to identify how small the chromosome chunks needed to be to run the job in an efficient amount of time and using an efficient amount of cpus. For this Shobana and I split up the troubleshooting tasks. I ran all tests using the merged bam file for all samples and she made a bam.list file with the names of all the bam files for freebayes to run through when calling SNPs. Once we got our initial code set up, we realized a few things: the jobs were running extremely slow, with low efficiency and would have big jumps in the amount of data that was being written to the output file. 

Initial code for my merged file tests below:
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

We had many tests time out and in some cases the computation memory was maxing out. We started our tests on the short partition and only ran jobs for max of 24 hours with 0 mem called so that we just used all memory available on that node. When we saw how slow it was running, we moved our tests to the long partition (example code above) where we tested our runs on nodes with 128 cpus. We still were finding that things were running slow and inefficiently and stopped to reanalyze our code and the output files we were getting. This is when we realized we misunderstood how freebayes-parallel worked. In the documentation for freebayes-parallel they provide an addition script called fasta_generate_regions.py which splits the reference into smaller chunks to run the program on. However, we initially believed that this cuts the reference up into smaller chunks and we needed to still have the --region flag with the regions we wanted it to analyze. So for example this was the regions file we input:

```
NC_044051.1:1-100000
NC_044051.1:100001-200000
```

The way we believed this worked was that we were telling freebayes that we wanted to call SNPs on these two 100kb windows of the NC_044051.1 chromosome, but we wanted it to split it between cpus by comparing to the reference in 25000 chunks using this code ```<(fasta_generate_regions.py ref.fa.fai 25000)```. This was incorrect and we lost time troubleshooting our incorrect code. What we were doing here was essentially giving it two regions flags and it would default to the first region flag which was spliting the entire reference into 25kb chunks. This is why we were timing out and why our memory usage was skyrocketing. 

Once we figured that out, we did the meaningful troubleshooting which was identifying how small to make the chunks and which nodes to test on. We first realized the merged file was much faster to run than the bam.list, which is what we had anticipated. Then we ran two 100 kb chunks of a single chromosome on the 128 cpu node and found it had low cpu efficiency (1.4%), but finished in 55 minutes using 1 GB of memory. We then moved to running freebayes on an entire chromosome (instead of just two 100kb chunks) on the 128 cpu node. This resulted in about 75% efficiency, 50 GB of memory used and 5.5 hour run time. Although this was pretty good, Research Computing limits the number of jobs individuals can submit to the 128 cpu nodes and we wanted to see if we could get similar runtimes on the nodes with fewer cpus. So we then split the chromosome three different ways: 1) 25kb chunks, 2) 100kb chunks and 3) 5MB chunks and we ran these on the nodes with 128 cpus and 64 cpus for comparison. We found we were reaching about 90% efficiency with 9-13 hour runs using 64 cores and 90 GB of memory when we split the genome up into 100kb chunks. So this is how we proceeded and ran each chromosome. Running the 23 chromosomes across nodes with 64 cores resulted in a total runtime of approximately 36 hours (which is really great) with final files sized to be 90-120GB each. This timing will vary depending on cluster resource availability. 

A final example script for one chromsome is below:
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
example of the beginning and end of the regions file for this run:
```

NC_044048.1:1-100000
NC_044048.1:100001-200000
NC_044048.1:200001-300000
NC_044048.1:300001-400000
NC_044048.1:400001-500000
NC_044048.1:500001-600000
NC_044048.1:600001-700000
NC_044048.1:700001-800000
NC_044048.1:800001-900000
NC_044048.1:900001-1000000
NC_044048.1:1000001-1100000
NC_044048.1:1100001-1200000
NC_044048.1:1200001-1300000
NC_044048.1:1300001-1400000

.
.
.


NC_044048.1:29900001-30000000
NC_044048.1:30000001-30100000
NC_044048.1:30100001-30200000
NC_044048.1:30200001-30300000
NC_044048.1:30300001-30400000
NC_044048.1:30400001-30500000
NC_044048.1:30500001-30600000
NC_044048.1:30600001-30700000
NC_044048.1:30700001-30800000
NC_044048.1:30800001-30875876

```

R script for making regions files from GFF file:
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

The total troubleshooting for Freebayes took about 5 weeks in total with my first meeting with Research Computing about it on March 9th. We could run troubleshooting steps on the Lotterhos nodes, but we wanted to first test on nodes with higher number of available cpus to ensure that we could run the program as efficiently as possible. The Lotterhos nodes each have 32 cpus available on them and no time limits. Although we wouldn't be limited on the time constraints, we would be limited in the available cpus. After our troubleshooting, we realized this would take each chromosome a significantly longer amount of time per chromosome (20-30 hours per chromosome as opposed to 9-12 hours) and because there are only two nodes available the total run time for all chromosomes would be 25 (average hr) * (23 chromosomes / 2 nodes) = 12 days of run time. By using the short partitiion on the cluster, we had access to many more nodes and nodes that had 64 nodes which was the most optimal number of cpus according to our earlier tests. 

