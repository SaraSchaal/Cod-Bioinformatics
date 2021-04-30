# Checking alignment

## First step is identifying your regions to calculate coverage
For my work, I wanted to calculate coverage per chromosome, so I downloaded the GFF file from NCBI for the Atlantic cod genome I used for alignment: https://www.ncbi.nlm.nih.gov/genome/?term=txid8049[orgn]

```
cod.gff <- read.table("src/alignment/GCF_902167405.1_gadMor3.0_genomic.gff", sep = "\t", quote = "")
scaffolds.gff <- cod.gff[cod.gff$V3 == "region",]
chrom.gff <- scaffolds.gff[1:23,]
write.table(chrom.gff, "src/alignment/GCF_902167405.1_gadMor3.0_genomic_chroms.gff", row.names = FALSE, 
            sep = "\t", col.names = FALSE, quote = FALSE)
```


## Second step is calculating coverage  

```
#!/bin/bash
#SBATCH --job-name=Pop1_16216_bedCov
#SBATCH --mem=50Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=bedtools_coverage/clustOut/Pop1_16216bedCov.%j.out
#SBATCH --error=bedtools_coverage/clustOut/Pop1_16216bedCov.%j.err
module load lotterhos/2020-08-24
source activate lotterhos-py38
bedtools coverage -a Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic_chroms.gff -b samtools_sortedBam_Out/Pop1_16216aln.sorted.bam -sorted -d > bedtools_coverage/Pop1_16216.coverageCalcDflag.txt 
```
```-a``` reference genome  

```-b``` sorted bam file  

```-d``` give per base coverage  


## Third subset for just the chromosome data and the columns of interest in the output file (reduces file sizes from 111GB to 17GB this part can easily be piped so you aren't creating that large intermediate file)  
```
#!/bin/bash
#SBATCH --job-name=Pop1_16216_alnCheck
#SBATCH --mem=2Gb
#SBATCH --mail-user=schaal.s@northeastern.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=lotterhos
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=bedtools_coverage/clustOut/Pop1_16216awk.%j.out
#SBATCH --error=bedtools_coverage/clustOut/Pop1_16216awk.%j.err
awk -F"\t" '$1~/NC*/' bedtools_coverage/Pop1_16216.coverageCalcDflag.txt | awk '{print $1,$11,$12}' > bedtools_coverage/Pop1_16216.coverageCalcChr.txt
```

## Fourth run the output file through the R script for calculating averages for different window sizes and plot results in ggplot   
```

### PROCESS COVERAGE DATA TO VISUALIZE ACROSS THE GENOME ###
## Sara M. Schaal

## LOAD LIBRARIES
library(ggplot2)
library(data.table)

## USER INPUTS
setwd("/scratch/schaal.s/CodGenomes") # SET YOUR WORKING DIRECTORY
samps <- c("Pop1_16216", "Pop1_17291", "Pop4_17236", "Pop5_17278", "Pop6_18017") # INPUT SAMPLES OF INTEREST
colors <- c("steelblue2", "chartreuse3", "orchid2", "firebrick1", "goldenrod2") # INPUT COLORS 
n <- 10000 # INPUT SIZE OF WINDOWS TO CALCULATE AVERAGE COVERAGE

## PROCESS DATA 
# first step through each sample
for(i in 1:length(samps)){
  df <- fread(paste0("bedtools_coverage/", samps[i], ".coverageCalcChr.txt"), 
                                                      sep =" ", quote = "")
  colnames(df) <- c("chrom", "base", "coverage")
  chroms <- unique(df[,1])
  df.sample.data <- NULL
  # step through each chromosome and break into the increment chunks you set with n
  for(j in 1:length(chroms)){
    df.chrom <- df[df$chrom == chroms[j], ]
    extra <- nrow(df.chrom) %% n
    # this next part is taking our windows set by n and giving each window a number id
    # this grouped dataframe is made by binding the original df.chrom with grouping values for the window size you want 
    # then finding what the last grouping value would be using because it will most like not be an even increment of your n
    # for example if you divide a chromosome by your n and get 3000.3 then you can easily input the first 3000
    # grouping values in with rep (middle part of the following cbind function) then the last grouping value will be 3001 which
    # you get by rounding using ceiling in the last part of this cbind
    grouped <- cbind(df.chrom, c(rep(1:(nrow(df.chrom)/n), each = n), rep(ceiling(nrow(df.chrom)/n), extra)))
    colnames(grouped)[4] <- "grouping"
    # now bind this chromosomes data in the full dataframe
    df.sample.data <- rbind(df.sample.data, grouped)
  }
  # finally take your new dataframe and calculate the average coverage for your increments using
  # this new grouping variable and the chromosome
  df.covAve <- aggregate(coverage~grouping + chrom, data = df.sample.data, FUN = mean)
 

  #### PLOTTING ####
  pdf(paste0("figures/", samps[i], "MaxcoveragePlot", n, ".pdf"), height= 15, width=15)

  ## I make two plots because there will undoubtly be some loci that have really high coverage which makes
  # it hard to see all the spread of the majority of the data. This first graph is set to the max coverage
  # found in the data frame and the second plot is setting your y limit to a more reasonable value for your
  # data. For me 100 x was good but feel free to change to what is appropriate for your data.

  print(ggplot(data = df.covAve, aes(x = grouping, y = coverage)) +
    geom_point(col = colors[i], alpha = 0.5) +
    facet_wrap(~chrom) +
    labs(y = "Coverage", x = paste0("location every ", n, " bases"), 
         title = paste0(samps[i], "Genome Coverage up to Max Coverage")) +
    ylim(0, max(df.covAve$coverage)) +
    xlim(0, max(df.covAve$grouping)) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), 
          legend.position = "none"))
  dev.off()
  
  pdf(paste0("figures/", samps[i], "100XcoveragePlot", n, ".pdf"), height= 15, width=15)
  
  print(ggplot(data = df.covAve, aes(x = grouping, y = coverage)) +
          geom_point(col = colors[i], alpha = 0.5) +
          facet_wrap(~chrom) +
          labs(y = "Coverage", x = paste0("location every ", n, " bases"),
               title = paste0(samps[i], "Genome Coverage up to 100X")) +
          ylim(0, 100) +
          xlim(0, max(df.covAve$grouping)) +
          theme_bw() + 
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"), 
                legend.position = "none"))
  dev.off()
  
  
}
```

