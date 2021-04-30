sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/samtoolsFilter/submissionFiles/", sampleNames$ID[i], "_submitFilter.sh", sep="")))
  
  writeLines(c("#!/bin/bash",
               "#SBATCH --job-name=samtoolsFilter",
               "#SBATCH --mem=750Mb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=short",
               "#SBATCH --time=0:30:00",
               "#SBATCH --cpus-per-task=2",
               "#SBATCH --output=../samtools_filter_Out/clustOut/samtoolsFilter.%j.out",
               "#SBATCH --error=../samtools_filter_Out/clustOut/samtoolsFilter.%j.err",
               "source activate lotterhos_utils_sara",
               paste0("samtools view -@2 -h -q 10 -F 0x100 -F 0x400 ../picard_Out/", sampleNames$ID[i], 
                      "aln.sorted.md.bam | mawk '$6 !~ /[1-9][0-9].[SH]/'| samtools view -@2 -b > ../samtools_filter_Out/", 
                      sampleNames$ID[i], ".f.bam", sep = ""),
               "wait",
               paste0("picard AddOrReplaceReadGroups I=../samtools_filter_Out/", sampleNames$ID[i], 
                      ".f.bam O=../labeled_bam_Out/", sampleNames$ID[i], ".f.rg.bam RGID=",  
                      sampleNames$ID[i], " RGSM=",  sampleNames$ID[i], 
                      " RGPL=Illumina RGLB=lib1 RGPU=unit1")
               
  ), fileConn)
}
