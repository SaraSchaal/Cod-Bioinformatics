sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/samtoolsFilter/submissionFiles/", sampleNames$ID[i], "_submitFilter.sh", sep="")))
  
  partition <- if(i < 300){
                 "#SBATCH --partition=lotterhos"
               } else if(i >= 46 & i < 150){
                 "#SBATCH --partition=short"
               } else {
                 "#SBATCH --partition=express"
               }
              
  writeLines(c("#!/bin/bash",
               "#SBATCH --job-name=picardMarkDup",
               "#SBATCH --mem=5Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               partition,
               "#SBATCH --time=1:00:00",
               "#SBATCH --cpus-per-task=2",
               "#SBATCH --output=samtools_filter_Out/clustOut/samtoolsFilter.%j.out",
               "#SBATCH --error=samtools_filter_Out/clustOut/samtoolsFilter.%j.err",
               "source activate lotterhos_utils_sara",
               paste0("samtools view -@32 -h -q 10 -F 0x100 -F 0x400 ", sampleNames$ID[i], "RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 32 -b", sep = "")
               
  ), fileConn)
 
}
