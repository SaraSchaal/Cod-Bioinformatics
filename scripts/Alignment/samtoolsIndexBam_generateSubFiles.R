#######################################
## Create Submission Files for FastP ##
#######################################

sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){
  
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/samtools/stats/", sampleNames$ID[i], "_submitStatsCalc.sh", sep="")))
  
  writeLines(c("#!/bin/bash",
               paste("#SBATCH --job-name=",sampleNames$ID[i],"_alnStatsCalc",sep=""),
               "#SBATCH --mem=20Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=4:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=samtools_indexedBam_Out/clustOut/",sampleNames$ID[i],".%j.out"),
               paste0("#SBATCH --error=samtools_indexedBam_Out/clustOut/",sampleNames$ID[i],".%j.err"),
               "module load samtools/1.9",
               paste0("samtools stats samtools_sortedBam_Out/", sampleNames$ID[i],
                      "aln.sorted.bam > samtools_alnCheck_Out/sam_stats/",
                      sampleNames$ID[i], "aln.sorted.bam.stats", sep = "")
               
  ), fileConn)
  
  ##system(paste("sbatch submissionFiles/",filename))
}


