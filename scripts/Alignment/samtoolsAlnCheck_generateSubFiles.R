#######################################
## Create Submission Files for FastP ##
#######################################

sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){
  
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/alignment/submissionFiles/", sampleNames$ID[i], "_submitAlnCheck.sh", sep="")))
  
  writeLines(c("#!/bin/bash",
               paste("#SBATCH --job-name=",sampleNames$ID[i],"_alnCheck",sep=""),
               "#SBATCH --mem=2Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=4:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=samtools_alnCheck_Out/clustOut/",sampleNames$ID[i],".%j.out"),
               paste0("#SBATCH --error=samtools_alnCheck_Out/clustOut/",sampleNames$ID[i],".%j.err"),
               "module load samtools/1.9",
               paste0("samtools view -Sbt BWA_genome/GCF_902167405.1_gadMor3.0_genomic.fna BWA_Out/",
                      sampleNames$ID[i],"aln.sam | samtools flagstat - > samtools_alnCheck_Out/",sampleNames$ID[i],
                      "alnCheck.txt",
                      sep = "")
               
  ), fileConn)
  
  ##system(paste("sbatch submissionFiles/",filename))
}


