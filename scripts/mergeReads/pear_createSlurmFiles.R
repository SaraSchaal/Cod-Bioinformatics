##############################
#### Submit jobs for PEAR ####
##############################

# Read in Sample Name table
sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){

  filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("scripts/mergeReads/submissionFiles/", sampleNames$ID[i], "_submitPear.sh", sep="")))

  writeLines(c("#!/bin/bash",
               paste("#SBATCH --job-name=",sampleNames$ID[i],"_submitPear.txt",sep=""),
               "#SBATCH --mem=1Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=8:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=Pear_submissionFiles/outFiles/",sampleNames$ID[i],".%j.out"),
               paste0("#SBATCH --error=Pear_submissionFiles/outFiles/",sampleNames$ID[i],".%j.err"),
               "module load lotterhos/2020-08-24",
               "source activate lotterhos-py38",
               paste0("pear -f FastP_Out/",sampleNames$ID[i],".R1.fq.gz -r FastP_Out/",sampleNames$ID[i],
                      ".R2.fq.gz -o Pear_Out/",sampleNames$ID[i],".R1.fq.gz -n 50", sep = "")

  ), fileConn)
}
