#######################################
## Create Submission Files for FastP ##
#######################################

sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){
  
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/BWA/submissionFiles/", sampleNames$ID[i], "_submitBWA.sh", sep="")))
  
  writeLines(c("#!/bin/bash",
               paste("#SBATCH --job-name=",sampleNames$ID[i],"_submitBWA",sep=""),
               "#SBATCH --mem=3Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=15:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=BWA_out/clustOut/",sampleNames$ID[i],".%j.out"),
               paste0("#SBATCH --error=BWA_out/clustOut/",sampleNames$ID[i],".%j.err"),
               "module load lotterhos/2020-08-24",
               "source activate lotterhos-py38",
               paste0("bwa mem BWA_genome/GCF_902167405.1_gadMor3.0_genomic.fna FastP_Out/",
                      sampleNames$ID[i],".R1.fq.gz FastP_Out/",sampleNames$ID[i],
                      ".R2.fq.gz  > BWA_Out/",sampleNames$ID[i],
                      "aln.sam -O 5 -B 3",
                      sep = "")
               
  ), fileConn)
  
  ##system(paste("sbatch submissionFiles/",filename))
}

## Bash For-loop 

# files = $(ls *bash)
# echo $files #Check to see if all files are accounted for
# for file in $files; do sbatch $file; done

