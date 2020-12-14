#######################################
## Create Submission Files for FastP ##
#######################################

sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){
  
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/FastP/submissionFiles/", sampleNames$ID[i], "_submitFastP.sh", sep="")))
  
  writeLines(c("#!/bin/bash",
               paste("#SBATCH --job-name=",sampleNames$ID[i],"_submitFastP.txt",sep=""),
               "#SBATCH --mem=2Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=1:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --output=FastP_Out/jobSum/clustOut/",sampleNames$ID[i],".%j.out"),
               paste0("#SBATCH --error=FastP_Out/jobSum/clustOut/",sampleNames$ID[i],".%j.err"),
               "module load lotterhos/2020-08-24",
               "source activate lotterhos-py38",
               paste0("fastp --in1 Stacks_Out/",sampleNames$ID[i],".1.fq.gz --in2 Stacks_Out/",sampleNames$ID[i],
                      ".2.fq.gz --out1 FastP_Out/",sampleNames$ID[i],".R1.fq.gz --out2 FastP_Out/",sampleNames$ID[i],
                      ".R2.fq.gz -q 15 -u 50 --trim_front1 1 --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 -j FastP_Out/jobSum/",
                      sampleNames$ID[i],".fp.json -h FastP_Out/jobSum/",sampleNames$ID[i],
                      ".fp.html &> FastP_Out/jobSum/",sampleNames$ID[i],".fp.trim.log", sep = "")
               
  ), fileConn)
  
  ##system(paste("sbatch submissionFiles/",filename))
}

## Bash For-loop 

# files = $(ls *bash)
# echo $files #Check to see if all files are accounted for
# for file in $files; do sbatch $file; done

