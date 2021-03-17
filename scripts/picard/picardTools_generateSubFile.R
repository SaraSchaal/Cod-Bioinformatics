sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/picard/submissionFiles/", sampleNames$ID[i], "_submitPicard.sh", sep="")))
  
  partition <- if(i < 150){
                 "#SBATCH --partition=lotterhos"
               } else {
                 "#SBATCH --partition=short"
               } 
                 #else {
               #   "#SBATCH --partition=short"
               # }
              
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=picardMarkDup", sampleNames$ID[i]),
               "#SBATCH --mem=10Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               partition,
               "#SBATCH --time=4:00:00",
               "#SBATCH --cpus-per-task=1",
               "#SBATCH --output=../picard_Out/clustOut/picard.%j.out",
               paste0("#SBATCH --error=../picard_Out/clustOut/picard.%j.", sampleNames$ID[i], ".err"),
               "source activate lotterhos_utils_sara",
               paste0("picard MarkDuplicates I=../samtools_sortedBam_Out/", sampleNames$ID[i],
                      "aln.sorted.bam O=../picard_Out/", sampleNames$ID[i], "aln.sorted.md.bam M=../picard_Out/logFiles/", 
                      sampleNames$ID[i], "md_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TAGGING_POLICY=OpticalOnly &> ../picard_Out/logFiles/",
                      sampleNames$ID[i], "aln.sorted.md.log", sep = "")
               
  ), fileConn)
 
}
