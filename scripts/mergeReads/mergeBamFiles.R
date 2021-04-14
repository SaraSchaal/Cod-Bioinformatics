sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"


  fileConn <- file(print("src/samtoolsFilter/submissionFiles/mergeBams.sh"))
  # inputs <- NULL
  # for(i in 1:nrow(sampleNames)){
  #   inputs <- paste0(inputs, paste0("./labeled_bam_Out/", sampleNames$ID[i], ".f.rg.bam "))
  # }
  writeLines(c("#!/bin/bash",
               "#SBATCH --job-name=samtoolsFilter",
               "#SBATCH --mem=1Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=10:00:00",
               "#SBATCH --cpus-per-task=32",
               "#SBATCH --output=../samtools_merge_Out/clustOut/samtoolsMerge.%j.out",
               "#SBATCH --error=../samtools_merge_Out/clustOut/samtoolsMerge.%j.err",
               "source activate lotterhos_utils_sara",
               "samtools merge -@64 -f filtered.merged.bam -b bam.list",
               
  ), fileConn)

