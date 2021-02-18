sampleNames <- read.csv("src/SampleNames.txt", header = FALSE)
colnames(sampleNames) <- "ID"

for(i in 1:nrow(sampleNames)){
  
}
  
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print("src/bedtools/submitBedtoolsCoverageCalc.sh")
  
  writeLines(c("#!/bin/bash",
               "#SBATCH --job-name=bedtoolsCovCalc",
               "#SBATCH --mem=50Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=4:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               "#SBATCH --output=bedtools/bedtools_CovCalc.%j.out",
               "#SBATCH --error=bedtools/bedtools_CovCalc.%j.err",
               "module load lotterhos/08-24-2020",
               "source activate lotterhos-py38",
               paste0("bedtools coverage , sep = """))
               
  ), fileConn)
  
  ##system(paste("sbatch submissionFiles/",filename))
}
