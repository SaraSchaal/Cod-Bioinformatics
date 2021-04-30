sampleNames <- read.csv("src/VCFnames.txt", header = FALSE)
colnames(sampleNames) <- "ID"
chromNum <- seq(from = 44048, to = 44070, by = 1)
chromNames <- paste0("NC_0", chromNum, ".1")
chromVCFfiles <- paste0("VarCall_freebayes-par.chrom_", chromNames, ".vcf" )

#VarCall_freebayes-par.chrom_NC_044060.1.vcf
for(i in 1:length(chromNames)){
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/bcftoolsFilterVCF/submissionFiles/", chromNames[i], "_submitFilter.sh", sep="")))
  
  writeLines(c("#!/bin/bash",
               "#SBATCH --job-name=bcftoolsFilter",
               "#SBATCH --mem=100Mb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=short",
               "#SBATCH --time=1:00:00",
               "#SBATCH --cpus-per-task=1",
               "#SBATCH --output=../clustOut/bcftoolsFilter.%j.out",
               "#SBATCH --error=../clustOut/bcftoolsFilter.%j.err",
               "source activate lotterhos_utils_sara",
               paste0("bcftools view -O v -o ../VarCall_", chromNames[i], ".f.vcf --exclude-types indels --max-alleles 2 -i \"MAF > 0.05 & QUAL > 30 & F_MISSING < 0.01 & FORMAT/DP > 5 & FORMAT/DP < 70\"  ../../10_freebayes/outFiles/", 
                      chromVCFfiles[i], " &> ../logFiles/", chromNames[i], "-BCFtools.log", sep = "")
         
               
  ), fileConn)
}
