######################
## Cluster For Loop ##
######################

# Read in Sample Name table
sampleNames <- read.table("src/SampleNames.txt")
colnames(sampleNames) <- "sampID"

# Create and submit job for each row
for(i in 1:nrow(sampleNames)){
 filename <- sampleNames$sampID[i] 
 fileConn<-file(print(paste(sampleNames$sampID[i],".sh",sep="")))
  writeLines(c("#!/bin/bash",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --job-name=",filename,".txt"),
               "#SBATCH --mem=10Mb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=lotterhos",
               "#SBATCH --time=24:00:00",
               "#SBATCH -N 1",
               "#SBATCH -n 1",
               paste0("#SBATCH --output=",sampleNames$sampID[i],".output"),
               paste0("#SBATCH --error=",sampleNames$sampID[i],".error"),
               "module load lotterhos/2019-11-15",
               paste0("fastp --in1 ", filename, ".1.fg.gz --in2 ", filename, ".2.fq.gz ",  
                      "--out1 ", filename, ".1.trimmed.fq.gz --out2 ", filename, ".2.trimmed.fq.gz --cut_front --cut_tail", 
                      " --cut_window_size 5 --cut_mean_quality 15 --correction $TW -q 15 -u 50" )
  ), fileConn)
  system(paste("sbatch *.sh")) # Creates bash, submits, but only runs first file. See below bash loop to run manually
}

## Bash For-loop 

# files = $(ls *bash)
# echo $files #Check to see if all files are accounted for
# for file in $files; do sbatch $file; done

