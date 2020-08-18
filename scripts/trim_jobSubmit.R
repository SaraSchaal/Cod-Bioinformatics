######################
## Cluster For Loop ##
######################

# Read in Parameter table
df =read.csv("df.csv")

# Create and submit job for each row
for(i in 1:length(df$row)){
 filename <- df$row[i] 
 #filename <- file(print(paste(df$row[i],".bash",sep="")))
 #fileConn<-file()
 fileConn<-file(print(paste(df$row[i],".bash",sep="")))
  writeLines(c("#!/bin/bash",
               "#SBATCH --nodes=1",
               "#SBATCH --tasks-per-node=1",
               paste0("#SBATCH --job-name=",filename,".txt"),
               "#SBATCH --mem=10Mb",
               "#SBATCH --mail-user=m.albecker@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --partition=short",
               "#SBATCH --time=00:05:00",
               "#SBATCH -N 1",
               "#SBATCH -n 1",
               paste0("#SBATCH --output=",df$row[i],".output"),
               paste0("#SBATCH --error=",df$row[i],".error"),
               "module load lotterhos/2019-11-15",
               paste0("Rscript --vanilla Testcode.R ",df$row[i]," ",df$replicate[i]," ", # Each parameter as tailing argument
                      df$delta_env[i]," ",df$delta_gen[i]," ",
                      df$sample_size[i]," ",df$n_genotypes[i]," ",
                      df$std_dev[i]," ",df$interaction[i])()
  ), fileConn)
  system(paste("sbatch *.sh")) # Creates bash, submits, but only runs first file. See below bash loop to run manually
}

## Bash For-loop 

# files = $(ls *bash)
# echo $files #Check to see if all files are accounted for
# for file in $files; do sbatch $file; done

