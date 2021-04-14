scaffolds <- read.table("src/alignment/GCF_902167405.1_gadMor3.0_assembly_report.txt", fill =TRUE)
chroms <- scaffolds[1:23,c(3,7,10)]
colnames(chroms) <- c("chromNum", "chromName", "chromLength")

for(i in 1:nrow(chroms)){
  #filename <- as.character(sampleNames$ID[i])
  fileConn <- file(print(paste("src/freebayes/submissionFiles/", chroms$chromName[i], "_submitFreebayes100kb.sh", sep="")))
  
  constraint <- if(i < 13){
                    "#SBATCH --constraint=zen2"
                } else {
                    "#SBATCH --constraint=cascadelake"
                }
  writeLines(c("#!/bin/bash",
               paste0("#SBATCH --job-name=chr", chroms$chromNum[i],"freebayes" ),
               "#SBATCH --partition=short",
               "#SBATCH --mem=150Gb",
               "#SBATCH --mail-user=schaal.s@northeastern.edu",
               "#SBATCH --mail-type=FAIL",
               "#SBATCH --time=24:00:00",
               "#SBATCH --cpus-per-task=64",
               constraint,
               paste0("#SBATCH --output=clustOut/", chroms$chromName[i] ,"_Freebayes.%j.out"),
               paste0("#SBATCH --error=clustOut/", chroms$chromName[i] ,"_Freebayes.%j.err"),
               "POPFILE=/scratch/schaal.s/CodGenomes/10_freebayes/poplist.txt",
               "REF=/scratch/schaal.s/CodGenomes/Cod_genome_data/GCF_902167405.1_gadMor3.0_genomic.fna",
               "source ~/miniconda3/bin/activate /work/rc/s.sekar/miniconda/envs/lotterhos_variantCallers",
               paste0("freebayes-parallel regionsFiles/", chroms$chromName[i], "_100kbRegions.txt 128 -f ${REF} -b ../labeled_bam_Out/mergedBam_n128_all_lot.bam --populations ${POPFILE} -m 5 -q 5 -E 3 --min-repeat-entropy 1 -n 10 -F 0.1 >> outFiles/VarCall_freebayes-par.chrom_", chroms$chromName[i], ".vcf")
            
  ), fileConn)
}
