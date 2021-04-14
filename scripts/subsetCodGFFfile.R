## Subset GFF file

cod.gff <- read.table("src/alignment/GCF_902167405.1_gadMor3.0_genomic.gff", sep = "\t", quote = "")
scaffolds.gff <- cod.gff[cod.gff$V3 == "region",]
chrom.gff <- scaffolds.gff[1:23,]
write.table(chrom.gff, "src/alignment/GCF_902167405.1_gadMor3.0_genomic_chroms.gff", row.names = FALSE, 
            sep = "\t", col.names = FALSE, quote = FALSE)
head(cod.gff)[1:8]
levels(cod.gff$V3)
levels(cod.gff$V2)
levels(cod.gff$V1)


## split linkage group 2 up
#NC_044051.1:1-100000
#NC_044051.1:100001-200000
#43798135
scaffolds <- read.table("src/alignment/GCF_902167405.1_gadMor3.0_assembly_report.txt", fill =TRUE)
chroms <- scaffolds[1:23,c(3,7,10)]
colnames(chroms) <- c("chromNum", "chromName", "chromLength")
options(scipen = 999)
for(i in 1:nrow(chroms)){
  chr.final.loc <- as.numeric(chroms$chromLength[i])
  n <- 100000
  final.sect <-  floor(chr.final.loc/n) * n
  chr.by100kb <- seq(from = 1, to = final.sect+1, by = n)
  chr.by100kb.ends <- c((chr.by100kb[1:(length(chr.by100kb)-1)])+(n-1), chr.final.loc)
  section.names <- paste0(chroms$chromName[i], ":", chr.by100kb, "-", chr.by100kb.ends)
  df.chr.100kb <- as.matrix(section.names)
  write.table(df.chr.100kb, paste0("src/alignment/", chroms$chromName[i], "_100kbRegions.txt"), row.names = FALSE, 
              sep = "\t", col.names = FALSE, quote = FALSE)
}





