## Subset GFF file

cod.gff <- read.table("src/alignment/GCF_902167405.1_gadMor3.0_genomic.gff", sep = "\t", quote = "")
chrom.gff <- cod.gff[cod.gff$V3 == "region",]
write.table(chrom.gff, "src/alignment/GCF_902167405.1_gadMor3.0_genomic_scaff_contigs.gff", row.names = FALSE, 
            sep = "\t", col.names = FALSE, quote = FALSE)
head(cod.gff)[1:8]
levels(cod.gff$V3)
levels(cod.gff$V2)
levels(cod.gff$V1)
