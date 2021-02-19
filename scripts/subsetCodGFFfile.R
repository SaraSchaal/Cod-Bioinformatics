## Subset GFF file

cod.gff <- read.table("src/alignment/GCF_902167405.1_gadMor3.0_genomic.gff", sep = "\t", quote = "")
chrom.gff <- cod.gff[cod.gff$V3 == "region",]
write.table(chrom.gff, "src/alignment/GCF_902167405.1_gadMor3.0_genomic_scaff_contigs.gff", row.names = FALSE, 
            sep = "\t", col.names = FALSE, quote = FALSE)
head(cod.gff)[1:8]
levels(cod.gff$V3)
levels(cod.gff$V2)
levels(cod.gff$V1)


## split linkage group 2 up
chr2.final.loc <- 28732775
chr2.by10k <- seq(from = 0, to = 28732775, by = 10000)
chr2.by10kb.start <- c(chr2.by10k, chr2.final.loc)
chr2.by10kb.end <- chr2.by10kb+9999
section.names <- paste0("NC_044049.1.sec", seq(from = 1, to = length(chr2.by10kb), 
                                               by = 1))

df.chr2.10kb <- as.matrix(cbind(section.names, chr2.by10kb.start, chr2.by10kb.end))
write.table(df.chr2.10kb, "src/alignment/GCF_902167405.1_gadMor3.0_chr2by10kb.gff", row.names = FALSE, 
            sep = "\t", col.names = FALSE, quote = FALSE)

