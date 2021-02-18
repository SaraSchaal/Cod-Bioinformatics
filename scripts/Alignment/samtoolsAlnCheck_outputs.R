### Analyzing Alignment Output ###
folder <- "src/alignment/alignCheckFiles/multiqc_data/"
alignResults <- read.csv(paste(folder, "multiqc_samtools_flagstat.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
head(alignResults)

alignResults$megaBases_seq <- (alignResults$paired_in_sequencing_passed * 300) / 10^6
alignResults$megaBases_mapped <- (alignResults$properly_paired_passed * 300) / 10^6
alignResults$coverage_seq <- alignResults$megaBases_seq / 613/2
alignResults$coverage_mapped <- alignResults$megaBases_mapped /613 /2

par(mfrow = c(2,1))
hist(alignResults$coverage_mapped, breaks = 12, main = "Properly Paired Genome Coverage", 
     xlab = "coverage (X)", col = "cornflowerblue")
av.cov <- round(mean(alignResults$coverage_mapped[which(!is.na(alignResults$coverage_mapped))]), 2)
sd.cov <- round(sd(alignResults$coverage_mapped[which(!is.na(alignResults$coverage_mapped))]), 2)
text(25, 60, labels = paste("mean = ", av.cov,"\n sd = ", sd.cov))

hist(alignResults$coverage_seq, breaks = 12, main = "Paired in Sequencing Genome Coverage", 
     xlab = "coverage (X)", col = "navyblue")
av.cov <- round(mean(alignResults$coverage_seq[which(!is.na(alignResults$coverage_seq))]), 2)
sd.cov <- round(sd(alignResults$coverage_seq[which(!is.na(alignResults$coverage_seq))]), 2)
text(25, 60, labels = paste("mean = ", av.cov,"\n sd = ", sd.cov))
