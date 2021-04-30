snpcov100kb <- read.table("src/freebayes/coverage.txt", header=FALSE, sep= "\t")
snpcov50kb <- read.table("src/freebayes/coverage50kb.txt", header=FALSE, sep= "\t")
snpcov150kb <- read.table("src/freebayes/coverage150kb.txt", header=FALSE, sep= "\t")
colnames(snpcov100kb) <- c("chrom", "windowStart", "windowEnd", "coverage")
colnames(snpcov50kb) <- c("chrom", "windowStart", "windowEnd", "coverage")
colnames(snpcov150kb) <- c("chrom", "windowStart", "windowEnd", "coverage")
sum(snpcov100kb$coverage) # should be 2766167
sum(snpcov50kb$coverage)
sum(snpcov150kb$coverage)
#snpcov100kb$perBase <- snpcov100kb$coverage/snpcov100kb$windowEnd
chromNum <- seq(from = 44048, to = 44070, by = 1)
chromNames <- paste0("NC_0", chromNum, ".1")

library(ggplot2)


ggplot(data = snpcov50kb, aes(x = windowEnd, y = coverage)) +
  geom_point(col = "cornflowerblue", alpha = 0.5) +
  facet_wrap(~chrom) +
  labs(y = "Coverage", x = "location every 50kb bases", 
       title = "Genome Coverage of SNPs every 50kb") +
  ylim(0, max(snpcov50kb$coverage)) +
  xlim(0, max(snpcov50kb$windowEnd)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

ggplot(data = snpcov100kb, aes(x = windowEnd, y = coverage)) +
  geom_point(col = "cornflowerblue", alpha = 0.5) +
  facet_wrap(~chrom) +
  labs(y = "Coverage", x = "location every 100kb bases", 
       title = "Genome Coverage of SNPs every 100kb") +
  ylim(0, max(snpcov100kb$coverage)) +
  xlim(0, max(snpcov100kb$windowEnd)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        legend.position = "none")

ggplot(data = snpcov150kb, aes(x = windowEnd, y = coverage)) +
  geom_point(col = "cornflowerblue", alpha = 0.5) +
  facet_wrap(~chrom) +
  labs(y = "Coverage", x = "location every 150kb bases", 
       title = "Genome Coverage of SNPs every 150 kb") +
  ylim(0, max(snpcov150kb$coverage)) +
  xlim(0, max(snpcov150kb$windowEnd)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        legend.position = "none")
