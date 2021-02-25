## Process 
setwd("/scratch/schaal.s/CodGenomes")
library(ggplot2)
samps <- c("Pop1_16216", "Pop1_17291", "Pop4_17236", "Pop5_17278", "Pop6_18017")
colors <- c("steelblue2", "chartreuse3", "orchid2", "firebrick1", "goldenrod2")
n <- 10000
full.data <- NULL
for(i in 1:length(samps)){
  df <- read.table(paste0("bedtools_coverage/", samps[i], ".coverageCalcChr.txt"), 
                                                      sep =" ", quote = "")
  colnames(df) <- c("chrom", "base", "coverage")
  df.sample.data <- NULL
  for(j in 1:length(unique(df[,1]))){
    chroms <- unique(df[,1])
    df.chrom <- df[df$chrom == chroms[j], ]
    print(n)
    extra <- nrow(df.chrom) - length(rep(1:(nrow(df.chrom)/n), each = n))
    grouped <- cbind(df.chrom, c(rep(1:(nrow(df.chrom)/n), each = n), rep(ceiling(nrow(df.chrom)/n), extra)))
    colnames(grouped)[4] <- "grouping"
    df.sample.data <- rbind(df.sample.data, grouped)
  }
  df.covAve <- aggregate(coverage~grouping + chrom, data = df.sample.data, FUN = mean)
  print(head(df.covAve))
  print(dim(df.covAve))
 
  pdf(paste0("figures/", samps[i], "MaxcoveragePlot", n, ".pdf"), height= 15, width=15)
   
  print(ggplot(data = df.covAve, aes(x = grouping, y = coverage)) +
    geom_point(col = colors[i], alpha = 0.5) +
    facet_wrap(~chrom) +
    labs(y = "Coverage", x = paste0("location every ", n, " bases"), 
         title = paste0(samps[i], "Genome Coverage up to Max Coverage")) +
    ylim(0, max(df.covAve$coverage)) +
    xlim(0, max(df.covAve$grouping)) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), 
          legend.position = "none"))
  dev.off()
  
  pdf(paste0("figures/", samps[i], "100XcoveragePlot", n, ".pdf"), height= 15, width=15)
  
  print(ggplot(data = df.covAve, aes(x = grouping, y = coverage)) +
          geom_point(col = colors[i], alpha = 0.5) +
          facet_wrap(~chrom) +
          labs(y = "Coverage", x = paste0("location every ", n, " bases"),
               title = paste0(samps[i], "Genome Coverage up to 100X")) +
          ylim(0, 100) +
          xlim(0, max(df.covAve$grouping)) +
          theme_bw() + 
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"), 
                legend.position = "none"))
  dev.off()
  
  
}

