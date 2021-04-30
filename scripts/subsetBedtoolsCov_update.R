## Process 
setwd("/scratch/schaal.s/CodGenomes")
#install.packages("packrat")
packrat::init("/scratch/schaal.s/CodGenomes/")
setwd("/scratch/schaal.s/CodGenomes/packrat")
Sys.setenv(MAKEFLAGS = "-j4") # if you requested 4 cores, uncomment this line to set your system to use them

library(packrat)
library(dplyr)
library(ggplot2)
library(data.table)
samps <- c("Pop1_16216", "Pop1_17291", "Pop4_17236", "Pop5_17278", "Pop6_18017")
colors <- c("steelblue2", "chartreuse3", "orchid2", "firebrick1", "goldenrod2")
n <- 10000

for(i in 1:length(samps)){
  df <- fread(paste0("bedtools_coverage/", samps[i], ".coverageCalcChr.txt"), 
                                                      sep =" ", quote = "")
  colnames(df) <- c("chrom", "base", "coverage")
  df.sample.data <- NULL
  chroms <- unique(df[,1])
  head(df)
  for(j in 1:length(chroms)){
    df.chrom <- df[df$chrom == chroms[j], ] #see if filter is faster
    #extra <- nrow(df.chrom) - length(rep(1:(nrow(df.chrom)/n), each = n))
    extra <- nrow(df.chrom) %% n
    # group basically is taking our windows set by n and giving each window a number id
    # this grouped dataframe is made by binding the original df.chrom with grouping values for the window size you want 
    # then finding what the last grouping value would be using because it will most like not be an even increment of your n
    # for example if you divide a chromosome by your n and get 3000.3 then you can easily input the first 3000
    # grouping values in with rep (middle part of the following cbind function) then the last grouping value will be 3001 which
    # you get by rounding using ceiling in the last part of this cbind
    grouped <- cbind(df.chrom, c(rep(1:(nrow(df.chrom)/n), each = n), rep(ceiling(nrow(df.chrom)/n), extra))) 
    colnames(grouped)[4] <- "grouping"
    df.sample.data <- rbind(df.sample.data, grouped) # add row
  }
  # Now you are all set to aggregate and get your average coverage per increment set by n
  df.covAve <- aggregate(coverage~grouping + chrom, data = df.sample.data, FUN = mean)

  png(paste0("figures/", samps[i], "MaxcoveragePlot", n, ".pdf"), height= 15, width=15)
   
  ## I make two plots because there will undoubtly be some loci that have really high coverage which makes
  # it hard to see all the spread of the majority of the data. This first graph is set to the max coverage
  # found in the data frame and the second plot is setting your y limit to a more reasonable value for your
  # data. For me 100 x was good but feel free to change to what is appropriate for your data.
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
  
  png(paste0("figures/", samps[i], "100XcoveragePlot", n, ".pdf"), height= 15, width=15)
  
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

