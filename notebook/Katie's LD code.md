Katie's LD code

```
#/bin/bash
​
cd $WD
​
VCF='/scratch/schaal.s/CodGenomes/11_SNPfiltering/merged.f.99ind.MAF05.vcf'
echo $VCF
​
RESDIR='/Users/lotterhos/Google Drive/katie_research/___Eastern Oyster Genome/Population Structure/large_outputs/ld_2020/AtlanticEastern/'
​
echo $RESDIR'geno_ld_window_499950-500000_AtlanticEastern'
​
###
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 499990 --ld-window-bp 500000 --out $RESDIR'geno_ld_window_499950-500000_AtlanticEastern'
tail $RESDIR'geno_ld_window_499950-500000_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_499950-500000_AtlanticEastern.geno.ld'
​
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 99990 --ld-window-bp 100000 --out $RESDIR'geno_ld_window_99990-100000_AtlanticEastern'
tail $RESDIR'geno_ld_window_99990-100000_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_99990-100000_AtlanticEastern.geno.ld'
​
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 49990 --ld-window-bp 50000 --out $RESDIR'geno_ld_window_49990-50000_AtlanticEastern'
tail $RESDIR'geno_ld_window_49990-50000_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_49990-50000_AtlanticEastern.geno.ld'
​
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 9990 --ld-window-bp 10000 --out $RESDIR'geno_ld_window_9990-10000_AtlanticEastern'
tail $RESDIR'geno_ld_window_9990-10000_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_9990-10000_AtlanticEastern.geno.ld'
​
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 4990 --ld-window-bp 5000 --out $RESDIR'geno_ld_window_4990-5000_AtlanticEastern'
tail $RESDIR'geno_ld_window_4990-5000_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_4990-5000_AtlanticEastern.geno.ld'
​
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 490 --ld-window-bp 500 --out $RESDIR'geno_ld_window_490-500_AtlanticEastern'
tail $RESDIR'geno_ld_window_490-500_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_490-500_AtlanticEastern.geno.ld'
​
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 190 --ld-window-bp 200 --out $RESDIR'geno_ld_window_190-200_AtlanticEastern'
tail $RESDIR'geno_ld_window_190-200_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_190-200_AtlanticEastern.geno.ld'
​
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 50 --ld-window-bp 50 --out $RESDIR'geno_ld_window_50-50_AtlanticEastern'
tail $RESDIR'geno_ld_window_50-50_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_50-50_AtlanticEastern.geno.ld'
​
​
​
## Capture inflection points
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 990 --ld-window-bp 1000 --out $RESDIR'geno_ld_window_990-1000_AtlanticEastern'
tail $RESDIR'geno_ld_window_990-1000_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_990-1000_AtlanticEastern.geno.ld'
​
vcftools --gzvcf $VCF --keep AtlanticEastern.pop --geno-r2 --ld-window-bp-min 2490 --ld-window-bp 2500 --out $RESDIR'geno_ld_window_2490-2500_AtlanticEastern'
tail $RESDIR'geno_ld_window_2490-2500_AtlanticEastern.geno.ld'
wc -l $RESDIR'geno_ld_window_2490-2500_AtlanticEastern.geno.ld'

```


R script
```
###############################################################################
#
# File      : LDAnalysisScript.R 
# History   : 10/21/2018  Created by K Bodie Weedop (KBW)
#
###############################################################################
​
###############################################################################
#
# This script works with the various LD analysis files given by VCFtools when
# running --geno-r2 on the original VCF file
# (Combined.SNP.TRSdp5g1FnDNAmaf052alleles.vcf.gz). 
#
###############################################################################
​
​
load.ld.data <- function (path = NULL) {
    if (!dir.exists(path)) {
        stop("ERROR: The path that you have provided is not a directory")
    } else {
        ld.files <- paste(path, list.files(path), sep = "")
    }
​
    ld.data.50bp <- read.csv(ld.files[which(grepl("_50-50", ld.files, fixed=TRUE) == TRUE)], 
                             sep="\t",
                             header=TRUE,
                             stringsAsFactors=FALSE)
​
    ld.data.200.250bp <- read.csv(ld.files[which(grepl("_200-250", ld.files, fixed=TRUE) == TRUE)], 
                                  sep="\t",
                                  header=TRUE,
                                  stringsAsFactors=FALSE)
​
    ld.data.450.500bp <- read.csv(ld.files[which(grepl("_450-500", ld.files, fixed=TRUE) == TRUE)], 
                                  sep="\t",
                                  header=TRUE,
                                  stringsAsFactors=FALSE)
​
    ld.data.4500.5000bp <- read.csv(ld.files[which(grepl("_4500-5000", ld.files, fixed=TRUE) == TRUE)], 
                                    sep="\t",
                                    header=TRUE,
                                    stringsAsFactors=FALSE)
​
    ld.data.9500.10kbp <- read.csv(ld.files[which(grepl("_9500-10000", ld.files, fixed=TRUE) == TRUE)], 
                                    sep="\t",
                                    header=TRUE,
                                    stringsAsFactors=FALSE)
​
    ld.data.49500.50kbp <- read.csv(ld.files[which(grepl("_49500-50000", ld.files, fixed=TRUE) == TRUE)], 
                                    sep="\t",
                                    header=TRUE,
                                    stringsAsFactors=FALSE)
​
    ld.data.99500.100kbp <- read.csv(ld.files[which(grepl("_99500-100000", ld.files, fixed=TRUE) == TRUE)], 
                                     sep="\t",
                                     header=TRUE,
                                     stringsAsFactors=FALSE)
​
    ld.data.499500.500kbp <- read.csv(ld.files[which(grepl("500000", ld.files, fixed=TRUE) == TRUE)], 
                                      sep="\t",
                                      header=TRUE,
                                      stringsAsFactors=FALSE)
​
    ld.vars <- list("50" = ld.data.50bp,
                    "200" = ld.data.200.250bp, 
                    "450" = ld.data.450.500bp, 
                    "4500" = ld.data.4500.5000bp,
                    "9500" = ld.data.9500.10kbp,
                    "49500" = ld.data.49500.50kbp,
                    "99500" = ld.data.99500.100kbp,
                    "499500" = ld.data.499500.500kbp)
    return(ld.vars)
}
​
# get the LD decay plots for all of the chromosomes in the ld.datasets
​
decay.plots <- function (path = NULL, file.suffix = NULL, plot.title = NULL) {
    ld.datasets <- load.ld.data(path = path)
​
    CHR <- c("NC_035780.1", "NC_035781.1", "NC_035782.1", "NC_035783.1", 
             "NC_035784.1", "NC_035785.1", "NC_035786.1", "NC_035787.1", 
             "NC_035788.1", "NC_035789.1")
    
    col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                  '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                  '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
                  '#ffffff', '#000000')
    
    chr.num <- length(unique(ld.datasets[[1]]$CHR))
    for (i in length(ld.datasets)) {
        if (length(unique(ld.datasets[[i]]$CHR)) != chr.num) {
            print(names(ld.datasets[i]))
            #stop("The number of chromosomes in the datasets do not match")
        } else {
​
        }
    }
    
    chr.index <- rep(NA, chr.num)
    for (i in 1:chr.num){
        chr.index[i] <- substr(unique(ld.datasets[[1]]$CHR)[i], start = 9, stop = 9)
    }
    chr.index <- sort(chr.index)
​
    ld.stats <- data.frame(matrix(NA, length(ld.datasets), chr.num))
    ld.error <- data.frame(matrix(NA, length(ld.datasets), chr.num))
    ld.stats[,1] <- as.numeric(names(ld.datasets))
    ld.error[,1] <- as.numeric(names(ld.datasets))
    for (i in 1:length(ld.datasets)) {
        ld.datasets[[i]]$R.2 <- as.numeric(ld.datasets[[i]]$R.2)
        for (j in 1:chr.num) {
            ld.stats[i, j+1] <- summary(ld.datasets[[i]]$R.2[which(ld.datasets[[i]]$CHR == paste("NC_03578", chr.index[j], ".1", sep=""))])[4]
            ld.error[i, j+1] <- sd(na.omit(ld.datasets[[i]]$R.2[which(ld.datasets[[i]]$CHR == paste("NC_03578", chr.index[j], ".1", sep=""))])) / sqrt(length(na.omit(ld.datasets[[i]]$R.2[which(ld.datasets[[i]]$CHR == paste("NC_03578", chr.index[j], ".1", sep=""))])))
        }
    }
​
    png(paste("figures/1LD_analysis/ldDecayPlot_", file.suffix, ".png", sep=""),
        width = 8,
        height = 8,
        units = "in",
        res = 300)
    par(bty="l")
    for (i in 2:ncol(ld.stats)) {
        if (i == 2){
            plot(jitter(ld.stats[,i]) ~ ld.stats[, 1],
                 pch = i,
                 col = col_vector[i],
                 cex = 2,
                 xlab="Distance (bp)", 
                 ylab="LD (r^2) +/- Std. Error",
                 ylim=c(0, 0.24),
                 type = "b",
                 log = "x",
                 main=paste("LD Decay:", plot.title, sep=" "))
                 arrows(ld.stats[,1], ld.stats[,i]-ld.error[,i], ld.stats[,1], ld.stats[,i]+ld.error[,i], col=col_vector[i], length=0.05, angle=90, code=3)
        } else {
            points(jitter(ld.stats[,i]) ~ ld.stats[, 1],
                   pch = i,
                   col = col_vector[i],
                   cex = 2,
                   type = "b")
            arrows(ld.stats[,1], ld.stats[,i]-ld.error[,i], ld.stats[,1], ld.stats[,i]+ld.error[,i], col=col_vector[i], length=0.05, angle=90, code=3)
        }
        legend("topright", legend = CHR, pch = 2:(length(CHR)+1), col = col_vector[2:(length(CHR)+1)])
    }
    dev.off()
}
​
# Get plot for all populations
#decay.plots(path="data/large_data/ldAnalysisData/allpops/",
#            file.suffix = "allpops",
#            plot.title = "All Populations")
#
## # Get plot for all populations when excluding LM
#decay.plots(path="data/large_data/ldAnalysisData/excluding_LM/",
#            file.suffix = "excluding_LM",
#            plot.title = "Excluding LM")
## 
## # Get plot for wild Atlantic populations
#decay.plots(path="data/large_data/ldAnalysisData/excluding_selection/",
#            file.suffix = "excluding_selection",
#            plot.title = "Atlantic Wild Populations")
#
## Get plot for selection Atlantic populations
# decay.plots(path="data/large_data/ldAnalysisData/excluding_wild/",
#             file.suffix = "excluding_wild",
#             plot.title = "Atlantic Selection Populations")
​
# Get plot for selection Atlantic populations
 decay.plots(path="data/large_data/ldAnalysisData/unrelated/",
             file.suffix = "unrelated",
             plot.title = "Unrelated Individuals")

```