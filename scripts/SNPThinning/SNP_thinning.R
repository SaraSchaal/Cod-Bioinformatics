#### SNP Thinning for LD & Identifying Outliers ####
## Sara M. Schaal

  packages_needed <-  c("bigsnpr", "bigstatsr", "vcfR", "xgboost", "dplyr", 
                        "ggplot2", "VGAM", "hexbin", "viridisLite", "pcadapt")

  ## Install packages that aren't installed already
  for (i in 1:length(packages_needed)){
    if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i], repos = "http://cran.us.r-project.org")}
  }

  ## Load each library
  for (i in 1:length(packages_needed)){
    library( packages_needed[i], character.only = TRUE)
  }
  
  ## Plink download
  plink2 <- download_plink2(AVX2 = FALSE)

  ## Input and Output folders
  folderIn <- "./src/bedfiles/"
  folderOut <- "./"

  ## filter for missing data (only need to run this once to make the files)
  # chrom 1
  #snp_plinkQC(plink2, paste0(folderIn, "bedfiles2/VarCall_NC_044048.1.f"),file.type = "--bfile", 
            #maf = 0.05, geno = 0, "--allow-extra-chr", prefix.out = paste0(folderIn, "bedfiles2/VarCall_NC_044048.1.f_QC"))
  
  # all data
  #snp_plinkQC(plink2, paste0(folderIn, "merged.f.99ind.MAF05"),file.type = "--bfile", 
              #maf = 0.05, geno = 0, "--allow-extra-chr", prefix.out = paste0(folderIn, "merged.f.99ind.MAF05_QC"))
  ## one sample dropped out due to missing data
  
  ## Reads in bed file made with plink above 
  bed.file.filt <- snp_readBed(paste0(folderIn, "merged.f.99ind.MAF05_QC.bed"))
  bigSNP <- snp_attach(bed.file.filt)

  ## pull out genotypes into a matrix and chromosomes into a vector
  G_full <- bigSNP$genotypes
  pos_full <- bigSNP$map$physical.pos
  G_full[1:100]
  n <- nrow(G_full)
  m <- ncol(G_full)
  table(G_full[1:m])
  sum(is.na(G_full[1:m])) # should have no NAs 

  # Store just G matrix
  G_mat <- G_full[1:n, 1:m]
  colnames(G_mat) <- bigSNP$map$physical.pos
  rownames(G_mat) <- bigSNP$fam$family.ID

  G_coded <- add_code256(big_copy(bigSNP$genotypes,
                                    type = "raw"),
                            code=bigsnpr:::CODE_012)

  # get chromosome values
  chrom_full <- bigSNP$map$chromosome
  chromNumbers <- seq(from = 44048, to = 44070, by = 1) 
  chromNames <- data.frame(chromID = paste0("NC_0", chromNumbers, ".1"), chromNum = 1:23)
  for(j in 1:nrow(chromNames)){
    chrom_full[chrom_full == chromNames$chromID[j]] <- j
  }
  chrom_full <- as.numeric(chrom_full)


  # remove SNPs that are in long range LD
  lrLD <- snp_autoSVD(G=G_full, infos.chr = chrom_full,
                        infos.pos = pos_full, size = 1)
      # set the window size increase window size to get closer to 200K
  str(lrLD)
  lrLD.ind <- attr(lrLD, "subset") # this contains new set of pruned SNPs
  length(lrLD.ind) 
  length(G_full[lrLD.ind])
  G_thin <- G_full[,lrLD.ind]
  colnames(G_thin) <- bigSNP$map$physical.pos[lrLD.ind]
  rownames(G_thin) <- bigSNP$fam$family.ID
  lrLD.tb <- attr(lrLD, "lrldr")
  sample(lrLD.ind)

  #df.pcs.full$pop <- factor(df.pcs.full$pop, levels = pop.names)

######################################################################################################
#### subset SNPs to 200K ####

  set.seed(48)
  subset200K <- sample(lrLD.ind, size = 200000)
  set.seed(NULL)
  #head(sample(lrLD.ind, size = 200000))
  subset200Korder <- subset200K[order(subset200K)]
  G_sub <- G_full[,subset200Korder]
  chrom_sub <- chrom_full[subset200Korder]
  pos_sub <- pos_full[subset200Korder]
  
  random200K <- list(G = G_sub, 
                     Pos = pos_sub,
                     Chr = chrom_sub,
                     Pop.ID = bigSNP$fam$family.ID,
                     Sample.ID = bigSNP$fam$sample.ID)
  Gsub_coded <- add_code256(big_copy(random200K$G,
                                     type = "raw"),
                            code=bigsnpr:::CODE_012)
  saveRDS(random200K, paste(folderOut, "random200KthinSNPs.rds"))
  
  lrLD_sub <- snp_autoSVD(G=Gsub_coded, infos.chr = random200K$Chr,
                      infos.pos = random200K$Pos, size = 1)
  lrLD_sub.tb <- attr(lrLD_sub, "lrldr")
  rownames(lrLD_sub$u) <- bigSNP$fam$sample.ID
  par(mar=c(4,4,1,1))
  plot(lrLD_sub$d, xlab="PC", ylab= "Variation explained", bty="l")
  
  df.pcs.sub <- data.frame(lrLD_sub$u, pop = bigSNP$fam$family.ID)
  colnames(df.pcs.sub) <- c(paste0("PC", 1:10), "pop")
  pop.names <- c("Mass.Winter", "Mass.Red", "Mass.Spring", "Cashes.Olive", "Cashes.Red", 
                 "Ice.SWOff", "Ice.SWNear", "Ice.NEOff", "Ice.NENear")
  for(i in 1:length(pop.names)){
    #df.pcs.full$pop[df.pcs.full$pop == paste0("Pop", i)] <- pop.names[i]
    #df.pcs.chrom1$pop[df.pcs.chrom1$pop == paste0("Pop", i)] <- pop.names[i]
    df.pcs.sub$pop[df.pcs.sub$pop == paste0("Pop", i)] <- pop.names[i]
  }
  df.pcs.sub$pop <- factor(df.pcs.sub$pop, levels = pop.names)
  
  ggplot(df.pcs.sub, aes(x = PC1, y = PC2, color = pop)) + 
    geom_point(alpha = 0.75) +
    scale_color_manual(name = "Population", 
                       values = viridis(length(pop.names))) + 
    labs(title = "All Samples Full Genome") +
    theme_bw()
  
  
  G_sub_GOM <- G_full[1:137,subset200Korder]
  chrom_sub_GOM <- chrom_full[subset200Korder]
  pos_sub_GOM <- pos_full[subset200Korder]
  
  random200K_GOM <- list(G = G_sub_GOM, 
                     Pos = pos_sub_GOM,
                     Chr = chrom_sub_GOM,
                     Pop.ID = bigSNP$fam$family.ID[1:137],
                     Sample.ID = bigSNP$fam$sample.ID[1:137])
  GsubGOM_coded <- add_code256(big_copy(random200K_GOM$G,
                                     type = "raw"),
                            code=bigsnpr:::CODE_012)
  saveRDS(random200K_GOM, paste(folderOut, "random200KthinSNPs_GOM.rds"))
  
  lrLD_sub_GOM <- snp_autoSVD(G=GsubGOM_coded, infos.chr = random200K_GOM$Chr,
                          infos.pos = random200K_GOM$Pos, size = 1)
  
  df.pcs.sub_GOM <- data.frame(lrLD_sub_GOM$u, pop = random200K_GOM$Pop.ID)
  colnames(df.pcs.sub_GOM) <- c(paste0("PC", 1:10), "pop")
  pop.names <- c("Mass.Winter", "Mass.Red", "Mass.Spring", "Cashes.Olive", "Cashes.Red", 
                 "Ice.SWOff", "Ice.SWNear", "Ice.NEOff", "Ice.NENear")
  for(i in 1:length(pop.names)){
    df.pcs.sub_GOM$pop[df.pcs.sub_GOM$pop == paste0("Pop", i)] <- pop.names[i]
  }
  df.pcs.sub_GOM$pop <- factor(df.pcs.sub_GOM$pop, levels = pop.names)
  
  ggplot(df.pcs.sub_GOM, aes(x = PC1, y = PC2, color = pop)) + 
    geom_point(alpha = 0.75) +
    scale_color_manual(name = "Population", 
                       values = viridis(length(pop.names))) + 
    labs(title = "All Samples Full Genome") +
    theme_bw()
  
  
  G_sub_ICE <- G_full[138:295,subset200Korder]
  chrom_sub_ICE <- chrom_full[subset200Korder]
  pos_sub_ICE <- pos_full[subset200Korder]
  
  random200K_ICE <- list(G = G_sub_ICE, 
                         Pos = pos_sub_ICE,
                         Chr = chrom_sub_ICE,
                         Pop.ID = bigSNP$fam$family.ID[138:295],
                         Sample.ID = bigSNP$fam$sample.ID[138:295])
  GsubICE_coded <- add_code256(big_copy(random200K_ICE$G,
                                        type = "raw"),
                               code=bigsnpr:::CODE_012)
  saveRDS(random200K_ICE, paste(folderOut, "random200KthinSNPs_ICE.rds"))
  
  lrLD_sub_ICE <- snp_autoSVD(G=GsubICE_coded, infos.chr = random200K_ICE$Chr,
                              infos.pos = random200K_ICE$Pos, size = 1)
  
  
  df.pcs.sub_ICE <- data.frame(lrLD_sub_ICE$u, pop = random200K_ICE$Pop.ID)
  colnames(df.pcs.sub_ICE) <- c(paste0("PC", 1:10), "pop")
  pop.names <- c("Mass.Winter", "Mass.Red", "Mass.Spring", "Cashes.Olive", "Cashes.Red", 
                 "Ice.SWOff", "Ice.SWNear", "Ice.NEOff", "Ice.NENear")
  for(i in 1:length(pop.names)){
    df.pcs.sub_ICE$pop[df.pcs.sub_ICE$pop == paste0("Pop", i)] <- pop.names[i]
  }
  df.pcs.sub_ICE$pop <- factor(df.pcs.sub_ICE$pop, levels = pop.names)
  
  ggplot(df.pcs.sub_ICE, aes(x = PC1, y = PC2, color = pop)) + 
    geom_point(alpha = 0.75) +
    scale_color_manual(name = "Population", 
                       values = viridis(length(pop.names))[6:9]) + 
    labs(title = "All Samples Full Genome") +
    theme_bw()
  

#### end subset 200K thinned snps
######################################################################################################    


######################################################################################################    
#### start outlier detection

## pcadapt
#[1:nrow(G_coded), 1:ncol(G_coded)]
#dist <- snp_pcadapt(G_thin, U.row = df.pcs.sub[, 1:10])

# global outliers analysis
  # full dataset run  pcadapt
    df.global.outliers <- read.pcadapt(t(G_mat))
    df.global <- as.data.frame(pos_full)
    colnames(df.global) <- "position"
    global.outliers <- pcadapt(df.global.outliers, K = 5)
    df.global$pca_ALL_PC1 <- global.outliers$loadings[,1]
    df.global$pca_ALL_PC2 <- global.outliers$loadings[,2]
    plot(global.outliers, option = "screeplot")
  
    folderOut <- "./Figures/Outliers/"
    png(paste0(folderOut, "pcadapt_globalOutliers.png"), width = 15, height = 7, units = 'in', res = 300)
      plot(global.outliers)
    dev.off()
  
  # thinned SNPs run pcadapt
    df.global.outliers.thin <- read.pcadapt(t(G_thin))
    global.outliers.thin.pcadapt <- pcadapt(df.global.outliers.thin)
    df.global$pca_PRUNED_PC1 <- NA 
    df.global$pca_PRUNED_PC2 <- NA
    df.global$pca_PRUNED_PC1[lrLD.ind] <- global.outliers.thin.pcadapt$loadings[,1]
    df.global$pca_PRUNED_PC2[lrLD.ind] <- global.outliers.thin.pcadapt$loadings[,2] 
    
    
  # outlier stats 
    # all data
    df.global$pcadapt_ALL_chisq <- as.numeric(global.outliers$chi2.stat)
    df.global$pcadapt_ALL_log10p <- -log10(global.outliers$pvalues)
    
    # thinned SNPs
    outliers <- snp_gc(snp_pcadapt(G_coded, U.row = lrLD$u[,1]))
    df.global$pcadapt_PRUNED_log10p <- -predict(outliers,log10=T)
    df.global$pcadapt_PRUNED_pvalue <- 10^(-df.global$pcadapt_PRUNED_log10p)
    df.global$qvalue <- qvalue(df.global$pcadapt_PRUNED_pvalue)$qvalues  
    df.global$pcadapt_outlier <- ifelse(df.global$qvalue > 0.01, FALSE, TRUE)
    df.global$pcadapt_outlier <- as.factor(df.global$pcadapt_outlier)
    
    
    # plot
    max_value_log10 <- max(!is.na(df.global$pcadapt_PRUNED_log10p))
    pcadapt.log10p <- ggplot(df.global, aes(x = position, y = pcadapt_PRUNED_log10p))  +
    
      geom_point(data = df.global[!is.na(df.global$pcadapt_outlier),], aes(color = pcadapt_outlier, shape = pcadapt_outlier))+
      scale_color_manual(values = c("black", "red")) +
      scale_shape_manual(values = c(19, 1)) +
      theme_classic() +
      theme(panel.background = element_blank(), 
            strip.background = element_rect(colour = "white", fill = "grey92"),
            text = element_text(size = 11)) +
      labs(title = "OutFLANK",
           y = "-log10(p-values)",
           x = " ") + 
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
    
    
    
#pop.names <- c("Mass.Winter", "Mass.Red", "Mass.Spring", "Cashes.Olive", "Cashes.Red", 
               #"Ice.SWOff", "Ice.SWNear", "Ice.NEOff", "Ice.NENear")

# Iceland cod outlier analysis
  G_mat_ice <- subset(G_mat,  subset = rownames(G_mat) == "Pop6" | rownames(G_mat) == "Pop7" |
                        rownames(G_mat) == "Pop8" | rownames(G_mat) == "Pop9")
  G_thin_ice <- subset(G_thin,  subset = rownames(G_thin) == "Pop6" | rownames(G_thin) == "Pop7" |
                       rownames(G_thin) == "Pop8" | rownames(G_thin) == "Pop9")
  df.ice.outliers <- read.pcadapt(t(G_thin_ice))
  ice.outliers <- pcadapt(df.ice.outliers, K = 2)
  plot(ice.outliers, option = "screeplot")
  png(paste0(folderOut, "pcadapt_iceOutliers.png"), width = 15, height = 7, units = 'in', res = 300)
    plot(ice.outliers)
  dev.off()
  
# outliers GOM offshore/nearshore
  G_mat_GOM <- subset(G_mat,  subset = rownames(G_mat) == "Pop1" | rownames(G_mat) == "Pop2" |
                         rownames(G_mat) == "Pop3" | rownames(G_mat) == "Pop4"| rownames(G_mat) == "Pop5")
  G_thin_GOM <- subset(G_thin,  subset = rownames(G_thin) == "Pop1" | rownames(G_thin) == "Pop2" |
                       rownames(G_thin) == "Pop3" | rownames(G_thin) == "Pop4"| rownames(G_thin) == "Pop5")
  df.GOM.outliers <- read.pcadapt(t(G_thin_GOM))
  GOM.outliers <- pcadapt(df.GOM.outliers, K = 2)
  plot(GOM.outliers, option = "screeplot")
  png(paste0(folderOut, "pcadapt_GOMOutliers.png"), width = 15, height = 7, units = 'in', res = 300)
    plot(GOM.outliers)
  dev.off()

# outliers iceland offshore/nearshore
G_thin_coast_GOM <- subset(G_thin,  subset = rownames(G_thin) == "Pop1" | rownames(G_thin) == "Pop2" |
                       rownames(G_thin) == "Pop3" | rownames(G_thin) == "Pop4"| rownames(G_thin) == "Pop5" | 
                         rownames(G_thin) == "Pop7"|rownames(G_thin) == "Pop9")
df.GOM.outliers <- pcadapt(df.GOM.outliers)
plot(GOM.outliers, option = "screeplot")

## outflank
  # ensure there are no NAs
  sum(is.na(G_thin)) # checks out; no NAs

# global outliers 
  FSTs_full <- MakeDiploidFSTMat(G_mat, locusNames = colnames(G_mat), popNames = rownames(G_mat))
  FSTs <- MakeDiploidFSTMat(G_thin, locusNames = colnames(G_thin), popNames = rownames(G_thin))
  global_out <- OutFLANK(FSTs, NumberOfSamples=295, qthreshold = 0.05, Hmin = 0.1)
  str(global_out)
  png(paste0(folderOut, "outflank_globalOutliers_FstFreq.png"), width = 15, height = 7, units = 'in', res = 300)
    OutFLANKResultsPlotter(global_out, withOutliers = TRUE,
                           NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           FALSE, RightZoomFraction = 0.05, titletext = NULL)
  dev.off()
  df.global.outliers <- pOutlierFinderChiSqNoCorr(FSTs_full, Fstbar = global_out$FSTNoCorrbar, 
                            dfInferred = global_out$dfInferred, qthreshold = 0.05, Hmin=0.1)
  df.global.outliers.thin <- pOutlierFinderChiSqNoCorr(FSTs, Fstbar = global_out$FSTNoCorrbar, 
                                                  dfInferred = global_out$dfInferred, qthreshold = 0.05, Hmin=0.1)
  global.outliers <- df.global.outliers$OutlierFlag==TRUE
  global.outliers.thin <- df.global.outliers.thin$OutlierFlag==TRUE
  png(paste0(folderOut, "outflank_globalOutliers.png"), width = 15, height = 7, units = 'in', res = 300)
    plot(df.global.outliers$LocusName[df.global.outliers$He>0.1], df.global.outliers$FST[df.global.outliers$He>0.1],
         xlab="Position", ylab="FST", col=rgb(0,0,0,0.2), main = "Global Outliers")
    points(df.global.outliers$LocusName[global.outliers], df.global.outliers$FST[global.outliers], col="red", pch=20)  
  dev.off()
  
  png(paste0(folderOut, "outflank_globalOutliersThin.png"), width = 15, height = 7, units = 'in', res = 300)
    plot(df.global.outliers.thin$LocusName[df.global.outliers.thin$He>0.1], df.global.outliers.thin$FST[df.global.outliers.thin$He>0.1],
         xlab="Position", ylab="FST", col=rgb(0,0,0,0.2), main = "Global Outliers")
    points(df.global.outliers.thin$LocusName[global.outliers.thin], df.global.outliers.thin$FST[global.outliers.thin], col="red", pch=20)  
  dev.off()
  
# outliers iceland offshore/nearshore
  FSTs_mat_ice <- MakeDiploidFSTMat(G_mat_ice, locusNames = colnames(G_mat_ice), popNames = rownames(G_mat_ice))
  FSTs_ice <- MakeDiploidFSTMat(G_thin_ice, locusNames = colnames(G_thin_ice), popNames = rownames(G_thin_ice))
  ice_out <- OutFLANK(FSTs_ice, NumberOfSamples=158, qthreshold = 0.05, Hmin = 0.1)
  str(ice_out)
  
  png(paste0(folderOut, "outflank_iceOutliers_FstFreq.png"), width = 15, height = 7, units = 'in', res = 300)
    OutFLANKResultsPlotter(ice_out, withOutliers = TRUE,
                           NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           FALSE, RightZoomFraction = 0.05, titletext = NULL)
  dev.off()
  
  df.ice.outliers <- pOutlierFinderChiSqNoCorr(FSTs_mat_ice, Fstbar = ice_out$FSTNoCorrbar, 
                                               dfInferred = ice_out$dfInferred, qthreshold = 0.05, Hmin=0.1)
  df.ice.outliers.thin <- pOutlierFinderChiSqNoCorr(FSTs_ice, Fstbar = ice_out$FSTNoCorrbar, 
                                                    dfInferred = ice_out$dfInferred, qthreshold = 0.05, Hmin=0.1)
  ice.outliers <- df.ice.outliers$OutlierFlag==TRUE
  ice.outliers.thin <- df.ice.outliers.thin$OutlierFlag==TRUE
  png(paste0(folderOut, "outflank_iceOutliers.png"), width = 15, height = 7, units = 'in', res = 300)
    plot(df.ice.outliers$LocusName[df.ice.outliers$He>0.1], df.ice.outliers$FST[df.ice.outliers$He>0.1],
         xlab="Position", ylab="FST", col=rgb(0,0,0,0.2), main = "Icelandic cod outliers")
    points(df.ice.outliers$LocusName[ice.outliers], df.ice.outliers$FST[ice.outliers], col="red", pch=20)  
  dev.off()
  
  png(paste0(folderOut, "outflank_iceOutliersThin.png"), width = 15, height = 7, units = 'in', res = 300)
    plot(df.ice.outliers.thin$LocusName[df.ice.outliers.thin$He>0.1], df.ice.outliers.thin$FST[df.ice.outliers.thin$He>0.1],
        xlab="Position", ylab="FST", col=rgb(0,0,0,0.2), main = "Icelandic cod outliers")
    points(df.ice.outliers.thin$LocusName[ice.outliers.thin], df.ice.outliers.thin$FST[ice.outliers.thin], col="red", pch=20)  
  dev.off()
  
# outliers GOM offshore/nearshore
  FSTs_mat_GOM <- MakeDiploidFSTMat(G_mat_GOM, locusNames = colnames(G_mat_GOM), popNames = rownames(G_mat_GOM))
  FSTs_GOM <- MakeDiploidFSTMat(G_thin_GOM, locusNames = colnames(G_thin_GOM), popNames = rownames(G_thin_GOM))
  GOM_out <- OutFLANK(FSTs_GOM, NumberOfSamples=137, qthreshold = 0.05, Hmin = 0.1)
  str(GOM_out)
  
  png(paste0(folderOut, "outflank_GOMOutliers_FstFreq.png"), width = 15, height = 7, units = 'in', res = 300)
    OutFLANKResultsPlotter(GOM_out, withOutliers = TRUE,
                           NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           FALSE, RightZoomFraction = 0.05, titletext = NULL)
  dev.off()
  
  df.GOM.outliers <- pOutlierFinderChiSqNoCorr(FSTs_mat_GOM, Fstbar = GOM_out$FSTNoCorrbar, 
                                               dfInferred = GOM_out$dfInferred, qthreshold = 0.05, Hmin=0.1)
  df.GOM.outliers.thin <- pOutlierFinderChiSqNoCorr(FSTs_GOM, Fstbar = GOM_out$FSTNoCorrbar, 
                                                    dfInferred = GOM_out$dfInferred, qthreshold = 0.05, Hmin=0.1)
  GOM.outliers <- df.GOM.outliers$OutlierFlag==TRUE
  GOM.outliers.thin <- df.GOM.outliers.thin$OutlierFlag==TRUE
  png(paste0(folderOut, "outflank_GOMOutliers.png"), width = 15, height = 7, units = 'in', res = 300)
    plot(df.GOM.outliers$LocusName[df.GOM.outliers$He>0.1], df.GOM.outliers$FST[df.GOM.outliers$He>0.1],
         xlab="Position", ylab="FST", col=rgb(0,0,0,0.2), main = "Gulf of Maine Outliers")
    points(df.GOM.outliers$LocusName[GOM.outliers], df.GOM.outliers$FST[GOM.outliers], col="red", pch=20)  
  dev.off()
  
  png(paste0(folderOut, "outflank_GOMOutliersThin.png"), width = 15, height = 7, units = 'in', res = 300)
    plot(df.GOM.outliers.thin$LocusName[df.GOM.outliers.thin$He>0.1], df.GOM.outliers.thin$FST[df.GOM.outliers.thin$He>0.1],
        xlab="Position", ylab="FST", col=rgb(0,0,0,0.2), main = "Gulf of Maine Outliers")
    points(df.GOM.outliers.thin$LocusName[GOM.outliers.thin], df.GOM.outliers.thin$FST[GOM.outliers.thin], col="red", pch=20)  
  dev.off()
  
  
#### end outlier detection
######################################################################################################    

######################################################################################################
## start Chrom 1
bed.file.filt.chrom1 <- snp_readBed(paste0(folderIn,"bedfiles2/VarCall_NC_044048.1.f_QC.bed"))
bigSNP.chrom1 <- snp_attach(bed.file.filt.chrom1)
G_chrom1 <- bigSNP.chrom1$genotypes
pos_chrom1 <- bigSNP.chrom1$map$physical.pos

G_chrom1[1:100]
n <- nrow(G_chrom1)
m <- ncol(G_chrom1)
table(G_chrom1[1:m])
sum(is.na(G_chrom1[1:m]))

chrom_1 <- bigSNP.chrom1$map$chromosome
for(j in 1:nrow(chromNames)){
  chrom_1[chrom_1 == chromNames$chromID[j]] <- j
}
chrom_1 <- as.numeric(chrom_1)

lrLD.chrom1 <- snp_autoSVD(G=G_chrom1, infos.chr = chrom_1,
                           infos.pos = pos_chrom1, size = 1)

str(lrLD.chrom1)
attr(lrLD.chrom1, "subset") # this contains new set of pruned SNPs
lrLD.ind.chrom1 <- attr(lrLD.chrom1, "subset") # this contains new set of pruned SNPs

length(lrLD.ind.chrom1) 
length(G_chrom1[lrLD.ind.chrom1])
lrLD.tb.chrom1 <- attr(lrLD.chrom1, "lrldr")

## Output lrLD regions
if(nrow(lrLD.tb.chrom1) != 0){
  write.table(folderOut, lrLD.tb.chrom1, row.names = FALSE,  col.names = FALSE)
}

#loadings pc1
plot(lrLD.ind.chrom1, lrLD.chrom1$v[,1], col = rgb(0,0,0,0.1))
#loadings pc2
plot(lrLD.ind.chrom1, lrLD.chrom1$v[,2], col = rgb(0,0,0,0.1))

# subset of high loading snps
chrom1.indexHighLoad <- lrLD.ind.chrom1[abs(lrLD.chrom1$v[,1]) > 0.01]

predict(chrom1.indexHighLoad)

new.pc.chrom1 <- predict(lrLD.chrom1)
df.pcs.chrom1 <- data.frame(new.pc.chrom1, pop = bigSNP.chrom1$fam$family.ID)
colnames(df.pcs.chrom1) <- c(paste0("PC", 1:10), "pop")
write.table(df.pcs.chrom1, "pcs_Chrom1.txt", row.names = FALSE, col.names = FALSE)
for(i in 1:length(pop.names)){
  df.pcs.chrom1$pop[df.pcs.chrom1$pop == paste0("Pop", i)] <- pop.names[i]
}
df.pcs.chrom1$pop <- factor(df.pcs.chrom1$pop, levels = pop.names)
                        
  
# plot chrom1
png(paste0(folderOut, "AllData_Chrom1_PC1_2.png"),  width = 7, height = 7, units = 'in', res = 300)
ggplot(df.pcs.chrom1, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) +
  scale_color_manual(name = "Population", 
                     values = viridis(length(pop.names))) + 
  labs(title = "All Samples Chrom 1") +
  theme_bw()
dev.off()

# separate by geographic location
ice.pops.chrom1 <- df.pcs.chrom1[df.pcs.chrom1$pop == "Ice.SWOff" | df.pcs.chrom1$pop ==  "Ice.SWNear" | df.pcs.chrom1$pop == "Ice.NEOff" | df.pcs.chrom1$pop == "Ice.NENear" ,]
us.pops.chrom1 <- df.pcs.chrom1[df.pcs.chrom1$pop == "Mass.Winter" | df.pcs.chrom1$pop ==  "Mass.Red" | df.pcs.chrom1$pop ==  "Mass.Spring" | df.pcs.chrom1$pop ==  "Cashes.Olive" | df.pcs.chrom1$pop ==  "Cashes.Red",]

png(paste0(folderOut, "IcelandicCod_chrom1_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(ice.pops.chrom1, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = viridis(length(pop.names))[6:9]) +
  labs(title = "Icelandic Cod Chrom1") + 
  theme_bw()
dev.off()

png(paste0(folderOut, "USCod_chrom1_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(us.pops.chrom1, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = viridis(length(pop.names))[1:5]) +
  labs(title = "Gulf of Maine Cod Chrom 1") + 
  theme_bw()
dev.off()


######
loc.subset <- sample(lrLD.ind.chrom1, 1500, replace = FALSE)

G_sub <- G_chrom1[,loc.subset]
POPNAMEID <-  df.pcs.chrom1$pop
heatmap(G_sub, cexCol = 0.3,  useRaster=TRUE,
        scale="none",
        Colv = NA, 
        #Rowv = NA,
        labRow = POPNAMEID, xlab="Mutation ID" ) 

#Colv = NA keeps the position stable instead of clustering with the position as well

lrLD_pos <- pos_chrom1[lrLD.ind.chrom1]
df.LDpos <- data.frame(lrLD.ind.chrom1 = lrLD.ind.chrom1, lrLD_pos = lrLD_pos)
df.invPos <- df.LDpos[df.LDpos$lrLD_pos > 11276210 & df.LDpos$lrLD_pos < 28335440,]
df.nonInvPos <- df.LDpos[df.LDpos$lrLD_pos <= 11276210 | df.LDpos$lrLD_pos >= 28335440,]

# add an 
plot(lrLD.ind.chrom1, lrLD.chrom1$v[,1]) #, col = rgb(0,0,0,0.1))
plot(lrLD.ind.chrom1, lrLD.chrom1$v[,2], col = rgb(0,0,0,0.1))

boxplot(df.pcs.chrom1$PC1~df.pcs.chrom1$pop)


oPC1 <- order(df.pcs.chrom1$PC1) #make sure this is the INDEXes of the individuals in the right order from lowest to highest PC loading
heatmap(G_chrom1[oPC1, df.invPos$lrLD.ind.chrom1], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
         middle="white"))
oPC2 <- order(df.pcs.chrom1$PC2)
heatmap(G_chrom1[oPC2, df.invPos$lrLD.ind.chrom1], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC2],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

heatmap(G_chrom1[oPC1, df.nonInvPos$lrLD.ind.chrom1], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))
heatmap(G_chrom1[oPC2, df.nonInvPos$lrLD.ind.chrom1], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC2],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                      middle="white"))

heatmap(G_chrom1[oPC1, chrom1.indexHighLoad], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

#subset for high absolute loadings



### end chrom1
########################################################################################################

########################################################################################################
### start chrom 2
snp_plinkQC(plink2, paste0(folderIn, "VarCall_NC_044049.1.f"),file.type = "--bfile", 
            maf = 0.05, geno = 0, "--allow-extra-chr", prefix.out = paste0(folderIn, "VarCall_NC_044049.1.f_QC"))
bed.file.filt.chrom2 <- snp_readBed(paste0(folderIn,"VarCall_NC_044049.1.f_QC.bed"))
bigSNP.chrom2 <- snp_attach(bed.file.filt.chrom2)
G_chrom2 <- bigSNP.chrom2$genotypes
pos_chrom2 <- bigSNP.chrom2$map$physical.pos

G_chrom2[1:100]
n <- nrow(G_chrom2)
m <- ncol(G_chrom2)
table(G_chrom2[1:m])
sum(is.na(G_chrom2[1:m]))
#884 10260927
chrom_2 <- bigSNP.chrom2$map$chromosome
for(j in 1:nrow(chromNames)){
  chrom_2[chrom_2 == chromNames$chromID[j]] <- j
}
chrom_2 <- as.numeric(chrom_2)

lrLD.chrom2 <- snp_autoSVD(G=G_chrom2, infos.chr = chrom_2,
                           infos.pos = pos_chrom2, size = 1)

str(lrLD.chrom2)
lrLD.ind.chrom2 <- attr(lrLD.chrom2, "subset") # this contains new set of pruned SNPs

length(lrLD.ind.chrom2) 
length(G_chrom2[lrLD.ind.chrom2])
lrLD.tb.chrom2 <- attr(lrLD.chrom2, "lrldr")

#loadings pc1
plot(lrLD.ind.chrom2, lrLD.chrom2$v[,1], col = rgb(0,0,0,0.1))
#loadings pc2
plot(lrLD.ind.chrom2, lrLD.chrom2$v[,2], col = rgb(0,0,0,0.1))

new.pc.chrom2 <- predict(lrLD.chrom2)
df.pcs.chrom2 <- data.frame(new.pc.chrom2, pop = bigSNP.chrom2$fam$family.ID)
colnames(df.pcs.chrom2) <- c(paste0("PC", 1:10), "pop")
write.table(df.pcs.chrom2, "pcs_Chrom2.txt", row.names = FALSE, col.names = FALSE)

for(i in 1:length(pop.names)){
  df.pcs.chrom2$pop[df.pcs.chrom2$pop  == paste0("Pop", i)] <- pop.names[i]
}
df.pcs.chrom2$pop <- factor(df.pcs.chrom2$pop, levels = pop.names)

# plot chrom1
png(paste0(folderOut, "AllData_Chrom2_PC1_2.png"),  width = 7, height = 7, units = 'in', res = 300)
ggplot(df.pcs.chrom2, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) +
  scale_color_manual(name = "Population", 
                     values = viridis(length(pop.names))) + 
  labs(title = "All Samples Chrom 2") +
  theme_bw()
dev.off()

# separate by geographic location
ice.pops.chrom2 <- df.pcs.chrom2[df.pcs.chrom2$pop == "Ice.SWOff" | df.pcs.chrom2$pop ==  "Ice.SWNear" | df.pcs.chrom2$pop == "Ice.NEOff" | df.pcs.chrom2$pop == "Ice.NENear" ,]
us.pops.chrom2 <- df.pcs.chrom2[df.pcs.chrom2$pop == "Mass.Winter" | df.pcs.chrom2$pop ==  "Mass.Red" | df.pcs.chrom2$pop ==  "Mass.Spring" | df.pcs.chrom2$pop ==  "Cashes.Olive" | df.pcs.chrom2$pop ==  "Cashes.Red",]

png(paste0(folderOut, "IcelandicCod_chrom2_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(ice.pops.chrom2, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = c("navyblue", "firebrick","cornflowerblue", "red")) +
  labs(title = "Icelandic Cod Chrom 2") + 
  theme_bw()
dev.off()

png(paste0(folderOut, "USCod_chrom2_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(us.pops.chrom2, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = c("navyblue", "firebrick", "goldenrod", "cornflowerblue", "red")) +
  labs(title = "Gulf of Maine Cod Chrom 2") + 
  theme_bw()
dev.off()

head(lrLD.ind.chrom2[abs(lrLD.chrom2$v[,1]) > 0.01])
lrLD_pos_c2 <- pos_chrom2[lrLD.ind.chrom2]
df.LDpos_c2 <- data.frame(lrLD.ind.chrom2 = lrLD.ind.chrom2, lrLD_pos = lrLD_pos_c2)
df.invPos_c2 <- df.LDpos_c2[df.LDpos_c2$lrLD_pos > 1673 & df.LDpos_c2$lrLD_pos < 23148,]
df.nonInvPos_c2 <- df.LDpos_c2[df.LDpos_c2$lrLD_pos <= 1673 | df.LDpos_c2$lrLD_pos >= 23148,]

oPC1_c2 <- order(df.pcs.chrom2$PC1) #make sure this is the INDEXes of the individuals in the right order from lowest to highest PC loading

heatmap(G_chrom2[oPC1_c2, df.invPos_c2$lrLD.ind.chrom2], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1_c2],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

oPC2_c2 <- order(df.pcs.chrom2$PC2)
heatmap(G_chrom2[oPC2, df.invPos_c2$lrLD.ind.chrom2], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC2_c2],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

heatmap(G_chrom2[oPC1_c2, df.nonInvPos_c2$lrLD.ind.chrom2], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1_c2],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

#### end chrom 2
########################################################################################################


########################################################################################################
#### start Chrom 7 
snp_plinkQC(plink2, paste0(folderIn, "VarCall_NC_044054.1.f"),file.type = "--bfile", 
            maf = 0.05, geno = 0, "--allow-extra-chr", prefix.out = paste0(folderIn, "VarCall_NC_044054.1.f_QC"))
bed.file.filt.chrom8 <- snp_readBed(paste0(folderIn,"VarCall_NC_044054.1.f_QC.bed"))
bigSNP.chrom7 <- snp_attach(bed.file.filt.chrom8)
G_chrom7 <- bigSNP.chrom7$genotypes
pos_chrom7 <- bigSNP.chrom7$map$physical.pos

G_chrom7[1:100]
n <- nrow(G_chrom7)
m <- ncol(G_chrom7)
table(G_chrom7[1:m])
sum(is.na(G_chrom7[1:m]))

chrom_7 <- bigSNP.chrom7$map$chromosome
for(j in 1:nrow(chromNames)){
  chrom_7[chrom_7 == chromNames$chromID[j]] <- j
}
chrom_7 <- as.numeric(chrom_7)

lrLD.chrom7 <- snp_autoSVD(G=G_chrom7, infos.chr = chrom_7,
                           infos.pos = pos_chrom7, size = 1)

str(lrLD.chrom7)
lrLD.ind.chrom7 <- attr(lrLD.chrom7, "subset") # this contains new set of pruned SNPs

length(lrLD.ind.chrom7) 
length(G_chrom7[lrLD.ind.chrom7])
lrLD.tb.chrom7 <- attr(lrLD.chrom7, "lrldr")

#loadings pc1
plot(lrLD.ind.chrom7, lrLD.chrom7$v[,1], col = rgb(0,0,0,0.1))
#loadings pc2
plot(lrLD.ind.chrom7, lrLD.chrom7$v[,2], col = rgb(0,0,0,0.1))

new.pc.chrom7 <- predict(lrLD.chrom7)
df.pcs.chrom7 <- data.frame(new.pc.chrom7, pop = bigSNP.chrom7$fam$family.ID)
colnames(df.pcs.chrom7) <- c(paste0("PC", 1:10), "pop")
write.table(df.pcs.chrom7, "pcs_Chrom7.txt", row.names = FALSE, col.names = FALSE)

for(i in 1:length(pop.names)){
  df.pcs.chrom7$pop[df.pcs.chrom7$pop  == paste0("Pop", i)] <- pop.names[i]
}
df.pcs.chrom7$pop <- factor(df.pcs.chrom7$pop, levels = pop.names)

# plot chrom1
png(paste0(folderOut, "AllData_Chrom7_PC1_2.png"),  width = 7, height = 7, units = 'in', res = 300)
ggplot(df.pcs.chrom7, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) +
  scale_color_manual(name = "Population", 
                     values = viridis(length(pop.names))) + 
  labs(title = "All Samples Chrom 7") +
  theme_bw()
dev.off()

# separate by geographic location
ice.pops.chrom7 <- df.pcs.chrom7[df.pcs.chrom7$pop == "Ice.SWOff" | df.pcs.chrom7$pop ==  "Ice.SWNear" | df.pcs.chrom7$pop == "Ice.NEOff" | df.pcs.chrom7$pop == "Ice.NENear" ,]
us.pops.chrom7 <- df.pcs.chrom7[df.pcs.chrom7$pop == "Mass.Winter" | df.pcs.chrom7$pop ==  "Mass.Red" | df.pcs.chrom7$pop ==  "Mass.Spring" | df.pcs.chrom7$pop ==  "Cashes.Olive" | df.pcs.chrom7$pop ==  "Cashes.Red",]

png(paste0(folderOut, "IcelandicCod_chrom7_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(ice.pops.chrom7, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = c("navyblue", "firebrick","cornflowerblue", "red")) +
  labs(title = "Icelandic Cod Chrom 7") + 
  theme_bw()
dev.off()

png(paste0(folderOut, "USCod_chrom7_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(us.pops.chrom7, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = c("navyblue", "firebrick", "goldenrod", "cornflowerblue", "red")) +
  labs(title = "Gulf of Maine Cod Chrom 7") + 
  theme_bw()
dev.off()

tail(lrLD.ind.chrom7[abs(lrLD.chrom7$v[,1]) > 0.01])
lrLD_pos_c7 <- pos_chrom7[lrLD.ind.chrom7]
df.LDpos_c7 <- data.frame(lrLD.ind.chrom7 = lrLD.ind.chrom7, lrLD_pos = lrLD_pos_c7)
df.invPos_c7 <- df.LDpos_c7[df.LDpos_c7$lrLD.ind.chrom7 > 52273 & df.LDpos_c7$lrLD.ind.chrom7 < 100353,]
df.nonInvPos_c7 <- df.LDpos_c7[df.LDpos_c7$lrLD.ind.chrom7 <= 52273 | df.LDpos_c7$lrLD.ind.chrom7 >= 100353,]

## inverted region
oPC1_c7 <- order(df.pcs.chrom7$PC1) #make sure this is the INDEXes of the individuals in the right order from lowest to highest PC loading
heatmap(G_chrom7[oPC1_c7, df.invPos_c7$lrLD.ind.chrom7], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1_c7],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

oPC2_c7 <- order(df.pcs.chrom7$PC2)
heatmap(G_chrom7[oPC2, df.invPos_c7$lrLD.ind.chrom7], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC2_c7],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

## Non inverted region
heatmap(G_chrom7[oPC1_c7, df.nonInvPos_c7$lrLD.ind.chrom7], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1_c7],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

heatmap(G_chrom7[oPC2, df.nonInvPos_c7$lrLD.ind.chrom7], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC2_c7],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))
# end chrom 7
########################################################################################################

########################################################################################################
#### start Chrom 12 
snp_plinkQC(plink2, paste0(folderIn, "VarCall_NC_044059.1.f"),file.type = "--bfile", 
            maf = 0.05, geno = 0, "--allow-extra-chr", prefix.out = paste0(folderIn, "VarCall_NC_044059.1.f_QC"))
bed.file.filt.chrom12 <- snp_readBed(paste0(folderIn,"VarCall_NC_044059.1.f_QC.bed"))
bigSNP.chrom12 <- snp_attach(bed.file.filt.chrom12)
G_chrom12 <- bigSNP.chrom12$genotypes
pos_chrom12 <- bigSNP.chrom12$map$physical.pos

G_chrom12[1:100]
n <- nrow(G_chrom12)
m <- ncol(G_chrom12)
table(G_chrom12[1:m])
sum(is.na(G_chrom12[1:m]))

chrom_12 <- bigSNP.chrom12$map$chromosome
for(j in 1:nrow(chromNames)){
  chrom_12[chrom_12 == chromNames$chromID[j]] <- j
}
chrom_12 <- as.numeric(chrom_12)

lrLD.chrom12 <- snp_autoSVD(G=G_chrom12, infos.chr = chrom_12,
                           infos.pos = pos_chrom12, size = 1)

str(lrLD.chrom12)
lrLD.ind.chrom12 <- attr(lrLD.chrom12, "subset") # this contains new set of pruned SNPs

length(lrLD.ind.chrom12) 
length(G_chrom12[lrLD.ind.chrom12])
lrLD.tb.chrom12 <- attr(lrLD.chrom12, "lrldr")

#loadings pc1
plot(lrLD.ind.chrom12, lrLD.chrom12$v[,1], col = rgb(0,0,0,0.1))
#loadings pc2
plot(lrLD.ind.chrom12, lrLD.chrom12$v[,2], col = rgb(0,0,0,0.1))

new.pc.chrom12 <- predict(lrLD.chrom12)
df.pcs.chrom12 <- data.frame(new.pc.chrom12, pop = bigSNP.chrom12$fam$family.ID)
colnames(df.pcs.chrom12) <- c(paste0("PC", 1:10), "pop")
write.table(df.pcs.chrom12, "pcs_Chrom12.txt", row.names = FALSE, col.names = FALSE)

for(i in 1:length(pop.names)){
  df.pcs.chrom12$pop[df.pcs.chrom12$pop  == paste0("Pop", i)] <- pop.names[i]
}
df.pcs.chrom12$pop <- factor(df.pcs.chrom12$pop, levels = pop.names)

# plot chrom1
png(paste0(folderOut, "AllData_Chrom12_PC1_2.png"),  width = 7, height = 7, units = 'in', res = 300)
ggplot(df.pcs.chrom12, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) +
  scale_color_manual(name = "Population", 
                     values = viridis(length(pop.names))) + 
  labs(title = "All Samples Chrom 12") +
  theme_bw()
dev.off()

# separate by geographic location
ice.pops.chrom12 <- df.pcs.chrom12[df.pcs.chrom12$pop == "Ice.SWOff" | df.pcs.chrom12$pop ==  "Ice.SWNear" | df.pcs.chrom12$pop == "Ice.NEOff" | df.pcs.chrom12$pop == "Ice.NENear" ,]
us.pops.chrom12 <- df.pcs.chrom12[df.pcs.chrom12$pop == "Mass.Winter" | df.pcs.chrom12$pop ==  "Mass.Red" | df.pcs.chrom12$pop ==  "Mass.Spring" | df.pcs.chrom12$pop ==  "Cashes.Olive" | df.pcs.chrom12$pop ==  "Cashes.Red",]

png(paste0(folderOut, "IcelandicCod_chrom12_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(ice.pops.chrom12, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = c("navyblue", "firebrick","cornflowerblue", "red")) +
  labs(title = "Icelandic Cod Chrom 12") + 
  theme_bw()
dev.off()

png(paste0(folderOut, "USCod_chrom12_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(us.pops.chrom12, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = c("navyblue", "firebrick", "goldenrod", "cornflowerblue", "red")) +
  labs(title = "Gulf of Maine Cod Chrom 12") + 
  theme_bw()
dev.off()

head(lrLD.ind.chrom12[abs(lrLD.chrom12$v[,1]) > 0.013], n=100)
lrLD_pos_c12 <- pos_chrom12[lrLD.ind.chrom12]
df.LDpos_c12 <- data.frame(lrLD.ind.chrom12 = lrLD.ind.chrom12, lrLD_pos = lrLD_pos_c12)
df.invPos_c12 <- df.LDpos_c12[df.LDpos_c12$lrLD.ind.chrom12 > 3993 & df.LDpos_c12$lrLD.ind.chrom12 <  46959,]
df.nonInvPos_c12 <- df.LDpos_c12[df.LDpos_c12$lrLD.ind.chrom12 <= 3993 | df.LDpos_c12$lrLD.ind.chrom12 >=  46959,]

## inverted region
oPC1_c12 <- order(df.pcs.chrom12$PC1) #make sure this is the INDEXes of the individuals in the right order from lowest to highest PC loading
heatmap(G_chrom12[oPC1_c12, df.invPos_c12$lrLD.ind.chrom12], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1_c12],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

oPC2_c12 <- order(df.pcs.chrom12$PC2)
heatmap(G_chrom12[oPC2, df.invPos_c12$lrLD.ind.chrom12], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC2_c12],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

## Non inverted region
heatmap(G_chrom12[oPC1_c12, df.nonInvPos_c12$lrLD.ind.chrom12], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC1_c12],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))

heatmap(G_chrom12[oPC2, df.nonInvPos_c12$lrLD.ind.chrom12], 
        Colv = NA, Rowv = NA,
        labRow = POPNAMEID[oPC2_c12],
        cexCol = 0.3, cexRow = 0.25, useRaster=TRUE,
        col=two.colors(3, start = "cornflowerblue", end="red", 
                       middle="white"))


# end chrom 12
########################################################################################################



########################################################################################################
# OLD CODE

#bed.file <- snp_readBed(paste0(folderIn, "VarCall_NC_044048.1.f.bed"))
#bigSNPobj <- snp_attach(bed.file)

pos <- bigSNPobj$map$physical.pos
n <- nrow(G_copy)
m <- ncol(G_copy)
table(G_copy[1:m])
sum(is.na(G_copy[1:m]))

## impute missing data so we can use in other functions
indNA <- which(is.na(G_copy[1:m])) # get indices of missing data don't actually need this
system.time(G_NoMiss <- snp_fastImpute(G, chrom)) # compute fastImpute
G_copy$code256 <- bigsnpr:::CODE_IMPUTE_PRED # add values that are estimated from imputation


# remove SNPs that are in LD and only keep the one in the region with highest MAF
# clump.ind <- snp_clumping(G_copy, infos.chr = chrom)
# length(clump.ind)
# length(G_copy[clump.ind])
# length(chrom[clump.ind])
# length(pos[clump.ind])
  
new.pc.full <- predict(lrLD)
df.pcs.full <- data.frame(new.pc.full, pop = bigSNP$fam$family.ID)
colnames(df.pcs.full) <- c(paste0("PC", 1:10), "pop")
write.table(df.pcs.full, "pcs_ALLDATA.txt", row.names = FALSE, col.names = FALSE)

# plot all data
png(paste0(folderOut, "AllData_FullGenome_PC1_2.png"),  width = 7, height = 7, units = 'in', res = 300)
ggplot(df.pcs.full, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) +
  scale_color_manual(name = "Population", 
                     values = viridis(length(pop.names))) + 
  labs(title = "All Samples Full Genome") +
  theme_bw()
dev.off()

# separate by geographic location
ice.pops.full <- df.pcs.full[df.pcs.full$pop == "Ice.SWOff" | df.pcs.full$pop ==  "Ice.SWNear" | df.pcs.full$pop == "Ice.NEOff" | df.pcs.full$pop == "Ice.NENear" ,]
us.pops.full <- df.pcs.full[df.pcs.full$pop == "Mass.Winter" | df.pcs.full$pop ==  "Mass.Red" | df.pcs.full$pop ==  "Mass.Spring" | df.pcs.full$pop ==  "Cashes.Olive" | df.pcs.full$pop ==  "Cashes.Red",]

png(paste0(folderOut, "IcelandicCod_FullGenome_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(ice.pops.full, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = viridis(length(pop.names))[6:9]) +
  labs(title = "Icelandic Cod Full Genome") + 
  theme_bw()
dev.off()

png(paste0(folderOut, "USCod_FullGenome_PC1_2.png"), width = 7, height = 7, units = 'in', res = 300)
ggplot(us.pops.full, aes(x = PC1, y = PC2, color = pop)) + 
  geom_point(alpha = 0.75) + 
  scale_color_manual(name = "Population", values = viridis(length(pop.names))[1:5]) +
  labs(title = "Gulf of Maine Cod Full Genome") + 
  theme_bw()
dev.off()