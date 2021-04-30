### Plot using samtools stats

Pop6_18017 <- read.table("samtools_alnCheck_Out/sam_stats/Pop6_18017aln.sorted.bam.stats.cov.txt")
Pop5_17278 <- read.table("samtools_alnCheck_Out/sam_stats/Pop5_17278aln.sorted.bam.stats.cov.txt")
Pop1_17291 <- read.table("samtools_alnCheck_Out/sam_stats/Pop1_17291aln.sorted.bam.stats.cov.txt")
Pop4_17236 <- read.table("samtools_alnCheck_Out/sam_stats/Pop4_17236aln.sorted.bam.stats.cov.txt")
samp.names <- c("Pop6_18017", "Pop5_17278", "Pop1_17291", "Pop4_17236")
par(mfrow = c(4,1))
Pop6_18017 <- Pop6_18017[1:50,3:4]
  colnames(Pop6_18017 ) <- c("coverage", "Frequency")
  barplot(Frequency~coverage, data = Pop6_18017 , main = "Pop6_18017",
          col = "cornflowerblue")
  
Pop5_17278 <- Pop5_17278[1:50,3:4]
  colnames(Pop5_17278) <- c("coverage", "Frequency")
  barplot(Frequency~coverage, data = Pop5_17278, main = "Pop5_17278",
          col = "cadetblue3")

Pop1_17291 <- Pop1_17291[1:50,3:4]
  colnames(Pop1_17291) <- c("coverage", "Frequency")
  barplot(Frequency~coverage, data = Pop1_17291, main = "Pop1_17291",
         col = "cadetblue1" )
  
Pop4_17236 <- Pop4_17236[1:50,3:4]
  colnames(Pop4_17236) <- c("coverage", "Frequency")
  barplot(Frequency~coverage, data = Pop4_17236, main = "Pop4_17236",
          col = "lightcyan")
  

