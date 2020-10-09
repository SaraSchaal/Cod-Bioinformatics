IDs <- read.csv("src/SampleNames.csv") 
nrow(IDs)
ncol(IDs)
nrow(unique(IDs))
duplicated(IDs)
IDs <- IDs[1:300,2]

write.table(IDs[IDs != "Pop8_18102" & 
                IDs != "Pop6_18016" &
                IDs != "Pop1_16237" &
                IDs != "Pop1_16230"], 
                "src/SampleNames.txt", row.names = FALSE, col.names = FALSE)


