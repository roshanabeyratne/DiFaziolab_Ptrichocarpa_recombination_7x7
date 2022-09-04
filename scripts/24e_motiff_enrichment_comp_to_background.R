### METADATA: 
# Written by CRA
# March 17, 2022
# Extract positive and negative sequences for motiff enrichment analysis
### 


# setting up working environment
setwd("/Volumes/GoogleDrive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7")
rm(list=ls())
library(tidyverse)
source(file = "./Stettler_14_linkage_mapping/new_scripts/CreateBinSpanCoordinates_function.R")
input_path<- "./Stettler_14_linkage_mapping/statistics/"
#

# subset CO locations to size <= 10k
CO_detailed <- read.csv(file=paste0(input_path, 'co.detailed.df.csv'), header=T, stringsAsFactors=F)
filtered_COs <- CO_detailed[CO_detailed$CO_SIZE <= 10000 & CO_detailed$CO_SIZE > 0, ]
hist(filtered_COs$CO_SIZE)
min(filtered_COs$CO_SIZE)
write.csv(filtered_COs, file=paste0(input_path, 'filtered_CO_le_10k_df.csv'), quote=F, row.names=F)

# create a dataframe with genomic windows of size 10k
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
chrom.size.df <- read.csv(file = paste0(input_path, "chrom.size.df_Stetv2.csv"), header=T, stringsAsFactors=F)
win.size <- 10000
total.bin.span.coord.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr01"
  bin.span.coord.df.temp <- CreateBinSpanCoordinates(chrom.name = chrom.name, 
                                                     bin.size = win.size, 
                                                     chrom.length = chrom.size.df[chrom.size.df$CHROM == chrom.name, "LENGTH"], 
                                                     shift.size = win.size)
  bin.span.coord.df.temp <- bin.span.coord.df.temp[-1,]
  bin.span.coord.df.temp <- bin.span.coord.df.temp[-(nrow(bin.span.coord.df.temp)),]
  bin.span.coord.df.temp2 <- cbind(rep(chrom.name, nrow(bin.span.coord.df.temp)), bin.span.coord.df.temp)
  total.bin.span.coord.df <- rbind(total.bin.span.coord.df, bin.span.coord.df.temp2)
}
total.bin.span.coord.df <- as.data.frame(total.bin.span.coord.df, stringsAsFactors = FALSE)
total.bin.span.coord.df[,2:3] <- apply(total.bin.span.coord.df[,2:3], 2, as.numeric)
colnames(total.bin.span.coord.df) <- c("CHROM", "START", "END")

# Find matching sized genomic locations

# randomly select genomic windows
neg_con_COs <- sample_n(total.bin.span.coord.df, size=nrow(filtered_COs))

# down-size randomly selected 10k window to the corresponding CO_size
neg_con_COs2 <- NULL
for (i in 1:nrow(filtered_COs)) {
  co_siz <- filtered_COs$CO_SIZE[i]
  chrom <- neg_con_COs$CHROM[i]
  new_start <- neg_con_COs$START[i]
  new_end <- neg_con_COs$START[i] + (co_siz-1)
  neg_con_COs2 <- rbind(neg_con_COs2, c(chrom, new_start, new_end))
}
neg_con_COs2 <- as.data.frame(neg_con_COs2, stringsAsFactors=F)
colnames(neg_con_COs2) <- colnames(neg_con_COs)
# hist(as.numeric(neg_con_COs2$START)-as.numeric(neg_con_COs2$END))
# head(neg_con_COs2)
write.csv(neg_con_COs2, file=paste0(input_path, 'neg_con_CO_le_10k_df10.csv'), quote=F, row.names=F)



