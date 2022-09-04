#!/usr/bin/R
setwd("/Users/cabeyrat/Google Drive File Stream/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: January 15, 2021
# USAGE: interactive
# PURPOSE: Track the number of informative markers (Het in focal parent) in each 30kbps size bins 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
output.path.phased.markers <- "./dwnld_from_spruce/phased_markers/"
statistics.output.path <- "./statistics/hotspot_analysis/30kb/"
########################################################################################################################################################################################################
source(file = "./new_scripts/CreateBinSpanCoordinates_function.R")
########################################################################################################################################################################################################
# read-in the file with Chromosome size calculated for Nisqually-V_3
# chrom.size.df <- read.csv(file = paste0(statistics.input.path, "chrom.size.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# # modifying coordinates of Chr19 for Stettler
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr11"] <- 19000000
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr16"] <- 15494368
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr19"] <- 16500000

# read-in the file with Chromosome size calculated for Stet-14_V2
chrom.size.df <- read.csv(file = paste0(statistics.input.path, "chrom.size.df_Stetv2.csv"), header = TRUE, stringsAsFactors = FALSE)



# Assign window start and end points based on window size
win.size <- 30000 # This is the final window size statistical analysis done on.
total.bin.span.coord.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr01"
  bin.span.coord.df.temp <- CreateBinSpanCoordinates(chrom.name = chrom.name, bin.size = win.size, chrom.length = chrom.size.df[chrom.size.df$CHROM == chrom.name, "LENGTH"], shift.size = win.size)
  bin.span.coord.df.temp2 <- cbind(rep(chrom.name, nrow(bin.span.coord.df.temp)), bin.span.coord.df.temp)
  total.bin.span.coord.df <- rbind(total.bin.span.coord.df, bin.span.coord.df.temp2)
}
total.bin.span.coord.df <- as.data.frame(total.bin.span.coord.df, stringsAsFactors = FALSE)
total.bin.span.coord.df[,2:3] <- apply(total.bin.span.coord.df[,2:3], 2, as.numeric)
colnames(total.bin.span.coord.df) <- c("CHROM", "START", "END")


count.all <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  count.per.chrom <- NULL 
  for (par.name in par.name.vec) { # par.name <- par.name.vec[1]
     
     # The set of Het markers that were phased using oneMap R-package.
     OM.phased.df <- read.csv(file = paste0(output.path.phased.markers, chrom.name, "_phased_markers_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
     # subset the windows for the chromosome 
     tmp.coord.tab <- total.bin.span.coord.df[total.bin.span.coord.df$CHROM == chrom.name, ]
     tmp.count <- apply(tmp.coord.tab, 1, function(x, y=OM.phased.df$POS){ start <- as.numeric(x[2]); end <- as.numeric(x[3]); count <- sum(y>=start & y<=end); return(count)})
     # sum(tmp.count)
     count.per.chrom <- cbind(count.per.chrom, tmp.count)
  }
  count.all <- rbind(count.all, count.per.chrom)
  cat(chrom.name, " ")
}

if (nrow(count.all) == nrow(total.bin.span.coord.df)) {
  tmp.df <- cbind(total.bin.span.coord.df, count.all)
  colnames(tmp.df) <- c(colnames(total.bin.span.coord.df), par.name.vec)
  write.csv(tmp.df, file=paste0(statistics.input.path, "info_mrkr_count_bin_size_30kb.csv"), quote=FALSE, row.names=FALSE)
} else {
  stop("Error in the dimenstions of count.all !")
}


# plotting the median marker number 
total.count.per.win <- apply(tmp.df[,-(1:3)], 1, median)
jpeg(file = paste0(statistics.output.path, "Median_markrs_per_bin.jpg"), width = 1000, height = 500, quality = 5000)
plot(total.count.per.win, col='grey', type='l', pch=16, cex=0.25, ylab='Median marker count per 30kb bin')
abline(h=mean(total.count.per.win))
legend(1, 160, legend = c("Genomewide expected marker count per bin"), col = c("black"), lty = 1, lwd = 2, cex = 0.8)
dev.off()
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################