#!/usr/bin/R
setwd("/Users/cabeyrat/Google Drive File Stream/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: January 06, 2021
# USAGE: Input file include a detailed CO count and informative Het marker count per halfsib family per 30 kb window.
# PURPOSE: Analyze whether there is a sprurious association between number of informative markers flanking the window including the index window and the observed number of COs
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
statistics.output.path <- "./statistics/hotspot_analysis/30kb/"
########################################################################################################################################################################################################

idx.to.rmv <- c(12002, 12044) # These are the Chr17 indexes identified as Stet-14 assembly errors
all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_30kb.csv"), header = TRUE, stringsAsFactors = FALSE )
colnames(all.par.df.2) <- colnames(all.par.df.2) %>% str_replace("^X", "")
all.par.df.2 <- all.par.df.2[-idx.to.rmv, ]
info.mrkr.count.df <- read.csv(file=paste0(statistics.input.path, "info_mrkr_count_bin_size_30kb.csv"), stringsAsFactors=FALSE )
info.mrkr.count.df <- info.mrkr.count.df[-idx.to.rmv, ]

# Take the mean CO counts and mean informative marker count
CO.count <- apply(all.par.df.2[-(1:3)], 1, mean)
MRKR.count <-  apply(info.mrkr.count.df[-(1:3)], 1, mean)

if (length(CO.count) == length(MRKR.count)) {
  mrkr.vs.co.df <- cbind(info.mrkr.count.df[,1:3], AVG_MRKR_COUNT=MRKR.count,AVG_CO_COUNT=CO.count)
} else {
  stop("error in CO.count or info.mrkr.count.df")
}



# Obtain the CO count for idx window and informative marker count for idx window and flanking windows. 
flnk.size <- 2
detailed.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  subset.mrkr.vs.co.df <- mrkr.vs.co.df[mrkr.vs.co.df$CHROM==chrom.name,]
  feasible.idx <- c((1+flnk.size):(nrow(subset.mrkr.vs.co.df)-flnk.size))
  for (idx in feasible.idx) { # idx <- 1000
    co.count <- subset.mrkr.vs.co.df$AVG_CO_COUNT[idx]
    mrkr.idxs <- c((idx-flnk.size):(idx-1), idx, (idx+1):(idx+flnk.size))
    mrkr.count <- mean(subset.mrkr.vs.co.df$AVG_MRKR_COUNT[mrkr.idxs])
    return.vec <- c(subset.mrkr.vs.co.df$CHROM[idx], subset.mrkr.vs.co.df$START[idx], subset.mrkr.vs.co.df$END[idx], mrkr.count, co.count)
    detailed.df <- rbind(detailed.df, return.vec)
  }
  cat(chrom.name, "\n")
}
detailed.df <- as.data.frame(detailed.df, stringsAsFactors=FALSE)
colnames(detailed.df) <- c(colnames(mrkr.vs.co.df)[1:3], "AVG_MRKR_COUNT", "AVG_CO_COUNT")

plot(AVG_CO_COUNT ~ AVG_MRKR_COUNT, data=detailed.df)
idx.outliers <- which(detailed.df$AVG_MRKR_COUNT<=14 & detailed.df$AVG_CO_COUNT>0)
outlier.df <- detailed.df[idx.outliers, ]

hotspot.cands <- read.csv(file = paste0(statistics.output.path, "30kb_hotspots_bonferronni_cutoff.csv"), header=TRUE, stringsAsFactors=FALSE)
tmp.idx <- which(hotspot.cands$CHROM %in% outlier.df$CHROM & hotspot.cands$START %in% outlier.df$START)
hotspot.cands[tmp.idx, ]


# library(ggplot2)
# ggplot(detailed.df, aes(AVG_MRKR_COUNT, AVG_CO_COUNT)) + geom_bin2d(bins=100) + scale_fill_viridis_c(breaks = c(100, 500, 1500, 2500, 4000)) + theme(axis.text.x=element_text(angle=50) ) + 
#   scale_y_discrete(breaks=waiver()) + scale_x_discrete(breaks=waiver())
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################