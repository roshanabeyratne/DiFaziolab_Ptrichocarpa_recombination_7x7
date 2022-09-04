#!/usr/bin/R
setwd("/Users/cabeyrat/Google Drive File Stream/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: October 24, 2020; Modified: December 16, 2020
# USAGE: Input file cross-over events observed for each focal parent for each window for all the chromosomes. The input is produced with script 14b_... 

# Following is the input file format used 
# CHROM  START    END 1863 1909 1950 2048 2066 2283 4593 2365 2393 2515 2572 2683 6909 7073
# 1 Chr01      1   2699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# 2 Chr01   2700  32699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# 3 Chr01  32700  62699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# 4 Chr01  62700  92699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# 5 Chr01  92700 122699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# 6 Chr01 122700 152699    0    0    0    0    0    0    0    0    0    0    0    0    0    0

# PURPOSE: Check whether a few individuals are driving the CO count difference between sexes.
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)

########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
permutation.input.path <- "./statistics/heterochiasmy_analysis/wmtsa/permutation_tabs/"
# statistics.output.path <- "./statistics/hotspot_analysis/250kb/"
statistics.output.path <- "./statistics/heterochiasmy_analysis/"
sliding.avg.output.path <- "./statistics/heterochiasmy_analysis/sliding_average/"
########################################################################################################################################################################################################
# Read-in the data file which has CO-counts per half-sib family by window. The window size size was selected as 30kb as the median Co-region approximates 30kb. 
all.par.df.2 <- read.csv(file=paste0(statistics.input.path, "co.bin_size_30kb.csv"), header=TRUE, stringsAsFactors=FALSE )
colnames(all.par.df.2)[4:17] <- par.name.vec
# checking whether the total number of cross-over counts are the same as observed
sum(apply(all.par.df.2[,4:17], 1, sum)) # 38846. This checks out



# Binning the  permuted and actual CO-counts into window sizes of each scale slid by 30kb each.
for (chrom.name in chrom.name.vec[1]) {
  # Test analysis is carried-out for Chr01 for now until a analysis pipeline is established.
  # chrom.name <- "Chr19"
  test.df <- all.par.df.2[all.par.df.2$CHROM==chrom.name, -(1:3)]
  cos.vec.fem <- apply(test.df[, 1:7], 1 , sum)
  cos.vec.male <- apply(test.df[, 8:14], 1 , sum)
  cos.vec.sex.diff <- (cos.vec.fem-cos.vec.male)
  
  permut.fem.df <- read.csv(paste0(permutation.input.path, chrom.name, "_fem.permut_df.csv"), stringsAsFactors=FALSE, header=TRUE )
  permut.male.df <- read.csv(paste0(permutation.input.path, chrom.name, "_male.permut_df.csv"), stringsAsFactors=FALSE, header=TRUE )
  permut.sex.diff.df <- (permut.fem.df-permut.male.df) 
  
  
  
  # Following is a function that divides a given chromosome into non overlapping windows. 
  CreateBinSpanCoordinates <- function(chrom.name, bin.size, chrom.length, shift.size) {
    bin.spacing <- seq(from = bin.size, to = chrom.length, by = shift.size)
    num.pos <- length(seq(from = bin.size, to = chrom.length, by = shift.size))
    
    snipping.size <- chrom.length - bin.spacing[num.pos]
    if (snipping.size %% 2 == 0) {
      snipping.size.one.end = snipping.size / 2
    } else {
      snipping.size.one.end = (snipping.size + 1) / 2
    }
    
    end.sites <- bin.spacing + snipping.size.one.end
    start.sites <- (end.sites - bin.size) + 1
    start.end.sites <- cbind(start.sites, end.sites)
    if (start.end.sites[1,1] > 1) {
      start.end.sites <- rbind(c(1, (start.end.sites[1,1] - 1)), start.end.sites)
    } 
    if (start.end.sites[nrow(start.end.sites),2] < chrom.length) {
      start.end.sites <- rbind(start.end.sites, c((start.end.sites[nrow(start.end.sites),2] + 1), chrom.length))
    }
    return(start.end.sites)
  }
  
  
  
  for (scale in 1:9) {
    # scale<-5
    # based on the wavelet analysis the appropriate window size for heterochiasmy analysis is assigned as (2^6)*30kb = 1920kb; Assigning window start and end points based on this new window size.
    selected.scale <- (2^scale)
    bin.span.coord.df <- as.data.frame(CreateBinSpanCoordinates(chrom.name=chrom.name, bin.size=selected.scale, chrom.length=(nrow(test.df)) , shift.size=1 ), stringsAsFactors=FALSE )
    colnames(bin.span.coord.df) <- c("START", "END")
    
    # dataframes with number of cross-overs and number of individuals used for each parent for subsequent hotspot analysis 
    num.ind.df <- read.csv(paste0(statistics.input.path, "number_of_individuals_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
    num.COs.df <- read.csv(paste0(statistics.input.path, "number_of_crossovers_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
    colnames(num.ind.df) <- colnames(num.ind.df) %>% str_replace("^X", "")
    colnames(num.COs.df) <- colnames(num.COs.df) %>% str_replace("^X", "")
    
    # binning the CO-counts for each half-sib family.
    binned.cos.by.par.df <- NULL
    for (par.name in par.name.vec) {
      tmp.permut <- apply(bin.span.coord.df, 1, y=as.numeric(test.df[,par.name]), function(x, y){start.pos <- x[1]; end.pos <- x[2]; tmp.count <- sum(y[start.pos:end.pos]); return(tmp.count)})
      binned.cos.by.par.df <- cbind(binned.cos.by.par.df, tmp.permut)
    }
    binned.cos.by.par.df <- as.data.frame(binned.cos.by.par.df, stringsAsFactors=FALSE)
    colnames(binned.cos.by.par.df) <- par.name.vec
    write.csv(binned.cos.by.par.df, file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs_by_par.csv"), quote=FALSE, row.names=FALSE)
    binned.cos.by.par.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs_by_par.csv"), stringsAsFactors=FALSE, header=TRUE)
    
    
    # for actual female and male groups
    binned.cos.fem <- apply(bin.span.coord.df, 1, y=cos.vec.fem, function(x, y){start.pos <- x[1]; end.pos <- x[2]; tmp.count <- sum(y[start.pos:end.pos]); return(tmp.count)})
    binned.cos.male <- apply(bin.span.coord.df, 1, y=cos.vec.male, function(x, y){start.pos <- x[1]; end.pos <- x[2]; tmp.count <- sum(y[start.pos:end.pos]); return(tmp.count)})
    actual.binned.COs <- cbind(binned.cos.fem, binned.cos.male)
    write.csv(actual.binned.COs, file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs.csv"), quote=FALSE, row.names=FALSE)
    actual.binned.COs <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs.csv"), stringsAsFactors=FALSE, header=TRUE)
    
    
    # for permuted female and male groups
    binned.permut.cos.fem.df <- NULL
    for (iter in 1:ncol(permut.fem.df)) {
      tmp.permut.fem <- apply(bin.span.coord.df, 1, y=as.numeric(permut.fem.df[,iter]), function(x, y){start.pos <- x[1]; end.pos <- x[2]; tmp.count <- sum(y[start.pos:end.pos]); return(tmp.count)})
      binned.permut.cos.fem.df <- cbind(binned.permut.cos.fem.df, tmp.permut.fem)
      cat(iter, " ")
    }
    binned.permut.cos.male.df <- NULL
    for (iter in 1:ncol(permut.male.df)) {
      tmp.permut.male <- apply(bin.span.coord.df, 1, y=as.numeric(permut.male.df[,iter]), function(x, y){start.pos <- x[1]; end.pos <- x[2]; tmp.count <- sum(y[start.pos:end.pos]); return(tmp.count)})
      binned.permut.cos.male.df <- cbind(binned.permut.cos.male.df, tmp.permut.male)
      cat(iter, " ")
    }
    write.csv(binned.permut.cos.fem.df, file=paste0(sliding.avg.output.path, chrom.name, "_", scale,"_permuted_binned_COs_FEM.csv"), quote=FALSE, row.names=FALSE)
    write.csv(binned.permut.cos.male.df, file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_MALE.csv"), quote=FALSE, row.names=FALSE)
    binned.permut.cos.fem.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
    binned.permut.cos.male.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
  }
  cat(chrom.name, "\n")
}


for (chrom.name in chrom.name.vec[1]) {
  jpeg(filename = paste0(sliding.avg.output.path, chrom.name, "_sliding_window_average.jpg"), width = 2000, height = 1000, quality = 150)
  op <- par(mfrow=c(4,3), mar=c(2,1,1,1))
  for (scale in 1:8) {
    binned.cos.by.par.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs_by_par.csv"), stringsAsFactors=FALSE, header=TRUE)
    colnames(binned.cos.by.par.df) <- par.name.vec
    actual.binned.COs <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs.csv"), stringsAsFactors=FALSE, header=TRUE)
    binned.permut.cos.fem.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
    binned.permut.cos.male.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
    
    
    # viasualize the total CO-count trends for permuted data-set and the two sexes
    # plot the permuted means for each sex 
    for (iter in 1:ncol(binned.permut.cos.male.df)) {
      if (iter == 1) {
        max.y.lim <- max(binned.permut.cos.male.df[, iter]/829)+0.1
        plot(binned.permut.cos.male.df[, iter]/829, col='grey', ylim=c(0, max.y.lim), main=paste0(chrom.name, " d", scale), xlab="overlapping window idx", ylab="mean CO-count")
      } else {
        lines(binned.permut.cos.male.df[, iter]/829, col='grey')
      }
    }
    for (iter in 1:ncol(binned.permut.cos.fem.df)) {lines(binned.permut.cos.fem.df[, iter]/829, col='grey')}
    
    # plot individual half-sib family counts per parent as well as the mean for each SEX
    for (par.name in par.name.vec[1:14]) { # par.name <- par.name.vec[1]
      colr <- ifelse(par.name %in% par.name.vec[1:7], yes=rgb(255, 0, 0, max=255, alpha=80, names="red50"), no=rgb(0, 0, 255, max=255, alpha=80, names="blue50"))
      if (par.name == par.name.vec[1]) {
        lines(binned.cos.by.par.df[, par.name] / num.ind.df[which(chrom.name == chrom.name.vec), par.name], col=colr)
      } else {
        lines(binned.cos.by.par.df[, par.name] / num.ind.df[which(chrom.name == chrom.name.vec), par.name], col=colr)
      }
    }
    lines(actual.binned.COs[,1]/829, type='l', col='red', lwd=1.5)
    lines(actual.binned.COs[,2]/829, col='blue', lwd=1.5)
  }
  dev.off()
}
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################
# WIP









# Differences between sexes
sex.diff.actual <- (binned.cos.fem - binned.cos.male)
sex.diff.permuted.df <-( binned.permut.cos.fem.df - binned.permut.cos.male.df)

p.val.vec <- NULL
for (win.idx in 1:nrow(sex.diff.permuted.df)) {
  tmp.p <- sum(abs(sex.diff.permuted.df[win.idx, ]) > abs(sex.diff.actual[win.idx])) / length(sex.diff.permuted.df[win.idx, ])
  p.val.vec <- c(p.val.vec, tmp.p)
}

# visualize sex difference
for (iter in 1:ncol(sex.diff.permuted.df)) {
  if (iter == 1) {
    plot(sex.diff.permuted.df[, iter], ylim=c(-100, 100), col='grey', main=chrom.name, xlab="", ylab="CO-count difference")
  } else {
    lines(sex.diff.permuted.df[, iter], col='grey')
  }
}
lines(sex.diff.actual, type='l', col='black', lwd=3)
lines(-log10(p.val.vec)*25, col=rgb(0, 255, 0, max=255, alpha=80, names="green50"), lwd=1.5)

# for (p in 1:length(p.val.vec)) {
# }
par(mfrow=c(1,1))












# just checking whether the expected pattern is observed for cross-over counts over each chromosome for each parent
library(RColorBrewer)
# colors <- c(brewer.pal(7, "Set3"), brewer.pal(7, "Set2"))
colors <- c('red', 'blue')
par(mfrow = c(1,2))
for (chrom.name in chrom.name.vec[1:2]) { # chrom.name <- chrom.name.vec[1]
  # jpeg(file = paste0(statistics.output.path, chrom.name, "_observed_p_value_on_Poisson_NULL.jpg"), width = 16, height = 9, unit = "in", res = 600)
  # pdf(file = paste0(statistics.output.path, chrom.name, "_observed_p_value_on_Poisson_NULL.pdf"), width = 16, height = 9, onefile = TRUE)
  # par(mfrow = c(2, 7))
  # colors <- c(brewer.pal(7, "Set3"), brewer.pal(7, "Set2"))
  
  
  # par.name.vec[c(1:12, 14)]
  for (par.name in par.name.vec[1:14]) { # par.name <- par.name.vec[1]
    if (par.name %in% par.name.vec[1:7]) {
      colr = 'red'
    } else {
      colr = 'blue'
    }
    if (par.name == par.name.vec[1]) {
      plot((all.par.df.2[all.par.df.2$CHROM == chrom.name, par.name] / num.ind.df[which(chrom.name == chrom.name.vec), par.name]), typ = "l", col = colr, main = chrom.name, xlab = "sliding_window_(not the win idx)",
           ylab = "number of cross-overs", ylim=c(0,0.5))
    } else {
      lines((all.par.df.2[all.par.df.2$CHROM == chrom.name, par.name] / num.ind.df[which(chrom.name == chrom.name.vec), par.name]), col = colr,
            # main = par.name, xlab = "window.idx",
            ylab = "number of cross-overs")
    }
  }
  
  # The mean of individual sexes
  female.vec <- apply((all.par.df.2[all.par.df.2$CHROM == chrom.name, par.name.vec[1:7]]), 1, sum) / 829
  male.ve <- apply((all.par.df.2[all.par.df.2$CHROM == chrom.name, par.name.vec[8:14]]), 1, sum) / 829 
  
  lines(female.vec, typ = "l", col = 'red', 
        # main = par.name, xlab = "window.idx", ylab = "number of cross-overs", ylim=c(0,0.5), 
        lwd=3)
  lines(male.ve, typ = "l", col = 'blue', lwd=3)
  
  
 
}
dev.off()
par(mfrow = c(1,1))
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################