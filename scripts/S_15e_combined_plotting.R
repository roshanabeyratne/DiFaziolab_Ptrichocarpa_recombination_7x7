# !/usr/bin/R
setwd("/Users/cabeyrat/Google Drive File Stream/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, & DMS for the recombination landscape analysis for 7x7 cross
# DATE: December 18, 2020
# USAGE: Input file include a detailed list of all the cross-over events observed for each focal parent for each offspring for all the chromosomes. The input is produced with script 22b_... 

# Following is the input file format used 
# CHROM  START    END 1863 1909 1950 2048 2066 2283 4593 2365 2393 2515 2572 2683 6909 7073
# Chr01      1   2699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01   2700  32699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01  32700  62699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01  62700  92699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01  92700 122699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01 122700 152699    0    0    0    0    0    0    0    0    0    0    0    0    0    0

# PURPOSE: as supplementary to S_15c and S_15d scripts ... this script produces the figures for the paper.

########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
# library("Rwave");
# library("wavethresh")
# library("haarfisz")
# library("DDHFm")
########################################################################################################################################################################################################


statistics.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/"
permutation.input.path <- "./statistics/heterochiasmy_analysis/wmtsa/"
diff.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/diff_btwn_sexes/"
wv.var.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/wv_var/"
wav.coffs.plot.path <- "./statistics/heterochiasmy_analysis/wmtsa/wav_coffs/"

par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
sliding.avg.output.path <- "./statistics/heterochiasmy_analysis/sliding_average/"
wv.var.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/wv_var/"
DWT.coffs.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/DWT_coffs/"
output.path <- "./statistics/heterochiasmy_analysis/combined_plots/"
########################################################################################################################################################################################################
# Produce the Figure for paper as advised by Steve. 
for (scale in 5) { # scale <- 6
  jpeg(filename = paste0(output.path, "d_", scale, "_sliding_window_average.jpg"), width = 2000, height = 1000, quality = 150)
  op <- par(mfrow=c(4,5), mar=c(2,1,1,1))
  
  for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr19"
    # loading-in required input files
    binned.cos.by.par.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs_by_par.csv"), stringsAsFactors=FALSE, header=TRUE)
    colnames(binned.cos.by.par.df) <- par.name.vec
    actual.binned.COs <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs.csv"), stringsAsFactors=FALSE, header=TRUE)
    binned.permut.cos.fem.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
    binned.permut.cos.male.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
    num.ind.df <- read.csv(paste0(statistics.input.path, "number_of_individuals_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
    colnames(num.ind.df) <- colnames(num.ind.df) %>% str_replace("^X", "")
    
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
  par(op)
  dev.off()
}



# plot the NULL distribution of DWT variance partition vs. actual SEX assignments
for (d in 3:7) { # d <- 5
  jpeg(filename = paste0(output.path, "d_", d, "_WV_variance.jpg"), width = 2000, height = 1000, quality = 150)
  op <- par(mfrow=c(4,5), mar=c(2,1,1,1))
  
  for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
    # read-in the input files
    fem.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
    male.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
    diff.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_DIFF.csv"), stringsAsFactors=FALSE, header=TRUE)
    wv.var.part.df <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_WV_var_for_scales.csv"), stringsAsFactors=FALSE, header=TRUE)
    
    fem.p.val <- sum(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d])> wv.var.part.df[1,d]) / length(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]))
    male.p.val <- sum(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d])> wv.var.part.df[2,d]) / length(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]))
    emp.p.val <- round(1-abs(fem.p.val-male.p.val), digits=2)
    
    hist(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]), breaks=100, main=paste0(chrom.name, " ", colnames(fem.wv.var.part.permut)[d], "; empirical p-value=", emp.p.val), 
         xlab="partitioned variance at the given scale")
    abline(v=wv.var.part.df[1,d], col='red', lwd=3, ylim=c(0, 50))
    abline(v=wv.var.part.df[2,d], col='blue', lwd=3, ylim=c(0, 50))
    
  }
  par(op)
  dev.off()
}



# Plot DWT coefficients for singl scale for all the chromosomes
for (scale in 6) {
  jpeg(filename = paste0(output.path, "d_", scale, "_DWT_coefficients.jpg"), width = 2000, height = 1000, quality = 150)
  op <- par(mfrow=c(4,5), mar=c(2,1,1,1))
  for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[2]
    
    # read-in the input files
    all.par.df.2 <- read.csv(file=paste0(statistics.input.path, "co.bin_size_30kb.csv"), header=TRUE, stringsAsFactors=FALSE )
    colnames(all.par.df.2)[4:17] <- par.name.vec
    test.df <- all.par.df.2[all.par.df.2$CHROM==chrom.name, -(1:3)]
    cos.vec.fem <- apply(test.df[, 1:7], 1 , sum)
    cos.vec.male <- apply(test.df[, 8:14], 1 , sum)
    
    
    permuted.fem.coffs.df <- read.csv(file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
    permuted.male.coffs.df <- read.csv(file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
    permuted.sex.diff.coffs.df <- read.csv(file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_DIFF.csv"), stringsAsFactors=FALSE, header=TRUE)
    actual.DWT.coffs <- read.csv(file=paste0(DWT.coffs.output.path, chrom.name, "_DWT_coffs.csv"), stringsAsFactors=FALSE, header=TRUE)
    
    num.levels <- as.integer(floor(logb(length(cos.vec.fem),base = 2)))+1
    num.coffs <- nrow(actual.DWT.coffs)/num.levels
    
    start.corrd <- (scale*num.coffs)-(num.coffs-1)
    end.coord <- (scale*num.coffs)
    
    for (iter in 1:ncol(permuted.fem.coffs.df)) {
      if (iter == 1) {
        max.y.lim <- max(permuted.fem.coffs.df[start.corrd:end.coord, iter])+1
        min.y.lim <- min(permuted.fem.coffs.df[start.corrd:end.coord, iter])-1
        plot(permuted.fem.coffs.df[start.corrd:end.coord, iter], typ='l', col='grey', ylim=c(min.y.lim, max.y.lim), main=paste0(chrom.name, " d", scale), 
             xlab="overlapping window idx", ylab="mean CO-count")
      } else {
        lines(permuted.fem.coffs.df[start.corrd:end.coord, iter], col='grey')
      }
    }
    for (iter in 1:ncol(permuted.male.coffs.df)) {lines(permuted.male.coffs.df[start.corrd:end.coord, iter], col='grey')}
    
    # plot individual half-sib family counts per parent as well as the mean for each SEX
    lines(actual.DWT.coffs[start.corrd:end.coord, 1], type='l', col='red', lwd=1)
    lines(actual.DWT.coffs[start.corrd:end.coord, 2], col='blue', lwd=1)
  }
  par(op)
  dev.off()
}
