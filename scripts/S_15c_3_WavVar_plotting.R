# !/usr/bin/R

# Written by CRA, & DMS for the recombination landscape analysis for 7x7 cross using CWT
# DATE: January 28, 2022
# PURPOSE: Plot wavelet variance partitioning for selected chromosomes at all scales 1 through 7 and produce PDF

setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
rm(list = ls())
library(wmtsa) 


# Set variabele and input output paths
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
statistics.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/"
permutation.input.path <- "./statistics/heterochiasmy_analysis/wmtsa/"
permutation.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/permutation_tabs/"
DWT.coffs.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/DWT_coffs/"
diff.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/diff_btwn_sexes/"
wv.var.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/wv_var/"
wav.coffs.plot.path <- "./statistics/heterochiasmy_analysis/wmtsa/wav_coffs/"
sliding.avg.output.path <- "./statistics/heterochiasmy_analysis/sliding_average/"
cwt.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/CWT/"



pdf(file=paste0(wv.var.output.path, "Chr1_19_8_11_WV_variance.pdf"), width=10, height=16)
layout(matrix(1:28, nrow=7, byrow=F))
# layout.show(n=28)
female_p_val_df <- NULL
male_p_val_df <- NULL
for (chrom.name in chrom.name.vec[c(1,14,11,13)]) { 
  fem_male_p_val <- NULL
  for (d in 1:7) {
    
    fem.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
    male.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
    diff.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_DIFF.csv"), stringsAsFactors=FALSE, header=TRUE)
    wv.var.part.df <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_WV_var_for_scales.csv"), stringsAsFactors=FALSE, header=TRUE)
    
    fem.p.val <- sum(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d])> wv.var.part.df[1,d]) / length(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]))
    male.p.val <- sum(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d])> wv.var.part.df[2,d]) / length(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]))
    fem_male_p_val <- rbind(fem_male_p_val, c(fem.p.val, male.p.val))
    
    
    hist(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]), breaks=100, main=paste0(chrom.name, " ", colnames(fem.wv.var.part.permut)[d], "; fem.p=", fem.p.val, " male.p.val=", male.p.val),
         xlab="partitioned variance at the given scale",
         ylim=c(0,80),
         lty='blank', col='grey')
    abline(v=wv.var.part.df[1,d], col='red', lwd=3, ylim=c(0, 50))
    abline(v=wv.var.part.df[2,d], col='blue', lwd=3, ylim=c(0, 50))
  }
  female_p_val_df <- cbind(female_p_val_df, fem_male_p_val[,1])
  male_p_val_df <- cbind(male_p_val_df, fem_male_p_val[,2])
}
dev.off()
