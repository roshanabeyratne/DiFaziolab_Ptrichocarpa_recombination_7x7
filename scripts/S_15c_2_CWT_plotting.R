# !/usr/bin/R

# Written by CRA, & DMS for the recombination landscape analysis for 7x7 cross using CWT
# DATE: January 28, 2022
# PURPOSE: Continuous wavelet transformation of CO data and produce PDFs

# Set path and variables
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
rm(list = ls())
library(wmtsa) # This package has certain requiremtns (ifultools) that are not available on CRAN and it is impossible to install afresh. Therefore I am switching to wavethresh R-package
library(dplyr)
library(RColorBrewer)
# display.brewer.all()



par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- paste0("Chr",sprintf('%02d',1:19))
lg.names <- c("LG-I", "LG-II", "LG-III", "LG-IV", "LG-V", "LG-VI", "LG-VII", "LG-VIII", "LG-IX", "LG-X", "LG-XI", "LG-XII", "LG-XIII", "LG-XIV",
              "LG-XV", "LG-XVI", "LG-XVII", "LG-XVIII", "LG-XIX")
statistics.input.path <- "./statistics/"
cwt.output.path <- "./statistics/heterochiasmy_analysis/wmtsa/CWT/"


# Read-in CO-counts per half-sib family for each non-overlapping window of size 30kb. 
all.par.df.2 <- read.csv(file=paste0(statistics.input.path, "co.bin_size_30kb_Stet_v2.csv"), header=TRUE, stringsAsFactors=FALSE )
colnames(all.par.df.2)[4:17] <- par.name.vec

# checking whether the total number of cross-over counts are the same as observed
sum(apply(all.par.df.2[,4:17], 1, sum)) # 38846. This checks out

# svg(file=paste0(cwt.output.path, "/all_CWT_2.svg"), width=8, height=10)
pdf(file=paste0(cwt.output.path, "/all_CWT_14_19.pdf"), width=8, height=11)
# layout(matrix(c(1,2,3,4,3,4,
#                 5,6,7,8,7,8,
#                 9,10,11,12,11,12), ncol=2, byrow=T))
layout(matrix(c(1,2,3,4,3,4,
                5,6,7,8,7,8,
                9,10,11,12,11,12,
                13,14,15,16,15,16,
                17,18,19,20,19,20,
                21,22,23,24,23,24), ncol=2, byrow=T))
# layout(matrix(c(1,2,3,4,3,4), ncol=2, byrow=T)
# layout.show(n=40)
for (chrom.name in chrom.name.vec[14:19]) {
  # chrom.name <- chrom.name.vec[1]
  lg_name <- lg.names[which(chrom.name.vec==chrom.name)]
  test.df <- all.par.df.2[all.par.df.2$CHROM==chrom.name, -(1:3)]
  cos.vec.fem <- apply(test.df[, 1:7], 1 , sum)
  cos.vec.male <- apply(test.df[, 8:14], 1 , sum)
  cos.vec.sex.diff <- (cos.vec.fem-cos.vec.male) # this is not a required parameter for downstream analysis.
  # males have 298 mre CO events than females in out dataset.
  
  
  # The CWT transformation and plotting of the coefficients using a 'Mexican-Hat wavelet'
  wavCWT.obj.fem <- wavCWT(cos.vec.fem, scale.range=deltat(cos.vec.fem)*c(1, length(cos.vec.fem)), n.scale=1000, wavelet="gaussian2", shift=5, variance=1)
  wavCWT.obj.male <- wavCWT(cos.vec.male, scale.range=deltat(cos.vec.male)*c(1, length(cos.vec.male)), n.scale=1000, wavelet="gaussian2", shift=5, variance=1)
  wavCWT.obj.sex.diff <- wavCWT(cos.vec.sex.diff, scale.range=deltat(cos.vec.male)*c(1, length(cos.vec.male)), n.scale=1000, wavelet="gaussian2", shift=5, variance=1)
  
  # obtain CWT coefficients
  cwt.fem.coffs <- as.matrix(wavCWT.obj.fem)
  cwt.male.coffs <- as.matrix(wavCWT.obj.male)
  cwt.sex.diff.coffs <- as.matrix(wavCWT.obj.sex.diff)
  
  # CWT.coffs.fem<-apply(as.matrix(wavCWT.obj.fem), 2, norm)
  # CWT.coffs.male<-apply(as.matrix(wavCWT.obj.fem), 2, norm)
  
  min.fem <- min(abs(as.numeric(cwt.fem.coffs)))
  min.male <- min(abs(as.numeric(cwt.male.coffs)))
  
  # plot(wavCWT.obj.male, col=blue.colors)
  # image(1:nrow(cwt.fem.coffs), log(1:ncol(cwt.fem.coffs), base=2), cwt.fem.coffs^0.5)
  # image(1:nrow(cwt.sex.diff.coffs), log(1:ncol(cwt.sex.diff.coffs), base=2), cwt.sex.diff.coffs)
  blue.colors <- brewer.pal(8, "Blues")
  
  op <- par(mar=c(0, 4.1, 2.1, 2.1), cex=0.5)
  y_lim_max <- (max(c(cos.vec.fem, cos.vec.male)))
  plot(cos.vec.fem, type='l', xlab=NULL, ylab='CO count', xaxt='n', main=paste0(lg_name, ' Female'), ylim=c(0, y_lim_max))
  plot(cos.vec.male, type='l', xlab=NULL, ylab='CO count',xaxt='n', main=paste0(lg_name, ' Male'), ylim=c(0, y_lim_max))
  par(op)
  # cwt.fem.coffs_t <- apply(cwt.fem.coffs, 2, rev)
  # cwt.male.coffs_t <- apply(cwt.male.coffs, 2, rev)
  # dev.new()
  op <- par(mar=c(5.1, 4.1, 0, 2.1), cex=0.5)
  f<-t(abs(cwt.fem.coffs)^0.25)
  f <- f[nrow(f):1,]
  f<-t(f)
  image(x=1:nrow(cwt.fem.coffs), y=-log(ncol(cwt.fem.coffs):1, base=2), z=f, col=blue.colors, yaxt='n', ylab='',xlab='30kb genomic window idx')
  labs=(c('30kb', '60kb', '120kb', '240kb', '480kb', '960kb', '1.9Mb', '3.85Mb', '7.68MB'))
  axis(side=2, at=-c(0:8), labels=labs, las=2)
  m<-t(abs(cwt.male.coffs)^0.25)
  m <- m[nrow(m):1,]
  m<-t(m)
  image(x=1:nrow(cwt.male.coffs), y=-log(ncol(cwt.male.coffs):1, base=2), z=m, col=blue.colors, yaxt='n', ylab='',xlab='30kb genomic window idx')
  # image(x=1:nrow(cwt.fem.coffs), y=log(1:ncol(cwt.fem.coffs), base=2), z=abs(cwt.fem.coffs)^0.25, col=blue.colors, yaxt='n', ylab='',xlab='30kb genomic window idx')
  # image(x=1:nrow(cwt.male.coffs), y=log(1:ncol(cwt.male.coffs), base=2), z=abs(cwt.male.coffs)^0.25, col=blue.colors, yaxt='n', ylab='',xlab='30kb genomic window idx')

  axis(side=2, at=-(c(0:8)), labels=labs, las=2)
  par(op)
}
dev.off()









