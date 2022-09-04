# !/usr/bin/R
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, & DMS for the recombination landscape analysis for 7x7 cross
# DATE: October 30, 2020; modified: December 02, 2020
# USAGE: Input file include a detailed list of all the cross-over events observed for each focal parent for each offspring for all the chromosomes. The input is produced with script 22b_... 

# Following is the input file format used 
# CHROM  START    END 1863 1909 1950 2048 2066 2283 4593 2365 2393 2515 2572 2683 6909 7073
# Chr01      1   2699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01   2700  32699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01  32700  62699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01  62700  92699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01  92700 122699    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# Chr01 122700 152699    0    0    0    0    0    0    0    0    0    0    0    0    0    0

# PURPOSE: Analyze COs across chromosomes with a continuous wavelet decomposition to visualize localized signals with our cross-over data. This analysis also include a discrete wavelet transform for 
# conucting statistical tests and effect analysis on multiple different scales within a chromosome that are orthogonal to each other. 

########################################################################################################################################################################################################
rm(list = ls())
library(wmtsa) # This package has certain requiremtns (ifultools) that are not available on CRAN and it is impossible to install afresh. Therefore I am switching to wavethresh R-package
library(dplyr)
library(stringr)
# library("Rwave");
# library("wavethresh")
# library("haarfisz")
# library("DDHFm")
########################################################################################################################################################################################################
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
########################################################################################################################################################################################################
# This is a code snippet for demonstrating how DWT works. 
# linchirp <- make.signal("bumps", n=1024)
# plot(wavShift(wavMODWT(linchirp, wavelet="s8", n.levels=10, keep.series=TRUE)))
# eda.plot((wavMODWT(linchirp, wavelet="s8", n.levels=10, keep.series=TRUE)))
# 
# # Different types of signals produced with wmtsa R-package
# op <- par(mfrow=c(1,1))
# nms <- c("blocks", "linchirp", "mishmash1", "bumps")
# z <- lapply(nms, make.signal)
# ifultools::stackPlot(x=seq(1024),y=z, ylab=nms)
# par(op)
########################################################################################################################################################################################################

# Read-in CO-counts per half-sib family for each non-overlapping window of size 30kb. 
# all.par.df.2 <- read.csv(file=paste0(statistics.input.path, "co.bin_size_30kb.csv"), header=TRUE, stringsAsFactors=FALSE )
all.par.df.2 <- read.csv(file=paste0(statistics.input.path, "co.bin_size_30kb_Stet_v2.csv"), header=TRUE, stringsAsFactors=FALSE )
colnames(all.par.df.2)[4:17] <- par.name.vec

# checking whether the total number of cross-over counts are the same as observed
sum(apply(all.par.df.2[,4:17], 1, sum)) # 38846. This checks out



# Mean cross-over count per individual at the whole chromosome level is higher for males in almost all the chromosomes. 
# The next step is to investigate whether there are difference in CO-counts between sexes at finer scales in a range of different scales.
for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr11"
  
  # test.df <- all.par.df.2[all.par.df.2$CHROM==chrom.name, -(1:3)]
  # cos.vec.fem <- apply(test.df[, 1:7], 1 , sum)
  # cos.vec.male <- apply(test.df[, 8:14], 1 , sum)
  # cos.vec.sex.diff <- (cos.vec.fem-cos.vec.male) # this is not a required parameter for downstream analysis. 
  # # males have 298 mre CO events than females in out dataset.
  
  
  # A tests for a signal with uneven mass at a certain region. This is only to illustrate how detection at different scales change with DWT 
  # set.seed(1960)
  # test.vec <- rnorm(1024)
  # test.vec <- rpois(1024, 3)
  # test.vec[100:125] <-  test.vec[100:125] + 10 # rchisq(length(100:150), df=1)
  # test.vec[500:525] <-  test.vec[500:525] + 10 # rchisq(length(100:150), df=1)
  # plot(test.vec, type='l')
  # cos.vec.sex.diff <- test.vec

  
  # PERMUTATION PIPELINE
  
  # # Permute SEX assignment without disturbing localization on Chromosome. 
  # set.seed(10)
  # permut.fem.df <- NULL
  # permut.male.df <- NULL
  # for (iter in 1:1000) {
  #   test.df.permut <- as.data.frame(t(apply(test.df, 1, function(x){sample(x, size=length(x), replace=FALSE)})), stringsAsFactors=FALSE )
  #   temp.fem.permut <- apply(test.df.permut[, 1:7], 1 , sum)
  #   temp.male.permut <- apply(test.df.permut[, 8:14], 1 , sum)
  # 
  #   permut.fem.df <- cbind(permut.fem.df, temp.fem.permut)
  #   permut.male.df <- cbind(permut.male.df, temp.male.permut)
  #   cat(iter, " ")
  # }
  # permut.sex.diff.df <- (permut.fem.df-permut.male.df)

  # # Permutation: sex assignment is permuted for parents in the half-sib family. USE this
  # set.seed(10)
  # permut.fem.df <- NULL
  # permut.male.df <- NULL
  # for (iter in 1:1000) {
  #   new.col.order <- sample(1:14, replace=FALSE)
  #   test.df.permut <- test.df[,new.col.order]
  #   temp.fem.permut <- apply(test.df.permut[, 1:7], 1 , sum)
  #   temp.male.permut <- apply(test.df.permut[, 8:14], 1 , sum)
  # 
  #   permut.fem.df <- cbind(permut.fem.df, temp.fem.permut)
  #   permut.male.df <- cbind(permut.male.df, temp.male.permut)
  #   cat(iter, " ")
  # }
  # colnames(permut.fem.df) <- paste0("P_", 1:1000)
  # colnames(permut.male.df) <- paste0("P_", 1:1000)
  # permut.sex.diff.df <- (permut.fem.df-permut.male.df)
  # write.csv(permut.fem.df, file=paste0(permutation.output.path, chrom.name, "_fem.permut_df.csv"), quote=FALSE, row.names=FALSE)
  # write.csv(permut.male.df, file=paste0(permutation.output.path, chrom.name, "_male.permut_df.csv"), quote=FALSE, row.names=FALSE)
  permut.fem.df <- read.csv(paste0(permutation.output.path, chrom.name, "_fem.permut_df.csv"), stringsAsFactors=FALSE, header=TRUE )
  permut.male.df <- read.csv(paste0(permutation.output.path, chrom.name, "_male.permut_df.csv"), stringsAsFactors=FALSE, header=TRUE )
  permut.sex.diff.df <- (permut.fem.df-permut.male.df) 
  
  
  
  # Pearson's product-moment correlation between combined male and female data-sets on the scale of 30kb window size
  # cor.test(cos.vec.fem, cos.vec.male) # (r^2=0.2798432; p<2.2e-16)
  # permuted.sex.corr.vec <- NULL
  # for (i in 1:1000) { # i<-1
  #   temp.cor <- cor.test(permut.fem.df[,i], permut.male.df[,i])
  #   permuted.sex.corr.vec <- c(permuted.sex.corr.vec, temp.cor$estimate)
  # }
  # hist(permuted.sex.corr.vec, breaks=10, main="R^2 of vectors when sex is permuted", xlab="R^2") # (r^2=0.3115163; p<2.2e-16)
  # sum(permuted.sex.corr.vec < 0.2798432)/length(permuted.sex.corr.vec) # p=0.007; for observing a value as extreme
  # NOTE: There seems to be more correlation between groups of individuals than the when they are assigned to actual male and female groups. 
  
  
  
  # Following code-block is not required for analysis with wmtsa R-package. padding of datasets is required for other wavelet analysis pipe-lines, so that total length of time-series 
  # is an integer power of base 2. 
  # since downstream analysis requires that the total size of the signal to be a power of 2. here I am adding NA padding based on the vector size
  # vec.size<-length(cos.vec.male)
  # exp.value<-ceiling(log(length(cos.vec.male), base = 2))
  # desired.vec.size<-2^exp.value
  # upstrm.pad.size<-ceiling((desired.vec.size-vec.size)/2)
  # dwnstrm.pad.size<-floor((desired.vec.size-vec.size)/2)
  # 
  # if(vec.size+upstrm.pad.size+dwnstrm.pad.size!=desired.vec.size) {
  #   stop("error in the desired length of the final vector to be analyzed!")
  # } else {
  #   cos.vec.male<-c(rep(0,upstrm.pad.size), cos.vec.male, rep(0,dwnstrm.pad.size))
  #   cos.vec.fem<-c(rep(0,upstrm.pad.size), cos.vec.fem, rep(0,dwnstrm.pad.size))
  # }
  
  
  
  # # Discrete-wavelet transformation using wmtsa R-package. 
  # # Produce DWT on total CO-count data for each sex, using a Haar wavelet. 
  # wavMODWT.obj.fem <- wavMODWT(cos.vec.fem, wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))), 
  #                              position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())
  # wavMODWT.obj.male <- wavMODWT(cos.vec.male, wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))), 
  #                               position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())
  # wavMODWT.obj.sex.diff <- wavMODWT(cos.vec.sex.diff, wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.sex.diff),base = 2))), 
  #                                   position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())
  
  # plotting just the orthogonal breakdown of the signal 
  # plot(wavMODWT.obj.fem)
  # plot(wavMODWT.obj.male)
  # plot(wavMODWT.obj.sex.diff)
  # plot the signal break-down as wella s variance at each scale
  # eda.plot(wavMODWT.obj.fem)
  # eda.plot(wavMODWT.obj.male)
  # eda.plot(wavMODWT.obj.sex.diff)
  # fem.coffs <- as.matrix(wavMODWT.obj.fem)
  # male.coffs <- as.matrix(wavMODWT.obj.male)
  # sex.diff.coffs <- as.matrix(wavMODWT.obj.sex.diff)
  # actual.DWT.coffs <- cbind(fem.coffs, male.coffs, sex.diff.coffs)
  # write.csv(actual.DWT.coffs, file=paste0(DWT.coffs.output.path, chrom.name, "_DWT_coffs.csv"), quote=FALSE, row.names=FALSE)
  actual.DWT.coffs <- read.csv(file=paste0(DWT.coffs.output.path, chrom.name, "_DWT_coffs.csv"), stringsAsFactors=FALSE, header=TRUE)
  
  
  
  # DWT for SEX permuted artificial data-set that was created. The assignment of SEX in this data-set is just symbolic.
  # permuted.fem.coffs.df <- NULL
  # permuted.male.coffs.df <- NULL
  # permuted.sex.diff.coffs.df <- NULL
  # for (iter in 1:1000) {
  #   wavMODWT.obj.fem.permut <- wavMODWT(permut.fem.df[,iter], wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))),
  #                                       position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())
  #   wavMODWT.obj.male.permut <- wavMODWT(permut.male.df[,iter], wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))),
  #                                        position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())
  #   wavMODWT.obj.sex.diff.permut <- wavMODWT(permut.sex.diff.df[, iter], wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))),
  #                                            position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())
  # 
  #   permuted.fem.coffs.df <- cbind(permuted.fem.coffs.df, as.matrix(wavMODWT.obj.fem.permut))
  #   permuted.male.coffs.df <- cbind(permuted.male.coffs.df, as.matrix(wavMODWT.obj.male.permut))
  #   permuted.sex.diff.coffs.df <- cbind(permuted.sex.diff.coffs.df, as.matrix(wavMODWT.obj.sex.diff.permut))
  #   cat(iter, " ")
  # }
  # write.csv(permuted.fem.coffs.df, file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_FEM.csv"), quote=FALSE, row.names=FALSE)
  # write.csv(permuted.male.coffs.df, file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_MALE.csv"), quote=FALSE, row.names=FALSE)
  # write.csv(permuted.sex.diff.coffs.df, file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_DIFF.csv"), quote=FALSE, row.names=FALSE)
  permuted.fem.coffs.df <- read.csv(file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
  permuted.male.coffs.df <- read.csv(file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
  # permuted.sex.diff.coffs.df <- read.csv(file=paste0(DWT.coffs.output.path, chrom.name, "_permuted_DWT_coffs_DIFF.csv"), stringsAsFactors=FALSE, header=TRUE)
  
  
  
  # plot the DWT wavelet coefficients for the difference. Not sure whether this block of code is the right analysis. Will be purged after review.
  # op <- par(mfrow=c(1,3))
  # hist(cos.vec.sex.diff, breaks=20, main="sex difference", xlab="CWT coefficients", col='lightgrey')
  # hist(as.vector(permuted.sex.diff.coffs.df), breaks=20, main="permuted sex difference", xlab="CWT coefficients", col='darkgrey')
  # par(op)
  # dev.off()

  
  
  # Extracting DWT-coefficients for each scale.
  num.levels <- as.integer(floor(logb(length(cos.vec.fem),base = 2)))+1
  num.coffs <- length(fem.coffs)/num.levels
  # d1.coffs<-which(grepl("d1\\(\\d+\\)", row.names(fem.coffs)))
  # d2.coffs<-which(grepl("d2\\(\\d+\\)", row.names(fem.coffs)))
  # plot DWT wavelet coefficients at each scale
  jpeg(filename = paste0(wav.coffs.plot.path, chrom.name, "_DWT_coefficients.jpg"), width = 2000, height = 1000, quality = 150)
  op <- par(mfrow=c(3,4))
  for (i in 1:num.levels) {
    start.corrd <- (i*num.coffs)-(num.coffs-1)
    end.coord <- (i*num.coffs)
    
    for (iter in 1:ncol(permuted.fem.coffs.df)) {
      if (iter == 1) {
        max.y.lim <- max(permuted.fem.coffs.df[start.corrd:end.coord, iter])+1
        min.y.lim <- min(permuted.fem.coffs.df[start.corrd:end.coord, iter])-1
        plot(permuted.fem.coffs.df[start.corrd:end.coord, iter], typ='l', col='grey', ylim=c(min.y.lim, max.y.lim), main=paste0(chrom.name, " d", i), 
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
  
  
  
  # calculating the correlation between male and female coefficients at each scale, for actual and one of the instances of permutation. and doing the same for the permuted vector.
  # However, in general there is high correlation between two random sets of individuals at each scale owing to highly shared recombination landscape. it is difficult to parse SEX based 
  # differences with this method. 
  # permut.fem.coffs <- as.vector(permuted.fem.coffs.df[,1])
  # permut.male.coffs <- as.vector(permuted.male.coffs.df[,1])
  # 
  # corr.df <- NULL
  # permu.corr.df <- NULL
  # for (i in 1:num.levels) {
  #   start.corrd <- (i*num.coffs)-(num.coffs-1)
  #   end.coord <- (i*num.coffs)
  #   
  #   temp.corr <- cor.test(fem.coffs[start.corrd:end.coord], male.coffs[start.corrd:end.coord])
  #   corr.df <- rbind(corr.df, c(temp.corr$estimate, temp.corr$conf.int[1], temp.corr$conf.int[2], temp.corr$p.value))
  #   
  #   temp.corr.permut <- cor.test(permut.fem.coffs[start.corrd:end.coord], permut.male.coffs[start.corrd:end.coord])
  #   permu.corr.df <- rbind(permu.corr.df, c(temp.corr.permut$estimate, temp.corr.permut$conf.int[1], temp.corr.permut$conf.int[2], temp.corr.permut$p.value))
  #   
  # }
  # corr.df <- as.data.frame(cbind(c(paste0("d", 1:10), "s10"), corr.df), stringsAsFactors=FALSE)
  # permu.corr.df <- as.data.frame(cbind(c(paste0("d", 1:10), "s10"), permu.corr.df), tringsAsFactors=FALSE)
  # 
  # colnames(corr.df) <- c("SCALE", "R_SQR", "LOW_CONFINT", "HI_CONFINT", "P_VAL")
  # colnames(permu.corr.df) <- c("SCALE", "R_SQR", "LOW_CONFINT", "HI_CONFINT", "P_VAL")
  # 
  # for (j in 2:ncol(corr.df)) {corr.df[,j] <- as.numeric(corr.df[,j])}
  # for (j in 2:ncol(permu.corr.df)) {permu.corr.df[,j] <- as.numeric(permu.corr.df[,j])}
  # 
  # # plot the correlation coefficients
  # plot(1:11, corr.df$R_SQR, pch=22, col='blue', xaxt = 'n', ylim=c(0, 1), xlab="scale", ylab="coefficient of determination", main="Correlation between males & females", frame.plot=FALSE)
  # axis(side=1, at=1:11, labels=corr.df$SCALE)
  # arrows(1:11, corr.df$LOW_CONFINT, 1:11, corr.df$HI_CONFINT, code=3, angle=90, length=0.1)
  
  
  
  # Proportion of variance at each sclae of the DWT transformed signal can be calculated using the DWT coefficients. The amount of variance at each level as a proportion of the total can then be
  # calculated. In this analysis the varinace is always highest at smaller scales. However, in order to understand which scale is most responsive to the SEX permutation 
  # (i.e at which scale are the SEX difference in the signal more prominent), we obtained the the NULL distribution of variance partitioning decoupled from SEX. SEX assignment is permuted and 
  # therefore the sex based clusters are just nominal. 
  fem.wv.var.part.permut <- NULL
  male.wv.var.part.permut <- NULL
  diff.wv.var.part.permut <- NULL
  for (iter in 1:1000) {
    tmp.obj.fem <- wavVar(permut.fem.df[,iter], xform="modwt", wavelet="haar", n.levels=NULL,
                        position=list(from=1,by=1,units=character()), units=character(),
                        documentation=character(), sdf=NULL, sdfargs=NULL,
                        sampling.interval=1, n.fft=1024)
    fem.wv.var.part.permut <- rbind(fem.wv.var.part.permut, tmp.obj.fem[[2]]$unbiased)


    tmp.obj.male <- wavVar(permut.male.df[,iter], xform="modwt", wavelet="haar", n.levels=NULL,
                        position=list(from=1,by=1,units=character()), units=character(),
                        documentation=character(), sdf=NULL, sdfargs=NULL,
                        sampling.interval=1, n.fft=1024)
    male.wv.var.part.permut <- rbind(male.wv.var.part.permut, tmp.obj.male[[2]]$unbiased)

    tmp.obj.diff <- wavVar((permut.fem.df[,iter]-permut.male.df[,iter]), xform="modwt", wavelet="haar", n.levels=NULL,
                           position=list(from=1,by=1,units=character()), units=character(),
                           documentation=character(), sdf=NULL, sdfargs=NULL,
                           sampling.interval=1, n.fft=1024)
    diff.wv.var.part.permut <- rbind(diff.wv.var.part.permut, tmp.obj.diff[[2]]$unbiased)
    cat(iter, " ")
  }
  write.csv(fem.wv.var.part.permut, file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_FEM.csv"), quote=FALSE, row.names=FALSE)
  write.csv(male.wv.var.part.permut, file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_MALE.csv"), quote=FALSE, row.names=FALSE)
  write.csv(diff.wv.var.part.permut, file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_DIFF.csv"), quote=FALSE, row.names=FALSE)
  fem.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
  male.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
  diff.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_DIFF.csv"), stringsAsFactors=FALSE, header=TRUE)
  
  
  
  # DWT-variance partitioning for observed data of FEMALE and MALEs
  wav.var.obj.fem <- wavVar(cos.vec.fem, xform="modwt", wavelet="haar", n.levels=NULL,
                        position=list(from=1,by=1,units=character()), units=character(),
                        documentation=character(), sdf=NULL, sdfargs=NULL,
                        sampling.interval=1, n.fft=1024)
  wav.var.obj.male <- wavVar(cos.vec.male, xform="modwt", wavelet="haar", n.levels=NULL,
                            position=list(from=1,by=1,units=character()), units=character(),
                            documentation=character(), sdf=NULL, sdfargs=NULL,
                            sampling.interval=1, n.fft=1024)
  wav.var.obj.diff <- wavVar((cos.vec.fem-cos.vec.male), xform="modwt", wavelet="haar", n.levels=NULL,
                             position=list(from=1,by=1,units=character()), units=character(),
                             documentation=character(), sdf=NULL, sdfargs=NULL,
                             sampling.interval=1, n.fft=1024)
  wv.var.part.df <- as.data.frame(rbind(wav.var.obj.fem[[2]]$unbiased, wav.var.obj.male[[2]]$unbiased, wav.var.obj.diff[[2]]$unbiased), stringsAsFactors=FALSE)
  write.csv(wv.var.part.df, file=paste0(wv.var.output.path, chrom.name, "_WV_var_for_scales.csv"), quote=FALSE, row.names=FALSE)
  wv.var.part.df <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_WV_var_for_scales.csv"), stringsAsFactors=FALSE, header=TRUE)
  
  # plot the NULL distribution of variance vs. actual
  jpeg(filename = paste0(wv.var.output.path, chrom.name, "_WV_variance.jpg"), width = 2000, height = 1000, quality = 150)
  op <- par(mfrow=c(4,3), mar=c(5.1, 4.1, 4.1, 2.1))
  for (d in 1:ncol(wv.var.part.df)) {
    fem.p.val <- sum(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d])> wv.var.part.df[1,d]) / length(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]))
    male.p.val <- sum(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d])> wv.var.part.df[2,d]) / length(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]))

    hist(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]), breaks=100, main=paste0(chrom.name, " ", colnames(fem.wv.var.part.permut)[d], "; fem.p=", fem.p.val, " male.p.val=", male.p.val),
         xlab="partitioned variance at the given scale")
    abline(v=wv.var.part.df[1,d], col='red', lwd=3, ylim=c(0, 50))
    abline(v=wv.var.part.df[2,d], col='blue', lwd=3, ylim=c(0, 50))
  }
  par(op)
  dev.off()
  
  # # Plot the emperical NULL vs. actual DWT warince at each scale for the difference of total CO counts between sexes
  # jpeg(filename = paste0(diff.output.path, chrom.name, "_WV_var_for_diff.jpg"), width = 2000, height = 1000, quality = 150)
  # op <- par(mfrow=c(4,3), mar=c(5.1, 4.1, 4.1, 2.1))
  # for (d in 1:ncol(wv.var.part.df)) {
  #   diff.p.val <- sum(diff.wv.var.part.permut[,d] > wv.var.part.df[3,d]) / length(diff.wv.var.part.permut[,d])
  #   hist(diff.wv.var.part.permut[,d], breaks=100, main=paste0(chrom.name, " ", colnames(fem.wv.var.part.permut)[d], "; diff.p=", diff.p.val), 
  #        xlab="partitioned variance at the given scale")
  #   abline(v=wv.var.part.df[3,d], col='black', lwd=3)
  # }
  # par(op)
  # dev.off()

  cat(chrom.name, "\n")
}



# plot the NULL distribution of DWT variance partition vs. actual SEX assignments. For each scale for all chromosomes. one p-value generated
for (d in 1:7) { # d <- 5
  jpeg(filename = paste0(wv.var.output.path, "d_", d, "_WV_variance.jpg"), width = 2000, height = 1000, quality = 150)
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




# This is the plotting of all the chromosomes at all scales 1 through 7 to observe patterns
# jpeg(filename = paste0(wv.var.output.path, "Chr01_10", "_WV_variance.jpg"), width = 2000, height = 1000, quality = 150)
jpeg(filename = paste0(wv.var.output.path, "Chr11_19", "_WV_variance.jpg"), width = 2000, height = 1000, quality = 150)
op <- par(mfrow=c(10,7), mar=c(1,1,3.1,1))
# for (chrom.name in chrom.name.vec[1:10]) {
female_p_val_df <- NULL
male_p_val_df <- NULL
for (chrom.name in chrom.name.vec[11:19]) { 
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
         xlab="partitioned variance at the given scale", lty='blank', col='grey')
    abline(v=wv.var.part.df[1,d], col='red', lwd=3, ylim=c(0, 50))
    abline(v=wv.var.part.df[2,d], col='blue', lwd=3, ylim=c(0, 50))
  }
  female_p_val_df <- cbind(female_p_val_df, fem_male_p_val[,1])
  male_p_val_df <- cbind(male_p_val_df, fem_male_p_val[,2])
}
par(op)
dev.off()

# Checking which scale is best
apply(male_p_val_df, 1, function(x){sum(x<=0.10 | x>=0.9)})
apply(female_p_val_df, 1, function(x){sum(x<=0.10 | x>=0.9)})
test_threshs <- c(0.4, 0.3, 0.2, 0.1, 0.05, 0.025)
for (p in test_threshs) {
  lty_p <- as.numeric(which(test_threshs==p))
  low_p <- p
  high_p <- 1-low_p
  tmp_vec <- (apply(male_p_val_df, 1, function(x){sum(x<=low_p | x>=high_p)}) + apply(female_p_val_df, 1, function(x){sum(x<=low_p | x>=high_p)}))/2 # decided to go with this
  if (p==test_threshs[1]) {
    plot(tmp_vec, type='l', ylim=c(0, 20), lty=lty_p, xlab='scale (d)', ylab='number of chromosomes showing either male or female bias')
  } else {
    lines(tmp_vec, lty=lty_p)
  }
}
# 
# apply(male_p_val_df, 1, function(x){sum(x<=0.10 | x>=0.9)}) + apply(female_p_val_df, 1, function(x){sum(x<=0.10 | x>=0.9)}) # decided to go with this
# apply(male_p_val_df, 1, function(x){sum(x<=0.10 | x>=0.9)}) + apply(female_p_val_df, 1, function(x){sum(x<=0.10 | x>=0.9)}) # decided to go with this
# apply(male_p_val_df, 1, function(x){sum(x<=0.05 | x>=0.95)}) + apply(female_p_val_df, 1, function(x){sum(x<=0.05 | x>=0.95)})
# apply(male_p_val_df, 1, function(x){sum(x<=0.025 | x>=0.975)}) + apply(female_p_val_df, 1, function(x){sum(x<=0.025 | x>=0.975)})



# # Just checking the differences in variance
# op <- par(mfrow=c(10,7), mar=c(1,1,3.1,1))
# # for (chrom.name in chrom.name.vec[1:10]) {
# for (chrom.name in chrom.name.vec[11:19]) {
#   for (d in 1:7) {
# 
#     fem.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
#     male.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
#     diff.wv.var.part.permut <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_permuted_WV_var_DIFF.csv"), stringsAsFactors=FALSE, header=TRUE)
#     wv.var.part.df <- read.csv(file=paste0(wv.var.output.path, chrom.name, "_WV_var_for_scales.csv"), stringsAsFactors=FALSE, header=TRUE)
# 
#     p.val <- sum(diff.wv.var.part.permut[,d]> wv.var.part.df[3,d]) / length(diff.wv.var.part.permut[,d])
# 
#     # hist((sqrt(male.wv.var.part.permut[,d]) - sqrt(fem.wv.var.part.permut[,d])), breaks=100, main=paste0(chrom.name, " ", colnames(fem.wv.var.part.permut)[d], "; fem.p=", fem.p.val, " male.p.val=", male.p.val),
#     #      xlab="partitioned variance at the given scale", lty='blank', col='grey')
#     hist(diff.wv.var.part.permut[,d], breaks=100, main=paste0(chrom.name, " ", colnames(diff.wv.var.part.permut)[d], "; p=", p.val),
#          xlab="partitioned variance at the given scale", lty='blank', col='grey')
#     abline(v=wv.var.part.df[3,d], col='red', lwd=3, ylim=c(0, 50))
#   }
# }
par(op)

########################################################################################################################################################################################################
# END of DWT analysis 
########################################################################################################################################################################################################
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

library(RColorBrewer)
display.brewer.all()

jpeg(filename = paste0(cwt.output.path, chrom.name, "_CWT_male.jpg"), width = 2000, height = 1000, quality = 150)
# plot(wavCWT.obj.male)
# image(1:nrow(cwt.fem.coffs), log(1:ncol(cwt.fem.coffs), base=2), cwt.fem.coffs^0.5)
# image(1:nrow(cwt.sex.diff.coffs), log(1:ncol(cwt.sex.diff.coffs), base=2), cwt.sex.diff.coffs)
blue.colors <- brewer.pal(8, "Blues")

dev.new()
op <- par(mfrow=c(1,2))
image(1:nrow(cwt.fem.coffs), log(1:ncol(cwt.fem.coffs), base=2), abs(cwt.fem.coffs)^0.25, col=blue.colors)
image(1:nrow(cwt.male.coffs), log(1:ncol(cwt.male.coffs), base=2), abs(cwt.male.coffs)^0.25, col=blue.colors)
par(op)
dev.off()

op <- par(mfrow=c(1,2))
plot(cos.vec.fem, type='l')
plot(cos.vec.male, type='l')
par(op)












# Animation of scale based variation of wavelet coefficients
for (iter in 1:ncol(cwt.fem.coffs)) { plot(cwt.fem.coffs[, iter], main=paste0("log scale =", log(iter)), col='red'); points(cwt.male.coffs[, iter], col='blue')}


# define a normalizing function
# norm <- function(data) {
#   norm <- (data-mean(data))/sd(data)
#   return(norm)
# }
# cwt.sex.diff.coffs.normed <- t(apply(cwt.sex.diff.coffs, 1, norm))



#  CWT for SEX permuted vectors for 100 iterations
permuted.cwt.fem.coffs.list <- list()
permuted.cwt.male.coffs.list <- list()
permuted.cwt.sex.diff.coffs.list <- list()

scale=100
for (iter in 1:100) { # iter<-1
 
  wavCWT.obj.fem.permut <- wavCWT(permut.fem.df[,iter], scale.range=deltat(cos.vec.fem)*c(1, length(cos.vec.fem)), n.scale=1000, wavelet="haar", shift=5, variance=1)
  wavCWT.obj.male.permut <- wavCWT(permut.male.df[,iter], scale.range=deltat(cos.vec.male)*c(1, length(cos.vec.male)), n.scale=1000, wavelet="haar", shift=5, variance=1)
  # wavCWT.obj.sex.diff.permut <- wavCWT(permut.sex.diff.df[,iter], scale.range=deltat(cos.vec.male)*c(1, length(cos.vec.male)), n.scale=1000, wavelet="gaussian2", shift=5, variance=1)
  
  permut.cwt.fem.coffs <- as.matrix(wavCWT.obj.fem.permut)
  permut.cwt.male.coffs <- as.matrix(wavCWT.obj.male.permut)
  # permut.sex.diff.coffs <- as.matrix(wavCWT.obj.sex.diff.permut)
  
  # permuted.cwt.fem.coffs.list[[iter]] <- as.matrix(wavCWT.obj.fem.permut)
  # permuted.cwt.male.coffs.list[[iter]] <- as.matrix(wavCWT.obj.male.permut)
  # permuted.cwt.sex.diff.coffs.list[[iter]] <- as.matrix(wavCWT.obj.sex.diff.permut)

  # permuted.cwt.male.coffs.df[iter] <- cbind(permuted.cwt.male.coffs.df, as.vector(as.matrix(wavCWT.obj.male.permut)))
  # if (iter==1){plot(permut.cwt.fem.coffs[,scale], type='l'); lines(permut.cwt.male.coffs[,scale], type='l')} else {lines(permut.cwt.fem.coffs[,scale], type='l'); lines(permut.cwt.male.coffs[,scale], type='l')}
  if (iter==1){plot(permut.cwt.fem.coffs[,scale], type='l'); lines((permut.cwt.male.coffs[,scale] + (min(permut.cwt.fem.coffs[,scale]) - min(permut.cwt.male.coffs[,scale]))), type='l')} else {lines(permut.cwt.fem.coffs[,scale], type='l'); lines((permut.cwt.male.coffs[,scale] + (min(permut.cwt.fem.coffs[,scale]) - min(permut.cwt.male.coffs[,scale]))), type='l')}
  
  # if (iter==1){plot(permut.sex.diff.coffs[,scale], type='l')} else {lines(permut.sex.diff.coffs[,scale], type='l')}
  
  cat(iter, " ")
  
}

points(cwt.fem.coffs[,scale], pch=16, cex=0.5, col='red')
lines(cwt.fem.coffs[,scale], pch=16, cex=0.5, col='red')

# points(permut.sex.diff.coffs[,scale], pch=16, cex=0.5, col='blue')
# lines(permut.sex.diff.coffs[,scale], pch=16, cex=0.5, col='blue')

# points((cwt.male.coffs[,scale] + (min(cwt.fem.coffs[,scale]) - min(cwt.male.coffs[,scale]))), pch=16, cex=0.5, col='blue')
# lines((cwt.male.coffs[,scale] + (min(cwt.fem.coffs[,scale]) - min(cwt.male.coffs[,scale]))), pch=16, cex=0.5, col='blue')


points((cwt.male.coffs[,scale] + (min(cwt.fem.coffs[,scale]) - min(cwt.male.coffs[,scale]))), pch=16, cex=0.5, col='blue')
lines((cwt.male.coffs[,scale] + (min(cwt.fem.coffs[,scale]) - min(cwt.male.coffs[,scale]))), pch=16, cex=0.5, col='blue')

# dim(permuted.fem.coffs.df)


# sliding window for CO counts female vs. male
CreateBinSpanCoordinates <- function(bin.size, chrom.length, shift.size) {
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



bin.span.coord.df <- CreateBinSpanCoordinates(33, length(cos.vec.fem), 3)
colnames(bin.span.coord.df) <- c("start.sites", "end.sites")

binned.fem <- NULL
binned.male <- NULL
for (bin in 1:nrow(bin.span.coord.df)) { # bin <- 1
  
  start.loc <- bin.span.coord.df[bin, "start.sites"]
  end.loc <-  bin.span.coord.df[bin, "end.sites"]
  
  fem.bin.sum <- sum(cos.vec.fem[start.loc:end.loc])
  binned.fem <- c(binned.fem, fem.bin.sum)
  
  male.bin.sum <- sum(cos.vec.male[start.loc:end.loc])
  binned.male <- c(binned.male, male.bin.sum)
  
}
plot(binned.male, type="l", col='blue')
lines(binned.fem, type="l", col='red')




# obtain a matrix with variance
var.mat <- apply(cwt.fem.coffs, 2, function(x){y <- var(x); z <- rep(y, length(x)); return(z)})
var.mat.2 <- apply(cwt.fem.coffs.2, 2, function(x){y <- var(x); z <- rep(y, length(x)); return(z)})

image(1:nrow(var.mat), log(1:ncol(var.mat)), log(var.mat))












master.mat <- matrix(rep(NA, 1685*474), 1685, 474, byrow=TRUE)
for (i in 1:nrow(master.mat)) { #i <- 1
  for (j in 1:ncol(master.mat)) { #j <- 2
    coff.to.eval <- abs(cwt.sex.diff.coffs[i, j])
    
    temp.coff.vec <- NULL
    for (iter in 1:100) { # iter <- 1
      temp.mat <- permuted.cwt.sex.diff.coffs.list[[iter]]
      temp.coff <- abs(temp.mat[i, j])
      temp.coff.vec <- c(temp.coff.vec, temp.coff)
    }
    upper.tail.p <- 1-(sum(temp.coff.vec <= coff.to.eval) / length(temp.coff.vec))
    master.mat[i, j] <- upper.tail.p
  }
  cat(i, " ")
}

image(1:nrow(cwt.fem.coffs), log(1:ncol(cwt.fem.coffs), base=2), cwt.fem.coffs, col=blue)
write.csv(master.mat, file=paste0(statistics.output.path, "master_mat.csv"), quote=FALSE, row.names=FALSE, col.names=FALSE )
image(cwt.sex.diff.coffs)

########################################################################################################################################################################################################
# This is the second path to the same analysis as immediately above by obtaining the difference between males and females after the CWT transformation is carried out.
post.cwt.sex.diff <- (cwt.fem.coffs-cwt.male.coffs)
image(1:nrow(post.cwt.sex.diff), log(1:ncol(post.cwt.sex.diff), base=2), post.cwt.sex.diff)

post.cwt.sex.diff.permuted <- matrix(rep(NA, 1685*474), 1685, 474, byrow=TRUE)

for (i in 1:nrow(post.cwt.sex.diff.permuted)) { #i <- 1
  for (j in 1:ncol(post.cwt.sex.diff.permuted)) { #j <- 2
    
    coff.to.eval <- abs(post.cwt.sex.diff[i, j])
    
    temp.coff.vec <- NULL
    for (iter in 1:100) { # iter <- 1
      temp.mat <- (permuted.cwt.fem.coffs.list[[iter]] - permuted.cwt.male.coffs.list[[iter]])
      temp.coff <- abs(temp.mat[i, j])
      temp.coff.vec <- c(temp.coff.vec, temp.coff)
    }
    
    upper.tail.p <- 1-(sum(temp.coff.vec <= coff.to.eval) / length(temp.coff.vec))
    post.cwt.sex.diff.permuted[i, j] <- upper.tail.p
  }
  cat(i, " ")
}


########################################################################################################################################################################################################
# As David advised make histograms of the SEX-permuted and non-permuted wavelet coefficients

par(mfrow = c(1,2))
mx <- max(as.vector(wavCWT.obj.fem))
mn <- min(as.vector(wavCWT.obj.fem))
hist(as.vector(wavCWT.obj.fem), col=rgb(173, 216, 230, max=255, alpha=100, names="lt.blue"), xlim = c(mn, mx), breaks = seq(mn-1, mx+1, 0.5), xlab="wavelet coefficients:female", main=NULL)
mx <- max(as.vector(wavCWT.obj.male))
mn <- min(as.vector(wavCWT.obj.male))
hist(as.vector(wavCWT.obj.male), col=rgb(255, 192, 203, max=255, alpha=100, names ="lt.pink"), xlim = c(-15, 15), breaks = seq(mn-1, mx+1, 0.5), xlab="wavelet coefficients:male", main=NULL)
dev.off()


observed.diff.sexes <- as.vector(wavCWT.obj.fem)-as.vector(wavCWT.obj.male)
permuted.diff.sexes <- as.vector(apply((permuted.fem.coffs.df - permuted.male.coffs.df), 1, mean))

par(mfrow = c(1,1))
hist(permuted.diff.sexes, col=rgb(255, 192, 203, max=255, alpha=80, names ="lt.pink"), xlim=c(-15, 15), ylim=c(0,500000), breaks = seq(-15,15,0.5), main=NULL)
hist(observed.diff.sexes, col=rgb(173, 216, 230, max=255, alpha=80, names="lt.blue"), xlim=c(-15, 15), breaks = seq(-15,15,0.5), xlab="wavelet coefficients", main="Diference_between_sexes", add=TRUE)

par(mfrow = c(1,1))
dev.off()

# Caluculate the quantiles
diff.vec.with.all.iterations <- as.vector(permuted.fem.coffs.df - permuted.male.coffs.df)
quantile(diff.vec.with.all.iterations, c(0.05, 0.1, 0.5, 0.95, 0.99))


plotting results from the CWT
blue <- c(seq(0, 0.9, length=100))
red <- c(seq(0, 0.9, length=100))
col <- c(rgb(blue, blue, 1), rgb(1, 1, 1), rev(rgb(1, red, red)))

# define a normalizing function
# norm <- function(data) {
#   norm <- (data-mean(data))/sd(data)
#   return(norm)
# }

# m<-matrix(1:4, 4, 1)
# layout(m, height=c(2,3,2,3))
# par(xaxs="i")
# par(mar=c(0, 2, 0, 2))
# 
# plot(cos.vec.fem, type="l", xaxt="n")
# image(as.matrix(wavCWT.obj.fem)[nrow(wavCWT.obj.fem):1,], col=col, axes=F)
# 
plot(cos.vec.male, type="l", xaxt="n")
image(CWT.coffs.male, col=col, axes=F)
# 
# par(mfrow=c(1,1))


# plot(wavCWT.obj.fem, series=TRUE, col=col)

########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################





















########################################################################################################################################################################################################
# This section is still work in progress
cor.test(fem.coffs[15166:16850], male.coffs[15166:16850])
cor.test(fem.coffs[16851:length(fem.coffs)], male.coffs[16851:length(fem.coffs)])
cor.test(fem.coffs[13480:15168], male.coffs[13480:15168])




par(mfrow=c(2,2))

plot(cos.vec.male, type = "l")
plot(cos.vec.fem, type = "l")
plot(haarfisz.cos.vec.male, type = "l")
plot(haarfisz.cos.vec.male, type = "l")


# CW transformation of the signals
output1 <- Mod(DOG(haarfisz.cos.vec.male, 9, nvoice=10, twoD=TRUE, plot=F, moment=2));
output2 <- Mod(DOG(haarfisz.cos.vec.fem, 9, nvoice=10, twoD=TRUE, plot=F, moment=2));

# define a normalizing function
norm <- function(data) {
  norm <- (data-mean(data))/sd(data)
  return(norm)
}

output1 <- apply(output1,2,norm)
output2 <- apply(output2,2,norm)

# plotting the wavelet transform
m<-matrix(1:4,4,1)
layout(m, height=c(2,3,2,3))
par(xaxs="i")
par(mar=c(0,2,0,2))

plot(haarfisz.cos.vec.male,type="l", xaxt="n");
image(output1,col=col,axes=F);
plot(haarfisz.cos.vec.fem,type="l", xaxt="n");
image(output2,col=col,axes=F);
dev.off()











for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr01"
  

  
  # This is a discrete wavelet transform of the data
 
  wavDWPT.obj.fem <- wavMODWT(cos.vec.fem, wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))))
  wavVar.fem <-  wavVar(wavDWPT.obj.fem, xform="modwt", wavelet="haar", n.levels=NULL)
  
  wavVar.fem<- wavVar(wavDWPT.obj.fem, xform="modwt", wavelet="haar", n.levels=NULL,
         position=list(from=1,by=1,units=character()), units=character(),
         documentation=character(), sdf=NULL, sdfargs=NULL,
         sampling.interval=1, n.fft=1024)
  
  
   wavDWPT.obj.fem <- wavDWPT(cos.vec.fem, wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))), 
                             position = list(from = 1, by = 1, units = character()), units=character(), title.data=character(), documentation=character())
 
  
  
  
   plot(wavDWPT.obj.fem)
  summary(wavDWPT.obj.fem)
  
  wavDWPT.obj.fem <- wavDWT(x.fem, wavelet="haar", n.levels=as.integer(floor(logb(length(x.fem),base = 2))), keep.series = TRUE)
  plot(wavShift(wavDWPT.obj.fem))
  eda.plot(wavDWPT.obj.fem)
  summary(wavDWPT.obj.fem)
  
  
  
  
  
  
  x.fem <- (cos.vec.fem - mean(cos.vec.fem)) / sqrt(var(cos.vec.fem)) # data needs to be normalized before applying the analysis 
  x.male <- (cos.vec.male - mean(cos.vec.male)) / sqrt(var(cos.vec.male)) # data needs to be normalized before applying the analysis 
  x.diff <- ((x.fem - x.male) - mean((x.fem - x.male))) / sqrt(var((x.fem - x.male)))
  
  # commparing to a white-noise model
  # this is a contrasting example 
  # D.statistic(rnorm(1024))
  # D.statistic(cumsum(rnorm(1024)))
  D.statistic(x.male)
  D.statistic(x.fem)
  
  D.table(n.sample = c(127, 130), significance = c(0.1, 0.05, 0.01), lookup = TRUE, n.realization = 10000, n.repetition = 3, tolerance = 1e-6)
  D.table(n.sample = c(length(x.fem)), significance = c(0.1, 0.05, 0.01), lookup = TRUE, n.realization = 10000, n.repetition = 3, tolerance = 1e-6)
  
  
  
  wavCWT.obj.fem <- wavCWT(x.fem, scale.range = deltat(x.fem) * c(1, length(x.fem)), n.scale = 100, wavelet = "gaussian2", shift = 5, variance = 1)
  wavCWT.obj.male <- wavCWT(x.male, scale.range = deltat(x.male) * c(1, length(x.male)), n.scale = 100, wavelet = "gaussian2", shift = 5, variance = 1)
  wavCWT.obj.diff <- wavCWT(x.diff, scale.range = deltat(x.diff) * c(1, length(x.diff)), n.scale = 100, wavelet = "gaussian2", shift = 5, variance = 1)
  
  
  fem.mat <- as.matrix(wavCWT.obj.fem)
  male.mat <- as.matrix(wavCWT.obj.male)
  
  corr.coffs.vec <- NULL
  for (i in 1 : ncol(fem.mat)) {
    corr.coffs <- cor(fem.mat[,i],male.mat[,i] )
    corr.coffs.vec <- c(corr.coffs.vec, corr.coffs)
  }
  plot(corr.coffs.vec, type = "l")
  
  
  plot(fem.mat[,78], male.mat[,78])
  
  
  
  
  
  
  
  
  jpeg(file = paste0(statistics.output.path, "sex_diff_wav_plot_", chrom.name, ".jpg"), width = 16, height = 9, unit = "in", res = 600)
  plot(wavCWT.obj.diff, series = T, main = paste0(chrom.name, "_Sex_difference"))
  dev.off()
  
  # plot(wavCWT.obj.diff, series = T, main = paste0(chrom.name, "_Sires"))
  
  jpeg(file = paste0(statistics.output.path, "wav_plot_", chrom.name, ".jpg"))
  for (i in 1:2) {
    if (i == 1) {
      ifultools::splitplot(1,2,1)
      plot(wavCWT.obj.fem, series = T, main = paste0(chrom.name, "_Dams"), add = TRUE)
    } else {
      ifultools::splitplot(1,2,2)
      plot(wavCWT.obj.male, series = T, main = paste0(chrom.name, "_Sires"), add = TRUE)
    }
  }
  dev.off()
  
  
}




par(mfrow = c(2,2))
hist(as.vector(wavCWT.obj.fem), col = rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue"), xlim = c(-8, 8), breaks = seq(-8,8,0.5), xlab = "wavelet coefficients:female", main = NULL)
hist(as.vector(wavCWT.obj.male), col = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink"), xlim = c(-8, 8), breaks = seq(-8,8,0.5), xlab = "wavelet coefficients:male", main = NULL)
hist(as.vector(wavCWT.obj.diff), col = 'blue', xlim = c(-8, 8), breaks = seq(-8,8,0.5), xlab = "wavelet coefficient difference between sexes", main = NULL)

wav.coffs.fem <- as.matrix(wavCWT.obj.fem)
wav.coffs.male <- as.matrix(wavCWT.obj.male)
diff.wav.coffs <- wav.coffs.fem - wav.coffs.male

########################################################################################################################################################################################################