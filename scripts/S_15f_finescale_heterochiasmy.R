#!/usr/bin/R
# setwd("/group/difazio/populus/gatk-7x7-Stet14/")
setwd("/Users/cabeyrat/Google Drive File Stream/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS, for the recombination landscape analysis for 7x7 cross
# DATE: December 23, 2020
# USAGE: Input files include average CO counts by parent and chromosome as well as a dataframe for chromosome physical size and genetc map size. 

# Following is the input file format used 
# CHROM  PAR   IND CO_SIZE DWNSTRM_FLNK UPSTRM_FLNK
# Chr01 1863 24708  103140      4675047     4778187
# Chr01 1863 24708   30914     26411972    26442886
# Chr01 1863 24708   61845     34445658    34507503
# Chr01 1863 24708    5612     46505278    46510890
# Chr01 1863 24709  164382       726848      891230
# Chr01 1863 24709  153359      2513166     2666525

# PURPOSE: Investigate SEX based CO count differences in the script 15b.. as it relates to windowed analysis (win.size=1.2Mb)
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
library(lme4)
require(car)
library(MASS)
# library(lmtest)
# library(glmmTMB)

########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
statistics.output.path <- "./statistics/gross_stat_plots/"
sliding.avg.output.path <- "./statistics/heterochiasmy_analysis/sliding_average/"

########################################################################################################################################################################################################
# p.val.dist <- NULL
for (chrom.name in chrom.name.vec) {  # chrom.name <- chrom.name.vec[2]
  jpeg(filename = paste0(sliding.avg.output.path, chrom.name, "_sliding_window_p_values.jpg"), width = 2000, height = 1000, quality = 150)
  op <- par(mfrow=c(3,3))
  for (scale in 1:8) {
    
    # read-in the data fils
    binned.cos.by.par.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs_by_par.csv"), stringsAsFactors=FALSE, header=TRUE)
    num.ind.df <- read.csv(paste0(statistics.input.path, "number_of_individuals_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
    num.COs.df <- read.csv(paste0(statistics.input.path, "number_of_crossovers_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
    colnames(binned.cos.by.par.df) <- colnames(binned.cos.by.par.df) %>% str_replace("^X", "")
    colnames(num.ind.df) <- colnames(num.ind.df) %>% str_replace("^X", "")
    colnames(num.COs.df) <- colnames(num.COs.df) %>% str_replace("^X", "")
    
    
    # produce dataframe for linear model analysis
    num.overlap.win <- nrow(binned.cos.by.par.df)
    reshaped.co.count.df <- NULL
    for (par.name in par.name.vec) {
      sex <- ifelse(par.name %in% par.name.vec[1:7], "F", "M")
      reshaped.co.count.df <- rbind(reshaped.co.count.df, cbind(c(1:num.overlap.win),
                                                                rep(par.name, num.overlap.win),
                                                                rep(num.ind.df[which(chrom.name.vec==chrom.name),par.name], num.overlap.win),
                                                                rep(sex, num.overlap.win),
                                                                binned.cos.by.par.df[,par.name]))
    }
    reshaped.co.count.df <- as.data.frame(reshaped.co.count.df, stringsAsFactors=FALSE)
    colnames(reshaped.co.count.df) <- c("OV.WIN.IDX", "PAR", "FAM_SIZE", "SEX", "CO_COUNT")
    reshaped.co.count.df$OV.WIN.IDX <- as.factor(reshaped.co.count.df$OV.WIN.IDX)
    reshaped.co.count.df$PAR <- as.factor(reshaped.co.count.df$PAR)
    reshaped.co.count.df$SEX <- as.factor(reshaped.co.count.df$SEX)
    reshaped.co.count.df$FAM_SIZE <- as.numeric(reshaped.co.count.df$FAM_SIZE)
    reshaped.co.count.df$CO_COUNT <- as.numeric(reshaped.co.count.df$CO_COUNT)
    
    # library(fitdistrplus)
    # fnb <- fitdist(reshaped.co.count.df$CO_COUNT, distr = "nbinom", method = 'mle')
    # fp <- fitdist(reshaped.co.count.df$CO_COUNT, distr = "pois", method = 'mle')
    # fg <- fitdist(reshaped.co.count.df$CO_COUNT, distr = "geom", method = 'mle')
    # 
    # op <- par(mfrow = c(2, 2))
    # plot.legend <- c("nbinom", "pois", "geom")
    # denscomp(list(fnb, fp, fg), legendtext = plot.legend)
    # qqcomp(list(fnb, fp, fg), legendtext = plot.legend)
    # plot.legend <- c("nbinom", "pois", "geom")
    # cdfcomp(list(fnb, fp, fg), legendtext = plot.legend)
    # ppcomp(list(fnb, fp, fg), legendtext = plot.legend)
    # par(op)
    # 
    # gof <- gofstat(list(fnb, fp, fg), fitnames = c("nbinom", "pois", "geom"))
    # 
    # glm.nb <- MASS::glm.nb(CO_COUNT ~ SEX + OV.WIN.IDX + FAM_SIZE, link = "log", data = reshaped.co.count.df)
    # sum.glm.nb <- summary(glm.nb, cor = FALSE)
    # win.coff.names <- names(glm.nb$coefficients)[-c(1:2, length(glm.nb$coefficients))]
    # win.idx.order <- order(as.numeric(sapply(win.coff.names, function(x){str_replace(x, pattern='OV.WIN.IDX', replacement="")})), decreasing=FALSE)
    # win.idx.coffs <- sum.glm.nb$coefficients[-c(1:2, length(glm.nb$coefficients)),]
    # win.idx.coffs <- win.idx.coffs[win.idx.order, ]
    # plot(-log10(win.idx.coffs[, 4]))
    # op <- par(mfrow=c(2,2))
    # plot(glm.nb)
    # par(op)
    
    
    # Chi-square test for FEM vs. MALE groups with permutation derived p-values
    # binned.permut.cos.fem.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_FEM.csv"), stringsAsFactors=FALSE, header=TRUE)
    # binned.permut.cos.male.df <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_permuted_binned_COs_MALE.csv"), stringsAsFactors=FALSE, header=TRUE)
    # actual.binned.COs <- read.csv(file=paste0(sliding.avg.output.path, chrom.name, "_", scale, "_binned_COs.csv"), stringsAsFactors=FALSE, header=TRUE)
    # 
    # actual.ch.sq.stat <- NULL
    # for (i in 1:nrow(actual.binned.COs)) {
    #   if (all(c(actual.binned.COs[i,1], actual.binned.COs[i,2]) == 0)) {
    #     null.ch.sq.stat <- c(null.ch.sq.stat, NA)
    #   } else {
    #     chi.mod <- chisq.test(c(actual.binned.COs[i,1], actual.binned.COs[i,2]))
    #     actual.ch.sq.stat <- c(actual.ch.sq.stat, chi.mod$statistic)
    #   }
    # }
    # 
    # null.ch.sq.stat <- NULL
    # for (iter in 1:ncol(binned.permut.cos.fem.df)) {
    #   tmp.ch.sq.stat <- NULL
    #   for (i in 1:nrow(binned.permut.cos.fem.df)) {
    #     if (all(c(binned.permut.cos.fem.df[i,iter], binned.permut.cos.male.df[i,iter]) == 0)) {
    #       tmp.ch.sq.stat <- c(tmp.ch.sq.stat, NA)
    #     } else {
    #       tmp.chi.mod <- chisq.test(c(binned.permut.cos.fem.df[i,iter], binned.permut.cos.male.df[i,iter]))
    #       tmp.ch.sq.stat <- c(tmp.ch.sq.stat, tmp.chi.mod$statistic)
    #     }
    #   }
    #   null.ch.sq.stat <- cbind(null.ch.sq.stat, tmp.ch.sq.stat)
    #   cat(iter, " ")
    # }
    # 
    # 
    # empirical.p.val <- NULL
    # for (i in 1:nrow(actual.binned.COs)) {
    #   if (is.na(actual.ch.sq.stat[i])) {
    #     empirical.p.val <- c(empirical.p.val, NA)
    #   } else {
    #     empirical.p.val <- c(empirical.p.val, sum(actual.ch.sq.stat[i] > null.ch.sq.stat[i,])/ncol(null.ch.sq.stat))
    #   }
    #   
    # }
    # sum(is.na(empirical.p.val))
    # plot(-log10(empirical.p.val), type='l')
    
    
    
    # Running a Poisson GLM for each window
    effect.size.vec <- NULL
    p.value.vec <- NULL
    for (i in 1:nrow(binned.cos.by.par.df)) {
      subset.reshaped.co.count.df <- reshaped.co.count.df[reshaped.co.count.df$OV.WIN.IDX == i, ]
      # tmp.glm.nb <- glm.nb <- MASS::glm.nb(CO_COUNT ~ SEX + FAM_SIZE, link = "log", data = subset.reshaped.co.count.df)
      # sum.tmp.glm.nb <- summary(tmp.glm.nb)
      tmp.glm.poiss <- stats::glm(CO_COUNT ~ SEX + FAM_SIZE, data = subset.reshaped.co.count.df, family = poisson(link = "log")) 
      sum.tmp.glm.poiss <- summary(tmp.glm.poiss)
      # AIC(tmp.glm.poiss, tmp.glm.nb)
      effet.size <- sum.tmp.glm.poiss$coefficients["SEXM", 1]
      p.value <- sum.tmp.glm.poiss$coefficients["SEXM", 4]
      effect.size.vec <- c(effect.size.vec, effet.size)
      p.value.vec <- c(p.value.vec, p.value)
      
    }
    # plot(effect.size.vec, type='o', pch=16, cex=0.5)
    plot(-log10(p.value.vec), pch=16, cex=0.5, col='black',type='o', main=paste0(chrom.name, " d", scale), ylim=c(0, 8))
    # abline(h=4.5)
    # hist(p.value.vec, breaks=100)
    # sum(p.value.vec < 0.05/1430)
    cat(chrom.name, " ", scale, "\n")
  }
  par(op)
  dev.off()
}

# hist(p.val.dist, breaks=100)






