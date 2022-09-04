#!/usr/bin/R
########################################################################################################################################################################################################

# PURPOSE: Investigate SEX based CO count differences in the script 15b.. as it relates to windowed analysis (win.size=960kb)
# METADATA: Written by CRA & DMS, for the recombination landscape analysis for 7x7 cross
# DATE: December 23, 2020; modified Fabruary 18, 2021
# USAGE: Input files include average CO counts by parent and chromosome as well as a dataframe for chromosome physical size and genetc map size. 

# head(binned.cos.by.par.df)
# CHROM   START     END 1863 1909 1950 2048 2066 2283 4593 2365 2393 2515 2572 2683 6909 7073 WIN.IDX
# Chr01       1  359287    0    0    0    3    0    0    1    0    0    0    0    0    0    0       1
# Chr01  359288 1319287    3    6   12   10    4    5    2    3    8    5    5    6    3    6       2
# Chr01 1319288 2279287    5    5    4   11    9    5    5   13   10    6    3    7    8   11       3
# Chr01 2279288 3239287    7    3    6    5    2    2    8    9   10    8   10    4    4    9       4
# Chr01 3239288 4199287   10    8   11    7    7   13    4    4    8   17    4    4    4   10       5
# Chr01 4199288 5159287    9    9   11   13   11   13    6   16   14   10   11    9   15   11       6
 

# head(gen.feat.df)
# CHROM   START     END   pct_at   pct_gc  num_A  num_C  num_G  num_T num_N num_oth seq_len COPIA GYPSY SIMP_REPEAT GENE CUM_CO_COUNT log_CUM_CO_COUNT WIN.IDX
# Chr01       1  359287 0.638862 0.361138 115219  64623  65129 114315     0       0  359286    40    60         116   36            4         0.607455       1
# Chr01  359288 1319287 0.660368 0.339632 320504 163073 162973 313449     0       0  959999    28   101         314  117           78         1.892373       2
# Chr01 1319288 2279287 0.658806 0.341194 317480 165571 161975 314973     0       0  959999    36    56         361  120          102         2.008813       3
# Chr01 2279288 3239287 0.649534 0.350466 310243 168812 167635 313309     0       0  959999    65    92         275  102           87         1.939769       4
# Chr01 3239288 4199287 0.655112 0.344888 317039 165199 165893 311868     0       0  959999    52    96         304   93          111         2.045519       5
# Chr01 4199288 5159287 0.667103 0.332897 321875 160434 159147 318543     0       0  959999    32    37         380  108          158         2.198795       6

########################################################################################################################################################################################################

# set working directory and load the packages used in the script
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
rm(list = ls())
library(dplyr)
library(stringr)
library(lme4)
library(MASS)
library(qvalue)
library(exactci)
# library(lmtest)
# library(glmmTMB)
# 
########################################################################################################################################################################################################

# set main input and output paths
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
statistics.output.path <- "./statistics/gross_stat_plots/"
sliding.avg.output.path <- "./statistics/heterochiasmy_analysis/sliding_average/"
gen.feat.output.path <- "./genomic_features/analysis_output/"


########################################################################################################################################################################################################

# read-in the data files
# binned.cos.by.par.df <- read.csv(file = paste0(statistics.input.path, "co.bin_size_1.92Mb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
# binned.cos.by.par.df <- read.csv(file = paste0(statistics.input.path, "co.bin_size_3.84Mb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
binned.cos.by.par.df <- read.csv(file = paste0(statistics.input.path, "co.bin_size_960kb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
colnames(binned.cos.by.par.df) <- colnames(binned.cos.by.par.df) %>% str_replace("^X", "")
binned.cos.by.par.df <- cbind(binned.cos.by.par.df, WIN.IDX=seq(1, nrow(binned.cos.by.par.df), by=1))




# removing the first and the last window for each LG as these are smaller than the assigned window size
binned.cos.by.par.df.2 <- NULL
for (chrom.name in chrom.name.vec) { 
  tmp.subset <- binned.cos.by.par.df[binned.cos.by.par.df$CHROM == chrom.name, ]
  idx.to.rmv <- c(1, nrow(tmp.subset))
  tmp.subset <- tmp.subset[-idx.to.rmv, ]
  binned.cos.by.par.df.2 <- rbind(binned.cos.by.par.df.2, tmp.subset)
}
# binned.cos.by.par.df.2[binned.cos.by.par.df.2$CHROM == 'Chr01' & binned.cos.by.par.df.2$START > 27239288, '6909'] <- NA # assign NA to Chr01 of 6909 from 27239288 to the end due to lack of markers.
binned.cos.by.par.df.2[binned.cos.by.par.df.2$CHROM == 'Chr01' & binned.cos.by.par.df.2$START > 29054288, '6909'] <- NA # assign NA to Chr01 of 6909 from 27239288 to the end due to lack of markers.
binned.cos.by.par.df.3 <- cbind(binned.cos.by.par.df.2, 
                                CUM_CO=apply(binned.cos.by.par.df.2[,4:17], 1, sum, na.rm=TRUE), 
                                FEM_CO=apply(binned.cos.by.par.df.2[,4:10], 1, sum, na.rm=TRUE),
                                MALE_CO=apply(binned.cos.by.par.df.2[,11:17], 1, sum, na.rm=TRUE),
                                AVG_CO=apply(binned.cos.by.par.df.2[,11:17], 1, mean, na.rm=TRUE),
                                VAR_CO=apply(binned.cos.by.par.df.2[,11:17], 1, var, na.rm=TRUE)
                                )

dev.new()
op <- par(mfrow=c(4,5))
for (chrom.name in chrom.name.vec) {
  subset.temp.df <- binned.cos.by.par.df.3[binned.cos.by.par.df.3$CHROM == chrom.name, ]
  plot(subset.temp.df$WIN.IDX, subset.temp.df$AVG_CO, type='h', ylab='mean CO count', xlab='non-overlap win.idx (960kb)', main=chrom.name)
}
par(op)

# remove windows that overalp with centromeres and telomeres
binned.cos.by.par.df.4 <- binned.cos.by.par.df.3[binned.cos.by.par.df.3$AVG_CO > 4, ]
nrow(binned.cos.by.par.df.4)

# check the distribuution of CO counts
co.count.dist.gw <- as.integer(as.matrix(binned.cos.by.par.df.4[4:17]))
dev.new()
op <- par(mfrow=c(2,1))
hist(co.count.dist.gw, breaks=100, xlab='CO count in non-overlapping window (960kbps)', main='genomewide windows CO-count histogram') # The zero counts are inflated in this distribution.
hist(rpois(n=length(co.count.dist.gw), lambda=mean(co.count.dist.gw, na.rm=T)), col='blue', breaks=100, xlab='CO count in non-overlapping window (960kbps)', main='Theoretical poisson')
par(op)
qqplot(co.count.dist.gw, rpois(n=length(co.count.dist.gw), lambda=mean(co.count.dist.gw, na.rm=T)))
abline(0,1)

# read-in the total COs and family size data for parents and chromosomes.
num.ind.df <- read.csv(paste0(statistics.input.path, "number_of_individuals_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
num.COs.df <- read.csv(paste0(statistics.input.path, "number_of_crossovers_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
colnames(num.ind.df) <- colnames(num.ind.df) %>% str_replace("^X", "")
colnames(num.COs.df) <- colnames(num.COs.df) %>% str_replace("^X", "")



# Calculating the index of dispersion for each window (var/mean)
binned.cos.by.par.df.4 <- cbind(binned.cos.by.par.df.4, IOD=binned.cos.by.par.df.4$VAR_CO/binned.cos.by.par.df.4$AVG_CO) # calculating index of dispersion
head(binned.cos.by.par.df.4)

# histograms of variance, mean and IOD in windows
dev.new()
op <- par(mfrow=c(2,4))
plot(binned.cos.by.par.df.4$AVG_CO, binned.cos.by.par.df.4$VAR_CO, xlab='average CO counts per window', ylab='variance CO counts per window', main=NULL)
abline(0,1)
hist(binned.cos.by.par.df.4$AVG_CO, breaks=100, xlab='Averag CO distribution', main=NULL)
hist(binned.cos.by.par.df.4$VAR_CO, breaks=100, xlab='CO count variance within windows', main=NULL)
hist(binned.cos.by.par.df.4$IOD, breaks=100, xlab='indx of dispersion (var/mean)', main=NULL)
plot(binned.cos.by.par.df.4$WIN.IDX, binned.cos.by.par.df.4$AVG_CO, ylab='mean(CO count)', xlab='genome-wide window idx.(960kb)', main=NULL)
plot(binned.cos.by.par.df.4$WIN.IDX, binned.cos.by.par.df.4$VAR_CO, ylab='var(CO count)', xlab='genome-wide window idx.(960kb)', main=NULL)
plot(binned.cos.by.par.df.4$WIN.IDX, binned.cos.by.par.df.4$IOD, ylab='Index of Dispersion (var/mean)', xlab='genome-wide window idx.(960kb)', main=NULL)
par(op)


# check the longitudinal patterns of IOD
dev.new()
plot(binned.cos.by.par.df.4$WIN.IDX, binned.cos.by.par.df.4$IOD, main='pre-filteration based on IOD') # same as above produced plot

# histogram of IOD
dev.new()
hist(binned.cos.by.par.df.4$IOD, breaks=100)

quantile(binned.cos.by.par.df.4$IOD, c(0.1, 0.9))
quantile(binned.cos.by.par.df.4$IOD, c(0.05, 0.95))

# remove the most extreme/outliers 5% and 95% of values
iod.low.hi.cutoffs <- quantile(binned.cos.by.par.df.4$IOD, c(0.05, 0.95)) # setting low and hi index of dispersion cutoffs. For null Poisson iod=1 iod.low.hi.cutoffs <- c(0.33, 3) 
windows.to.rmv.idx <- which(binned.cos.by.par.df.4$IOD <= iod.low.hi.cutoffs[1] | binned.cos.by.par.df.4$IOD >= iod.low.hi.cutoffs[2])
length(windows.to.rmv.idx)


# remove windows based on a IOD cutoffs
binned.cos.by.par.df.5 <- binned.cos.by.par.df.4[binned.cos.by.par.df.4$IOD > iod.low.hi.cutoffs[1] &  binned.cos.by.par.df.4$IOD < iod.low.hi.cutoffs[2], ]
nrow(binned.cos.by.par.df.5)
dev.new()
plot(binned.cos.by.par.df.5$WIN.IDX, binned.cos.by.par.df.5$IOD, main='post-filtered based on IOD')


# visualize genoewide trends after filtering based on IOD
dev.new()
op <- par(mfrow=c(2,4))
plot(binned.cos.by.par.df.5$AVG_CO, binned.cos.by.par.df.5$VAR_CO, xlab='average CO counts per window', ylab='variance CO counts per window', main=NULL)
abline(0,1)
hist(binned.cos.by.par.df.5$AVG_CO, breaks=100, xlab='Averag CO distribution', main=NULL)
hist(binned.cos.by.par.df.5$VAR_CO, breaks=100, xlab='CO count variance within windows', main=NULL)
hist(binned.cos.by.par.df.5$IOD, breaks=100, xlab='indx of dispersion (var/mean)', main=NULL)
plot(binned.cos.by.par.df.5$WIN.IDX, binned.cos.by.par.df.5$AVG_CO, ylab='mean(CO count)', xlab='genome-wide window idx.(960kb)', main=NULL)
plot(binned.cos.by.par.df.5$WIN.IDX, binned.cos.by.par.df.5$VAR_CO, ylab='var(CO count)', xlab='genome-wide window idx.(960kb)', main=NULL)
plot(binned.cos.by.par.df.5$WIN.IDX, binned.cos.by.par.df.5$IOD, ylab='Index of Dispersion (var/mean)', xlab='genome-wide window idx.(960kb)', main=NULL)
par(op)



# op <- par(mfrow=c(4,5))
# for (chrom.name in chrom.name.vec) {
#   subset.temp.df <- binned.cos.by.par.df.5[binned.cos.by.par.df.5$CHROM == chrom.name, ]
#   plot(subset.temp.df$WIN.IDX, apply(subset.temp.df[,4:17], 1, mean, na.rm=TRUE), type='l', ylab='mean CO count', xlab='non-overlap win.idx (960kb)', main=chrom.name)
# }
# par(op)



# Check the correlation between males and females for CO counts in windows
dev.new()
plot(binned.cos.by.par.df.5$FEM_CO, binned.cos.by.par.df.5$MALE_CO, xlab='Female cumulative CO count', ylab='Male cumulative CO count')
abline(0,1)
cum.cor.test <- cor.test(binned.cos.by.par.df.5$FEM_CO, binned.cos.by.par.df.5$MALE_CO, method = 'spearman') # very high correlation observed between male and female cumulative counts across the genome.
text(80, 40, labels=paste0('R^2= ', round(cum.cor.test$estimate, digits = 3), ' p-value= ', round(cum.cor.test$p.value, digits=4)))

# KS-test for the males vs. females
# ks.test(x=binned.cos.by.par.df.5$FEM_CO, y=binned.cos.by.par.df.5$MALE_CO)
# ks.test(x=binned.cos.by.par.df.5$FEM_CO, y=binned.cos.by.par.df.5$MALE_CO, alternative="greater")
# centering and scaling the male and female vectors so that they can be compared without bias. 
# methods: A KS-test was carried out for the centered and scaled vectors of male and female CO counts in 960 kb windows 
# with centromeres and chromosome ends removed. Index of dispersion outliers were also removed.
ks.test(x=scale(binned.cos.by.par.df.5$FEM_CO, center=T, scale=T), y=scale(binned.cos.by.par.df.5$MALE_CO, center=T, scale=T))
ks.test(x=scale(binned.cos.by.par.df.5$FEM_CO, center=T, scale=F), y=scale(binned.cos.by.par.df.5$MALE_CO, center=T, scale=F))
# NOTE: p-value does not change whether scale T or F


# check the distribution of observed vs. null CO counts genomewide
observed.counts <- table(as.integer(as.matrix(binned.cos.by.par.df.5[4:17]), na.rm=TRUE))
random.counts <- table(rpois(lambda=mean(as.integer(as.matrix(binned.cos.by.par.df.5[4:17])), na.rm=TRUE), length(as.integer(as.matrix(binned.cos.by.par.df.5[4:17])))))
op <- par(mfrow=c(2,1))
plot(observed.counts)
plot(random.counts, xlim=c(0,30))
par(op)
qqplot(observed.counts, as.vector(random.counts))
abline(0,1)


# dataframe for lm() analysis.
num.overlap.win <- nrow(binned.cos.by.par.df.5)
reshaped.co.count.df <- NULL
for (par.name in par.name.vec) { # par.name <- par.name.vec[1]
  sex <- ifelse(par.name %in% par.name.vec[1:7], "F", "M")
  reshaped.co.count.df <- rbind(reshaped.co.count.df,
                                cbind(binned.cos.by.par.df.5$WIN.IDX,
                                      binned.cos.by.par.df.5[,c("CHROM", "START", "END")],
                                      rep(par.name, num.overlap.win),
                                      rep(num.ind.df[which(chrom.name.vec==chrom.name),par.name], num.overlap.win),
                                      rep(sex, num.overlap.win),
                                      binned.cos.by.par.df.5[, par.name]))
}
reshaped.co.count.df <- as.data.frame(reshaped.co.count.df, stringsAsFactors=FALSE)
colnames(reshaped.co.count.df) <- c("WIN.IDX", "CHROM", "START", "END", "PAR", "FAM_SIZE", "SEX", "CO_COUNT")
reshaped.co.count.df$WIN.IDX <- as.factor(reshaped.co.count.df$WIN.IDX)
reshaped.co.count.df$PAR <- as.factor(reshaped.co.count.df$PAR)
reshaped.co.count.df$SEX <- as.factor(reshaped.co.count.df$SEX)
reshaped.co.count.df$FAM_SIZE <- as.numeric(reshaped.co.count.df$FAM_SIZE)
reshaped.co.count.df$CO_COUNT <- as.numeric(reshaped.co.count.df$CO_COUNT)
# reshaped.co.count.df$CO_COUNT <- sample(x=reshaped.co.count.df$CO_COUNT, size=length(reshaped.co.count.df$CO_COUNT), replace=FALSE)
# reshaped.co.count.df <- reshaped.co.count.df[reshaped.co.count.df$CO_COUNT != 0, ]
# reshaped.co.count.df <- reshaped.co.count.df[!reshaped.co.count.df$CHROM %in% c(chrom.name.vec[c(1,4,12)]), ]
head(reshaped.co.count.df)




# Poisson.test of homogeneity for each window to compare cumulative CO counts between females and males. This will help track down chromosome scale CO rate differences to a finer 960kbps scale.
# Running permutations for each window to obtain the NULl distribution of p-values from Poisson.test() to check whether the assumptions are correct.
permut.p.value.vec <- NULL
mean.vec <- NULL
var.vec <- NULL
for (i in binned.cos.by.par.df.5$WIN.IDX) {
  subset.reshaped.co.count.df <- reshaped.co.count.df[reshaped.co.count.df$WIN.IDX == i, ]
  
  # remove parents with NA counts from the windowed analysis
  if (any(is.na(subset.reshaped.co.count.df$CO_COUNT))) {
    subset.reshaped.co.count.df <- subset.reshaped.co.count.df[-(which(is.na(subset.reshaped.co.count.df$CO_COUNT))), ]
  }
  
  iter.p.val.vec <- NULL
  for (iter in 1:100) {
    
    iter.subset.reshaped.co.count.df <- subset.reshaped.co.count.df
    iter.subset.reshaped.co.count.df$SEX <- sample(x=iter.subset.reshaped.co.count.df$SEX, size=length(iter.subset.reshaped.co.count.df$SEX), replace=FALSE)
    
    fem.vec <- iter.subset.reshaped.co.count.df$CO_COUNT[iter.subset.reshaped.co.count.df$SEX == 'F']
    male.vec <- iter.subset.reshaped.co.count.df$CO_COUNT[iter.subset.reshaped.co.count.df$SEX == 'M']
    # Downsampling if counts are different
    fem.count <-  ifelse(length(fem.vec) > length(male.vec), sum(sample(x=fem.vec, size=length(male.vec), replace=FALSE)), sum(fem.vec))
    male.count <- sum(male.vec)
    
    e.pois.mod <- poisson.exact(c(fem.count, male.count), T = 1, r = 1, alternative = "two.sided", tsmethod='central', plot=FALSE, midp=TRUE)
    # e.pois.mod <- poisson.test(c(fem.count, male.count), T = 1, r = 1, alternative = "two.sided")
    mid.p <- e.pois.mod$p.value
    iter.p.val.vec <- c(iter.p.val.vec, mid.p)

  }
  
  var.vec <- c(var.vec,  var(iter.p.val.vec))
  mean.vec <- c(mean.vec, mean(iter.p.val.vec))
  permut.p.value.vec <- c(permut.p.value.vec, iter.p.val.vec)
  
  cat(i, '/', binned.cos.by.par.df.5$WIN.IDX[length(binned.cos.by.par.df.5$WIN.IDX)], "\n")

}

# check for any trends visually after SEX permutation to identify as outliers
op <- par(mfrow=c(2,2))
hist(permut.p.value.vec, breaks=1000)
# length(mean.vec)
# length(var.vec)
# nrow(binned.cos.by.par.df.5)
plot(mean.vec, var.vec)
plot(binned.cos.by.par.df.5$WIN.IDX, mean.vec)
plot(binned.cos.by.par.df.5$WIN.IDX, var.vec)
par(op)
# none of the windows can be identified as outliers.



# Following is a code-block to check and confirm that the p-value distribution is the same as observed above when random Poisson vectors 
# are used with lambda deisgnated uniquely for each window
null.count <- binned.cos.by.par.df.5$FEM_CO
null.p.value.vec <- NULL

iter.p.val.vec <- NULL
for (l in null.count) {
  for (i in 1:100) {
    fem.count <- rpois(1, lambda = null.count[l])
    male.count <- rpois(1, lambda = null.count[l])
    
    e.pois.mod <- poisson.exact(c(fem.count, male.count), T = 1, r = 1, alternative = "two.sided", tsmethod='central', plot=FALSE, midp=TRUE)
    # e.pois.mod <- poisson.test(c(fem.count, male.count), T = 1, r = 1, alternative = "two.sided")
    mid.p <- e.pois.mod$p.value
    
    iter.p.val.vec <- c(iter.p.val.vec, mid.p)
    
  }
  
  null.p.value.vec <- c(null.p.value.vec, iter.p.val.vec)
}
hist(null.p.value.vec, breaks=1000)



# comparing the two NULL distributions.
op <- par(mfrow=c(2,1))
hist(null.p.value.vec, breaks=1000)
hist(permut.p.value.vec, breaks=1000)
par(op)





# # Since the p-value distribution of Poisson test is a bit unwieldy carrying out a lm() analysis to identify which windows to exclude from analysis
# permut.effect.size.vec <- NULL
# permut.p.value.vec <- NULL
# permut.r.sq.vec <- NULL
# mean.vec <- NULL
# var.vec <- NULL
# for (i in unique(reshaped.co.count.df$WIN.IDX)) {
#   subset.reshaped.co.count.df <- reshaped.co.count.df[reshaped.co.count.df$WIN.IDX == i, ]
#   
#   # remove individuals with NA from the dataset entirely
#   if (any(is.na(subset.reshaped.co.count.df$CO_COUNT))) {
#     subset.reshaped.co.count.df <- subset.reshaped.co.count.df[-(which(is.na(subset.reshaped.co.count.df$CO_COUNT))), ]
#   }
#   
#   if (sum(subset.reshaped.co.count.df$CO_COUNT) == 0) {
#     mean.vec <- c(mean.vec, NA)
#     var.vec <- c(var.vec, NA)
#     next
#     # stop("all counts are Zero for the window ! \n")
#     
#   } else {
#     
#     iter.p.val.vec <- NULL
#     iter.effect.vec <- NULL
#     iter.r.square.vec <- NULL
#     
#     for (iter in 1:100) {
#       iter.subset.reshaped.co.count.df <- subset.reshaped.co.count.df
#       iter.subset.reshaped.co.count.df$SEX <- sample(iter.subset.reshaped.co.count.df$SEX, length(iter.subset.reshaped.co.count.df$SEX), replace=FALSE)
#       # iter.tmp.glm.poiss <- stats::glm(CO_COUNT ~ SEX + scale(FAM_SIZE), data = iter.subset.reshaped.co.count.df, family = poisson(link = "log"))
#       iter.tmp.glm.poiss <- lm(log(CO_COUNT + 1) ~ SEX + scale(FAM_SIZE), data = iter.subset.reshaped.co.count.df)
#       iter.sum.tmp.glm.poiss <- summary(iter.tmp.glm.poiss)
#       
#       
#       iter.r.sq <- iter.sum.tmp.glm.poiss$r.squared
#       iter.effet.size <- iter.sum.tmp.glm.poiss$coefficients["SEXM", 1]
#       iter.p.value <- iter.sum.tmp.glm.poiss$coefficients["SEXM", 4]
#       
#       iter.effect.vec <- c(iter.effect.vec, iter.effet.size)
#       iter.p.val.vec <- c(iter.p.val.vec, iter.p.value)
#       iter.r.square.vec <- c(iter.r.square.vec, iter.r.sq)
#     }
#     
#     var.vec <- c(var.vec,  var(iter.p.val.vec))
#     mean.vec <- c(mean.vec, mean(iter.p.val.vec))
#     
#     permut.effect.size.vec <- c(permut.effect.size.vec, iter.effect.vec)
#     permut.p.value.vec <- c(permut.p.value.vec, iter.p.val.vec)
#     permut.r.sq.vec <- c(permut.r.sq.vec, iter.r.square.vec)
#   }
# }
# # hist(permut.p.value.vec, breaks=100)
# length(mean.vec)
# length(var.vec)
# nrow(binned.cos.by.par.df.3)
# 
# plot(mean.vec, var.vec)
# # temp.idx <- which((var.vec > 0.1 | var.vec < 0.06) | (mean.vec > 0.55 | mean.vec < 0.45))
# temp.idx <- which(var.vec < 0.06 | var.vec > 0.1)
# length(temp.idx)
# binned.cos.by.par.df.3[temp.idx, ]






# carry out a Poisson.test() for the rest to check for fine-scale heterochiasmy.
set.seed(1019)
p.value.vec <- NULL
rate.ratio.vec <- NULL
count.diff.vec <- NULL
null.lambda.vec <- NULL
for (i in unique(reshaped.co.count.df$WIN.IDX)) {
  subset.reshaped.co.count.df <- reshaped.co.count.df[reshaped.co.count.df$WIN.IDX == i, ]
  
  # remove individuals with NA from the dataset entirely. 
  if (any(is.na(subset.reshaped.co.count.df$CO_COUNT))) {
    subset.reshaped.co.count.df <- subset.reshaped.co.count.df[-(which(is.na(subset.reshaped.co.count.df$CO_COUNT))), ]
  }
  
  fem.vec <- subset.reshaped.co.count.df$CO_COUNT[subset.reshaped.co.count.df$SEX == 'F']
  male.vec <- subset.reshaped.co.count.df$CO_COUNT[subset.reshaped.co.count.df$SEX == 'M']
  fem.count <-  ifelse(length(fem.vec) > length(male.vec), sum(sample(x=fem.vec, size=length(male.vec), replace=FALSE)), sum(fem.vec))
  male.count <- sum(male.vec)
  count.diff <- fem.count-male.count
  # fem.count <-  ifelse(length(fem.vec) > length(male.vec), round(mean(sample(x=fem.vec, size=length(male.vec)), replace=FALSE)), round(mean(fem.vec))) 
  # male.count <- round(mean(male.vec))
  
  e.pois.mod <- poisson.exact(c(fem.count, male.count), T = 1, r = 1, alternative = "two.sided", tsmethod='central', plot=FALSE, midp=TRUE)
  mid.p <- e.pois.mod$p.value
  tmp.ratio <- e.pois.mod$estimate
  rate.ratio <- tmp.ratio
  # rate.ratio <- ifelse(tmp.ratio < 1, -(1/tmp.ratio), tmp.ratio)
  null.lambda <- e.pois.mod$parameter
  
  p.value.vec <- c(p.value.vec, mid.p)
  rate.ratio.vec <- c(rate.ratio.vec, rate.ratio)
  count.diff.vec <- c(count.diff.vec, count.diff)
  null.lambda.vec <- c(null.lambda.vec, null.lambda)
}

hist(p.value.vec, breaks=1000)
if (nrow(binned.cos.by.par.df.5) == length(null.lambda.vec) & nrow(binned.cos.by.par.df.5) == length(rate.ratio.vec)) {
  binned.cos.by.par.df.6 <- cbind(binned.cos.by.par.df.5, RATE_RATIO=rate.ratio.vec, NULL_LAMBDA=null.lambda.vec)
} else {
  stop("ERROR !!! B")
}




# Storey-Tibshirani procudure for FDR 
qvalue.obj <- qvalue(p=p.value.vec)
hist(qvalue.obj$pvalues, breaks=100)
plot(qvalue.obj$qvalues, qvalue.obj$pvalues)
# hist(qvalue.obj)
# plot(qvalue.obj)
sum(qvalue.obj$qvalues < 0.2)
fdr.rates <- seq(0, 1, by = 0.001)
num.sig.windows <- sapply(fdr.rates, function(y, x = qvalue.obj$qvalues){
  z <- sum(x < y, na.rm=TRUE)
  return(z)
})
plot(num.sig.windows ~ fdr.rates, main = 'Number of windows with heterochiasmy vs. FDR', xlab = 'FDR (q-value)', ylab = 'Number of windows significant', type = 'l', lwd = 2)



# Benjamini-Hochberg procedure for FDR
bh.p.value.vec <- p.adjust(p.value.vec, method='BH', n=length(p.value.vec))
# sum(bh.p.value.vec < 0.2)
hist(bh.p.value.vec, breaks=100)
plot(bh.p.value.vec, p.value.vec)

fdr.rates <- seq(0, 1, by = 0.001)
num.sig.windows <- sapply(fdr.rates, function(y, x = bh.p.value.vec){
  z <- sum(x < y, na.rm=TRUE)
  return(z)
})
plot(num.sig.windows ~ fdr.rates, main = 'Number of windows with heterochiasmy vs. FDR', xlab = 'FDR (q-value)', ylab = 'Number of windows significant', type = 'l', lwd = 2)

# Setting the FDR cutoff and selecting windows with recombination bias.
binned.cos.by.par.df.5[which(bh.p.value.vec < 0.25),]
nrow(binned.cos.by.par.df.5[which(bh.p.value.vec < 0.25),]) #16


# plotting different FDR values (TS & BH) against each other to understand the relationship between them.
plot(bh.p.value.vec, qvalue.obj$qvalue)
# sum(b.h.fdr.win.idx. %in% s.t.fdr.win.idx.)
# binned.cos.by.par.df.5[binned.cos.by.par.df.5$WIN.IDX %in% b.h.fdr.win.idx., ]



# Storey-Tibshirani FDR cutoffs
if (length(binned.cos.by.par.df.5$WIN.IDX) == length(qvalue.obj$qvalues)) {
  s.t.fdr.win.idx. <- binned.cos.by.par.df.5$WIN.IDX[which(qvalue.obj$qvalues < 0.25)]
  binned.cos.by.par.df.7 <- cbind(binned.cos.by.par.df.6, ST.FDR=qvalue.obj$qvalues)
} else {
  stop("ERROR!!!  A" )
}

# Benjamini-Hochberg cutoffs
if (length(binned.cos.by.par.df.5$WIN.IDX) == length(bh.p.value.vec)) {
  b.h.fdr.win.idx <- binned.cos.by.par.df.5$WIN.IDX[which(bh.p.value.vec < 0.25)]
  # length(b.h.fdr.win.idx )
  binned.cos.by.par.df.8 <- cbind(binned.cos.by.par.df.7, BH.FDR=bh.p.value.vec)
  
} else {
  stop("ERROR!!!  A" )
}
head(binned.cos.by.par.df.8)
nrow(binned.cos.by.par.df.8)
# write.csv(binned.cos.by.par.df.8, paste0(statistics.input.path, "genomewide_heterochiasmy_final_df.csv"), quote = FALSE, row.names = FALSE)
binned.cos.by.par.df.8 <- read.csv(file = paste0(statistics.input.path, "genomewide_heterochiasmy_final_df.csv"), header = TRUE, stringsAsFactors = FALSE )
colnames(binned.cos.by.par.df.8) <- colnames(binned.cos.by.par.df.8) %>% str_replace("^X", '')




# Following commented chunk moved to a separate script for plotting (Script S_17b_2)
# # plotting the GLM estimates for SEX effect as well as the CO counts for each sex to look at the pattern. 
# # jpeg(file = paste0(statistics.output.path, "genomewide_fine_scale_heterochiasmy.jpg"), width = 2000, height = 1000, quality = 5000)
# pdf(file = paste0(statistics.output.path, "genomewide_fine_scale_heterochiasmy.pdf"), width = 20, height = 10)
# # layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
# #                 rep(20, 19),
# #                 rep(21, 19), 
# #                 rep(22, 19)),nrow=4, byrow=TRUE))
# # layout.show(n=22)
# layout(matrix(c(1:19, rep(20, 19), rep(21, 19)),nrow=3, byrow=T))
# # layout.show(n=21)
# 
# chrom.size.df <- read.csv(file = paste0(statistics.input.path, "chrom.size.df_Stetv2.csv"), header = TRUE, stringsAsFactors = FALSE)
# lg.names <- c("LG-I", "LG-II", "LG-III", "LG-IV", "LG-V", "LG-VI", "LG-VII", "LG-VIII", "LG-IX", "LG-X", "LG-XI", "LG-XII", "LG-XIII", "LG-XIV",
#               "LG-XV", "LG-XVI", "LG-XVII", "LG-XVIII", "LG-XIX")
# # In order to apply the linear models the data needs to be arranged in a different table format.
# chrom.reshaped.co.count.df <- NULL
# for (j in 1:ncol(num.COs.df)) {chrom.reshaped.co.count.df <- rbind(chrom.reshaped.co.count.df, cbind(chrom.size.df, rep(colnames(num.COs.df)[j]), num.COs.df[,j], num.ind.df[,j]))}
# colnames(chrom.reshaped.co.count.df) <- c(colnames(chrom.size.df), "PAR", "CO_COUNT", "FAM_SIZE")
# chrom.reshaped.co.count.df$SEX <- as.factor(rep(c("F", "M"), each = 133))
# chrom.reshaped.co.count.df$CHROM <- as.factor(chrom.reshaped.co.count.df$CHROM)
# # sum(chrom.reshaped.co.count.df$CO_COUNT) # 38846
# chrom.reshaped.co.count.df$CO_RATE_PER_MEIOSIS <- chrom.reshaped.co.count.df$CO_COUNT / chrom.reshaped.co.count.df$FAM_SIZE
# 
# op<-par(mar=c(4.1,3.8,4.1,0))
# col.vec <- c("red", "blue")
# for (chrom.name in chrom.name.vec) {
#   y_lab<-ifelse(chrom.name==chrom.name.vec[1], "CO rate per meiosis", '')
#   boxplot(chrom.reshaped.co.count.df$CO_RATE_PER_MEIOSIS[chrom.reshaped.co.count.df$CHROM == chrom.name] ~ chrom.reshaped.co.count.df$SEX[chrom.reshaped.co.count.df$CHROM == chrom.name], 
#           main = lg.names[chrom.name.vec==chrom.name], ylab=y_lab, xlab='',
#           frame.plot=F,
#           col = col.vec)
# }
# par(op)
# # op <- par(mfrow=c(3, 1), mar=c(5.1, 5.1, 0.5, 2.1), cex.lab=1.2, cex.axis=1.2)
# 
# # # 1st plot
# # op <- par(mar=c(5.1, 5.1, 0, 2.1), cex.lab=1.2, cex.axis=1.2)
# # plot(binned.cos.by.par.df.8$WIN.IDX, -log10(binned.cos.by.par.df.8$BH.FDR), type='p', pch=16, cex=0.75, 
# #      # xlab="genomic window index",
# #      xlab="",
# #      ylab="-log10(q-value)", col='white', ylim=c(0, 1.5))
# # for (chrom.name in chrom.name.vec) {
# #   col.chrom <- ifelse(which(chrom.name.vec %in% chrom.name)%%2 == 1, 'grey48', 'grey25')
# #   subset.binned.df <- binned.cos.by.par.df.8[binned.cos.by.par.df.8$CHROM == chrom.name, ]
# #   points(subset.binned.df$WIN.IDX, -log10(subset.binned.df$BH.FDR), pch=16, cex=0.7, col=col.chrom)
# #   abline(v=subset.binned.df$WIN.IDX[nrow(subset.binned.df)] ,col='grey', lty=2)
# #   # text(x=round(mean(subset.binned.df$WIN.IDX)), y=1.4, labels=chrom.name)
# # }
# # # points(b.h.fdr.win.idx, -log10(binned.cos.by.par.df.8$BH.FDR[binned.cos.by.par.df.8$WIN.IDX %in% b.h.fdr.win.idx]), col='black')
# # abline(h=-log10(0.2), col='black', lty=2)
# # par(op)
# 
# 
# # 2nd plot
# # par(new=TRUE)
# # plot(binned.cos.by.par.df.8$WIN.IDX, -log10(binned.cos.by.par.df.8$RATE_RATIO), type='h', col='lightgrey', axes=FALSE, bty="n", xlab='', ylab='')
# op <- par(mar=c(5.1, 5.1, 0, 2.1), cex.lab=2, cex.axis=2, xaxs='i')
# plot(binned.cos.by.par.df.8$WIN.IDX, -log10(binned.cos.by.par.df.8$RATE_RATIO), type='h', col='grey57', 
#      # xlab='genomic window index', 
#      xlim=c(0, binned.cos.by.par.df.8$WIN.IDX[length(binned.cos.by.par.df.8$WIN.IDX)]+1),
#      xlab='',
#      ylab='-log10(Fem CO / Male CO)')
# 
# lines(b.h.fdr.win.idx, -log10(binned.cos.by.par.df.8$RATE_RATIO[binned.cos.by.par.df.8$WIN.IDX %in% b.h.fdr.win.idx]), col='blue', type='h', lwd=2)
# lines(b.h.fdr.win.idx[5], -log10((binned.cos.by.par.df.8$RATE_RATIO[binned.cos.by.par.df.8$WIN.IDX %in% b.h.fdr.win.idx])[6]), col='red', type='h', lwd=2)
# 
# # points(binned.cos.by.par.df.8$WIN.IDX[binned.cos.by.par.df.8$WIN.IDX %in% b.h.fdr.win.idx], -log10(binned.cos.by.par.df.8$RATE_RATIO[binned.cos.by.par.df.8$WIN.IDX %in% b.h.fdr.win.idx]), pch=24, cex=0.5, col='red')
# abline(h=0, lty=2, col='grey')
# for (chrom.name in chrom.name.vec) {
#   subset.binned.df <- binned.cos.by.par.df.8[binned.cos.by.par.df.8$CHROM == chrom.name, ]
#   abline(v=subset.binned.df$WIN.IDX[nrow(subset.binned.df)] ,col='grey', lty=2)
# }
# # axis(side=4, at=pretty(range(-log10(binned.cos.by.par.df.8$RATE_RATIO), na.rm=TRUE)))
# # mtext("-log10(q-value)", side=2, line=3, cex=1)
# # mtext("-log10(Fem CO / Male CO)", side=4, line=3, cex=1)
# par(op)
# 
# 
# 
# # 3rd plot
# op <- par(mar=c(5.1, 5.1, 0, 2.1), cex.lab=2, cex.axis=2, xaxs='i')
# for (par.name in par.name.vec[1:14]) { # par.name <- par.name.vec[1]
#   colr <- ifelse(par.name %in% par.name.vec[1:7], yes=rgb(255, 0, 0, max=255, alpha=80, names="red50"), no=rgb(0, 0, 255, max=255, alpha=80, names="blue50"))
#   colr <- 'white'
#   
#   if (par.name == par.name.vec[1]) {
#     plot(binned.cos.by.par.df.8$WIN.IDX, binned.cos.by.par.df.8[, par.name], col=colr, ylim=c(0, 30), type='l',
#          xlab='genomic window index', ylab='cross-over count')
#     points(b.h.fdr.win.idx, rep(0, length(b.h.fdr.win.idx)), col='black', pch=17)
#   } else {
#     lines(binned.cos.by.par.df.8$WIN.IDX, binned.cos.by.par.df.8[, par.name], col=colr)
#   }
# }
# 
# 
# for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
#   
#   subset.binned.df <- binned.cos.by.par.df.8[binned.cos.by.par.df.8$CHROM == chrom.name, ]
#   splitAt <- function(x, pos) {unname(split(x, cumsum(seq_along(x) %in% pos)))} 
#   defined.gap <- 3
#   win.indices <- subset.binned.df$WIN.IDX
#   idx.diff.vec <- diff(win.indices)
#   diff.idx <- which(idx.diff.vec > defined.gap)
#   gap.sep.list <- splitAt(win.indices, (diff.idx + 1))
#   
#   for (par.name in par.name.vec) { # par.name <- par.name.vec[1]
#     colr <- ifelse(par.name %in% par.name.vec[1:7], yes=rgb(255, 0, 0, max=255, alpha=80, names="red50"), no=rgb(0, 0, 255, max=255, alpha=80, names="blue50"))
#     
#     for (l in 1:length(gap.sep.list)) {
#       segmnt <- gap.sep.list[[l]]
#       tmp.idx <- which(binned.cos.by.par.df.8$WIN.IDX %in% segmnt)
#       lines(binned.cos.by.par.df.8$WIN.IDX[tmp.idx], binned.cos.by.par.df.8[tmp.idx, par.name], col=colr, ylim=c(0, 30), type='l', xlab='genomic window index', ylab='cross-over count')
#     }
#   }
#   
#   abline(v=subset.binned.df$WIN.IDX[nrow(subset.binned.df)] ,col='grey', lty=2)
#   
#   for (l in 1:length(gap.sep.list)) {
#     segmnt <- gap.sep.list[[l]]
#     tmp.idx <- which(binned.cos.by.par.df.8$WIN.IDX %in% segmnt)
#     lines(binned.cos.by.par.df.8$WIN.IDX[tmp.idx], apply(binned.cos.by.par.df.8[tmp.idx, par.name.vec[1:7]], 1, mean, na.rm=TRUE), col='red', lwd=1.5)
#     lines(binned.cos.by.par.df.8$WIN.IDX[tmp.idx], apply(binned.cos.by.par.df.8[tmp.idx, par.name.vec[8:14]], 1, mean, na.rm=TRUE), col='blue', lwd=1.5)
#   }
# }
# par(op)
# dev.off()









# 
# for (par.name in par.name.vec[1:14]) { # par.name <- par.name.vec[1]
#   colr <- ifelse(par.name %in% par.name.vec[1:7], yes=rgb(255, 0, 0, max=255, alpha=80, names="red50"), no=rgb(0, 0, 255, max=255, alpha=80, names="blue50"))
#   for (chrom.name in chrom.name.vec) {
#     
#   }
#   
#   if (par.name == par.name.vec[1]) {
#     for (chrom.name in chrom.name.vec) {
#      
#       
#       for (segmnt in gap.sep.list) {
#         segmnt <- unlist(segmnt)
#         plot(binned.cos.by.par.df.8$WIN.IDX, binned.cos.by.par.df.8[segmnt, par.name], col=colr, ylim=c(0, 30), type='l', xlab='genomic window index', ylab='cross-over count')
#       }
#       
#     }
# 
#     # plot(binned.cos.by.par.df.8$WIN.IDX, binned.cos.by.par.df.8[, par.name], col=colr, ylim=c(0, 30), type='l', xlab='genomic window index', ylab='cross-over count')
#   } else {
#     
#     lines(binned.cos.by.par.df.8$WIN.IDX, binned.cos.by.par.df.8[, par.name], col=colr)
#   }
# }
# 
# 
# 
# 
# lines(binned.cos.by.par.df.8$WIN.IDX, apply(binned.cos.by.par.df.8[, par.name.vec[1:7]], 1, mean, na.rm=TRUE), col='red', lwd=2)
# lines(binned.cos.by.par.df.8$WIN.IDX, apply(binned.cos.by.par.df.8[, par.name.vec[8:14]], 1, mean, na.rm=TRUE), col='blue', lwd=2)
# 
# # select.idx <- which(win.based.hetchim.df.2$WIN.IDX %in% win.based.hetchim.df.3$WIN.IDX)
# # points(win.based.hetchim.df.2$WIN.IDX[select.idx], apply(win.based.hetchim.df.2[select.idx, par.name.vec[8:14]], 1, mean))
# par(op)





# extracting windows for which males have significantly higher CO counts compared to females
male.hi.wins <- which(binned.cos.by.par.df.8$WIN.IDX %in% b.h.fdr.win.idx & binned.cos.by.par.df.8$RATE_RATIO < 1)

# windows for which male-female CO rate difference is not significant
avg.wins <- which((binned.cos.by.par.df.8$BH.FDR > 0.4 & binned.cos.by.par.df.8$BH.FDR < 0.8) & binned.cos.by.par.df.8$RATE_RATIO < 1)
# avg.wins <- which((!binned.cos.by.par.df.8$WIN.IDX %in% b.h.fdr.win.idx) & binned.cos.by.par.df.8$RATE_RATIO < 1)

# make a dataframe for analaysis
tmp.1.df <- cbind(binned.cos.by.par.df.8[male.hi.wins, ], CATEG=rep("HI", length(male.hi.wins)))
tmp.2.df <- cbind(binned.cos.by.par.df.8[avg.wins, ], CATEG=rep("AVG", length(avg.wins)))
gen.feat.selected.df <- rbind(tmp.1.df, tmp.2.df)
gen.feat.selected.df <- gen.feat.selected.df[order(gen.feat.selected.df$WIN.IDX), ]

# read-in the genomic features for these windows
win.size <- '960kb'
gen.feat.df <- read.csv(file = paste0(gen.feat.output.path, win.size, "_genomic_features_Stet_v2_raw.csv"), header = TRUE, stringsAsFactors = FALSE )
colnames(gen.feat.df)[which(colnames(gen.feat.df) == 'WIN_IDX')] <- 'WIN.IDX'
subset.gen.feat.df <- gen.feat.df[gen.feat.df$WIN.IDX %in% gen.feat.selected.df$WIN.IDX, ]

# produce a common data table with all the data together.
if (all(gen.feat.selected.df$START == subset.gen.feat.df$START)) {
  combined.df <- cbind(gen.feat.selected.df, subset.gen.feat.df[, 4:18])
} else {
  stop("ERROR !!! C")
}


# visualize the difference in genomic features using boxplots
pdf(file = paste0(statistics.output.path, "genomic_correlates_and_male_bias.pdf"), width = 4, height =8)
layout(matrix(rep(1:5, each=2), nrow=5, byrow=T))
# layout.show(n=5)
op <- par(mar=c(4.1, 3.1, 2.1, 2.1))
# boxplot(AVG_CO ~ CATEG, data=combined.df, xlab="Male CO count differential", ylab="Average CO count", col=c("blue", "grey"), border="black")
# stripchart(AVG_CO ~ CATEG, vertical = TRUE, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
# boxplot(FEM_CO ~ CATEG, data=combined.df, xlab="Male CO count differential", ylab="Fem cumulative COs", col=c("blue", "grey"), border="black")
# stripchart(FEM_CO ~ CATEG, vertical = TRUE, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
# boxplot(MALE_CO ~ CATEG, data=combined.df, xlab="Male CO count differential", ylab="Male cumulative COs", col=c("blue", "grey"), border="black")
# stripchart(MALE_CO ~ CATEG, vertical = TRUE, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
# op <- par(mfrow=c(2, 3))
box_names<- c('male biased', 'background')
colrs <- c("blue", "grey")
boxplot(GENE ~ CATEG, data=combined.df, ylab="Gene count", col=colrs, border="black", names=box_names, las=2, ylim=c(0,400), horizontal = T)
stripchart(GENE ~ CATEG, vertical = F, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
boxplot(COPIA ~ CATEG, data=combined.df, ylab="Copia count", col=colrs, border="black", names=box_names, las=2, ylim=c(0,400), horizontal = T)
stripchart(COPIA ~ CATEG, vertical = F, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
boxplot(GYPSY ~ CATEG, data=combined.df, ylab="Gypsy count", col=colrs, border="black", names=box_names, las=2, ylim=c(0,400), horizontal = T)
stripchart(GYPSY ~ CATEG, vertical = F, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
boxplot(SIMP_REPEAT ~ CATEG, data=combined.df, ylab="Simple repeats count", col=colrs, border="black", names=box_names, las=2, ylim=c(0,400), horizontal = T)
stripchart(SIMP_REPEAT ~ CATEG, vertical = F, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
# boxplot(pct_at ~ CATEG, data=combined.df, ylab="pct AT", col=colrs, border="black", names=box_names, las=2, horizontal = T)
# stripchart(pct_at ~ CATEG, vertical = F, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
boxplot(pct_gc ~ CATEG, data=combined.df, ylab="pct GC", col=colrs, border="black", names=box_names, las=2, horizontal = T)
stripchart(pct_gc ~ CATEG, vertical = F, data = combined.df, method = "jitter", add = TRUE, pch = 20, col =c('black'))
par(op)
dev.off()



# compare the means of the two categories using a t-test
t.test(pct_at ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(pct_gc ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(GENE ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(COPIA ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(GYPSY ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(SIMP_REPEAT ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)

t.test(log10(GENE + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(log10(COPIA + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(log10(GYPSY + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(log10(SIMP_REPEAT + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)




# compare the means of the two categories using a non parametric test (Wilcoxon-Mann-Whitney test)
wilcox.test(formula=pct_at ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=pct_gc ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=GENE ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=COPIA ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=GYPSY ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=SIMP_REPEAT ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
# 
# wilcox.test(formula=log10(GENE + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
# wilcox.test(formula=log10(COPIA + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
# wilcox.test(formula=log10(GYPSY + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
# wilcox.test(formula=log10(SIMP_REPEAT + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################