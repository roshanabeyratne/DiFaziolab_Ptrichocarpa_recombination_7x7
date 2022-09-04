#!/usr/bin/R
setwd("/Users/cabeyrat/Google Drive File Stream/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping/statistics")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS, for the recombination landscape analysis for 7x7 cross
# DATE: January 15, 2021
# USAGE: Input files include average CO counts by parent and chromosome as well as a dataframe for chromosome physical size and genetc map size. 

# Following is the input file format used with average cross-over counts per half-sib family per chromosome
#          V1        V2        V3        V4        V5        V6        V7       V8        V9       V10       V11       V12       V13       V14
# 1  2.4900000 2.0677966 2.9022556 2.8225806 2.2015504 2.7203390 2.2336449 3.338583 3.5137615 3.1777778 2.5877193 2.6796117 1.4521739 3.1507937
# 2  1.5600000 1.4237288 1.6917293 1.6290323 1.2093023 1.6271186 1.2429907 1.897638 1.7339450 2.0814815 1.4210526 1.4466019 1.4347826 1.7539683
# 3  1.2600000 0.9406780 1.5413534 1.3790323 1.0387597 1.4661017 1.0467290 1.850394 1.6972477 1.3851852 1.3859649 1.1553398 1.2173913 1.5158730
# 4  1.5500000 1.0932203 1.4736842 1.3306452 1.2480620 1.4237288 1.1495327 1.559055 1.6238532 1.3703704 1.3421053 1.3689320 1.1652174 1.4365079
# 5  1.6400000 1.1440678 1.7218045 1.6451613 1.4806202 1.7542373 1.2149533 1.755906 1.7889908 1.5777778 1.4385965 1.4368932 1.4347826 1.5952381
# 6  1.6700000 1.3813559 1.9473684 1.8629032 1.3643411 1.9491525 1.3738318 2.314961 1.9908257 1.8962963 1.6929825 1.5922330 1.7043478 1.7380952
# 7  1.0400000 0.7627119 0.9774436 1.0725806 0.8604651 0.9067797 0.8971963 1.015748 1.0642202 1.0814815 0.8245614 0.8543689 0.9217391 1.0158730
# 8  1.3400000 1.1525424 1.4887218 1.2741935 0.9457364 1.0338983 0.9626168 1.496063 1.4128440 1.2666667 1.1929825 1.0873786 1.0869565 1.3333333
# 9  1.0400000 0.8220339 1.0526316 1.0161290 0.7286822 0.9406780 0.7383178 1.133858 1.0275229 0.9777778 0.9473684 0.9514563 0.8869565 0.9841270
# 10 1.4400000 0.9745763 1.4360902 1.3306452 1.0465116 1.5508475 1.0747664 1.511811 1.7431193 1.5777778 1.2719298 1.2135922 1.2521739 1.7142857
# 11 1.1400000 0.9067797 1.1654135 0.9274194 0.8759690 1.1186441 0.8130841 1.236220 1.0550459 1.1111111 1.0087719 0.8252427 0.9304348 1.1190476
# 12 1.0100000 0.7542373 0.9248120 0.8870968 0.8527132 1.1694915 0.7196262 1.094488 1.0825688 0.9777778 0.9736842 0.8737864 0.8347826 1.0079365
# 13 1.1500000 0.9237288 1.0902256 0.9516129 0.8604651 1.1355932 0.8130841 1.110236 1.1651376 1.0962963 1.0964912 0.8932039 0.9561404 1.1428571
# 14 0.9200000 0.8135593 1.1353383 1.0967742 0.8062016 0.9322034 0.8130841 1.301587 1.4220183 1.2296296 0.9824561 0.8543689 1.1304348 1.0555556
# 15 0.8787879 0.7542373 0.9473684 0.9516129 0.8604651 0.8813559 0.8317757 1.118110 0.9449541 1.1703704 0.9115044 0.7864078 0.7739130 0.9682540
# 16 0.9600000 0.6440678 0.9548872 0.8709677 0.8139535 1.0932203 0.7102804 1.039370 1.0917431 1.0518519 0.9385965 0.9611650 0.8695652 1.0952381
# 17 0.9100000 0.6779661 1.0601504 0.9596774 0.6666667 0.9322034 0.7663551 1.078740 1.1559633 1.0223881 0.8771930 0.8155340 0.9122807 1.0634921
# 18 1.1500000 0.8305085 1.0601504 1.0806452 0.7906977 1.0169492 0.8504673 1.149606 1.2110092 1.0074074 0.9649123 0.9514563 0.8434783 0.9444444
# 19 0.6000000 0.6186441 0.8571429 0.9032258 0.6821705 0.9406780 0.6542056 1.110236 1.0458716 0.7925926 0.9122807 0.7184466 0.8347826 0.9841270

# PURPOSE: Model average CO-counts given the chromosome LENGTH and SEX for all the halfsib families
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
library(car)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
input.path <- "./statistics/"
statistics.output.path <- "./statistics/gross_stat_plots/"
########################################################################################################################################################################################################
# The dataframe with Chromosome sizes in base-pairs and genetic-map size
# chrom.size.df <- read.csv(file = paste0(input.path, "chrom.size.df.csv"), header = TRUE, stringsAsFactors = FALSE)
chrom.size.df <- read.csv(file = paste0(input.path, "chrom.size.df_Stetv2.csv"), header = TRUE, stringsAsFactors = FALSE)
# detailed CO dataframe
co.count.df <-  read.csv(paste0(input.path, "co.count.df.csv"), stringsAsFactors = FALSE, header = TRUE)

# These are the standardized average counts per half-sib family. The E[CO counts per individual]  are multiplied by sqrt(IND) and devided by the sqrt(E[variance of CO counts]), as per CLT.
# std.co.count.df <-  read.csv(paste0(statistics.output.path, "std.co.count.df.csv"), stringsAsFactors = FALSE, header = TRUE)
# std.co.count.df <-  read.csv(paste0("std.co.count.df.csv"), stringsAsFactors = FALSE, header = TRUE)
# Use this if std.co.count.df is the dataset of choice.
# avg.count.df <- std.co.count.df

# Use this if regular co.count.df is the dataset of choice. 
avg.count.df <- co.count.df


colnames(avg.count.df) <- par.name.vec
row.names(avg.count.df) <- chrom.name.vec


# Analyzing CO counts averaged by sex
avg.count.df.by.sex <- cbind(c(apply(avg.count.df[1:7], 1, mean), apply(avg.count.df[8:14], 1, mean)), (rep(c("F", "M"), each = 19)), rep(as.numeric(chrom.size.df$LENGTH), 2))
avg.count.df.by.sex <- as.data.frame(avg.count.df.by.sex, stringsAsFactors = FALSE)
colnames(avg.count.df.by.sex) <- c("AVG_CO_COUNT", "SEX", "LENGTH")
avg.count.df.by.sex$AVG_CO_COUNT <- as.numeric(avg.count.df.by.sex$AVG_CO_COUNT)
avg.count.df.by.sex$LENGTH <- as.numeric(avg.count.df.by.sex$LENGTH)
avg.count.df.by.sex$SEX <- as.factor(avg.count.df.by.sex$SEX)


# General linear model for chromosome length as the covariate. Y ~ xb + e
lm.avg <- lm(AVG_CO_COUNT~LENGTH, data = avg.count.df.by.sex)
summary(lm.avg)

# check for distribution of normality of residuals for downstream statistical analysis
hist(residuals(lm.avg), breaks = 10)
shapiro.test(residuals(lm.avg)) 
# data:  residuals(lm.avg)
#  W = 0.98017, p-value = 0.7238
plot(fitted(lm.avg), residuals(lm.avg))

# levene's test to establish homogeneity of variance across groups of SEX to be tested.
avg.count.df.by.sex$resids <- residuals(lm.avg)
leveneTest(resids ~ SEX, avg.count.df.by.sex)

# t-test for groups of SEX after regressing LENGTH of chromosome out
t.test(residuals(lm.avg)[1:19], residuals(lm.avg)[20:38], paired = TRUE, conf.level = 0.95)


# plotting the data 
plot(avg.count.df.by.sex$LENGTH, avg.count.df.by.sex$AVG_CO_COUNT, col = 'white', xlab = "Chromosome size (bp)", ylab = "Average Count of Cross-overs per individual in half-sib family")
points(avg.count.df.by.sex$LENGTH[avg.count.df.by.sex$SEX == "F"], avg.count.df.by.sex$AVG_CO_COUNT[avg.count.df.by.sex$SEX == "F"], col = 'red', pch = 19)
points(avg.count.df.by.sex$LENGTH[avg.count.df.by.sex$SEX == "M"], avg.count.df.by.sex$AVG_CO_COUNT[avg.count.df.by.sex$SEX == "M"], col = 'blue', pch = 19)



# Analyzing CO counts per individual across half-sib families using standardized data per CLT.
# reshape the co.count.df to a form that can be used for plotting as (x,y) coordinates
reshaped.co.count.df <- NULL
for (j in 1:ncol(avg.count.df)) {
  reshaped.co.count.df <- rbind(reshaped.co.count.df, cbind(chrom.size.df, rep(colnames(avg.count.df)[j]), avg.count.df[,j]))
}
colnames(reshaped.co.count.df) <- c(colnames(chrom.size.df), "PAR", "AVG_CO_COUNT")
reshaped.co.count.df$SEX <- as.factor(rep(c("F", "M"), each = 133))

# if (any(reshaped.co.count.df$AVG_CO_COUNT == 0)) {
#   stop("zeros exist!")
# } else {
#   reshaped.co.count.df$AVG_CO_COUNT <- log(reshaped.co.count.df$AVG_CO_COUNT)
# }


hist(reshaped.co.count.df$AVG_CO_COUNT, breaks = 150)
shapiro.test(reshaped.co.count.df$AVG_CO_COUNT)
lm.mod <- lm(AVG_CO_COUNT ~ SEX + LENGTH, data = reshaped.co.count.df)
# lm.mod <- lm(AVG_CO_COUNT ~ SEX, data = reshaped.co.count.df[reshaped.co.count.df$CHROM == chrom.name.vec[3],])
summary(lm.mod)

# Call:
#   lm(formula = AVG_CO_COUNT ~ SEX + LENGTH, data = reshaped.co.count.df)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.4967 -1.2302 -0.0224  1.1820  5.8818 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 9.121e+00  3.402e-01  26.808  < 2e-16 ***
#   SEXM        1.219e+00  2.338e-01   5.213 3.77e-07 ***
#   LENGTH      2.514e-07  1.421e-08  17.699  < 2e-16 ***
#   ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.906 on 263 degrees of freedom
# Multiple R-squared:  0.5641,	Adjusted R-squared:  0.5608 
# F-statistic: 170.2 on 2 and 263 DF,  p-value: < 2.2e-16
# 

shapiro.test(residuals(lm.mod))
# data:  residuals(lm.mod)
# W = 0.99167, p-value = 0.1379
plot(fitted(lm.mod), residuals(lm.mod))
hist(residuals(lm.mod))


# Checking each chromosome for which chromosome is driving differences in sex
# reshaped.co.count.df.per.chrom <- reshaped.co.count.df[reshaped.co.count.df$CHROM == chrom.name.vec[3],]
lm.mod.length.only <- lm(AVG_CO_COUNT ~ LENGTH, data = reshaped.co.count.df)
summary(lm.mod.length.only)
shapiro.test(residuals(lm.mod.length.only))
reshaped.co.count.df$LENGTH_FIT <- residuals(lm.mod.length.only)

sex.mod.lm <- lm(LENGTH_FIT ~ SEX, data = reshaped.co.count.df)
summary(sex.mod.lm)

# Fitting the model per chromosome to identify whether evidence exists in favor of sex-bias.
sex.per.chrom.mod <- lm(lm(LENGTH_FIT ~ SEX, data = reshaped.co.count.df[reshaped.co.count.df$CHROM == chrom.name.vec[10], ]))
summary(sex.per.chrom.mod)
# chromosomes 3, 10 & 19 show interesting trends of sex bias in CO counts. 




# plotting the data: regression of average CO counts across halfsib families on chromosome length
jpeg(filename = paste0(statistics.output.path, "CO_count_vs_chromosome_size.jpg"), width = 480, height = 480, quality = 500)
plot(reshaped.co.count.df$LENGTH, reshaped.co.count.df$AVG_CO_COUNT, pch = 19, cex = 0.5, col = 'black', xlab = "chromosome size (bp)", ylab = "mean count of cross-overs across half-sib families")
co.v.size.mod <- lm(AVG_CO_COUNT ~ LENGTH + SEX, data = reshaped.co.count.df) # define the linear model
abline(co.v.size.mod, lty = 5, lwd = 2, col = 'blue') # plot the regression line
# define confidence intervals based on the linear regression model
CI <- as.data.frame(predict(co.v.size.mod, interval = 'confidence'))
# add the original x-values to the dataframe
CI$length <- reshaped.co.count.df$LENGTH
# make the lwr and upr confidence interval lines
lines(CI$length[1:19], round(CI$lwr[1:19], 4), lty = 1, lwd = 2, col = 'grey')
lines(CI$length[1:19], round(CI$upr[1:19], 4), lty = 1, lwd = 2, col = 'grey')
mtext(paste0("R^2 = ", round(summary(co.v.size.mod)$adj.r.squared, 3), ";", "p-value = 2.2e-16"))
dev.off()


# The p-value for the model was obtained by the following command.
cor.test(reshaped.co.count.df$LENGTH, reshaped.co.count.df$AVG_CO_COUNT)


# caluculatig the genome-wide avg. CO count per Mb
reshaped.co.count.df$AVG_CO_PER_MB <- reshaped.co.count.df$AVG_CO_COUNT / (reshaped.co.count.df$LENGTH) / 1*(10^-6)

jpeg(filename = paste0(statistics.output.path, "CO_count_by_parent_boxplot.jpg"), width = 800, height = 450, quality = 500)
col.vec <- rep(c("red", "blue") , each = 7)
boxplot(AVG_CO_PER_MB ~ PAR, data = reshaped.co.count.df, xlab = "Half-sib family identity by focal parent", ylab = "Genoe-wide average cross-over count per Mb", col = col.vec)
dev.off()

########################################################################################################################################################################################################
# experimenting with different models
summary(sex.included.mod)
boxplot(AVG_CO_COUNT ~ SEX, data = reshaped.co.count.df, xlab = "Half-sib family identity by focal parent", ylab = "Mean cross-over count")

# test the two fit models 
anova(sex.included.mod, co.v.size.mod)
# Model 1: AVG_CO_COUNT ~ LENGTH + SEX
# Model 2: AVG_CO_COUNT ~ LENGTH
# Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1    263 12.071                                  
# 2    264 13.228 -1   -1.1572 25.213 9.477e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

residuals.length <- residuals(co.v.size.mod)
boxplot(residuals.length ~ reshaped.co.count.df$SEX, xlab = "Sex of the focal parent", ylab = "Residuals from Mean CO-count ~ physical length")
leveneTest(residuals.length ~ reshaped.co.count.df$SEX)
t.test(residuals.length ~ reshaped.co.count.df$SEX)
reshaped.co.count.df2 <- reshaped.co.count.df[reshaped.co.count.df$CHROM == chrom.name.vec[3],]
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################