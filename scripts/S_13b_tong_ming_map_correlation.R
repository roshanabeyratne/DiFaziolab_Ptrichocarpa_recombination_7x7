#!/usr/bin/R
setwd("Google\ Drive/Shared\ drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS, for the recombination landscape analysis for 7x7 cross
# DATE: September 24, 2020
# PURPOSE: Unerstand the average CO-region size
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
input.path <- "./statistics/"
statistics.output.path <- "./statistics/gross_stat_plots/"
########################################################################################################################################################################################################

num.vec <- c(11,8,36,1,18,7,18,
               18,17,17,19,12,19,16,
               18,19,18,18,11,32,19,
               20,15,17,18,18,18,17,
               18,17,16,18,15,18,24,
               17,20,19,25,13,4,18,
               21,13,11,14,16,14,15)
hist(num.vec, breaks = 14, xlab = "Number of offspring", main = NULL, col = 'grey')
plot(density(num.vec))
var(num.vec)


tong.ming <- c(270.58, 173.73, 142.07, 140.46, 179.27, 147.41, 96.31, 160.62, 113.89, 137.40, 103.03, 78.41, 83.31, 139.03, 93.49, 107.46, 104.04, 108.02, 101.3)
female.map <- c(264.9, 156.1, 126.1, 132.3, 162.2, 175.8, 89.8, 114.5, 94.1, 130.8, 92, 87.2, 93.6, 92.1, 87.9, 86.3, 87.1, 98.4, 68.3)
male.map <- c(304.6, 172.6, 137.8, 137.1, 157.8, 173.9, 100.8, 126.7, 96.9, 151.2, 103.7, 97.4, 108.9, 113.1, 93.6, 103.2, 91.9, 96.6, 90.4)


pdf(file=paste0(statistics.output.path, "Tong_Ming_comparison.pdf"), width=8, height=8)
plot(tong.ming, male.map, typ = "p", col = 'blue', pch = 19, 
     xlab='Genetic map length (cM) T.M. Yin et al., 2004',
     ylab='Genetic map length (cM) 7x7 map averages')
points(tong.ming, female.map, typ = "p", col = 'red', pch = 19)
lm.tong <- lm(tong.ming~ female.map)
print(cor.test(tong.ming, female.map, method = "pearson"))  
print(cor.test(tong.ming, male.map, method = "pearson"))  

abline(lm(tong.ming~ female.map), col = 'red', lty=2)
abline(lm(tong.ming~ male.map), col = 'blue', lty=2)
dev.off()
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################