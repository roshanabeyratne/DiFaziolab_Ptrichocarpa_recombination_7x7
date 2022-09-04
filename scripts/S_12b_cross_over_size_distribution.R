#!/usr/bin/R
setwd("/Users/cabeyrat/Google Drive File Stream/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS, for the recombination landscape analysis for 7x7 cross
# DATE: January 15, 2021
# USAGE: Input file include a detailed list of all the cross-over events observed for each focal parent for each offspring for all the chromosomes. The input is produced with script 22b_... 

# Following is the input file format used 
# CHROM  PAR   IND CO_SIZE DWNSTRM_FLNK UPSTRM_FLNK
# Chr01 1863 24708  103140      4675047     4778187
# Chr01 1863 24708   30914     26411972    26442886
# Chr01 1863 24708   61845     34445658    34507503
# Chr01 1863 24708    5612     46505278    46510890
# Chr01 1863 24709  164382       726848      891230
# Chr01 1863 24709  153359      2513166     2666525

# PURPOSE: Unerstand the average CO-region size
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
input.path <- "./statistics/"
statistics.output.path <- "./statistics/gross_stat_plots/"
########################################################################################################################################################################################################
# Read-n the detailed CO dataframe
CO.df <- read.csv(paste0(input.path, "co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)

# this code block is to identify the breakdown of individuals in full-sib families for the Table S1 in the manuscript
ind_829 <- unique(CO.df$IND)
parentage_v4 <- read.table(file='../phenotypes/CRA/7x7_parentage_v4.txt', stringsAsFactors=F, header=T)
parentage_829 <- parentage_v4[parentage_v4$OldClone %in% ind_829, ]
table(parentage_829$matched_father, parentage_829$matched_mother)
sum(colSums(table(parentage_829$matched_father, parentage_829$matched_mother)))


cleaned.co.size <- CO.df$CO_SIZE[-(which(CO.df$CO_SIZE == 0))]

max.x <- max(cleaned.co.size)
# jpeg(filename = paste0(statistics.output.path, "CO_region_size_distribution.jpg"), width = 480, height = 480, quality = 500)
pdf(file=paste0(statistics.output.path, "CO_region_size_distribution.pdf"), width = 5, height = 5)
hist(cleaned.co.size, breaks = 10000, 
     xlab = "Size of CO-region (bps)", ylab = "frequency", main = NULL
     ,xlim = c(0, max.x)
     )
text(x = c(3000000), y = c(200), labels = (paste0("mean CO size = ", round(mean(CO.df$CO_SIZE))/1000, "Kb")), pos = 3)
text(x = c(3000000), y = c(300), labels = (paste0("median CO size = ", round(median(CO.df$CO_SIZE))/1000, "Kb")), pos = 3)
dev.off()

pdf(file=paste0(statistics.output.path, "CO_region_size_distribution_zoom.pdf"), width = 5, height = 5)
hist(cleaned.co.size, breaks = 10000, 
     xlab = "Size of CO-region (bps)", ylab = "frequency", main = NULL
     ,xlim = c(0, 200000)
)
# text(x = c(3000000), y = c(200), labels = (paste0("mean CO size = ", round(mean(CO.df$CO_SIZE))/1000, "Kb")), pos = 3)
# text(x = c(3000000), y = c(300), labels = (paste0("median CO size = ", round(median(CO.df$CO_SIZE))/1000, "Kb")), pos = 3)
dev.off()


quantile(cleaned.co.size, c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
# 25%       50%       75%       90% 
#   13722.25  30715.50  55600.50 103432.80 
mean(cleaned.co.size)
median(cleaned.co.size)
max(cleaned.co.size, na.rm = T)
# [1] 5338316
#########################################################################################################################################################################################################
# END
#########################################################################################################################################################################################################