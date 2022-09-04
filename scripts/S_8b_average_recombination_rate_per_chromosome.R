#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS, for the recombination landscape analysis for 7x7 cross
# DATE: September 16, 2020
# USAGE: Input files include genetic maps and table with physical vs. genetic map sizes. The marker names contain the chromosome and marker putative physical position. 

# Following is the input file / genetic map format; 
# 1 Chr06_417767 0
# 1 Chr06_5256 9.99999999868873e-06
# 1 Chr06_685931 0.72899955945013
# 1 Chr06_472891 0.729009559450128
# 1 Chr06_811815 2.18926150399718
# 1 Chr06_804902 2.18927150399718
# 1 Chr06_811828 2.18928150399718
# 1 Chr06_1049955 2.91925669384554
# 1 Chr06_1037834 2.91926669384554
# 1 Chr06_1297578 3.64925286278189
# 1 Chr06_1139116 3.64926286278189
# 1 Chr06_1071037 3.64927286278189

# ollowing is the input file chrom.size.df with physical vs. genetic map sizes
# CHROM   LENGTH GEN_DIST
# 1  Chr01 50495398      350
# 2  Chr02 25263042      210
# 3  Chr03 25000000      200
# 4  Chr04 24267058      175
# 5  Chr05 25890711      200
# 6  Chr06 27912132      250
# 7  Chr07 15610920      120
# 8  Chr08 19465468      175
# 9  Chr09 12948749      125
# 10 Chr10 22580539      175
# 11 Chr11 18501278      125
# 12 Chr12 15760353      150
# 13 Chr13 16320724      125
# 14 Chr14 18920901      150
# 15 Chr15 15278584      125
# 16 Chr16 14494368      125
# 17 Chr17 16080365      125
# 18 Chr18 16958307      150
# 19 Chr19 15942152      150

# PURPOSE: create a boxplot that displays the average recombination rate for each chromosome. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
# args <- commandArgs(TRUE)
# chrom.name <- args[1]
# par.name <- args[2]
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

genetic.map.input.path <- "./onemap/bin_markers/genetic_maps/"
rippled.genetic.map.input.path <- "./onemap/bin_markers/genetic_maps/rippled_maps/"
########################################################################################################################################################################################################
rec.rate.df  <- NULL
gen.map.stat.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  tmp.gen.map.size <- NULL
  tmp.phys.map.size <- NULL
  tmp.num.markers <- NULL
  
  for (par.name in par.name.vec) { # par.name <- "1863"
    
    if (file.exists(paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"))) {
      gen.map <- read.table(file = paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    } else {
      gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    }
    
    string <- paste0(chrom.name, "_")
    physical.distance <- gen.map[,2] %>% str_replace(string, "") %>% as.numeric()
    genetic.distance <- gen.map[,3] %>% as.numeric()
    
    total.phys.dist <- (physical.distance[length(physical.distance)] - physical.distance[1]) / 1000000 # obtain the total physical distance over which the markers spread
    total.gen.dist <- genetic.distance[length(genetic.distance) ]# obtain the total genetic map size for the markers
    total.markers <- length(genetic.distance)
    
    rec.rate.df <- rbind(rec.rate.df, c(chrom.name, par.name, (total.gen.dist/total.phys.dist)))
    tmp.gen.map.size <- c(tmp.gen.map.size, total.gen.dist)
    tmp.phys.map.size <- c(tmp.phys.map.size, total.phys.dist)
    tmp.num.markers <- c(tmp.num.markers, total.markers)
    
  }
  gen.map.stat.df <- rbind(gen.map.stat.df, c(median(tmp.gen.map.size[1:7]), median(tmp.gen.map.size[8:14]), median(tmp.num.markers[1:7]), median(tmp.num.markers[8:14]), median(tmp.phys.map.size)))
}

rec.rate.df <- as.data.frame(rec.rate.df, stringsAsFactors = FALSE)
colnames(rec.rate.df) <- c("CHROM", "PAR", "REC_RATE")
rec.rate.df[,"REC_RATE"] <- as.numeric(rec.rate.df[,"REC_RATE"])



gen.map.stat.df <- apply(gen.map.stat.df, 2, as.numeric)
gen.map.stat.df <- apply(gen.map.stat.df, c(1,2), function(x){round(x, digits = 1)})
gen.map.stat.df <-  as.data.frame(gen.map.stat.df, stringsAsFactors = FALSE)

boxplot(REC_RATE ~ CHROM, data = rec.rate.df, xlab = "Chromosome", ylab = "Average recombination rate (cM/Mb)")
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################