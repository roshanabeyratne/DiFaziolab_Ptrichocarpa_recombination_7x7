#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: August 30, 2020
# USAGE: Input files include genetic maps. The marker names contain the chromosome and marker putative physical position. 

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

# PURPOSE: This is to identify the distribution of the distance between neighboring markers (as per Stet-14 physical distance), included in the genetic map. 
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
genmap.stats.output.path <- "./onemap/bin_markers/genetic_maps/stats/"
########################################################################################################################################################################################################
neighbor.dist.phys <- NULL
neighbor.dist.gen <- NULL


for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  for (par.name in par.name.vec) { # par.name <- "1863"
    
    if (file.exists(paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"))) {
      gen.map <- read.table(file = paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    } else {
      gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    }
    
    # gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    string <- paste0(chrom.name, "_")
    physical.distance <- gen.map[,2] %>% str_replace(string, "") %>% as.numeric()
    genetic.distance <- gen.map[,3] %>% as.numeric()
    
    # This is a test code block to check whether the marker order in the genetic map conforms to the physical order.
    if (all(as.numeric(physical.distance) == sort(as.numeric(physical.distance), decreasing = FALSE))) {
      neighbor.dist.phys <- c(neighbor.dist.phys, diff(physical.distance))
    } else {
      neighbor.dist.phys <- c(neighbor.dist.phys, diff(as.numeric(physical.distance), decreasing = FALSE))
      cat("error in order do ripple_seq", " ", chrom.name, " ", par.name, "\n")
    }
    
    neighbor.dist.gen <- c(neighbor.dist.gen, diff(genetic.distance))
  }
}

median(neighbor.dist.phys) # 126109
mean(neighbor.dist.phys) # 207220
median(neighbor.dist.gen) # 1.031314
mean(neighbor.dist.gen) # 1.283817

# jpeg(filename = paste0(genmap.stats.output.path, "/histogram_of_nearest_neighbor_distance.jpg"), width = 1000, height = 500, quality=1500)
pdf(file=paste0(genmap.stats.output.path, "/histogram_of_nearest_neighbor_distance.pdf"), width = 8, height = 5)
op <- par(mfrow=c(1,2))
hist(neighbor.dist.phys, breaks = 10000, xlab = "Physical distance to the nearest neighboring marker (bps)", ylab = "frequency", 
     main = "Histogram of nearest neighbor marker physical distance")
hist(neighbor.dist.gen, breaks = 1000, xlab = "Genetic distance to the nearest neighboring marker (bps)", ylab = "frequency", 
     main = "Histogram of nearest neighbor marker genetic-map distance")
par(op)
dev.off()

# output for ripple_seq adjustment
# error in order do ripple_seq   Chr01   2283 
# error in order do ripple_seq   Chr01   2365 
# error in order do ripple_seq   Chr01   2393 
# error in order do ripple_seq   Chr01   2515 
# error in order do ripple_seq   Chr01   6909 
# error in order do ripple_seq   Chr01   7073 
# error in order do ripple_seq   Chr02   2283 
# error in order do ripple_seq   Chr02   2683 
# error in order do ripple_seq   Chr02   7073 
# error in order do ripple_seq   Chr03   7073 
# error in order do ripple_seq   Chr04   1863 
# error in order do ripple_seq   Chr04   2283 
# error in order do ripple_seq   Chr05   1863 
# error in order do ripple_seq   Chr05   2283 
# error in order do ripple_seq   Chr05   2393 
# error in order do ripple_seq   Chr06   2283 
# error in order do ripple_seq   Chr06   2572 
# error in order do ripple_seq   Chr07   1863 
# error in order do ripple_seq   Chr07   2283 
# error in order do ripple_seq   Chr07   6909 
# error in order do ripple_seq   Chr09   2365 
# error in order do ripple_seq   Chr09   7073 
# error in order do ripple_seq   Chr10   2365 
# error in order do ripple_seq   Chr11   2048 
# error in order do ripple_seq   Chr11   2283 
# error in order do ripple_seq   Chr11   2365 
# error in order do ripple_seq   Chr11   7073 
# error in order do ripple_seq   Chr12   2393 
# error in order do ripple_seq   Chr12   7073 
# error in order do ripple_seq   Chr13   2048 
# error in order do ripple_seq   Chr14   2283 
# error in order do ripple_seq   Chr14   7073 
# error in order do ripple_seq   Chr15   2048 
# error in order do ripple_seq   Chr15   2283 
# error in order do ripple_seq   Chr15   2365 
# error in order do ripple_seq   Chr16   2048 
# error in order do ripple_seq   Chr16   2066 
# error in order do ripple_seq   Chr16   6909 
# error in order do ripple_seq   Chr17   2283 
# error in order do ripple_seq   Chr17   2365 
# error in order do ripple_seq   Chr17   2393 
# error in order do ripple_seq   Chr17   2683 
# error in order do ripple_seq   Chr17   6909 
# error in order do ripple_seq   Chr17   7073 
# error in order do ripple_seq   Chr18   4593 
# error in order do ripple_seq   Chr19   1863 
# error in order do ripple_seq   Chr19   2283 
# error in order do ripple_seq   Chr19   2572

# Note of September 10, 2020:
# After correcting above maps by ripple_seq function as implementded in S_7b script. Following maps still has markers that do not conform to the physical order. 
# error in order do ripple_seq   Chr17   2393 
# error in order do ripple_seq   Chr17   6909 
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################