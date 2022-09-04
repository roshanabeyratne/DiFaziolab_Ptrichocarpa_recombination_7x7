#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: Sep 6, 2020
# USAGE: Input files include genetic maps. The marker names contain the chromosome and marker putative physical position. 
# # Following is the input file / genetic map format; 
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

# This is the list of genetic maps which require manual adjustment of markers based on the ripple_seq output 
## output for ripple_seq adjustment
# error in order do ripple_seq   Chr01   2365 
# error in order do ripple_seq   Chr01   2393 # done
# error in order do ripple_seq   Chr01   6909 # done
# error in order do ripple_seq   Chr02   2048 # done
# error in order do ripple_seq   Chr02   2283 # done
# error in order do ripple_seq   Chr02   2683 # done
# error in order do ripple_seq   Chr02   7073 # done 
# error in order do ripple_seq   Chr03   7073 # done
# error in order do ripple_seq   Chr04   1863 # done 
# error in order do ripple_seq   Chr04   2283 # done
# error in order do ripple_seq   Chr04   7073 # done 
# error in order do ripple_seq   Chr05   2283 # done 
# error in order do ripple_seq   Chr06   1863 # done 
# error in order do ripple_seq   Chr06   2283 # done 
# error in order do ripple_seq   Chr06   2572 # done
# error in order do ripple_seq   Chr07   1863 # done
# error in order do ripple_seq   Chr07   4593 # done 
# error in order do ripple_seq   Chr07   6909 # done 
# error in order do ripple_seq   Chr09   2365 # done 
# error in order do ripple_seq   Chr09   7073 # done 
# error in order do ripple_seq   Chr10   2365 # done 
# error in order do ripple_seq   Chr10   2572 # done 
# error in order do ripple_seq   Chr11   2048 # done 
# error in order do ripple_seq   Chr11   2283 # done 
# error in order do ripple_seq   Chr11   7073 # done 
# error in order do ripple_seq   Chr12   7073 # done 
# error in order do ripple_seq   Chr13   2393 # done 
# error in order do ripple_seq   Chr13   2572 # done 
# error in order do ripple_seq   Chr14   2283 # done
# error in order do ripple_seq   Chr15   2283 # done 
# error in order do ripple_seq   Chr16   2048 # done 
# error in order do ripple_seq   Chr16   2066 # done 
# error in order do ripple_seq   Chr16   2515 # done 
# error in order do ripple_seq   Chr16   6909 # done 
# error in order do ripple_seq   Chr17   1909 # done
# error in order do ripple_seq   Chr17   2365 # done
# error in order do ripple_seq   Chr17   2393 # not-changed
# error in order do ripple_seq   Chr17   2683 # done
# error in order do ripple_seq   Chr17   6909 # not-changed
# error in order do ripple_seq   Chr17   7073 # done 
# error in order do ripple_seq   Chr18   7073 # done
# error in order do ripple_seq   Chr19   1909 # done
# error in order do ripple_seq   Chr19   2572 # done

# PURPOSE: Rectify incorrectly ordered markers (if a marker order does not conform to the physical order and has an alternative marker order than what is initially suggested by order_seq funciton, 
# it is defined as placed incorrectly) using ripple_seq function in Onemap R-package.
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
library(onemap)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

rippled.genetic.map.output.path <- "./onemap/bin_markers/genetic_maps/rippled_maps/"
########################################################################################################################################################################################################

for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  for (par.name in par.name.vec) { # par.name <- "2365"
    load(paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/binned_db_raw/", chrom.name, "_binned_db_raw_object_", par.name,".RData"))
    load(paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/twopts/", chrom.name, "_twopts_object_", par.name,".RData"))
    load(paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/LG1_rec/", chrom.name, "_full_genetic_map_object_", par.name,".RData"))
    
    LG.rec.map <- make_seq(LG1.rec, "safe")
    diff(LG1.rec$ord$seq.num)
    ripple_seq(LG.rec.map, ws = 4, LOD = 3, tol = 10E-5)
    
    temp.map <- make_seq(twopts, sort(LG1.rec$ord$seq.num, decreasing = FALSE))
    temp.map2 <- map(temp.map)
    LG1.rec$ord
    
    diff(temp.map2$seq.num)
    write_map(temp.map2, file.out = paste0(rippled.genetic.map.output.path, chrom.name, "_", par.name, "_genetic.map"))
    
  }
}
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################