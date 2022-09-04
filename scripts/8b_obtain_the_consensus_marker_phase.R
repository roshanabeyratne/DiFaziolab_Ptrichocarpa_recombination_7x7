#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: Created: June 17, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name. 

# LD based haplotypes are obtained from the first four columns of the following file format; 
# CHROM,POS,hap1,hap2,24757,24758,24759,24760,24763,24765,24772,24773,24777,24778,24781,24782,24783,24784,24785,24793,24794,26705,25431,25436,25437,25438,25439,25443,25448,25452,25453,25454,25455,25460,25465,25470,25472,25476,25483,25322,25325,25332,25334,25336,25339,25503,25504,25505,25506,25507,25509,25510,25511,25512,25513,25514,26979,24802,24804,24805,24806,24809,24810,24811,24812,24817,24818,24823,24825,24828,24830,24831,24833,24834,26998,25515,25516,25517,25518,25519,25520,25522,25523,25526,25527,25533,25535,24837,24839,24840,24841,24842,24844,24847,24848,24849,24850,24851,24852,24853,24859,24864,24865,24867,24868,25530,25531,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492,25493,25495,25498,25500
# Chr15,  173833,0,1,NA,1,1,1,1,1,1,NA,1,NA,1,2,1,NA,1,1,1,1,NA,NA,2,2,1,1,1,1,2,1,2,1,2,1,1,1,NA,1,2,NA,2,1,2,NA,NA,2,NA,NA,2,NA,1,NA,2,1,2,2,NA,2,1,NA,1,2,1,1,1,2,1,1,2,1,2,1,2,1,2,2,2,1,NA,2,1,1,1,2,2,2,1,2,1,2,1,2,2,2,2,2,2,1,NA,2,1,1,NA,1,1,1,1,1,1,1,1,2,1,1,2,NA,2,2,2,2,2,1
# Chr15,  174711,0,1,2,1,1,1,1,1,1,NA,1,2,1,2,1,2,1,2,1,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,2,1,2,NA,2,1,2,NA,NA,2,NA,NA,NA,2,1,NA,2,1,2,2,2,2,1,2,1,2,1,1,1,2,1,1,2,1,2,1,NA,2,NA,2,2,1,NA,NA,NA,1,1,1,2,2,1,1,1,2,1,NA,2,2,2,2,NA,2,2,2,1,1,2,1,2,2,1,1,1,2,1,1,1,NA,2,NA,2,NA,2,2,NA,2
# Chr15,  178446,1,0,2,1,2,NA,NA,1,1,NA,NA,2,1,2,1,2,NA,2,1,1,2,2,2,2,1,NA,2,1,2,1,2,1,2,2,1,2,2,1,2,2,2,NA,2,2,2,2,2,2,1,2,1,2,2,1,2,2,2,2,1,2,1,2,1,1,1,2,1,1,1,1,2,1,2,2,1,2,2,1,2,1,2,1,1,2,2,2,1,2,1,2,2,2,2,2,2,1,2,2,2,2,1,1,2,1,2,NA,1,1,1,NA,1,2,NA,NA,NA,NA,NA,NA,NA,NA,NA,2

# Onemap based haplotypes are obtained from the following file format; 
# CHROM,POS,HAP1,HAP2
# Chr15,173833,0,1
# Chr15,174711,0,1
# Chr15,178446,1,0
# Chr15,189430,1,0
# Chr15,195830,0,1
# Chr15,211967,1,0
# Chr15,226300,0,1
# Chr15,242915,0,1
# Chr15,242997,0,1

# PURPOSE: Compare marker phasing using the two methods (onemap based and adjacent LD based), to obtain a consensus marker set and phases for the focal parent.
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
args <- commandArgs(TRUE)
chrom.name <- args[1]
par.name <- args[2]
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
input.path.LD.phase <- "./haplotypes/by_halfsib_family/"
input.path.onemap.phase <- "./onemap/phase_parents/phased_markers/"
output.path.consensus <- "./consensus_phase/"
#######################################################################################################################################################################################################

# for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr03"
  # for (par.name in par.name.vec) { # par.name <- "2515"
    
    # reading in the focal parental marker phase obtained using the adjacent-marker LD method 
    phase.LD <- read.csv(file = paste0(input.path.LD.phase, chrom.name, "_halfsibfamily_haplotypes_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    
    # reading in the focal parental marker phase obtained using all marker clustering method using onemap:R-package
    phase.onemap <- read.csv(file = paste0(input.path.onemap.phase, chrom.name, "_phased_markers_", par.name, ".csv"), stringsAsFactors = FALSE, header =TRUE)
    
    # merge the two data-sets uisng the CHROM and POS columns
    common.phase.tab <- merge(phase.LD[,1:4], phase.onemap, by = c("CHROM", "POS"))
    common.phase.tab$POS <- sapply(common.phase.tab$POS, as.numeric)
    ordrd.idxs <- order(common.phase.tab$POS, decreasing = FALSE)
    common.phase.tab.2 <- common.phase.tab[ordrd.idxs, ]
    
    # sum(phase.LD$POS %in% phase.onemap$POS)
    # sum(common.phase.tab.2$hap1 == common.phase.tab.2$HAP2)
    
    # incorporate markers that were excluded from the adjacent-marker LD method
    for (i in 1:nrow(phase.LD)) {
      pos <- phase.LD$POS[i]
      if (pos %in% common.phase.tab.2$POS) {
        next
      } else {
        common.phase.tab.2 <- rbind(common.phase.tab.2, c(chrom.name, pos, phase.LD$hap1[i], phase.LD$hap2[i], NA, NA))
      }
    }
    
    # incorporate markers that were excluded from the onemap based method
    for (i in 1:nrow(phase.onemap)) {
      pos <- phase.onemap$POS[i]
      if (pos %in% common.phase.tab.2$POS) {
        next
      } else {
        common.phase.tab.2 <- rbind(common.phase.tab.2, c(chrom.name, pos, NA, NA, phase.onemap$HAP1[i], phase.onemap$HAP2[i]))
      }
    }
    
    # which(is.na(common.phase.tab.2$HAP1))
    colnames(common.phase.tab.2) <- c(colnames(common.phase.tab.2[1:2]), "LD.hap.1", "LD.hap.2", "OM.hap.1", "OM.hap.2")
    ordered.idxs.2 <- order(common.phase.tab.2$POS, decreasing = FALSE)
    common.phase.tab.3 <- common.phase.tab.2[ordered.idxs.2, ]
    write.csv(common.phase.tab.3, file = paste0(output.path.consensus, chrom.name, "_consensus_phase_", par.name, ".csv"), row.names = FALSE, quote = FALSE)
    cat(paste0(output.path.consensus, chrom.name, "_consensus_phase_", par.name, ".csv"), " Successfully completed!", "\n")
  # }
# }
########################################################################################################################################################################################################
# END
# Following is the output file format compiles phased focal parent phases side by side; 
    # CHROM,POS,LD.hap.1,LD.hap.2,OM.hap.1,OM.hap.2
    # Chr16,10000291,0,1,0,1
    # Chr16,10000322,0,1,0,1
    # Chr16,10000348,0,1,0,1
    # Chr16,10000443,0,1,0,1
    # Chr16,10000520,0,1,0,1
    # Chr16,10000640,0,1,0,1
    # Chr16,10000665,0,1,0,1
    # Chr16,10000725,0,1,0,1
    # Chr16,10000735,0,1,0,1