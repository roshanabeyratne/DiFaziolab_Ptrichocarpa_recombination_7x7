#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: March 23, 2020; Modified - May 30, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name. Input files include haplotype contribution for the focal parent for a fulsib family.

# This is the intput format for the focal parent haplotypes
# CHROM,POS,REF,ALT,1863,2365,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553
# Chr15,8393,C,T,1,1,NA,NA,NA,NA,NA,1,NA,0,0,NA,0,NA
# Chr15,8421,C,T,1,1,NA,NA,NA,NA,NA,1,NA,0,0,0,0,NA
# Chr15,8426,T,A,1,1,NA,NA,NA,NA,NA,1,0,0,0,NA,0,NA
# Chr15,8805,G,A,1,1,NA,0,NA,NA,NA,1,NA,0,0,0,0,NA
# Chr15,8951,A,T,0,1,0,0,0,0,0,NA,0,0,0,0,0,NA
# Chr15,9019,C,T,1,1,NA,0,NA,NA,NA,1,NA,NA,0,NA,0,NA
# Chr15,9178,T,C,1,1,NA,0,NA,0,0,1,NA,NA,0,NA,0,NA
# Chr15,10452,C,T,1,1,NA,NA,NA,NA,NA,1,NA,NA,0,0,0,NA
# Chr15,10571,A,G,1,1,NA,0,0,NA,NA,1,NA,NA,0,NA,0,NA

# PURPOSE: This followup to the script number 4b where each focal parent's alllele contribution in a full-sib family are combined across all the offspring for a given focal parent across all 
# the families which it is a part of. This is the template for the onemap input. 
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
input.path.haplotypes <- "./haplotypes/by_fulsib_family/"
output.path.combined.haplo <- "./haplotypes/by_parent/"
########################################################################################################################################################################################################
# for (chrom.name in chrom.name.vec[1:19]) { # chrom.name <- chrom.name.vec[1]
  # for (par.name in par.name.vec[1:14]) { # par.name <- "2283"
    
    if (any(par.name %in% par.name.vec[1:7])) {
      other.par.name.vec <- par.name.vec[8:14]
      focal.column <- 5
      common.tab <- read.csv(file = paste0(input.path.haplotypes, chrom.name, "_fulsib_haplotypes_", par.name, "X2365_2365.csv"), stringsAsFactors = FALSE, header = TRUE)
      common.tab <- common.tab[, c(1:4, focal.column)]
      
      for (other.par.name in other.par.name.vec) {
        # other.par.name <- "2572"
        temp.df <- read.csv(file = paste0(input.path.haplotypes, chrom.name, "_fulsib_haplotypes_", par.name, "X", other.par.name, "_", par.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)
        temp.df2 <- as.data.frame(temp.df[, -(1:6)], stringsAsFactors = FALSE)
        colnames(temp.df2) <- colnames(temp.df)[-(1:6)] # Counters issues for fulsib famillies with one offspring
        common.tab <- cbind(common.tab, temp.df2)
      }
      
    } else if (any(par.name %in% par.name.vec[8:14])) {
      other.par.name.vec <- par.name.vec[1:7]
      focal.column <- 6
      common.tab <- read.csv(file = paste0(input.path.haplotypes, chrom.name, "_fulsib_haplotypes_1863X", par.name, "_1863.csv"), stringsAsFactors = FALSE, header = TRUE)
      common.tab <- common.tab[, c(1:4, focal.column)]
      
      for (other.par.name in other.par.name.vec) {
        temp.df <- read.csv(file = paste0(input.path.haplotypes, chrom.name, "_fulsib_haplotypes_", other.par.name, "X", par.name, "_", par.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)
        temp.df2 <- as.data.frame(temp.df[, -(1:6)], stringsAsFactors = FALSE)
        colnames(temp.df2) <- colnames(temp.df)[-(1:6)] # Counters issues for fulsib famillies with one offspring
        common.tab <- cbind(common.tab, temp.df2)
      }
    }
    
    column.name.order.1 <- c(colnames(common.tab)[1:4], paste0("X", par.name))
    column.name.order.2 <- colnames(common.tab)[!(colnames(common.tab) %in% column.name.order.1)]
    
    if ((length(c(column.name.order.1, column.name.order.2))) == ncol(common.tab)) {
      
      parent.haplo.tab <- common.tab[, c(column.name.order.1, column.name.order.2)]
      colnames(parent.haplo.tab) <- colnames(parent.haplo.tab) %>% str_replace(".GT", "") %>% str_replace("^X", "")
      write.csv(parent.haplo.tab, file = paste0(output.path.combined.haplo, chrom.name, "_haplotypes_", par.name, ".csv"), quote = FALSE, row.names = FALSE)
      cat(chrom.name, "\t", par.name, "\t", "Successfully completed!", "\n")
    } else {
      stop("error with number of columns")
    }
  # }
# }
########################################################################################################################################################################################################
#END
    # This is the output format for focal parent haplotype contribution in all its offspring.
    # CHROM,POS,REF,ALT,1863,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553,26548,26549,26551,26552,26554,26555,26556,26557,24716,24717,24718,24719,24720,24721,24722,24723,24724,24725,24726,24727,24728,24729,24730,24731,24732,24734,26564,26565,26566,26567,26570,26572,26573,26574,26575,26577,26578,26581,26584,26587,26589,26591,26592,26593,25998,26598,26600,26602,26603,26605,26606,26607,26608,26610,26611,26612,26613,26615,26616,26620,26623,26624,26625,26798,24736,24740,24745,24747,26631,26635,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677
    # Chr15,8393,C,T,1,NA,NA,NA,NA,NA,1,NA,0,0,NA,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,1,NA,0,1,1,0,1,1,NA,0,0,1,0,NA,0,1,NA,NA,0,NA,NA,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,0,0,1,1,0,0,NA,NA,0,0,0,1,0,NA,NA,NA,NA,NA,0,1,1,1,NA,0,0
    # Chr15,8421,C,T,1,NA,NA,NA,NA,NA,1,NA,0,0,0,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,NA,NA,0,1,0,NA,0,1,NA,NA,0,0,1,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,1,0,1,1,0,0,1,NA,0,0,0,1,0,NA,NA,NA,NA,NA,0,1,1,NA,NA,NA,0
    # Chr15,8426,T,A,1,NA,NA,NA,NA,NA,1,0,0,0,NA,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,NA,NA,0,1,0,NA,0,1,NA,NA,0,NA,1,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,1,0,1,1,0,0,1,NA,0,0,0,1,0,NA,NA,NA,NA,NA,NA,1,1,NA,NA,0,0
    # Chr15,8805,G,A,1,NA,0,NA,NA,NA,1,NA,0,0,0,0,NA,0,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,0,NA,0,NA,0,0,0,1,NA,0,0,0,1,NA,NA,0,1,1,1,1,0,1,1,0,NA,1,1,1,0,0,1,NA,0,0,0,1,0,1,1,0,0,1,0,0,0,0,1,NA,NA,NA,0,NA,NA,NA,1,1,NA,NA,NA,0