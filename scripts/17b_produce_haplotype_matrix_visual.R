#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: July 23, 2020
# USAGE: This script is run nteractively on spruce. Input files include smothed-out bakground haplotype contribution from the focal parent

# Following is the input file format used to produce the pdfs; 
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1

# PURPOSE: Produce a visual representation of focal parent haplotypes in all offspring.
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
input.path <- "./resolved_haps/by_parent/tmp/"
input.path.corrected <- "./resolved_haps/by_parent/"
input.path.OM.based <- "./haplotypes/OM_based/"
output.path.jpgs <- "./resolved_haps/by_parent/tmp/jpgs/"
output.path.corrected.jpgs <- "./resolved_haps/by_parent/jpgs/"
########################################################################################################################################################################################################
# capabilities() # just to check whether R has X11 capabilities to run image command
# This loop is for the non-corrected portion of the output of script 16b_...
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[2]
  for (par.name in par.name.vec) { # par.name <- "1950"
    
    # read-in the phased marker haplotypes for offspring
    phased.offspring.tab <- read.csv(file = paste0(input.path.OM.based, chrom.name, "_phased_offspring_OM_",  par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    par.haps.tab <- read.csv(paste0(input.path, chrom.name, "_tmp_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
    
    # checking if any of the offspring have NAs
    tmp.logic <- apply(par.haps.tab[,-(1:2)], 2, function(x){
      tmp <- any(is.na(x))
      return(tmp)
    })
    
    if (any(tmp.logic)) {
      stop("offsprings contain NA values at markers!")
    } else {
      
      # obtain only the overlapping markers
      par.haps.tab <- par.haps.tab[par.haps.tab[,2] %in% phased.offspring.tab[,2], ]
      temp.tab <- as.matrix(par.haps.tab[, -(1:2)])
      # unique(as.vector(temp.tab))
      
      jpeg(filename = paste0(output.path.jpgs, chrom.name, "_tmp_resolved_background_haplotypes_all_offspring_", par.name,".jpg"))
      # pdf(file = paste0(output.path.pdfs, chrom.name, "_tmp_resolved_background_haplotypes_all_offspring_", par.name,".pdf"))
      image(temp.tab, xlab = "markers", ylab = "offspring", col=c("yellow", "blue", "red"))
      dev.off()
      cat(paste0(output.path.jpgs, chrom.name, "_tmp_resolved_background_haplotypes_all_offspring_", par.name,".jpg", " Successfully Completed!", "\n"))
    }
  }
}
########################################################################################################################################################################################################

# This loop is for the corrected portion of the output of script 16b_...
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[2]
  for (par.name in par.name.vec) { # par.name <- "1950"
    
    # read-in the phased marker haplotypes for offspring
    phased.offspring.tab <- read.csv(file = paste0(input.path.OM.based, chrom.name, "_phased_offspring_OM_",  par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    par.haps.tab <- read.csv(paste0(input.path.corrected, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
    
    # checking if any of the offspring have NAs
    tmp.logic <- apply(par.haps.tab[,-(1:2)], 2, function(x){
      tmp <- any(is.na(x))
      return(tmp)
    })
    
    if (any(tmp.logic)) {
      stop("offsprings contain NA values at markers!")
    } else {
      
      # obtain only the overlapping markers
      par.haps.tab <- par.haps.tab[par.haps.tab[,2] %in% phased.offspring.tab[,2], ]
      temp.tab <- as.matrix(par.haps.tab[, -(1:2)])
      # unique(as.vector(temp.tab))
      
      
      jpeg(filename = paste0(output.path.corrected.jpgs, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".jpg"))
      # pdf(file = paste0(output.path.pdfs, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".pdf"))
      image(temp.tab, xlab = "markers", ylab = "offspring", col=c("yellow", "blue", "red"))
      dev.off()
      cat(paste0(output.path.corrected.jpgs, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".jpg", " Successfully Completed!", "\n"))
    }
  }
}
########################################################################################################################################################################################################
# END