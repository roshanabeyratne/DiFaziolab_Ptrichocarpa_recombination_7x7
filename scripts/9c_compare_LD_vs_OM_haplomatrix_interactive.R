#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: Created: July 25, 2020
# USAGE: This script needs to be run interactively on spruce as image() function in base-R does not work as a batch script on spruce.
#
# Following is the input format for this script which is the output fof script 9b. (key: 0-missing, 1-hap_1, 2-hap_2)
# CHROM,POS,LD.hap.1,LD.hap.2,OM.hap.1,OM.hap.2,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553,24757,24758,24759,24760,24763,24765,24772,24773,24777,24778,24781,24782,24783
# Chr02,6617,0,1,0,1,1,2,2,2,1,1,2,2,1,2,2,1,2,1,2,2,1,1,1,2,2,1,1,1,1,2,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,0,1,1,1,0,2,1,2,2,2,1,2,1,1,2,1,2,2,1,2,2,2,2,1,1,1,1,2,2,1,1,1,2,2,2,2,1,1,2,1,1,2,1,1,2,2,1,2

# PURPOSE: This script takes the output of the script 9b for both LD and OM methods and makes focal parent haplotype matrix for all of its offspring.
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

input.path.LD.based <- "./haplotypes/LD_based/"
input.path.OM.based <- "./haplotypes/OM_based/"

output.path.jpegs <- "./consensus_phase/jpegs/"
########################################################################################################################################################################################################

for (chrom.name in chrom.name.vec) {# chrom.name <- "Chr01"
  for (par.name in par.name.vec) {# par.name <- "1909"
    for (pipline in c("LD", "OM")) { # pipline <- "LD"
      if (pipline == "LD") {
        phased.offspring.path <- input.path.LD.based
      } else if (pipline == "OM") {
        phased.offspring.path <- input.path.OM.based
      } else {
        stop("pipiline not identified for outputpath!")
      }
      
      if (file.exists(paste0(phased.offspring.path, chrom.name, "_phased_offspring_", pipline, "_",  par.name, ".csv"))) {
        phased.offspring.tab <- read.csv(file = paste0(phased.offspring.path, chrom.name, "_phased_offspring_", pipline, "_",  par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
        temp.tab <- as.matrix(phased.offspring.tab[, -(1:6)])
        jpeg(filename = paste0(output.path.jpegs, chrom.name, "_haplotypes_", pipline, "_", par.name, ".jpg"))
        image(temp.tab, xlab = "markers", ylab = "offspring")
        dev.off()
        cat(paste0(output.path.jpegs, chrom.name, "_haplotypes_", pipline, "_", par.name, ".jpg"), " Successfully completed!", "\n")
      } else {
        stop("file does not exist")
      } 
    }
  }
}
########################################################################################################################################################################################################
# END