#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: July 29, 2020
# USAGE: This script is run nteractively on spruce. Input files include smothed-out bakground haplotype contribution from the focal parent

# Following is the input file format used that has focal parent haplotype identification (hap_1 or hap_2) for all of its offspring;  
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1

# PURPOSE: Remove the offspring identified as having high missingness and higher number of focal parent haplotype switches. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
# args <- commandArgs(TRUE)
# chrom.name <- args[1]
# par.name <- args[2]
########################################################################################################################################################################################################
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
input.path <- "./resolved_haps/by_parent/"
output.path <- "./resolved_haps/by_parent/curated/"
########################################################################################################################################################################################################
# following script carries out manual curation of the file and saves it in a different location. This is done using a script to keep a record. 
for (chrom.name in chrom.name.vec) {
  par.name <- "1909"
  # if (chrom.name == "Chr13") {
  #   offspring.names <- c("25334", "24859", "25526")
  # } else {
    offspring.names <- c("25334", "24859")
  # }
  par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
  colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
  
  included.cols <- !colnames(par.hap.tab) %in% offspring.names
  par.hap.tab <- par.hap.tab[, included.cols]
  write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
  cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
}
#####################################
for (chrom.name in chrom.name.vec) {
  par.name <- "1950"
  offspring.names <- c("24892", "24954", "24957")
  par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
  colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
  
  included.cols <- !colnames(par.hap.tab) %in% offspring.names
  par.hap.tab <- par.hap.tab[, included.cols]
  write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
  cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
}
#####################################
for (chrom.name in chrom.name.vec) {
  par.name <- "2048"
  offspring.names <- c("26410")
  par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
  colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
  
  included.cols <- !colnames(par.hap.tab) %in% offspring.names
  par.hap.tab <- par.hap.tab[, included.cols]
  write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
  cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
}
#####################################
for (chrom.name in chrom.name.vec) {
  par.name <- "2515"
  if (chrom.name == "Chr17") {
    offspring.names <- c("25334", "24892", "26564")
  } else {
    offspring.names <- c("25334", "24892")
  }
  par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
  colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
  
  included.cols <- !colnames(par.hap.tab) %in% offspring.names
  par.hap.tab <- par.hap.tab[, included.cols]
  write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
  cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
}
#####################################
for (chrom.name in chrom.name.vec) {
  par.name <- "7073"
  offspring.names <- c("24954", "24957")
  par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
  colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
  
  included.cols <- !colnames(par.hap.tab) %in% offspring.names
  par.hap.tab <- par.hap.tab[, included.cols]
  write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
  cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
}
#####################################
for (chrom.name in chrom.name.vec) {
  par.name <- "6909"
  if (chrom.name == "Chr13") {
    offspring.names <- c("24859", "25653") 
  } else if (chrom.name == "Chr17") {
    offspring.names <- c("24859", "26265") 
  } else {
    offspring.names <- c("24859")
  }
  par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
  colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
  
  included.cols <- !colnames(par.hap.tab) %in% offspring.names
  par.hap.tab <- par.hap.tab[, included.cols]
  write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
  cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
}
#####################################
for (chrom.name in chrom.name.vec) {
  par.name <- "2683"
  offspring.names <- c("26410")
  par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
  colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
  
  included.cols <- !colnames(par.hap.tab) %in% offspring.names
  par.hap.tab <- par.hap.tab[, included.cols]
  write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
  cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
}
########################################################################################################################################################################################################
# Following are the parents and offspring for which specific individuals need to be removed for a specific chromosome owing to high missingness leading to high switches.
chrom.name <- "Chr14"
par.name <- "2365"
offspring.names <- c("26705")
par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
included.cols <- !colnames(par.hap.tab) %in% offspring.names
par.hap.tab <- par.hap.tab[, included.cols]
write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
#####################################
chrom.name <- "Chr15"
par.name <- "1863"
offspring.names <- c("24715")
par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
included.cols <- !colnames(par.hap.tab) %in% offspring.names
par.hap.tab <- par.hap.tab[, included.cols]
write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
#####################################
chrom.name <- "Chr15"
par.name <- "2572"
offspring.names <- c("26151")
par.hap.tab <- read.csv(file = paste0(input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
colnames(par.hap.tab) <- colnames(par.hap.tab) %>% str_replace("^X", "")
included.cols <- !colnames(par.hap.tab) %in% offspring.names
par.hap.tab <- par.hap.tab[, included.cols]
write.csv(par.hap.tab, file = paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
cat(paste0(output.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
########################################################################################################################################################################################################
# END