#!/usr/bin/R
setwd('/gpfs/group/difazio/populus/gatk-7x7-Stet14')
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: Dec 2, 2021
# USAGE: Input files include smoothed-out bakground haplotype contribution from the focal parent. Some individuals were removed from relevant parents due to high missingness and switches with script
# 19b. Therefore the path to such curated files is "./resolved_haps/by_parent/curated/". If a parent is not corrected such as described before, the focal parent hap contribution is at
# "./resolved_haps/by_parent/"

# Following is the input file format used that has focal parent haplotype identification (hap_1 or hap_2) for all of its offspring; This version of the background haplotype contains entries for which
# minimum haplotype block has been limited to 500Kb. (refer end of script 16b_ for details).
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1

# PURPOSE: Output focal parental haplotype contribution in offspring for all half-sib families
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

hap.input.path <- "./resolved_haps/by_parent/"
curated.hap.input.path <- "./resolved_haps/by_parent/curated/"
qtl_haplotype_path <- "./haplotypes_for_qtl_analysis/"
########################################################################################################################################################################################################
CO.df <- NULL
co.count.df <- NULL
std.co.count.df <- NULL

for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  mean.CO.count.per.chrom <- NULL
  std.CO.count.per.chrom <- NULL
  for (par.name in par.name.vec) { # par.name <- "1863"
    
    # Read-in background haplotype tables for each focal parent per chromosome. If a curated file with certain offspring removed exists then I read from that location.
    if (file.exists(paste0(curated.hap.input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"))) {
      haplo.tab <- read.csv(file = paste0(curated.hap.input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
    } else {
      haplo.tab <- read.csv(file = paste0(hap.input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
    }
    colnames(haplo.tab) <- colnames(haplo.tab) %>% str_replace("^X", "")
    
    write.csv(haplo.tab, file = paste0(qtl_haplotype_path, chrom.name, "_parental_haplotypes_for_halfsib_fam_", par.name, ".csv"), quote = FALSE, row.names = FALSE)
    cat(chrom.name, ": ", par.name, "\n")
  }
}
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################