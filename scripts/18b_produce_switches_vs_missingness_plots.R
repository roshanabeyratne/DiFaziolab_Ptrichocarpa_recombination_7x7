#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: July 23, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name. Input files include smothed-out bakground haplotype contribution from the focal parent

# Following is the input file format used to produce the pdfs; 
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1

# Following is the input format with haplotype contribution for each offspring based on a trio-analysis as produced by script 9b_.... contains entries 0,1 or 2 ('0' assigned when haplotypes not certain)
# CHROM,POS,LD.hap.1,LD.hap.2,OM.hap.1,OM.hap.2,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553,24757,24758,24759,24760,24763,24765,24772,24773,24777,24778,24781,24782,24783
# Chr02,6617,0,1,0,1,1,2,2,2,1,1,2,2,1,2,2,1,2,1,2,2,1,1,1,2,2,1,1,1,1,2,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,0,1,1,1,0,2,1,2,2,2,1,2,1,1,2,1,2,2,1,2,2,2,2,1,1,1,1,2,2,1,1,1,2,2,2,2,1,1,2,1,1,2,1,1,2,2,1,2


# PURPOSE: Produce a visual representation of focal parent haplotypes in all offspring.
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
input.path <- "./resolved_haps/by_parent/tmp/"
input.path.OM.based <- "./haplotypes/OM_based/"
output.path <- "./haplotypes/switch_miss/"
########################################################################################################################################################################################################
# capabilities()
for (chrom.name in chrom.name.vec) {# chrom.name <- "Chr02"
  for (par.name in par.name.vec) {# par.name <- "1909"
    if (file.exists(paste0(input.path, chrom.name, "_tmp_resolved_background_haplotypes_all_offspring_", par.name,".csv"))) {
      par.haps.tab <- read.csv(paste0(input.path, chrom.name, "_tmp_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
      colnames(par.haps.tab) <- colnames(par.haps.tab) %>% str_replace("^X", "")
      
      phased.offspring.tab <- read.csv(file = paste0(input.path.OM.based, chrom.name, "_phased_offspring_OM_",  par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
      colnames(phased.offspring.tab) <- colnames(phased.offspring.tab) %>% str_replace("^X", "")
      
      # subset the dataframe with parental haplotypes to for which data is available
      par.haps.tab <- par.haps.tab[par.haps.tab[,2] %in% phased.offspring.tab[,2], ]
      phased.offspring.tab <- phased.offspring.tab[phased.offspring.tab[,2] %in% par.haps.tab[,2], ]
      
      
      tmp.switch.tab <- par.haps.tab[, -(1:2)]
      tmp.miss.tab <- phased.offspring.tab[, -(1:6)]
      
      # if (all(dim(tmp.switch.tab) == dim(tmp.miss.tab))) {
      #   for (off.col in 1:ncol(tmp.switch.tab)) { # off.col <- 1
      #     bg.hap <- tmp.switch.tab[,off.col]
      #     
      #     # Imputing zero-islands that are flannked by the same haplotype.
      #     splitAt <- function(x, pos) {unname(split(x, cumsum(seq_along(x) %in% pos)))} # function splits avector of marker haplotype assignments to contiguous stretches of a given haplotype.
      #     max.dist.1s <- 1
      #     inconclusive.idx <- which(bg.hap == 0)
      #     idx.diff.vec <- diff(inconclusive.idx)
      #     diff.idx <- which(idx.diff.vec > max.dist.1s)
      #     zero.island.list <- splitAt(inconclusive.idx, (diff.idx + 1))
      #     
      #     
      #     if (length(zero.island.list) > 0) { # check whether undeiceded haplotype blocks exist
      #       for (iter in 1:length(zero.island.list)) { # iter <- 1
      #         mrkr.pos <- zero.island.list[[iter]]
      #         left.pos <- mrkr.pos[1]
      #         right.pos <- mrkr.pos[length(mrkr.pos)]
      #         
      #         # if the zero.island.list is at the begining of the chromosome then haplotype to the immediately downstream used to impute
      #         if (left.pos == 1) {
      #           bg.hap[mrkr.pos] <- rep(bg.hap[right.pos + 1], length(mrkr.pos))
      #           
      #           # if the zero.island.list is at the end of the chromosome then haplotype to the immediately upstream used to impute
      #         } else if (right.pos == length(bg.hap)) {
      #           bg.hap[mrkr.pos] <- rep(bg.hap[left.pos - 1], length(mrkr.pos))
      #           
      #           # if the zero.island.list is flanked by the same haplotype the zeros are imputed with the flanking haplotype
      #         } else if ((bg.hap[(left.pos - 1)] == bg.hap[(right.pos + 1)])) {
      #           bg.hap[mrkr.pos] <- rep(bg.hap[left.pos - 1], length(mrkr.pos))
      #           
      #         } else {
      #           next
      #         }
      #       }
      #     } else {
      #       bg.hap <- bg.hap
      #     }
      #     
      #     # assign the corrected vector
      #     # plot(bg.hap, pch = 16, cex = 0.2, col = "red")
      #     tmp.switch.tab[,off.col] <- bg.hap
      #   }
      # } else {
      #   stop("errror with the two main dataframes!")
      # }
      
      # Count the number of haplotype switches by means of counting number of blocks
      switch.vec <- apply(tmp.switch.tab, 2, function(x){
        tmp <- diff(x)
        tmp.diff <- sum(tmp > 0)
        return(tmp.diff)
      })
      
      miss.vec <- apply(tmp.miss.tab, 2, function(x){
        tmp <- sum(x == 0, na.rm = TRUE)
        return(tmp)
      })
      
      
      offspring.names <- colnames(tmp.switch.tab) %>% str_replace("^X", "")
      
      jpeg(filename = paste0(output.path, chrom.name, "_switch_miss_offspring_", par.name,".jpg"), width = 1000, height = 1000)
      if (max(switch.vec) < 10) {
        plot(miss.vec, switch.vec, pch = 16, cex = 1, col = 'grey', ylim = c(0, 10), main = par.name)
      } else {
        plot(miss.vec, switch.vec, pch = 16, cex = 1, col = 'grey', main = par.name)
      }
      text(miss.vec, switch.vec, labels = offspring.names, cex = 1, offset = 0.5, col = 'black')
      dev.off()
      cat(paste0(output.path, chrom.name, "_switch_miss_offspring_", par.name,".jpg", " Successfully Completed!", "\n"))
    } else {
      stop("cannot find file containing background haploypes!")
    }
  }
}
########################################################################################################################################################################################################
# END