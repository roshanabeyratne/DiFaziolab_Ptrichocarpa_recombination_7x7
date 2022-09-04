#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: July 22, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name. Input files include smothed-out bakground haplotype contribution from the focal parent and 
# another input file that identifies all the offspring for the given focal parent.

# This is the main intput with actual and smoothed-background haplotypes produced for a given offspring with script 15b; 
# CHROM   POS REF ALT 2048 2515 hap1_2048 hap2_2048 hap1_2515 hap2_2515 25663 2048_hap_combo 2048_background 2515_hap_combo 2515_background
# Chr02  6617   T   C    0    2         0         0         1         1     1              0               0              0               0
# Chr02  8424   A   G    0    2         0         0         1         1     1              0               0              0               0
# Chr02 12411   C   G    0    2         0         0         1         1     1              0               0              0               0
# Chr02 13186   T   G    0    2         0         0         1         1     1              0               0              0               0
# Chr02 14242   G   A    0    2         0         0         1         1     1              0               0              0               0
# Chr02 14663   G   A    0    2         0         0         1         1     1              0               0              0               0

# The other input format identifies all the offspring for a given focal parent; 
# This is the output format of script 5b_ for focal parent haplotype contribution in all its offspring.
# CHROM,POS,REF,ALT,1863,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553,26548,26549,26551,26552,26554,26555,26556,26557,24716,24717,24718,24719,24720,24721,24722,24723,24724
# Chr15,8393,C,T,1,NA,NA,NA,NA,NA,1,NA,0,0,NA,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,1,NA,0,1,1,0,1,1,NA,0,0,1,0,NA,0,1,NA,NA,0,NA,NA,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,0,0, 1
# Chr15,8421,C,T,1,NA,NA,NA,NA,NA,1,NA,0,0,0,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,NA,NA,0,1,0,NA,0,1,NA,NA,0,0,1,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,1,0,1,1
# Chr15,8426,T,A,1,NA,NA,NA,NA,NA,1,0,0,0,NA,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,NA,NA,0,1,0,NA,0,1,NA,NA,0,NA,1,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,1,0,1,1
# Chr15,8805,G,A,1,NA,0,NA,NA,NA,1,NA,0,0,0,0,NA,0,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,0,NA,0,NA,0,0,0,1,NA,0,0,0,1,NA,NA,0,1,1,1,1,0,1,1,0,NA,1,1,1,0,0,1,NA,0,0,0,1,0,1,1,0,0,1

# PURPOSE: This combines all offspring background-hapotype data for a focal parent into one data table. Also imputes missing marker data for a given offspring provided the missing marker is flanked by 
# the same haplotype. If not the markers are coded as containing missing information for the offspring.  
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
input.path.all_offpring <- "./haplotypes/by_parent/"
input.path.offspring.background <- "./resolved_haps/"
output.path.haplo.background <- "./resolved_haps/by_parent/"
tmp.output.path.haplo.background <- "./resolved_haps/by_parent/tmp/"
########################################################################################################################################################################################################
# for (chrom.name in chrom.name.vec[1:19]) { # chrom.name <- chrom.name.vec[2]
#   for (par.name in par.name.vec[1:14]) { # par.name <- "2365"  
    temp.tab <- read.csv(file = paste0(input.path.all_offpring, chrom.name, "_haplotypes_", par.name, ".csv"), header  = TRUE, stringsAsFactors = FALSE)
    colnames(temp.tab) <- colnames(temp.tab) %>% str_replace("^X", "")
    all.off.names <- colnames(temp.tab)[-(1:5)]
    
    
    # following code block merges the hapotype information from the focal parent into one table based on CHROM and POS information. 
    for (offspring.name in all.off.names) { # offspring.name <- "25208"
      if (offspring.name == all.off.names[1]) {
        # tmp.1 <- read.csv(file = paste0(input.path.offspring.background, chrom.name, "_resolved_haplotypes_both_parents_", offspring.name,".csv"), header  = TRUE, stringsAsFactors = FALSE)
        # tmp.2 <- read.csv(file = paste0(input.path.offspring.background, chrom.name, "_resolved_haplotypes_both_parents_", offspring.name,".csv"), header  = TRUE, stringsAsFactors = FALSE)
        
        temp.off.info.tab <- read.csv(file = paste0(input.path.offspring.background, chrom.name, "_resolved_haplotypes_both_parents_", offspring.name,".csv"), header  = TRUE, stringsAsFactors = FALSE)
        colnames(temp.off.info.tab) <- colnames(temp.off.info.tab) %>% str_replace("^X", "")
        off.info.tab <- cbind(temp.off.info.tab[1:2], temp.off.info.tab[,paste0(par.name, "_background")])
        colnames(off.info.tab)[ncol(off.info.tab)] <- offspring.name
        
      } else {
        temp.off.info.tab <- read.csv(file = paste0(input.path.offspring.background, chrom.name, "_resolved_haplotypes_both_parents_", offspring.name,".csv"), header  = TRUE, stringsAsFactors = FALSE)
        colnames(temp.off.info.tab) <- colnames(temp.off.info.tab) %>% str_replace("^X", "")
        temp.off.info.tab <- cbind(temp.off.info.tab[1:2], temp.off.info.tab[,paste0(par.name, "_background")])
        colnames(temp.off.info.tab)[ncol(temp.off.info.tab)] <- offspring.name
        
        off.info.tab <- merge(off.info.tab, temp.off.info.tab, all = TRUE, by = c("CHROM", "POS"))
        colnames(off.info.tab)[ncol(off.info.tab)] <- offspring.name
      }
    }
    # at this stage any given cell would have foloowing values
    # 0 - hap undecided; 1 - hap_1; 2 - hap_2; NA - the marker position has missing information for the given offspring
    
    # I amking a copy of the off.info.tab as a place-holder for later use. This will be explained below closer to the actual code-block.
    tmp.off.info.tab <- off.info.tab
    
    # with the following code if marker infromation for a given position is not available then this is imputed based on flanking markers. It is assumed that the physical position of markers indicate 
    # the right order of markers. 
    for (offspring.name in all.off.names) { # offspring.name <- '25208'
      bg.hap <- off.info.tab[,offspring.name]
      
      # first imputing the NAs
      splitAt <- function(x, pos) {unname(split(x, cumsum(seq_along(x) %in% pos)))} # function splits avector of marker haplotype assignments to contiguous stretches of a given haplotype.
      max.dist.1s <- 1
      inconclusive.idx <- which(is.na(bg.hap))
      idx.diff.vec <- diff(inconclusive.idx)
      diff.idx <- which(idx.diff.vec > max.dist.1s)
      NA.island.list <- splitAt(inconclusive.idx, (diff.idx + 1))
      
      if (length(NA.island.list) > 0) { # checking whether NAs even exist before running imputation
        for (iter in 1:length(NA.island.list)) { # iter <- 1
          mrkr.pos <- NA.island.list[[iter]]
          left.pos <- mrkr.pos[1]
          right.pos <- mrkr.pos[length(mrkr.pos)]
          
          # if the NA.island is at the begining of the chromosome then haplotype to the immediately downstream used to impute
          if (left.pos == 1) {
            bg.hap[mrkr.pos] <- rep(bg.hap[right.pos + 1], length(mrkr.pos))
            
            # if the NA.island is at the end of the chromosome then haplotype to the immediately upstream used to impute
          } else if (right.pos == length(bg.hap)) {
            bg.hap[mrkr.pos] <- rep(bg.hap[left.pos - 1], length(mrkr.pos))
            
            # if the NA.island is flanked by the same haplotype then flanking haplotype used to impute NAs
          } else if (bg.hap[(left.pos - 1)] == bg.hap[(right.pos + 1)]) {
            bg.hap[mrkr.pos] <- rep(bg.hap[left.pos - 1], length(mrkr.pos))
            
            # if the NA.islands are not flanked by the same haplotype then zero is assigned
          } else if (bg.hap[(left.pos - 1)] != bg.hap[(right.pos + 1)]) {
            bg.hap[mrkr.pos] <- rep(0, length(mrkr.pos))
          } else {
            stop("haplotype imputation error!")
          }
        }
      } else {
        next
      }
      
      
      # second imputing the Zeros at the ends of chromosomes and zero-islands that are flannked by the same haplotype.
      splitAt <- function(x, pos) {unname(split(x, cumsum(seq_along(x) %in% pos)))} # function splits avector of marker haplotype assignments to contiguous stretches of a given haplotype.
      max.dist.1s <- 1
      inconclusive.idx <- which(bg.hap == 0)
      idx.diff.vec <- diff(inconclusive.idx)
      diff.idx <- which(idx.diff.vec > max.dist.1s)
      zero.island.list <- splitAt(inconclusive.idx, (diff.idx + 1))
      
      if (length(zero.island.list) > 0) { # check whether undeiceded haplotype blocks exist
        for (iter in 1:length(zero.island.list)) { # iter <- 1
          mrkr.pos <- zero.island.list[[iter]]
          left.pos <- mrkr.pos[1]
          right.pos <- mrkr.pos[length(mrkr.pos)]
          
          # if the zero.island.list is at the begining of the chromosome then haplotype to the immediately downstream used to impute
          if (left.pos == 1) {
            bg.hap[mrkr.pos] <- rep(bg.hap[right.pos + 1], length(mrkr.pos))
            
            # if the zero.island.list is at the end of the chromosome then haplotype to the immediately upstream used to impute
          } else if (right.pos == length(bg.hap)) {
            bg.hap[mrkr.pos] <- rep(bg.hap[left.pos - 1], length(mrkr.pos))
            
            # if the zero.island.list is flanked by the same haplotype the zeros are imputed with the flanking haplotype
          } else if ((bg.hap[(left.pos - 1)] == bg.hap[(right.pos + 1)])) {
            bg.hap[mrkr.pos] <- rep(bg.hap[left.pos - 1], length(mrkr.pos))
            
          } else {
            next
          }
        }
      } else {
        bg.hap <- bg.hap
      }
      
      # This is a place-holder for the corrected haplotype w/o the final correction. It is useful to make this output at this point since this allows us to figure-out offspring that are 
      # sprurious matches to putative parents. If an offspring is a low match (based on relatedness metrics) they would show higher amount of small haplotype blocks. These blocks will be removed if 
      # the last correction is carried out. However, the last correction is useful to get rid of small blocks of haplotypes that are artifacts of allele-dropout and genotype-error as explained later.
      tmp.bg.hap <- bg.hap
      # plot(bg.hap, pch = 16, cex = 0.2, col = "red")
      
      
      # As noticed and commented on script 15b_..., at times algorithm that produces haplotypes.vec, can mistakenly construe that there is a double haplotype switch in areas in the genome where both 
      # focal parents are Het for considered markers and the offspring observed as HomRef or HomAlt. Higher allele drop-out in areas like this may lead the algorithm to derive a short double switch 
      # in focal parent haplotype. This can be corrected by imposing a minimum size limit for a given haplotype block. However, this may only be suitable for CO event analysis and not GC analysis. 
      # Here I am changing stretches with very small haplotype blocks due to above explained errors as zeros. The zeros are not imputed by flanking haplotypes though. 
      
      min.physical.bp.size <- 500000 # This is the minimum allowed size for a haplotype block in bps to be identified as legitimte
      for (haplo in c(1,2)) { # haplo <- 1
        splitAt <- function(x, pos) {unname(split(x, cumsum(seq_along(x) %in% pos)))} # function splits a vector of marker haplotype assignments to contiguous stretches of a given haplotype.
        max.dist.1s <- 1
        inconclusive.idx <- which(bg.hap == haplo)
        idx.diff.vec <- diff(inconclusive.idx)
        diff.idx <- which(idx.diff.vec > max.dist.1s)
        small.hap.block.list <- splitAt(inconclusive.idx, (diff.idx + 1))
        
        if (length(small.hap.block.list) > 0) { # run the next block only if there are breaks in haplotype
          for (b in 1:length(small.hap.block.list)) { # b <- 1
            block <- unlist(small.hap.block.list[b])
            
            block.start <- block[1]
            block.end <- block[length(block)]
            
            # Trying to grasp the physical distance from end-to-end of the small haplotype chunk
            start.bp.pos <- off.info.tab$POS[block.start]
            end.bp.pos <- off.info.tab$POS[block.end]
            physical.dist.hap.block <- end.bp.pos - start.bp.pos
            
            
            # following IF-block takes care of the sistuation where a small haplotype block smaller than that specified as the min cutoff occurs at the very edges of the chromosome.
            if ((block.start == 1 ) | (block.end == length(bg.hap))) {
              next
            } else if (physical.dist.hap.block < min.physical.bp.size) {
              
              # The haplotype block is removed only if the minimum critera as described above are not met.
              bg.hap[block.start:block.end] <- c(rep(0, length(block.start:block.end)))
            }
          }
        }
      }
      # plot(bg.hap, pch = 16, cex = 0.2, col = "red")
      
      
      # assign the imputed haplotypes after the final correction based on min hap block size
      if ((length(bg.hap) == length(off.info.tab[,offspring.name])) & all(!is.na(bg.hap))) {
        off.info.tab[,offspring.name] <- bg.hap
        
      } else {
        stop("imputation error with offspring background haplotype after final correcion")
      }
      
      
      # This is the imputed haplotypes without the final correction as descibed in the above code-block.
      if ((length(tmp.bg.hap) == length(tmp.off.info.tab[,offspring.name])) & all(!is.na(tmp.bg.hap))) {
        tmp.off.info.tab[,offspring.name] <- tmp.bg.hap
        
      } else {
        stop("imputation error with offspring background haplotype w/o correction")
      }
      
      
    }
    write.csv(off.info.tab, file = paste0(output.path.haplo.background, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
    cat(paste0(output.path.haplo.background, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))
    
    
    write.csv(tmp.off.info.tab, file = paste0(output.path.haplo.background, chrom.name, "_tmp_resolved_background_haplotypes_all_offspring_", par.name,".csv"), quote = FALSE, row.names = FALSE)
    cat(paste0(tmp.output.path.haplo.background, chrom.name, "_tmp_resolved_background_haplotypes_all_offspring_", par.name,".csv", " Successfully completed!", "\n"))

#   }
# }
########################################################################################################################################################################################################
# END
# The output contains focal parent haplotype identification (hap_1 or hap_2) for all of its offspring.
# Following is the output file format; 
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
########################################################################################################################################################################################################