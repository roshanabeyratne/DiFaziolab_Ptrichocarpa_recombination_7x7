#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: Created: June 19, 2020; Modified July 17, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name. Input files include haplotype contribution for the focal parent for all of its offspring.
# as well as the two different phasing outputs derived using adjacent marker LD and Onemap based clustering.

# focal parent haplotype contribution containing input file format as follows; 
# This is the output format for focal parent haplotype contribution in all its offspring.
# CHROM,POS,REF,ALT,1863,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553,26548,26549,26551,26552,26554,26555,26556,26557,24716,24717,24718,24719,24720,24721,24722,24723,24724
# Chr15,8393,C,T,1,NA,NA,NA,NA,NA,1,NA,0,0,NA,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,1,NA,0,1,1,0,1,1,NA,0,0,1,0,NA,0,1,NA,NA,0,NA,NA,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,0,0,1
# Chr15,8421,C,T,1,NA,NA,NA,NA,NA,1,NA,0,0,0,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,NA,NA,0,1,0,NA,0,1,NA,NA,0,0,1,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,1,0,1
# Chr15,8426,T,A,1,NA,NA,NA,NA,NA,1,0,0,0,NA,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,NA,NA,0,1,0,NA,0,1,NA,NA,0,NA,1,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,1,0,1

# Following is the input file format with two different phasing outcomes for the focal parent side by side; 
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

# PURPOSE: Derive haplotype contribution of focal parent in each of the pffspring and plot the focal parent haplotype contribution in each offspring for all offspring in the halfsib family.
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
chrom.name.vec <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

input.path.consensus <- "./consensus_phase/"
input.path.combined.haplo <- "./haplotypes/by_parent/"

output.path.jpegs <- "./consensus_phase/jpegs/"
output.path.LD.based <- "./haplotypes/LD_based/"
output.path.OM.based <- "./haplotypes/OM_based/"
########################################################################################################################################################################################################

# for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr04"
#   for (par.name in par.name.vec) { # par.name <- "6909"
    
    
    # Read-in the consensus phase of the TRUTH dataset markers. Since the marker order based on physical position is jumbled, I am correcting it as well.
    parent.phase.tab <- read.csv(file = paste0(input.path.consensus, chrom.name, "_consensus_phase_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    ordrd.idx <- order(parent.phase.tab$POS, decreasing = FALSE)
    parent.phase.tab <- parent.phase.tab[ordrd.idx, ]
    
    # read-in the haplotype contribution of each focal parent for all of its offspring. This information is compiled using script 4b and 5b. 
    parent.haplo.contrb.tab <- read.csv(file = paste0(input.path.combined.haplo, chrom.name, "_haplotypes_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    parent.haplo.contrb.tab.het <- parent.haplo.contrb.tab[(parent.haplo.contrb.tab$POS %in% parent.phase.tab$POS), ]
    colnames(parent.haplo.contrb.tab.het) <- colnames(parent.haplo.contrb.tab.het) %>% str_replace("^X", "")
    
    
    if (all(parent.phase.tab$POS == parent.haplo.contrb.tab.het$POS)) { # making sure that the marker order is right between the two datasets before cbind
      
      tmp.combined.tab <- cbind(parent.phase.tab, parent.haplo.contrb.tab.het[6:ncol(parent.haplo.contrb.tab.het)])
      tmp.combined.tab <- tmp.combined.tab[!is.na(tmp.combined.tab$OM.hap.1), ] # obtain markers that have a complete dataset outputted from the OM pipeline. 
      
      
      # This code section is implemented to loop over haplotypes derived based on the adjacent marker LD and Onemap methods
      for (pipline in c("LD", "OM")) { #  pipline <- "LD"
        
        # identify whether a given allele contributed by focal parent to the offspring comes rom HAP1 or HAP2
        temp.LD <- t(apply(tmp.combined.tab[,-(1:2)], 1, function(x, pl = pipline){ # The first two columns (CHROM and POS) are ommitted from the input to this function.
          LD.hap.1 <- as.numeric(x[1])
          LD.hap.2 <- as.numeric(x[2])
          OM.hap.1 <- as.numeric(x[3])
          OM.hap.2 <- as.numeric(x[4])
          off.haps <- as.numeric(x[5:length(x)])
          
          if (pl == "LD") {
            conversion.key <- c(LD.hap.1, LD.hap.2)
            conversion.value <- c(1, 2)
          } else if (pl == "OM") {
            conversion.key <- c(OM.hap.1, OM.hap.2)
            conversion.value <- c(1, 2)
          } else {
            stop("pipiline not specified for haplotypes!")
          }
          
          if (any(is.na(conversion.key))) { # This takes care of instances where the parental haplotype is missing - NA
            hap.vec <- rep(0, length(off.haps))
          } else {
            hap.vec <- (sapply(off.haps, function(y, key = conversion.key, value = conversion.value){
              if (is.na(y)) { # if the allle contribution of focal parent is not determined in an offspring at the marker it is assigned a zero. 
                hap <- 0
              } else {
                hap <- value[which(key == y)]
              }
              return(hap)
            }))
          }
          return(hap.vec)
        }))
        
        
        if ((nrow(temp.LD) == nrow(tmp.combined.tab)) & (ncol(temp.LD) == (ncol(tmp.combined.tab) - 6))) {
          
          # output a colored image that identifies focal parent haplotype in the offspring. This could be HAP1 - 1, HAP2 - 2 or missing - 0
          # pdf(paste0(output.path.jpegs, chrom.name, "_haplotypes_", pipline, "_", par.name, ".pdf"))
          # image(temp.LD, xlab = "markers", ylab = "offspring")
          # dev.off()
          # cat(paste0(output.path.jpegs, chrom.name, "_haplotypes_", pipline, "_", par.name, ".pdf"), " Successfully completed!")
          
          # specify the output paths based on the pipeline used (either LD or OM)
          if (pipline == "LD") {
            phased.offspring.path <- output.path.LD.based
          } else if (pipline == "OM") {
            phased.offspring.path <- output.path.OM.based
          } else {
            stop("pipiline not identified for outputpath!")
          }
          
          
          phased.offspring.tab <- cbind(tmp.combined.tab[, 1:6], as.data.frame(temp.LD, stringsAsFactors = FALSE))
          for (j in 1:ncol(phased.offspring.tab)) {
            phased.offspring.tab[,j] <- unlist( phased.offspring.tab[,j])
          }
          
          if (ncol(phased.offspring.tab) == ncol(tmp.combined.tab)) {
            colnames(phased.offspring.tab) <- colnames(tmp.combined.tab)
            write.csv(phased.offspring.tab, file = paste0(phased.offspring.path, chrom.name, "_phased_offspring_", pipline, "_",  par.name, ".csv"), row.names = FALSE, quote = FALSE)
            cat(paste0(phased.offspring.path, chrom.name, "_phased_offspring_", par.name, ".csv"), " Successfully completed!", "\n")
            
          } else {
            stop("columns not equal between tmp.combined.tab and phased.offspring.tab")
          }
          
        } else {
          stop("error in the dimensions of temp.LD table which may relate to the number of offspring for the focal parent!")
        }
      }
    } else {
      stop("parent.phase.tab and parent.haplo.contrb.tab.het do not have same number of markers or the marker order is not contiguous!")
    }
    
#   }
# }
########################################################################################################################################################################################################
# END 
########################################################################################################################################################################################################
# Following is the output format produced by above script that contains either 0,1 or 2
# CHROM,POS,LD.hap.1,LD.hap.2,OM.hap.1,OM.hap.2,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553,24757,24758,24759,24760,24763,24765,24772,24773,24777,24778,24781,24782,24783
# Chr02,6617,0,1,0,1,1,2,2,2,1,1,2,2,1,2,2,1,2,1,2,2,1,1,1,2,2,1,1,1,1,2,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,0,1,1,1,0,2,1,2,2,2,1,2,1,1,2,1,2,2,1,2,2,2,2,1,1,1,1,2,2,1,1,1,2,2,2,2,1,1,2,1,1,2,1,1,2,2,1,2
########################################################################################################################################################################################################
