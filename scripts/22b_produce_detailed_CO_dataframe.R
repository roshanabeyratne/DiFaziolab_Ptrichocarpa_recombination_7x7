#!/usr/bin/R
# setwd("/group/difazio/populus/gatk-7x7-Stet14/")
setwd('/gpfs/group/difazio/populus/gatk-7x7-Stet14')
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: Sep 18, 2020
# USAGE: Input files include smoothed-out bakground haplotype contribution from the focal parent. Some individuals were removed from relevant parents due to high missingness and switches with script
# 19b. Therefore the path to such curated files is "./resolved_haps/by_parent/curated/". If a parent is not corrected such as described before, the focal parent hap contribution is at
# "./resolved_haps/by_parent/"

# Following is the input file format used that has focal parent haplotype identification (hap_1 or hap_2) for all of its offspring; This version of the background haplotype contains entries for which
# minimum haplotype block has been limited to 500Kb. (refer end of script 16b_ for details).
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1

# PURPOSE: produce a detailed account of all the cross-over locations for all individuals for their focal parent haplotype. Also produce a dataframe wit summary statistics. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

hap.input.path <- "./resolved_haps/by_parent/"
curated.hap.input.path <- "./resolved_haps/by_parent/curated/"
statistics.output.path <- "./recombination_statistics/"
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
    
    # Identify loci with higher number of unresolved haplotypes and exclude them from mapping as these may lead to errors in maps (jumps, and order of markers)
    sum.mrkr.info <- apply(haplo.tab[,-(1:2)], 1, function(x){
      tmp <- (sum(x == 0) / length(x)) * 100
      return(tmp)
    })
    # cut.off <- 1 # >= 5% of individuals for this given loci has unresolved focal parental haplotype.
    tolerated.missingness <- 0.05
    num.offspring <- ncol(haplo.tab) - 2
    unrslvd.hap.outliers <- which(sum.mrkr.info > round(num.offspring * tolerated.missingness))
    
    # remove these loci from the original haplotype tab if loci to be removed are > 0
    if (sum(unrslvd.hap.outliers) > 0) {
      haplo.tab <- haplo.tab[-(unrslvd.hap.outliers), ]
    } else {
      haplo.tab <- haplo.tab
    }
    
    
    # Creating a full data table with CO information
    raw.CO.count <- NULL
    for (i in 3:ncol(haplo.tab)) { # i <- 9
      
      ind.name <- colnames(haplo.tab)[i]
      
      # creates a list of CO regions
      splitAt <- function(x, pos) {unname(split(x, cumsum(seq_along(x) %in% pos)))} # function splits a vector of marker haplotype assignments to contiguous stretches of a given haplotype.
      max.dist.1s <- 1
      zero.idx <- which(haplo.tab[,i] == 0)
      idx.diff.vec <- diff(zero.idx)
      diff.idx <- which(idx.diff.vec > max.dist.1s)
      zero.block.list <- splitAt(zero.idx, (diff.idx + 1))
      
      
      # Check if there are CO-events for an individual
      if (length(zero.block.list) > 0) { 
        for (b in 1:length(zero.block.list)) { # b <- 3
          
          # collecting information regarding the CO region
          zero.region <- as.numeric(unlist(zero.block.list[b]))
          upstrm.flnk <- haplo.tab$POS[zero.region[1] - 1]
          dwnstrm.flnk <- haplo.tab$POS[zero.region[length(zero.region)] + 1]
          CO.region.size <- (dwnstrm.flnk - upstrm.flnk)
          
          # incorporate data into the main data table
          CO.df <- rbind(CO.df, c(chrom.name, par.name, ind.name, CO.region.size, upstrm.flnk, dwnstrm.flnk))
        }
        
        # collect data for the secondary table collecting summary stats. This is not absolutely essential.
        tmp <- diff(zero.idx)
        raw.CO.count <- c(raw.CO.count, (sum(tmp > 1) + 1))
        
      } else {
        CO.df <- rbind(CO.df, c(chrom.name, par.name, ind.name, 0, NA, NA))
        raw.CO.count <- c(raw.CO.count, 0)
      }
    }
    mean.CO.count <- sum(raw.CO.count) / length(3:ncol(haplo.tab))
    variance.CO.count <- var(raw.CO.count)
    std.CO.count <- sqrt(length(raw.CO.count)) * (mean.CO.count / sqrt(variance.CO.count)) # https://llc.stat.purdue.edu/2014/41600/notes/prob1804.pdf
    
    mean.CO.count.per.chrom <- c(mean.CO.count.per.chrom, mean.CO.count)
    std.CO.count.per.chrom <- c(std.CO.count.per.chrom, std.CO.count)
    
  }
  co.count.df <- rbind(co.count.df, mean.CO.count.per.chrom)
  std.co.count.df  <- rbind(std.co.count.df, std.CO.count.per.chrom)
}


write.csv(co.count.df, file = paste0(statistics.output.path, "co.count.df.csv"), quote = FALSE, row.names = FALSE)
write.csv(std.co.count.df, file = paste0(statistics.output.path, "std.co.count.df.csv"), quote = FALSE, row.names = FALSE)
# co.count.df <-  read.csv(paste0(statistics.output.path, "co.count.df.csv"), stringsAsFactors = FALSE, header = TRUE)
# barplot(apply(co.count.df, 2, sum), main = 'Average CO per parent for the whole genome', xlab = par.name.vec)
# barplot(apply(co.count.df, 1, mean), main = 'Average cross-over counts by chromosome')

CO.df <- as.data.frame(CO.df, stringsAsFactors = FALSE)
CO.df[,4:6] <- apply(CO.df[,4:6], 2, as.numeric)
colnames(CO.df) <- c("CHROM", "PAR", "IND", "CO_SIZE", "UPSTRM_FLNK", "DWNSTRM_FLNK")
write.csv(CO.df, file = paste0(statistics.output.path, "co.detailed.df.csv"), quote = FALSE, row.names = FALSE)
# CO.df <- read.csv(paste0(statistics.output.path, "co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# CO.df <- read.csv(paste0("./co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################