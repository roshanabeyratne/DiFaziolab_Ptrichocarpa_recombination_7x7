#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: March 16, 2020, Modified: June 14, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name. Input files include haplotype contribution for the focal parent. 

# This is the output format for focal parent haplotype contribution in all its offspring.
# CHROM,POS,REF,ALT,1863,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553,26548,26549,26551,26552,26554,26555,26556,26557,24716,24717,24718,24719,24720,24721,24722,24723,24724,24725,24726,24727,24728,24729,24730,24731,24732,24734,26564,26565,26566,26567,26570,26572,26573,26574,26575,26577,26578,26581,26584,26587,26589,26591,26592,26593,25998,26598,26600,26602,26603,26605,26606,26607,26608,26610,26611,26612,26613,26615,26616,26620,26623,26624,26625,26798,24736,24740,24745,24747,26631,26635,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677
# Chr15,8393,C,T,1,NA,NA,NA,NA,NA,1,NA,0,0,NA,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,1,NA,0,1,1,0,1,1,NA,0,0,1,0,NA,0,1,NA,NA,0,NA,NA,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,0,0,1,1,0,0,NA,NA,0,0,0,1,0,NA,NA,NA,NA,NA,0,1,1,1,NA,0,0
# Chr15,8421,C,T,1,NA,NA,NA,NA,NA,1,NA,0,0,0,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,NA,NA,0,1,0,NA,0,1,NA,NA,0,0,1,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,1,0,1,1,0,0,1,NA,0,0,0,1,0,NA,NA,NA,NA,NA,0,1,1,NA,NA,NA,0
# Chr15,8426,T,A,1,NA,NA,NA,NA,NA,1,0,0,0,NA,0,NA,NA,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,NA,NA,0,1,0,NA,0,1,NA,NA,0,NA,1,NA,NA,0,1,1,1,1,0,1,1,NA,0,1,1,1,0,0,1,0,NA,0,0,1,0,1,1,0,0,1,NA,0,0,0,1,0,NA,NA,NA,NA,NA,NA,1,1,NA,NA,0,0
# Chr15,8805,G,A,1,NA,0,NA,NA,NA,1,NA,0,0,0,0,NA,0,1,0,0,1,1,1,1,NA,NA,0,0,NA,0,1,0,1,1,0,0,NA,NA,0,1,1,0,1,1,0,NA,0,NA,0,0,0,1,NA,0,0,0,1,NA,NA,0,1,1,1,1,0,1,1,0,NA,1,1,1,0,0,1,NA,0,0,0,1,0,1,1,0,0,1,0,0,0,0,1,NA,NA,NA,0,NA,NA,NA,1,1,NA,NA,NA,0

# PURPOSE: Phase heterozygous parental markers based on LD between adjacent markers. Physical position of markers as aligned to the Stettler-14 reference genome is assumed accurate for this analysis. 
# All offspring for a given focal parent are used and focal parent's allele contribution towards an offspring is derived using trio genotype information for a given marker.
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
input.path.combined.haplo <- "./haplotypes/by_parent/"
output.path.parent.haplo <- "./haplotypes/parent_haplotypes/"
output.path.hs.fam <- "./haplotypes/by_halfsib_family/"
########################################################################################################################################################################################################

# for (chrom.name in chrom.name.vec[1:19]) { # chrom.name <- chrom.name.vec[3]
  # for (par.name in par.name.vec[1:14]) { # par.name <- par.name.vec[8]
    parent.haplo.tab <- read.csv(file = paste0(input.path.combined.haplo, chrom.name, "_haplotypes_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    colnames(parent.haplo.tab) <- colnames(parent.haplo.tab) %>% str_replace(".GT", "") %>% str_replace("^X", "")

    
    par.het.haplo.tab <- parent.haplo.tab[parent.haplo.tab[,par.name] == 1, ] # Selecting markers for which focal parent is heterozygous
    par.het.haplo.tab <- par.het.haplo.tab[!is.na(par.het.haplo.tab$POS), ] # 20841 markers for par.name = "2283
    
    # hist(par.het.haplo.tab$POS, breaks = 100000)
    # par.het.haplo.tab <- par.het.haplo.tab[-c(3264, 14012), ]
    
    
    # identify indices of markers for which information content in offspring is low. Following conditions would make focal parent allele contribution in offspring unreliable. Therefore NA is assigned
    # offspring genotype is missing or the other parent of the offspring has a missing genotype.
    # offspring genotype is not consistent with parental genotyes (mendelian violation).
    less.informative.marker.count <- apply(par.het.haplo.tab[,6:ncol(par.het.haplo.tab)], 1, function(x){ # offspring information start from column-6 onwards.
      count.Nas <- sum(is.na(x))
      return(count.Nas)
    })
    # hist(less.informative.marker.count, breaks = 100)
    less.informative.marker.idx <- which(less.informative.marker.count > (ncol(par.het.haplo.tab) * 0.25)) # removing markers that have more than 25% of offspring with missing information (NA)
    par.het.haplo.tab <- par.het.haplo.tab[-less.informative.marker.idx, ]
    
    
    ########################################################################################################### LD first round START
    # following section calculates adjacent marker LD. This analysis assumes that the physical marker positions of the Stet-14 referene genome accurately reflects physical positions of the parents of 
    # the 7x7 cross. i.e. adjacent markers are considered to be the ones for which the physical position is closest. This is the first round of LD calculations between markers.
    ld.df <- NULL
    for (j in 6:ncol(par.het.haplo.tab)) {  # j <- 13
     
      # making vectors with adjacent marker pairs for offspring. offspring columns start at column 6
      tmp.vec1 <-  par.het.haplo.tab[1:(nrow(par.het.haplo.tab) - 1), j]
      tmp.vec2 <-  par.het.haplo.tab[2:nrow(par.het.haplo.tab), j]
      tmp.ld.vec <- paste0(tmp.vec1, tmp.vec2)
      ld.df <- cbind(ld.df, tmp.ld.vec)
    }
    ld.df <- as.data.frame(ld.df, stringsAsFactors = FALSE)
    
    
    LD.tab <- as.data.frame(t(apply(ld.df, 1, function(x){ # calculate the LD between adjacent markers.
      geno.vec <- as.character(x)
      X11 <- sum(geno.vec == "11", na.rm = TRUE)
      X00 <- sum(geno.vec == "00", na.rm = TRUE)
      X10 <- sum(geno.vec == "10", na.rm = TRUE)
      X01 <- sum(geno.vec == "01", na.rm = TRUE)
      toal.informative.loci <- (X11 + X00 + X10 + X01) # this is the total number of informative offspring per given loci.
      LD <- ((X11 * X00) - (X10 * X01)) / (toal.informative.loci ^ 2) # As per equation in, Genetics of Populations - Philip W. Hedrick, 4th Edition (page 528).
      return.vec <- c(LD, toal.informative.loci)
      return(return.vec)
    })), stringsAsFactors = FALSE)
    
    
    position.ids <- par.het.haplo.tab$POS[1:(nrow(par.het.haplo.tab) - 1)] # physical position of the upstream marker
    pos.diffs <- (par.het.haplo.tab$POS[2:nrow(par.het.haplo.tab)] - par.het.haplo.tab$POS[1:(nrow(par.het.haplo.tab) - 1)]) # the distance between the markers for which LD is calculated
    LD.tab <- cbind(position.ids, pos.diffs, LD.tab)
    colnames(LD.tab) <- c("POS_ID", "POS_DIFF", "LD", "INF_INDS")
    write.csv(LD.tab, file = paste0(output.path.parent.haplo, chrom.name, "_parent_haplotypes_", par.name, ".csv"), row.names = FALSE, quote = FALSE)
    
    
    ########################################################################################################### LD second round START
    # I have decided to set 0.1 as the cutoff for adjacent LD. If LD values obtained for a given marker with its upstream partner and downstream partner are both lower than the cutoff, this marker is 
    # is onsidered erroneous. Furthermore, if for a given loci pair LD value is considered a tolerable estimate if the minimum number of offspring used in the LD calculation is >= 20.
    # larger sections of markers that show low LD to the flanking regions (structural anomalies of the focal parent) are not treated in this script.
    ld.cut.off <- 0.1
    min.num.iform.offspring <- 20
    
    ids.to.exclude <- NULL
    for (pos in 2:nrow(LD.tab)) {
      ld.pos.prev <-  abs(LD.tab$LD[pos - 1])
      ld.pos.curr <- abs(LD.tab$LD[pos])
      
      if ((ld.pos.prev < ld.cut.off) & (ld.pos.curr < ld.cut.off)) {
        ids.to.exclude <- c(ids.to.exclude, LD.tab$POS_ID[pos])
      }
    }
    
    idx.to.exclude.2 <- which(is.na(LD.tab$LD) | (LD.tab$INF_INDS <= min.num.iform.offspring) | LD.tab$POS_ID %in% ids.to.exclude) #(ncol(ld.df)*min.num.iform.offspring)
    pos.ids.to.exclude <- LD.tab$POS_ID[idx.to.exclude.2]
    par.het.haplo.tab2 <- par.het.haplo.tab[!(par.het.haplo.tab$POS %in% pos.ids.to.exclude), ]
    
    
    # LD between adjacent markers are calculated again, after some of the markers are removed from further phasing. This is section is the same as the first. 
    ld.df2 <- NULL
    for (j in 6:ncol(par.het.haplo.tab2)) { # making vectors with adjacent marker pairs for offspring. offspring columns start at column 6
      # j <- 13
      tmp.vec1 <-  par.het.haplo.tab2[1:(nrow(par.het.haplo.tab2) - 1), j]
      tmp.vec2 <-  par.het.haplo.tab2[2:nrow(par.het.haplo.tab2), j]
      tmp.ld.vec <- paste0(tmp.vec1, tmp.vec2)
      ld.df2 <- cbind(ld.df2, tmp.ld.vec)
    }
    ld.df2 <- as.data.frame(ld.df2, stringsAsFactors = FALSE)
    
    LD.tab2 <- as.data.frame(t(apply(ld.df2, 1, function(x){ # calculate the LD between adjacent markers
      geno.vec <- as.character(x)
      X11 <- sum(geno.vec == "11", na.rm = TRUE)
      X00 <- sum(geno.vec == "00", na.rm = TRUE)
      X10 <- sum(geno.vec == "10", na.rm = TRUE)
      X01 <- sum(geno.vec == "01", na.rm = TRUE)
      toal.informative.loci <- (X11 + X00 + X10 + X01)
      LD <- ((X11 * X00) - (X10 * X01)) / (toal.informative.loci ^ 2) # As per equation in, Genetics of Populations - Philip W. Hedrick, 4th Edition (page 528)
      return.vec <- c(LD, toal.informative.loci)
      return(return.vec)
    })), stringsAsFactors = FALSE)
    
    position.ids <- par.het.haplo.tab2$POS[1:(nrow(par.het.haplo.tab2) - 1)]
    pos.diffs <- (par.het.haplo.tab2$POS[2:nrow(par.het.haplo.tab2)] - par.het.haplo.tab2$POS[1:(nrow(par.het.haplo.tab2) - 1)])
    LD.tab2 <- cbind(position.ids, pos.diffs, LD.tab2)
    colnames(LD.tab2) <- c("POS_ID", "POS_DIFF", "LD", "INF_INDS")
    write.csv(LD.tab2, file = paste0(output.path.parent.haplo, chrom.name, "_parent_haplotypes_2_", par.name, ".csv"), row.names = FALSE, quote = FALSE)
    
    
    ########################################################################################################### Focal parent phasing based on adjacent marker LD
    # following code block uses the adjecent marker LD to identify the phase of the focal parent. Since low LD markers are removed ideally with above scripts, here I have considered any positive LD 
    # value to indicate markers in coupled phase and repulsion if the adjacent marker LD is negative. Here hap1 is considered arbitrarily to start with 0. 
    hap1 <- 0
    starting.hap <- 0
    for (idx in 1:nrow(LD.tab2)) {
      
      if (LD.tab2$LD[idx] < 0) {
        adjusted.haplotype <- abs(starting.hap - 1)
        hap1 <- c(hap1, adjusted.haplotype)
        starting.hap <- adjusted.haplotype
      
        } else if (LD.tab2$LD[idx] >= 0) {
        adjusted.haplotype <- starting.hap
        hap1 <- c(hap1, adjusted.haplotype)
        starting.hap <- adjusted.haplotype
      }
      
      hap2 <- (rep(1, length(hap1)) - hap1)
    }
    
    
    halfsib.family.haplo.tab <- cbind(par.het.haplo.tab2[, 1:2], hap1, hap2, par.het.haplo.tab2[, 6:(ncol(par.het.haplo.tab2))])
    halfsib.family.haplo.tab <- as.data.frame(halfsib.family.haplo.tab, stringsAsFactors = FALSE)
    
    # following code block identifies to which focal parental haplotype thier focal parental contribution fall.
    offspring.haplo.df <- t(apply(halfsib.family.haplo.tab, 1, function(x) {
      x <- as.character(x)
      chrom <- x[1]
      pos <- x[2]
      hap1 <- as.character(x[3])
      hap2 <- as.character(x[4])
      
      haplo.vec <- x[5:length(x)]
      offspring.haplo.per.pos <- sapply(haplo.vec, function(y, h1 = hap1, h2 = hap2) {
        y <- as.numeric(y)
        if (is.na(y)) {
          temp.var <- NA
        } else {
          if (y == h1) {
            temp.var <- 1
          } else if (y == h2) {
            temp.var <- 2
          }
        }
        return(temp.var)
        })

      full.line.per.pos <- as.character(c(chrom, pos, hap1, hap2, offspring.haplo.per.pos))
      return(full.line.per.pos)
      }))
    
    
    offspring.haplo.df <- as.data.frame(offspring.haplo.df, stringsAsFactors = FALSE)
    colnames(offspring.haplo.df) <- colnames(halfsib.family.haplo.tab)
    write.csv(offspring.haplo.df, file = paste0(output.path.hs.fam, chrom.name, "_halfsibfamily_haplotypes_", par.name, ".csv"), row.names = FALSE, quote = FALSE)
    cat(paste0(output.path.hs.fam, chrom.name, "_halfsibfamily_haplotypes_", par.name, ".csv"), " Successfully completed!", "\n")
  # }
# }
########################################################################################################################################################################################################
###END
    # This is the first output file format of the focal parent haplotypes
    # POS_ID,POS_DIFF,LD,INF_INDS
    # 20263,4399,-0.216581446311176,74
    # 24662,28060,-0.202546296296296,72
    # 52722,34543,0.16764061358656,74
    # 87265,155,0.184441197954711,74
    # 87420,1416,0.170349131388092,77
    # 88836,6361,0.178545507330527,83
    # 95197,4046,0.196014277215943,82
    # 99243,1789,0.189866666666667,75
    # 101032,10668,0.195158566335147,73
    
    
    
    # This is the second output format of the parental haplotypes in the focal parent's offspring
    # CHROM,POS,hap1,hap2,24757,24758,24759,24760,24763,24765,24772,24773,24777,24778,24781,24782,24783,24784,24785,24793,24794,26705,25431,25436,25437,25438,25439,25443,25448,25452,25453,25454,25455,25460,25465,25470,25472,25476,25483,25322,25325,25332,25334,25336,25339,25503,25504,25505,25506,25507,25509,25510,25511,25512,25513,25514,26979,24802,24804,24805,24806,24809,24810,24811,24812,24817,24818,24823,24825,24828,24830,24831,24833,24834,26998,25515,25516,25517,25518,25519,25520,25522,25523,25526,25527,25533,25535,24837,24839,24840,24841,24842,24844,24847,24848,24849,24850,24851,24852,24853,24859,24864,24865,24867,24868,25530,25531,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492,25493,25495,25498,25500
    # Chr15,  173833,0,1,NA,1,1,1,1,1,1,NA,1,NA,1,2,1,NA,1,1,1,1,NA,NA,2,2,1,1,1,1,2,1,2,1,2,1,1,1,NA,1,2,NA,2,1,2,NA,NA,2,NA,NA,2,NA,1,NA,2,1,2,2,NA,2,1,NA,1,2,1,1,1,2,1,1,2,1,2,1,2,1,2,2,2,1,NA,2,1,1,1,2,2,2,1,2,1,2,1,2,2,2,2,2,2,1,NA,2,1,1,NA,1,1,1,1,1,1,1,1,2,1,1,2,NA,2,2,2,2,2,1
    # Chr15,  174711,0,1,2,1,1,1,1,1,1,NA,1,2,1,2,1,2,1,2,1,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,2,1,2,NA,2,1,2,NA,NA,2,NA,NA,NA,2,1,NA,2,1,2,2,2,2,1,2,1,2,1,1,1,2,1,1,2,1,2,1,NA,2,NA,2,2,1,NA,NA,NA,1,1,1,2,2,1,1,1,2,1,NA,2,2,2,2,NA,2,2,2,1,1,2,1,2,2,1,1,1,2,1,1,1,NA,2,NA,2,NA,2,2,NA,2
    # Chr15,  178446,1,0,2,1,2,NA,NA,1,1,NA,NA,2,1,2,1,2,NA,2,1,1,2,2,2,2,1,NA,2,1,2,1,2,1,2,2,1,2,2,1,2,2,2,NA,2,2,2,2,2,2,1,2,1,2,2,1,2,2,2,2,1,2,1,2,1,1,1,2,1,1,1,1,2,1,2,2,1,2,2,1,2,1,2,1,1,2,2,2,1,2,1,2,2,2,2,2,2,1,2,2,2,2,1,1,2,1,2,NA,1,1,1,NA,1,2,NA,NA,NA,NA,NA,NA,NA,NA,NA,2
########################################################################################################################################################################################################
########################################################################################################### test code block
# testing.df <- read.csv(file = paste0("./haplotypes/by_halfsib_family/", chrom.name, "_halfsibfamily_haplotypes_", par.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)
# positions.to.exclude <- par.het.haplo.tab$POS[c(3264, 14012)]
# haplotype.first <- read.csv(file = paste0("./haplotypes/by_halfsib_family/", chrom.name, "_halfsibfamily_haplotypes_2_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
# haplotype.first <- haplotype.first[, 1:4]
# haplotype.second <- read.csv(file = paste0("./haplotypes/by_halfsib_family/", chrom.name, "_halfsibfamily_haplotypes_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
# haplotype.second <- haplotype.second[, 1:4]
# haplotype.second.subset <- haplotype.second[haplotype.second$POS %in% haplotype.first$POS, ] 
# all(haplotype.second.subset$hap2 == haplotype.first$hap1)
#######################################################################################################################################################################################################
