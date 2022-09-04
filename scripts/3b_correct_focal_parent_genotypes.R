#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
#######################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: June 08, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name 
# Input file format is in; 
# example input file format for the fulsib family

# CHROM, POS, REF, ALT, 4593, 2393, 26241, 26244, 26246, 26249, 26250, 26251, 26252, 26253, 26274, 26278, 26293, 26302, 26304
# Chr05, 7196, A, C, 1, 2, 1, NA, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1
# Chr05, 8867, A, G, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1
# Chr05, 9407, A, C, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 0, 1

# PURPOSE:This code attempts to correct the genotype errors of parents based on trio information based on a Liklihood based method for the whole halfsib family. 
# This is identified as the "Dobzhansky scoring procedure" - REF: https://www.genetics.org/content/genetics/119/2/465.full.pdf
#######################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
#######################################################################################################################################################################################################
args <- commandArgs(TRUE)
chrom.name <- args[1]
par.name <- args[2]
#######################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
# parentage.file.name <- "./lists/7x7_parentage_v4.txt"
parentage.file.name <- "./scripts_lists/7x7_parentage_v4.txt"
parantage_v4.tab <- read.table(file = parentage.file.name, stringsAsFactors = FALSE, header = TRUE)
input.path <- "./genotypes/"
output.path <- "./genotypes/parents/"
#######################################################################################################################################################################################################

# for (chrom.name in chrom.name.vec[1:19]) { # chrom.name <- chrom.name.vec[2]
  # for (par.name in par.name.vec[1:14]) { # par.name <- "2515"
    
    # following code block iterates through 7 full-sib famillies for a given focal parent and makes a genotype file for the whole halfsib family including all the parents. 
    if (any(par.name %in% par.name.vec[1:7])) {
      
      other.par.names.vec <- par.name.vec[8:14]
      halfsib.fam.geno.tab <- NULL
      for (other.name in other.par.names.vec[1:7]) { # other.name <- "6909"
        digit.tab <- as.matrix(read.csv(file = paste0(input.path, chrom.name, "_fulsib_genotypes_", par.name, "X", other.name, ".csv"), stringsAsFactors = FALSE, header = TRUE) )
        # cat(nrow(digit.tab), "\n")
        uniq.col.names <-  colnames(digit.tab)[!(colnames(digit.tab) %in% colnames(halfsib.fam.geno.tab))]
        halfsib.fam.geno.tab <- cbind(halfsib.fam.geno.tab, digit.tab[, uniq.col.names])
        
      }
    } else {
      
      other.par.names.vec <- par.name.vec[1:7]
      halfsib.fam.geno.tab <- NULL
      for (other.name in other.par.names.vec[1:7]) { # other.name <- "1863"
        digit.tab <- as.matrix(read.csv(file = paste0(input.path, chrom.name, "_fulsib_genotypes_", other.name, "X", par.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)) 
        # cat(nrow(digit.tab), "\n")
        uniq.col.names <-  colnames(digit.tab)[!(colnames(digit.tab) %in% colnames(halfsib.fam.geno.tab))]
        halfsib.fam.geno.tab <- cbind(halfsib.fam.geno.tab, digit.tab[, uniq.col.names])
        
      }
    }
    
    colnames(halfsib.fam.geno.tab) <- colnames(halfsib.fam.geno.tab) %>% str_replace("X", "")
    parents.in.hs.fam <- par.name.vec[par.name.vec %in% c(par.name, other.par.names.vec)]
    offspring.names <- colnames(halfsib.fam.geno.tab)[ !(colnames(halfsib.fam.geno.tab) %in% c(parents.in.hs.fam, colnames(halfsib.fam.geno.tab)[1:4])) ]
    # col.name.order <- c(colnames(halfsib.fam.geno.tab)[1:4], parents.in.hs.fam,  offspring.names)
    col.name.order <- colnames(halfsib.fam.geno.tab)
    halfsib.fam.geno.tab <- as.data.frame(halfsib.fam.geno.tab[, col.name.order], stringsAsFactors = FALSE)
    rm(digit.tab)
    ####################################################################################################################################################################################################
    
    # This rate of error was calculated for the offspring in an earlier script - allele_drop_out.R
    geno.error.rate <- 0.1
    
    # Defining probability scores for each trio genotypes. If any of the three individuals have a missing genotype, they are not included in the calculations. 
    prob.scores.0 <- c(1, geno.error.rate, geno.error.rate, 0.5, 0.5, geno.error.rate, geno.error.rate, 1, geno.error.rate)
    par.geno.mat.0 <- matrix(prob.scores.0, nrow = 3, ncol = 3, byrow = TRUE)
    rownames(par.geno.mat.0) <- 0:2
    colnames(par.geno.mat.0) <- 0:2
    
    # par.geno.mat.0
    #       0    1    2
    # 0   1.00 0.15 0.15
    # 1   0.50 0.50 0.15
    # 2   0.15 1.00 0.15
    
    prob.scores.1 <- c(0.5, 0.5, geno.error.rate, 0.25, 0.5, 0.25, geno.error.rate, 0.5, 0.5)
    par.geno.mat.1 <- matrix(prob.scores.1, nrow = 3, ncol = 3, byrow = TRUE)
    rownames(par.geno.mat.1) <- 0:2
    colnames(par.geno.mat.1) <- 0:2
    
    # par.geno.mat.1
    #       0   1    2
    # 0   0.50 0.5 0.15
    # 1   0.25 0.5 0.25
    # 2   0.15 0.5 0.50
     
    prob.scores.2 <- c(geno.error.rate, 1, geno.error.rate, geno.error.rate, 0.5, 0.5, geno.error.rate, geno.error.rate, 1)
    par.geno.mat.2 <- matrix(prob.scores.2, nrow = 3, ncol = 3, byrow = TRUE)
    rownames(par.geno.mat.2) <- 0:2
    colnames(par.geno.mat.2) <- 0:2
    
    
    # following are the weighting matrices. rows refer other.parent's genotypes and columns an offspring's genotypes
    # par.geno.mat.2
    #       0    1    2
    # 0   0.15 1.00 0.15
    # 1   0.15 0.50 0.50
    # 2   0.15 0.15 1.00
    
    # par.1.expected.geno.vec.2 <- par.1.expected.geno.vec
    focal.par.geno.vec <- NULL
    for (i in 1:nrow(halfsib.fam.geno.tab)) { # i <- 35 
      geno.0 <- 1
      geno.1 <- 1
      geno.2 <- 1
      
     focal.par.observed.geno <- as.character(as.integer(halfsib.fam.geno.tab[i, par.name]))
      for (offspring in offspring.names) { # offspring <- "25866"
        
        off.geno <- as.character(as.integer(halfsib.fam.geno.tab[i, offspring]))
        
        if (!is.na(off.geno)) {
          tmp.par <- parantage_v4.tab[which(parantage_v4.tab$OldClone == offspring), 3:4]
          other.par.geno <- as.character(as.integer(halfsib.fam.geno.tab[i, tmp.par[tmp.par != par.name]]))
          
          if (!is.na(other.par.geno)) {
            
            geno.0 <- geno.0 * par.geno.mat.0[other.par.geno, off.geno]
            geno.1 <- geno.1 * par.geno.mat.1[other.par.geno, off.geno]
            geno.2 <- geno.2 * par.geno.mat.2[other.par.geno, off.geno]
            
          } else {
            # if the other parent's genotype is missing then the offspring is not included in the likelihood calclations
            next
          }
        } else {
          # if the offspring's genotype is missing then this offspring is not included in the likelihood calculation
          next
        }
      }
      
      data.likelihood <- (c(geno.0, geno.1, geno.2))
      focal.par.expected.geno <- as.character(c(0:2)[which.max(data.likelihood)])
      
      if ((!is.na(focal.par.observed.geno)) & (focal.par.expected.geno != focal.par.observed.geno)) {
        
        max.likelihood <- max(data.likelihood)
        obs.geno.likelihood <- data.likelihood[as.integer(focal.par.observed.geno) + 1]
        log.like.ratio <- log(max.likelihood / obs.geno.likelihood) # This is a note to self. https://stephens999.github.io/fiveMinuteStats/likelihood_ratio_simple_models.html #https://www.g3journal.org/content/7/5/1393
       
        if (log.like.ratio > 5) { # focal parent genotype is changed only if the expected genotype is 100 times more likely thatn the observed
          focal.par.geno <-  focal.par.expected.geno 
        } else {
          focal.par.geno <- focal.par.observed.geno
        }
      } else if (is.na(focal.par.observed.geno)) {
        focal.par.geno <- focal.par.expected.geno
      } else {
        focal.par.geno <- focal.par.observed.geno
      }
      focal.par.geno.vec <- c(focal.par.geno.vec, focal.par.geno)
    }
    # sum(focal.par.geno.vec != as.integer(halfsib.fam.geno.tab[1:500, par.name]), na.rm = TRUE)
    # which(focal.par.geno.vec != as.integer(halfsib.fam.geno.tab[1:500, par.name]))
    
    obsrvd.v.exp <- as.data.frame(cbind(halfsib.fam.geno.tab[, 1:4], halfsib.fam.geno.tab[, par.name], focal.par.geno.vec), stringsAsFactrs = FALSE)
    colnames(obsrvd.v.exp) <- c(colnames(halfsib.fam.geno.tab)[1:5], "CORRECTED_GENOTYPE")
    write.csv(obsrvd.v.exp, file = paste0(output.path, chrom.name, "_corrected_genotypes_", par.name, ".csv"), quote = FALSE, row.names = FALSE)
    
  # } # par.name bracket
# } # chrom.name bracket
#######################################################################################################################################################################################################
# END
# Following is an example of the output format; 
#     CHROM,POS,REF,ALT,1863,CORRECTED_GENOTYPE
#     Chr12,    7909,C,A, 0,0
#     Chr12,    7914,G,T, 0,0
#     Chr12,   14391,C,T, 1,1
#     Chr12,   14660,C,G, 2,2
#     Chr12,   14676,A,T, 1,1
#     Chr12,   14721,T,A, 1,1
#     Chr12,   15155,G,T, 1,1
#     Chr12,   15321,T,C, 2,2
#     Chr12,   15977,T,C, 1,1