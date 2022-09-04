#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: June 30, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name, mother.name and father.name 
# input file format is as follows, 

# CHROM	POS	REF	ALT	1863.GT	1909.GT	1950.GT	2048.GT	2066.GT	2283.GT	2365.GT	2393.GT	24708.GT	24709.GT	24710.GT	24711.GT	24712.GT
# Chr01	7940	T	 C	C/C	    C/C	    C/C	    T/C	    C/C	    T/C	    T/T	    C/C	    T/T	      T/T	      T/T       T/T       T/T

# PURPOSE:This code converts the GATK VariantsToTable format to 012 format for the truth data-set and for individuals in a ful-sib family. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
args <- commandArgs(TRUE)
chrom.name <- args[1]
mother.name <- args[2]
father.name <- args[3]
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
parentage.file.name <- "./scripts_lists/7x7_parentage_v4.txt"
parantage_v4.tab <- read.table(file = parentage.file.name, stringsAsFactors = FALSE, header = TRUE)
output.path.012 <- "./genotypes/"
########################################################################################################################################################################################################
# for (chrom.name in chrom.name.vec[2]) { # chrom.name <- chrom.name.vec[2]
  
  # reading the SNP truth-dataset in GATK VatriantToTable format for each chromosome
  snp.file.name <- paste0("./truth_data/7x7_Stet14_SNP_truthdataset_", chrom.name, ".table") 
  tmp.variant.table <- read.table(file = snp.file.name, stringsAsFactors = FALSE, header = TRUE)
  colnames(tmp.variant.table) <- colnames(tmp.variant.table) %>% str_replace(".GT", "") %>% str_replace("^X", "") # removing GATK suffixes and 'X' in front of the clone-names when they are read to R
  
  
  # for (mother.name in par.name.vec[1:7]) { # mother.name <- "2283"
    # for (father.name in par.name.vec[8:14]) { # father.name <- "7073"
      
      
      # Identify offspring clone names/individuals for the full-sib family given the most recent parentage assignement and then subset the truth-dataset for the full-sib family.  trio information.
      offspring.idx <- which((parantage_v4.tab$matched_mother == mother.name) & (parantage_v4.tab$matched_father == father.name))
      tmp.offspring.names <- as.character(parantage_v4.tab$OldClone[offspring.idx]) # old clone names used in the 7x7 trial are used here
      offspring.names <- sort(colnames(tmp.variant.table)[which(colnames(tmp.variant.table) %in% tmp.offspring.names)])
      tmp.subset.var.tab <- tmp.variant.table[, c(colnames(tmp.variant.table)[1:4], mother.name, father.name, offspring.names)]
      
      
      # Following is a function to convert genotypes in G/C format into 012 format.
      tmp.tab <- as.data.frame(t(apply(tmp.subset.var.tab, 1, function(x) {
        x <- as.character(x)
        ref.allele <- x[3]
        alt.allele <- x[4]
        hom.ref.geno <- paste0(ref.allele, "/", ref.allele)
        hom.var.geno <- paste0(alt.allele, "/", alt.allele)
        het.geno <- paste0(ref.allele, "/", alt.allele)
        missing.geno <- paste0(".", "/", ".")
        template.vec <- c(hom.ref.geno, het.geno, hom.var.geno, missing.geno)
        convertion.vec <- c(0, 1, 2, NA) # 0-HomRef, 1-Het, 2-HomAlt. The missing genotypes denoted as "./." are assigned NA
        
        bin.vec <- as.vector(sapply(x[-(1:4)], function(y, template = template.vec, conversion = convertion.vec) {
          genotype <- y
          binary <- conversion[which(template %in% genotype)]
          return(binary)
        }))
        return(bin.vec)
      })))
      
      
      if (nrow(tmp.tab) < ncol(tmp.tab)) { # this IF block takes care of complications in function due to one offspring for a given fulsib family
        tmp.tab <- as.data.frame(t(tmp.tab), stringsAsFactors = FALSE)
      }
      

      digit.tab <- cbind(tmp.subset.var.tab[, 1:4], tmp.tab)
      if (ncol(digit.tab) == ncol(tmp.subset.var.tab)) {
        colnames(digit.tab) <- c(colnames(tmp.subset.var.tab))
        write.csv(digit.tab, file = paste0(output.path.012, chrom.name, "_fulsib_genotypes_", mother.name, "X", father.name, ".csv"), quote = FALSE, row.names = FALSE)
        cat("Successfully completed processing..... ", paste0(output.path.012, chrom.name, "_fulsib_genotypes_", mother.name, "X", father.name, ".csv"), "\n")
      } else {
        stop("Dimensions do not match between digit.tab and tmp.subset.var.tab")
      }
      
#     } # father.name
#   } # mother.name
# } # chrom.name
########################################################################################################################################################################################################
# END
      # example output file format for the fulsib family
      # CHROM, POS, REF, ALT, 4593, 2393, 26241, 26244, 26246, 26249, 26250, 26251, 26252, 26253, 26274, 26278, 26293, 26302, 26304
      # Chr05, 7196, A, C, 1, 2, 1, NA, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1
      # Chr05, 8867, A, G, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1
      # Chr05, 9407, A, C, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 0, 1
      # Chr05, 10957, G, A, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      # Chr05, 11004, A, G, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 2, 0, 1
      # Chr05, 14286, G, A, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1
      # Chr05, 15015, G, A, 0, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1
      # Chr05, 17807, T, A, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1
      # Chr05, 18034, C, G, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 2, 0, 1