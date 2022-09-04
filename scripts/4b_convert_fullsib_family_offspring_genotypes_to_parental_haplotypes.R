#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: March 23, 2020; Modified - June 12, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name, mother.name and father.name. input files include truthdata-set and the corrcted parent genotypes.

# SNP truthdata-set format;
# CHROM	POS	REF	ALT	1863.GT	1909.GT	1950.GT	2048.GT	2066.GT	2283.GT	2365.GT	2393.GT	24708.GT	24709.GT	24710.GT	24711.GT	24712.GT
# Chr01	7940	T	 C	C/C	    C/C	    C/C	    T/C	    C/C	    T/C	    T/T	    C/C	    T/T	      T/T	      T/T       T/T       T/T

# The corrected genotype format for the parents; 
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

# PURPOSE:This code substitutes parental genotypes with the corrected version from an earlier script. Then extracts full-sib family inforamation for a given focal parent from the SNP truth-dataset
# for each chromosome. It first converts the GATK VariantsToTable format to 012 format and then uses the trio info to identify the focal parental allle contribution in offspring in each of it's 
# full-sib families. 
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
par.geno.input.path <- "./genotypes/parents/"

output.path.012 <- "./genotypes/012/"
output.path.haplotypes <- "./haplotypes/by_fulsib_family/"
########################################################################################################################################################################################################

# for (chrom.name in chrom.name.vec[1:19]) { # chrom.name <- chrom.name.vec[2]
  
  # reading the SNP truth-dataset in GATK VatriantToTable format for each chromosome
  snp.file.name <- paste0("./truth_data/7x7_Stet14_SNP_truthdataset_", chrom.name, ".table") 
  # snp.file.name <- paste0("./truthdataset/7x7_Stet14_SNP_truthdataset_", chrom.name, ".table") 
  tmp.variant.table <- read.table(file = snp.file.name, stringsAsFactors = FALSE, header = TRUE)
  # tmp.variant.table <- read.table(file = paste0("./lists/7x7_Stet14_bi_SNP_VQSR_", chrom.name, ".table"), stringsAsFactors = FALSE, header = TRUE) # reading the SNP VQSR for each chromosome
  # tmp.variant.table <- read.table(file = paste0("./lists/SNP_truthdata_set_for_2283_Chr02.table"), stringsAsFactors = FALSE, header = TRUE) # just for checking 2283 pseudo testcross markers
  colnames(tmp.variant.table) <- colnames(tmp.variant.table) %>% str_replace(".GT", "") %>% str_replace("^X", "") # removing GATK suffixes and 'X' in front of the clone-names when they are read to R
  
  
  
  # substitute the corrected parental genotypes to the table.
  for (par.name in par.name.vec) { # par.name <- "2283"
    tmp.corrected.geno <- read.csv(file = paste0(par.geno.input.path, chrom.name, "_corrected_genotypes_", par.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)
    if (all(tmp.corrected.geno$POS == tmp.variant.table$POS)) {
      par.name.idx <- which(colnames(tmp.variant.table) == par.name)
      tmp.variant.table[, par.name.idx] <- tmp.corrected.geno$CORRECTED_GENOTYPE
    } else {
      stop(paste0("Corrected genotype position IDs do not match with the SNP truthset for parent: ", par.name))
    }
  }
  rm(tmp.corrected.geno)
  
  
  # for (mother.name in par.name.vec[1:7]) { # mother.name <- "1863"
    # for (father.name in par.name.vec[8:14]) { # father.name <- "2572"
      
      
      # Identify offspring clone names/individuals for the full-sib family given the most recent parentage assignement and then subset the truth-dataset for the full-sib family.  trio information.
      offspring.idx <- which((parantage_v4.tab$matched_mother == mother.name) & (parantage_v4.tab$matched_father == father.name))
      tmp.offspring.names <- as.character(parantage_v4.tab$OldClone[offspring.idx]) # old clone names are used here
      offspring.names <- sort(colnames(tmp.variant.table)[which(colnames(tmp.variant.table) %in% tmp.offspring.names)])
      tmp.subset.var.tab.par.cols <- tmp.variant.table[, c(mother.name, father.name)]
      tmp.subset.var.tab <- tmp.variant.table[, c(colnames(tmp.variant.table)[1:4], offspring.names)]

      
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
      
      digit.tab <- cbind(tmp.subset.var.tab[, 1:4], tmp.subset.var.tab.par.cols, tmp.tab)
      colnames(digit.tab) <- c(colnames(tmp.subset.var.tab)[1:4], colnames(tmp.subset.var.tab.par.cols), colnames(tmp.subset.var.tab)[-(1:4)])
      write.csv(digit.tab, file = paste0(output.path.012, chrom.name, "_fulsib_genotypes_", mother.name, "X", father.name, ".csv"), quote = FALSE, row.names = FALSE)
      
      
      # Followng function identifies each focal parent allele contribution in offspring. The hap.code matrix below identified mendelian violations as NA as well as codes Het*Het crosses that result
      # Het offspring as offspring as well.
      for (focal.par in c(mother.name, father.name)) { # focal.par <- mother.name
        if (focal.par == mother.name) {
          focal.col <- 5
          other.col <- 6
        } else {
          focal.col <- 6
          other.col <- 5
        }
        
        tmp.tab <- as.data.frame(t(apply(digit.tab, 1, function(x, focal.idx = focal.col, other.idx = other.col){
          hap.code <- c("0-NA-NA-NA-NA", "0-0-NA-NA", "NA-0-NA-NA", "NA-NA-NA-NA",
                        "0-1-NA-NA", "0-NA-1-NA", "NA-0-1-NA", "NA-NA-NA-NA",
                        "NA-1-NA-NA", "NA-1-1-NA", "NA-NA-1-NA", "NA-NA-NA-NA",
                        "NA-NA-NA-NA", "NA-NA-NA-NA", "NA-NA-NA-NA", "NA-NA-NA-NA")
          hap.key.df <- as.data.frame(matrix(hap.code, nrow = 4, ncol = 4, byrow = TRUE), stringsAsFactors = FALSE)
          
          
          rows.to.process <- as.character(x)
          focal.geno <- as.numeric(rows.to.process[focal.idx])
          other.geno <- as.numeric(rows.to.process[other.idx])
          if (is.na(focal.geno)) {
            row.idx <- 4 # row index of the hap.key.df 
          } else {
            row.idx <- focal.geno + 1
          }
          if (is.na(other.geno)) {
            col.idx <- 4 # column index of the hap.key.df 
          } else {
            col.idx <- other.geno + 1
          }
          
          offspring.geno <- as.numeric(rows.to.process[7:length(x)])
          template.vec <- c(0, 1, 2, NA)
          hap.string <- as.numeric(unlist(str_split(hap.key.df[row.idx, col.idx], "-")))
          
          bin.vec <- as.vector(sapply(offspring.geno, function(y, template = template.vec, conversion = hap.string) {
            genotype <- y
            binary <- conversion[which(template %in% genotype)]
            return(binary)
          }))
          return(bin.vec)
        })), stringsAsFactors = FALSE)
        
        if (nrow(tmp.tab) < ncol(tmp.tab)) { # this IF block takes care of complications in function due to one offspring for a given fulsib family
          tmp.tab <- as.data.frame(t(tmp.tab), stringsAsFactors = FALSE)
        }
        
        haplo.tab <- cbind(digit.tab[, 1:6], tmp.tab)
        colnames(haplo.tab) <- colnames(digit.tab)
        write.csv(haplo.tab, file = paste0(output.path.haplotypes, chrom.name, "_fulsib_haplotypes_", mother.name, "X", father.name, "_", focal.par, ".csv"), quote = FALSE, row.names = FALSE)
        cat("Successfully completed processing..... ", chrom.name, "\t", paste0(output.path.haplotypes, chrom.name, "_fulsib_haplotypes_", mother.name, "X", father.name, "_", focal.par, ".csv"), "\n")
      }
      
    # } # mother.name bracket
  # } # father.name bracket
# } # chrom.name bracket
#######################################################################################################################################################################################################
# END
# This is the output format for the focal parent haplotypes
# CHROM,POS,REF,ALT,1863,2365,24708,24709,24710,24711,24712,24713,24714,24715,25426,25427,26550,26553
# Chr15,8393,C,T,1,1,NA,NA,NA,NA,NA,1,NA,0,0,NA,0,NA
# Chr15,8421,C,T,1,1,NA,NA,NA,NA,NA,1,NA,0,0,0,0,NA
# Chr15,8426,T,A,1,1,NA,NA,NA,NA,NA,1,0,0,0,NA,0,NA
# Chr15,8805,G,A,1,1,NA,0,NA,NA,NA,1,NA,0,0,0,0,NA
# Chr15,8951,A,T,0,1,0,0,0,0,0,NA,0,0,0,0,0,NA
# Chr15,9019,C,T,1,1,NA,0,NA,NA,NA,1,NA,NA,0,NA,0,NA
# Chr15,9178,T,C,1,1,NA,0,NA,0,0,1,NA,NA,0,NA,0,NA
# Chr15,10452,C,T,1,1,NA,NA,NA,NA,NA,1,NA,NA,0,0,0,NA
# Chr15,10571,A,G,1,1,NA,0,0,NA,NA,1,NA,NA,0,NA,0,NA