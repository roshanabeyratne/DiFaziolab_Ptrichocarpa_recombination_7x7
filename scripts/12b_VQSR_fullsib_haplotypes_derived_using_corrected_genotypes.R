#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: July 03, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name, mother.name and father.name. 

# The input file format for the fulsib family genotypes are is follows; 
# CHROM,POS,REF,ALT,1863,6909,24736,24740,24745,24747,26631,26635
# Chr15,1143,CA,C,NA,NA,NA,NA,NA,NA,NA,NA
# Chr15,1202,CCCCAAACCCAAA,C,NA,NA,NA,NA,NA,NA,NA,NA
# Chr15,1241,A,AC,NA,NA,NA,NA,NA,0,NA,NA
# Chr15,1245,T,TA,NA,NA,NA,NA,NA,0,NA,NA
# Chr15,1525,A,AC,NA,0,NA,NA,NA,NA,NA,NA
# Chr15,1645,CTAAACCCT,C,NA,0,NA,NA,NA,NA,NA,NA
# Chr15,1669,AAACCCT,A,NA,NA,NA,0,NA,0,NA,NA
# Chr15,1709,C,CA,0,NA,NA,NA,NA,NA,NA,NA
# Chr15,1710,AAACCCT,A,0,NA,NA,NA,NA,NA,NA,NA

# The input files for focal parent corrected genotype is as follows; 
# CHROM,POS,REF,ALT,1909,CORRECTED_GENOTYPE
# Chr15,    1143,CA,C,NA,0
# Chr15,    1202,CCCCAAACCCAAA,C,NA,0
# Chr15,    1241,A,AC,NA,0
# Chr15,    1245,T,TA,NA,0
# Chr15,    1525,A,AC, 0,0
# Chr15,    1645,CTAAACCCT,C, 0,0
# Chr15,    1669,AAACCCT,A, 0,0
# Chr15,    1709,C,CA, 2,0
# Chr15,    1710,AAACCCT,A, 0,0


# PURPOSE:This code substitutes parental genotypes with the corrected version from an earlier script 11b. Then extracts full-sib family inforamation for a given focal parent from the SNP and INDEL 
# VQSR-dataset for each chromosome. It uses the trio info to identify the focal parental allle contribution in offspring in each of it's full-sib families. 
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
chrom.name.vec <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

input.path.fs.fam.geno.VQSR <- "./genotypes/VQSR/"
input.path.corrected.geno <- "./genotypes/parents/VQSR/"

output.path.012 <- "./genotypes/012/VQSR/"
output.path.haplotypes <- "./haplotypes/by_fulsib_family/VQSR/"
########################################################################################################################################################################################################
# for (chrom.name in chrom.name.vec[1:19]) { # chrom.name <- chrom.name.vec[2]
  # for (mother.name in par.name.vec[1:7]) { # mother.name <- "1863"
    # for (father.name in par.name.vec[8:14]) { # father.name <- "2572"
      
      VQSR.all.geno <- read.csv(file = paste0(input.path.fs.fam.geno.VQSR, chrom.name, "_fulsib_genotypes_", mother.name, "X", father.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)
      colnames(VQSR.all.geno) <- colnames(VQSR.all.geno) %>% str_replace("X", "")
      
      corrected.geno.mother <- read.csv(file = paste0(input.path.corrected.geno, chrom.name, "_corrected_genotypes_", mother.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)
      colnames(corrected.geno.mother) <- colnames(corrected.geno.mother) %>% str_replace("X", "")
      
      corrected.geno.father <- read.csv(file = paste0(input.path.corrected.geno, chrom.name, "_corrected_genotypes_", father.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)
      colnames(corrected.geno.father) <- colnames(corrected.geno.father) %>% str_replace("X", "")
      
      if (all(corrected.geno.mother$POS == corrected.geno.father$POS)) {
        
        included.idx <- which(VQSR.all.geno$POS %in% corrected.geno.mother$POS)
        if (nrow(corrected.geno.mother) == length(included.idx)) {
          digit.tab <- VQSR.all.geno[included.idx, ]
          digit.tab[,mother.name] <- corrected.geno.mother$CORRECTED_GENOTYPE
          digit.tab[,father.name] <- corrected.geno.father$CORRECTED_GENOTYPE
          write.csv(digit.tab, file = paste0(output.path.012, chrom.name, "_fulsib_genotypes_", mother.name, "X", father.name, ".csv"), quote = FALSE, row.names = FALSE)  
        } else {
          stop("error in number of markers subsetted to digit.tab!")
        }
      } else {
        stop("parent corrected genotypes do not have equal number of markers")
      }
      
      
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
        
        # Sys.time()
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
        # Sys.time()
        
        if (nrow(tmp.tab) < ncol(tmp.tab)) { # this IF block takes care of complications in function due to one offspring for a given fulsib family
          tmp.tab <- as.data.frame(t(tmp.tab), stringsAsFactors = FALSE)
        }
        
        haplo.tab <- cbind(digit.tab[, 1:6], tmp.tab)
        colnames(haplo.tab) <- colnames(digit.tab)
        write.csv(haplo.tab, file = paste0(output.path.haplotypes, chrom.name, "_fulsib_haplotypes_", mother.name, "X", father.name, "_", focal.par, ".csv"), quote = FALSE, row.names = FALSE)
        cat("Successfully completed processing..... ", chrom.name, "\t", paste0(output.path.haplotypes, chrom.name, "_fulsib_haplotypes_", mother.name, "X", father.name, "_", focal.par, ".csv"), "\n")
      }
    # }
  # }
# }
########################################################################################################################################################################################################
# END
# following is the output format that identifies the haplotype contribution by the focal parents to their fulsib family based on trio information; 
# CHROM,POS,REF,ALT,2283,2365,25954,25955,25956,25957,25962,25964,25966,25967,25969,25971,25973,25975,25979,25983,25987,25988,25990,25991
# Chr02,1032,CCTAAACT,C,0,0,0,NA,NA,0,NA,NA,0,0,0,0,NA,NA,NA,0,0,NA,NA,0
# Chr02,1062,CCTAAAAA,C,0,0,0,NA,NA,0,0,NA,0,0,0,0,NA,0,0,0,0,0,NA,0
# Chr02,1089,G,T,0,0,0,0,NA,0,0,0,0,0,0,0,NA,0,0,0,0,0,0,0
# Chr02,1090,ACCCTAAACACCTAAAACCCTAAC,A,0,0,0,0,NA,0,0,0,0,0,0,0,NA,0,0,0,0,0,0,0