#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: March 16, 2020, Modified: June 04, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name.This is the input file format of haplotype contribution from the focal parent

# CHROM,POS,REF,ALT,2066,25376,25380,25382,25383,25384,25386,25902,26355,26357,26418,26422,26423,26424,26427,26439,26443,26446,26447,25820,25821,25822,25823,25825,25826,25828,25829,25830,25831,25832,25836,25838,25839,25840,25841,25842,25341,25343,25346,25363,25365,25369,25370,26366,26369,26373,26378,26380,26385,26386,26388,26392,26395,25843,25846,25849,25854,25858,25861,25862,25869,25870,25871,25872,25875,25878,25879,25882,25889,25891,25895,25912,25919,25920,25924,25943,25947,26310,26313,26317,26318,26329,26331,26333,26336,26344,25918,25927,26256,26259,26260,26261,26262,26263,26264,26265,26266,26267,26268,26269,26271,26348,26349,26350,26352,26353,24977,24981,24983,24986,24987,24988,26679,26693,26704,26708,26725,26728,26736,26740,26745,26752,26753,26754,26759,26761,26765,26768,26771,26778
# Chr15,8393,C,T,1,0,NA,NA,1,0,NA,1,1,0,NA,1,NA,NA,1,0,0,0,1,0,0,0,1,0,0,NA,1,1,1,NA,1,0,0,1,1,NA,NA,1,NA,0,NA,0,NA,NA,NA,1,0,1,0,NA,1,NA,NA,0,1,NA,NA,0,0,NA,0,0,1,0,1,0,1,NA,0,1,1,1,0,1,1,1,0,1,0,0,1,1,0,1,NA,0,0,0,1,1,1,1,0,1,0,0,1,1,0,1,1,1,1,1,1,1,1,0,NA,1,1,NA,1,0,1,NA,1,NA,1,0,0,1,NA,1,NA,1,1,1,1,NA
# Chr15,8421,C,T,1,0,NA,NA,1,0,NA,1,1,NA,NA,1,NA,NA,1,0,0,0,1,0,0,0,1,0,0,NA,1,1,1,NA,1,0,0,1,1,NA,NA,1,NA,0,NA,0,NA,NA,NA,1,0,1,0,NA,1,NA,NA,0,1,NA,NA,0,0,NA,0,0,1,0,1,0,1,NA,0,1,1,1,0,1,1,1,0,1,0,0,1,1,0,1,NA,0,0,0,1,1,1,1,NA,1,0,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,NA,1,0,1,NA,1,NA,1,0,0,1,NA,1,NA,1,1,1,1,NA

# PURPOSE: Phase heterozygous parental markers based on twopt method in Onemap. Markers are coded in both coupled and reulsion phases. Description needs to be completed fo this section
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
input.path.combined.haplo <- "./haplotypes/by_parent/"
# input.path.combined.haplo <- "../../../../../Desktop/"
########################################################################################################################################################################################################
total.NA.count <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  for (par.name in par.name.vec) { # par.name <- par.name.vec[2] 
    
    parent.haplo.tab <- read.csv(file = paste0(input.path.combined.haplo, chrom.name, "_haplotypes_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    colnames(parent.haplo.tab) <- colnames(parent.haplo.tab) %>% str_replace(".GT", "") %>% str_replace("^X", "")
    
    # Selecting markers for which focal parent is heterozygous
    par.het.haplo.tab <- parent.haplo.tab[parent.haplo.tab[,par.name] == 1, ] 
    par.het.haplo.tab <- par.het.haplo.tab[!is.na(par.het.haplo.tab$POS), ] # 20841 markers for par.name = "2283"
    
    # count the amount of missing data for loci
    NA.count.per.loci <- apply(par.het.haplo.tab[,6:ncol(par.het.haplo.tab)], 1, function(x){ # offspring information start from column-6 onwards.
      count.Nas <- sum(is.na(x)) / length(x)
      return(count.Nas)
    })
    
    total.NA.count <- c(total.NA.count, NA.count.per.loci)
  }
}

hist(total.NA.count, breaks = 1000)
median(total.NA.count)
mean(total.NA.count)
quantile(total.NA.count, c(.05, 0.75, 0.9))
# 5%        75%        90% 
# 0.08527132 0.26470588 0.31775701 
# 
# i have decided to go with a cut-off of 25% (This means that a missingness rate of 25% would be allowed for loci)
########################################################################################################################################################################################################
#END
########################################################################################################################################################################################################