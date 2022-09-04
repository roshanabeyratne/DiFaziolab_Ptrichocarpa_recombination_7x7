#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: August 31, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name.This is the input file format of haplotype contribution from the focal parent
# PURPOSE: Track the process of marker elimination up to the process of producing the genetic maps for eah focal parent and chromosome at each stage
# Following are the stages identified along with scripts that produce these outputs;
# 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

input.path.combined.haplo <- "./haplotypes/by_parent/"
output.path.phased.markers <- "./onemap/phase_parents/phased_markers/"
hap.input.path <- "./resolved_haps/by_parent/"
curated.hap.input.path <- "./resolved_haps/by_parent/curated/"
genmap.stats.output.path <- "./onemap/bin_markers/genetic_maps/stats/"
########################################################################################################################################################################################################
mrkr.pipeline.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  for (par.name in par.name.vec) { # par.name <- par.name.vec[1]
    
    # all markers in truthdata
    parent.haplo.tab <- read.csv(file = paste0(input.path.combined.haplo, chrom.name, "_haplotypes_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    colnames(parent.haplo.tab) <- colnames(parent.haplo.tab) %>% str_replace(".GT", "") %>% str_replace("^X", "")
    truth.mrkr.count <- nrow(parent.haplo.tab)
    
    # filtering in all the het markers
    par.het.haplo.tab <- parent.haplo.tab[parent.haplo.tab[,par.name] == 1, ] # Selecting markers for which focal parent is heterozygous
    par.het.haplo.tab <- par.het.haplo.tab[!is.na(par.het.haplo.tab$POS), ] 
    het.truth.mrkr.count <- nrow(par.het.haplo.tab)
    hom.truth.mrkr.count <- truth.mrkr.count - het.truth.mrkr.count 
    
    # apply a hard-filter: markers for which >25% of the data is missing are omitted
    less.informative.marker.count <- apply(par.het.haplo.tab[,6:ncol(par.het.haplo.tab)], 1, function(x){ # offspring information start from column-6 onwards.
      count.Nas <- sum(is.na(x))
      return(count.Nas)
    })
    # hist(less.informative.marker.count, breaks = 100)
    less.informative.marker.idx <- which(less.informative.marker.count > (ncol(par.het.haplo.tab) * 0.25)) # removing markers that have more than 25% of offspring with missing information (NA)
    par.het.haplo.tab <- par.het.haplo.tab[-less.informative.marker.idx, ]
    het.post.missing.filter.mrkr.count <- nrow(par.het.haplo.tab)
    
    # markers remaining after the Onemap based marker clustering and phasing.
    tmp.df <- read.csv(file = paste0(output.path.phased.markers, chrom.name, "_phased_markers_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    OM.phased.mrkr.count <- nrow(tmp.df)
    
    
    # At this stage the filtered het markers that are phased and Hom markers out of the truthdata set are combined for focal parents. 
    if (file.exists(paste0(curated.hap.input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"))) {
      haplo.tab <- read.csv(file = paste0(curated.hap.input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
    } else {
      haplo.tab <- read.csv(file = paste0(hap.input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
    }
    colnames(haplo.tab) <- colnames(haplo.tab) %>% str_replace("^X", "")
    hap.identified.mrkr.count <- nrow(haplo.tab) # this count is approximately equal to (hom.truth.mrkr.count + OM.phased.mrkr.count)
    
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
    mis.adjusted.hap.identified.mrkr.count <- nrow(haplo.tab)
    
    # Re-loading the LG1.rec object
    load(paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/LG1_rec/", chrom.name, "_full_genetic_map_object_", par.name,".RData"))
    post.binning.mrkr.count <- length(LG1.rec$ord.all$seq.num)
    post.mapping.safe.mrkr.count <- length(LG1.rec$ord$seq.num)
    
    
    mrkr.pipeline.df <- rbind(mrkr.pipeline.df, c(chrom.name, par.name, truth.mrkr.count, het.truth.mrkr.count, hom.truth.mrkr.count, het.post.missing.filter.mrkr.count, OM.phased.mrkr.count, 
                                                  hap.identified.mrkr.count, mis.adjusted.hap.identified.mrkr.count, post.binning.mrkr.count, post.mapping.safe.mrkr.count))
  }
}
colnames(mrkr.pipeline.df) <- c("CHROM", "PARENT", "TRUTH", "TRUTH_HET", "TRUTH_HOM", "TRUTH_MISS_FILTD", "OM_PHASED", "COMBINED", "OM_INPUT_POSTMISS_FILT", "POST_BINNED", "POST_MAPPED")
write.csv(mrkr.pipeline.df, file = paste0(genmap.stats.output.path, "mrkr_pipeline_table.csv"), row.names = FALSE, quote = FALSE)
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################