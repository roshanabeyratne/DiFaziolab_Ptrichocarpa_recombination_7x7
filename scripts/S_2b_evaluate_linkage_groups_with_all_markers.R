#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: August 5, 2020
# USAGE: Input files include smoothed-out bakground haplotype contribution from the focal parent. Some individuals were removed from relevant parents due to high missingness and switches with script
# 19b. Therefore the path to such curated files is "./resolved_haps/by_parent/curated/". If a parent is not corrected such as described before, the focal parent hap contribution is at 
# "./resolved_haps/by_parent/"

# Following is the input file format used that has focal parent haplotype identification (hap_1 or hap_2) for all of its offspring;  
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1

# PURPOSE: Remove the offspring identified as having high missingness and higher number of focal parent haplotype switches. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
library(onemap)
########################################################################################################################################################################################################
#  checking bash idxs created.
# parentName=( 1863 1909 1950 2048 2066 2283 4593 2365 2393 2515 2572 2683 6909 7073 )
# chromNames=( Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 )
# for i in $(seq 1 266); do Nchrom=19; Token=$((${i}-1)); parIdx=$((${Token}/${Nchrom})); chromIdx=$((${Token}%${Nchrom}); echo ${i} ${chromNames[${chromIdx}]} ${parentName[${parIdx}]}; done
# args <- commandArgs(TRUE)
# chrom.name <- args[1]
# par.name <- args[2]
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

hap.input.path <- "./resolved_haps/by_parent/"
onemap.format.input.path <- "./onemap/bin_markers/input_files/all_markers/"
genetic.map.output.path <- "./onemap/bin_markers/genetic_maps/all_markers/"
bin.info.output.path <- "./onemap/bin_markers/bin_info/all_markers/"
########################################################################################################################################################################################################
for (par.name  in par.name.vec[1]) {  # par.name <- "1863"
  haplo.tab.combo <- NULL
    for (chrom.name in chrom.name.vec[1:19]) { # chrom.name <- chrom.name.vec[1]
  
      haplo.tab <- read.csv(file = paste0(hap.input.path, chrom.name, "_resolved_background_haplotypes_all_offspring_", par.name,".csv"), header = TRUE, stringsAsFactors = FALSE)
      colnames(haplo.tab) <- colnames(haplo.tab) %>% str_replace("^X", "")
      
      # combine the haplo.tabs of all the chromosomes into one master table    
      haplo.tab.combo <- rbind(haplo.tab.combo, haplo.tab)
      cat(chrom.name, "\n")
    }
  
  # coding for Onemap backcross configuration. hap_1 is coded as 'a' and hap_2 as 'ab'. zeros are coded as missing data as these haplotypes weren't determined with earlier scripts. 
  onemap.format <- apply(haplo.tab.combo[-(1:2)], c(1,2), function(x){ 
    x <- as.integer(x)
    if (is.na(x)) {
      stop("focal parent haplotype cannot be NA!")
    } else {
      if (x == 0) {
        genotype <-  "-"
      } else if (x == 1) {
        genotype = "a"
      } else if (x == 2) {
        genotype <-  "ab"
      } else {
        stop("genotypes can be only coded with integers 0, 1 & 2")
      }
    }
    return(genotype)
  })
  onemap.format <- as.data.frame(onemap.format, stringsAsFactors = FALSE)
  
  # produce marker names all markers
  markers.names <- as.character(paste0("*", haplo.tab.combo$CHROM, "_", haplo.tab.combo$POS))
  
  # produce individual names. This is not the same as colnames as Onemap requires them to be I1, I2, I3 ... 
  ind.names <- paste0("I", 1:ncol(onemap.format))
  
  # compile metadata lines for the input file
  first.line <- noquote(paste("data", "type", "f2", "backcross", sep = " "))
  ind.names.line <- noquote(paste(ind.names,collapse = " "))
  metadata <- noquote(paste(length(ind.names), nrow(onemap.format), "0","0","0", sep = " "))
  
  # create a filehandle for the output and write elements iteratively
  onemap.file = paste0(onemap.format.input.path, par.name, "_marker_binning_all_contigs", ".raw")
  unlink(onemap.file)
  rm(FILE1)
  FILE1 <- file(onemap.file,"w")
  write(first.line, FILE1, append = TRUE)
  write(metadata, FILE1, append = TRUE)
  write(ind.names.line, FILE1, append = TRUE)
  
  for (i in 1:nrow(onemap.format)) { #i<- 1
    line <- noquote(paste(markers.names[i],"A.H",paste(onemap.format[i,],collapse = " "),sep = " "))
    write(line, FILE1, append = TRUE)
  }
  close(FILE1)
  
  
  # Read-in the Onemap input files
  db.raw <- read_onemap(inputfile = onemap.file)
  save(db.raw, file=paste0(onemap.format.input.path, par.name, "_all_mrkrs_db_raw.RData"))
  bins.indb.raw <- find_bins(db.raw, exact = TRUE)
  save(bins.indb.raw, file=paste0(onemap.format.input.path, par.name, "_all_mrkrs_bins_indb_raw.RData"))

  # bin markers that have exactly the same information so that compute time reduces
  if (length(bins.indb.raw$bins) < bins.indb.raw$info$n.mar) {
    binned.db.raw <- create_data_bins(db.raw, bins.indb.raw)
  } else {
    binned.db.raw <- db.raw
  }
  
  
  # writing the binned marker information separately so that it can be reconciled with the other focal parents.
  # lapply(bins.indb.raw[[1]], function(x) write.table( data.frame(x), paste0(bin.info.output.path, par.name, "_binned_markers.txt"), append = T, sep = ',' ))
  
  # calculate twopoint recombination fractions
  LOD.sug <- suggest_lod(binned.db.raw)
  twopts <- rf_2pts(binned.db.raw, LOD = LOD.sug, max.rf = 0.5)
  save(twopts, file=paste0(onemap.format.input.path, par.name, "_all_mrkrs_twopts.RData"))
  
  markers.for.binned.db.raw <- NULL
  markers.for.binned.db.raw <- make_seq(twopts, "all")
  
  linkage.groups <- group(markers.for.binned.db.raw, LOD = LOD.sug, max.rf = 0.5)
  save(linkage.groups, file=paste0(onemap.format.input.path, par.name, "_all_mrkrs_linkage_groups.RData"))
  # lapply(linkage.groups, function(x) write.table( data.frame(x), paste0(genetic.map.output.path, par.name, "_linkage_grouo_info.txt"), append = T, sep = ',' ))
}


# Checking whether using a lower value for 'max.rf' argument changes the number of LGs identified or markers clustered within LGs
# db.raw <- read_onemap(inputfile = onemap.file)
load('/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/input_files/all_markers/1863_all_mrkrs_db_raw.RData')
# bins.indb.raw <- find_bins(db.raw, exact = TRUE)
load('/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/input_files/all_markers/1863_all_mrkrs_bins_indb_raw.RData') 
# bin markers that have exactly the same information so that compute time reduces
if (length(bins.indb.raw$bins) < bins.indb.raw$info$n.mar) {
  binned.db.raw <- create_data_bins(db.raw, bins.indb.raw)
} else {
  binned.db.raw <- db.raw
}
# calculate twopoint recombination fractions
LOD.sug <- suggest_lod(binned.db.raw)
twopts <- rf_2pts(binned.db.raw, LOD = LOD.sug, max.rf = 0.4)
markers.for.binned.db.raw <- NULL
markers.for.binned.db.raw <- make_seq(twopts, "all")
linkage.groups <- group(markers.for.binned.db.raw, LOD = LOD.sug, max.rf = 0.4)

for (i in 1:19) { # i=1
  chrom.name <- chrom.name.vec[i]
  tmp_LG <- make_seq(linkage.groups, arg = i)
  idxs <- tmp_LG$seq.num
  all_markrs <- colnames(as.data.frame((binned.db.raw$geno), stringsAsFactors=F)[1,])
  all_markrs_chrom <- all_markrs %>% str_replace("_.*", '')
  cat(chrom.name, ":", all(all_markrs_chrom[idxs]==chrom.name), "\n")
}
##############################################################################################################################################################
# END