#!/usr/bin/R
setwd("/gpfs/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7")
# getwd()

# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: Aug 26, 2021
# USAGE: Input files include smoothed-out bakground haplotype contribution for 10kb genomic windows for all parents
# Following is the input file format used that has focal parent haplotype identification (hap_1 or hap_2) for all of its offspring;
# start_pos end_pos 1863_24708 1863_24709 1863_24710 1863_24711 1863_24712 1863_24713 1863_24714 
# 1         1    3968          B          B          B          A          A          B         
# 2      3969   13968          B          B          B          A          A          B          
# 3     13969   23968          B          B          B          A          A          B          
# 4     23969   33968          B          B          B          A          A          B          
# PURPOSE: Select loci with minimum unresolved haplotypes and produce genetic maps with marker ordering using Onemap-R package

rm(list = ls())
library(stringr)
library(onemap)
source(file = "./scripts_lists/Automated_Orderseq_function.R")
# source(file = "./Stettler_14_linkage_mapping/new_scripts/Automated_Orderseq_function.R")

#  checking bash idxs created.
# parentName=( 1863 1909 1950 2048 2066 2283 4593 2365 2393 2515 2572 2683 6909 7073 )
# chromNames=( Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 )
# for i in $(seq 1 266); do Nchrom=19; Token=$((${i}-1)); parIdx=$((${Token}/${Nchrom})); chromIdx=$((${Token}%${Nchrom}); echo ${i} ${chromNames[${chromIdx}]} ${parentName[${parIdx}]}; done
args <- commandArgs(TRUE)
chrom.name <- args[1]

# Set-up variables
chrom_name_vec <- paste0("Chr",sprintf('%02d',1:19))
par_name_vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", 
                  "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
onemap_format_input_path <- "./onemap/consensus_maps/onemap_inputs/"
# onemap_format_input_path <- "./Stettler_14_linkage_mapping/consensus_map/onemap_inputs/"
genetic.map.output.path <- "./onemap/consensus_maps/genetic_maps/"

# Obtain offspring haplotype configuration for all parents
# for (chrom.name in chrom_name_vec) {
  # chrom.name <- chrom_name_vec[19]
  hap_path <- "./onemap/consensus_maps/interval_haplotypes_combined/"
  # hap_path <- "./quantitative_genetic_analysis_7X7/qtl_analysis/interval_haplotypes_combined/"
  chrom_haps <- read.csv(file=paste0(hap_path, chrom.name, "_offspring_intvl_haplotypes_combined.csv"),
                         header=T, stringsAsFactors=F) 
  colnames(chrom_haps) <- colnames(chrom_haps) %>% str_replace("^X", "")
  # cat(dim(chrom_haps), "\n")
# }
# NOTE: All chromosomes do not have the same number of offspring. This is because some of the offspring were dopped
# from the analysis based on high number of switches between haplotypes to avoid map inflation.

# tmp_vec <- as.matrix(chrom_haps[,3:1660])
# unique(as.vector(tmp_vec)) # Only A, B and NA are possible values

# # NOTE:  This section is for changing the haplotype coding for just one parent (GW-1909).
# # This is to check whether arbitrary coding of haplotype makes a difference in genetic map size.
# # Combining just two parents to test whether arbitrary haplotype coding makes a difference
# chrom_haps <- chrom_haps[,1:220] # subsetting to parents GW-1863 and GW-1909
# tmp_chrom_haps <- chrom_haps[,103:220]
# tmp_chrom_haps <- apply(tmp_chrom_haps, c(1,2), function(x){
#   if (is.null(x)) {
#     stop("focal parent haplotype cannot be empty!")
#   } else {
#     if (is.na(x)) {
#       genotype = NA
#     } else if (x == 'A') {
#       genotype = "B"
#     } else if (x == 'B') {
#       genotype <-  "A"
#     } else {
#       stop("genotypes can be only coded with A, B & NA")
#     }
#   }
#   return(genotype)
# })
# tmp_chrom_haps <- as.data.frame(tmp_chrom_haps, stringsAsFactors = FALSE)
# chrom_haps[,103:220] <- tmp_chrom_haps

# Re-code the haplotypes to 'onemap' format
# coding for Onemap backcross configuration. hap_A is coded as 'a' and hap_2 as 'ab'. NAs are coded as '-'
onemap.format <- apply(chrom_haps[-(1:2)], c(1,2), function(x){
  if (is.null(x)) {
    stop("focal parent haplotype cannot be empty!")
  } else {
    if (is.na(x)) {
      genotype <-  "-"
    } else if (x == 'A') {
      genotype = "a"
    } else if (x == 'B') {
      genotype <-  "ab"
    } else {
      stop("genotypes can be only coded with A, B & NA")
    }
  }
  return(genotype)
})
onemap.format <- as.data.frame(onemap.format, stringsAsFactors = FALSE)

# produce marker names all markers
markers.names <- as.character(paste0("*", chrom.name, '_', round((chrom_haps$start_pos+chrom_haps$end_pos)/2, digits=0)))

# produce individual names. This is not the same as colnames as Onemap requires them to be I1, I2, I3 ...
ind.names <- paste0("I", 1:ncol(onemap.format))

# compile metadata lines for the input file
first.line <- noquote(paste("data", "type", "f2", "backcross", sep = " "))
ind.names.line <- noquote(paste(ind.names,collapse = " "))
metadata <- noquote(paste(length(ind.names), nrow(onemap.format), "0","0","0", sep = " "))

# create a filehandle for the output and write elements iteratively
onemap.file = paste0(onemap_format_input_path, chrom.name, "_marker_binning.raw")
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
bins.indb.raw <- find_bins(db.raw, exact = TRUE)

# bin markers that have exactly the same information so that compute time reduces
if (length(bins.indb.raw$bins) < bins.indb.raw$info$n.mar) {
  binned.db.raw <- create_data_bins(db.raw, bins.indb.raw)
} else {
  binned.db.raw <- db.raw
}

# calculate twopoint recombination fractions
LOD.sug <- suggest_lod(binned.db.raw) # can use LOD = 10 but get the same result
twopts <- rf_2pts(binned.db.raw, LOD = 10, max.rf = 0.5)

markers.for.binned.db.raw <- NULL
markers.for.binned.db.raw <- make_seq(twopts, "all")

linkage.groups <- group(markers.for.binned.db.raw, LOD = LOD.sug, max.rf = 0.5)
LG1 <- make_seq(linkage.groups, arg = 1)
# set the mapping funtion that converts recombination fraction to distance between markers
set_map_fun(type = "kosambi")

LG1.rec <- modified_order_seq(input.seq = LG1, n.init = 7, subset.search = "twopt", twopt.alg = "rec", THRES = 3, touchdown = TRUE)
save(LG1.rec, file = paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/consensus_maps/genetic_maps/", chrom.name, "_genetic_map.RData"))

# making a sequence with the ordered markers
LG.rec.map <- make_seq(LG1.rec, "safe")
ripple_seq(LG.rec.map, ws = 4, LOD = 3, tol = 10E-5)

# make genetic maps
LG1_rec.2 <- LG1.rec
LG1.rec <- LG1_rec.2$ord

write_map(LG1.rec, file.out=paste0(genetic.map.output.path, chrom.name, "_", "_genetic.map"))
cat(paste0(genetic.map.output.path, chrom.name, "_", "_genetic.map"), "Successfully completed!")

### END
