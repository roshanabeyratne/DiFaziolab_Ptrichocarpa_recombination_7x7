#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: July 31, 2020
# USAGE: Input files include smoothed-out bakground haplotype contribution from the focal parent. Some individuals were removed from relevant parents due to high missingness and switches with script
# 19b. Therefore the path to such curated files is "./resolved_haps/by_parent/curated/". If a parent is not corrected such as described before, the focal parent hap contribution is at
# "./resolved_haps/by_parent/"

# Following is the input file format used that has focal parent haplotype identification (hap_1 or hap_2) for all of its offspring;
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1

# PURPOSE: Select loci with minimum unresolved haplotypes and produce genetic maps with marker ordering using Onemap-R package
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
library(onemap)

source(file = "./scripts_lists/Automated_Orderseq_function.R")
########################################################################################################################################################################################################
#  checking bash idxs created.
# parentName=( 1863 1909 1950 2048 2066 2283 4593 2365 2393 2515 2572 2683 6909 7073 )
# chromNames=( Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 )
# for i in $(seq 1 266); do Nchrom=19; Token=$((${i}-1)); parIdx=$((${Token}/${Nchrom})); chromIdx=$((${Token}%${Nchrom}); echo ${i} ${chromNames[${chromIdx}]} ${parentName[${parIdx}]}; done
args <- commandArgs(TRUE)
chrom.name <- args[1]
par.name <- args[2]
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

hap.input.path <- "./resolved_haps/by_parent/"
curated.hap.input.path <- "./resolved_haps/by_parent/curated/"
onemap.format.input.path <- "./onemap/bin_markers/input_files/"

bin.info.output.path <- "./onemap/bin_markers/bin_info/"
rf.graph.output.path <- "./onemap/bin_markers/rf_graphs/"
genetic.map.output.path <- "./onemap/bin_markers/genetic_maps/"
########################################################################################################################################################################################################
# for (chrom.name in chrom.name.vec[1]) { # chrom.name <- chrom.name.vec[1]
#   for (par.name in par.name.vec[1]) { # par.name <- "2365"

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
    

    # coding for Onemap backcross configuration. hap_1 is coded as 'a' and hap_2 as 'ab'. zeros are coded as missing data as these haplotypes weren't determined with earlier scripts.
    onemap.format <- apply(haplo.tab[-(1:2)], c(1,2), function(x){
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
    markers.names <- as.character(paste0("*", haplo.tab$CHROM, "_", haplo.tab$POS))

    # produce individual names. This is not the same as colnames as Onemap requires them to be I1, I2, I3 ...
    ind.names <- paste0("I", 1:ncol(onemap.format))

    # compile metadata lines for the input file
    first.line <- noquote(paste("data", "type", "f2", "backcross", sep = " "))
    ind.names.line <- noquote(paste(ind.names,collapse = " "))
    metadata <- noquote(paste(length(ind.names), nrow(onemap.format), "0","0","0", sep = " "))

    # create a filehandle for the output and write elements iteratively
    onemap.file = paste0(onemap.format.input.path, chrom.name, "_marker_binning_", par.name,".raw")
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
    save(binned.db.raw, file = paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/binned_db_raw/", chrom.name, "_binned_db_raw_object_", par.name,".RData"))
    # load(paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/binned_db_raw/", chrom.name, "_binned_db_raw_object_", par.name,".RData"))

    # writing the binned marker information separately so that it can be reconciled with the other focal parents.
    bin.file <- paste0(bin.info.output.path, chrom.name, "_", par.name, "_binned_markers.txt")
    unlink(bin.file)
    lapply(bins.indb.raw[[1]], function(x) write.table( data.frame(x), bin.file, append = T, sep = ',' ))

    # calculate twopoint recombination fractions
    LOD.sug <- suggest_lod(binned.db.raw) # can use LOD = 10 but get the same result
    twopts <- rf_2pts(binned.db.raw, LOD = 10, max.rf = 0.5)
    save(twopts, file = paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/twopts/", chrom.name, "_twopts_object_", par.name,".RData"))
    # load(paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/twopts/", chrom.name, "_twopts_object_", par.name,".RData"))
    
    markers.for.binned.db.raw <- NULL
    markers.for.binned.db.raw <- make_seq(twopts, "all")
    # marker_type(markers.for.binned.db.raw) # This displays the marker names and the sequence number in Onemap object file that was read-in.

    linkage.groups <- group(markers.for.binned.db.raw, LOD = LOD.sug, max.rf = 0.5)
    if (linkage.groups$n.groups > 1) {
      stop("Number of LGs are greater than one!")
    } else {
      LG1 <- make_seq(linkage.groups, arg = 1)
      # marker_type(LG1) # This displays the marker names and the sequence number in Onemap object file that was read-in.

      # set the mapping funtion that converts recombination fraction to distance between markers
      set_map_fun(type = "kosambi")

      # order markers based on the 'twopt' method in Onemap
      # LG1.rec <- order_seq(input.seq = LG1, n.init = 7, subset.search = "twopt", twopt.alg = "rec", THRES = 3, touchdown = TRUE)
      LG1.rec <- modified_order_seq(input.seq = LG1, n.init = 7, subset.search = "twopt", twopt.alg = "rec", THRES = 3, touchdown = TRUE)
      save(LG1.rec, file = paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/LG1_rec/", chrom.name, "_full_genetic_map_object_", par.name,".RData"))
      # load(paste0("/gpfs/group/difazio/populus/gatk-7x7-Stet14/onemap/bin_markers/genetic_maps/RData/LG1_rec/", chrom.name, "_full_genetic_map_object_", par.name,".RData"))

      # making a sequence with the ordered markers
      LG.rec.map <- make_seq(LG1.rec, "safe")
      ripple_seq(LG.rec.map, ws = 4, LOD = 3, tol = 10E-5)
      # temp.std.out <- capture.output(ripple_seq(LG.rec.map, ws = 4, LOD = 3, tol = 10E-5))
      # temp.map <- make_seq(twopts, c(LG1.rec$ord$seq.num[1:111], LG1.rec$ord$seq.num[c(113, 112, 114)]))
      # temp.map <- make_seq(twopts, sort(LG1.rec$ord$seq.num, decreasing = FALSE))
      #  temp.map2 <- map(temp.map)

      # Produce an rf /LOD graph which enables detection of mapping errors in markers
      pdf(file = paste0(rf.graph.output.path, chrom.name, "_rf_LOD_plot_", par.name,".pdf"))
      rf_graph_table(LG.rec.map, inter = FALSE)
      dev.off()

      # make genetic maps
      LG1_rec.2 <- LG1.rec
      LG1.rec <- LG1_rec.2$ord
      # LG1.rec <- LG1_rec.2$ord.all
      # LG1.rec <- input.seq2
      write_map(LG1.rec, file.out = paste0(genetic.map.output.path, chrom.name, "_", par.name, "_genetic.map"))
      cat(paste0(genetic.map.output.path, chrom.name, "_", par.name, "_genetic.map"), "Successfully completed!")
    }

#   }
# }
########################################################################################################################################################################################################
# END