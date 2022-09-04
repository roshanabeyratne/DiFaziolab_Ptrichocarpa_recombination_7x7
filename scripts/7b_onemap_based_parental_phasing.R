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
library(onemap)
source(file = "./scripts_lists/MakeOnemapInput_function.R")
# source(file = "./new_scripts/MakeOnemapInput_function.R")
########################################################################################################################################################################################################
args <- commandArgs(TRUE)
chrom.name <- args[1]
par.name <- args[2]
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
input.path.combined.haplo <- "./haplotypes/by_parent/"
# input.path.combined.haplo <- "../../../../../Desktop/"
tmp.input.path <- "./onemap/phase_parents/input_files/"
tmp.output.path <- "./onemap/phase_parents/output_files/"
output.path.phased.markers <- "./onemap/phase_parents/phased_markers/"
#######################################################################################################################################################################################################
# metadata.df <- NULL
# for (chrom.name in chrom.name.vec[1:19]) { # chrom.name <- chrom.name.vec[1]
  # for (par.name in par.name.vec[1:14]) { # par.name <- par.name.vec[1]
    
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
    ####################################################################################################################################################################################################
    
    # the maximum number of markers in the intensive set of markers to be clustered are decided first. The number of intesnsive groups is calculated afterwards. Onemap clustering is then carried-out
    # for each intensive group first. 
    max.intensive <- 200
    intensive.group.num <- round(nrow(par.het.haplo.tab) / max.intensive)
    intensive.group.df <- NULL
    critical.length <- 205
    for (b in 1:intensive.group.num) { # b <- 1
      seq.tmp <- seq(from = b, to = nrow(par.het.haplo.tab), by = intensive.group.num)
      # if (nrow(par.het.haplo.tab) %in% seq.tmp) {
      #   stop()
      # }
      seq.length <- length(seq.tmp)
      # # cat(seq.length, "\n")
      
      # following while loop assures that the seq.tmp vector has the max.ntensive number of elements. This is required in order to efficiently use cbind function.
      while (seq.length != critical.length) {
        seq.tmp <- c(seq.tmp[1], seq.tmp)
        seq.length <- length(seq.tmp)
      }
      intensive.group.df <- cbind(intensive.group.df, seq.tmp)
    } # intensive group bracket
    intensive.group.df <- as.data.frame(intensive.group.df)
    colnames(intensive.group.df) <- c(1:intensive.group.num)
    
    # #  this is a test code-block to see whether there are marker position duplications 
    # for (i in 1:ncol(intensive.group.df)) {
    #   for (j in 1:ncol(intensive.group.df)) {
    #     if (any(intensive.group.df[,i] %in% intensive.group.df[,j])) {
    #       if (i != j) {
    #         cat(i, " ", j, "\n")
    #       }
    #     }
    #   }
    # }
    # there were no marker duplications among intensive groups. 
    ####################################################################################################################################################################################################
    
    
    # A random subset of markers from each intensive group of markers are selected to do a defined number of reconciliatory groups.
    reconc.group.num <- 5
    sample.num.per.grp <- 10
    reconc.group.df <- NULL
    
    for (l in 1:reconc.group.num) { # i <- 1
      tmp.recon <- NULL
      for (i in 1:ncol(intensive.group.df)) { # i <- 5
        sample.tmp <- sample(unique(intensive.group.df[, i]), sample.num.per.grp, replace = FALSE)
        tmp.recon <- c(tmp.recon, sample.tmp)
      }
      tmp.recon <- sort(tmp.recon, decreasing = FALSE)
      reconc.group.df <- cbind(reconc.group.df, tmp.recon)
    }
    reconc.group.df <- as.data.frame(reconc.group.df)
    colnames(reconc.group.df) <- LETTERS[1:reconc.group.num]
    # which(reconc.group.df[,1] %in% reconc.group.df[,2]) # reconciliatory groups may share markers between each other
    
    # This is a test code-block to check whether each reconc group has any marker ositions duplicated
    # for (j in 1:ncol(reconc.group.df)) {
    #   if (length(reconc.group.df[, j]) == length(unique(reconc.group.df[,j]))) {
    #     cat(TRUE, "\n")
    #   }
    # }
    
    
    for (l in LETTERS[1:reconc.group.num]) { # l <- "A"
      # produce a temporary table to format Onemap input files. 
      tmp.tab <- par.het.haplo.tab[reconc.group.df[, l], ]
      
      # Creating input files in onemap format for markers in both coupled and repulsion. 
      onemap.file <- paste0(tmp.input.path, chrom.name, "_tmp_reconc_", par.name, "_", l, ".raw")
      MakeOnemapInput(allele.tab = tmp.tab, output.file.name = onemap.file)
      
      
      # phasing markers in the intensive group using Onemap:R-package. The markers in the same linkage group are the phased haplotype for the focal parent.
      db.raw <- read_onemap(inputfile = onemap.file)
      LOD.sug <- suggest_lod(db.raw)
      # LOD.sug <- 12
      twopts <- rf_2pts(db.raw, LOD = LOD.sug, max.rf = 0.5)
      markers.for.binned.db.raw <- NULL
      markers.for.binned.db.raw <- make_seq(twopts, "all")
      
      linkage.groups <- group(markers.for.binned.db.raw, LOD = LOD.sug, max.rf = 0.5)
      if (linkage.groups$n.groups < 2) {
        stop("only one linkage group produced")
      } else {
        lg.vec <- NULL
        for (lg in 1:linkage.groups$n.groups) {
          LG <-  make_seq(linkage.groups, arg = lg)
          lg.length <- length(LG$seq.num)
          lg.vec <- c(lg.vec, lg.length)
        }
        
        lg.1.idx <- order(lg.vec, decreasing = TRUE)[1] # obtain the largest and the second largest linkage groups 
        lg.2.idx <- order(lg.vec, decreasing = TRUE)[2]
        
        LG.1 <-  make_seq(linkage.groups, arg = lg.1.idx)
        LG.2 <-  make_seq(linkage.groups, arg = lg.2.idx)
        
        LG.1.marker.names <- as.character(linkage.groups$marnames[LG.1$seq.num])
        LG.2.marker.names <- linkage.groups$marnames[LG.2$seq.num]
        # length(LG.1.marker.names)
        
        write.table(LG.1.marker.names, file = paste0(tmp.output.path, chrom.name, "_reconc_", par.name, "_", l, "_", 1, ".txt"), row.names = FALSE, quote = FALSE)
        write.table(LG.2.marker.names, file = paste0(tmp.output.path, chrom.name, "_reconc_", par.name, "_", l, "_", 2, ".txt"), row.names = FALSE, quote = FALSE)
      }
    } # reconc group bracket
    ###################################################################################################################################################################################################
    
    # put together the reconc groups using intensive group-1. 
    # produce a temporary table to format Onemap input files for intensive group -1 
    tmp.tab <- par.het.haplo.tab[unique(intensive.group.df[, 1]), ]
    
    # Creating input files in onemap format for markers in both coupled and repulsion.
    onemap.file <- paste0(tmp.input.path, chrom.name, "_tmp_intensive_1", par.name, ".raw")
    MakeOnemapInput(allele.tab = tmp.tab, output.file.name = onemap.file)
    
    
    # phasing markers in the intensive group using Onemap:R-package. The markers in the same linkage group are the phased haplotype for the focal parent.
    db.raw <- read_onemap(inputfile = onemap.file)
    LOD.sug <- suggest_lod(db.raw)
    # LOD.sug <- 12
    twopts <- rf_2pts(db.raw, LOD = LOD.sug, max.rf = 0.5)
    markers.for.binned.db.raw <- NULL
    markers.for.binned.db.raw <- make_seq(twopts, "all")
    
    linkage.groups <- group(markers.for.binned.db.raw, LOD = LOD.sug, max.rf = 0.5)
    if (linkage.groups$n.groups < 2) {
      stop("only one linkage group produced in first intensive group")
    } else {
      lg.vec <- NULL
      for (lg in 1:linkage.groups$n.groups) {
        LG <-  make_seq(linkage.groups, arg = lg)
        lg.length <- length(LG$seq.num)
        lg.vec <- c(lg.vec, lg.length)
      }
      
      lg.1.idx <- order(lg.vec, decreasing = TRUE)[1] # obtain the largest and the second largest linkage groups 
      lg.2.idx <- order(lg.vec, decreasing = TRUE)[2]
      
      LG.1 <-  make_seq(linkage.groups, arg = lg.1.idx)
      LG.2 <-  make_seq(linkage.groups, arg = lg.2.idx)
      
      LG.1.marker.names <- linkage.groups$marnames[LG.1$seq.num]
      LG.2.marker.names <- linkage.groups$marnames[LG.2$seq.num]
      # length(LG.1.marker.names)
    }
    # tmp.set <- LG.2.marker.names %>% str_replace("_c", "") %>%  str_replace("_r", "") 
    # length(tmp.set) == length(unique(tmp.set))
    
    
    reconcile.collective <- NULL
    for (l in LETTERS[1:reconc.group.num]) { # l <- LETTERS[1]
      test.markers.1 <- read.table(file = paste0(tmp.output.path, chrom.name, "_reconc_", par.name, "_", l, "_", 1, ".txt"), stringsAsFactors = FALSE, header = TRUE)
      test.markers.1 <- apply(test.markers.1, 1, as.character)
      
      test.markers.2 <- read.table(file = paste0(tmp.output.path, chrom.name, "_reconc_", par.name, "_", l, "_", 2, ".txt"), stringsAsFactors = FALSE, header = TRUE)
      test.markers.2 <- apply(test.markers.2, 1, as.character)
      
      # cat(length(test.markers.1), "\n")
      
      if ((sum(LG.1.marker.names %in% test.markers.1) > 3) & (sum(LG.1.marker.names %in% test.markers.2) > 3)) {
        stop("intensive group markers in both reconcile groups!")
      } else if (sum(LG.1.marker.names %in% test.markers.1) > 2) {
        reconcile.collective <- c(reconcile.collective, test.markers.1)
      } else if (sum(LG.1.marker.names %in% test.markers.2) > 2) {
        reconcile.collective <- c(reconcile.collective, test.markers.2)
      } else {
        stop("intensive group markers are not in any of the reconcile linkage groups!")
      }
    }
    reconcile.collective.resolved <- unique(reconcile.collective)
    
    # following code block tests whether there are marker duplications within the reconcile.collective.resolved set
    # length(reconcile.collective.resolved)
    # tmp.set <- reconcile.collective.resolved %>% str_replace("_c", "") %>%  str_replace("_r", "") 
    # length(tmp.set) == length(unique(tmp.set))
    ###################################################################################################################################################################################################
    
    
    # Now I am putting the rest of the intensive groups together using the reconcile.collective.resolved
    intensive.collective <- NULL
    # for (n in c(1:ncol(intensive.group.df))[-c(30, 112)]) { # this is for "1950" and Chr01
    for (n in 1:ncol(intensive.group.df)) { # n <- 90
      
      tmp.tab <- par.het.haplo.tab[unique(intensive.group.df[, n]), ]
      
      # Creating input files in onemap format for markers in both coupled and repulsion.      
      onemap.file <- paste0(tmp.input.path, chrom.name, "_tmp_intensive_", n, "_", par.name, ".raw")
      MakeOnemapInput(allele.tab = tmp.tab, output.file.name = onemap.file)
      
      
      # phasing markers in the intensive group using Onemap:R-package. The markers in the same linkage group are the phased haplotype for the focal parent.
      db.raw <- read_onemap(inputfile = onemap.file)
      LOD.sug <- suggest_lod(db.raw)
      # LOD.sug <- 8
      twopts <- rf_2pts(db.raw, LOD = LOD.sug, max.rf = 0.5)
      markers.for.binned.db.raw <- NULL
      markers.for.binned.db.raw <- make_seq(twopts, "all")
      linkage.groups <- group(markers.for.binned.db.raw, LOD = LOD.sug, max.rf = 0.5)
      
      if (linkage.groups$n.groups < 2) {
        next(paste0("only one linkage group produced in the intensive group ", n))
      } else {
        lg.vec <- NULL
        for (lg in 1:linkage.groups$n.groups) {
          LG <-  make_seq(linkage.groups, arg = lg)
          lg.length <- length(LG$seq.num)
          lg.vec <- c(lg.vec, lg.length)
        }
        
        lg.1.idx <- order(lg.vec, decreasing = TRUE)[1] # obtain the largest and the second largest linkage groups 
        lg.2.idx <- order(lg.vec, decreasing = TRUE)[2]
        
        LG.1 <-  make_seq(linkage.groups, arg = lg.1.idx)
        LG.2 <-  make_seq(linkage.groups, arg = lg.2.idx)
        
        LG.1.marker.names <- linkage.groups$marnames[LG.1$seq.num]
        LG.2.marker.names <- linkage.groups$marnames[LG.2$seq.num]
        # length(LG.1.marker.names)
      }
      
      #  This section checks for markers in both confiigurations. 
      tmp.marker.names.1.dup  <- LG.1.marker.names %>% str_replace("_c", "") %>%  str_replace("_r", "") %>% duplicated
      tmp.marker.names.2.dup <- LG.2.marker.names %>% str_replace("_c", "") %>%  str_replace("_r", "") %>% duplicated
      if (any(tmp.marker.names.1.dup) | any(tmp.marker.names.2.dup)) {
        next("marker duplicated!")
      }
      
      
      if ((sum(LG.1.marker.names %in% reconcile.collective.resolved) > 4) & (sum(LG.2.marker.names %in% reconcile.collective.resolved) > 4)) {
        stop("intensive group markers in both reconcile groups!")
      } else if (sum(LG.1.marker.names %in% reconcile.collective.resolved) > 4 ) {
        intensive.collective <- c(intensive.collective, LG.1.marker.names)
      } else if (sum(LG.2.marker.names %in% reconcile.collective.resolved) > 4 ) {
        intensive.collective <- c(intensive.collective, LG.2.marker.names)
      } else {
        stop("intensive group markers do not meet threshold level of match!")
      }
    } # intensive.collective bracket
    # length(unique(intensive.collective))
    
    
    # This is a test code block to see whether there are any duplications of markers or same marker available in both phases.
    test.set <- intensive.collective %>% str_replace("_c", "") %>%  str_replace("_r", "") 
    # length(test.set)
    if (length(unique(test.set)) == length(test.set)) {
      # following code block converts the markers to haplotypes for the focal parent
      phase.df <- NULL
      for (m in 1:length(intensive.collective)) { # convert to 01 format
        
        tmp.1 <- unlist(str_split(intensive.collective[m], "_"))
        pos <- as.numeric(tmp.1[1])
        phase <- tmp.1[2]
        
        if (phase == "c") {
          tmp.2 <- c(0, 1)
        } else if (phase == "r") {
          tmp.2 <- c(1, 0)
        } else {
          stop("phase error!")
        }
        
        phase.df <- rbind(phase.df, c(chrom.name, pos, tmp.2))   
      }
      phase.df <- as.data.frame(phase.df, stringsAsFactors = FALSE)
      colnames(phase.df) <- c("CHROM", "POS", "HAP1", "HAP2")
      
      # order the markers based on their physical position.
      order.idx <- order(as.integer(phase.df$POS), decreasing = FALSE)
      phase.df.2 <- phase.df[order.idx, ]
      
      write.csv(phase.df.2, file = paste0(output.path.phased.markers, chrom.name, "_phased_markers_", par.name, ".csv"), row.names = FALSE, quote = FALSE)
      cat(paste0(output.path.phased.markers, chrom.name, "_phased_markers_", par.name, ".csv"), " Successfully completed!", "\n")
      # metadata.df <- rbind(metadata.df, c(chrom.name, par.name, nrow(par.het.haplo.tab), length(intensive.collective)))
    } else {
      stop("final marker set has duplications")
    }
  # } # par.name bracket
# } # chrom.name bracket
# colnames(metadata.df) <- c("CHROM", "PAR", "NUM_MRKRS_LD_METHOD", "NUM_MRKRS_ONEMAP_CLSTR")
# write.csv(metadata.df, file = paste0(output.path.phased.markers, row.names = FALSE, quote = FALSE)
########################################################################################################################################################################################################
# END
# The output of the script contains phased haplotype for the focal parent as follows;
    # CHROM,POS,HAP1,HAP2
    # Chr15,173833,0,1
    # Chr15,174711,0,1
    # Chr15,178446,1,0
    # Chr15,189430,1,0
    # Chr15,195830,0,1
    # Chr15,211967,1,0
    # Chr15,226300,0,1
    # Chr15,242915,0,1
    # Chr15,242997,0,1
########################################################################################################################################################################################################
# This code block checks the phasing accuracy of the LD based method. 
# phase.df.LD <- read.csv(file = paste0("./haplotypes/by_halfsib_family/", chrom.name, "_halfsibfamily_haplotypes_", par.name, ".csv"), stringsAsFactors = FALSE, header = TRUE)
# phase.df.LD <- phase.df.LD[, 1:4]
# common.idx <- which(phase.df.LD$POS %in% phase.df.2$POS)
# phase.df.LD.2 <- phase.df.LD[common.idx, ]
