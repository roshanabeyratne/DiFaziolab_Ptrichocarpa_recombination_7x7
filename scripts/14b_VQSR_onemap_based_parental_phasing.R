#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: July 14, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and par.name.

# This is the input file format of haplotype contribution from the focal parent
# CHROM,POS,REF,ALT,2066,25376,25380,25382,25383,25384,25386,25902,26355,26357,26418,26422,26423,26424,26427,26439,26443,26446,26447,25820,25821,25822,25823,25825,25826,25828,25829,25830,25831,25832,
# Chr15,8393,C,T,1,0,NA,NA,1,0,NA,1,1,0,NA,1,NA,NA,1,0,0,0,1,0,0,0,1,0,0,NA,1,1,1,NA,1,0,0,1,1,NA,NA,1,NA,0,NA,0,NA,NA,NA,1,0,1,0,NA,1,NA,NA,0,1,NA,NA,0,0,NA,0,0,1,0,1,0,1,NA,0,1,1,1,0,1,1,1,0,1,0,0,
# Chr15,8421,C,T,1,0,NA,NA,1,0,NA,1,1,NA,NA,1,NA,NA,1,0,0,0,1,0,0,0,1,0,0,NA,1,1,1,NA,1,0,0,1,1,NA,NA,1,NA,0,NA,0,NA,NA,NA,1,0,1,0,NA,1,NA,NA,0,1,NA,NA,0,0,NA,0,0,1,0,1,0,1,NA,0,1,1,1,0,1,1,1,0,1,0,0

# This is the input file format for the already phased TRUTH dataset markers from an earlier script (Script 7b_onemap_based_parental_phasing.R). file has phased haplotype for the focal parent;
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

# PURPOSE: Phase heterozygous parental markers from the [VQSR + Truth] datasets based on group::Onemap() function. This function identifies marker linkage groups, using results from two-point
# (pairwise) analysis and the _transitive_ property of linkage.twopt method in Onemap. Markers are coded in both coupled (0|1) and repulsion (1|0) phases. 
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

TRUTH.input.path.combined.haplo <- "./haplotypes/by_parent/"
VQSR.input.path.combined.haplo <- "./haplotypes/by_parent/VQSR/"
input.path.phased.truth.markers <- "./onemap/phase_parents/phased_markers/"

tmp.input.path <- "./onemap/phase_parents/input_files/VQSR/"
tmp.output.path <- "./onemap/phase_parents/output_files/VQSR"
output.path.phased.markers <- "./onemap/phase_parents/phased_markers/VQSR/"
########################################################################################################################################################################################################

# for (chrom.name in chrom.name.vec[1:19]) {# chrom.name <- chrom.name.vec[2]
#   for (par.name in par.name.vec[1:14]) {# par.name <- par.name.vec[8]
    
    
    # read-in haplotpes from truth-dataset and the phased markers.
    TRUTH.parent.haplo.tab <- read.csv(file = paste0(TRUTH.input.path.combined.haplo, chrom.name, "_haplotypes_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    colnames(TRUTH.parent.haplo.tab) <- colnames(TRUTH.parent.haplo.tab) %>% str_replace(".GT", "") %>% str_replace("^X", "")
    
    # read-in the phased markers from the output of an earlier script (Script 7b_onemap_based_parental_phasing.R)
    truth.phase.df <- read.csv(file = paste0(input.path.phased.truth.markers, chrom.name, "_phased_markers_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    TRUTH.parent.haplo.tab.het <- TRUTH.parent.haplo.tab[TRUTH.parent.haplo.tab$POS %in% truth.phase.df$POS, ]
    if (!all(TRUTH.parent.haplo.tab.het[,par.name] == 1)) {
      stop("phased truthdataset markers are not all heterozygous!")
    }
    
    
    # read-in haplotpes from  VQSR-dataset markers. This dataset contains bi-allelic SNPs and INDELs
    VQSR.parent.haplo.tab <- read.csv(file = paste0(VQSR.input.path.combined.haplo, chrom.name, "_haplotypes_", par.name, ".csv"), header = TRUE, stringsAsFactors = FALSE)
    colnames(VQSR.parent.haplo.tab) <- colnames(VQSR.parent.haplo.tab) %>% str_replace(".GT", "") %>% str_replace("^X", "")
    
    # Selecting markers for which focal parent is heterozygous in the VQSR markers
    VQSR.parent.haplo.tab.het <- VQSR.parent.haplo.tab[VQSR.parent.haplo.tab[,par.name] == 1, ]
    VQSR.parent.haplo.tab.het <- VQSR.parent.haplo.tab.het[!is.na(VQSR.parent.haplo.tab.het$POS), ]
    
    if (all(VQSR.parent.haplo.tab.het[,par.name] == 1)) { # check if all the subset markers are Het in the focal parent.
      if (all(!(VQSR.parent.haplo.tab.het$POS %in% TRUTH.parent.haplo.tab.het$POS))) { # check that none of the markers in the VQSR set is duplicated in the TRUTH set 
        
        
        # following code block identifies indices of markers for which information content in offspring is low. 
        # for clarification, following are the conditions that would make focal parent allele contribution in offspring unreliable. Therefore NA is assigned for the aoffspring.
        # offspring genotype is missing or the other parent of the offspring has a missing genotype.
        # offspring genotype is not consistent with the expected genotype given parental genotyes (mendelian violation).
        less.informative.marker.count <- apply(VQSR.parent.haplo.tab.het[,6:ncol(VQSR.parent.haplo.tab.het)], 1, function(x){ # offspring information start from column-6 onwards.
          count.Nas <- sum(is.na(x))
          return(count.Nas)
        })
        # hist(less.informative.marker.count, breaks = 100)
        
        # removing markers that have more than 25% of offspring with missing information (NA)
        num.offspring <- (ncol(VQSR.parent.haplo.tab.het) - 5)
        less.informative.marker.idx <- which(less.informative.marker.count > (num.offspring * 0.25)) 
        VQSR.parent.haplo.tab.het <- VQSR.parent.haplo.tab.het[-less.informative.marker.idx, ]
      } else {
        stop("markers are duplicated in the VQSR and TRUTH sets!")
      }
    } else {
      stop("not all markers in VQSR.parent.haplo.tab.het are heterozygus!")
    }
    
    ####################################################################################################################################################################################################
    # In the following section  VQSR markers are clustered in batches with a subset of already phased TRUTH dataset stretched along the chromosome.
    # the maximum number of markers in a batch of VQSR markers are specified first. These markers should be mostly evenly stretched along the chromosome. Here these batches are called intensive sets. 
    # The number of intensive groups is calculated afterwards. 
    max.intensive <- 200
    intensive.group.num <- round(nrow(VQSR.parent.haplo.tab.het) / max.intensive)
    critical.length <- 205
    
    intensive.group.df <- NULL
    for (b in 1:intensive.group.num) { # b <- 1
      seq.tmp <- seq(from = b, to = nrow(VQSR.parent.haplo.tab.het), by = intensive.group.num)
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
    # A common subset of TRUTH dataset markers need to be added to the VQSR intensive groups. The phase of this subset is already known. Following code block specifies this subset. 
    max.phased.truth.mrkrs <- 300  # the size of the subset is arbitrarily assigned as 300 markers.
    seq.truth <- unique(floor(seq(from = 1, to = nrow(TRUTH.parent.haplo.tab.het), length.out = max.phased.truth.mrkrs)))
    reconcile.tab.het <- TRUTH.parent.haplo.tab.het[seq.truth, ] # subset rconcilation marker haplotypes
    recon.mrkr.phase.tab <- truth.phase.df[truth.phase.df$POS %in% reconcile.tab.het$POS, ]
    
    # identify the phase of subset of the reconc markers ( I have used HAP1 as the selected phase).
    recon.mrkr.phase <- as.vector(apply(recon.mrkr.phase.tab, 1, function(x){
      chrom <- x[1]
      pos <- as.numeric(x[2])
      phase <- (x[3])
      if (phase == 0) {
        mrkr <- paste0(pos, "_c")
      } else if (phase == 1) {
        mrkr <- paste0(pos, "_r")
      } else {
        stop("error in marker phase input of the truthset!")
      }
      return(mrkr)
    }))
    
    ####################################################################################################################################################################################################
    # Following clusters the intensive group markers iteratively along with the subset of phased TRUTH dataset markers. twopt method of Onemap-R package is used to calculate pairwise recombination 
    # fractions. Onemap::group function is used identify linkage groups. At the final step, out of the two largest linkage groups, the linkage group that corresponds HAP 1 of the TRUTH subset is used 
    # as the scaffold of choice arbitrarily to identify HAP 1 of VQSR marker haplotype configuarion.  
    intensive.collective <- NULL
    # for (n in c(1:ncol(intensive.group.df))[-c(30, 112)]) { # this is for "1950" and Chr01
    for (n in 1:ncol(intensive.group.df)) { # n <- 80 #  this is the number of intensive groups for the VQSR dataset.
      
      tmp.tab <- rbind(VQSR.parent.haplo.tab.het[unique(intensive.group.df[, n]), ], reconcile.tab.het) # combine VQSR subset and TRUTH subsets into one dataframe
      ordrd.idx <- order(tmp.tab$POS, decreasing = FALSE)
      tmp.tab <- tmp.tab[ordrd.idx, ]
      
      
      # Creating input files in onemap format for markers in both coupled and repulsion.      
      onemap.file <- paste0(tmp.input.path, chrom.name, "_tmp_intensive_", n, "_", par.name, ".raw")
      MakeOnemapInput(allele.tab = tmp.tab, output.file.name = onemap.file) # function sourced at the begining of the script. Makes all the input files in Onemap format. 
      
      
      # phasing markers in the intensive group using Onemap:R-package. The markers in the same linkage group are the phased haplotype for the focal parent.
      db.raw <- read_onemap(inputfile = onemap.file)
      LOD.sug <- suggest_lod(db.raw)
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
      
      #  This section checks whether markers in coupled and repulsion configurations are both clustered in the same linkage group. 
      tmp.marker.names.1.dup  <- LG.1.marker.names %>% str_replace("_c", "") %>%  str_replace("_r", "") %>% duplicated
      tmp.marker.names.2.dup <- LG.2.marker.names %>% str_replace("_c", "") %>%  str_replace("_r", "") %>% duplicated
      if (any(tmp.marker.names.1.dup) | any(tmp.marker.names.2.dup)) {
        next("marker duplicated!")
      } else {
        if ((sum(LG.1.marker.names %in% recon.mrkr.phase) > 200 ) & (sum(LG.2.marker.names %in% recon.mrkr.phase) > 200 )) { # 200/300 TRUTH subset markers are required in a linkage group to proceed.
          stop("intensive group markers in both reconcile groups!")
        } else if (sum(LG.1.marker.names %in% recon.mrkr.phase) > 200 ) {
          # the TRUTH subset markers are removed from the linkage group
          tmp.mrkrs.1  <- LG.1.marker.names %>% str_replace("_c", "") %>%  str_replace("_r", "") %>% as.numeric()
          mrkrs.to.rm <- which(tmp.mrkrs.1 %in% truth.phase.df$POS)
          mrkrs.to.add <- LG.1.marker.names[-(mrkrs.to.rm)]
          
          intensive.collective <- c(intensive.collective, mrkrs.to.add)
        } else if (sum(LG.2.marker.names %in% recon.mrkr.phase) > 200 ) {
          # the TRUTH subset markers are removed from the linkage group
          tmp.mrkrs.1  <- LG.2.marker.names %>% str_replace("_c", "") %>%  str_replace("_r", "") %>% as.numeric()
          mrkrs.to.rm <- which(tmp.mrkrs.1 %in% truth.phase.df$POS)
          mrkrs.to.add <- LG.2.marker.names[-(mrkrs.to.rm)]
          
          intensive.collective <- c(intensive.collective, mrkrs.to.add)
        } else {
          stop("intensive group markers do not meet threshold level of match!")
        }
      }
      # sum(LG.2.marker.names %in% recon.mrkr.phase)
    } # intensive.collective bracket
    # following code block tests whether there are marker duplications within the reconcile.collective.resolved set
    # tmp.set <- intensive.collective %>% str_replace("_c", "") %>%  str_replace("_r", "") 
    # length(tmp.set) == length(unique(tmp.set))
    
    ########################################################################################################################################################################################################
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
#   } # par.name bracket
# } # chrom.name bracket
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
