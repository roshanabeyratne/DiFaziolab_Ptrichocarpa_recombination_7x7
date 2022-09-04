#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("J:/Stettler_14_linkage_mapping")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA and DMS, for the recombination landscape analysis for 7x7 cross
# DATE: May 12, 2020; Modified July 17, 2020
# USAGE: This script is used along with a bash submission script that specifies the chrom.name and mother.name and father.name. Input files include phased markers of both parents in the following format; 
# CHROM   POS HAP1 HAP2
# 1 Chr02 18950    0    1
# 2 Chr02 42845    1    0
# 3 Chr02 49350    1    0
# 4 Chr02 54307    1    0
# 5 Chr02 55886    1    0
# 6 Chr02 59320    1    0

# as well as an input file that contains all the genotypes for a given full-sib family. The parent genotypese are the corrected genotypes based on the likelihood based correction. 
# CHROM   POS REF ALT 2048 2515 25663 25666 25667 25670 25672 25673 25674 25677 25679 25684 25685 25688 25691 25692 25697 25699 25700 25702
# 1 Chr02  6617   T   C    0    2     1     1     1     1     1     1     1     1     1     1     1    NA     2     1     1     1     1     1
# 2 Chr02  8424   A   G    0    2     1     1     1     1     1     1     1     1     1     1     1     1    NA     0     1     0     1     1
# 3 Chr02 12411   C   G    0    2     1     1     1     1     1     1     1     1     0     1     1     1     1     1     1     1     1     1
# 4 Chr02 13186   T   G    0    2     1     1     0     1     0     1     1     1     1     1     0     0    NA     0     1     1     1     1
# 5 Chr02 14242   G   A    0    2     1     1     1     1     1     1     1     1     1     1     1     0     1     1     1     1     0     1
# 6 Chr02 14663   G   A    0    2     1     1     1     1     1     1    NA     1     1     1     1     1     1     1     1     1     1     1

# PURPOSE: This script identifies parental haplotypes across the chromosome for a given offspring individual with high confidence using a sliding window approach. This analysis leaves out any window 
# that cannot be confidently assigned to either haplotype as zero/unknown for further downstream analysis. Furthermore, the script calls functions that smooth identified haplotypes as well as resolve
# cross-over locations to a narrow margin whcih then enables downstream analysis. The output format is per offspring.
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
source(file = "./scripts_lists/SmoothHaplotypes_function.R")
# source(file = "./new_scripts/SmoothHaplotypes_function.R")
source(file = "./scripts_lists/ObtainConsensus_function.R")
# source(file = "./new_scripts/ObtainConsensus_function.R")
source(file = "./scripts_lists/ResolveCrossovers_function.R")
# source(file = "./new_scripts/ResolveCrossovers_function.R")
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073")
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

OM.phased.mrkr.input.path <- "./onemap/phase_parents/phased_markers/"
# OM.phased.mrkr.input.path <- "./test_folder/"
genotype.input.path.012 <- "./genotypes/012/"
# genotype.input.path.012 <- "./test_folder/"
haps.resolved.output.path <- "./resolved_haps/"
########################################################################################################################################################################################################
# for (chrom.name in chrom.name.vec[2]) {# chrom.name <- chrom.name.vec[2]
#   for (mother.name in par.name.vec[1:7]) {# mother.name <- par.name.vec[3]
#     for (father.name in par.name.vec[8:14]) {# father.name <- par.name.vec[9]
      
      # read the OM phased Het markers for mother and father
      informative.mother.mrkrs <- read.csv(file = paste0(OM.phased.mrkr.input.path, chrom.name, "_phased_markers_", mother.name, ".csv"), header  = TRUE, stringsAsFactors = FALSE)
      informative.father.mrkrs <- read.csv(file = paste0(OM.phased.mrkr.input.path, chrom.name, "_phased_markers_", father.name, ".csv"), header  = TRUE, stringsAsFactors = FALSE)
      
      # read all the genotypes for fulsib family
      all.genos.tab <- read.csv(file = paste0(genotype.input.path.012, chrom.name, "_fulsib_genotypes_", mother.name, "X", father.name, ".csv"), header  = TRUE, stringsAsFactors = FALSE)
      colnames(all.genos.tab) <- colnames(all.genos.tab) %>% str_replace("^X", "")
      
      if (all(informative.mother.mrkrs$POS %in% all.genos.tab$POS) & all(informative.father.mrkrs$POS %in% all.genos.tab$POS)) { # all phased markers in mother & father  available in all.genos.tab
        
        # merge all.genos.tab with the phased mother haplotypes and reorganizing the table
        tmp.merge.1 <- merge(all.genos.tab, informative.mother.mrkrs, all.x = TRUE, by = c("CHROM", "POS"))
        tmp.idx <- ncol(tmp.merge.1)
        
        # the following code block is required in the case a full-sib family has one offspring and the subset of the genotable is a vector rather than a dataframe in R
        if (!(is.data.frame(tmp.merge.1[, 7:(tmp.idx - 2)]))) {
          offspring.data.table <- as.data.frame(tmp.merge.1[, 7:(tmp.idx - 2)])
          colnames(offspring.data.table) <- colnames(all.genos.tab)[ncol(all.genos.tab)]
        } else {
          offspring.data.table <- as.data.frame(tmp.merge.1[, 7:(tmp.idx - 2)])
        }
        
        tmp.merge.1 <- cbind(tmp.merge.1[,1:6], tmp.merge.1[,c((tmp.idx - 1), tmp.idx)], offspring.data.table)
        colnames(tmp.merge.1)[7:8] <- c(paste0("hap1_", mother.name), paste0("hap2_", mother.name))
        
        # merge all.genos.tab with the phased father haplotypes and reorganizing the table
        tmp.merge.2 <- merge(tmp.merge.1, informative.father.mrkrs, all.x = TRUE, by = c("CHROM", "POS"))
        tmp.idx <- ncol(tmp.merge.2)
        
        # the following code block is required in the case a full-sib family has one offspring and the subset of the genotable is a vector rather than a dataframe in R
        if (!(is.data.frame(tmp.merge.2[, 9:(tmp.idx - 2)]))) {
          offspring.data.table <- as.data.frame(tmp.merge.2[, 9:(tmp.idx - 2)])
          colnames(offspring.data.table) <- colnames(tmp.merge.1)[ncol(tmp.merge.1)]
        } else {
          offspring.data.table <- tmp.merge.2[, 9:(tmp.idx - 2)]
        }
        tmp.merge.2 <- cbind(tmp.merge.2[,1:8], tmp.merge.2[,c((tmp.idx - 1), tmp.idx)], offspring.data.table)
        colnames(tmp.merge.2)[9:10] <- c(paste0("hap1_", father.name), paste0("hap2_", father.name))
        
        # now tmp.merge.2 has all phased Het markers as well as HomRef and HomAlt markers for both parents
        if ((nrow(tmp.merge.1) == nrow(tmp.merge.2)) & (nrow(tmp.merge.1) == nrow(all.genos.tab))) {
          
          # assign all Hom positions in mother its relevant haplotypes. This is trivial since allles in haplotypes are uninformative.
          hom.ref.idx.mother <- which(tmp.merge.2[,mother.name] == 0 | tmp.merge.2[,mother.name] == 2)
          tmp.merge.2[hom.ref.idx.mother, 7:8] <- tmp.merge.2[hom.ref.idx.mother, mother.name] / 2 
          
          # assign all Hom positions in father its relevant haplotypes. This is trivial since allles in haplotypes are uninformative.
          hom.ref.idx.father <- which(tmp.merge.2[,father.name] == 0 | tmp.merge.2[,father.name] == 2)
          tmp.merge.2[hom.ref.idx.father, 9:10] <- tmp.merge.2[hom.ref.idx.father, father.name] / 2 
          
          # identify and exclude marker positions for which haplotypes are not informative either in mother or father. This is the case if the focal parent genotype is Het and the marker position was 
          # removed from the list of markers that was phased for that given parent.
          idx.to.exclude <- unique(c(which(is.na(tmp.merge.2[,7])), which(is.na(tmp.merge.2[,9])))) # hap1_ columns of both parents are looked at to identfy these marker positions.
          tmp.fs.tab <- tmp.merge.2[-idx.to.exclude, ]
          ordrd.idx <- order(tmp.fs.tab$POS, decreasing = FALSE)
          tmp.fs.tab <- tmp.fs.tab[ordrd.idx, ]
         
        } else {
          stop("merging tables was not successful!")
        }

      } else {
        stop("phased markers of mother or father missing from all.genos.tab!")
      }
      
      ##################################################################################################################################################################################################
      # following code block processes the offspring iteratively
      off.spr.count <- 0
      for (j in 11:ncol(tmp.fs.tab)) {# j <- 17 offspring columns start at 11th column.
        
        # markers with missing genotype in offspring are excluded. A missing genotype observation in offspring is not informative at this stage to identify the background parental haplotype. 
        off.info.tab <- tmp.fs.tab[!is.na(tmp.fs.tab[, j]), c(1:10, j)] # This is the offspring table
        consensus.bin.tab <- off.info.tab # This is a temporary place holder into which final result will be appended.
        
        temp_out <- NULL
        
        for (focal.parent in c(mother.name, father.name)) { # focal.parent <- mother.name
          num.het.mrkrs <- 20 # number of Het markers (these are the haplotype informative markers) that should be included in a window for the given focal parent
          slide.by <- 10
          break.points.1 <- NULL
          break.points.2 <- NULL
          count <- 0
          lag.count <- slide.by
          
          # Defining the limits of the windows. Each window contains 20 het markers in the focal parent and they overlap by 10.
          for (pos in 1:nrow(off.info.tab)) {
            # print(pos)
            if (off.info.tab[pos, focal.parent] == 1) {
              count <- count + 1
              lag.count <- lag.count + 1
            }
            if (count == num.het.mrkrs) {
              break.points.1 <- c(break.points.1, pos)
              count <- 0
            }
            if (lag.count == num.het.mrkrs) {
              break.points.2 <- c(break.points.2, pos)
              lag.count <- 0
            }
          }
          
          bin.span.tmp.1 <- NULL
          for (i in 1:length(break.points.1)) {
            if (i == 1) {
              tmp.s <- 1
            } else {
              tmp.s <- break.points.1[i - 1] + 1
            }
            tmp.e <- break.points.1[i]
            bin.span.tmp.1 <- rbind(bin.span.tmp.1, c(tmp.s, tmp.e))
          }
          
          
          bin.span.tmp.2 <- NULL
          for (i in 1:length(break.points.2)) {
            if (i == 1) {
              next
            } else {
              tmp.s <- break.points.2[i - 1] + 1
            }
            tmp.e <- break.points.2[i]
            bin.span.tmp.2 <- rbind(bin.span.tmp.2, c(tmp.s, tmp.e))
          }
          
          # This is a diagnostic block to check if all the windows hav 20 Het markers.
          checkpoint.vec <- NULL
          for (z in c(1,2)) { # z <- 1
            bin.span.tmp <- get(paste0("bin.span.tmp.",z))
            coll.vec <- NULL
            for (i in 1:nrow(bin.span.tmp)) {
              num.het <- sum(off.info.tab[bin.span.tmp[i,1]:bin.span.tmp[i,2], focal.parent] == 1)
              if (num.het != num.het.mrkrs) {
                stop("error in window with num.het.mrkrs!")
              } else {
                coll.vec <- c(coll.vec, num.het)
              }
            }
            checkpoint.vec <- c(checkpoint.vec, all(coll.vec == num.het.mrkrs))
          }
          
          if (all(checkpoint.vec)) {
            
            # Building a dataframe with ovralapping window cordinates in terms of start and end indices.
            bin.span.tmp <- rbind(bin.span.tmp.1, bin.span.tmp.2)
            bin.span.tmp[,1] <- sort(bin.span.tmp[,1], decreasing = FALSE)
            bin.span.tmp[,2] <- sort(bin.span.tmp[,2], decreasing = FALSE)
            bin.span.tab <- as.data.frame(bin.span.tmp) 
            colnames(bin.span.tab) <- c("start.sites", "end.sites")
            
          } else {
            stop("not all windows have 20 HEt markers!")
          }
          
          
          
          haplotypes.vec <- NULL
          allele.drop.rate <- 0.1535688 # estimated in a previous script in /Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping/scripts/allele_drop_out.R
          
          # Given a number of heterozygous markers in the parent half of them should be Het in a given offspring. This is the case because irrespective of the other parent's genotype the expectation 
          # of a Het genotype in the offspring is one-half. Therefore allele dropout would only affect half of the offpring markers. In addition, allele dropout could result in a genotype that is 
          # either 0 or 2, so half of the errors will match with the alternate haplotype, propagating the observation twice so that error would be observed twice. (mismatch for one haplotype is also a
          # match to the other haplotype widening the gap even more). This is why we multiply by 1.5 (1.0 + 0.5) expectation of Hom offspring is 0.5 and the allele.drop.rate calculated for these offspring
          cut.off <- round(num.het.mrkrs - (num.het.mrkrs * allele.drop.rate * 0.5) * 1.5)
          
          
          for (i in 1:nrow(bin.span.tab)) { # i <- 2 # specifies the window iteration number
            
            large.window.tab <- off.info.tab[bin.span.tab$start.sites[i]:bin.span.tab$end.sites[i], ] 
            large.window.tab <- large.window.tab[large.window.tab[,focal.parent] == 1,] # making a subset table with info for the window with 20  Het markers
            mother.hap.1 <- large.window.tab[1:nrow(large.window.tab), 7] # mother's first haplotype
            mother.hap.2 <- large.window.tab[1:nrow(large.window.tab), 8] # mother's second haplotype
            father.hap.1 <- large.window.tab[1:nrow(large.window.tab), 9] # father's first haplotype
            father.hap.2 <- large.window.tab[1:nrow(large.window.tab), 10] # father's second haplotype
            
            # these are the four possible haplotype combinations 
            d1s1 <- mother.hap.1 + father.hap.1 
            d1s2 <- mother.hap.1 + father.hap.2
            d2s1 <- mother.hap.2 + father.hap.1
            d2s2 <- mother.hap.2 + father.hap.2
            off.geno <- large.window.tab[1:nrow(large.window.tab), ncol(large.window.tab)] # offspring's genotype
            
            d1s1.mismatch <- sum((d1s1 != off.geno), na.rm = TRUE)
            d1s2.mismatch <- sum(d1s2 != off.geno, na.rm = TRUE)
            d2s1.mismatch <- sum(d2s1 != off.geno, na.rm = TRUE)
            d2s2.mismatch <- sum(d2s2 != off.geno, na.rm = TRUE)
            
            if (focal.parent == mother.name) {# for the mother 
              mismatch.min.hap.1 <- min(c(d1s1.mismatch, d1s2.mismatch))
              mismatch.min.hap.2 <- min(c(d2s1.mismatch, d2s2.mismatch))
            } else if (focal.parent == father.name) { # for father
              mismatch.min.hap.1 <- min(c(d1s1.mismatch, d2s1.mismatch))
              mismatch.min.hap.2 <- min(c(d1s2.mismatch, d2s2.mismatch))
            } else {
              stop("focal parent not defined in the iteration")
            }
            
            mismatch.min.hap <- which.min(c(mismatch.min.hap.1, mismatch.min.hap.2)) # which haplotype (hap1 or hap 2) contains the minimum number of mismatches to the observed offspring genotype
            min.amnt <- min(c(mismatch.min.hap.1, mismatch.min.hap.2)) # the number of mismatches out of a max of num.het.mrkrs
            mismatch.diff <- abs(mismatch.min.hap.1 - mismatch.min.hap.2) # this is the pivotal quantity that indicates whether there is a clear haplotype winner.
            
            # if (mismatch.diff == 0) { # test code block
            #   stop("mismatch.diff is zero!!")
            # }
            
            if (mismatch.diff >= cut.off) {
              haplotypes <- mismatch.min.hap
            } else {
              haplotypes <- 0
            }
            haplotypes.vec <- c(haplotypes.vec, haplotypes)
            
            temp_out <- rbind(temp_out,c(mismatch.min.hap.1,mismatch.min.hap.2, mismatch.diff,haplotypes))
          } # bin span end bracket
          # plot(haplotypes.vec, pch = 16, cex = 0.5, col = 'red', ylab = "haplotype", xlab = "window.index", main = "")
          
          ##############################################################################################################################################################################################
          # This code block calls functions that smooth haplotypes identified with high confidence using large windows. Smoothing is carried out to identify the background haplotype for the offspring. 
          # Clear information of the background haplotype for the offspring is required for ;
          # 1) produce genetic linkage maps for the focal parent
          # 2) identify possible gene conversion events 
          smoothed.windows <- SmoothHaplotypes(hap.vec = haplotypes.vec, max.buffer = 5)
          # plot(smoothed.windows, pch = 16, cex = 0.2, col = "red")

          # We noticed that, at times algorithm that produces haplotypes.vec above, can mistakenly construe that there is a double haplotype switch in areas in the genome where both focal parents are 
          # Het for considered markers and the offspring observed as HomRef or HomAlt. Higher allele drop-out in areas like this may lead the algorithm to derive a short double switch in focal parent
          # haplotype. This can be corrected by imposing a minimum size limit for a given haplotype block. However, this may only be suitable for CO event analysis and not GC analysis. 
          # This will be tackled in the next script 16b_.....
          
          smoothed.haps <- ObtainConsensus(hap.vec = smoothed.windows, data.tab = off.info.tab, coord.tab = bin.span.tab)
          # plot(smoothed.haps, pch = 16, cex = 0.2, col = "black")
          
          if (all(smoothed.haps == 1) | all(smoothed.haps == 2)) { # If crossover events are not present in the background haplotypes then ResolveCrossovers function is not implemented.
            resolved.background.haplotypes <- smoothed.haps
          } else {
            if (focal.parent == mother.name) {
              resolved.background.haplotypes <- ResolveCrossovers(background.haps = smoothed.haps, data.tab = off.info.tab, focal.par = mother.name, other.par = father.name, max.buffer = 12)
            } else {
              resolved.background.haplotypes <- ResolveCrossovers(background.haps = smoothed.haps, data.tab = off.info.tab, focal.par = father.name, other.par = mother.name, max.buffer = 12)
            }
          }
          # plot(resolved.background.haplotypes, pch = 16, cex = 0.2, col = "blue")
          
          
          # Assign the large window based haplotypes to markers without any smoothing. 
          consensus.haps <- ObtainConsensus(hap.vec = haplotypes.vec, data.tab = off.info.tab, coord.tab = bin.span.tab)
          # plot(consensus.haps, pch = 16, cex = 0.2, col = "red")
          consensus.bin.tab <- cbind(consensus.bin.tab, consensus.haps, resolved.background.haplotypes)
        } # focal.parent bracket
        
        # The final consensus.bin.tab produced by above code block that iterates over both focal parents has background haplotypes as smoothed as well as actual haplotypes observed that assist in 
        # identifying gene-conversion events
        colnames(consensus.bin.tab) <- c(colnames(off.info.tab), paste0(mother.name, "_hap_combo"), paste0(mother.name, "_background"), paste0(father.name, "_hap_combo"), paste0(father.name, "_background"))
        # plot(consensus.bin.tab[, 11], pch = 16, cex = 0.25, col = 'black')
        # plot(consensus.bin.tab[, 12], pch = 16, cex = 0.25, col = 'red')
        # plot(consensus.bin.tab[, 13], pch = 16, cex = 0.25, col = 'blue')
        # plot(consensus.bin.tab[, 14], pch = 16, cex = 0.25, col = 'red')
        # plot(consensus.bin.tab[, 15], pch = 16, cex = 0.25, col = 'blue')
        
        offspring.name <- colnames(consensus.bin.tab)[11]
        write.csv(consensus.bin.tab, file = paste0(haps.resolved.output.path, chrom.name, "_resolved_haplotypes_both_parents_", offspring.name,".csv"), quote = FALSE, row.names = FALSE)
        cat(paste0(haps.resolved.output.path, chrom.name, "_resolved_haplotypes_both_parents_", offspring.name,".csv"), " Successfully completed!", "\n")
        off.spr.count <- off.spr.count + 1
      }# offspring bracket
      if (off.spr.count == (ncol(all.genos.tab) - 6)) {
        cat(paste0("all offspring Successfully completed!", "\n"))
      } else {
        stop("mismatch in the number of offspring processed vs. actual available in the full-sib family!")
      }
      
#     }# father bracket
#   }# mother bracket
# }# chrom bracket
# consensus.bin.tab.tmp <- consensus.bin.tab
# consensus.bin.tab <- consensus.bin.tab.tmp

# temp_out1<-temp_out
# temp_out2<-temp_out
# 
# sum(temp_out1[,4 ]==0)
# sum(temp_out2[,4 ]==0)
# 
# hist(temp_out1[,1 ],breaks=100)
# hist(temp_out2[,1 ],breaks=100)
# hist(temp_out2[,3 ],breaks=100)
# ######################################################################################################################################################################################################
# END
# The output of the script contains background haplotypes as smoothed as well as actual haplotypes for a given offspring as folllows; 
# CHROM   POS REF ALT 2048 2515 hap1_2048 hap2_2048 hap1_2515 hap2_2515 25663 2048_hap_combo 2048_background 2515_hap_combo 2515_background
# Chr02  6617   T   C    0    2         0         0         1         1     1              0               0              0               0
# Chr02  8424   A   G    0    2         0         0         1         1     1              0               0              0               0
# Chr02 12411   C   G    0    2         0         0         1         1     1              0               0              0               0
# Chr02 13186   T   G    0    2         0         0         1         1     1              0               0              0               0
# Chr02 14242   G   A    0    2         0         0         1         1     1              0               0              0               0
# Chr02 14663   G   A    0    2         0         0         1         1     1              0               0              0               0
########################################################################################################################################################################################################