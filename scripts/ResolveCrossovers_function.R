#########################################################################################################################################################################################################
# METADATA: Written by CRA and DMS, for the recombination landscape analysis for 7x7 cross
# DATE: May 14, 2020

# USAGE: This function script will be called by a script that identfies background parental haplotypes based on phased parental markers and offspring genotype observations. ex: 15b_big_window_both_parents_DMS.R.  
# function takes three arguments as input. 
# background.haps = smoothed.haps # This is the output of the ObtainConsensus() function called in the main script before current function is called.
# data.tab = off.info.tab # Following is an example of the information carried by this table; 
# CHROM   POS REF ALT 2048 2515 hap1_2048 hap2_2048 hap1_2515 hap2_2515 25663
# 1 Chr02  6617   T   C    0    2         0         0         1         1     1
# 2 Chr02  8424   A   G    0    2         0         0         1         1     1
# 3 Chr02 12411   C   G    0    2         0         0         1         1     1
# 4 Chr02 13186   T   G    0    2         0         0         1         1     1
# 5 Chr02 14242   G   A    0    2         0         0         1         1     1
# 6 Chr02 14663   G   A    0    2         0         0         1         1     1

# focal.par = mother.name
# other.par = father.name
# max.buffer = 12 # This is the maximum number of marker to be used to identify the crossover locations. 

# PURPOSE: Given background haplotype assignment using the SmoothHaplotypes function and ObtainConsensus function together, the begining and end of chromosomes as well as vicinity of the cross-overs 
# will contain markers that cannot be assigned to either haplotype with confidence. This script tries to delimit cross-overs to a much narrower delimitation.

# FUNCTTION DESCRIPTION: Resolves the cross-over locations to a much narrow delimitation.
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
#########################################################################################################################################################################################################
ResolveCrossovers <- function(background.haps, data.tab, focal.par, other.par, max.buffer){
  resolved.hap <- background.haps
  splitAt <- function(x, pos) {unname(split(x, cumsum(seq_along(x) %in% pos)))} # function splits avector of marker haplotype assignments to contiguous stretches of a given haplotype.
  max.dist.1s <- 1
  inconclusive.idx <- which(background.haps == 0)
  idx.diff.vec <- diff(inconclusive.idx)
  diff.idx <- which(idx.diff.vec > max.dist.1s)
  zero.island.list <- splitAt(inconclusive.idx, (diff.idx + 1))
  
  for (num in 1:length(zero.island.list)) { # num <- 1
    zero.island <- unlist(zero.island.list[num])
    zero.island.coord <- zero.island[1]:zero.island[length(zero.island)]
    if (any(zero.island == 1) | any(zero.island == length(background.haps))) { # if the stretch of undefined haplotypes are in the begining or the end of the chromosome they are not processed further.
      next("the zero island processed is at the begining or end of the chromosome!")
    } else {
      tmp.tab <- data.tab[zero.island[1]:zero.island[length(zero.island)], ]
      
      # the first code block tries to identify single marker haploype by considering other parent's genotype. 
      hap.assignment <- NULL
      for (i in 1:nrow(tmp.tab)) {
        focal.par.geno <- as.integer(tmp.tab[i, focal.par])
        if (focal.par.geno == 0 | focal.par.geno == 2) {
          hap <- 0
        } else {
          other.par.geno <- as.integer(tmp.tab[i, other.par])
          offspring.geno <- as.integer(tmp.tab[i, 11])
          hap.1 <- as.integer(tmp.tab[i, paste0("hap1_", focal.par)])
          hap.2 <- as.integer(tmp.tab[i, paste0("hap2_", focal.par)])
          
          if (other.par.geno == 0) {
            allele <- ifelse(offspring.geno == 2, NA,  offspring.geno)
          }
          if (other.par.geno == 2) {
            allele <- ifelse(offspring.geno == 0, NA, offspring.geno - 1)
          }
          if (other.par.geno == 1) {
            if (offspring.geno == 0) {
              allele <- 0
            } else if (offspring.geno == 2) {
              allele <- 1
            } else {
              allele <- NA
            }
          }
          
          if (!(is.na(allele))) {
            hap <- which(c(hap.1, hap.2) == allele)
          } else {
            hap <- 0
          }
        }
        hap.assignment <- c(hap.assignment, hap)
      }
      # plot(hap.assignment)
      
      # following code block carries-out a smoothing opertaion on the identified haplotypes for the single markers. This is similar but not the same as SmoothHaplotypes function in a different script.
      fwd.rev.single.mrkr.haps <- NULL
      for (direction in c("FWD", "REV")) {
        buffer <- c(0, 0)
        # max.buffer <- 12
        if (direction == "FWD") {
          iter.seq <- c(1:length(hap.assignment))
        } else {
          iter.seq <- c(length(hap.assignment):1)
        }
        smooth.hap.vec <- NULL
        for (win.num in iter.seq) {
          hap <- hap.assignment[win.num] # identify single marker hap assignment based on first code block
          if (hap != 0) {
            if (buffer[hap] < max.buffer) { # when haplotype assignment is 1 or 2, then that haplotype's buffer is incremented by one point and other reduced by one. 
                                            # The highest allowed for the buffer is defined in the function. 
              buffer[hap] <- buffer[hap] + 1 # increase buffer score of the haplotype
            }
            if (buffer[3 - hap] > 0) { # decrease one from the other haplotype's buffer score
              buffer[3 - hap] <- buffer[3 - hap]  - 1
            }
          }
          if (buffer[1] == buffer[2]) { # when buffer score for neither haplotype is prominent, zero is assigned as the haplotype
            impute.hap <- 0
          } else {
            impute.hap <- which.max(buffer) # assign haplotype with the highest buffer score
          }
          smooth.hap.vec <- c(smooth.hap.vec, impute.hap)
        }
        if (direction == "REV") {
          smooth.hap.vec <- rev(smooth.hap.vec) # correcting for the orientation of the reverse direction smoothing
        }
        fwd.rev.single.mrkr.haps <- as.data.frame(cbind(fwd.rev.single.mrkr.haps, smooth.hap.vec))
      }
      consensus.hap <- apply(fwd.rev.single.mrkr.haps, 1, function(x){ # taking the consensus of the FWD and REV direction smoothing outputs. Zero is assigned if the FWD and REV assignments do not agree.
        ifelse(x[1] == x[2], x[1], 0)
      })
    }
    # plot(consensus.hap)
    resolved.hap[zero.island.coord] <- consensus.hap
  }
  return(resolved.hap)
}
# resolved.background.haplotype <- ResolveCrossOvers(background.haps = temp, data.tab = off.info.tab, focal.par = mother.name, other.par = father.name, max.buffer = 12)
# plot(background.haps, pch = 16, cex = 0.2, col = "black")
# plot(resolved.background.haplotype, pch = 16, cex = 0.2, col = "black")
# #########################################################################################################################################################################################################
# END