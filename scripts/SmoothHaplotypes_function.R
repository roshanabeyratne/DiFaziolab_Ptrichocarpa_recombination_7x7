#########################################################################################################################################################################################################
# METADATA: Written by CRA and DMS, for the recombination landscape analysis for 7x7 cross
# DATE: May 14, 2020

# USAGE: This function script will be called by a script that identfies background parental haplotypes based on phased parental markers and offspring genotype observations. ex: 15b_big_window_both_parents_DMS.R.  
# function takes two arguments as input. 
# hap.vec = haplotypes.vec # This is the background haplotype based on a sliding large window smoothing 
# max.buffer = 5 # This is the maximum number of such windows to consider in each direction of smoothing.

# PURPOSE: takes in a vector of observed haplotype assignments with noise and tries to identify contiguous stretches of the same haplotype. Initial haplotype assignments can be 0,1 or 2. 
# This function can be used to identify the background parental haplotypes across the chromosome for a given offspring individual and delimits the cross-over locations. 

# FUNCTTION DESCRIPTION: The vector is processed in both directions (FWD & REV), to identify contiguos stretches of possible background haplotypes for each offspring. At each progression 
# for a given direction, the leading haplotype (haplotype with the highest score) will be assigned to the window.
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
#########################################################################################################################################################################################################
SmoothHaplotypes <- function(hap.vec, max.buffer){
  fwd.rev.smoothed.haps <- NULL
  for (direction in c("FWD", "REV")) {
    buffer <- c(0, 0) # initial state
    # max.buffer <- 5 # this can be specified in the function as 'max.buffer'
    
    if (direction == "FWD") {
      iter.seq <- c(1:length(hap.vec))
    } else {
      iter.seq <- c(length(hap.vec):1)
    }
    
    smooth.hap.vec <- NULL
    for (win.num in iter.seq) {
      hap <- hap.vec[win.num] # identifies initial hap assignment to the large bin (1 or 2 are informative and 0 is not)
      if (hap != 0) {
        if (buffer[hap] < max.buffer) { # when the large window considered has a haplotype assignment 1 or 2, then buffer-score for that haplotype is incremented by one point and alt reduced by one. 
                                        # The highest allowed for the buffer is 5 and lowest is 0. 
          buffer[hap] <- buffer[hap] + 1 # increase buffer score of the haplotype identified in the large window
        }
        if (buffer[3 - hap] > 0) { # decrease one from the other haplotype only if the buffer-score is greater than zero.
          buffer[3 - hap] <- buffer[3 - hap]  - 1
        }
      }
      if (buffer[1] == buffer[2]) { # in the case where buffer-scores are equal the haplotype is not clear
        impute.hap <- 0
      } else {
        impute.hap <- which.max(buffer) # assign haplotype with the highest buffer score
      }
      smooth.hap.vec <- c(smooth.hap.vec, impute.hap)
    }
    if (direction == "REV") {
      smooth.hap.vec <- rev(smooth.hap.vec) # correcting for the orientation of the reverse direction smoothing
    }
    fwd.rev.smoothed.haps <- as.data.frame(cbind(fwd.rev.smoothed.haps, smooth.hap.vec))
  }
  # plot(fwd.rev.smoothed.haps[,1], col = "black", pch = 16, cex = 0.25)
  # plot(fwd.rev.smoothed.haps[,2], col = "blue", pch = 16, cex = 0.25)
  # plot(haplotypes.vec, col = "red", pch = 16, cex = 0.25)
  # sum(fwd.rev.smoothed.haps[, 1] != fwd.rev.smoothed.haps[, 2])
  # which(fwd.rev.smoothed.haps[, 1] != fwd.rev.smoothed.haps[, 2])
  
  smoothed.haps <- apply(fwd.rev.smoothed.haps, 1, function(x){ # taking the consensus of the FWD and REV direction smoothing outputs. Zero is assigned if the FWD and REV assignments do not agree.
    ifelse(x[1] == x[2], x[1], 0)
  })
  # sum(smoothed.haps == 0)
  return(smoothed.haps)
}
#########################################################################################################################################################################################################
# END