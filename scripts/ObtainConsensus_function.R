#########################################################################################################################################################################################################
# METADATA: Written by CRA and DMS, for the recombination landscape analysis for 7x7 cross
# DATE: May 14, 2020

# USAGE: This function script will be called by a script that identfies background parental haplotypes based on phased parental markers and offspring genotype observations. ex: 15b_big_window_both_parents_DMS.R.  
# function takes three arguments as input. 
# hap.vec = haplotypes.vec # This is the background haplotype based on a sliding large window smoothing 
# data.tab = off.info.tab # Following is an example of the information carried by this table; 
# CHROM   POS REF ALT 2048 2515 hap1_2048 hap2_2048 hap1_2515 hap2_2515 25663
# 1 Chr02  6617   T   C    0    2         0         0         1         1     1
# 2 Chr02  8424   A   G    0    2         0         0         1         1     1
# 3 Chr02 12411   C   G    0    2         0         0         1         1     1
# 4 Chr02 13186   T   G    0    2         0         0         1         1     1
# 5 Chr02 14242   G   A    0    2         0         0         1         1     1
# 6 Chr02 14663   G   A    0    2         0         0         1         1     1

# coord.tab = bin.span.tab # This carries the start and end site coordinates for each window in following format; 
# start.sites end.sites
# 1           1       313
# 2         220       355
# 3         314       385
# 4         356       429
# 5         386       454
# 6         430       465

# PURPOSE: Obtaining the consensus of the sliding windows. Giving values to the markers based on the overlapping windows to which they belong. 
# Most markers (except at the start and end of chromosome), are in two windows (windows overlap), so they will have two values, one for the odd windows and one for the even

# FUNCTTION DESCRIPTION: obtains the consensus haplotype for the sliding window analysis.
# key for obtaining the consensus of odd & even windows  
# 0 0 -> 0
# 0 1 -> 1
# 0 2 -> 2
# 1 1 -> 1
# 2 2 -> 2
# 1 2 -> 0
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
#########################################################################################################################################################################################################
ObtainConsensus <- function(hap.vec, data.tab, coord.tab){
  odd.seq <- seq(from = 1, to = length(hap.vec), by = 2) # all the odd numbered windows
  even.seq <- seq(from = 2, to = length(hap.vec), by = 2) # all the even numbered windows
  even.bin.tab <- data.tab
  even.bin.tab[, 11] <- rep(0, nrow(even.bin.tab))
  odd.bin.tab <- even.bin.tab
  
  # assign all markers as per odd-numbered windows
  for (k in odd.seq) {
    # k<- 1
    start <- coord.tab$start.sites[k]
    end <- coord.tab$end.sites[k]
    odd.bin.tab[start:end, 11] <- rep(hap.vec[k], length(odd.bin.tab[start:end, 11]))
  }
  # plot(odd.bin.tab[,11], pch = 16, cex = 0.25, col = 'black')
  
  # assign all markers as per even-numbered windows
  for (k in even.seq) {
    # k<- 2
    start <- coord.tab$start.sites[k]
    end <- coord.tab$end.sites[k]
    even.bin.tab[start:end, 11] <- rep(hap.vec[k], length(even.bin.tab[start:end, 11]))
  }
  # plot(even.bin.tab[, 11], pch = 16, cex = 0.25, col = 'blue')
  
  # takes the consensus of odd and even windows
  consensus <- apply(cbind(odd.bin.tab[, 11], even.bin.tab[, 11]), 1, function(x){
    tmp <- ifelse(x[1] == x[2], x[1], (x[1] + x[2]))
    return(tmp)
  })
  if (sum(consensus == 3) > 0) { 
    consensus[which(consensus == 3)] <- 0 
  }
  return(consensus)
}
#########################################################################################################################################################################################################
# END