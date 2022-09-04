#########################################################################################################################################################################################################
# METADATA: Written by CRA and DMS, for the recombination landscape analysis for 7x7 cross
# DATE: December 27, 2020

# USAGE: This function will be called by a script that needs to create a bin.span.coordinate.df  
# function takes four arguments as input. 
# chrom.name
# bin.size - window size in bps
# chrom.length - total contig size
# shift.size - shift scale in bps
# data.tab = off.info.tab # Following is an example of the information carried by this table; 

# PURPOSE: Create overlppin/non-overlapping window sequence for each chromosome
# FUNCTTION DESCRIPTION: Given a chromosome size in bps and a window size as well as a shift size creates windows and their indices.
#########################################################################################################################################################################################################
# Following is a function that divides a given chromosome into non overlapping windows. 
CreateBinSpanCoordinates <- function(chrom.name, bin.size, chrom.length, shift.size) {
  bin.spacing <- seq(from = bin.size, to = chrom.length, by = shift.size)
  num.pos <- length(seq(from = bin.size, to = chrom.length, by = shift.size))
  
  snipping.size <- chrom.length - bin.spacing[num.pos]
  if (snipping.size %% 2 == 0) {
    snipping.size.one.end = snipping.size / 2
  } else {
    snipping.size.one.end = (snipping.size + 1) / 2
  }
  
  end.sites <- bin.spacing + snipping.size.one.end
  start.sites <- (end.sites - bin.size) + 1
  start.end.sites <- cbind(start.sites, end.sites)
  if (start.end.sites[1,1] > 1) {
    start.end.sites <- rbind(c(1, (start.end.sites[1,1] - 1)), start.end.sites)
  } 
  if (start.end.sites[nrow(start.end.sites),2] < chrom.length) {
    start.end.sites <- rbind(start.end.sites, c((start.end.sites[nrow(start.end.sites),2] + 1), chrom.length))
  }
  return(start.end.sites)
}
# #########################################################################################################################################################################################################
# END