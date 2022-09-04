# METADATA: Written by CRA and DMS, for the recombination landscape analysis for 7x7 cross

# DATE: June 05, 2020

# FUNCTTION DESCRIPTION: # Creates input files in onemap format for markers in both coupled and repulsion, given an allele table for offspring and a filename for the .raw file to be outputted.  

# PURPOSE: Given a table with offspring and their allele contribution from the focal parent, this function creates the necessary input files for Onemap-R-package. The function needs this input file
# in 0,1,NA format. 0-Ref, 1-Alt, NA-missing/cannot be called. Following is the input file format

# CHROM    POS REF ALT 7073 24748 24749 24752 24753 24756 26643 26649
# Chr19 527031   T   C    1     1     0    NA     0     1     0     1
# Chr19 558957   A   C    1     1     1     0     0     1    NA     1
# Chr19 561141   T   C    1     1     1     0    NA     1     0     1
# Chr19 562431   C   G    1     1     1    NA     0     1     0     1
# Chr19 660818   T   C    1    NA     1     0     0     1     0     1
# Chr19 755105   G   A    1     1     1     0     0     1     0     1

# REQUIREMENTS: dplyr , stringr R-libraries are required.
########################################################################################################################################################################################################

MakeOnemapInput <- function(allele.tab, output.file.name) {
  
  library(dplyr)
  library(stringr)
  
  # coding for the coupled configuration, where Ref allele is coded as a and Alt allele as ab
  onemap.format.1 <- apply(allele.tab[6:ncol(allele.tab)], c(1,2), function(x){ 
    x <- as.integer(x)
    if (is.na(x)) {
      genotype <- "-"
    } else {
      if (x == 0) {
        genotype <-  "a"
      } else if (x == 1) {
        genotype = "ab"
      } else {
        stop("error in genotype format!")
      }
    }
    return(genotype)
  })
  onemap.format.1 <- as.data.frame(onemap.format.1, stringsAsFactors = FALSE)
  colnames(onemap.format.1) <- colnames(onemap.format.1) %>% str_replace("^X", "") # gets rid of the unnecessary X prefixes, (if there are any), when reading colnames that are numbers 
  
  
  # coding for the repulsion configuration, where Ref allele is coded as ab and Alt allele as a
  onemap.format.2 <- apply(allele.tab[6:ncol(allele.tab)], c(1,2), function(x){
    x <- as.integer(x)
    if (is.na(x)) {
      genotype <- "-"
    } else {
      if (x == 0) {
        genotype <-  "ab"
      } else if (x == 1) {
        genotype = "a"
      } else {
        stop("error in genotype format!")
      }
    }
    return(genotype)
  })
  onemap.format.2 <- as.data.frame(onemap.format.2, stringsAsFactors = FALSE)
  colnames(onemap.format.2) <- colnames(onemap.format.2) %>% str_replace("^X", "")
  
  
  # produce marker names for coupled and repulsion markers
  markers.cup <- NULL
  markers.rep <- NULL
  for (i in 1:nrow(allele.tab)) {
    #i<- 1
    markers.cup <- c(markers.cup, paste0("*",allele.tab$POS[i],"_c"))
    markers.rep <- c(markers.rep, paste0("*",allele.tab$POS[i],"_r"))
  }
  
  
  # produce individual names. This is not the same as colnames as Onemap requires them to be I1, I2, I3 ... 
  ind.names <- NULL
  for (i in 1:ncol(onemap.format.1)) {
    #i<- 1
    temp1 <- paste0("I",i)
    ind.names <- c(ind.names,temp1)
  }
  
  
  # compile metadata lines for the input file
  first.line <- noquote(paste("data", "type", "f2", "backcross", sep = " "))
  ind.names.line <- noquote(paste(ind.names,collapse = " "))
  metadata <- noquote(paste(length(ind.names), nrow(onemap.format.1)*2, "0","0","0", sep = " "))
  
  
  # create a filehandle for the output and write elements iteratively
  onemap.file = output.file.name
  unlink(onemap.file)
  rm(FILE1)
  
  FILE1 <- file(onemap.file,"w")
  write(first.line, FILE1, append = TRUE)
  write(metadata, FILE1, append = TRUE)
  write(ind.names.line, FILE1, append = TRUE)
  
  for (i in 1:nrow(onemap.format.1)) {
    #i<- 1
    line.cup <- noquote(paste(markers.cup[i],"A.H",paste(onemap.format.1[i,],collapse = " "),sep = " "))
    write(line.cup, FILE1, append = TRUE)
    line.rep <- noquote(paste(markers.rep[i],"A.H",paste(onemap.format.2[i,],collapse = " "),sep = " "))
    write(line.rep, FILE1, append = TRUE)
  }
  
  close(FILE1)
}
########################################################################################################################################################################################################
# END



 













