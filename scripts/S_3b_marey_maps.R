#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: August 1, 2020
# USAGE: Input files include genetic maps. The marker names contain the chromosome and marker putative physical position. 

# Following is the input file / genetic map format; 
# 1 Chr06_417767 0
# 1 Chr06_5256 9.99999999868873e-06
# 1 Chr06_685931 0.72899955945013
# 1 Chr06_472891 0.729009559450128
# 1 Chr06_811815 2.18926150399718
# 1 Chr06_804902 2.18927150399718
# 1 Chr06_811828 2.18928150399718
# 1 Chr06_1049955 2.91925669384554
# 1 Chr06_1037834 2.91926669384554
# 1 Chr06_1297578 3.64925286278189
# 1 Chr06_1139116 3.64926286278189
# 1 Chr06_1071037 3.64927286278189

# PURPOSE: Plot marrey maps with physical vs genetic map distance. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
# args <- commandArgs(TRUE)
# chrom.name <- args[1]
# par.name <- args[2]
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

genetic.map.input.path <- "./onemap/bin_markers/genetic_maps/"
rippled.genetic.map.input.path <- "./onemap/bin_markers/genetic_maps/rippled_maps/"
marey.plot.output.path <- "./onemap/bin_markers/genetic_maps/marey/"
########################################################################################################################################################################################################
mrkr.count <- NULL
map.size <- NULL
lg.names <- c("LG-I", "LG-II", "LG-III", "LG-IV", "LG-V", "LG-VI", "LG-VII", "LG-VIII", "LG-IX", "LG-X", "LG-XI", "LG-XII", "LG-XIII", "LG-XIV",
              "LG-XV", "LG-XVI", "LG-XVII", "LG-XVIII", "LG-XIX")

for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  for (par.name in par.name.vec) { # par.name <- "2365"
    
    if (file.exists(paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"))) {
      gen.map <- read.table(file = paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    } else {
      gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    }
    
    string <- paste0(chrom.name, "_")
    physical.distance <- gen.map[,2] %>% str_replace(string, "") %>% as.numeric()
    genetic.distance <- gen.map[,3] %>% as.numeric()
    
    mrkr.count <- c(mrkr.count, length(genetic.distance)) # This is the number of markers included in the genetic map
    map.size <- c(map.size, genetic.distance[length(genetic.distance)]) # This is the max size of the genetic map for the last marker in the map.
    
    jpeg(filename = paste0(marey.plot.output.path, chrom.name, "_marey_plot_", par.name,".jpg"), width = 1000, height = 1000)
    plot((physical.distance/1000000), genetic.distance, pch = 16, cex = 1, xlab = "Physical map position (Mb)", ylab = "Genetic map position (cM)", main = paste0(chrom.name, "_", par.name))
    points((physical.distance[183]/1000000), genetic.distance[183], pch = 16, cex = 1, col = 'blue')
    dev.off()
    cat(paste0(marey.plot.output.path, chrom.name, "_marey_plot_", par.name,".jpg", " Successfully Completed!", "\n"))
  }
}

write.table(cbind(mrkr.count, map.size), file = paste0(genetic.map.input.path, "mrkr_count_vs_genmap_size.txt"), quote = FALSE, row.names = FALSE)
stat.tab <- read.table(file = paste0(genetic.map.input.path, "mrkr_count_vs_genmap_size.txt"), header = FALSE, stringsAsFactors = FALSE)
plot(stat.tab[,1], stat.tab[,2], xlab  = "Number of markers in the genetic map", ylab = "Genetic map size (cM)")
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################