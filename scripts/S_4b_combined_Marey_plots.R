#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: August 28, 2020
# USAGE: Input files include genetic maps and table with physical vs. genetic map sizes. The marker names contain the chromosome and marker putative physical position. 

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

# ollowing is the input file chrom.size.df with physical vs. genetic map sizes
# CHROM   LENGTH GEN_DIST
# 1  Chr01 50495398      350
# 2  Chr02 25263042      210
# 3  Chr03 25000000      200
# 4  Chr04 24267058      175
# 5  Chr05 25890711      200
# 6  Chr06 27912132      250
# 7  Chr07 15610920      120
# 8  Chr08 19465468      175
# 9  Chr09 12948749      125
# 10 Chr10 22580539      175
# 11 Chr11 18501278      125
# 12 Chr12 15760353      150
# 13 Chr13 16320724      125
# 14 Chr14 18920901      150
# 15 Chr15 15278584      125
# 16 Chr16 14494368      125
# 17 Chr17 16080365      125
# 18 Chr18 16958307      150
# 19 Chr19 15942152      150

# PURPOSE: Plot a variety of marey map styles with physical vs genetic map distance. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
library(RColorBrewer)
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
proposed.plots.output.path <- "./onemap/bin_markers/genetic_maps/stats/"
########################################################################################################################################################################################################
# Following dataframe has chromosome size in Approximated physical and genetic distance measures
chrom.size.df <- read.csv(file = "./scripts_lists/chrom.size.df_Stetv2.csv", header = TRUE, stringsAsFactors = FALSE)

# add a column with the predicted centromere location to chrom.size.df based on https://doi.org/10.3389/fgene.2019.00487 
chrom.size.df$CENTRMR.LOC <- c(18745000, 18930000, 5845000, 13245000, 14600000, 15895000, 7125000, 15740000, 810000, 5225000, 9585000, 7850000, 8955000,16610000, 6235000, 8475000, 7350500, 7125000, 7710000)
chrom.size.df$CENTRMR.LFT <- c(12550001, 16680001, 4380001, 11770001, 12960001, 12510001, 5620001, 14670001, 1, 4220001, 8390001, 6370001, 8400001, 14290001, 4480001, 7470001, 5700001, 6130001, 4970001)
chrom.size.df$CENTRMR.RIGHT <- c(24940000, 21180000, 7310000, 14720000, 16240000, 19280000, 8630000, 16800000, 1620000, 6230000, 10780000, 9330000,9510000, 18930000, 7990000, 9480000, 8990000, 8120000, 10450000) 

lg.names <- c("LG-I", "LG-II", "LG-III", "LG-IV", "LG-V", "LG-VI", "LG-VII", "LG-VIII", "LG-IX", "LG-X", "LG-XI", "LG-XII", "LG-XIII", "LG-XIV",
              "LG-XV", "LG-XVI", "LG-XVII", "LG-XVIII", "LG-XIX")
head(chrom.size.df)
# jpeg(filename = paste0(proposed.plots.output.path, "Chr01_marey_plot.jpg"), width = 480, height = 480, quality = 150)

# jpeg(filename = paste0(proposed.plots.output.path, "combined_marey_plot.jpg"), width = 2000, height = 1000, quality = 150)
# pdf(file = paste0(proposed.plots.output.path, "Chr01_marey_plot.pdf"), width = 16, height = 12)
pdf(file = paste0(proposed.plots.output.path, "combined_marey_plot2.pdf"), width = 21, height = 25)
colors <- c(brewer.pal(7, "Set1"), brewer.pal(7, "Dark2"))
op <- par(mfrow = c(5,4), mar=c(5.1, 4.1, 4.1, 2.1), cex.lab=1.75, cex.axis=1.5, cex.sub=1.0, cex.main=2)
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  
  chrom.length <- chrom.size.df$LENGTH[which(chrom.size.df$CHROM == chrom.name)]
  gen.map.size <-  chrom.size.df$GEN_DIST[which(chrom.size.df$CHROM == chrom.name)]
  plot((seq(from = 1, to = chrom.length, length.out = 10)/1000000), seq(from = 1, to = gen.map.size, length.out = 10), pch = 16, cex = 0.5, col = 'white', 
       xlab = "Physical map position (Mb)",
       ylab = "Genetic map position (cM)",
       # xlab = NULL, 
       # ylab = NULL, 
       main = lg.names[which(chrom.name.vec %in% chrom.name)])
  
  upstream.margin <- (chrom.size.df$CENTRMR.LFT[chrom.size.df$CHROM == chrom.name])/1000000
  downstream.margin <- (chrom.size.df$CENTRMR.RIGHT[chrom.size.df$CHROM == chrom.name])/1000000
  rect(xleft = upstream.margin, xright = downstream.margin, ybottom = 0, ytop = gen.map.size, density = 1000, col = "light grey")
  centrmr.loc <- chrom.size.df$CENTRMR.LOC[chrom.size.df$CHROM == chrom.name]
  abline(v = centrmr.loc/1000000)

  
  for (par.name in par.name.vec) { # par.name <- "6909"
    
    if (file.exists(paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"))) {
      gen.map <- read.table(file = paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    } else {
      gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    }
    
    string <- paste0(chrom.name, "_")
    physical.distance <- gen.map[,2] %>% str_replace(string, "") %>% as.numeric()
    genetic.distance <- gen.map[,3] %>% as.numeric()
    
    
    # jpeg(filename = paste0(marey.plot.output.path, chrom.name, "_marey_plot_", par.name,".jpg"), width = 1000, height = 1000)
    if (any(par.name %in% par.name.vec[1:7])) {
      # lines((physical.distance/1000000), genetic.distance, pch = 16, cex = 0.5, col = 'red')
      # points((physical.distance/1000000), genetic.distance, pch = 16, cex = 0.5, col = 'red')
      # lines((physical.distance/1000000), genetic.distance, pch = 16, cex = 0.5, col = colors[which(par.name.vec == par.name)])
      points((physical.distance/1000000), genetic.distance, pch = 16, cex = 0.5, col = colors[which(par.name.vec == par.name)])
    } else {
      # lines((physical.distance/1000000), genetic.distance, pch = 16, cex = 0.5, col = 'blue')
      # points((physical.distance/1000000), genetic.distance, pch = 16, cex = 0.5, col = 'blue')
      # lines((physical.distance/1000000), genetic.distance, pch = 16, cex = 0.5, col = colors[which(par.name.vec == par.name)])
      points((physical.distance/1000000), genetic.distance, pch = 16, cex = 0.5, col = colors[which(par.name.vec == par.name)])
    }
  
    # dev.off()
    # cat(paste0(marey.plot.output.path, chrom.name, "_marey_plot_", par.name,".jpg", " Successfully Completed!", "\n"))
  }
}
# plot((seq(from = 1, to = chrom.length, length.out = 10)/1000000), seq(from = 1, to = gen.map.size, length.out = 10), pch = 16, cex = 0.5, col = 'white', 
#      xlab = NULL, ylab = NULL, main = NULL, axes=FALSE, ann=FALSE)
# legend(5, 150, legend=par.name.vec, col=colors, pch=16, cex=0.9, bty='n')

plot(1:10, 1:10, pch = 16, cex = 0.5, col = 'white', 
     xlab = NULL, ylab = NULL, main = NULL, axes=FALSE, ann=FALSE)
x_y_coords <- cbind(y=seq(from=2, to=8, by=1), x=rep(c(3, 7), each=7))
#points(x_y_coords[,2], x_y_coords[,1], pch=16, col=colors)
nudge=0.25
text(x_y_coords[,2]+nudge, x_y_coords[,1], labels= par.name.vec, col=colors, cex=2)

par(op)
dev.off()


########################################################################################################################################################################################################
# Create an inmage that shows all chromosomes in one plot for each parent
# NOTE TO SELF: to complete this code write a block that corrects for the increase in genetic map
for (par.name in par.name.vec[10]) {
  # par.name <- "2515"
  # jpeg(file = paste0("./onemap/all_chrom_plots/",par.name,"_all_chromosomes.jpg"), width = 10, height = 10, unit = "in", res = 600)
  # tiff(file = paste0(v,"_HS_for_all_chromosomes","v2.tiff"),units="in", width=10, height=10, compression = 'lzw', res = 300)
  par(mai = c(0,0,0,0))
  #plot(1:15,ty="n",axes=F,xlab="",ylab="")
  #text(6,20,"\\VE",vfont=c("sans serif","bold"),cex=2)
  #text(16,20,"\\MA",vfont=c("sans serif","bold"),cex=2)
  par(mfcol = c(10,2), oma = c(5.0, 4.5, 4.5, 4.5), mai = c(0.05,0.05,0.05,0.05))
  par(mfg = c(1,1))
  par(cex.axis = 1.0)
  # par(mar = c(1, 1, 1, 1),  mai=c(0.3,0.3,0.3,0.3))#, oma = c(4, 4, 0.5, 0.5))
  
  for (chrom.name in chrom.name.vec[1:9]) {
    if (file.exists(paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"))) {
      gen.map <- read.table(file = paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    } else {
      gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    }
    
    str <- paste0(chrom.name,"_")
    gen.map$V2 <- as.numeric(gen.map$V2 %>% stringr::str_replace(pattern = str, replacement = ""))/1000000
    gen.map$V3 <- as.numeric(gen.map$V3)
    plot(gen.map$V2,gen.map$V3, main = NULL, ylim = c(0,350), xaxt = "n", yaxt = "n", type = "p", pch = 16,cex = 0.5, col = "black", xlim = c(0,50))
    axis(2, at = c(0,150,300), labels = c(0, 150, 300))
    text(2,300, chrom.name, cex  = 2.0)
    if (chrom.name == "Chr05") {
      mtext("Genetic Distance (cM)", 2, 3, cex = 1.5)
    }
  }
  #dev.off()
  
  for (chrom.name in chrom.name.vec[10]) {
    
    if (file.exists(paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"))) {
      gen.map <- read.table(file = paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    } else {
      gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    }
    
    str <- paste0(chrom.name,"_")
    gen.map$V2 <- as.numeric(gen.map$V2 %>% stringr::str_replace(pattern = str, replacement = ""))/1000000
    gen.map$V3 <- as.numeric(gen.map$V3)
    plot(gen.map$V2, gen.map$V3, main = NULL, ylim = c(0,350), xaxt = "n", yaxt = "n", type = "p", pch = 16, cex = 0.5, col = "black", xlim = c(0,50))
    text(2, 300, chrom.name, cex = 2.0)
    axis(1, at = (0:5)*10, labels = (0:5)*10)
    axis(2, at = c(0,150,300), labels = c(0,150,300))
    mtext("Physical position (million bp)", 1, 3, at = 50 ,cex = 1.5)
  }
  
  for (chrom.name in chrom.name.vec[11:18]) {
    
    if (file.exists(paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"))) {
      gen.map <- read.table(file = paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    } else {
      gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    }
    
    str <- paste0(chrom.name,"_")
    gen.map$V2 <- as.numeric(gen.map$V2 %>% stringr::str_replace(pattern = str, replacement = ""))/1000000
    gen.map$V3 <- as.numeric(gen.map$V3)
    plot(gen.map$V2,gen.map$V3, main = NULL, ylim = c(0,350), xaxt = "n", yaxt = "n", type = "p", pch = 16, cex = 0.5, col = "black", xlim = c(0,50))
    axis(4, at = c(0,150,300), labels = c(0,150,300))
    text(2, 300, chrom.name, cex = 2.0)
  }
  
  for (chrom.name in chrom.name.vec[19]) {
   
    if (file.exists(paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"))) {
      gen.map <- read.table(file = paste0(rippled.genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    } else {
      gen.map <- read.table(file = paste0(genetic.map.input.path, chrom.name, "_", par.name, "_genetic.map"), header = FALSE, stringsAsFactors = FALSE)
    }
    
    str <- paste0(chrom.name,"_")
    gen.map$V2 <- as.numeric(gen.map$V2 %>% stringr::str_replace(pattern = str, replacement = ""))/1000000
    gen.map$V3 <- as.numeric(gen.map$V3)
    plot(gen.map$V2,gen.map$V3, main = NULL, ylim = c(0,350), xaxt = "n", yaxt = "n", type = "p", pch = 16, cex = 0.5, col = "black", xlim = c(0,50))
    text(2, 300, chrom.name, cex = 2.0)
    axis(1, at = (0:5)*10, labels = (0:5)*10)
    axis(4, at = c(0,150,300), labels = c(0,150,300))
  }
  # dev.off()
}
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################