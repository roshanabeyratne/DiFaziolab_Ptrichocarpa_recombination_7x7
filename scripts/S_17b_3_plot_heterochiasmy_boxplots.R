# set working directory and load the packages used in the script
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
rm(list = ls())
library(dplyr)
library(stringr)
library(lme4)
library(MASS)
library(qvalue)
library(exactci)


# set main input and output paths
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
statistics.output.path <- "./statistics/gross_stat_plots/"
sliding.avg.output.path <- "./statistics/heterochiasmy_analysis/sliding_average/"
gen.feat.output.path <- "./genomic_features/analysis_output/"
lg.names <- c("LG-I", "LG-II", "LG-III", "LG-IV", "LG-V", "LG-VI", "LG-VII", "LG-VIII", "LG-IX", "LG-X", "LG-XI", "LG-XII", "LG-XIII", "LG-XIV",
              "LG-XV", "LG-XVI", "LG-XVII", "LG-XVIII", "LG-XIX")

######  
# chromosome box-plots
chrom.size.df <- read.csv(file = paste0(statistics.input.path, "chrom.size.df_Stetv2.csv"), header = TRUE, stringsAsFactors = FALSE)
num.ind.df <- read.csv(paste0(statistics.input.path, "number_of_individuals_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
num.COs.df <- read.csv(paste0(statistics.input.path, "number_of_crossovers_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
colnames(num.ind.df) <- colnames(num.ind.df) %>% str_replace("^X", "")
colnames(num.COs.df) <- colnames(num.COs.df) %>% str_replace("^X", "")
lg.names <- c("LG-I", "LG-II", "LG-III", "LG-IV", "LG-V", "LG-VI", "LG-VII", "LG-VIII", "LG-IX", "LG-X", "LG-XI", "LG-XII", "LG-XIII", "LG-XIV",
              "LG-XV", "LG-XVI", "LG-XVII", "LG-XVIII", "LG-XIX")

# In order to apply the linear models the data needs to be arranged in a different table format.
chrom.reshaped.co.count.df <- NULL
for (j in 1:ncol(num.COs.df)) {chrom.reshaped.co.count.df <- rbind(chrom.reshaped.co.count.df, cbind(chrom.size.df, rep(colnames(num.COs.df)[j]), num.COs.df[,j], num.ind.df[,j]))}
colnames(chrom.reshaped.co.count.df) <- c(colnames(chrom.size.df), "PAR", "CO_COUNT", "FAM_SIZE")
chrom.reshaped.co.count.df$SEX <- as.factor(rep(c("F", "M"), each = 133))
chrom.reshaped.co.count.df$CHROM <- as.factor(chrom.reshaped.co.count.df$CHROM)
# sum(chrom.reshaped.co.count.df$CO_COUNT) # 38846
chrom.reshaped.co.count.df$CO_RATE_PER_MEIOSIS <- chrom.reshaped.co.count.df$CO_COUNT / chrom.reshaped.co.count.df$FAM_SIZE

pdf(file = paste0(statistics.output.path, "genomewide_fine_scale_heterochiasmy_box.pdf"), width=4, height=16)
layout(matrix(1:19, nrow=19, byrow=T))
# layout.show(n=19)
op<-par(mar=c(2.1,2.8,2.1,0))
col.vec <- c("red", "blue")
for (chrom.name in chrom.name.vec) {
  y_lab<-ifelse(chrom.name==chrom.name.vec[1], "CO rate per meiosis", '')
  boxplot(chrom.reshaped.co.count.df$CO_RATE_PER_MEIOSIS[chrom.reshaped.co.count.df$CHROM == chrom.name] ~ chrom.reshaped.co.count.df$SEX[chrom.reshaped.co.count.df$CHROM == chrom.name],
          ylim=c(0,4),
          horizontal=T, 
          main = lg.names[chrom.name.vec==chrom.name], ylab=y_lab, xlab='',
          frame.plot=T,
          col = col.vec)
  # stripchart(chrom.reshaped.co.count.df$CO_RATE_PER_MEIOSIS[chrom.reshaped.co.count.df$CHROM == chrom.name] ~ chrom.reshaped.co.count.df$SEX[chrom.reshaped.co.count.df$CHROM == chrom.name],
  #            vertical = F, method = "jitter", add = TRUE, pch = 20, col =c('black'))
}
par(op)
dev.off()
### END 