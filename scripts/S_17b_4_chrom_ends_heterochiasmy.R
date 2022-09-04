#!/usr/bin/R

# Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: February 9, 2020
# Purpose: Test the hypothesis that recombinaiton rate is higher in males closer to the chromosome ends

# Set the environment and read inputs
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
rm(list = ls())
library(dplyr)
library(stringr)



par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
input.path <- "./statistics/"
statistics.output.path <- "./statistics/gross_stat_plots/"

chrom.size.df <- read.csv(file = paste0(input.path, "chrom.size.df_Stetv2.csv"), header = TRUE, stringsAsFactors = FALSE)


# CO.df <- read.csv(paste0(input.path, "./co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# # any((CO.df$CO_SIZE == 0))
# # Define chromosome ends as 1Mb upstream and 1Mb downstream
# block_size <- 1000000
# tmp_vec <- NULL
# for (i in 1:nrow(chrom.size.df)) {
#   upstrm_block <- c(1, block_size)
#   dwnstrm_block <- c(chrom.size.df$LENGTH[i]-block_size, chrom.size.df$LENGTH[i])
#   tmp_vec <- rbind(tmp_vec, c(upstrm_block, dwnstrm_block))
# }
# chrom.size.df2 <- cbind(chrom.size.df, tmp_vec)
# colnames(chrom.size.df2)[4:7] <- c('up_start', 'up_end', 'dwn_start', 'dwn_end')


# Read-in CO count data for 7x7 cross
binned.cos.by.par.df <- read.csv(file = paste0(input.path, "co.bin_size_960kb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
colnames(binned.cos.by.par.df) <- colnames(binned.cos.by.par.df) %>% str_replace("^X", "")
head(binned.cos.by.par.df)

# Assign window indces by chromosome
win.idx <- NULL
for (chrom_name in chrom.name.vec) {
  tmp <- binned.cos.by.par.df[binned.cos.by.par.df$CHROM==chrom_name, ]
  win.idx <- c(win.idx, 1:nrow(tmp))
}
binned.cos.by.par.df <- cbind(binned.cos.by.par.df, WIN.IDX=win.idx)
head(binned.cos.by.par.df)

# Evaluate whether there are differences in CO counts between males and females
# at the ends of chromsomes
range <- 5 # assigning this many consecutive windows as belonging to the distal end 
p_values <- NULL
for (chrom_name in chrom.name.vec) {
  subset_COs <- binned.cos.by.par.df[binned.cos.by.par.df$CHROM==chrom_name, ]
  upstrm_COs <- subset_COs[1:range, 4:17] %>% colSums()
  dwnstrm_COs <- subset_COs[(nrow(subset_COs)-range):nrow(subset_COs), 4:17] %>% colSums()
  
  upstrm_test <- t.test(x=upstrm_COs[1:7], y=upstrm_COs[8:14], mu=0)
  dwnstrm_test <- t.test(x=dwnstrm_COs[1:7], y=dwnstrm_COs[8:14], mu=0)
  p_values <- rbind(p_values, c(chrom_name, upstrm_test$p.value, dwnstrm_test$p.value))
  
  # upstrm_test <- wilcox.test(x=upstrm_COs[1:7], y=upstrm_COs[8:14], mu=0)
  # dwnstrm_test <- wilcox.test(x=dwnstrm_COs[1:7], y=dwnstrm_COs[8:14], mu=0)
  # p_values <- rbind(p_values, c(chrom_name, upstrm_test$p.value, dwnstrm_test$p.value))
}
p_values <- as.data.frame(p_values, stringsAsFactors = F)
colnames(p_values) <- c('chrom_name', 'upstream_end_p_value', 'dwnstream_end_p_value')
p_values$upstream_end_p_value <- as.numeric(p_values$upstream_end_p_value)
p_values$dwnstream_end_p_value <- as.numeric(p_values$dwnstream_end_p_value)
str(p_values)
print(p_values)

# Check if the window idx can predict the COs
binned.cos.by.par.df <- cbind(binned.cos.by.par.df, cum_CO_fem=apply(binned.cos.by.par.df[,4:10], 1, sum))
binned.cos.by.par.df <- cbind(binned.cos.by.par.df, cum_CO_male=apply(binned.cos.by.par.df[,11:17], 1, sum))
str(binned.cos.by.par.df)

lm_mod_fem <- lm(cum_CO_fem~WIN.IDX, data=binned.cos.by.par.df)
summary(lm_mod_fem)
lm_mod_male <- lm(cum_CO_male~WIN.IDX, data=binned.cos.by.par.df)
summary(lm_mod_male)

plot(binned.cos.by.par.df$WIN.IDX, binned.cos.by.par.df$cum_CO_fem, xlab='window idx within chromosome', ylab='cumulative CO count (female)')
plot(binned.cos.by.par.df$WIN.IDX, binned.cos.by.par.df$cum_CO_male, xlab='window idx within chromosome', ylab='cumulative CO count (male)')



# read-in the centromere locations from Nisqually
# add a column with the predicted centromere location to chrom.size.df based on https://doi.org/10.3389/fgene.2019.00487 
chrom.size.df$CENTRMR.LOC <- c(18745000, 18930000, 5845000, 13245000, 14600000, 15895000, 7125000, 15740000, 810000, 5225000, 9585000, 7850000, 8955000,16610000, 6235000, 8475000, 7350500, 7125000, 7710000)
chrom.size.df$CENTRMR.LFT <- c(12550001, 16680001, 4380001, 11770001, 12960001, 12510001, 5620001, 14670001, 1, 4220001, 8390001, 6370001, 8400001, 14290001, 4480001, 7470001, 5700001, 6130001, 4970001)
chrom.size.df$CENTRMR.RIGHT <- c(24940000, 21180000, 7310000, 14720000, 16240000, 19280000, 8630000, 16800000, 1620000, 6230000, 10780000, 9330000,9510000, 18930000, 7990000, 9480000, 8990000, 8120000, 10450000) 
head(chrom.size.df)

# add a column with distance to centromere
cntrmr.loc <- NULL
for (chrom_name in chrom.name.vec) {
  tmp <- binned.cos.by.par.df[binned.cos.by.par.df$CHROM==chrom_name, ]
  cntrmr.loc <- c(cntrmr.loc, rep(chrom.size.df$CENTRMR.LOC[chrom.size.df$CHROM==chrom_name], nrow(tmp)))
}
binned.cos.by.par.df <- cbind(binned.cos.by.par.df, dist_to_cntmr=abs(binned.cos.by.par.df$START-cntrmr.loc))
str(binned.cos.by.par.df)
head(binned.cos.by.par.df)

# Evaluate whether distance to centromere has an effect
lm_mod_fem2 <- lm(cum_CO_fem~dist_to_cntmr, data=binned.cos.by.par.df)
summary(lm_mod_fem2)

lm_mod_male2 <- lm(cum_CO_male~dist_to_cntmr, data=binned.cos.by.par.df)
summary(lm_mod_male2)



# Evaluate differential effect of distance to centromere between males and female
# new dataframe
colnames(binned.cos.by.par.df)[19:20] <- 'cum_COs'
tmp <- rbind(binned.cos.by.par.df[,c(1:3, 18:19, 21)], binned.cos.by.par.df[,c(1:3, 18, 20:21)])
tmp <- cbind(tmp, sex=rep(c('F', 'M'), each=nrow(binned.cos.by.par.df)))
str(tmp)
head(tmp)

lm_mod3 <- lm(cum_COs~sex*dist_to_cntmr, data=tmp)
lm_mod4 <- lm(cum_COs~sex+dist_to_cntmr, data=tmp)
lm_mod5 <- lm(cum_COs~dist_to_cntmr, data=tmp)
lm_mod6 <- lm(cum_COs~sex, data=tmp)

AIC(lm_mod3, lm_mod4, lm_mod5, lm_mod6)

summary(lm_mod3)
summary(lm_mod4)
summary(lm_mod5)
summary(lm_mod6)



# checking whether there is a significant CO count difference between sexes near the centromere
# Identify windows that are 1Mb from the centromere
range_size <- 1000000
tmp2 <- tmp[tmp$dist_to_cntmr <= range_size, ]
head(tmp2)

lm_mod7 <- lm(cum_COs~sex*dist_to_cntmr, data=tmp2)
summary(lm_mod7)

# Conducting additional post-hoc tests 
t.test(tmp2$cum_COs[tmp2$sex=='M'], tmp2$cum_COs[tmp2$sex=='F'], mu=0)
wilcox.test(tmp2$cum_COs[tmp2$sex=='M'], tmp2$cum_COs[tmp2$sex=='F'], mu=0)

