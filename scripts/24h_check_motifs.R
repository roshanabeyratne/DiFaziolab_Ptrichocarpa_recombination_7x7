### METADATA: 
# Written by CRA
# March 17, 2022
# Check for enriched motifs from output files downloaded from spruce
### 


# setting up working environment
setwd("/Volumes/GoogleDrive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7")
rm(list=ls())
library(tidyverse)
input_path<- "./Stettler_14_linkage_mapping/motif_analysis/comp_to_background/"
input_dir_seq <- c('le_10kb', paste0('le_10kb', 2:10))
input_dir_seq2 <- paste0('le_10kb_neg_neg', 3:10)

#

# read-in the input motif tables to compare
le_10kb1_df <- read.csv(file=paste0(input_path, input_dir_seq[1], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb2_df <- read.csv(file=paste0(input_path, input_dir_seq[2], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb3_df <- read.csv(file=paste0(input_path, input_dir_seq[3], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb4_df <- read.csv(file=paste0(input_path, input_dir_seq[4], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb5_df <- read.csv(file=paste0(input_path, input_dir_seq[5], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb6_df <- read.csv(file=paste0(input_path, input_dir_seq[6], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb7_df <- read.csv(file=paste0(input_path, input_dir_seq[7], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb8_df <- read.csv(file=paste0(input_path, input_dir_seq[8], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb9_df <- read.csv(file=paste0(input_path, input_dir_seq[9], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb10_df <- read.csv(file=paste0(input_path, input_dir_seq[10], '/motif_table.txt'), header=F, stringsAsFactors=F)


# read-in the negative-negative motif tables to compare
le_10kb3_neg_neg_df <- read.csv(file=paste0(input_path, input_dir_seq2[1], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb4_neg_neg_df <- read.csv(file=paste0(input_path, input_dir_seq2[2], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb5_neg_neg_df <- read.csv(file=paste0(input_path, input_dir_seq2[3], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb6_neg_neg_df <- read.csv(file=paste0(input_path, input_dir_seq2[4], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb7_neg_neg_df <- read.csv(file=paste0(input_path, input_dir_seq2[5], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb8_neg_neg_df <- read.csv(file=paste0(input_path, input_dir_seq2[6], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb9_neg_neg_df <- read.csv(file=paste0(input_path, input_dir_seq2[7], '/motif_table.txt'), header=F, stringsAsFactors=F)
le_10kb10_neg_neg_df <- read.csv(file=paste0(input_path, input_dir_seq2[8], '/motif_table.txt'), header=F, stringsAsFactors=F)



# Change the above read tables to proper columns
for (i in 1:10) {
  temp <- get(paste0('le_10kb', i, '_df'))
  tmp1 <- as.data.frame(t(apply(temp, 1, function(x){x<- as.character(x); y<- x %>% str_split("\t") %>% unlist(); return(y)})), stringsAsFactors=F)
  tmp2 <- as.data.frame(t(apply(tmp1, 1, function(x){x<- as.character(x[1]); y<- x %>% str_split("-") %>% unlist(); return(y)})), stringsAsFactors=F)
  tmp_tab <- cbind(tmp2, tmp1[,2])
  colnames(tmp_tab) <- c('idx', 'motif', 'p_val')
  assign(paste0('le_10kb', i, '_df'), tmp_tab)
}

# Change the above read tables to proper columns
for (i in 3:10) {
  temp <- get(paste0('le_10kb', i, '_neg_neg_df'))
  tmp1 <- as.data.frame(t(apply(temp, 1, function(x){x<- as.character(x); y<- x %>% str_split("\t") %>% unlist(); return(y)})), stringsAsFactors=F)
  tmp2 <- as.data.frame(t(apply(tmp1, 1, function(x){x<- as.character(x[1]); y<- x %>% str_split("-") %>% unlist(); return(y)})), stringsAsFactors=F)
  tmp_tab <- cbind(tmp2, tmp1[,2])
  colnames(tmp_tab) <- c('idx', 'motif', 'p_val')
  assign(paste0('le_10kb', i, '_neg_neg_df'), tmp_tab)
}


# trying to implement a fuzzy match.













tmp_tab <- le_10kb3_df
tmp_tab[which(tmp_tab[,2] %in% le_10kb1_df[,2]),]
tmp_tab[which(tmp_tab[,2] %in% le_10kb2_df[,2]),]
tmp_tab[which(tmp_tab[,2] %in% le_10kb4_df[,2]),]
tmp_tab[which(tmp_tab[,2] %in% le_10kb5_df[,2]),]
tmp_tab[which(tmp_tab[,2] %in% le_10kb6_df[,2]),]
tmp_tab[which(tmp_tab[,2] %in% le_10kb7_df[,2]),]
tmp_tab[which(tmp_tab[,2] %in% le_10kb8_df[,2]),]
tmp_tab[which(tmp_tab[,2] %in% le_10kb9_df[,2]),]
tmp_tab[which(tmp_tab[,2] %in% le_10kb10_df[,2]),]


