# setting up working environment
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7")
rm(list=ls())
library(tidyverse)


jf_15mers <- read.table(file='Stettler_14_linkage_mapping/motif_analysis/jellyFish/Stet14_v2_15mer_enrichment.txt', 
                        header=F, stringsAsFactors = F, nrows = 146787282)
colnames(jf_15mers) <- c('counts', '15_mer_seq')
total_counts <- sum(jf_15mers$counts)
hist(jf_15mers$counts) 

count_cut_off <- 1
jf_15mers_cut_off_plus <- jf_15mers[jf_15mers$counts > count_cut_off, ]
nrow(jf_15mers_cut_off_plus)
hist(jf_15mers_cut_off_plus$counts, breaks=100) 

# configure MEME motif file
meme_motifs <- read.table(file='Stettler_14_linkage_mapping/motif_analysis/BLAST_analysis/Final/motif_table_p_only.txt', 
                          header=F, stringsAsFactors = F)
meme_motifs <- meme_motifs[,-2]
temp <- as.data.frame(t(sapply(meme_motifs[,1], function(x){y <- unlist(str_split(x, '-')); return(y)})), stringsAsFactors=F, row.names = F)
meme_motifs <- cbind(temp, meme_motifs[,2])
colnames(meme_motifs) <- c('qseqid_idx', 'motif', 'p-value')
str(meme_motifs)

# configure BLAST output
blast_out <- read.table(file='Stettler_14_linkage_mapping/motif_analysis/BLAST_analysis/Final/blast.out_p_only.txt', 
                          header=F, stringsAsFactors = F)
colnames(blast_out) <- c('qseqid', 'qlen', 'sseqid', 'slen', 'mismatch', 'gaps', 'length', 'pident', 'bitscore', 'evalue', 'qseq', 'sseq')
temp <- as.data.frame(t(sapply(blast_out[,3], function(x){y <- unlist(str_split(x, '_')); return(y)})), stringsAsFactors=F, row.names = F)
temp <- cbind(temp, unite(temp[,1:5], sseqid2))
temp$V7 <- as.integer(temp$V7)
str(temp)
blast_out <- data.frame(cbind(blast_out, temp[,c(8,5,7)]))
str(blast_out)

# obtain the query seq idx
temp <- as.data.frame(t(sapply(blast_out[,1], function(x){y <- unlist(str_split(x, '_')); return(y)})), stringsAsFactors=F, row.names = F)
blast_out <- data.frame(cbind(as.numeric(temp[,4]), blast_out), stringsAsFactors = F)
colnames(blast_out) <- c('qseqid_idx', 'qseqid', 'qlen', 'sseqid', 'slen', 'mismatch', 'gaps', 'length', 'pident', 'bitscore', 'evalue', 'qseq', 'sseq', 'sseqid2', 'sseq_idx', 'counts')
blast_out$counts <- as.integer(blast_out$counts)
str(blast_out)



# calculating percentiles for BLAST hits
num_jf_15mers_cut_off_plus <- nrow(jf_15mers_cut_off_plus)
percentiles <- NULL
for (i in 1:nrow(blast_out)) {
  if (blast_out$counts[i] < count_cut_off) {
    percentiles <- c(percentiles, NA)
  } else{
    perct <-(1 - (sum(jf_15mers_cut_off_plus$counts > blast_out$counts[i]) / num_jf_15mers_cut_off_plus)) *100 
    percentiles <- c(percentiles, perct)
  }
}
blast_out <- cbind(blast_out, pct=percentiles)


# 
tmp_sseq <- NULL
for (i in 1:nrow(blast_out)) {
  tmp_sseq_idx <- as.integer(blast_out$sseq_idx[i])
  tmp_sseq <- c(tmp_sseq, as.character(jf_15mers$`15_mer_seq`[tmp_sseq_idx]))
}
blast_out2 <- cbind(blast_out, sseq_full=tmp_sseq)
meme_motifs$qseqid_idx <- as.integer(meme_motifs$qseqid_idx)

# merging BLAST output information and the MEME enriched DNA sequence information
meme_motifs2 <- merge(meme_motifs, blast_out2, by = 'qseqid_idx', all.x=T)
meme_motifs4 <- meme_motifs2[,c('qseqid_idx', 'motif', 'qlen', 'p-value', 'sseq_full', 'slen', 'counts', 
                                'pct', 'mismatch', 'gaps', 'length', 'pident', 'bitscore', 'evalue')]
motif_length <- sapply(meme_motifs4$motif, nchar)
meme_motifs4$qlen <- motif_length
write.table(meme_motifs4, file='Stettler_14_linkage_mapping/motif_analysis/BLAST_analysis/Final/meme_motif_sup_tab_2.txt', quote=F, row.names=F)



# BLAST on existing motifs in publications

# collection of previous motifs
pub_motifs <- read.table(file='Stettler_14_linkage_mapping/motif_analysis/BLAST_analysis/Final/pub_motifs.txt', header=F, stringsAsFactors=F )
colnames(pub_motifs) <- c('pub_idx', 'motif', 'motif_len', 'author', 'sprecies')

# Hits to pubs
blast_out_pub_motifs <- read.table(file='Stettler_14_linkage_mapping/motif_analysis/BLAST_analysis/Final/blast.out_pub_motifs.txt', header=T, stringsAsFactors=F )
query_seq_ids <- as.integer(blast_out_pub_motifs$qseqid %>% str_replace('enriched_motif_seq_', ''))
query_seqs <- meme_motifs$motif[query_seq_ids]

pub_motif_idxs <- as.integer(gsub("[^0-9.-]", "", blast_out_pub_motifs$sseqid ))
pub <- as.character(gsub("[^a-z | A-Z]", "", blast_out_pub_motifs$sseqid ))
pub_motif_seqs <- pub_motifs$motif[pub_motif_idxs]

blast_out_pub_motifs_2 <- cbind(blast_out_pub_motifs, query_seqs, query_seq_ids, pub_motif_seqs, pub_motif_idxs)
blast_out_pub_motifs_2 <- blast_out_pub_motifs_2[,c('qseqid', 'query_seq_ids', 'query_seqs', 'qlen',
                                                    'sseqid', 'pub_motif_idxs', 'pub_motif_seqs', 'slen', 
                                                    "mismatch", "gaps", "length", "pident", "bitscore", 'evalue')]
write.table(blast_out_pub_motifs_2, file='Stettler_14_linkage_mapping/motif_analysis/BLAST_analysis/Final/blast.out_pub_motifs_aug.txt', quote=F, row.names=F)

set.seed(012)
rand_vec <- rnorm(100)
rand_cos<- cos(pi/4) * rand_vec
rand_sin <- sin(pi/4) * rand_vec
var(rand_vec)
var(rand_cos) + var(rand_sin) + (2*var(rand_cos, rand_sin))
var(rand_cos + rand_sin)
