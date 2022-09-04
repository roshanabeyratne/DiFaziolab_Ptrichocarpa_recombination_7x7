#!/usr/bin/R
########################################################################################################################################################################################################
# PURPOSE: Compare pairwise LD estimates based on 220 individuals in the core natural population of 7x7 parents 
# and cumulative CO counts estimated based on the 7x7 mapping population in 1Mb windows. 
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: June 5, 2021
# USAGE: Input files as follows 

#  head(CO.df)
# CHROM  PAR   IND CO_SIZE DWNSTRM_FLNK UPSTRM_FLNK
# Chr01 1863 24708  103140      4675047     4778187
# Chr01 1863 24708   30914     26411972    26442886
# Chr01 1863 24708   61845     34445658    34507503
# Chr01 1863 24708    5612     46505278    46510890
# Chr01 1863 24709  164382       726848      891230
# Chr01 1863 24709  153359      2513166     2666525

# head(chrom.size.df)
# CHROM   LENGTH GEN_DIST
# Chr01 50495398      350
# Chr02 25263042      210
# Chr03 25000000      200
# Chr04 24267058      175
# Chr05 25890711      200
# Chr06 27912132      250

# head(ld.tab)
# CHR POS1  POS2 N_INDV      R.2
# Chr01 7940  8325    205 0.903499
# Chr01 7940 12928    211 0.924399
# Chr01 7940 13875    211 0.936510
# Chr01 7940 14436    211 0.269450
# Chr01 7940 14952    209 0.411868
# Chr01 7940 16690    211 0.391945

########################################################################################################################################################################################################
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
rm(list = ls())
library(dplyr)
library(stringr)
library(ggplot2)
source(file = "./new_scripts/CreateBinSpanCoordinates_function.R")

########################################################################################################################################################################################################
# set main input and output paths
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
statistics.input.path <- "./statistics/"
# statistics.output.path <- "./statistics/hotspot_analysis/250kb/"
statistics.output.path <- "./statistics/hotspot_analysis/30kb/"
input.path.genomic.features <- "./genomic_features/output/"

########################################################################################################################################################################################################
# The dataframe containing the respective chromosome sizes for Nisq_V3. 
# chrom.size.df <- read.csv(file = paste0(statistics.input.path, "chrom.size.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# # modifying coordinates of Chr19 for Stettler
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr11"] <- 19000000
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr16"] <- 15494368
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr19"] <- 16500000

# new chromosome size obtained from Ran's folder for Stet-14_V2
# Chr01    49678573 
# Chr02    25309097 
# Chr03    23634591 
# Chr04    23180857 
# Chr05    23930588 
# Chr06    26840282 
# Chr07    14830512 
# Chr08    20251161 
# Chr09    12986220 
# Chr10    22611514 
# Chr11    18762988 
# Chr12    15038488 
# Chr13    15829923 
# Chr14    18880830 
# Chr15    15059279 
# Chr16    14646863 
# Chr17    16018856 
# Chr18    15951547 
# Chr19    16527936 
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr01"] <- 49678573
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr02"] <- 25309097
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr03"] <- 23634591
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr04"] <- 23180857
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr05"] <- 23930588
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr06"] <- 26840282
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr07"] <- 14830512
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr08"] <- 20251161
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr09"] <- 12986220
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr10"] <- 22611514
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr11"] <- 18762988
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr12"] <- 15038488
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr13"] <- 15829923
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr14"] <- 18880830
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr15"] <- 15059279
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr16"] <- 14646863
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr17"] <- 16018856
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr18"] <- 15951547
# chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr19"] <- 16527936
# write.csv(chrom.size.df, file = paste0(statistics.input.path, "chrom.size.df_Stetv2.csv"), quote=FALSE,  row.names = FALSE)
chrom.size.df <- read.csv(file = paste0(statistics.input.path, "chrom.size.df_Stetv2.csv"), header = TRUE, stringsAsFactors = FALSE)

# The dataframe with detailed cross-over events produced as per script 22b_...
CO.df <- read.csv(paste0(statistics.input.path, "co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# sum(is.na(CO.df$CO_SIZE))
# sum(CO.df$CO_SIZE == 0)
# sum(CO.df$CO_SIZE != 0)
# sum(CO.df$CO_SIZE == 0) + sum(CO.df$CO_SIZE != 0) == nrow(CO.df)

# check whether same set of offspring are present for both male and female groups
# all(unique(CO.df$IND[CO.df$CHROM=="Chr01" & CO.df$PAR %in% par.name.vec[1:7]]) %in% unique(CO.df$IND[CO.df$CHROM=="Chr01" & CO.df$PAR %in% par.name.vec[8:14]])) 

# Assign window start and end points based on window size
# win.size <- 1000000
# win.size <- 250000

# win.size <- 3840000 # (scale = 7)
# win.size <- 1920000 # (scale = 6)
win.size <- 960000 # (scale = 5)
# win.size <- 480000 # (scale = 4)
# win.size <- 240000 # (scale = 3)
# win.size <- 120000 # (scale = 2)
# win.size <- 60000 # (scale = 1)
# win.size <- 30000
# win.size <- 20000 
# win.size <- 10000 
# total.bin.span.coord.df <- NULL
# for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr01"
#   bin.span.coord.df.temp <- CreateBinSpanCoordinates(chrom.name = chrom.name, bin.size = win.size, chrom.length = chrom.size.df[chrom.size.df$CHROM == chrom.name, "LENGTH"], shift.size = win.size)
#   bin.span.coord.df.temp2 <- cbind(rep(chrom.name, nrow(bin.span.coord.df.temp)), bin.span.coord.df.temp)
#   total.bin.span.coord.df <- rbind(total.bin.span.coord.df, bin.span.coord.df.temp2)
# }
# total.bin.span.coord.df <- as.data.frame(total.bin.span.coord.df, stringsAsFactors = FALSE)
# total.bin.span.coord.df[,2:3] <- apply(total.bin.span.coord.df[,2:3], 2, as.numeric)
# colnames(total.bin.span.coord.df) <- c("CHROM", "START", "END")


# Follwoing bed files are written for use with the intention of extracting genomic features from Stet-14 in a later script
# write.table(total.bin.span.coord.df, file = paste0(statistics.input.path, "30kb_window.bed"), quote = FALSE, row.names = FALSE, col.names=FALSE )
# write.table(total.bin.span.coord.df, file = paste0(statistics.input.path, "250kb_window.bed"), quote = FALSE, row.names = FALSE, col.names=FALSE )
# write.table(total.bin.span.coord.df, file = paste0(statistics.input.path, "1Mb_window.bed"), quote = FALSE, row.names = FALSE, col.names=FALSE )
# write.table(total.bin.span.coord.df, file = paste0(statistics.input.path, "960kb_window.bed"), quote = FALSE, row.names = FALSE, col.names=FALSE )


# Assign the COs to windows based on the the CO-region overlap. Creating a dataframe with CO counts in each overlapping windows.
# Each CO was assigned to a single window. This window was selected as the midpoint of the CO-region.
# CO.df <- read.csv(paste0(statistics.input.path, "co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# CO.df <- CO.df[!is.na(CO.df$DWNSTRM_FLNK) & CO.df$CO_SIZE <= win.size, ]


# all.par.df <- NULL
# for (par.name in par.name.vec) { # par.name <- par.name.vec[9]
#   cat(date(), "\n")
#   all.chroms.vec <- NULL
#   for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[18]
# 
#     tmp.co.tab <- CO.df[(CO.df$CHROM == chrom.name & CO.df$PAR == par.name), ] # subset to the half-sib fam
#     tmp.co.tab <- tmp.co.tab[tmp.co.tab$CO_SIZE != 0, ] # removing individuals that do not have an observed CO
# 
#     # if (nrow(tmp.co.tab) == 0) { # test code to capture entries with no data
#     #   stop()
#     # }
# 
#     tmp.coord.tab <- total.bin.span.coord.df[total.bin.span.coord.df$CHROM == chrom.name, ]  # subset the windows for the chromosome
#     # tmp.vec <- rep(0, nrow(total.bin.span.coord.df))
# 
#     sum.vec <- rep(0, nrow(tmp.coord.tab))
#     if (nrow(tmp.co.tab) > 0) {
#       for (i in 1:nrow(tmp.co.tab)) { # i <- 6
#         bp.coords <- c(tmp.co.tab[i,"DWNSTRM_FLNK"], tmp.co.tab[i,"UPSTRM_FLNK"])
#         co.region <- seq(from = (bp.coords[1] + 1), to = (bp.coords[2] - 1), by = 1)
# 
#         # selecting the mid-point of the CO
#         if ((length(co.region) %% 2) == 0) {
#           cand.idx <- c( (length(co.region) %/% 2), ((length(co.region) %/% 2) + 1) )
#           random.idx <- sample(cand.idx, 1) # randomly selecting an index out of the two possible.
#           mid.point <- co.region[random.idx]
#         } else {
#           mid.point <- co.region[ceiling((length(co.region) / 2))]
#         }
# 
#         # Following code block is if CO region is assigned to more than one window without considering the midpoint
#         # tmp.vec <- apply(tmp.coord.tab[2:3], 1, function(x, y = bp.coords){
#         #   test.1 <- ifelse(((x[1] > y[1] & x[1] < y[2]) | (x[2] >= y[1] & x[2] < y[2])), 1, 0)
#         #   test.2 <- ifelse((y[1] > x[1] & y[2] < x[2]) & (y[2] > x[1] & y[2] < x[2]) , 1, 0)
#         #   test <- test.1 + test.2
#         #   return(test)
#         # })
# 
#         # Considering the CO-midpoint assign the CO to a window
#         tmp.vec <- apply(tmp.coord.tab[2:3], 1, function(x, y = mid.point){
#           x <- as.numeric(x)
#           test <- ifelse(((y >= x[1]) & (y <= x[2])), 1, 0)
#           return(test)
#         })
# 
#         if (sum(tmp.vec) == 1) {
#           sum.vec <- sum.vec + tmp.vec
#         } else {
#           stop("ERROR with tmp.vec!")
#         }
#       }
#     } else {
#       sum.vec <- rep(0, nrow(tmp.coord.tab))
#     }
# 
#     all.chroms.vec <- c(all.chroms.vec, sum.vec)
#   }
# 
#   cat(date(), "\n")
#   all.par.df <- cbind(all.par.df, all.chroms.vec)
#   cat(par.name, " done!", "\n")
# }
# colnames(all.par.df) <- par.name.vec
# all.par.df.2 <- cbind(total.bin.span.coord.df, all.par.df)
# sum(apply(all.par.df.2[-(1:3)], 1, sum)) # Just making sure that the total num of CO counts match after the processing above into windows; should be 38,846

# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_10kb.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_30kb.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_250kb.csv"), quote = FALSE, row.names = FALSE)
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_10kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_30kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_250kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# 
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_3.84Mb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_1.92Mb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_960kb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_480kb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_240kb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_120kb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_60kb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_30kb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_250kb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_1Mb_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_30kb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_960kb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_1Mb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )

colnames(all.par.df.2) <- colnames(all.par.df.2) %>% str_replace("^X", "") # or # colnames(all.par.df.2)[4:17] <- par.name.vec
all.par.df.2[(all.par.df.2$CHROM == 'Chr01' & all.par.df.2$START > 29054288), '6909'] <- NA
sum(apply(all.par.df.2[-(1:3)], 1, sum, na.rm=TRUE), na.rm=TRUE)
hist(apply(all.par.df.2[-(1:3)], 1, sum), breaks=100, xlab='cumulative COs per window', main='Histogram: CO count')

# # CO distribution landscape per chromosome
# LD_v_CO_df <- NULL
# for (chrom.name in chrom.name.vec) {
#   # chrom.name <- chrom.name.vec[1]
#   cat(chrom.name, "\n")
#   ld.tab <- read.table(file=paste0("./LD_analysis/", chrom.name, "_ld-window-bp_10000.geno.ld"), header=T, stringsAsFactors = F) # pairwise LD file per chromosome
#   subset.tmp <- all.par.df.2[all.par.df.2$CHROM == chrom.name, ]
#   # head(ld.tab)
#   
#   mean_ld <- NULL
#   for (i in 1:nrow(subset.tmp)) {
#     start_pos <- subset.tmp$START[i]
#     end_pos <- subset.tmp$END[i]
#     temp_ld_tab <- ld.tab[(ld.tab$POS1 >= start_pos & ld.tab$POS1 <= end_pos) & (ld.tab$POS2 >= start_pos & ld.tab$POS1 <= end_pos), ]
#     mean_ld <- c(mean_ld, mean(temp_ld_tab$R.2, na.rm = T))
#   }
#   subset.tmp <- cbind(subset.tmp, cum_CO=apply(subset.tmp[-(1:3)], 1, sum, na.rm=TRUE), pw_LD=mean_ld)
#   LD_v_CO_df <- rbind(LD_v_CO_df, subset.tmp)
# }
# write.table(LD_v_CO_df, file='./LD_analysis/LD_v_CO_df.txt', quote=F, row.names=F)
LD_v_CO_df <- read.table(file='./LD_analysis/LD_v_CO_df.txt', header=T, stringsAsFactors=F)
colnames(LD_v_CO_df) <- colnames(LD_v_CO_df) %>% str_replace("^X", '')


# # permutation test for correlation parameter
# # Permutations are indifferent to chromsomes and carries out a genomewide shuffling
# pmuted_cum_CO <- LD_v_CO_df[,1:3]
# for (i in 1:9999) {
#   tmp_pmut <- sample(LD_v_CO_df$cum_CO, length(LD_v_CO_df$cum_CO), replace=F)
#   pmuted_cum_CO <- cbind(pmuted_cum_CO, tmp_pmut)
#   cat(i, "\n")
# }
# colnames(pmuted_cum_CO)[4:10002] <- paste0('pmut_CO_', 1:9999)
# write.table(pmuted_cum_CO, file='./LD_analysis/pmuted_cum_CO.txt', quote=F, row.names=F)
pmuted_cum_CO <- read.table(file='./LD_analysis/pmuted_cum_CO.txt', header=T, stringsAsFactors=F)


# Check for correlation
temp_cor <- cor.test(LD_v_CO_df$cum_CO, LD_v_CO_df$pw_LD, na.rm=T, method='spearman')
stat_cor <-  temp_cor$estimate
permut_cor_vec <- c(stat_cor)
for (j in 4:10002) {
  temp_cor <- cor.test(pmuted_cum_CO[,j], LD_v_CO_df$pw_LD, na.rm=T, method='spearman')
  permut_cor_vec <- c(permut_cor_vec, temp_cor$estimate)
  cat(j, "\n")
}
# sum(is.na(permut_cor_vec))
hist(permut_cor_vec, main='Histogram the permuted dataset', xlab='correlation coefficient (r-Pearson)')
# One-tailed test for LD and CO-count correlation
p_value_cor <- sum(abs(permut_cor_vec) >= abs(stat_cor)) / length(permut_cor_vec)


# plot cumulative-CO vs. pairwise LD
dev.new(width=10, height=8)
op <- par(mfrow=c(4,5))
for (chrom.name in chrom.name.vec) {
  subset.tmp <- LD_v_CO_df[LD_v_CO_df$CHROM == chrom.name,]
  plot(subset.tmp$START, subset.tmp$cum_CO, type='l', lwd=2, col='red', 
       # xlab=paste(chrom.name, 'win.size:', win.size), ylab='cumulative CO count'
       xlab="", ylab="")
  par(new=TRUE)
  plot(subset.tmp$START, subset.tmp$pw_LD, type='l', lwd=2, col='blue',
       xaxt="n", yaxt="n",
       xlab="", ylab="")
  axis(side=4)
  # mtext(text=paste(chrom.name, 'window start position - size(', win.size, 'bps)'), side=1, line=3, cex=0.75)
  mtext(text=chrom.name, side=1, line=3, cex=0.75)
  mtext(text=paste('cumulative CO count'), side=2, line=2, cex=0.75)
  mtext(text=paste('average pairwise LD'), side=4, line=2, cex=0.75)
  
}
par(op)


LD_v_CO_df <- cbind(LD_v_CO_df, win.idx=1:nrow(LD_v_CO_df))

op <- par(mfrow=c(1,1))
plot(LD_v_CO_df$win.idx, LD_v_CO_df$cum_CO, type='h', col='white', 
     ylim=c(-200,200),
     xlab='genomic window idx (1Mb)', ylab='pairwise LD/ cumulative CO count')
for (chrom.name in chrom.name.vec) {
  subset.tmp <- LD_v_CO_df[LD_v_CO_df$CHROM == chrom.name,]
  lines(subset.tmp$win.idx, subset.tmp$cum_CO, type='h', col='orange', lwd=2)
}
# par(new=TRUE)
lines(LD_v_CO_df$win.idx, -200*LD_v_CO_df$pw_LD, type='h', col='white')
for (chrom.name in chrom.name.vec) {
  subset.tmp <- LD_v_CO_df[LD_v_CO_df$CHROM == chrom.name,]
  lines(subset.tmp$win.idx, -200*subset.tmp$pw_LD, type='h', col='green', lwd=2)
}
par(op)

plot(LD_v_CO_df$pw_LD, LD_v_CO_df$cum_CO, xlab='Average pair-wise LD within 960 kb genomic windows', 
     ylab='cumulative CO count within 960 kb genomic windows')
lm_mod <- lm(cum_CO ~ pw_LD, data=LD_v_CO_df)
sum_lm_mod <- summary(lm_mod)
abline(lm_mod, lty=2, col='blue')
text(x=0.8, y=200, labels=paste0('r^2=', round(stat_cor, digits=3)))
# Some of the points with low LD and low CO count need to be investigated. 
LD_v_CO_df[LD_v_CO_df$pw_LD < 0.3 & LD_v_CO_df$cum_CO < 25,]  # Identify points with low LD and low CO counts
# These appear to be starting and ending windows that are smaler than the regular genomic windows thus carry
# low Co counts.




# Remove terminal windows and re-plot the relationship between LD and cumulative CO counts
trimmed_df <- NULL
for (chrom.name in chrom.name.vec) {
  subset.tmp <- LD_v_CO_df[LD_v_CO_df$CHROM == chrom.name,]
  subset.tmp <- subset.tmp[-c(1, nrow(subset.tmp)), ]
  trimmed_df <- rbind(trimmed_df, subset.tmp)
}

pdf(file="./statistics/gross_stat_plots/COs_vs_LD_960kb.pdf")
stat_cor<-cor.test(trimmed_df$pw_LD, trimmed_df$cum_CO)
plot(trimmed_df$pw_LD, trimmed_df$cum_CO,
     pch=16,
     col='grey',
     xlab='Average pair-wise LD within 960 kb genomic windows', 
     ylab='cumulative CO count within 960 kb genomic windows')
abline(lm_mod, lty=2, col='blue')
text(x=0.8, y=200, labels=paste0('r=', round(stat_cor$estimate, digits=3)))
dev.off()

trimmed_df[trimmed_df$pw_LD < 0.3 & trimmed_df$cum_CO < 40,]
summary(lm(trimmed_df$pw_LD~trimmed_df$cum_CO))
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################