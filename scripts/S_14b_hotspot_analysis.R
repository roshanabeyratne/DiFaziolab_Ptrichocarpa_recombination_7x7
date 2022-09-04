#!/usr/bin/R
########################################################################################################################################################################################################
# PURPOSE: Analyze HOT / COLD regions of the genome for recombination.
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: January 15, 2021
# USAGE: Input file include a detailed list of all the cross-over events observed for each focal parent for each offspring for all the chromosomes. The input is produced with script 22b_... 

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

########################################################################################################################################################################################################

# set working directory, load the packages and source functions used in the script
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
# win.size <- 960000 # (scale = 5)
# win.size <- 480000 # (scale = 4)
# win.size <- 240000 # (scale = 3)
# win.size <- 120000 # (scale = 2)
# win.size <- 60000 # (scale = 1)
win.size <- 30000
# win.size <- 20000 
# win.size <- 10000 
total.bin.span.coord.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr01"
  bin.span.coord.df.temp <- CreateBinSpanCoordinates(chrom.name = chrom.name, bin.size = win.size, chrom.length = chrom.size.df[chrom.size.df$CHROM == chrom.name, "LENGTH"], shift.size = win.size)
  bin.span.coord.df.temp2 <- cbind(rep(chrom.name, nrow(bin.span.coord.df.temp)), bin.span.coord.df.temp)
  total.bin.span.coord.df <- rbind(total.bin.span.coord.df, bin.span.coord.df.temp2)
}
total.bin.span.coord.df <- as.data.frame(total.bin.span.coord.df, stringsAsFactors = FALSE)
total.bin.span.coord.df[,2:3] <- apply(total.bin.span.coord.df[,2:3], 2, as.numeric)
colnames(total.bin.span.coord.df) <- c("CHROM", "START", "END")


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

# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_10kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_30kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_250kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# # all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_1Mb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )

all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_30kb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
colnames(all.par.df.2) <- colnames(all.par.df.2) %>% str_replace("^X", "") # or # colnames(all.par.df.2)[4:17] <- par.name.vec
all.par.df.2[(all.par.df.2$CHROM == 'Chr01' & all.par.df.2$START > 29054288), '6909'] <- NA
sum(apply(all.par.df.2[-(1:3)], 1, sum, na.rm=TRUE), na.rm=TRUE)
hist(apply(all.par.df.2[-(1:3)], 1, sum), breaks=100, xlab='cumulative COs per window', main='Histogram: CO count')



# CO distribution landscape per chromosome
op <- par(mfrow=c(4,5))
for (chrom.name in chrom.name.vec) {
  subset.tmp <- all.par.df.2[all.par.df.2$CHROM==chrom.name, ]
  plot(subset.tmp$START, apply(subset.tmp[-(1:3)], 1, sum, na.rm=TRUE), type='l', xlab=paste(chrom.name, 'win.size:', win.size), ylab='cumulative CO count')
}
par(op)



# Construct dataframes with number of cross-overs and number of individuals used for each parent for subsequent hotspot analysis 
num.COs.df <-  NULL
num.ind.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  tmp.num <- NULL
  tmp.num.2 <- NULL
  for (par.name in par.name.vec) { # par.name <- par.name.vec[10]
    
    tmp.tab <- CO.df[CO.df$CHROM == chrom.name & CO.df$PAR == par.name, ]
    tmp.num <- c(tmp.num, length(unique(tmp.tab$IND)))
    tmp.num.2 <- c(tmp.num.2, sum(tmp.tab$CO_SIZE != 0))
  }
  num.ind.df <- rbind(num.ind.df, tmp.num)
  num.COs.df <- rbind(num.COs.df, tmp.num.2)
}

num.ind.df <- as.data.frame(num.ind.df, stringsAsFactors = FALSE)
num.COs.df <- as.data.frame(num.COs.df, stringsAsFactors = FALSE)
colnames(num.ind.df) <- par.name.vec
colnames(num.COs.df) <- par.name.vec


# check whether the total number of COs are equal in both dataframes. Should deiifer if the windows are overlapping in the all.par.df.2 dataframe
sum(apply(all.par.df.2[,-(1:3)], 1, sum, na.rm=TRUE)) == sum(apply(num.COs.df, 1, sum)) 
fem.cum.COs <- sum(apply(all.par.df.2[,4:10], 1, sum, na.rm=TRUE)) 
male.cum.COs <- sum(apply(all.par.df.2[,11:17], 1, sum, na.rm=TRUE)) 
poisson.test(x=c(round(fem.cum.COs/7) , round(male.cum.COs/7)), alternative="two.sided")




# parameterizing the NULL model that the probability of observing a certain number of cross-overs within each window can be modeled 
# with a Poisson distribution for the whole genome at once
tot.num.cos <- sum(apply(num.COs.df, 2, sum, na.rm=TRUE))
cos.per.window <- apply(all.par.df.2[,-(1:3)], 1, sum, na.rm=TRUE) # hist(cos.per.window) # sum(cos.per.window)
num.windows <- nrow(all.par.df.2)

lambda <- tot.num.cos / num.windows # genomewide average 
iod <- var(cos.per.window) / lambda # index of dispersion

# This is a code block to check whether accounting for NAs applied to 6909 parent makes a difference to the genome-wide average lambda
# mean.cos.per.window <- apply(all.par.df.2[,-(1:3)], 1, mean, na.rm=TRUE) # hist(cos.per.window) # sum(cos.per.window)
# gw.cos.avg.per.win <- mean(mean.cos.per.window*14)
# It does not make a difference!



sum(apply(all.par.df.2[,-(1:3)], 1, mean, na.rm=TRUE) > (5*lambda)) # definition of hotspot according the maize paper. However, in our analysis this dows not yield anuy 30kb window
obs.p.vals <- ppois(cos.per.window, lambda = lambda, lower.tail = FALSE) # p-values for obs. based on Poisson-NULL parameterized by lambda above
bon.f.cut.off <- (0.05 / num.windows) # FDR methods cannot be used here as the p-values under the NULL are not uniformly distributed for count data with these small values for counts.





# Visualize the distribution of COs over windows in comparison with a POisson NULL distribution parameterized with the genomewide average CO counts per window
set.seed(1234)
random.poiss.vec <- rpois(length(cos.per.window), lambda = lambda)
# jpeg(file = paste0(statistics.output.path, "observed_vs_Poisson_null_250kb.jpg"), width = 16, height = 9, unit = "in", res = 600)
# op <- par(mfrow=c(2,1))
pdf(file = paste0(statistics.output.path, "observed_vs_Poisson_null_250kb.pdf"), width = 8, height = 8)
hist(cos.per.window,  col = rgb(173, 216, 230, max = 255, alpha = 80, names = "lt.blue"), 
     breaks = seq(0, max(cos.per.window), 1),
     border=F,
     # probability=T, 
     xlim = c(0, 32), ylim = c(0, 6000),
     main = "Histogram: observed", xlab = "CO-count per window cumulative")
hist(random.poiss.vec,col = rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink"), 
     breaks = seq(0,max(cos.per.window),1),
     border=F,
     # probability=T, 
     add=T,
     xlim = c(0, 32), ylim = c(0, 6000),
     main = paste0("Histogram: Poisson_null", " (LAMBDA= ", round(lambda, digits=3), ")"), xlab = "CO-count per window cumulative")
dev.off()
par(op)
# qqplot(random.poiss.vec, cos.per.window)
# abline(0,1)
# dev.off()

# jpeg(file = paste0(statistics.output.path, "observed_p_value_on_Poisson_NULL.jpg"), width = 16, height = 9, unit = "in", res = 600)
# p-value distribution under NULL for discrete distributions do not show a Uniform distribution. This applies only to continuous distributions.
# hist(obs.p.vals, breaks = 100, main = paste0("Histogram of observed p-values under the Poisson NULL (lambda = ", lambda, ")"))
# dev.off()



# # calculating the p-values based on permutations of COs observed within windows
# tot.num.cos <- sum(apply(num.COs.df, 2, sum))
# num.windows <- nrow(total.bin.span.coord.df)
# lambda <- (tot.num.cos / num.windows)
# p.vals.empirical <- NULL
# for (i in 1:100) {
#   permuted.num.COs.df <- apply(all.par.df.2[,4:17], 2, function(x){
#     y <- sample(x, length(x), replace = FALSE)
#     return(y)
#   })
#   
#   cos.per.window <- apply(permuted.num.COs.df, 1, sum)
#   tmp.vals <- ppois(cos.per.window, lambda = lambda, lower.tail = FALSE)
#   p.vals.empirical <- c(p.vals.empirical, tmp.vals)
#   cat(i, " ")
# }
# # note to self: The p-value distribution for the NULL hypothesis is not UNIFORM for Poiss distributions to start with
# jpeg(file = paste0(statistics.output.path, "empirical_p_value_based_on_permuting_COs_observed_in_windows.jpg"), width = 16, height = 9, unit = "in", res = 600)
# hist(p.vals.empirical, breaks = 100, main = paste0("Histogram of emperical p-values under the Poisson NULL (lambda = ", lambda, ")"))
# dev.off()
# emp.p.cutt.off <- quantile(p.vals.empirical, c(0.05))




# plotting the p-value distribution for the whole genome with Bonferonni and empirical cut-offs. No false-positive screening is carried out yet. Will be done next
jpeg(file = paste0(statistics.output.path, "Hotspot_distribution_whole_genome_pre_filtering.jpg"), width = 1000, height = 500, quality = 5000)
plot(1:length(obs.p.vals), -1*(log10(obs.p.vals)), type="p", cex=0.5,
     xlab="non-overlap window index (size 30kbps)", 
     ylab="-log10(p-value)", 
     main=paste("whole-genome p-values (Poisson null):", expression(LAMBDA), "=", round(lambda, digits=1))
)
abline(h = -1*log10(bon.f.cut.off), lty=2, lwd=2, col='red')
# abline(h = -1 * log10(emp.p.cutt.off), lty = 2, lwd = 2, col = 'grey')
# legend(1, 30, legend = c("Bonferroni", "Empirical"), col = c("red", "grey"), lty = 2, lwd = 2, cex = 0.8)
# legend(1, 15, legend = c("Bonferroni cutt-off", "Empirical cut-off"), col = c("red", "grey"), lty = 2, lwd = 2, cex = 0.8)
legend(1, 15, legend = c("Bonferroni cutt-off"), col = c("red"), lty = 2, lwd = 2, cex = 0.8)
dev.off()




# We identified a re-arrangement artifact in Chr17 in Stet-14 we have removed these bins from further analysis.
total.bin.span.coord.df[which.max(-1*(log10(obs.p.vals))),]
tmp.CO.df <- CO.df[CO.df$CHROM=='Chr17',]
tmp.CO.df.2 <- tmp.CO.df[order(tmp.CO.df$DWNSTRM_FLNK, decreasing=FALSE), ]
idx.to.rmv <- which(-1*(log10(obs.p.vals)) > 15)
total.bin.span.coord.df[idx.to.rmv,]




# manually removing two windows in Chr17 identified as artifacts due to a re-arrangement in Chr17 of Stet-14 from all relevant data-sets 
idx.to.rmv.2 <- c(11708, 11749)
total.bin.span.coord.df <- total.bin.span.coord.df[-idx.to.rmv.2,]
all.par.df.3 <- all.par.df.2[-idx.to.rmv.2,]

info.mrkr.count.df <- read.csv(file=paste0(statistics.input.path, "info_mrkr_count_bin_size_30kb.csv"), stringsAsFactors=FALSE )
info.mrkr.count.df <- info.mrkr.count.df[-idx.to.rmv.2, ]

obs.p.vals <- obs.p.vals[-idx.to.rmv.2]
if (length(obs.p.vals) == nrow(all.par.df.3)) {
  hotspot.df <- cbind(all.par.df.3, p.value=obs.p.vals)
} else {
  stop("Dimension mismatch!")
}
hist(hotspot.df$p.value)




# identify how many windows surpass the Bonferroni cut.off and output the location. This will be useful for extracting sequences for DNA sequence motif analysis. 
idx.hotspots.bon <- which(hotspot.df$p.value < bon.f.cut.off)
length(idx.hotspots.bon)
# idx.hotspots.emp <- which(obs.p.vals < emp.p.cutt.off)
# lengp)th(idx.hotspots.em



# arbitrarily decided to select window indices for which the -log10(p-value) within the range 1.5 - 3.5. This is for the MEME based DNA seq motif analysis as negative controls.
hist(-log10(hotspot.df$p.value))
nc.tmp.idx <- which(-log10(hotspot.df$p.value) >= 1 & -log10(hotspot.df$p.value)  <= 5)
length(nc.tmp.idx)
write.csv(hotspot.df[nc.tmp.idx, ], file = paste0(statistics.output.path, "30kb_negative_control.csv"), quote = FALSE, row.names = FALSE)

# for (iter in 1:100) {
#   idx.neg.ctrl <- sample(idx.tmp, length(idx.hotspots.bon), replace=FALSE)
#   write.csv(all.par.df.3[idx.neg.ctrl, ], file = paste0(statistics.output.path, "30kb_negative_control_", iter, ".csv"), quote = FALSE, row.names = FALSE)
#   cat(iter, " ")
# }




# Below I am imposing two more stringent cut-off criteria in order to increase the confidence of the candidate CO hotspots.
# Following section is to filter out false-positives based on validating the higher avg.CO count observed for windows based on available informative markers OR  based on a few parents influencing
# the outcome.
# Investigate the trend between informative marker count and the number of observed CO counts
CO.count <- apply(all.par.df.3[,par.name.vec], 1, mean, na.rm=TRUE)
MRKR.count <-  apply(info.mrkr.count.df[-(1:3)], 1, mean, na.rm=TRUE)
if (length(CO.count) == length(MRKR.count)) {
  mrkr.vs.co.df <- cbind(info.mrkr.count.df[,1:3], AVG_MRKR_COUNT=MRKR.count,AVG_CO_COUNT=CO.count)
} else {
  stop("error in CO.count or info.mrkr.count.df")
}
# plot(CO.count, MRKR.count)
head(mrkr.vs.co.df)




# Obtain the avg.CO count for ech window and avg.informative marker count within the span of idx window and flanking windows. 
flnk.size <- 2
detailed.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
  subset.mrkr.vs.co.df <- mrkr.vs.co.df[mrkr.vs.co.df$CHROM==chrom.name,]
  feasible.idx <- c((1+flnk.size):(nrow(subset.mrkr.vs.co.df)-flnk.size))
  for (idx in feasible.idx) { # idx <- 1000
    co.count <- subset.mrkr.vs.co.df$AVG_CO_COUNT[idx]
    mrkr.idxs <- c((idx-flnk.size):(idx-1), idx, (idx+1):(idx+flnk.size))
    mrkr.count <- mean(subset.mrkr.vs.co.df$AVG_MRKR_COUNT[mrkr.idxs])
    return.vec <- c(subset.mrkr.vs.co.df$CHROM[idx], subset.mrkr.vs.co.df$START[idx], subset.mrkr.vs.co.df$END[idx], mrkr.count, co.count)
    detailed.df <- rbind(detailed.df, return.vec)
  }
  cat(chrom.name, "\n")
}
detailed.df <- as.data.frame(detailed.df, stringsAsFactors=FALSE)
colnames(detailed.df) <- c(colnames(mrkr.vs.co.df)[1:3], "AVG_MRKR_COUNT", "AVG_CO_COUNT")
detailed.df$START <- as.numeric(detailed.df$START)
detailed.df$END <- as.numeric(detailed.df$END)
detailed.df$AVG_MRKR_COUNT <- as.numeric(detailed.df$AVG_MRKR_COUNT)
detailed.df$AVG_CO_COUNT <- as.numeric(detailed.df$AVG_CO_COUNT)
str(detailed.df)
hist(detailed.df$AVG_CO_COUNT, breaks=50)



# plotting the avg.mrkr.count vs. the avg.co.count to identify outliers.
hist(detailed.df$AVG_MRKR_COUNT)
plot(detailed.df$AVG_MRKR_COUNT, detailed.df$AVG_CO_COUNT, main='plot to detect CO outliers', xlab='mean number of flanking markers', ylab='CO count')

# # to contrast with the average number of COs per window. Taking the distribution of COs per window for all individuals
# hist(as.numeric(as.matrix(all.par.df.3[-(1:3)])), breaks=50)
# sum(as.numeric(as.matrix(all.par.df.3[-(1:3)]))) / length(as.numeric(as.matrix(all.par.df.3[-(1:3)])))
# for (j in 4:17) {
#   cat(colnames(all.par.df.3)[j], ": " ,which(all.par.df.3[,j] > 6), "\n")
# }
# 
# all.par.df.3[c(373),]
# all.par.df.3[c(799, 5018, 7156),]
# all.par.df.3[c(516, 3699),]
# all.par.df.3[c(2380, 5027, 10854),]
# all.par.df.3[c(568, 10434, 11649),]




# # check the number of outlier windows without flanking markers but >0 CO-events
# CO.counts.vs.mrkr.req <- NULL
# for (i in 0:200) {
#   tmp.var <- sum((detailed.df$AVG_MRKR_COUNT >=i & detailed.df$AVG_MRKR_COUNT < i+1) & (detailed.df$AVG_CO_COUNT > 0 & detailed.df$AVG_CO_COUNT < 0.5))
#   CO.counts.vs.mrkr.req <- rbind(CO.counts.vs.mrkr.req, c(i, tmp.var))
# }
# plot(CO.counts.vs.mrkr.req[,1], CO.counts.vs.mrkr.req[,2], xlab='number of flanking markers', ylab='windows with CO count > 0')




# Identify the outliers and overlap with identified hotspots.
length(idx.hotspots.bon)
hotspot.candidates <- hotspot.df[idx.hotspots.bon, ]
hist(as.numeric(as.matrix(hotspot.candidates[-(1:3)])), breaks=50)
sum(as.numeric(as.matrix(hotspot.candidates[-(1:3)]))) / length(as.numeric(as.matrix(hotspot.candidates[-(1:3)])))

# First-critera to reduce False-positives
# These are the outliers identified with few individuals skewing the result
poten.outliers <- NULL
thresh.co.count <- 4
for (j in 4:17) {
  cat(colnames(hotspot.candidates)[j], ": " ,which(hotspot.candidates[,j] > thresh.co.count), "\n")
  poten.outliers <- c(poten.outliers, which(hotspot.candidates[,j] > thresh.co.count))
}
# hotspot.candidates[sort(poten.outliers),]
# nrow(hotspot.candidates[sort(poten.outliers),])



# evaluate the index of dispersion for the whole dataset of windows as an alternative to the above arbitrary method.
# evaluate index of dispersion for all windows across all parents in order to gauge the null iod pattern 
iod.vec <- apply(all.par.df.3[-(1:3)], 1, function(x){
  var.i <- var(x, na.rm=TRUE)
  mean.i <- mean(x, na.rm=TRUE)
  iod <- ifelse(mean.i > 0, var.i/mean.i, NA)
  return(iod)
})
hist(iod.vec)
iod.cut.off <- quantile(iod.vec, c(0.1, 0.9), na.rm=T)

# same analysis as above for candidate hotspots
iod.vec.cand <- apply(hotspot.candidates[4:17], 1, function(x){
  var.i <- var(x, na.rm=TRUE)
  mean.i <- mean(x, na.rm=TRUE)
  iod <- ifelse(mean.i > 0, var.i/mean.i, NA)
  return(iod)
})
iod.outlier.idx <- which(iod.vec.cand < iod.cut.off[1] | iod.vec.cand > iod.cut.off[2])
length(iod.outlier.idx) # 50




# Second criteria for reducing False-positives.
# These are the outliers identified with flanking marker count
idx.mrkr.count.outliers <- which((detailed.df$AVG_CO_COUNT > 0) & (detailed.df$AVG_MRKR_COUNT < 2)) # windows that contain at least one CO in a window and where the average marker count in 2 adjacent flanking 
# windows is less than 1 
mrkr.count.outliers.df <- detailed.df[idx.mrkr.count.outliers, ]
tmp.idx <- which((hotspot.candidates$CHROM %in% mrkr.count.outliers.df$CHROM) & (hotspot.candidates$START %in% mrkr.count.outliers.df$START))
hotspot.candidates[tmp.idx, ]



# merging the two-sets of outliers
# tmp.hotspot.rmv <- rbind(hotspot.candidates[sort(poten.outliers),], hotspot.candidates[tmp.idx, ]) # iod.outlier.idx
tmp.hotspot.rmv <- rbind(hotspot.candidates[sort(iod.outlier.idx),], hotspot.candidates[tmp.idx, ]) 
tmp.hotspot.rmv.df <- unique(tmp.hotspot.rmv)
tmp.hotspot.rmv.df
dim(tmp.hotspot.rmv.df)

confint.hotspot.candidates <- hotspot.candidates[!((hotspot.candidates$CHROM %in% tmp.hotspot.rmv.df$CHROM) & (hotspot.candidates$START %in% tmp.hotspot.rmv.df$START)),]
nrow(confint.hotspot.candidates)


# Identify the actual hotspots from the set
conf.hotspot.cand.idx <- which((hotspot.df$CHROM %in% confint.hotspot.candidates$CHROM) & (hotspot.df$START %in% confint.hotspot.candidates$START))
length(conf.hotspot.cand.idx) # 67


# hotspot distribution with false-positives marked in a different color
pdf(file = paste0(statistics.output.path, "Hotspot_distribution_whole_genome_post_filtering.pdf"), width = 15, height = 10)
# jpeg(file = paste0(statistics.output.path, "Hotspot_distribution_whole_genome_post_filtering.jpg"), width = 1000, height = 500, quality = 5000)
# op <- par(mfrow=c(2,1), mar=c(4.1, 4.1, 1, 2.1))
plot(1:length(hotspot.df$p.value), -1*(log10(hotspot.df$p.value)), type="p", cex=0.5, col='white',
     xaxt='n',
     # xlab="non-overlap window index (size 30kbps)",
     xlab='',
     ylab="-log10(p-value)", 
     main=paste("whole-genome p-values (Poisson null):", expression(LAMBDA), "=", round(lambda, digits=1))
)
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[2]
  col.pal <- ifelse(which(chrom.name.vec %in% chrom.name) %% 2 == 0, 'black', 'grey')
  win.idxs.for.chrom <- which(hotspot.df$CHROM == chrom.name)
  points(win.idxs.for.chrom, -1*(log10(hotspot.df$p.value))[win.idxs.for.chrom], col=col.pal, pch=19, cex=0.5)
}
points(conf.hotspot.cand.idx, -1*(log10(hotspot.df$p.value))[conf.hotspot.cand.idx], col='blue', cex=1, pch=23)
abline(h = -1*log10(bon.f.cut.off), lty=2, lwd=2, col='red')
legend(1, 15, legend = c("Bonferroni cutt-off"), col = c("red"), lty = 2, lwd = 2, cex = 0.8, box.lwd=0)
# Identify mid-point win.idx for each chromosome
axis.points <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[2]
  win.idxs.for.chrom <- which(hotspot.df$CHROM == chrom.name)
  mid.win.idxs.for.chrom <- round(mean(c(win.idxs.for.chrom[1], win.idxs.for.chrom[length(win.idxs.for.chrom)])))
  axis.points <- c(axis.points, mid.win.idxs.for.chrom)
}
# Draw new x-axis
axis(side=1, at=axis.points, labels=chrom.name.vec, las=2 )
# if (nrow(hotspot.df) == nrow(info.mrkr.count.df)) {
#   mean.mrkr.count <- apply(info.mrkr.count.df[-(1:3)], 1, mean)
#   plot(mean.mrkr.count, type='h', xlab="non-overlap window index (size 30kbps)", ylab='mean informative marker count per window', col='grey')
#   abline(h = mean(mean.mrkr.count), lty=2, lwd=2, col='white')
# } else {
#   stop("ERROR in dimesions!")
# }
par(op)
dev.off()


# # Same data as above visualized for one chromosome at a time 
# jpeg(file = paste0(statistics.output.path, "Hotspots_per_chromosome.jpg"), width = 16, height = 9, unit = "in", res = 600)
# op <- par(mfrow = c(2,10))
# for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr01"
# 
#   test.df <- hotspot.df[hotspot.df$CHROM == chrom.name, -(1:3)]
#   num.windows <- nrow(test.df)
# 
#   tot.num.cos <- sum(num.COs.df[which(chrom.name.vec == chrom.name), ])
#   lambda <- (tot.num.cos / num.windows)
# 
#   cos.per.window <- apply(test.df, 1, sum)
#   obs.p.vals <- ppois(cos.per.window, lambda = lambda, lower.tail = FALSE)
#   # hist(obs.p.vals, breaks = 100)
# 
#   plot(1:num.windows, -1 * (log10(test.df$p.value)), type = "p", pch = 19, cex = 0.25, xlab = "window.idx", ylab = "-log10(p-value)", main = chrom.name)
#   bon.f.cut.off <- (0.05 / num.windows)
#   abline(h = -1 * log10(bon.f.cut.off))
# }
# par(op)
# dev.off()


nrow(hotspot.df[conf.hotspot.cand.idx, ])
write.csv(hotspot.df[conf.hotspot.cand.idx, ], file = paste0(statistics.output.path, "30kb_hotspots_bonferronni_cutoff_post_filtering.csv"), quote = FALSE, row.names = FALSE)



# # Doing the same analysis as above for males and females separately
# pdf(file = paste0(statistics.output.path, "Hotspot_distribution_dam_v_sires.pdf"), width = 16, height = 9, onefile = TRUE)
all.par.df.4 <- all.par.df.3[!(all.par.df.3$CHROM %in% tmp.hotspot.rmv.df$CHROM & all.par.df.3$START %in% tmp.hotspot.rmv.df$START & all.par.df.3$END %in% tmp.hotspot.rmv.df$END), ]
cos.per.window.fem <- apply(all.par.df.4[,4:10], 1, sum, na.rm=TRUE)
cos.per.window.male <- apply(all.par.df.4[,11:17], 1, sum, na.rm=TRUE)

tot.num.cos.fem <- sum(cos.per.window.fem) 
tot.num.cos.male <- sum(cos.per.window.male) 

lambda.fem <- (tot.num.cos.fem /length(cos.per.window.fem))
lambda.male <- (tot.num.cos.male /length(cos.per.window.male))

obs.p.vals.fem <- ppois(cos.per.window.fem, lambda = lambda.fem, lower.tail = FALSE)
obs.p.vals.male <- ppois(cos.per.window.male, lambda = lambda.male, lower.tail = FALSE)

y.max <- max(c(-log10(obs.p.vals.fem), -log10(obs.p.vals.male)))

idx.hotspots.bon.fem <- which(obs.p.vals.fem < bon.f.cut.off)
idx.hotspots.bon.male <- which(obs.p.vals.male < bon.f.cut.off)

length(idx.hotspots.bon.fem) # 38
length(idx.hotspots.bon.male) # 24

sum(idx.hotspots.bon.fem %in% idx.hotspots.bon.male) # 0

sum(all.par.df.4$CHROM[idx.hotspots.bon.fem] %in% confint.hotspot.candidates$CHROM & all.par.df.4$START[idx.hotspots.bon.fem] %in% confint.hotspot.candidates$START) # 21
sum(all.par.df.4$CHROM[idx.hotspots.bon.male] %in% confint.hotspot.candidates$CHROM & all.par.df.4$START[idx.hotspots.bon.male] %in% confint.hotspot.candidates$START) # 11



# # Testing for correlation between male and female within hotspots
# plot(cos.per.window.fem[idx.hotspots.bon], cos.per.window.male[idx.hotspots.bon])
# cor.hot <- cor.test(cos.per.window.fem[idx.hotspots.bon], cos.per.window.male[idx.hotspots.bon], method='pearson')
# cor.test(cos.per.window.fem, cos.per.window.male,method='pearson')
# # since hotpots are more negatively correlated than when the whole genome is condsidered, I randomly sampled 127 locations with replacement to to obtain a confidence interval
# permut.corr <- numeric(length=1000)
# for (iter in 1:1000) {
#   tmp.idx <- sample(1:length(cos.per.window.fem), length(idx.hotspots.bon), replace=FALSE)
#   tmp.cor <- cor.test(cos.per.window.fem[tmp.idx], cos.per.window.male[tmp.idx], method='pearson')
#   permut.corr[iter] <- tmp.cor$estimate
# }
# range(permut.corr)
# quantile(permut.corr, c(0.025, 0.975))

# # Creating a cloud of points for all bins males vs. female counts 
# fem.male.tab <- as.data.frame(table(cos.per.window.fem, cos.per.window.male))
# colnames(fem.male.tab) <- c("FEM", "MALE", "value")
# jpeg(file = paste0(statistics.output.path, "fem_vs_male_CO_count_heatmap.jpg"), width = 1000, height = 1000, quality = 5000)
# # op <- par(c(1,2))
# ggplot(fem.male.tab, aes(FEM, MALE, col=value, fill=value, label=value)) +
#   geom_tile() +
#   geom_text(col="black") +
#   theme_minimal() +
#   scale_fill_gradient2(low="white", mid="yellow", high="red") +
#   scale_color_gradient2(low="white", mid="yellow", high="red")
# # par(op)
# dev.off()
# # produce same plot above for SEX permuted data.



# jpeg(file = paste0(statistics.output.path, "Hotspot_distribution_dam_v_sires_whole_genome.jpg"), width = 1000, height = 500, quality = 150)

pdf(file = paste0(statistics.output.path, "Hotspot_distribution_dam_v_sires_whole_genome.pdf"), width = 8, height = 8)
op<-par(cex=1)
ymin <- min(c(-1*(log10(obs.p.vals.fem)), 1*(log10(obs.p.vals.male))))
ymax <- max(c(-1*(log10(obs.p.vals.fem)), 1*(log10(obs.p.vals.male))))

plot(1:length(obs.p.vals.fem), -1*(log10(obs.p.vals.fem)), type="p", pch=16, cex=0.7, col='red',
     ylim=c(ymin, ymax),
     yaxt='n',
     # main=paste0("whole-genome p-values(poisson NULL: female-lambda=", round(lambda.fem, digits=1), ";", "male-lambda=", round(lambda.male, digits=1),")"),
     xlab="non-overlap window index (size 30kbps)", ylab="-log10(p-value)" 
)
axis(side=2, at=c(-10,-5,0,5,10), labels=c(10,5,0,5,10))

abline(h = -1*log10(bon.f.cut.off), lty=2, lwd=2, col='black')
points(1:length(obs.p.vals.male), log10(obs.p.vals.male), type="p", pch=16, cex=0.7, col='blue')
abline(h = log10(bon.f.cut.off), lty=2, lwd=2, col='black')
points(idx.hotspots.bon.fem, -1*log10(obs.p.vals.fem[idx.hotspots.bon.fem]), type="p", cex=0.8, col='red', pch=23)
points(idx.hotspots.bon.male, 1*log10(obs.p.vals.male[idx.hotspots.bon.male]), type="p", cex=0.8, col='blue', pch=23)
legend(1, -9, legend = c("Bonferroni cutt-off"), col = c('black'), lty = 2, lwd = 2, cex = 0.8, box.col='grey', box.lwd=0)
par(op)
dev.off()






###################################################################################
# Below I am analyzing the relationship of various genomic features with CO hotspots vs. randomly selected neutral windows
# read-in the genomic features file
win.size <- '30kb'
copia <- read.table(file=paste0(input.path.genomic.features, win.size, "_window_copia_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
gypsy <- read.table(file=paste0(input.path.genomic.features, win.size, "_window_gypsy_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
simp.rep <- read.table(file=paste0(input.path.genomic.features, win.size, "_window_Simple_repeat_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
genes <- read.table(file=paste0(input.path.genomic.features, win.size, "_window_gene_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
nuc.content <- read.table(file=paste0(input.path.genomic.features, win.size, "_window_NUC_content_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
colnames(nuc.content) <- c("CHROM", "START" ,"END", "pct_at", "pct_gc", "num_A", "num_C", 'num_G', "num_T", "num_N", "num_oth", "seq_len")
tail(nuc.content)
if (nrow(all.par.df.2) == nrow(nuc.content)) {
  tot.co.count <- apply(all.par.df.2[-(1:3)], 1, sum)
  gen.feat.df <- cbind(nuc.content, COPIA=copia[,4], GYPSY=gypsy[,4], SIMP_REPEAT=simp.rep[,4], GENE=genes[,4], CUM_CO_COUNT=tot.co.count, 
                       log_CUM_CO_COUNT=log10(tot.co.count+0.05), WIN_IDX=1:nrow(nuc.content))
} else {
  stop("ERROR in dimensions!")
}
head(gen.feat.df)
gen.feat.df <- cbind(gen.feat.df, phys.Mb=round(gen.feat.df$START/(10^7)))



# observe the general trend of counts for each genomic feature.
op <- par(mfrow=c(2,3))
for (i in c("COPIA", "GYPSY", "SIMP_REPEAT", "GENE", "pct_gc", "CUM_CO_COUNT")) {
  hist(gen.feat.df[,which(colnames(gen.feat.df) %in% i)])
}
par(op)
# deciding to go with a count of <= 20 as low and anything greater than that as high



gen.feat.df.2 <- gen.feat.df[which(!is.na(gen.feat.df$CUM_CO_COUNT)), ]
head(gen.feat.df.2)

# Visualize the trend of repeat elements and other genomic features in comparison to CO-count using a CIRCOS-plot
library("circlize")
circos.clear()
circos.par( track.height = 0.1, track.margin=c(0,0), gap.degree=3)
circos.initialize(sectors=gen.feat.df.2$CHROM, x=gen.feat.df.2$START)
circos.track(sectors=gen.feat.df.2$CHROM, x=gen.feat.df.2$START, y=gen.feat.df.2$CUM_CO_COUNT,
             panel.fun = function(x, y) {
               circos.lines(x, y)
               circos.axis(h='top', major.tick=TRUE, direction='outside', labels.facing='outside', major.at=10^7*(seq(0, 20, 0.5)), minor.ticks=8, labels.cex = 0.25)
               circos.yaxis(side='left', labels.cex=0.25 )
             })

circos.track(sectors=gen.feat.df.2$CHROM, x=gen.feat.df.2$START, y=gen.feat.df.2$GENE,
             panel.fun = function(x, y) {
               circos.lines(x, y, col='blue')
               circos.yaxis(side='left', labels.cex=0.25 )
             })
circos.track(sectors=gen.feat.df.2$CHROM, x=gen.feat.df.2$START, y=gen.feat.df.2$COPIA,
             panel.fun = function(x, y) {
               circos.lines(x, y, col='purple')
               circos.yaxis(side='left', labels.cex=0.25 )
             })
circos.track(sectors=gen.feat.df.2$CHROM, x=gen.feat.df.2$START, y=gen.feat.df.2$GYPSY,
             panel.fun = function(x, y) {
               circos.lines(x, y, col='orange')
               circos.yaxis(side='left', labels.cex=0.25 )
             })
circos.track(sectors=gen.feat.df.2$CHROM, x=gen.feat.df.2$START, y=gen.feat.df.2$SIMP_REPEAT,
             panel.fun = function(x, y) {
               circos.lines(x, y, col='green')
               circos.yaxis(side='left', labels.cex=0.25 )
             })
circos.info()
circos.clear()





negative.control.df <- (hotspot.df[nc.tmp.idx,])
candidate.hotspot.df <- hotspot.df[conf.hotspot.cand.idx,]
nrow(negative.control.df)
nrow(candidate.hotspot.df)

nc.idx.gen.feat <- which(gen.feat.df$CHROM %in% negative.control.df$CHROM & gen.feat.df$START %in%  negative.control.df$START & gen.feat.df$END %in%  negative.control.df$END)
hot.idx.gen.feat <- which(gen.feat.df$CHROM %in% candidate.hotspot.df$CHROM & gen.feat.df$START %in%  candidate.hotspot.df$START & gen.feat.df$END %in%  candidate.hotspot.df$END)
length(nc.idx.gen.feat)
length(hot.idx.gen.feat)

nc.combined.df <- cbind(negative.control.df, gen.feat.df[nc.idx.gen.feat, ], CATEG=rep("BG", length(nc.idx.gen.feat)))
candidate.combined.df <- cbind(candidate.hotspot.df, gen.feat.df[hot.idx.gen.feat, ], CATEG=rep("HOT", length(hot.idx.gen.feat)))
temp.df <- rbind(nc.combined.df, candidate.combined.df)
combined.df <- temp.df[order(temp.df$WIN_IDX), ]

# # checking
# head(combined.df)
# sum(combined.df$CATEG == 'HOT')


# visualize the difference in genomic features using boxplots
pdf(file = paste0(statistics.output.path, "hotspot_gw_features.pdf"), width=8, height=16)
op <- par(mfrow=c(6, 1), mar=c(3.1,4.1,2.1,2.1), cex=1)

# boxplot(pct_at ~ CATEG, data=combined.df, xlab=NULL, ylab="% AT", col=c("lightpink", "firebrick"), border="black", las=2, horizontal=T)
boxplot(pct_gc ~ CATEG, data=combined.df, xlab=NULL, ylab="% GC", col=c("lightpink", "firebrick"), border="black", las=2, horizontal=T)
boxplot(CUM_CO_COUNT ~ CATEG, data=combined.df, xlab=NULL, ylab="Average CO count", col=c("lightpink", "firebrick"), border="black", las=2, horizontal=T, ylim=c(0,45))
boxplot((GENE + 1) ~ CATEG, data=combined.df, xlab=NULL, ylab="Gene count", col=c("lightpink", "firebrick"), border="black", las=2, horizontal=T, ylim=c(0,45))
boxplot((COPIA + 1) ~ CATEG, data=combined.df, xlab=NULL, ylab="Copia count", col=c("lightpink", "firebrick"), border="black", las=2, horizontal=T, ylim=c(0,45))
boxplot((GYPSY + 1) ~ CATEG, data=combined.df, xlab=NULL, ylab="Gypsy count", col=c("lightpink", "firebrick"), border="black", las=2, horizontal=T, ylim=c(0,45))
boxplot((SIMP_REPEAT + 1) ~ CATEG, data=combined.df, xlab=NULL, ylab="Simple repeats count", col=c("lightpink", "firebrick"), border="black", las=2, horizontal=T, ylim=c(0,45))
par(op)
dev.off()

t.test(pct_at ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(pct_gc ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(log(GENE + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(log(COPIA + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(log(GYPSY + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
t.test(log(SIMP_REPEAT + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)

# compare the means of the two categories using a non parametric test (Wilcoxon-Mann-Whitney test)
wilcox.test(formula=CUM_CO_COUNT ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=pct_at ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=pct_gc ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=(GENE + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=(COPIA + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=(GYPSY + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)
wilcox.test(formula=(SIMP_REPEAT + 1) ~ CATEG, data=combined.df, alternatie='two.sided', var.equal = FALSE)


# # testing for association with genomic features
# # https://online.stat.psu.edu/stat500/lesson/8/8.1
# 
# # COPIA
# cut.off.gf <- quantile(gen.feat.df$COPIA, c(0.5, 0.75, 0.8, 0.9, 0.95))[3]
# num.hi.copia.in.hs <- sum(gen.feat.df$COPIA[hot.idx.gen.feat]<= cut.off.gf)
# num.lo.copia.in.hs <- sum(gen.feat.df$COPIA[hot.idx.gen.feat] > cut.off.gf)
# num.hi.copia.in.cs <- sum(gen.feat.df$COPIA[nc.idx.gen.feat]<= cut.off.gf)
# num.lo.copia.in.cs <- sum(gen.feat.df$COPIA[nc.idx.gen.feat] > cut.off.gf)
# copia.chi <- chisq.test(matrix(c(num.hi.copia.in.hs, num.lo.copia.in.hs, num.hi.copia.in.cs, num.lo.copia.in.cs), ncol=2, byrow=TRUE))
# 
# # GYPSY
# cut.off.gf <- quantile(gen.feat.df$GYPSY, c(0.5, 0.75, 0.8, 0.9, 0.95))[3]
# num.hi.gypsy.in.hs <- sum(gen.feat.df$GYPSY[hot.idx.gen.feat]<= cut.off.gf)
# num.lo.gypsy.in.hs <- sum(gen.feat.df$GYPSY[hot.idx.gen.feat] > cut.off.gf)
# num.hi.gypsy.in.cs <- sum(gen.feat.df$GYPSY[nc.idx.gen.feat]<= cut.off.gf)
# num.lo.gypsy.in.cs <- sum(gen.feat.df$GYPSY[nc.idx.gen.feat] > cut.off.gf)
# gypsy.chi <- chisq.test(matrix(c(num.hi.gypsy.in.hs, num.lo.gypsy.in.hs, num.hi.gypsy.in.cs, num.lo.gypsy.in.cs), ncol=2, byrow=TRUE))
# 
# # Simple Repeats
# cut.off.gf <- quantile(gen.feat.df$SIMP_REPEAT, c(0.5, 0.75, 0.8, 0.9, 0.95))[3]
# num.hi.sr.in.hs <- sum(gen.feat.df$SIMP_REPEAT[hot.idx.gen.feat]<= cut.off.gf)
# num.lo.sr.in.hs <- sum(gen.feat.df$SIMP_REPEAT[hot.idx.gen.feat] > cut.off.gf)
# num.hi.sr.in.cs <- sum(gen.feat.df$SIMP_REPEAT[nc.idx.gen.feat]<= cut.off.gf)
# num.lo.sr.in.cs <- sum(gen.feat.df$SIMP_REPEAT[nc.idx.gen.feat] > cut.off.gf)
# sr.chi <- chisq.test(matrix(c(num.hi.sr.in.hs, num.lo.sr.in.hs, num.hi.sr.in.cs, num.lo.sr.in.cs), ncol=2, byrow=TRUE))
# 
# 
# # Gene content
# cut.off.gf <- quantile(gen.feat.df$GENE, c(0.5, 0.75, 0.8, 0.9, 0.95))[3]
# num.hi.gene.in.hs <- sum(gen.feat.df$GENE[hot.idx.gen.feat]<= cut.off.gf)
# num.lo.gene.in.hs <- sum(gen.feat.df$GENE[hot.idx.gen.feat] > cut.off.gf)
# num.hi.gene.in.cs <- sum(gen.feat.df$GENE[nc.idx.gen.feat]<= cut.off.gf)
# num.lo.gene.in.cs <- sum(gen.feat.df$GENE[nc.idx.gen.feat] > cut.off.gf)
# gene.chi <- chisq.test(matrix(c(num.hi.gene.in.hs, num.lo.gene.in.hs, num.hi.gene.in.cs, num.lo.gene.in.cs), ncol=2, byrow=TRUE))
# 
# 
# # GC content
# cut.off.gf <- quantile(gen.feat.df$pct_gc, c(0.5, 0.75, 0.8, 0.9, 0.95))[3]
# num.hi.gc.in.hs <- sum(gen.feat.df$pct_gc[hot.idx.gen.feat]<= cut.off.gf)
# num.lo.gc.in.hs <- sum(gen.feat.df$pct_gc[hot.idx.gen.feat] > cut.off.gf)
# num.hi.gc.in.cs <- sum(gen.feat.df$pct_gc[nc.idx.gen.feat]<= cut.off.gf)
# num.lo.gc.in.cs <- sum(gen.feat.df$pct_gc[nc.idx.gen.feat] > cut.off.gf)
# gc.chi <- chisq.test(matrix(c(num.hi.gc.in.hs, num.lo.gc.in.hs, num.hi.gc.in.cs, num.lo.gc.in.cs), ncol=2, byrow=TRUE))
# 
# 
# 
# # AT content
# cut.off.at <- quantile(gen.feat.df$pct_at, c(0.5, 0.75, 0.8, 0.9, 0.95))[3]
# num.hi.at.in.hs <- sum(gen.feat.df$pct_at[hot.idx.gen.feat]<= cut.off.gf)
# num.lo.at.in.hs <- sum(gen.feat.df$pct_at[hot.idx.gen.feat] > cut.off.gf)
# num.hi.at.in.cs <- sum(gen.feat.df$pct_at[nc.idx.gen.feat]<= cut.off.gf)
# num.lo.at.in.cs <- sum(gen.feat.df$pct_at[nc.idx.gen.feat] > cut.off.gf)
# at.chi <- chisq.test(matrix(c(num.hi.gc.in.hs, num.lo.gc.in.hs, num.hi.gc.in.cs, num.lo.gc.in.cs), ncol=2, byrow=TRUE))
# 
# 


# 
# hotspot.total.repeats <- sum(gen.feat.df$SIMP_REPEAT[hot.idx.gen.feat], gen.feat.df$GYPSY[hot.idx.gen.feat], gen.feat.df$COPIA[hot.idx.gen.feat])
# hotspot.feat <- c(sum(gen.feat.df$SIMP_REPEAT[hot.idx.gen.feat]), sum(gen.feat.df$GYPSY[hot.idx.gen.feat]), sum(gen.feat.df$COPIA[hot.idx.gen.feat])) / hotspot.total.repeats
# 
# coldspot.total.repeats <- sum(gen.feat.df$SIMP_REPEAT[nc.idx.gen.feat], gen.feat.df$GYPSY[nc.idx.gen.feat], gen.feat.df$COPIA[nc.idx.gen.feat])
# coldspot.feat <- c(sum(gen.feat.df$SIMP_REPEAT[nc.idx.gen.feat]), sum(gen.feat.df$GYPSY[nc.idx.gen.feat]), sum(gen.feat.df$COPIA[nc.idx.gen.feat])) / coldspot.total.repeats
# 
# two.way.tab <- t(matrix(c(hotspot.feat, coldspot.feat), nrow=2, byrow=TRUE ))
# two.way.tab <- round(two.way.tab*100)
# chisq.test(two.way.tab)
# 
# 
# p.val.df <- NULL
# for (iter in 1:10000) {
#   
#   # nc.idx <- sample(nc.idx.gen.feat, length(hot.idx.gen.feat), replace=TRUE)
#   nc.idx <- sample(hot.idx.gen.feat, length(hot.idx.gen.feat), replace=TRUE)
#   
#   pc.idx <- sample(hot.idx.gen.feat, length(hot.idx.gen.feat), replace=TRUE)
#   # pc.idx <- sample(nc.idx.gen.feat, length(hot.idx.gen.feat), replace=TRUE)
#   
#   # p.val.copia <- chisq.test(c(sum(gen.feat.df$COPIA[nc.idx]), sum(gen.feat.df$COPIA[pc.idx])))$p.value
#   p.val.copia <- poisson.test(c(sum(gen.feat.df$COPIA[nc.idx]), sum(gen.feat.df$COPIA[pc.idx])), alternative="two.sided")$p.value
#   
#   # p.val.gypsy <- chisq.test(c(sum(gen.feat.df$GYPSY[nc.idx]), sum(gen.feat.df$GYPSY[pc.idx])))$p.value
#   p.val.gypsy <- poisson.test(c(sum(gen.feat.df$GYPSY[nc.idx]), sum(gen.feat.df$GYPSY[pc.idx])), alternative="two.sided")$p.value
#   
#   # p.val.simp.rep <- chisq.test(c(sum(gen.feat.df$SIMP_REPEAT[nc.idx]), sum(gen.feat.df$SIMP_REPEAT[pc.idx])))$p.value
#   p.val.simp.rep <- poisson.test(c(sum(gen.feat.df$SIMP_REPEAT[nc.idx]), sum(gen.feat.df$SIMP_REPEAT[pc.idx])), alternative="two.sided")$p.value
#   
#   # p.val.gene <- chisq.test(c(sum(gen.feat.df$GENE[nc.idx]), sum(gen.feat.df$GENE[pc.idx])))$p.value
#   p.val.gene <- poisson.test(c(sum(gen.feat.df$GENE[nc.idx]), sum(gen.feat.df$GENE[pc.idx])), alternative="two.sided")$p.value
#   
#   p.val.gc <- chisq.test(c(sum(gen.feat.df$pct_gc[nc.idx]), sum(gen.feat.df$pct_gc[pc.idx])))$p.value
#   # p.val.gc <- poisson.test(c(sum(gen.feat.df$pct_gc[nc.idx]), sum(gen.feat.df$pct_gc[pc.idx])), alternative="two.sided")$p.value
#   
#   p.val.df <- rbind(p.val.df, c(p.val.copia, p.val.gypsy, p.val.simp.rep, p.val.gene, p.val.gc))
# }
# # null.distro.p.val <- p.val.df
# observed.distro.p.val <- p.val.df
# 
# 
# 
# 
# op <- par(mfrow=c(2,5))
# for (i in 1:5) {
#   hist(null.distro.p.val[,i])
# }
# for (i in 1:5) {
#   hist(observed.distro.p.val[,i])
# }
# par(op)
# # There seems to be a higher association with simple repeats. 
# 
# 
# 
# 
# 
# p.val.vec <- NULL
# for (iter in 1:100) {
#   repeat.cont.df <- matrix(nrow=3, ncol=2, byrow=TRUE,  
#                            c(c(sum(sample(gen.feat.df$COPIA[nc.idx.gen.feat], length(hot.idx.gen.feat), replace=TRUE)), sum(gen.feat.df$COPIA[hot.idx.gen.feat])),
#                              c(sum(sample(gen.feat.df$GYPSY[nc.idx.gen.feat], length(hot.idx.gen.feat), replace=TRUE)), sum(gen.feat.df$GYPSY[hot.idx.gen.feat])),
#                              c(sum(sample(gen.feat.df$SIMP_REPEAT[nc.idx.gen.feat], length(hot.idx.gen.feat), replace=TRUE)), sum(gen.feat.df$SIMP_REPEAT[hot.idx.gen.feat])))
#                            
#                            # c(c(sum(gen.feat.df$COPIA[nc.idx.gen.feat])/length(nc.idx.gen.feat), sum(gen.feat.df$COPIA[hot.idx.gen.feat])/length(hot.idx.gen.feat)), 
#                            #   c(sum(gen.feat.df$GYPSY[nc.idx.gen.feat])/length(nc.idx.gen.feat), sum(gen.feat.df$GYPSY[hot.idx.gen.feat])/length(hot.idx.gen.feat)), 
#                            #   c(sum(gen.feat.df$SIMP_REPEAT[nc.idx.gen.feat])/length(nc.idx.gen.feat), sum(gen.feat.df$SIMP_REPEAT[hot.idx.gen.feat])/length(hot.idx.gen.feat)))
#                            
#   )
#   row.names(repeat.cont.df) <- c("COPIA", "GYPSY", "SIMP_REPEAT")
#   colnames(repeat.cont.df) <- c("COLD", "HOT")
#   p.val.vec <- c(p.val.vec, chisq.test(repeat.cont.df)$p.value)
# }
# hist(p.val.vec, breaks=10)
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################