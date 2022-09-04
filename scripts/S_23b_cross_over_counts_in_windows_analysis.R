#!/usr/bin/R
# setwd("/group/difazio/populus/gatk-7x7-Stet14/")
setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping/statistics")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: October 19, 2020
# USAGE: Input file include a detailed list of all the cross-over events observed for each focal parent for each offspring for all the chromosomes. The input is produced with script 22b_... 

# Following is the input file format used 
# CHROM  PAR   IND CO_SIZE DWNSTRM_FLNK UPSTRM_FLNK
# Chr01 1863 24708  103140      4675047     4778187
# Chr01 1863 24708   30914     26411972    26442886
# Chr01 1863 24708   61845     34445658    34507503
# Chr01 1863 24708    5612     46505278    46510890
# Chr01 1863 24709  164382       726848      891230
# Chr01 1863 24709  153359      2513166     2666525

# Following is the input file format of the respective chromosome size.  
# CHROM   LENGTH GEN_DIST
# Chr01 50495398      350
# Chr02 25263042      210
# Chr03 25000000      200
# Chr04 24267058      175
# Chr05 25890711      200
# Chr06 27912132      250

# PURPOSE: 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

statistics.output.path <- "./recombination_statistics/"
########################################################################################################################################################################################################
# The dataframe containing the respective chromosome sizes. 
chrom.size.df <- read.csv(file = "./scripts_lists/chrom.size.df.csv", header = TRUE, stringsAsFactors = FALSE)
# chrom.size.df <- read.csv(file = "./chrom.size.df.csv", header = TRUE, stringsAsFactors = FALSE)


# Read-in the dataframe with detailed cross-over events produced as per script 22b_...
# CO.df <- read.csv(paste0(statistics.output.path, "co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
CO.df <- read.csv(paste0("./co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)


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
  
  return(start.end.sites)
}




# Assign window start and end points based on window size
win.size <- 250000 # This is the final window size statistical analysis done on.
total.bin.span.coord.df <- NULL
for (chrom.name in chrom.name.vec) {
  bin.span.coord.df.temp <- CreateBinSpanCoordinates(chrom.name = chrom.name, bin.size = win.size, chrom.length = chrom.size.df[chrom.size.df$CHROM == chrom.name, "LENGTH"], shift.size = win.size)
  bin.span.coord.df.temp2 <- cbind(rep(chrom.name, nrow(bin.span.coord.df.temp)), bin.span.coord.df.temp)
  total.bin.span.coord.df <- rbind(total.bin.span.coord.df, bin.span.coord.df.temp2)
}

total.bin.span.coord.df <- as.data.frame(total.bin.span.coord.df, stringsAsFactors = FALSE)
total.bin.span.coord.df[,2:3] <- apply(total.bin.span.coord.df[,2:3], 2, as.numeric)
colnames(total.bin.span.coord.df) <- c("CHROM", "START", "END")



# Assign the COs to windows based on the the CO-region overlap. Creating a dataframe with CO counts in each overlapping windows
all.par.df <- NULL
for (par.name in par.name.vec) { # par.name <- par.name.vec[9]
  cat(date(), "\n")
  all.chroms.vec <- NULL
  for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[18]
    
    tmp.co.tab <- CO.df[(CO.df$CHROM == chrom.name & CO.df$PAR == par.name), ]
    tmp.co.tab <- tmp.co.tab[tmp.co.tab$CO_SIZE != 0, ]
    tmp.coord.tab <- total.bin.span.coord.df[total.bin.span.coord.df$CHROM == chrom.name, ]
    
    # tmp.vec <- rep(0, nrow(total.bin.span.coord.df))
    
    sum.vec <- rep(0, nrow(tmp.coord.tab))
    for (i in 1:nrow(tmp.co.tab)) { # i <- 2
      bp.coords <- c(tmp.co.tab[i,"LFT_FLNK"], tmp.co.tab[i,"RT_FLNK"])
      tmp.vec <- apply(tmp.coord.tab[2:3], 1, function(x, y = bp.coords){
        test.1 <- ifelse(((x[1] > y[1] & x[1] < y[2]) | (x[2] >= y[1] & x[2] < y[2])), 1, 0)
        test.2 <- ifelse((y[1] > x[1] & y[2] < x[2]) & (y[2] > x[1] & y[2] < x[2]) , 1, 0)
        test <- test.1 + test.2
        return(test)
      })
      sum.vec <- sum.vec + tmp.vec
    }
    all.chroms.vec <- c(all.chroms.vec, sum.vec)
  }
  cat(date(), "\n")
  all.par.df <- cbind(all.par.df, all.chroms.vec)
  cat(par.name, " done!", "\n")
}
colnames(all.par.df) <- par.name.vec
all.par.df.2 <- cbind(total.bin.span.coord.df, all.par.df)
write.csv(all.par.df.2, file = paste0(statistics.output.path, "co.bin_size_10kb.csv"), quote = FALSE, row.names = FALSE)

# all.par.df.2 <- read.csv(paste0(statistics.output.path, "co.bin_size_10kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(paste0(statistics.output.path, "co.bin_size_250kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = "./co.bin_size_10kb.csv", header = TRUE, stringsAsFactors = FALSE )
all.par.df.2 <- read.csv(file = "./co.bin_size_250kb.csv", header = TRUE, stringsAsFactors = FALSE )
colnames(all.par.df.2) <- colnames(all.par.df.2) %>% str_replace("^X", "")
# plot(apply(all.par.df.2[all.par.df.2$CHROM == "Chr01",4:17], 1, mean), typ = "l")


# Sample code block to plot the cross-over density
# library(RColorBrewer)
# colors <- c(brewer.pal(7, "Set3"), brewer.pal(7, "Set2"))
# for (par.name in par.name.vec[2]) { # par.name <- par.name.vec[1]
#   if (par.name == par.name.vec[1]) {
#     plot((all.par.df.2[all.par.df.2$CHROM == chrom.name.vec[4], paste0("X",par.name)]), typ = "l", col = colors[par.name.vec == par.name])
#   } else {
#     points((all.par.df.2[all.par.df.2$CHROM == chrom.name.vec[4], paste0("X",par.name)]), typ = "l", col = colors[par.name.vec == par.name])
#   }
# }




# Construct a dataframe for hotspot analysis 
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


tmp.df <- (all.par.df.2[all.par.df.2$CHROM == chrom.name.vec[1], 4:17])
plot(apply(tmp.df, 1, mean), typ = "l", col = 'black')
plot(tmp.df[,1], typ = "l", col = 'black')


all.par.df.3 <- cbind(seq(1, nrow(all.par.df.2), by = 1), all.par.df.2)
colnames(all.par.df.3) <- c("WIN_IDX", colnames(all.par.df.2))
reshaped.co.count.df <- NULL
for (par.name in par.name.vec) { # par.name <- "2515"
  for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr01"
    tmp.tab <- all.par.df.3[all.par.df.3$CHROM == chrom.name, c("WIN_IDX", "CHROM", "START", "END", par.name)]
    sex.assignment <-  ifelse(par.name %in% par.name.vec[1:7], "F", "M")
    num.inds <- num.ind.df[which(chrom.name %in% chrom.name.vec),par.name]
    comb.tmp.tab <- cbind(tmp.tab, rep(par.name, nrow(tmp.tab)), rep(sex.assignment, nrow(tmp.tab)), rep(num.inds, nrow(tmp.tab)))
    colnames(comb.tmp.tab) <- c(colnames(tmp.tab)[1:4], "CO_COUNT", "PAR", "SEX", "NUM_IND")
    
    reshaped.co.count.df <- rbind(reshaped.co.count.df, comb.tmp.tab)
  }
}
reshaped.co.count.df$SEX <- as.factor(reshaped.co.count.df$SEX)
reshaped.co.count.df$CHROM <- as.factor(reshaped.co.count.df$CHROM)
reshaped.co.count.df$PAR <- as.factor(reshaped.co.count.df$PAR)
reshaped.co.count.df$WIN_IDX <- as.factor(reshaped.co.count.df$WIN_IDX)

# head(reshaped.co.count.df)



reshaped.co.count.df.tmp <- reshaped.co.count.df[reshaped.co.count.df$CHROM == "Chr01", ]

# Applying a zero-inflated model
library(pscl)
f1 <- formula(CO_COUNT ~ WIN_IDX)
Zip1 <- zeroinfl(f1, dist = "negbin", link = "logit", data = reshaped.co.count.df.tmp)
hist(reshaped.co.count.df.tmp$CO_COUNT)










glm.poiss <- stats::glm(CO_COUNT ~ WIN_IDX, data = reshaped.co.count.df.tmp, family = poisson(link = "log")) 
glm.poiss <- stats::glm(CO_COUNT ~ SEX + NUM_IND + WIN_IDX, data = reshaped.co.count.df.tmp, family = poisson(link = "log")) 
summary(glm.poiss)
hist(glm.poiss$coefficients)
plot(glm.poiss$coefficients)
plot(glm.poiss$coefficients[glm.poiss$coefficients > -1])




glm.nb <- MASS::glm.nb(CO_COUNT ~ SEX + NUM_IND + WIN_IDX, PAR, link = "log", data = reshaped.co.count.df.tmp)
glm.nb <- MASS::glm.nb(CO_COUNT ~ SEX, NUM_IND + WIN_IDX, PAR, link = "log", data = reshaped.co.count.df.tmp)



summary(glm.nb)

plot(glm.poiss$coefficients[-(1:2)])
hist(glm.poiss$coefficients[-(1:2)])

plot(glm.nb$coefficients[-(1:2)])
hist(glm.nb$coefficients[-(1:2)])
100*(glm.nb$null.deviance - glm.nb$deviance) / glm.nb$null.deviance

library(lme4)
glmer.nb.fit <- lme4::glmer.nb(CO_COUNT ~  NUM_IND + WIN_IDX + (1|PAR), data = reshaped.co.count.df.tmp)







########################################################################################################################################################################################################
# Testing for model fit per window. 
library(lme4)
require(car)
library(MASS)
library(lmtest)

# identify rows (windows) that have all zero entries or an excess of zero entries.
rows.to.exclude <- which(apply(all.par.df.2[,4:17], 1, mean) < 0.3)
all.par.df.3 <- all.par.df.2[-rows.to.exclude, ]
colnames(all.par.df.3) <- colnames(all.par.df.3) %>% str_replace("^X", "")


glm.poiss.fit.vec <- NULL
glm.nb.fit.vec <- NULL
mod.compar.vec <- NULL

for (i in 1:nrow(all.par.df.3)) { # i <- 31
  sex.assignment <- rep(c(1,2), each = 7) # assign sex as a binary categorical factor
  
  num.ind.vec <- NULL
  for (par.name in par.name.vec) { # par.name <- "2365"
    num.ind.vec <- c(num.ind.vec, length(unique(CO.df[CO.df$CHROM == all.par.df.3$CHROM[i] & CO.df$PAR == par.name, "IND"])))
    
  }
  # sum(num.ind.vec)
  
  # subset dataframe per window
  stat.tab <- as.data.frame(cbind(as.vector(as.matrix(all.par.df.3[i,4:17])), num.ind.vec, sex.assignment))
  colnames(stat.tab) <- c("CO_COUNTS", "NUM_IND", "SEX")
  tmp.var <- factor(stat.tab$SEX)
  stat.tab$SEX <- tmp.var
  
  glm.poiss <- stats::glm(CO_COUNTS ~ SEX + NUM_IND, data = stat.tab, family = poisson(link = "log")) # POISSON
  glm.nb <- MASS::glm.nb(CO_COUNTS ~ SEX + NUM_IND, data = stat.tab) # NEGATIVE BINOMIAL; corrects the poisson model for overdispersion and zero inflation
  sum.glm.poiss <- summary(glm.poiss)
  sum.glm.nb <- summary(glm.nb)
  
  glm.poiss.fit.vec <- c(glm.poiss.fit.vec, pchisq(sum.glm.poiss$deviance, sum.glm.poiss$df.residual, lower.tail = FALSE)) # http://datavoreconsulting.com/programming-tips/count-data-glms-choosing-poisson-negative-binomial-zero-inflated-poisson/
  glm.nb.fit.vec <- c(glm.nb.fit.vec, pchisq(sum.glm.nb$deviance, sum.glm.nb$df.residual, lower.tail = FALSE))
  
  lr.stat <- lrtest(glm.nb, glm.poiss)$Chisq[2] # model comparison using a log-likelihood ratio (LR) test
  # statistic <- as.numeric(2*((logLik(glm.nb)) - (logLik(glm.poiss)))) # lrtest is an implementation of the same calculations mentioned in above lines: https://stats.stackexchange.com/questions/127505/compare-poisson-and-negative-binomial-regression-with-lr-test
  # pchisq(statistic, 1, lower.tail = FALSE)
  mod.compar.vec <- c(mod.compar.vec, lr.stat)
}



# This is a plot that visualizes the (CO count variance / CO count mean) i.e. coefficeint of dispersion (VMR) for each window. 
# Following is the expected value of VMR for several of the canonical distributions.  
# constant random variable	VMR = 0	not dispersed
# binomial distribution	0 < VMR < 1	under-dispersed
# Poisson distribution	VMR = 1	
# negative binomial distribution	VMR > 1	over-dispersed

setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping/statistics")
# all.par.df.2 <- read.csv(file = "./co.bin_size_250kb.csv", header = TRUE, stringsAsFactors = FALSE )

jpeg(filename = paste0(statistics.output.path, "variance_inflation.jpg"), width = 480, height = 480, quality = 500)
win.mean.vec <- apply(all.par.df.2[,4:17], 1, mean) # mean estimate for each window
win.var.vec <- apply(all.par.df.2[,4:17], 1, var) # variance estimate for window
# plotting over-dispersion
plot(win.mean.vec, win.var.vec, xlab = 'mean', ylab = 'variance')
var.inf.function <- lm(win.var.vec~win.mean.vec) # shows existence of variance inflation
sum.var.inf.function <- summary(var.inf.function)
mtext(paste0("R^2 = ", round(sum.var.inf.function$adj.r.squared, 4)))
abline(var.inf.function, col = 'blue')
lines(0:25, 0:25, col = 'black')
dev.off()




########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################   