#!/usr/bin/R
setwd("/group/difazio/populus/gatk-7x7-Stet14/")
# setwd("/Users/cabeyrat/ownCloud/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA, for the recombination landscape analysis for 7x7 cross
# DATE: September 20, 2020
# USAGE: Input files include smoothed-out bakground haplotype contribution from the focal parent. Some individuals were removed from relevant parents due to high missingness and switches with script
# 19b. Therefore the path to such curated files is "./resolved_haps/by_parent/curated/". If a parent is not corrected such as described before, the focal parent hap contribution is at
# "./resolved_haps/by_parent/"

# Following is the input file format used that has focal parent haplotype identification (hap_1 or hap_2) for all of its offspring;
# CHROM,POS,24748,24749,24752,24753,24756,26643,26649,26650,26652,26653,26654,26655,26657,26660,26664,26668,26675,26677,25289,25297,25301,25305,25307,25310,25311,25317,25486,25487,25490,25491,25492
# Chr02,6617,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1
# Chr02,8424,1,1,2,2,2,2,2,2,2,1,1,1,2,1,2,2,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,1,1,2,2,2,1,2,2,1,2,1,2,2,2,2,2,2,1,1,2,1,1,2,2,1,1,1,2,2,1,2,2,2,1,2,2,2,1,2,2,1,2,2,1,1

# PURPOSE: Select loci with minimum unresolved haplotypes and carry out statistics on the cross-overs identified in this script as switches from one haplotype to the other. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

hap.input.path <- "./resolved_haps/by_parent/"
curated.hap.input.path <- "./resolved_haps/by_parent/curated/"
statistics.output.path <- "./recombination_statistics/"
########################################################################################################################################################################################################
chrom.size.df <- read.csv(file = "./scripts_lists/chrom.size.df_Stetv2.csv", header = TRUE, stringsAsFactors = FALSE)
co.count.df <-  read.csv(paste0(statistics.output.path, "co.count.df.csv"), stringsAsFactors = FALSE, header = TRUE)
colnames(co.count.df) <- par.name.vec
row.names(co.count.df) <- chrom.name.vec

# reshape the co.count.df to a form that can be used for plotting as (x,y) coordinates
reshaped.co.count.df <- NULL
for (j in 1:ncol(co.count.df)) {
  reshaped.co.count.df <- rbind(reshaped.co.count.df, cbind(chrom.size.df, rep(colnames(co.count.df)[j]), co.count.df[,j]))
}
colnames(reshaped.co.count.df) <- c(colnames(chrom.size.df), "PAR", "AVG_CO_COUNT")


pdf(file=paste0(statistics.output.path, "CO_vs_chrom_size_vs.pdf"), width=8, height=8)
# regression of average CO counts across halfsib families on chromosome length
plot(reshaped.co.count.df$LENGTH, reshaped.co.count.df$AVG_CO_COUNT, pch = 19, cex=0.75, col = 'grey', xlab = "chromosome size (bp)", ylab = "mean count of cross-overs across half-sib families")
# define the linear model
co.v.size.mod <- lm(AVG_CO_COUNT ~ LENGTH, data = reshaped.co.count.df)
# pot the regression line
abline(co.v.size.mod, lty = 5, lwd = 2, col = 'blue')
# define confidence intervals based on the linear regression model
# CI <- as.data.frame(predict(co.v.size.mod, interval = 'confidence'))
# add the original x-values to the dataframe
# CI$length <- reshaped.co.count.df$LENGTH
# make the lwr and upr confidence interval lines
# lines(CI$length[1:19], round(CI$lwr[1:19], 4), lty = 1, lwd = 2, col = 'black')
# lines(CI$length[1:19], round(CI$upr[1:19], 4), lty = 1, lwd = 2, col = 'black')
mtext(paste0("R^2 = ", round(summary(co.v.size.mod)$adj.r.squared, 3)))
dev.off()


########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################

# barplot(apply(co.count.df, 2, sum), main = 'Average CO per parent for the whole genome', xlab = par.name.vec)
# barplot(apply(co.count.df, 1, mean), main = 'Average cross-over counts by chromosome')

CO.df <- read.csv(paste0(statistics.output.path, "co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# CO.df <- read.csv(paste0("./co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
########################################################################################################################################################################################################
hist(CO.df$CO_SIZE, breaks = 1000, xlab = "Size of CO-region (bps)", ylab = "frequency", main = NULL)
text(x = c(4000000), y = c(7000), labels = (paste0("mean CO size = ", round(mean(CO.df$CO_SIZE))/1000, "Kb")), pos = 3)
text(x = c(4000000), y = c(6500), labels = (paste0("median CO size = ", round(median(CO.df$CO_SIZE))/1000, "Kb")), pos = 3)

quantile(CO.df$CO_SIZE, c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
# 25%       50%       75%       90% 
#   13722.25  30715.50  55600.50 103432.80 
#   
max(CO.df$CO_SIZE, na.rm = T)
# [1] 5338316
# 

# following code block does not run to completion and gives error: Error: cannot allocate vector of size 735.2 Mb
# tmp.CO.df <- CO.df[CO.df$CHROM == chrom.name.vec[3] & CO.df$PAR == par.name.vec[14], ]
# tmp.CO.loci <- NULL
# for (i in 1:nrow(tmp.CO.df)) { # i <- 4
#   if (tmp.CO.df$CO_SIZE[i] != 0) {
#     tmp.CO.loci <- c(tmp.CO.loci, seq(from = tmp.CO.df$LFT_FLNK[i], to = tmp.CO.df$RT_FLNK[i], by = 1))
#   }
# }
# plot(density(tmp.CO.loci))
#   

########################################################################################################################################################################################################

# # creating a dataframe that has data in the form of an input to the boxplot funtion
# co.data.df <- NULL
# for (i in 1:nrow(co.count.df)) { # i <- 1
#   chrom.name.set <- rep(chrom.name.vec[i], ncol(co.count.df))
#   tmp.df <- cbind(CHROM=chrom.name.set, PAR=paste0(rep(c('F_','M_'),each=7), par.name.vec), CO.COUNT=co.count.df[i,])
#   co.data.df <- rbind(co.data.df, tmp.df)
# }
# 
# co.data.df <- as.data.frame(co.data.df, stringsAsFactors = FALSE)
# colnames(co.data.df) <- c("CHROM", "PAR", "CO.COUNT")
# co.data.df[,"CO.COUNT"] <- as.numeric(co.data.df[,"CO.COUNT"])
# write.csv(co.data.df, file = paste0(statistics.output.path, "co.data.df.csv"), quote = FALSE, row.names = FALSE)
# co.data.df <- read.csv(paste0(statistics.output.path, "co.data.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# 
# boxplot(CO.COUNT ~ CHROM, data = co.data.df, xlab = "Chromosome", ylab = "Mean cross-over count")
# boxplot(CO.COUNT ~ PAR, data = co.data.df, xlab = "Half-sib family identity by focal parent", ylab = "Mean cross-over count")

########################################################################################################################################################################################################
chrom.size.df <- read.csv(file = "./scripts_lists/chrom.size.df_Stetv2.csv", header = TRUE, stringsAsFactors = FALSE)

# regression of average CO counts across halfsib families on chromosome length
plot(chrom.size.df[,"LENGTH"], apply(co.count.df, 1, mean), type = "p", col = "blue", xlab = "chromosome size (bp)", ylab = "mean count of cross-overs across half-sib families")
co.size.mod <- lm(apply(co.count.df, 1, mean) ~ chrom.size.df[,"LENGTH"])
sum.co.size.mod <- summary(co.size.mod)
mtext(paste0("R^2 = ", sum.co.size.mod$adj.r.squared))
abline(co.size.mod)

# NOTES:
# Cross-over count has high correlation to the size of the chromosome. In certain species the number of COs do not depend on the chromosome length and rather shows a specific obligatory number of COs.
########################################################################################################################################################################################################
# script for CreateBinSpanCoordinates function
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

win.size <- 250000 # This is the final window size statistical analysis done on.
# assign bin start and end points based on bin size
total.bin.span.coord.df <- NULL
for (chrom.name in chrom.name.vec) {
  bin.span.coord.df.temp <- CreateBinSpanCoordinates(chrom.name = chrom.name, bin.size = win.size, chrom.length = chrom.size.df[chrom.size.df$CHROM == chrom.name, "LENGTH"], shift.size = win.size)
  bin.span.coord.df.temp2 <- cbind(rep(chrom.name, nrow(bin.span.coord.df.temp)), bin.span.coord.df.temp)
  total.bin.span.coord.df <- rbind(total.bin.span.coord.df, bin.span.coord.df.temp2)
}

total.bin.span.coord.df <- as.data.frame(total.bin.span.coord.df, stringsAsFactors = FALSE)
total.bin.span.coord.df[,2:3] <- apply(total.bin.span.coord.df[,2:3], 2, as.numeric)
colnames(total.bin.span.coord.df) <- c("CHROM", "START", "END")
########################################################################################################################################################################################################
# creating a dataframe with COs assigned to non overlapping windows
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
# all.par.df.2 <- read.csv(paste0(statistics.output.path, "co.bin_size_10kb.csv"), stringsAsFactors = FALSE, header = TRUE)

# all.par.df.2 <- read.csv(paste0(statistics.output.path, "co.bin_size_10kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(paste0(statistics.output.path, "co.bin_size_250kb.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = "./co.bin_size_10kb.csv", header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = "./co.bin_size_250kb.csv", header = TRUE, stringsAsFactors = FALSE )
# plot(apply(all.par.df.2[all.par.df.2$CHROM == "Chr01",4:17], 1, mean), typ = "l")

library(RColorBrewer)
colors <- c(brewer.pal(7, "Set3"), brewer.pal(7, "Set2"))
for (par.name in par.name.vec) { # par.name <- par.name.vec[1]
  if (par.name == par.name.vec[1]) {
    plot((all.par.df.2[all.par.df.2$CHROM == chrom.name.vec[4], paste0("X",par.name)]), typ = "l", col = colors[par.name.vec == par.name])
  } else {
    points((all.par.df.2[all.par.df.2$CHROM == chrom.name.vec[4], paste0("X",par.name)]), typ = "l", col = colors[par.name.vec == par.name])
  }
}
########################################################################################################################################################################################################
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

library(lme4)
require(car)
library(MASS)
library(lmtest)

# identify rows (windows) that have all zero entries or an excess of zero entries.
rows.to.exclude <- which(apply(all.par.df.2[,4:17], 1, mean) < 0.3)
all.par.df.3 <- all.par.df.2[-rows.to.exclude, ]

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

########################################################################################################################################################################################################
# difference between the number of windows that show misfit to NB vs poiss
barplot((c(sum(glm.nb.fit.vec < 0.05), sum(glm.poiss.fit.vec < 0.05))), main =  'Negative Binomial vs. Poisson', ylab = 'count of windows that deviate from model') 
misfit.rows <- which(glm.nb.fit.vec < 0.05) # these are the rows that do not fit the NB model with an alpha =< 0.05

null.distribution <- rchisq(length(mod.compar.vec),1) # Deviance between Poiss and NB models is distributed Chi-sq with df = 1
dev.off()
hist(mod.compar.vec, breaks = seq(0, 50, by = 0.5),
     main = expression(paste("Deviance between Negative Binomial and Poisson vs. ", chi^2, "df = 1")),
     xlab = 'Deviance',
     xlim = c(0,40), ylim = c(0,200), 
     col = rgb(1,0,0,0.5), border = 'black') # https://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r
hist(null.distribution, breaks = seq(0, 50, by = 0.5), 
     xlim = c(0,40), ylim = c(0,200), 
     col = rgb(0,0,1,0.5), border = 'white', 
     add = TRUE)
# max(null.distribution)

sum(glm.poiss.fit.vec > 0.05) # 153 # check the number of windows fitiing the model as per arbitrary threshold
sum(glm.nb.fit.vec > 0.05) # 352 # check the number of windows fitiing the model as per arbitrary threshold

lower.tail.cut.off <- quantile(glm.nb.fit.vec, 0.05)
windows.to.test <- which(glm.nb.fit.vec > lower.tail.cut.off)# specify the windows to test based on whether or not fits the NB PDF. 

########################################################################################################################################################################################################
# comparing the fit to two nested models for selected windows to evaluate whether adding Sex as an effect improves model fit
# glm1 -> CO.count ~ Num.ind
# glm2 -> CO.count ~ Num.ind + Sex


lr.stat.vec <- NULL 
sex.beta.vec <- NULL
sex.beta.p.val.vec <- NULL
intcpt.vec <- NULL
mu.vec <- NULL
msme.intcpt.vec <- NULL
num.ind.mod.vec <- NULL

for (i in windows.to.test) { # i <- 823
  sex.assignment <- rep(c(1, 2), each = 7)
  
  num.ind.vec <- NULL
  for (par.name in par.name.vec) { # par.name <- "2365"
    num.ind.vec <- c(num.ind.vec, length(unique(CO.df[CO.df$CHROM == all.par.df.3$CHROM[i] & CO.df$PAR == par.name, "IND"])))
    
  }
  
  
  stat.tab <- as.data.frame(cbind(as.vector(as.matrix(all.par.df.3[i, 4:17])), num.ind.vec, sex.assignment))
  colnames(stat.tab) <- c("CO_COUNTS", "NUM_IND", "SEX")
  tmp.var <- factor(stat.tab$SEX)
  stat.tab$SEX <- tmp.var
  
  glm.1 <- MASS::glm.nb(CO_COUNTS ~ NUM_IND, data = stat.tab)
  glm.2 <- MASS::glm.nb(CO_COUNTS ~ NUM_IND + SEX, data = stat.tab)
  lr.stat <- lrtest(glm.2, glm.1)$Chisq[2]
  lr.stat.vec <- c(lr.stat.vec, lr.stat)
  
  sum.glm.1 <- summary(glm.1)
  sum.glm.2 <- summary(glm.2)
  
  sex.p.val <- sum.glm.2$coefficients[3, 4]
  sex.beta.p.val.vec <- c(sex.beta.p.val.vec, sex.p.val)
  
  sex.beta <- sum.glm.2$coefficients[3,1]
  sex.beta.vec <- c(sex.beta.vec, sex.beta)
  
  intcpt <- exp(coef(glm.2)[1]) # extract intercept value for glm.2 model fit
  intcpt.vec <- c(intcpt.vec, intcpt)
  
  num.ind.mod <- sum.glm.1$coefficients[2,1]
  num.ind.mod.vec <- c(num.ind.mod.vec, num.ind.mod)
  
  nbinom <- MASS::fitdistr(stat.tab$CO_COUNTS, "negative binomial") # estimate the negative binomial parameters
  mu <- nbinom$estimate[[2]]
  mu.vec <- c(mu.vec, mu)
  # qqp(subset.genome.tab$CO_COUNT, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
  
}

plot(exp(num.ind.mod.vec), type = 'l', xlab = "window id", ylab = "exp(coef(NUM_IND))") # shows that the beta-coefficient for number of offspring varies by window, which is not ideal
# plot(intcpt.vec, type = 'l', ylim = c(0,50))
# plot(mu.vec, type = 'l')

##############################################################################################################################################################











