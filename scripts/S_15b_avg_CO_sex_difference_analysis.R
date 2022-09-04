#!/usr/bin/R
# setwd("/group/difazio/populus/gatk-7x7-Stet14/")
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS, for the recombination landscape analysis for 7x7 cross
# DATE: September 20, 2020
# USAGE: Input files include average CO counts by parent and chromosome as well as a dataframe for chromosome physical size and genetc map size. 

# Following is the input file format used 
# CHROM  PAR   IND CO_SIZE DWNSTRM_FLNK UPSTRM_FLNK
# Chr01 1863 24708  103140      4675047     4778187
# Chr01 1863 24708   30914     26411972    26442886
# Chr01 1863 24708   61845     34445658    34507503
# Chr01 1863 24708    5612     46505278    46510890
# Chr01 1863 24709  164382       726848      891230
# Chr01 1863 24709  153359      2513166     2666525

# PURPOSE: Model CO-count given the SEX, chromosome size and other covariates. 
########################################################################################################################################################################################################
rm(list = ls())
library(dplyr)
library(stringr)
library(lme4)
require(car)
library(MASS)
# library(lmtest)
# library(glmmTMB)

########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")
input.path <- "./statistics/"
statistics.output.path <- "./statistics/gross_stat_plots/"

########################################################################################################################################################################################################
# reading in the chromosome size data and the CO details
# chrom.size.df <- read.csv(file = "./scripts_lists/chrom.size.df.csv", header = TRUE, stringsAsFactors = FALSE)
# chrom.size.df <- read.csv(file = paste0(input.path, "chrom.size.df.csv"), header = TRUE, stringsAsFactors = FALSE)
chrom.size.df <- read.csv(file = paste0(input.path, "chrom.size.df_Stetv2.csv"), header = TRUE, stringsAsFactors = FALSE)
CO.df <- read.csv(paste0(input.path, "./co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# any((CO.df$CO_SIZE == 0))



# reading in the total number of CO counts and the family-size data for each half-sib family
# num.COs.df <-  NULL
# num.ind.df <- NULL
# for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[1]
#   tmp.num <- NULL
#   tmp.num.2 <- NULL
#   for (par.name in par.name.vec) { # par.name <- par.name.vec[10]
#     
#     tmp.tab <- CO.df[CO.df$CHROM == chrom.name & CO.df$PAR == par.name, ]
#     tmp.num <- c(tmp.num, length(unique(tmp.tab$IND)))
#     tmp.num.2 <- c(tmp.num.2, sum(tmp.tab$CO_SIZE != 0))
#   }
#   num.ind.df <- rbind(num.ind.df, tmp.num)
#   num.COs.df <- rbind(num.COs.df, tmp.num.2)
# }
# num.ind.df <- as.data.frame(num.ind.df, stringsAsFactors = FALSE)
# num.COs.df <- as.data.frame(num.COs.df, stringsAsFactors = FALSE)
# colnames(num.ind.df) <- par.name.vec
# colnames(num.COs.df) <- par.name.vec
# write.csv(num.ind.df, file = paste0(input.path, "number_of_individuals_per_par_and_chrom.csv"), quote = FALSE, row.names = FALSE)
num.ind.df <- read.csv(paste0(input.path, "number_of_individuals_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
# write.csv(num.COs.df, file = paste0(input.path, "number_of_crossovers_per_par_and_chrom.csv"), quote = FALSE, row.names = FALSE)
num.COs.df <- read.csv(paste0(input.path, "number_of_crossovers_per_par_and_chrom.csv"), header = TRUE, stringsAsFactors = FALSE)
colnames(num.ind.df) <- colnames(num.ind.df) %>% str_replace("^X", "")
colnames(num.COs.df) <- colnames(num.COs.df) %>% str_replace("^X", "")



# Check the total CO count differences between sexes
tot.COs <- as.data.frame(t(apply(num.COs.df, 1, function(x) {
  f.sum <- sum(x[1:7])
  m.sum <- sum(x[8:14])
  return(c(f.sum, m.sum))
})), stringsAsFactors=FALSE )
colnames(tot.COs) <- c("FEM", "MALE")
rownames(tot.COs) <- chrom.name.vec
sum(apply(tot.COs, 1, sum)) # 38846



# Chi-square test for each chromosome to check for equality of CO numbers since the total number of participating offspring are similar.
critical.p.val <- 0.05/length(chrom.name.vec)
chi.test.df <- NULL
for (chrom in 1:19){
  chi.test.sum <- chisq.test(c(tot.COs[chrom,"FEM"], tot.COs[chrom,"MALE"]))
  chi.test.df <- rbind(chi.test.df, c(chi.test.sum$statistic, chi.test.sum$p.value))
}
colnames(chi.test.df) <- c("CHI.STAT", "CHI.P.VAL")
chi.test.df <- cbind(CHROM=chrom.name.vec, chi.test.df)
# in every chromosome, Total male CO count is higher. Here the Total fam sizes for male and female groups are the same at 829 individuals
# CHROM   CHI.STAT            CHI.P.VAL
# "Chr01" "19.982898289829"   "7.81378945986676e-06"
# "Chr02" "11.7333333333333"  "0.000613905375006302"
# "Chr03" "14.9038718291055"  "0.000113130556138208"
# "Chr04" "2.28773168578994"  "0.130400373618151"
# "Chr05" "0.934604904632153" "0.333669195320745"
# "Chr06" "9.46291208791209"  "0.00209667947621377"
# "Chr07" "0.821292775665399" "0.36480245958835"
# "Chr08" "3.39526860522425"  "0.0653837132261563"
# "Chr09" "3.03055378739656"  "0.0817101234978999"
# "Chr10" "13.7831060272767"  "0.000205172933512095"
# "Chr11" "1.30478440637921"  "0.25334114264351"
# "Chr12" "2.70659833440103"  "0.0999338807608696"
# "Chr13" "2.63284457478006"  "0.104674029918986"
# "Chr14" "17.1600928074246"  "3.43578973142228e-05"
# "Chr15" "3.59789750328515"  "0.0578526944272907"
# "Chr16" "9.26640926640927"  "0.00233394823164995"
# "Chr17" "8.34030045721751"  "0.00387751847126799"
# "Chr18" "0.835265405735204" "0.360754387509099"
# "Chr19" "12.5533141210375"  "0.00039550424827411"


# In order to apply the linear models the data needs to be arranged in a different table format.
reshaped.co.count.df <- NULL
for (j in 1:ncol(num.COs.df)) {reshaped.co.count.df <- rbind(reshaped.co.count.df, cbind(chrom.size.df, rep(colnames(num.COs.df)[j]), num.COs.df[,j], num.ind.df[,j]))}
colnames(reshaped.co.count.df) <- c(colnames(chrom.size.df), "PAR", "CO_COUNT", "FAM_SIZE")
reshaped.co.count.df$SEX <- as.factor(rep(c("F", "M"), each = 133))
reshaped.co.count.df$CHROM <- as.factor(reshaped.co.count.df$CHROM)
sum(reshaped.co.count.df$CO_COUNT) # 38846
reshaped.co.count.df$CO_RATE_PER_MEIOSIS <- reshaped.co.count.df$CO_COUNT / reshaped.co.count.df$FAM_SIZE


jpeg(file = paste0(statistics.output.path, "CO_counts_by_SEX_per_CHROM.jpg"), width = 2000, height = 500, quality = 5000)

dev.new()
op <- par(mfrow = c(1,19))
col.vec <- c("red", "blue")
for (chrom.name in chrom.name.vec) {
  boxplot(reshaped.co.count.df$CO_RATE_PER_MEIOSIS[reshaped.co.count.df$CHROM == chrom.name] ~ reshaped.co.count.df$SEX[reshaped.co.count.df$CHROM == chrom.name], 
          main = chrom.name, ylab = "CO rate per meiosis",
          col = col.vec)
}
par(op)
dev.off()

##############################################################################################################################################################
# Analyze distribution of CO-count data
# install.packages('fitdistrplus', dependencies = TRUE) # https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
# library(fitdistrplus)
# # The emperical PDF an the CDF of the CO_COUNT per half-sib family
# jpeg(file = paste0(statistics.output.path, "Emperical_density_of_total_CO_events_per_hs_fam.jpg"), width = 16, height = 9, unit = "in", res = 600)
# op <- par(mfrow = c(2, 2))
# plotdist(reshaped.co.count.df$CO_COUNT, histo = TRUE, demp = TRUE)
# par(op)
# dev.off()
# 
# fnb <- fitdist(reshaped.co.count.df$CO_COUNT, distr = "nbinom", method = 'mle')
# fp <- fitdist(reshaped.co.count.df$CO_COUNT, distr = "pois", method = 'mle')
# fg <- fitdist(reshaped.co.count.df$CO_COUNT, distr = "geom", method = 'mle')
# 
# jpeg(file = paste0(statistics.output.path, "Emperical_vs_parametric_distributions_of_total_CO_events_per_hs_fam.jpg"), width = 16, height = 9, unit = "in", res = 600)
# op <- par(mfrow = c(1, 2))
# plot.legend <- c("nbinom", "pois", "geom")
# denscomp(list(fnb, fp, fg), legendtext = plot.legend)
# qqcomp(list(fnb, fp, fg), legendtext = plot.legend)
# par(op)
# dev.off()
# 
# jpeg(file = paste0(statistics.output.path, "Emperical_vs_parametric_CDFs_of_total_CO_events_per_hs_fam.jpg"), width = 16, height = 9, unit = "in", res = 600)
# op <- par(mfrow = c(1, 2))
# plot.legend <- c("nbinom", "pois", "geom")
# cdfcomp(list(fnb, fp, fg), legendtext = plot.legend)
# ppcomp(list(fnb, fp, fg), legendtext = plot.legend)
# par(op)
# dev.off()

# Goodness of fit statistics for the different theoretical models for model selection
# gof <- gofstat(list(fnb, fp, fg), fitnames = c("nbinom", "pois", "geom"))
# gof # lower the AIC, better the fit to model
##############################################################################################################################################################
# Following block of code is a manual curation of the dta for adjusting known errors in the geneic maps.
reshaped.co.count.df.mod <- reshaped.co.count.df[-229, ] # removing the known outlier in Chr01 GW-6909
idx.to.rmv <- which(reshaped.co.count.df.mod$CHROM==chrom.name.vec[17])# remove Chr17 data since structural re-arrangement was found on Stet-14 and may bias inference
reshaped.co.count.df.mod <- reshaped.co.count.df.mod[-idx.to.rmv,]


# Prediction model for CO_COUNT using covariates
lm.1 <- lm(CO_COUNT ~ scale(LENGTH) + scale(FAM_SIZE), data = reshaped.co.count.df.mod)
sum.lm.1 <- summary(lm.1)
plot(fitted(lm.1), resid(lm.1, type = "pearson"))
lm.2 <- lm(CO_COUNT ~ scale(LENGTH) + scale(FAM_SIZE) + SEX, data = reshaped.co.count.df.mod)
sum.lm.2 <- summary(lm.2)
plot(fitted(lm.2), resid(lm.2, type = "pearson"))


glm.poiss.1 <- stats::glm(CO_COUNT ~ SEX*CHROM + scale(FAM_SIZE), data = reshaped.co.count.df.mod, family = poisson(link = "log")) 
summary(glm.poiss.1)
glm.poiss.2 <- stats::glm(CO_COUNT ~ SEX*scale(LENGTH) + scale(FAM_SIZE), data = reshaped.co.count.df.mod, family = poisson(link = "log")) 
summary(glm.poiss.2)
glm.poiss.3 <- stats::glm(CO_COUNT ~ PAR*scale(LENGTH) + scale(FAM_SIZE) , data = reshaped.co.count.df.mod, family = poisson(link = "log")) 
summary(glm.poiss.3)
glmer.poiss.1 <- glmer(CO_COUNT ~ SEX*scale(LENGTH) + scale(FAM_SIZE) + (1|SEX:PAR), data = reshaped.co.count.df.mod, family=poisson(link = "log")) # testing the CHROM as a variable with LENGTH
summary(glmer.poiss.1, cor = FALSE)
plot(glmer.poiss.1)
glmer.poiss.2 <- glmer(CO_COUNT ~ SEX*CHROM + scale(FAM_SIZE) + (1|SEX:PAR), data = reshaped.co.count.df.mod, family=poisson(link = "log")) # testing the CHROM as a variable with LENGTH
sum.glmer.poiss.2 <- summary(glmer.poiss.2, cor = FALSE)
plot(glmer.poiss.2)
glmer.poiss.3 <- glmer(CO_COUNT ~ SEX*CHROM + scale(LENGTH) + scale(FAM_SIZE) + (1|SEX:PAR), data = reshaped.co.count.df.mod, family=poisson(link = "log")) # testing the CHROM as a variable with LENGTH
sum.glmer.poiss.3 <- summary(glmer.poiss.2, cor = FALSE)



jpeg(file = paste0(statistics.output.path, "Residuals_vs_fit_stat_model.jpg"), width = 1000, height = 500, quality = 5000)
plot(glmer.poiss.2, xlab="Fitted CO counts", ylab="Pearson residuals", main="Statistical model validation; glmer.poiss.2")
dev.off()


# jpeg(file = paste0(statistics.output.path, "Residuals_vs_covariates_stat_model.jpg"), width = 1000, height = 500, quality = 5000)
pdf(file = paste0(statistics.output.path, "Residuals_vs_covariates_stat_model.pdf"), width = 10, height = 10)
op <- par(mfrow=c(2,2), mar=c(5.1,4.1,1,2.1))
plot(reshaped.co.count.df.mod$CHROM, resid(sum.glmer.poiss.2, type='Pearson'), xlab='LG / Chromosome', ylab='Pearson-residuals', las=2)
plot(reshaped.co.count.df.mod$PAR, resid(sum.glmer.poiss.2, type='Pearson'), xlab='Half-sib family name', ylab='Pearson-residuals', las=2)
plot(reshaped.co.count.df.mod$LENGTH, resid(sum.glmer.poiss.2, type='Pearson'), xlab='Physical length of the Chromosome in base pairs', ylab='Pearson-residuals')
plot(reshaped.co.count.df.mod$SEX, resid(sum.glmer.poiss.2, type='Pearson'), xlab='Sex', ylab='Pearson-residuals')
par(op)
dev.off()

pdf(file = paste0(statistics.output.path, "fitted_vs_residuals.pdf"), width = 5, height = 5)
plot(glmer.poiss.2)
dev.off()


glm.nb.1 <- MASS::glm.nb(CO_COUNT ~ SEX*CHROM + scale(FAM_SIZE), link = "log", data = reshaped.co.count.df.mod)
summary(glm.nb.1, cor = FALSE)
glm.nb.2 <- MASS::glm.nb(CO_COUNT ~ SEX + CHROM + scale(FAM_SIZE), link = "log", data = reshaped.co.count.df.mod)
summary(glm.nb.2, cor = FALSE)
glm.nb.3 <- MASS::glm.nb(CO_COUNT ~ SEX*scale(LENGTH) + scale(FAM_SIZE), link = "log", data = reshaped.co.count.df.mod)
summary(glm.nb.3, cor = FALSE)
glm.nb.4 <- MASS::glm.nb(CO_COUNT ~ SEX + scale(LENGTH) + scale(FAM_SIZE), link = "log", data = reshaped.co.count.df.mod)
summary(glm.nb.4, cor = FALSE)
glmer.nb.1 <- glmer.nb(CO_COUNT ~ SEX + scale(LENGTH) + scale(FAM_SIZE) + (1|SEX:PAR), data = reshaped.co.count.df.mod) # testing the CHROM as a variable with LENGTH
summary(glmer.nb.1, cor = FALSE)
plot(glmer.nb.1)


# statistical model selection
AIC(lm.1, lm.2, glm.poiss.1, glm.poiss.2, glm.poiss.3, glmer.poiss.1, glmer.poiss.2, glm.nb.1, glm.nb.2, glm.nb.3, glm.nb.4, glmer.nb.1)
anova(glm.poiss.1, glm.poiss.2, glmer.poiss.1, glmer.poiss.2, glm.nb.1, glm.nb.2, glm.nb.3, glm.nb.4, glmer.nb.1)



# Since the best model is glmer.poiss.2 here I am checking whether chromsome size is an underlying covariate leading the effect.
glmer.estimates <- sum.glmer.poiss.2$coefficients[3:19,1]
glmer.estimates <- c(0, glmer.estimates)
chrom.size.vec <- chrom.size.df$LENGTH[-17] # chrom.size.vec <- chrom.size.df$LENGTH[-c(1, 17)]
chrom.mod <- lm(glmer.estimates ~ chrom.size.vec)
sum.chrom.mod <- summary(chrom.mod)


jpeg(file = paste0(statistics.output.path, "Chrom_effect_size_on_length_stat_model.jpg"), width = 1000, height = 500, quality = 5000)
op <- par(mfrow = c(1,1))
plot(glmer.estimates ~ chrom.size.vec, xlab='LG size in bps', ylab='GLMM(Poiss) estimates for CHROM effect', main=paste0('R^2=', round(sum.lm.length.estim$adj.r.squared, digits=2)))
text(x=chrom.size.vec, y=glmer.estimates, labels=chrom.name.vec[-17], offset=1, pos=2, cex=0.75)
abline(sum.chrom.mod$coefficients[1,1], sum.chrom.mod$coefficients[2,1], lty=2, col='blue')
par(op)
dev.off()
# Chromosome fixed effect is largely determined by the size of the chromosome as expected.

##############################################################################################################################################################
# # Model adequecy (overfitting test) validated using leave-one-out cross validation
# reshaped.co.count.df.mod.2 <- reshaped.co.count.df.mod#[reshaped.co.count.df.mod$CHROM != chrom.name.vec[1], ]
# pred.vec <- NULL
# obs.vec <- NULL
# for (par.name in par.name.vec) {
#   mod.data <- reshaped.co.count.df.mod.2[(reshaped.co.count.df.mod.2$PAR != par.name), ]
#   valid.data <- reshaped.co.count.df.mod.2[(reshaped.co.count.df.mod.2$PAR == par.name), ]
#  
#   tmp.glmer.poiss.2 <- glmer(CO_COUNT ~ SEX*CHROM + scale(FAM_SIZE) + (1|SEX:PAR), data = mod.data, family=poisson(link = "log"))
#   sum.tmp.glmer.poiss.2 <- summary(tmp.glmer.poiss.2, cor = FALSE)
#   
#   predictions <-  round(predict(object=tmp.glmer.poiss.2, newdata=valid.data, re.form=NA, type="response"))
#   observations <- valid.data$CO_COUNT
#   
#   pred.vec <- c(pred.vec, predictions)
#   obs.vec <- c(obs.vec, observations)
# }
# cor.result <- cor.test(obs.vec, pred.vec, method='pearson', conf.level=0.95)
# 
# # cor.result$conf.int[1:2]
# plot(pred.vec, obs.vec, xlab="predicted CO count", ylab='observed CO count', 
#      main=paste0("Leave-one-parent-out cross-validation -- ", "Pearson-R^2=", round(cor.result$estimate, digits=2), 
#                 "(95% CI: ",round(cor.result$conf.int[1] ,digits=2), "- " ,round(cor.result$conf.int[2] ,digits=2),  ") :p-value=", formatC(cor.result$p.value, format = "e", digits = 2)))

########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################






# MISC....
# # Poisson model
# glm.poiss <- stats::glm(CO_COUNT ~ SEX + FAM_SIZE + LENGTH, data = reshaped.co.count.df, family = poisson(link = "log")) 
# summary(glm.poiss)
# 100*(glm.poiss$null.deviance - glm.poiss$deviance) / glm.poiss$null.deviance # Pseudo R^2 - this is the equivalent to the R^2 in general linear model 
# drop1(glm.poiss, test = "Chi") # Model selection and testing explanatory variables one by one 
# anova(glm.poiss)
# theta <- glm.poiss$deviance/glm.poiss$df.residual # since in the Poisson model the over-dispersion parameter is assumed 1 and this needs to be tested before accepting these results. 
# # this needs to be in and around 1 to accept a Poisson model # There is overdispersion in the above model and the readons could be missing covariates or interactions, outliers in the response, 
# # variable, non-linear effects of covariates entered as linear terms in the systematic part of the model, and choice of the wrong link function
# 
# # quasi-Poisson model assumes a theta that is != 1
# glm.q.poiss <- stats::glm(CO_COUNT ~ SEX*FAM_SIZE + LENGTH, data = reshaped.co.count.df, family = "quasipoisson") 
# summary(glm.q.poiss)
# 100*(glm.q.poiss$null.deviance - glm.q.poiss$deviance) / glm.q.poiss$null.deviance # Pseudo R^2 - this is the equivalent to the R^2 in general linear model 
# drop1(glm.q.poiss, test = "F") # Model selection and testing explanatory variables one by one 
# # Model validation
# resids.p <- resid(glm.q.poiss, type = "pearson")
# resids.d <- resid(glm.q.poiss, type = "deviance")
# mu <- predict(glm.q.poiss, type = "response")
# errors <- reshaped.co.count.df$CO_COUNT - mu
# pearson.resids <- errors / sqrt(4.190437 * mu)
# op <- par(mfrow = c(2, 2))
# plot(x = mu, y = errors, main = "Response residuals")
# plot(x = mu, y = resids.p, main = "Pearson residuals")
# plot(x = mu, y = pearson.resids, main = "Pearson residuals")
# plot(x = mu, y = resids.d, main = "Deviance residuals")
# par(op)
# # dev.off()
# 
# # fitting a negative-binomial GLM 
# reshaped.co.count.df$SCALED_LENGTH <- scale(reshaped.co.count.df$LENGTH, center = TRUE)
# glm.nb <- MASS::glm.nb(CO_COUNT ~ SEX + FAM_SIZE + SCALED_LENGTH, link = "log", data = reshaped.co.count.df)
# summary(glm.nb, cor = FALSE)
# boxplot(resid(glm.nb) ~ reshaped.co.count.df$CHROM)
# op <- par(mfrow=c(2,2))
# plot(glm.nb)
# par(op)
# # removing the outlier 6909
# reshaped.co.count.df.mod <- reshaped.co.count.df[-229, ]
# idx.to.rmv <- which(reshaped.co.count.df.mod$CHROM==chrom.name.vec[17])
# reshaped.co.count.df.mod <- reshaped.co.count.df.mod[-idx.to.rmv,]
# # removing Chr17 data
# glm.nb <- MASS::glm.nb(CO_COUNT ~ SEX + FAM_SIZE + SCALED_LENGTH, link = "log", data = reshaped.co.count.df.mod)
# summary(glm.nb, cor = FALSE)
# # model selection processs
# # This drops one variable at a time from te model and compares it to the full model with a Chi-sq goodness of fit test with df = (difference in df between models)
# # <none> here in this output refers to the full model
# drop1(glm.nb, test = "Chi") 
# step(glm.nb)
# stepAIC(glm.nb)
# # With the above step all the variables used in the model were significant and none were dropped. 
# 100 * (glm.nb$null.deviance - glm.nb$deviance) / glm.nb$null.deviance # This is the closest metric to a R^2 we can get using a GLM model
# # reality check that there is no over-dispersion and model assumptions are valid.
# theta <- glm.nb$deviance / glm.nb$df.residual
# (theta)
# # Model validation
# resids.p <- resid(glm.nb, type = "pearson")
# resids.d <- resid(glm.nb, type = "deviance")
# mu <- predict(glm.nb, type = "response")
# errors <- reshaped.co.count.df$CO_COUNT - mu
# pearson.resids <- errors / sqrt(4.190437 * mu)
# op <- par(mfrow = c(1, 1))
# plot(x = mu, y = errors, main = "Response residuals")
# plot(x = mu, y = resids.p, main = "Pearson residuals")
# plot(x = mu, y = pearson.resids, main = "Pearson residuals")
# plot(x = mu, y = resids.d, main = "Deviance residuals")
# par(op)
# dev.off()
# 
# 
# 
# # with chromosome interaction
# glm.nb.1 <- MASS::glm.nb(CO_COUNT ~ SEX*CHROM + scale(FAM_SIZE), link = "log", data = reshaped.co.count.df.mod)
# summary(glm.nb.1, cor = FALSE)
# boxplot(resid(glm.nb) ~ reshaped.co.count.df.mod$CHROM)
# op <- par(mfrow = c(2,2))
# plot(glm.nb)
# par(op)
# 
# # fitting a negative-binomial GLM with PAR parameter instead of SEX
# glm.nb.2 <- MASS::glm.nb(CO_COUNT ~ PAR + scale(FAM_SIZE) + scale(LENGTH) + SEX, link = "log", data = reshaped.co.count.df.mod)
# summary(glm.nb.2, cor = FALSE)
# boxplot(resid(glm.nb.2) ~ reshaped.co.count.df.mod$CHROM)
# drop1(glm.nb.2)
# op <- par(mfrow = c(2,2))
# plot(glm.nb.2)
# par(op)
# 
# # fitting a negative-binomial GLM with CHROM*SEX
# glm.nb.interaction <- MASS::glm.nb(CO_COUNT ~ SEX * SCALED_LENGTH + FAM_SIZE, link = "log", data = reshaped.co.count.df.mod)
# summary(glm.nb.interaction, cor = FALSE)
# op <- par(mfrow = c(2,2))
# plot(glm.nb.interaction)
# par(op)
# boxplot(resid(glm.nb.interaction) ~ reshaped.co.count.df.mod$CHROM)
# # lrtest(glm.nb.2, glm.nb.interaction)
# 
# # fitting PAR as a mixed effect with a GLMM
# glm.nb.w.par <- glmer.nb(CO_COUNT ~ scale(LENGTH) + scale(FAM_SIZE) + SEX + (1|SEX:PAR), data = reshaped.co.count.df.mod) # testing the CHROM as a variable with LENGTH
# summary(glm.nb.w.par, cor = FALSE)
# glm.nb.w.par.poiss <- glmer(CO_COUNT ~ scale(LENGTH) + scale(FAM_SIZE) + SEX + (1|PAR), data = reshaped.co.count.df.mod, family=poisson(link = "log")) # testing the CHROM as a variable with LENGTH
# summary(glm.nb.w.par.poiss, cor = FALSE)
# 
# op <- par(mfrow = c(2,2))
# plot(glm.nb.w.par)
# plot(glm.nb.w.par.poiss)
# par(op)
# boxplot(resid(glm.nb.w.par) ~ reshaped.co.count.df$CHROM)
# 
# # fitting the same model as aboe after removing outlier (Chr01 in "6909")
# reshaped.co.count.df.2 <- reshaped.co.count.df[-(reshaped.co.count.df$CHROM == chrom.name.vec[1] & reshaped.co.count.df$PAR == par.name.vec[13]), ]
# # reshaped.co.count.df.2 <- reshaped.co.count.df[-(reshaped.co.count.df$CHROM == chrom.name.vec[1]), ]
# # reshaped.co.count.df.2 <- reshaped.co.count.df[-(reshaped.co.count.df$PAR == par.name.vec[13]), ]
# 
# glm.nb.w.par.2 <- glmer.nb(CO_COUNT ~ SCALED_LENGTH + FAM_SIZE + SEX + (1|PAR), data = reshaped.co.count.df.2)
# summary(glm.nb.w.par.2, cor = FALSE)
# boxplot(resid(glm.nb.w.par.2) ~ reshaped.co.count.df.2$CHROM)
# plot(glm.nb.w.par.2)
# 
# # Model comparison using AIC
# anova(glm.poiss, glm.q.poiss, glm.nb, glm.nb.1, glm.nb.2, glm.nb.interaction, glm.nb.w.par, glm.nb.w.par.2)
# AIC(glm.nb, glm.nb.1, glm.nb.2, glm.nb.interaction, glm.nb.w.par)
########################################################################################################################################################################################################