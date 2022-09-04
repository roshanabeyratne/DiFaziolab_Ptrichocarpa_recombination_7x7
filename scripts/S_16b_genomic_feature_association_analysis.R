#!/usr/bin/R
########################################################################################################################################################################################################
# METADATA: Written by CRA & DMS for the recombination landscape analysis for 7x7 cross
# DATE: January 14, 2021; last modified February 14 2021
# USAGE: Input file format include a genomic feature count per halfsib family per 30 kb window.
# PURPOSE: evaluate association between various genomic features and the correlation count
########################################################################################################################################################################################################
setwd("/Users/cabeyrat/Google Drive/Shared drives/DiFazioLab/BESC/BESC_DATA/7x7/Stettler_14_linkage_mapping")
rm(list = ls())
library(dplyr)
library(stringr)
library(lme4)
library(caret)
library(ggplot2)
library(gridExtra)
library(psych)

########################################################################################################################################################################################################
par.name.vec <- c("1863", "1909", "1950", "2048", "2066", "2283", "4593", "2365", "2393", "2515", "2572", "2683", "6909", "7073" )
chrom.name.vec <- c( "Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19")

input.path <- "./genomic_features/output/"
statistics.input.path <- "./statistics/"
output.path <- "./genomic_features/analysis_output/"

########################################################################################################################################################################################################
# win.size <- "250kb"
# win.size <- "1Mb"
win.size <- "960kb"

# read-in the input files with the features to analyze
copia <- read.table(file=paste0(input.path, win.size, "_window_copia_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
gypsy <- read.table(file=paste0(input.path, win.size, "_window_gypsy_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
simp.rep <- read.table(file=paste0(input.path, win.size, "_window_Simple_repeat_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
genes <- read.table(file=paste0(input.path, win.size, "_window_gene_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
nuc.content <- read.table(file=paste0(input.path, win.size, "_window_NUC_content_output_7x7"), stringsAsFactors=FALSE, header=FALSE)
colnames(nuc.content) <- c("CHROM", "START" ,"END", "pct_at", "pct_gc", "num_A", "num_C", 'num_G', "num_T", "num_N", "num_oth", "seq_len")
tail(nuc.content)



# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_30kb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_250kb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
# all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_1Mb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_960kb_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )
colnames(all.par.df.2) <- colnames(all.par.df.2) %>% str_replace("^X", "")

tot.co.count <- apply(all.par.df.2[-(1:3)], 1, sum)
gen.feat.df <- cbind(nuc.content, COPIA=copia[,4], GYPSY=gypsy[,4], SIMP_REPEAT=simp.rep[,4], GENE=genes[,4], CUM_CO_COUNT=tot.co.count, log_CUM_CO_COUNT=log10(tot.co.count+0.05), WIN_IDX=1:nrow(nuc.content))
# write.csv(gen.feat.df, file = paste0(output.path, win.size, "_genomic_features_Stet_v2_raw.csv"), quote = FALSE, row.names = FALSE)



# # Visualize genomic patterns of LTRs, GC content and gene content.
jpeg(file = paste0(output.path, "genomic_features_hist_plot_960kb.jpg"), width = 1000, height = 1000)
op <- par(mfrow=c(7,1))
plot( scale(CUM_CO_COUNT) ~ WIN_IDX, data = gen.feat.df, type='h')
plot( scale(pct_gc) ~ WIN_IDX, data = gen.feat.df, type='h')
plot( scale(pct_at) ~ WIN_IDX, data = gen.feat.df, type='h')
plot( scale(COPIA) ~ WIN_IDX, data = gen.feat.df, type='h')
plot( scale(GYPSY) ~ WIN_IDX, data = gen.feat.df, type='h')
plot( scale(SIMP_REPEAT) ~ WIN_IDX, data = gen.feat.df, type='h')
plot( scale(GENE) ~ WIN_IDX, data = gen.feat.df, type='h')
par(op)
dev.off()


# # Visualize genomic patterns of LTRs, GC content and gene content.
# p0 <- ggplot(gen.feat.df, aes(x=WIN_IDX, y=CUM_CO_COUNT)) + 
#   geom_area(fill='black') + xlab("") + 
#   coord_cartesian(ylim=c(0,200)) + 
#   ylab("CO count") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# p1 <- ggplot(gen.feat.df, aes(x=WIN_IDX, y=pct_gc)) + 
#   geom_area(fill='green') + xlab("") + 
#   coord_cartesian(ylim=c(0.3,0.35)) + 
#   ylab("percent GC content") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# p2 <- ggplot(gen.feat.df, aes(x=WIN_IDX, y=GENE)) + 
#   geom_area(fill='grey') + xlab("") + 
#   coord_cartesian(ylim=c(0,150)) + 
#   ylab("gene count") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# p3 <- ggplot(gen.feat.df, aes(x=WIN_IDX, y=COPIA)) + 
#   geom_area(fill='blue') + xlab("") + 
#   coord_cartesian(ylim=c(0,150)) + 
#   ylab("Copia count") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# p4 <- ggplot(gen.feat.df, aes(x=WIN_IDX, y=GYPSY)) + 
#   geom_area(fill='red') + xlab("") + 
#   coord_cartesian(ylim=c(0,350)) + 
#   ylab("Gypsy count") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# p5 <- ggplot(gen.feat.df, aes(x=WIN_IDX, y=SIMP_REPEAT)) + 
#   geom_area(fill=col5) + xlab("genome-wide window index (960 kbps non-overlapping)") + 
#   coord_cartesian(ylim=c(0,400)) + 
#   ylab("simple repeats count") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# grid.arrange(p0, p2, p5, p4, p1, p3, nrow=6)
# 


# pairs plot with R^2
# jpeg(file = paste0(output.path, "genomic_features_pairs_plot_960kb.jpg"), width = 1000, height = 1000)
pdf(file = paste0(output.path, "genomic_features_pairs_plot_960kb.pdf"), width = 12, height = 12)
pairs.panels(gen.feat.df[, -c(1:3, 6:12, 18,19)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = FALSE,  # show density plots
             scale =F,
             pch=1,
             col='grey',
             ellipses = TRUE # show correlation ellipses
)
dev.off()


# # permute the CO counts in genomic windows
# permuted.gen.feat.df <- gen.feat.df
# permuted.gen.feat.df$CUM_CO_COUNT <- sample(gen.feat.df$CUM_CO_COUNT, size=length(gen.feat.df$CUM_CO_COUNT), replace=TRUE)


hist(gen.feat.df$CUM_CO_COUNT, breaks=10)
hist(gen.feat.df$COPIA, breaks=10)
hist(gen.feat.df$GYPSY, breaks=10)
hist(gen.feat.df$SIMP_REPEAT, breaks=10)
hist(gen.feat.df$GENE, breaks=10)



# As per Steve/David, the centromere location and the other genomic correlates have higher multicolinearity in modelling. Therefore, following is an excercise to remove the centromeric locations from 
# the dataset.
op <- par(mfrow=c(4,5), mar=c(1,1,1,1))
for (chrom.name in chrom.name.vec) {
  select.idx <- which(gen.feat.df$CHROM == chrom.name)
  plot(gen.feat.df$WIN_IDX[select.idx], gen.feat.df$CUM_CO_COUNT[select.idx], type='l', xlab=chrom.name, ylab='CUM_CO_COUNT')
}
par(op)
# hist(gen.feat.df$CUM_CO_COUNT, breaks=1000)
# quantile(gen.feat.df$CUM_CO_COUNT, c(0.05))



# removing windows that have < 50 CUM_CO_COUNT to make a new dataset. This is done in order to remove the centromeres and the telomeres from the dataset. Near centromeres and telomeres
# the CO counts behave according to a different set of rules and this is where most of the LTRs are harbored.
filtered.idx <- NULL
for (chrom.name in chrom.name.vec) {
  tmp.idx <- which(gen.feat.df$CHROM == chrom.name & gen.feat.df$CUM_CO_COUNT >= 50)
  filtered.idx <- c(filtered.idx, tmp.idx)
}
gen.feat.df <- gen.feat.df[filtered.idx, ]
# write.csv(gen.feat.df, file = paste0(output.path, win.size, "_genomic_features_Stet_v2.csv"), quote = FALSE, row.names = FALSE)
# gen.feat.df <- read.csv(file = paste0(output.path, win.size, "_genomic_features_Stet_v2.csv"), header = TRUE, stringsAsFactors = FALSE )



# considering linear models
lm.mod.1 <- lm(CUM_CO_COUNT ~ 1 + scale(pct_at), data=gen.feat.df) # not siginificant
summary(lm.mod.1)
lm.mod.2 <- lm(CUM_CO_COUNT ~ 1 + scale(pct_gc), data=gen.feat.df) # significant
summary(lm.mod.2)
lm.mod.3 <- lm(CUM_CO_COUNT ~ 1 + scale(GENE), data=gen.feat.df) # highly significant
summary(lm.mod.3)
lm.mod.4 <- lm(CUM_CO_COUNT ~ 1 + scale(COPIA), data=gen.feat.df) # highly significant
summary(lm.mod.4)
lm.mod.5 <- lm(CUM_CO_COUNT ~ 1 + scale(GYPSY), data=gen.feat.df) # highly significant
summary(lm.mod.5)
lm.mod.6 <- lm(CUM_CO_COUNT ~ 1 + scale(SIMP_REPEAT), data=gen.feat.df) # highly significant
summary(lm.mod.6)
AIC(lm.mod.1, lm.mod.2, lm.mod.3, lm.mod.4, lm.mod.5, lm.mod.6)


lm.mod.7 <- lm(CUM_CO_COUNT ~ 1 + scale(pct_gc) + scale(GENE) + scale(COPIA) + scale(GYPSY) + scale(SIMP_REPEAT), data=gen.feat.df) # GENE and SIMP_REPEAT not significant
summary(lm.mod.7)
lm.mod.8 <- lm(CUM_CO_COUNT ~ 1 + scale(pct_gc) + scale(COPIA) + scale(GYPSY) + scale(SIMP_REPEAT), data=gen.feat.df) # SIMP_REPEAT not that significant
summary(lm.mod.8)
lm.mod.9 <- lm(CUM_CO_COUNT ~ 1 + scale(pct_gc) + scale(COPIA) + scale(GYPSY), data=gen.feat.df) # all significant
summary(lm.mod.9)


lm.mod.10 <- lm(CUM_CO_COUNT ~ 1 + scale(pct_gc) + scale(GENE) + scale(SIMP_REPEAT) + scale(COPIA), data=gen.feat.df) # all significant
summary(lm.mod.10)
lm.mod.11 <- lm(CUM_CO_COUNT ~ 1 + scale(pct_gc) + scale(COPIA) + scale(GYPSY), data=gen.feat.df) # all significant
summary(lm.mod.11)
lm.mod.12 <- lm(CUM_CO_COUNT ~ 1 + scale(GENE) + scale(COPIA) + scale(GYPSY), data=gen.feat.df) # all significant
summary(lm.mod.12)
AIC(lm.mod.1, lm.mod.2, lm.mod.3, lm.mod.4, lm.mod.5, lm.mod.6, lm.mod.7, lm.mod.8, lm.mod.9, lm.mod.10, lm.mod.11, lm.mod.12)

# drop1(lm.mod.7)
# stepAIC(lm.mod.7)


# Model diagnositcs
op <- par(mfrow=c(2,2))
# plot(lm.mod.8)
plot(lm.mod.7)
par(op)


# cross-validation for lm.mod.8 using caret package
set.seed(123) 
train.control <- trainControl(method = "repeatedcv", number = 5)
# Train the model
# model.form <- formula(CUM_CO_COUNT ~ pct_gc + COPIA + GYPSY + SIMP_REPEAT)
model.form <- formula(CUM_CO_COUNT ~pct_gc + GENE + COPIA + GYPSY + SIMP_REPEAT)
model.lm <- train(model.form, data=gen.feat.df, method = "lm", trControl = train.control)
# Summarize the results
print(model.lm)



# investigate whether genomic features could explain differences between CO counts
subset.par.2.df <- all.par.df.2[all.par.df.2$CHROM %in% gen.feat.df$CHROM & all.par.df.2$START %in% gen.feat.df$START & all.par.df.2$END %in% gen.feat.df$END, ]
if (nrow(subset.par.2.df) != nrow(gen.feat.df)) {
  stop("ERROR in subsetting!")
} else {
  fem.CO.count <- apply(subset.par.2.df[4:10], 1, sum)
  male.CO.count <- apply(subset.par.2.df[11:17], 1, sum)
  count.diff <- male.CO.count - fem.CO.count
  gen.feat.df <- cbind(gen.feat.df, "fem.CO.count"=fem.CO.count, "male.CO.count"=male.CO.count, diff=(male.CO.count - fem.CO.count))
}
head(gen.feat.df)
hist(gen.feat.df$diff)
shapiro.test(gen.feat.df$diff)

# considering linear models for the differnce of CO counts between sexes
lm.mod.13 <- lm(diff ~ 1 + scale(pct_at), data=gen.feat.df) # not siginificant
summary(lm.mod.13)
lm.mod.14 <- lm(diff ~ 1 + scale(pct_gc), data=gen.feat.df) # significant
summary(lm.mod.14)
lm.mod.15 <- lm(diff ~ 1 + scale(GENE), data=gen.feat.df) # highly significant
summary(lm.mod.15)
lm.mod.16 <- lm(diff ~ 1 + scale(COPIA), data=gen.feat.df) # highly significant
summary(lm.mod.16)
lm.mod.17 <- lm(diff ~ 1 + scale(GYPSY), data=gen.feat.df) # highly significant
summary(lm.mod.17)
lm.mod.18 <- lm(diff ~ 1 + scale(SIMP_REPEAT), data=gen.feat.df) # highly significant
summary(lm.mod.18)
lm.mod.19 <- lm(diff ~ 1 + scale(pct_gc) + scale(GENE) + scale(COPIA) + scale(GYPSY) + scale(SIMP_REPEAT), data=gen.feat.df) # GENE and SIMP_REPEAT not significant
summary(lm.mod.19)
lm.mod.20 <- lm(diff ~ 1 + scale(pct_gc) + scale(GENE), data=gen.feat.df) # SIMP_REPEAT not that significant
summary(lm.mod.20)
step(lm.mod.19)
AIC(lm.mod.13, lm.mod.14, lm.mod.15, lm.mod.16, lm.mod.17, lm.mod.18, lm.mod.19, lm.mod.20)


##############################################################################################################################################################
# Here I am fitting a GLM model rather than LMs
# glm.mod.1 <- stats::glm(CUM_CO_COUNT ~ pct_at, data=gen.feat.df, family = poisson(link = "log")) 
# summary(glm.mod.1)
# glm.mod.2 <- stats::glm(CUM_CO_COUNT ~ pct_gc, data=gen.feat.df, family = poisson(link = "log")) 
# summary(glm.mod.2)
# glm.mod.3 <- stats::glm(CUM_CO_COUNT ~ pct_at + GENE, data=gen.feat.df, family = poisson(link = "log")) 
# summary(glm.mod.3)
# glm.mod.4 <- stats::glm(CUM_CO_COUNT ~ pct_at + COPIA, data=gen.feat.df, family = poisson(link = "log")) 
# summary(glm.mod.4)
# glm.mod.4 <- stats::glm(CUM_CO_COUNT ~ pct_at + GENE + COPIA, data=gen.feat.df, family = poisson(link = "log")) 
# summary(glm.mod.4)
# glm.mod.5 <- stats::glm(CUM_CO_COUNT ~ pct_at + GENE + COPIA + GYPSY, data=gen.feat.df, family = poisson(link = "log")) 
# summary(glm.mod.5)
# glm.mod.6 <- stats::glm(CUM_CO_COUNT ~ pct_at + GENE + COPIA + GYPSY + SIMP_REPEAT, data=gen.feat.df, family = poisson(link = "log")) 
# summary(glm.mod.6)
# glm.mod.7 <- stats::glm(CUM_CO_COUNT ~ pct_at + GENE + COPIA + GYPSY + SIMP_REPEAT + CHROM, data=gen.feat.df, family = poisson(link = "log")) 
# sum.glm.mod.7 <- summary(glm.mod.7)
# 

# glmer.poiss.1 <- glmer(CUM_CO_COUNT ~ scale(pct_at) + scale(GENE) + scale(COPIA) + scale(GYPSY) + scale(SIMP_REPEAT) + (1|CHROM), data = gen.feat.df, family=poisson(link = "log")) # testing the CHROM as a variable with LENGTH
# sum.glmer.poiss.1 <- summary(glmer.poiss.1, cor = FALSE)
# plot(glmer.poiss.1)
# 
# AIC(glm.mod.1, glm.mod.2, glm.mod.3, glm.mod.4, glm.mod.5, glm.mod.6, glm.mod.7, glmer.poiss.1)

# op <- par(mfrow=c(2,2))
# plot(glm.mod.7)
# par(op)

# # checking all LM and GLM models together
# AIC(lm.mod.1, lm.mod.2, lm.mod.3, lm.mod.4, lm.mod.5, lm.mod.6, lm.mod.7, glm.mod.1, glm.mod.2, glm.mod.3, glm.mod.4, glm.mod.5, glm.mod.6, glm.mod.7, glmer.poiss.1)
# 
# # This is the closest metric to a R^2 we can get using a GLM model
# 100 * (glm.mod.7$null.deviance - glm.mod.7$deviance) / glm.mod.7$null.deviance 
# 
# 
# # cross-validation for the GLMs
# set.seed(123) 
# train.control <- trainControl(method = "cv", number = 5)
# # Train the model
# model.glm <- train(CUM_CO_COUNT ~ pct_at + GENE + COPIA + GYPSY + SIMP_REPEAT, data=gen.feat.df, method = "glm", trControl = train.control)
# # Summarize the results
# print(model.glm)
########################################################################################################################################################################################################
# END
########################################################################################################################################################################################################


# # manually enter centromere positions
# cent.pos.vec <- c(75, 80, 35, 55, 55, 60, 30, 60, 0, 25, 40, 30, 35, 60, 30, 35, 30, 30, 30)
# dist.frm.cent <- NULL
# for (chrom.name in chrom.name.vec) {
#   chrom.seg.idx <- seq(from=1, to=sum(gen.feat.df$CHROM == chrom.name), by=1)
#   tmp.dist.cen <- abs(chrom.seg.idx - cent.pos.vec[which(chrom.name.vec %in% chrom.name)])
#   dist.frm.cent <- c(dist.frm.cent, tmp.dist.cen)
# }
# if (length(dist.frm.cent) == nrow(gen.feat.df)) {
#   gen.feat.df <- cbind(gen.feat.df, DIST_FRM_CENT=dist.frm.cent)
# } else {
#   stop()
# }
###########################################################################################################################################################


# Making a dataset with sex as a factor
fem_df <- gen.feat.df[, c('pct_gc', 'COPIA', 'GYPSY', 'SIMP_REPEAT', 'GENE', 'WIN_IDX', 'fem.CO.count')]
male_df <- gen.feat.df[, c('pct_gc', 'COPIA', 'GYPSY', 'SIMP_REPEAT', 'GENE', 'WIN_IDX', 'male.CO.count')]
colnames(fem_df) <- c('pct_gc', 'COPIA', 'GYPSY', 'SIMP_REPEAT', 'GENE', 'WIN_IDX', 'cum_CO_count')
colnames(male_df) <- c('pct_gc', 'COPIA', 'GYPSY', 'SIMP_REPEAT', 'GENE', 'WIN_IDX', 'cum_CO_count')
dataset <- rbind(fem_df, male_df)
dataset$sex <- rep(c('F', 'M'), each=319)
head(dataset)



# considering linear models with sexe as a factor
lm.mod.21 <- lm(cum_CO_count ~ 1 + scale(pct_gc) + scale(GENE) + scale(COPIA) + scale(GYPSY) + scale(SIMP_REPEAT) + sex, data=dataset)
summary(lm.mod.21)
lm.mod.22 <- lm(cum_CO_count ~ 1 + scale(pct_gc) + scale(GENE) + scale(COPIA) + scale(GYPSY) + sex, data=dataset)
summary(lm.mod.22)
lm.mod.23 <- lm(cum_CO_count ~ 1 + scale(pct_gc) + scale(COPIA) + scale(GYPSY) + sex, data=dataset)
summary(lm.mod.23)
lm.mod.24 <- lm(cum_CO_count ~ 1 + scale(pct_gc) + scale(GENE) + scale(COPIA) + scale(GYPSY) + scale(SIMP_REPEAT), data=dataset)
summary(lm.mod.24)

AIC(lm.mod.21, lm.mod.22, lm.mod.23, lm.mod.24)

anova(lm.mod.21, lm.mod.22)
anova(lm.mod.21, lm.mod.23)
anova(lm.mod.21, lm.mod.24)


# cross-validation for lm.mod.21 using caret package
set.seed(123) 
train.control <- trainControl(method = "repeatedcv", number = 5)
# Train the model
model.form <- formula(cum_CO_count ~ 1 + scale(pct_gc) + scale(GENE) + scale(COPIA) + scale(GYPSY) + scale(SIMP_REPEAT) + sex)
model.form <- formula(cum_CO_count ~ 1 + scale(pct_gc) + scale(GENE) + scale(COPIA) + scale(GYPSY) + sex)
model.lm <- train(model.form, data=dataset, method = "lm", trControl = train.control)
# Summarize the results
print(model.lm)



###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
########################################################################################################################################################################################################