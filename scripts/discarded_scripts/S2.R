test.1 <- rnorm(1000)
test.2 <- test.1 + 4
test.1 <- rpois(1000, 1)
test.2 <- test.1 + 10
# test.2 <- rpois(1000, 1)
# test.2 <- c(test.2[1:300]+3, test.2[301:1000])
# plot(test.1, col='blue', type='l')
# plot(test.2, col='red', type='l')
ks.test(test.1, test.2)


wavMODWT.obj.test.1 <- wavMODWT(test.1, wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))), 
                                    position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())
wavMODWT.obj.test.2 <- wavMODWT(test.2, wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))), 
                                position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())
wav.coffs.1 <- as.matrix(wavMODWT.obj.test.1)
wav.coffs.2 <- as.matrix(wavMODWT.obj.test.2)



for (i in 1:num.levels) {
  start.corrd <- (i*num.coffs)-(num.coffs-1)
  end.coord <- (i*num.coffs)
  
  # tmp.sum.sq.fem <- sum(permuted.fem.coffs.df[start.corrd:end.coord, iter]^2)
  # tmp.sum.sq.male <- sum(permuted.male.coffs.df[start.corrd:end.coord, iter]^2)
  # # tmp.sum.sq.diff <- sum((permuted.fem.coffs.df[start.corrd:end.coord, iter] - permuted.male.coffs.df[start.corrd:end.coord, iter])^2)
  # tmp.sum.sq.diff <- sum(permuted.sex.diff.coffs.df[start.corrd:end.coord, iter]^2)
  # 
  # scale.sum.sq.fem <- c(scale.sum.sq.fem, tmp.sum.sq.fem)
  # scale.sum.sq.male <- c(scale.sum.sq.male, tmp.sum.sq.male)
  # scale.sum.sq.diff <- c(scale.sum.sq.diff, tmp.sum.sq.diff)
  
  ks.test.mod <- ks.test(wav.coffs.1[start.corrd:end.coord],wav.coffs.2[start.corrd:end.coord])
  # tmp.ks.vec <- c(tmp.ks.vec, ks.test.mod$statistic) 
  cat(ks.test.mod$p.value, "\n")
  
}




tmp.obj.1 <- wavVar(test.1, xform="modwt", wavelet="haar", n.levels=NULL,
                    position=list(from=1,by=1,units=character()), units=character(),
                    documentation=character(), sdf=NULL, sdfargs=NULL,
                    sampling.interval=1, n.fft=1024)
tmp.obj.2 <- wavVar(test.2, xform="modwt", wavelet="haar", n.levels=NULL,
                    position=list(from=1,by=1,units=character()), units=character(),
                    documentation=character(), sdf=NULL, sdfargs=NULL,
                    sampling.interval=1, n.fft=1024)



wavMODWT.obj.fem <- wavMODWT(cos.vec.fem-cos.vec.male, wavelet="haar", n.levels=as.integer(floor(logb(length(cos.vec.fem),base = 2))), 
                             position = list(from=1, by=1, units=character()), units=character(), title.data=character(), documentation=character())



# permuted.ks.dist.df <- NULL
# for (iter in 1:1000) {
#   tmp.ks.vec <- NULL
#   for (i in 1:num.levels) {
#     start.corrd <- (i*num.coffs)-(num.coffs-1)
#     end.coord <- (i*num.coffs)
#     
#     ks.test.mod <- ks.test(permuted.fem.coffs.df[start.corrd:end.coord, iter], permuted.male.coffs.df[start.corrd:end.coord, iter])
#     chisq.test.mod <- chisq.test(permuted.fem.coffs.df[start.corrd:end.coord, iter], permuted.male.coffs.df[start.corrd:end.coord, iter])
#     
#     tmp.ks.vec <- c(tmp.ks.vec, ks.test.mod$statistic) 
#     
#   }
#   permuted.ks.dist.df <- rbind(permuted.ks.dist.df, tmp.ks.vec)
#   cat(iter, " ")
# }
# 

# constant off-set of the same wave does not affect DWT wariance decomposition. Also if the individuals irrespective of their SEX do not show marked variation in CO-pattern along the chromosome, 
# variance partition analysis does not identify that. scale based wavelet coefficient comparisons between sexes can be used to identify these differences. 
# permuted.ks.dist.df <- NULL
# for (iter in 1:1000) {
#   tmp.ks.vec <- NULL
#   for (i in 1:num.levels) {
#     start.corrd <- (i*num.coffs)-(num.coffs-1)
#     end.coord <- (i*num.coffs)
#     
#     ks.test.mod <- ks.test(permuted.fem.coffs.df[start.corrd:end.coord, iter], permuted.male.coffs.df[start.corrd:end.coord, iter])
#     chisq.test.mod <- chisq.test(permuted.fem.coffs.df[start.corrd:end.coord, iter], permuted.male.coffs.df[start.corrd:end.coord, iter])
#     
#     tmp.ks.vec <- c(tmp.ks.vec, ks.test.mod$statistic) 
#     
#   }
#   permuted.ks.dist.df <- rbind(permuted.ks.dist.df, tmp.ks.vec)
#   cat(iter, " ")
# }
# 
# # actual KS-test 
# actual.ks.dist.vec <- NULL
# for (i in 1:num.levels) {
#   start.corrd <- (i*num.coffs)-(num.coffs-1)
#   end.coord <- (i*num.coffs)
#   
#   # print(length(unique(fem.coffs[start.corrd:end.coord])))
#   
#   ks.test.mod <- ks.test(fem.coffs[start.corrd:end.coord], male.coffs[start.corrd:end.coord])
#   actual.ks.dist.vec <- c(actual.ks.dist.vec, ks.test.mod$statistic)
#   chisq.test.mod <- chisq.test(fem.coffs[start.corrd:end.coord], male.coffs[start.corrd:end.coord])
# }
# 
# # plot the NULL KS distance vs. observed for each scale
# jpeg(filename = paste0(statistics.output.path, chrom.name, "_combined_DWT_scale_KS_test_stat_distribution.jpg"), width = 2000, height = 1000, quality = 150)
# par(mfrow=c(4,3))
# for (d in 1:ncol(permuted.ks.dist.df)) {
#   ks.p.val <- sum(permuted.ks.dist.df[,d] > actual.ks.dist.vec[d]) / nrow(permuted.ks.dist.df)
#   # male.p.val <- sum(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d])> wv.var.part.df[2,d]) / length(c(male.wv.var.part.permut[,d], fem.wv.var.part.permut[,d]))
#   
#   hist(permuted.ks.dist.df[,d], breaks=100, 
#        main=paste0(chrom.name, " d", d, "; KS_distance p-value=", ks.p.val), 
#        xlab="KS test statistic")
#   abline(v=actual.ks.dist.vec[d], col='black', lwd=3)
# }
# dev.off()


# # This is an alternate way to estimate variance partition at different scales using DWT coefficients when sex permuted. This is not preferred as EDOFs are not accurate in this manual derivation of
# DWT variance at each scale.
# permuted.power.df <- NULL
# permuted.power.diff.df <- NULL
# permuted.ks.D.val.df <- NULL
# 
# for (iter in 1:1000) {
#   
#   scale.sum.sq.fem <- NULL
#   scale.sum.sq.male <- NULL
#   scale.sum.sq.diff <- NULL
#   
#   tmp.ks.vec <- NULL
#   for (i in 1:num.levels) {
#     start.corrd <- (i*num.coffs)-(num.coffs-1)
#     end.coord <- (i*num.coffs)
#     
#     tmp.sum.sq.fem <- sum(permuted.fem.coffs.df[start.corrd:end.coord, iter]^2)
#     tmp.sum.sq.male <- sum(permuted.male.coffs.df[start.corrd:end.coord, iter]^2)
#     # tmp.sum.sq.diff <- sum((permuted.fem.coffs.df[start.corrd:end.coord, iter] - permuted.male.coffs.df[start.corrd:end.coord, iter])^2)
#     tmp.sum.sq.diff <- sum(permuted.sex.diff.coffs.df[start.corrd:end.coord, iter]^2)
# 
#     scale.sum.sq.fem <- c(scale.sum.sq.fem, tmp.sum.sq.fem)
#     scale.sum.sq.male <- c(scale.sum.sq.male, tmp.sum.sq.male)
#     scale.sum.sq.diff <- c(scale.sum.sq.diff, tmp.sum.sq.diff)
#     
#   }
#   permuted.ks.D.val.df <- rbind(permuted.ks.D.val.df, tmp.ks.vec)
#   
#   var.part.vec.fem <- scale.sum.sq.fem /sum(scale.sum.sq.fem)
#   var.part.vec.male <- scale.sum.sq.male /sum(scale.sum.sq.male)
#   var.part.vec.diff <- scale.sum.sq.diff /sum(scale.sum.sq.male)
# 
#   permuted.power.df <- rbind(permuted.power.df, var.part.vec.fem)
#   permuted.power.df <- rbind(permuted.power.df, var.part.vec.male)
#   permuted.power.diff.df <- rbind(permuted.power.diff.df, var.part.vec.diff)
#   
# }
# colnames(permuted.power.df) <- paste0("d_", 1:num.levels)
# colnames(permuted.power.diff.df) <- paste0("d_", 1:num.levels)
# # boxplot((permuted.power.df))
# 
# reshaped.df <- NULL
# for (scale in 1:num.levels) {
#   tmp <- as.numeric(permuted.power.df[, scale])
#   tmp.2 <- rep(paste0("d_", scale), length(tmp))
#   reshaped.df <- rbind(reshaped.df, cbind(tmp.2, tmp))
# }
# reshaped.df <- as.data.frame(reshaped.df, stringsAsFactors=FALSE)
# colnames(reshaped.df) <- c("SCALE", "VAR_PART")
# reshaped.df$VAR_PART <- as.numeric(reshaped.df$VAR_PART)
# 
# # plot(reshaped.df$VAR_PART ~ jitter(as.numeric(as.factor(reshaped.df$SCALE)), 1), pch=16, cex=0.1, xlab="", xaxt='n', ylab="partitioning of variance")
# # plot(reshaped.df$VAR_PART ~ as.factor(reshaped.df$SCALE), pch=16, cex=0.1, add=T, xlab=NULL, ylab=NULL)
# 
# # calculate actual variance partitioning estimates when the SEX assignments are accurate/actual.
# scale.sum.sq.fem <- NULL
# scale.sum.sq.male <- NULL
# scale.sum.sq.diff <- NULL
# 
# for (i in 1:num.levels) {
#   start.corrd <- (i*num.coffs)-(num.coffs-1)
#   end.coord <- (i*num.coffs)
#   
#   # print(length(unique(fem.coffs[start.corrd:end.coord])))
#   
#   fem.tmp.sum.sq <- sum(fem.coffs[start.corrd:end.coord]^2)
#   male.tmp.sum.sq <- sum(male.coffs[start.corrd:end.coord]^2)
#   # diff.tmp.sum.sq <- sum((fem.coffs[start.corrd:end.coord] - male.coffs[start.corrd:end.coord])^2)
#   diff.tmp.sum.sq <- sum(sex.diff.coffs[start.corrd:end.coord]^2)
#   
#   ks.test.mod <- ks.test(fem.coffs[start.corrd:end.coord], male.coffs[start.corrd:end.coord])
#   cat(ks.test.mod$statistic, "\n")
# 
#   
#   # scale.sum.sq.fem <- c(scale.sum.sq.fem, fem.tmp.sum.sq)
#   # scale.sum.sq.male <- c(scale.sum.sq.male, male.tmp.sum.sq)
#   # scale.sum.sq.diff <- c(scale.sum.sq.diff, diff.tmp.sum.sq)
#   
# }
# var.part.fem <- scale.sum.sq.fem /sum(scale.sum.sq.fem)
# var.part.male <- scale.sum.sq.male /sum(scale.sum.sq.male)
# var.part.diff <- scale.sum.sq.diff /sum(scale.sum.sq.diff)
# 
# # calculate the p-value for scale variance partition observations
# fem.p.vals <- NULL
# male.p.vals <- NULL
# diff.pvals <- NULL
# 
# for (scale in 1:num.levels) {
#   fem.p.vals <- c(fem.p.vals, (sum(as.vector(permuted.power.df[, scale]) > var.part.fem[scale]) / (length(permuted.power.df[, scale]) + 1)))
#   male.p.vals <- c(male.p.vals, (sum(as.vector(permuted.power.df[, scale]) > var.part.male[scale]) / (length(permuted.power.df[, scale]) + 1))) 
#   diff.pvals <- c(diff.pvals, (sum(as.vector(permuted.power.diff.df[, scale]) > var.part.diff[scale]) / (length(permuted.power.diff.df[, scale]) + 1))) 
#   
# }
# var.part.df <- cbind(SCALE=paste0("d_", 1:num.levels), FEM.P.VAL=round(fem.p.vals, 3), MALE.P.VAL=round(male.p.vals, 3), DIFF.P.VAL=round(diff.pvals, 3))
# write.csv(var.part.df, file=paste0(statistics.output.path, chrom.name, "_variance_partition_table.csv"), quote=FALSE, row.names=FALSE  )
# 
# # visualize the effect of variance partitioning by permuting sex. The rationale is that if at a given scale the sex effect has a major role to play, then the 
# jpeg(filename = paste0(statistics.output.path, chrom.name, "_combined_DWT_scale_variance_partitioning.jpg"), width = 2000, height = 1000, quality = 150)
# par(mfrow=c(2,6))
# for (scale in 1:num.levels) {
#   
#   # plot(density(permuted.power.df[,scale]), main=colnames(permuted.power.df)[scale])
#   # hist(permuted.power.df[,scale], breaks=100, main=colnames(permuted.power.df)[scale], xlab="partitioned variance at the given scale")
#   hist(permuted.power.diff.df[,scale], breaks=100, main=colnames(permuted.power.df)[scale], xlab="partitioned variance at the given scale")
#   
#   
#   # abline(v=var.part.fem[scale], col='red', lwd=3)
#   # abline(v=var.part.male[scale], col='blue', lwd=3)
#   abline(v=var.part.diff[scale], col='blue', lwd=3)
#   
#   
# }
# # par(mfrow=c(1,1))
# dev.off()




# 15d
# The difference between SEX averaged distributions could be due to chance. Checking whether the differences are due to chance using a KS test and comparing it to a NULL KS statistic
# derived from the permutation groups. 
# permuted.ks.dist.vec <- NULL
# for (iter in 1:1000) {
#   ks.test.mod <- ks.test(binned.permut.cos.fem.df[, iter], binned.permut.cos.male.df[, iter])
#   permuted.ks.dist.vec <- c(permuted.ks.dist.vec, ks.test.mod$statistic) 
#   cat(iter, " ")
# }
# actual.ks.dist <- ks.test(binned.cos.fem, binned.cos.male)$statistic
# ks.p.val <- sum(permuted.ks.dist.vec > actual.ks.dist) / length(permuted.ks.dist.vec)
# 
# hist(permuted.ks.dist.vec, breaks=100, main=paste0(chrom.name, " d",scale, " ;ks_stat p-value=", ks.p.val), xlab="ks test statistic")
# abline(v=actual.ks.dist, col='black', lwd=3)