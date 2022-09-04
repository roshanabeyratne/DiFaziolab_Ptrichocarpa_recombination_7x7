library(WaveletComp)

# Example 1
x = periodic.series(start.period = 50, length = 100, make.plot = TRUE)
class(x)
x = x + 0.2*rnorm(1000)  # add some noise
plot(x, xlim = c(0,200))

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7))
reconstruct(my.w, plot.waves = FALSE, lwd = c(1,2),
            legend.coords = "bottomleft", ylim = c(-1.8, 1.8))


# Example 2
x = periodic.series(start.period = 20, end.period = 100, length = 1000, make.plot = TRUE) # linearly increasing trend
x = x + 0.2*rnorm(1000)
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))


# Example 3
x1 <- periodic.series(start.period = 80, length = 1000)
x2 <- periodic.series(start.period = 30, length = 1000)
x <- x1 + x2 + 0.2*rnorm(1000)
plot((x))
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels") )

reconstruct(my.w, sel.period = 80, plot.waves = TRUE, lwd = c(1,2),
            legend.coords = "bottomleft")


# Example 4
x1 <- periodic.series(start.period = 100, length = 500)
x2 <- 1.2*periodic.series(start.period = 60, length = 500)
x <- c(x1, x2) + 0.3*rnorm(1000)
plot((x))

y1 <- periodic.series(start.period = 100, length = 1000)
y2 <- 1.2*periodic.series(start.period = 60, length = 1000)
y <- (y1 + y2)/2  + 0.3*rnorm(1000)
plot((y))


my.data <- data.frame(x = x, y = y)
my.wx <- analyze.wavelet(my.data, "x", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 16, upperPeriod = 256,
                         make.pval = TRUE, n.sim = 10)
my.wy <- analyze.wavelet(my.data, "y", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 16, upperPeriod = 256,
                         make.pval = TRUE, n.sim = 10)

maximum.level = 1.001*max(my.wx$Power.avg, my.wy$Power.avg)
wt.avg(my.wx, maximum.level = maximum.level)
wt.avg(my.wy, maximum.level = maximum.level)


# Example 5
x1 <- periodic.series(start.period = 100, length = 400, make.plot = T)
x2 <- 1.2*periodic.series(start.period = 50, length = 200, make.plot = T)
x  <- c(x1, x2, x1) + 0.2*rnorm(1000)
plot(x)

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        method = "white.noise",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 32, upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = "interval", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))
wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2))
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         useRaster = FALSE, max.contour.segments = 500)



# column<- 8
# 
tmp.1 <- apply(all.par.df.2[all.par.df.2$CHROM == "Chr01",4:10], 1, sum)
tmp.1 <- (tmp.1 - mean(tmp.1))/sqrt(var(tmp.1))
tmp.2 <- apply(all.par.df.2[all.par.df.2$CHROM == "Chr01",11:17], 1, sum)
tmp.2 <- (tmp.2 - mean(tmp.2))/sqrt(var(tmp.2))


tmp.1 <- ((all.par.df.2[all.par.df.2$CHROM == "Chr01",6] - mean(all.par.df.2[all.par.df.2$CHROM == "Chr01",6]))/ sqrt(var(all.par.df.2[all.par.df.2$CHROM == "Chr01",6])))
tmp.2 <- ((all.par.df.2[all.par.df.2$CHROM == "Chr01",7] - mean(all.par.df.2[all.par.df.2$CHROM == "Chr01",7]))/ sqrt(var(all.par.df.2[all.par.df.2$CHROM == "Chr01",7])))

my.data <- data.frame(x = tmp.1, y = tmp.2)

# my.data <- data.frame(x = all.par.df.2[all.par.df.2$CHROM == "Chr01",column] - mean(all.par.df.2[all.par.df.2$CHROM == "Chr01",column]))
my.w <- analyze.wavelet(my.data, "x",
                        method = "white.noise",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 2, upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.wy, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7))
maximum.level = 1.001*max(my.w$Power.avg, my.wy$Power.avg)
wt.avg(my.w, maximum.level = maximum.level)
  
par(mfrow=c(1,1))




my.wx <- analyze.wavelet(my.data, "x", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 256,
                         make.pval = TRUE, n.sim = 10)
my.wy <- analyze.wavelet(my.data, "y", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 256,
                         make.pval = TRUE, n.sim = 10)
wt.image(my.wx, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.image(my.wy, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7))

maximum.level = 1.001*max(my.wx$Power.avg, my.wy$Power.avg)
wt.avg(my.wx, maximum.level = maximum.level)
wt.avg(my.wy, maximum.level = maximum.level)














