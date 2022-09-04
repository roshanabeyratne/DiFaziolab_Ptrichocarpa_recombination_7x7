



# The dataframe containing the respective chromosome sizes. 
chrom.size.df <- read.csv(file = paste0(statistics.input.path, "chrom.size.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# modifying coordinates of Chr19 for Stettler
chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr11"] <- 19000000
chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr16"] <- 15494368
chrom.size.df$LENGTH[chrom.size.df$CHROM == "Chr19"] <- 16500000

# The dataframe with detailed cross-over events produced as per script 22b_...
CO.df <- read.csv(paste0(statistics.input.path, "co.detailed.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# sum(is.na(CO.df$CO_SIZE))
# sum(CO.df$CO_SIZE == 0)
# sum(CO.df$CO_SIZE != 0)
# sum(CO.df$CO_SIZE == 0) + sum(CO.df$CO_SIZE != 0) == nrow(CO.df)

# check whether same set of offspring are present for both male and female groups
# all(unique(CO.df$IND[CO.df$CHROM=="Chr01" & CO.df$PAR %in% par.name.vec[1:7]]) %in% unique(CO.df$IND[CO.df$CHROM=="Chr01" & CO.df$PAR %in% par.name.vec[8:14]])) 



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
  if (start.end.sites[1,1] > 1) {
    start.end.sites <- rbind(c(1, (start.end.sites[1,1] - 1)), start.end.sites)
  } 
  if (start.end.sites[nrow(start.end.sites),2] < chrom.length) {
    start.end.sites <- rbind(start.end.sites, c((start.end.sites[nrow(start.end.sites),2] + 1), chrom.length))
  }
  return(start.end.sites)
}




# Assign window start and end points based on window size
# win.size <- 250000 # This is the final window size statistical analysis done on.
# win.size <- 10000 # This is the final window size statistical analysis done on.
win.size <- 1920000 # This is the final window size statistical analysis done on.

total.bin.span.coord.df <- NULL
for (chrom.name in chrom.name.vec) { # chrom.name <- "Chr01"
  bin.span.coord.df.temp <- CreateBinSpanCoordinates(chrom.name = chrom.name, bin.size = win.size, chrom.length = chrom.size.df[chrom.size.df$CHROM == chrom.name, "LENGTH"], shift.size = win.size/(2^6))
  bin.span.coord.df.temp2 <- cbind(rep(chrom.name, nrow(bin.span.coord.df.temp)), bin.span.coord.df.temp)
  total.bin.span.coord.df <- rbind(total.bin.span.coord.df, bin.span.coord.df.temp2)
}

total.bin.span.coord.df <- as.data.frame(total.bin.span.coord.df, stringsAsFactors = FALSE)
total.bin.span.coord.df[,2:3] <- apply(total.bin.span.coord.df[,2:3], 2, as.numeric)
colnames(total.bin.span.coord.df) <- c("CHROM", "START", "END")



# Assign the COs to windows based on the the CO-region overlap. Creating a dataframe with CO counts in each overlapping windows.
# Each CO was assigned to a single window. This window was selected as the midpoint of the CO-region.
all.par.df <- NULL
for (par.name in par.name.vec) { # par.name <- par.name.vec[9]
  cat(date(), "\n")
  all.chroms.vec <- NULL
  for (chrom.name in chrom.name.vec) { # chrom.name <- chrom.name.vec[18]
    
    tmp.co.tab <- CO.df[(CO.df$CHROM == chrom.name & CO.df$PAR == par.name), ]
    
    # removing individuals that do not have an observed CO
    tmp.co.tab <- tmp.co.tab[tmp.co.tab$CO_SIZE != 0, ]
    
    # subset the windows for the chromosome 
    tmp.coord.tab <- total.bin.span.coord.df[total.bin.span.coord.df$CHROM == chrom.name, ]
    
    # tmp.vec <- rep(0, nrow(total.bin.span.coord.df))
    
    sum.vec <- rep(0, nrow(tmp.coord.tab))
    for (i in 1:nrow(tmp.co.tab)) { # i <- 6
      bp.coords <- c(tmp.co.tab[i,"DWNSTRM_FLNK"], tmp.co.tab[i,"UPSTRM_FLNK"])
      co.region <- seq(from = (bp.coords[1] + 1), to = (bp.coords[2] - 1), by = 1)
      
      # selecting the mid-point of the CO
      if ((length(co.region) %% 2) == 0) {
        cand.idx <- c( (length(co.region) %/% 2), ((length(co.region) %/% 2) + 1) )
        random.idx <- sample(cand.idx, 1) # randomly selecting an index out of the two possible.
        mid.point <- co.region[random.idx]
      } else {
        mid.point <- co.region[ceiling((length(co.region) / 2))]
      }
      
      # Following code block is if CO region is assigned to more than one window without considering the midpoint
      # tmp.vec <- apply(tmp.coord.tab[2:3], 1, function(x, y = bp.coords){
      #   test.1 <- ifelse(((x[1] > y[1] & x[1] < y[2]) | (x[2] >= y[1] & x[2] < y[2])), 1, 0)
      #   test.2 <- ifelse((y[1] > x[1] & y[2] < x[2]) & (y[2] > x[1] & y[2] < x[2]) , 1, 0)
      #   test <- test.1 + test.2
      #   return(test)
      # })
      
      # Considering the CO-midpoint assign the CO to a window
      tmp.vec <- apply(tmp.coord.tab[2:3], 1, function(x, y = mid.point){
        x <- as.numeric(x)
        test <- ifelse(((y >= x[1]) & (y <= x[2])), 1, 0)
        return(test)
      })
      
      # In this code since the windows are overlapping, CO event can be assigned to more than one window.
      # if (sum(tmp.vec) == 1) {
      #   sum.vec <- sum.vec + tmp.vec
      # } else {
      #   stop("ERROR with tmp.vec!")
      # }
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


write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_10kb.csv"), quote = FALSE, row.names = FALSE)
write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_30kb_2.csv"), quote = FALSE, row.names = FALSE)
write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_250kb.csv"), quote = FALSE, row.names = FALSE)
write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_2000kb.csv"), quote = FALSE, row.names = FALSE)
write.csv(all.par.df.2, file = paste0(statistics.input.path, "co.bin_size_1920kb.csv"), quote = FALSE, row.names = FALSE)

all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_10kb.csv"), header = TRUE, stringsAsFactors = FALSE )
all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_30kb_2.csv"), header = TRUE, stringsAsFactors = FALSE )
all.par.df.2 <- read.csv(file = paste0(statistics.input.path, "co.bin_size_250kb.csv"), header = TRUE, stringsAsFactors = FALSE )

colnames(all.par.df.2) <- colnames(all.par.df.2) %>% str_replace("^X", "")
# colnames(all.par.df.2)[4:17] <- par.name.vec



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





# just checking whether the expected pattern is observed for cross-over counts over each chromosome for each parent
library(RColorBrewer)
# colors <- c(brewer.pal(7, "Set3"), brewer.pal(7, "Set2"))
colors <- c('red', 'blue')
par(mfrow = c(1,2))
for (chrom.name in chrom.name.vec[1:2]) { # chrom.name <- chrom.name.vec[1]
  # jpeg(file = paste0(statistics.output.path, chrom.name, "_observed_p_value_on_Poisson_NULL.jpg"), width = 16, height = 9, unit = "in", res = 600)
  # pdf(file = paste0(statistics.output.path, chrom.name, "_observed_p_value_on_Poisson_NULL.pdf"), width = 16, height = 9, onefile = TRUE)
  # par(mfrow = c(2, 7))
  # colors <- c(brewer.pal(7, "Set3"), brewer.pal(7, "Set2"))
  
  
  # par.name.vec[c(1:12, 14)]
  for (par.name in par.name.vec[1:14]) { # par.name <- par.name.vec[1]
    if (par.name %in% par.name.vec[1:7]) {
      colr = 'red'
    } else {
      colr = 'blue'
    }
    if (par.name == par.name.vec[1]) {
      plot((all.par.df.2[all.par.df.2$CHROM == chrom.name, par.name] / num.ind.df[which(chrom.name == chrom.name.vec), par.name]), typ = "l", col = colr, main = chrom.name, xlab = "sliding_window_(not the win idx)",
           ylab = "number of cross-overs", ylim=c(0,0.5))
    } else {
      lines((all.par.df.2[all.par.df.2$CHROM == chrom.name, par.name] / num.ind.df[which(chrom.name == chrom.name.vec), par.name]), col = colr,
            # main = par.name, xlab = "window.idx",
            ylab = "number of cross-overs")
    }
  }
  
  # The mean of individual sexes
  female.vec <- apply((all.par.df.2[all.par.df.2$CHROM == chrom.name, par.name.vec[1:7]]), 1, sum) / 829
  male.ve <- apply((all.par.df.2[all.par.df.2$CHROM == chrom.name, par.name.vec[8:14]]), 1, sum) / 829 
  
  lines(female.vec, typ = "l", col = 'red', 
        # main = par.name, xlab = "window.idx", ylab = "number of cross-overs", ylim=c(0,0.5), 
        lwd=3)
  lines(male.ve, typ = "l", col = 'blue', lwd=3)
  
  
  
}
dev.off()
par(mfrow = c(1,1))
########################################################################################################################################################################################################
# END
###