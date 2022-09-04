modified_order_seq <- function (input.seq, n.init = 5, subset.search = c("twopt", "sample"), 
          subset.n.try = 30, subset.THRES = 3, twopt.alg = c("rec", 
                                                             "rcd", "ser", "ug"), THRES = 3, touchdown = FALSE, tol = 0.1) 
{
  if (!any(class(input.seq) == "sequence")) 
    stop(deparse(substitute(input.seq)), " is not an object of class 'sequence'")
  if (n.init < 2) 
    stop("'n.init' must be greater than or equal to 2")
  if (!is.logical(touchdown)) 
    stop("'touchdown' must be logical")
  if (!touchdown && THRES <= 1e-09) 
    stop("Threshold must be greater than 0 if 'touchdown' is FALSE")
  if (touchdown && THRES <= (1 + 1e-09)) 
    stop("Threshold must be greater than 1 if 'touchdown' is TRUE")
  if (length(input.seq$seq.num) <= n.init) {
    cat("   Length of sequence ", deparse(substitute(input.seq)), 
        " is less than n.init \n   Returning the best order using compare function:\n")
    ifelse(length(input.seq$seq.num) == 2, seq.ord <- map(input.seq, 
                                                          tol = 1e-04), seq.ord <- make_seq(compare(input.seq = input.seq, 
                                                                                                    tol = 1e-04), 1))
    seq.ord <- map(seq.ord, tol = 1e-04)
    structure(list(ord = seq.ord, mrk.unpos = NULL, LOD.unpos = NULL, 
                   THRES = THRES, ord.all = seq.ord, data.name = input.seq$data.name, 
                   twopt = input.seq$twopt), class = "order")
  }
  else {
    cross.type <- class(get(input.seq$data.name, pos = 1))[2]
    if (cross.type == "f2") 
      FLAG <- "f2"                                                      
    else if (cross.type == "backcross" || cross.type == "riself" || 
             cross.type == "risib") 
      FLAG <- "bc"                                                        #Run this line first then go to line 42
    else if (cross.type == "outcross") 
      FLAG <- "outcross"
    else stop("Invalid cross type\n")
    if (FLAG == "bc" || FLAG == "f2") {
      subset.search <- match.arg(subset.search)
      if (subset.search == "twopt") {
        cat("\nCross type: ", cross.type, "\nChoosing initial subset using 'two-point' approach\n")
        twopt.alg <- match.arg(twopt.alg)
        tpt.type <- switch(EXPR = twopt.alg, rec = {
          seq.rec <- record(input.seq = input.seq, tol = 0.1)$seq.num       # Run this line second then go to line 43
          seq.init <- seq.rec[unique(round(seq(from = 1, to = length(seq.rec), length.out = n.init)))] # Run this phrase third to see the marker names for initial order colnames(get(input.seq$data.name,pos = 1)$geno)[seq.init] # go to line 143
          for (iter in 1:10) { # iter <- 1
            seq.ord <- compare(input.seq = make_seq(get(input.seq$twopt),seq.init, twopt = input.seq$twopt), n.best = 50)
            input.seq2.temp <- make_seq(seq.ord, 1)   
            input.seq2.temp.idx <- input.seq2.temp$seq.num
            seq.init.temp.idx <- order(as.numeric(gsub(paste0(chrom.name,"_"), "", colnames(get(input.seq$data.name,pos = 1)$geno)[input.seq2.temp.idx])), decreasing = FALSE)
            seq.init.mod <- input.seq2.temp.idx[seq.init.temp.idx]
            if (all(input.seq2.temp$seq.num == seq.init.mod) | all(input.seq2.temp$seq.num == seq.init.mod[n.init:1])) {
              break
            } else {
              break.points <- unique(round(seq(from = 1, to = length(seq.rec), length.out = (n.init + 1))))
              seq.init <- NULL
              for (i in 1:n.init) { # i <- 1
                window.start <- break.points[i]
                if (i == n.init) {
                  window.end <- (break.points[i + 1])
                } else {
                  window.end <- (break.points[i + 1] - 1)
                }
                seq.init.iter <- sample(seq.rec[window.start:window.end], 1, replace = FALSE)
                seq.init <- c(seq.init, seq.init.iter)
              }
              # quarter.point <- floor(length(seq.rec) * 0.10)
              # seq.init <- sample(seq.rec[-c(1:quarter.point, length(seq.rec):(length(seq.rec)-quarter.point))], n.init, replace = FALSE)
              cat("iteration number...", " ", iter, "\n" )
            }
          }
        }, rcd = {
          seq.rcd <- rcd(input.seq = input.seq, tol = 0.1)$seq.num
          seq.init <- seq.rcd[unique(round(seq(from = 1, 
                                               to = length(seq.rcd), length.out = n.init)))]
        }, ser = {
          seq.ser <- seriation(input.seq = input.seq, 
                               tol = 0.1)$seq.num
          seq.init <- seq.ser[unique(round(seq(from = 1, 
                                               to = length(seq.ser), length.out = n.init)))]
        }, ug = {
          seq.ug <- ug(input.seq = input.seq, tol = 0.1)$seq.num
          seq.init <- seq.ug[unique(round(seq(from = 1, 
                                              to = length(seq.ug), length.out = n.init)))]
        })
        seq.rest <- input.seq$seq.num[-pmatch(seq.init, input.seq$seq.num)]          # Run this line fifth, go to line 62
        seq.mis <- apply(as.matrix(get(input.seq$data.name, pos = 1)$geno[, seq.rest]), 2, function(x) sum(x == 0)) # Run this line sixth, go to line 65
        names(seq.mis) <- colnames(get(input.seq$data.name, pos = 1)$geno)[seq.rest]                             #Run this line seventh, go to line 68. 
        if (FLAG == "bc") {
          rest.ord <- pmatch(names(seq.mis), colnames(get(input.seq$data.name, pos = 1)$geno))                   #Run this line eighth, go to line 70. 
          seq.work <- pmatch(c(seq.init, rest.ord), input.seq$seq.num)                      #Run this line ninth, go to line 149.
        }
        else if (FLAG == "f2") {
          seq.type <- get(input.seq$data.name, pos = 1)$segr.type.num[seq.rest]
          names(seq.type) <- names(seq.mis)
          tp.ord <- sort(seq.type)
          rest.ord <- c(sort(seq.mis[pmatch(names(tp.ord)[tp.ord == 
                                                            1], names(seq.mis))]), sort(seq.mis[pmatch(names(tp.ord)[tp.ord == 
                                                                                                                       4], names(seq.mis))]), sort(seq.mis[pmatch(names(tp.ord)[tp.ord == 
                                                                                                                                                                                  2], names(seq.mis))]), sort(seq.mis[pmatch(names(tp.ord)[tp.ord == 
                                                                                                                                                                                                                                             3], names(seq.mis))]))
          rest <- pmatch(names(rest.ord), colnames(get(input.seq$data.name, 
                                                       pos = 1)$geno))
          seq.work <- pmatch(c(seq.init, rest), input.seq$seq.num)
        }
        else stop("Invalid cross type\n")
      }
      else if (subset.search == "sample") {
        cat("\nCross type: ", cross.type, "\nChoosing initial subset using the 'sample' approach\n")
        LOD.test <- i <- 0
        while (abs(LOD.test) < abs(subset.THRES) && i < 
               subset.n.try) {
          smp.seq <- make_seq(get(input.seq$twopt), sample(input.seq$seq.num, 
                                                           size = n.init), twopt = input.seq$twopt)
          res.test <- compare(smp.seq)
          LOD.test <- res.test$best.ord.LOD[2]
          i < -i + 1
        }
        if (abs(LOD.test) >= abs(subset.THRES)) {
          seq.init <- res.test$best.ord[1, ]
          seq.rest <- input.seq$seq.num[-pmatch(seq.init, 
                                                input.seq$seq.num)]
          seq.mis <- apply(as.matrix(get(input.seq$data.name, 
                                         pos = 1)$geno[, seq.rest]), 2, function(x) sum(x == 
                                                                                          0))
          names(seq.mis) <- colnames(get(input.seq$data.name, 
                                         pos = 1)$geno)[seq.rest]
          if (FLAG == "bc") {
            rest.ord <- pmatch(names(seq.mis), colnames(get(input.seq$data.name, 
                                                            pos = 1)$geno))
            seq.work <- pmatch(c(seq.init, rest.ord), 
                               input.seq$seq.num)
          }
          else if (FLAG == "f2") {
            seq.type <- get(input.seq$data.name, pos = 1)$segr.type.num[seq.rest]
            names(seq.type) <- names(seq.mis)
            tp.ord <- sort(seq.type)
            rest.ord <- c(sort(seq.mis[pmatch(names(tp.ord)[tp.ord == 
                                                              1], names(seq.mis))]), sort(seq.mis[pmatch(names(tp.ord)[tp.ord == 
                                                                                                                         4], names(seq.mis))]), sort(seq.mis[pmatch(names(tp.ord)[tp.ord == 
                                                                                                                                                                                    2], names(seq.mis))]), sort(seq.mis[pmatch(names(tp.ord)[tp.ord == 
                                                                                                                                                                                                                                               3], names(seq.mis))]))
            rest <- pmatch(names(rest.ord), colnames(get(input.seq$data.name, 
                                                         pos = 1)$geno))
            seq.work <- pmatch(c(seq.init, rest), input.seq$seq.num)
          }
          else stop("Invalid cross type\n")
        }
        else stop("Cannot find any subset using 'subset.n.try'=", 
                  subset.n.try, " and 'subset.THRES'= ", subset.THRES, 
                  "\n")
      }
    }
    else if (FLAG == "outcross") {
      cat("\nCross type: outcross\nUsing segregation types of the markers to choose initial subset\n")
      segregation.types <- get(input.seq$data.name, pos = 1)$segr.type.num[input.seq$seq.num]
      if (sum(segregation.types == 7) > sum(segregation.types == 
                                            6)) 
        segregation.types[segregation.types == 6] <- 8
      seq.work <- order(segregation.types)
      seq.init <- input.seq$seq.num[seq.work[1:n.init]]
    }
    else stop("Invalid cross type")
    seq.ord <- compare(input.seq = make_seq(get(input.seq$twopt), seq.init, twopt = input.seq$twopt), n.best = 50)  #Run this line fourth to evaluate the initial marker ordering. 
                                                                                              #If the markers are not ordered right, then use seq.init[i]<-  to assigne new markers to order
                                                                                              #use seq.rec to see the marker order based on recombination fraction
                                                                                              #Then run line 143 again and iterate until right order obtained. check order with phrase in line 44
                                                                                              #Once done go back to line 60
                                                                                              
    input.seq2 <- make_seq(seq.ord, 1)                                                        #Run this line tenth all the way to line 205. But run it one line at a time
    cat("\n\nRunning try algorithm\n")
    for (i in (n.init + 1):length(input.seq$seq.num)) { # i <- 8
      seq.ord <- try_seq(input.seq2, input.seq$seq.num[seq.work[i]], tol = tol) 
      if (all(seq.ord$LOD[-which(seq.ord$LOD == max(seq.ord$LOD))[1]] < -THRES)) 
        input.seq2 <- make_seq(seq.ord, which.max(seq.ord$LOD)) 
    }
    mrk.unpos <- input.seq$seq.num[which(is.na(match(input.seq$seq.num, input.seq2$seq.num)))]
    LOD.unpos <- NULL
    cat("\nLOD threshold =", THRES, "\n\nPositioned markers:", input.seq2$seq.num, "\n\n")
    cat("Markers not placed on the map:", mrk.unpos, "\n")
    if (touchdown && length(mrk.unpos) > 0) {
      cat("\n\n\nTrying to map remaining markers with LOD threshold ", THRES - 1, "\n")
      for (i in mrk.unpos) { # i <- 2
        seq.ord <- try_seq(input.seq2, i, tol = tol)
        if (all(seq.ord$LOD[-which(seq.ord$LOD == max(seq.ord$LOD))[1]] < 
                (-THRES + 1))) 
          input.seq2 <- make_seq(seq.ord, which.max(seq.ord$LOD))
      }
      mrk.unpos <- input.seq$seq.num[which(is.na(match(input.seq$seq.num, input.seq2$seq.num)))]
      cat("\nLOD threshold =", THRES - 1, "\n\nPositioned markers:", input.seq2$seq.num, "\n\n")
      cat("Markers not placed on the map:", mrk.unpos, "\n")
    }
    if (length(mrk.unpos) > 0) {
      LOD.unpos <- matrix(NA, length(mrk.unpos), (length(input.seq2$seq.num) + 1))
      j <- 1
      cat("\n\nCalculating LOD-Scores\n")
      for (i in mrk.unpos) { # i <- 3
        LOD.unpos[j, ] <- try_seq(input.seq = input.seq2, mrk = i, tol = tol)$LOD
        j <- j + 1
      }
    }
    else mrk.unpos <- NULL
    input.seq3 <- input.seq2
    if (!is.null(mrk.unpos)) {
      cat("\n\nPlacing remaining marker(s) at most likely position\n")
      which.order <- order(apply(LOD.unpos, 1, function(x) max(x[-which(x == 0)[1]])))
      for (i in mrk.unpos[which.order]) { # i <- 68
        seq.ord <- try_seq(input.seq3, i, tol)
        input.seq3 <- make_seq(seq.ord, which(seq.ord$LOD == 0)[sample(sum(seq.ord$LOD == 0))[1]])
      }
    }
    cat("\nEstimating final genetic map using tol = 10E-5.\n\n")
    input.seq2 <- map(input.seq2, tol = 1e-04)
    input.seq3 <- map(input.seq3, tol = 1e-04)                                            #This is the last line of code to run. Finally run LG1_rec<- input.seq2
    structure(list(ord = input.seq2, mrk.unpos = mrk.unpos, 
                   LOD.unpos = LOD.unpos, THRES = THRES, ord.all = input.seq3, 
                   data.name = input.seq$data.name, twopt = input.seq$twopt), 
              class = "order")
  }
}