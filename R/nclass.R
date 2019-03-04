# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2016
# Last update: March 2019
# Version 2.1
# Licence GPL v3 


if (!isGeneric("nclass")) {
  setGeneric("nclass", function(x,th)
    standardGeneric("nclass"))
}

# select the optimum number of classes based on one-standard error rule

setMethod('nclass', signature(x='RasterLayer'), 
          function(x,th=0.005) {
            if (missing(th)) th <- 0.005
            if (canProcessInMemory(x,2)) {
              w <- which(!is.na(x[]))
              w <- w[sample(length(w),min(round(length(w) * 0.8),1e4))]
            } else {
              w <- min(round(ncell(x) * 0.8),1e4)
              w <- sampleRandom(x,w,cells=TRUE)[,1]
            }
            xx <- x[w]
            Loop <- TRUE
            i <-  2
            o <- c()
            while (Loop) {
              cc <- cor(xx,categorize(xx,i),method='spearman')
              if (cc >= (1-th) | i > 100) Loop <- FALSE
              o <- c(o,cc)
              i <- i + 1
            }
            i <- i-1
            se <- sd(o,na.rm = TRUE) / sqrt(i)
            which(o > (max(o,na.rm = TRUE) - se))[1] + 1
            # o <- o[2:i] - o[1:(i-1)]
            # i <- i-1
            # o <- o[1:(i-1)] - o[2:i]
            # w <- which(abs(o) <= th)
            # if (length(w) > 0) return(w[1]+2)
            # else return(i+1)
          }
)

#------------ 

setMethod('nclass', signature(x='numeric'), 
          function(x,th=0.005) {
            if (missing(th)) th <- 0.005
            x <- x[which(!is.na(x))]
            Loop <- TRUE
            i <-  2
            o <- c()
            while (Loop) {
              cc <- cor(x,categorize(x,i),method='spearman')
              if (cc >= (1-th) | i > 100) Loop <- FALSE
              o <- c(o,cc)
              i <- i + 1
            }
            i <- i-1
            se <- sd(o,na.rm = TRUE) / sqrt(i)
            which(o > (max(o,na.rm = TRUE) - se))[1] + 1
            #o <- o[2:i] - o[1:(i-1)]
            #i <- i-1
            #o <- o[1:(i-1)] - o[2:i]
            #w <- which(abs(o) <= th)
            #if (length(w) > 0) return(w[1]+2)
            #else return(i+1)
          }
)
