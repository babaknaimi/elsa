# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2014
# Version 1.0
# Licence GPL v3 

if (!isGeneric("dif2list")) {
  setGeneric("dif2list", function(x, pattern)
    standardGeneric("dif2list"))
}



setMethod('dif2list', signature(x='data.frameORmatrix'), 
          function(x, pattern) {
            
            x <- x[,1:2]
            
            .f <- function(code,d) {
              d <- t(apply(d,1,function(x) {abs(x - code)}))
              nc <- ncol(d)
              ss <- rep(0,nrow(d))
              for (i in 1:nrow(d)) {
                j <- 1
                while (j <= nc) {
                  if (d[i,j] != 0) {
                    ss[i] <-(nc - j + 1)
                    j <- nc + 1
                  } else j <- j + 1
                }
              }
              ss
            }
            
            if (missing(pattern)) {
              u <- unlist(strsplit(as.character(x[1,2]),''))
              pattern <- rep(1,length(u))
            }
            
            p <- list()
            o <- 1
            for (j in 1:length(pattern)) {
              p[[j]] <- c(o:(o+pattern[j]-1))
              o <- (j+pattern[j])
            }
            
            s <- sapply(x[,2],function(x) {strsplit(as.character(x),'')})
            if (!all(sapply(s,function(x) {length(x) == sum(pattern)}))) stop("the provided codes does not match the pattern or have inconsistency!")
            
            d <- data.frame(matrix(nrow=length(s),ncol=length(pattern)))
            for (i in 1:length(s)) {
              for (j in 1:length(pattern)) {
                d[i,j] <- as.numeric(paste(s[[i]][p[[j]]],collapse=''))
              }
            }
            gc <- x[,1]
            dT <- list()
            for (i in 1:length(gc)) {
              n <- .f(as.numeric(d[i,]),d)
              names(n) <- gc
              dT[[as.character(gc[i])]] <- n
            }
            dT
          }
)
