# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2016
# Last Update : March 2023
# Version 1.8
# Licence GPL v3 
#-----------



if (!isGeneric("categorize")) {
  setGeneric("categorize", function(x,nc,probs,...)
    standardGeneric("categorize"))
}	

setMethod('categorize', signature(x='RasterLayer'), 
          function(x,nc,probs,filename='',verbose=TRUE,...)  {
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(nc)) {
              nc <- nclass(x)
              if (verbose) cat(paste("the optimum number of class has been identified as ",nc,"!\n"))
            }
            
            if (missing(probs)) probs <- FALSE
            else if (is.null(probs) || (is.logical(probs) && !probs)) probs <- FALSE
            else {
              if (is.numeric(probs) && length(probs) == 2 && all(probs <= 1) && all(probs >= 0) && probs[2] > probs[1]) {
                probs <- probs
              } else {
                warning('probs is not appropriately specified, e.g. c(0.025,0.975); NULL is considered')
                probs <- FALSE
              }
            }
            #-----
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- cellStats(x,'range')
              if (is.numeric(probs)) {
                # the quantile is used to avoid the effect of outlier on binning!
                .rq <- quantile(x,probs=probs)
                n <- (.rq[2] - .rq[1])/ nc
                nc <- seq(.rq[1],.rq[2],n)
                nc[1] <- r[1]
              } else {
                n <- (r[2] - r[1])/ nc
                nc <- seq(r[1],r[2],n)
              }
              
              nc[1] <- nc[1] - n
              if (nc[length(nc)] < r[2]) nc[length(nc)] <- r[2]
            }
            out <- raster(x)
            #-----
            if (canProcessInMemory(out)) {
              out[] <- .Call('categorize', as.vector(x[]), as.vector(nc), PACKAGE='elsa')
              if (filename != '') out <- writeRaster(out, filename, ...)
            } else {
              out <- writeStart(out, filename,...)
              tr <- blockSize(out, minblocks=3)
              pb <- pbCreate(tr$n, label='categorize',...)
              
              for (i in 1:tr$n) {
                v <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
                v <- .Call('categorize', v, as.vector(nc), PACKAGE='elsa')
                out <- writeValues(out, v, 1)
                pbStep(pb)
              }
              out <- writeStop(out)      
              pbClose(pb)
            }
            return(out)
          }
)
            

# in the following, we assume the layers in RasterStackBrick are the same variables in different times,
# therefore, the categorization should use the same base (class-number) for all layers (i.e., range of classes are specified based on the rangle of all values/layers together)
setMethod('categorize', signature(x='RasterStackBrick'), 
          function(x,nc,probs,filename='',verbose=TRUE,...)  {
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(nc)) {
              nc <- nclass(x)
              if (verbose) cat(paste("the optimum number of class has been identified as ",nc,"!\n"))
            }
            
            
            if (missing(probs)) probs <- FALSE
            else if (is.null(probs) || (is.logical(probs) && !probs)) probs <- FALSE
            else {
              if (is.numeric(probs) && length(probs) == 2 && all(probs <= 1) && all(probs >= 0) && probs[2] > probs[1]) {
                probs <- probs
              } else {
                warning('probs is not appropriately specified, e.g. c(0.025,0.975); NULL is considered')
                probs <- FALSE
              }
            }
            #-----
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- cellStats(x,'range')
              
              if (nlayers(x) > 1) {
                rr <- list()
                for (i in 1:nlayers(x)) {
                  rr[[i]] <- r[,i]
                }
              } else rr <- list(r[,1])
              
              
              if (is.numeric(probs)) {
                # the quantile is used to avoid the effect of outliers on binning!
                .rq <- t(quantile(x,probs=probs))
                ncl <- list()
                for (i in 1:ncol(.rq)) {
                  n <- (.rq[2,i] - .rq[1,i]) / nc
                  ncl[[i]] <- seq(.rq[1,i],.rq[2,i],n)
                  ncl[[i]][1] <- rr[[i]][1] - n
                  ncl[[i]][length(ncl[[i]])] <- rr[[i]][1]
                }
                
              } else {
                ncl <- list()
                for (i in 1:length(rr)) {
                  n <- (rr[[i]][2] - rr[[i]][1]) / nc
                  ncl[[i]] <- seq(rr[[i]][1],rr[[i]][2],n)
                  ncl[[i]][1] <- rr[[i]][1] - n
                }
              }
            }  else {
              if (is.numeric(nc)) {
                ncl <- vector('list',nlayers(x))
                for (i in 1:nlayers(x)) ncl[[i]] <- nc
              } else if (is.list(nc)) {
                if (length(nc) == nlayers(x)) ncl <- nc
                else stop('The provided list in nc has a different length than the number of layers in x!')
              } else stop('nc should be either a single number (number of class) or a numeric vector (list for multiple layers) specifying the categorisation thresholds!')
              
            }
            
            if (nlayers(x) > 1) {
              out <- brick(x,values=FALSE)
              for (i in 1:nlayers(x)) {
                out[[i]][] <- .Call('categorize', as.vector(x[[i]][]), as.vector(ncl[[i]]), PACKAGE='elsa')
              }
            } else {
              out <- raster(x)
              out[] <- .Call('categorize', as.vector(x[]), as.vector(ncl[[1]]), PACKAGE='elsa')
            }
            #-----
            names(out) <- names(x)
            
            if (filename != '') {
              writeRaster(out,filename=filename,...)
            }
            
            return(out)
          }
)


#-----------

setMethod('categorize', signature(x='SpatRaster'), 
          function(x,nc,probs,filename='',verbose=TRUE,...)  {
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(nc)) {
              nc <- nclass(x[[1]])
              if (nlyr(x) > 1) {
                if (verbose) cat(paste("the optimum number of class has been identified as ",nc," (ONLY the first layer is considered)!\n"))
              } else {
                if (verbose) cat(paste("the optimum number of class has been identified as ",nc,"!\n"))
              }
            }
            
            
            if (missing(probs)) probs <- FALSE
            else if (is.null(probs) || (is.logical(probs) && !probs)) probs <- FALSE
            else {
              if (is.numeric(probs) && length(probs) == 2 && all(probs <= 1) && all(probs >= 0) && probs[2] > probs[1]) {
                probs <- probs
              } else {
                warning('probs is not appropriately specified, e.g. c(0.025,0.975); NULL is considered')
                probs <- FALSE
              }
            }
            #-----
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- global(x,'range',na.rm=TRUE)
              if (nlyr(x) > 1) {
                rr <- list()
                for (i in 1:nlyr(x)) {
                  rr[[i]] <- t(r)[,i]
                }
              } else rr <- list(t(r)[,1])
              
              #---
              if (is.numeric(probs)) {
                # the quantile is used to avoid the effect of outliers on binning!
                .rq <- t(global(x,fun=quantile,probs=probs,na.rm=TRUE))
                ncl <- list()
                for (i in 1:ncol(.rq)) {
                  n <- (.rq[2,i] - .rq[1,i]) / nc
                  ncl[[i]] <- seq(.rq[1,i],.rq[2,i],n)
                  ncl[[i]][1] <- rr[[i]][1] - n
                  ncl[[i]][length(ncl[[i]])] <- rr[[i]][1]
                }
                
              } else {
                ncl <- list()
                for (i in 1:length(rr)) {
                  n <- (rr[[i]][2] - rr[[i]][1]) / nc
                  ncl[[i]] <- seq(rr[[i]][1],rr[[i]][2],n)
                  ncl[[i]][1] <- rr[[i]][1] - n
                }
              }
              
              # nc[1] <- nc[1] - n
              # if (nc[length(nc)] < r[2]) nc[length(nc)] <- r[2]
            } else {
              if (is.numeric(nc)) {
                ncl <- vector('list',nlyr(x))
                for (i in 1:nlyr(x)) ncl[[i]] <- nc
              } else if (is.list(nc)) {
                if (length(nc) == nlyr(x)) ncl <- nc
                else stop('The provided list in nc has a different length than the number of layers in x!')
              } else stop('nc should be either a single number (number of class) or a numeric vector (list for multiple layers) specifying the categorisation thresholds!')
              
            }
            
            out <- rast(x)
            #-----
            
            for (i in 1:nlyr(x)) {
              out[[i]][] <- .Call('categorize', as.vector(x[[i]][]), as.vector(ncl[[i]]), PACKAGE='elsa')
            }
            names(out) <- names(x)
            
            if (filename != '') writeRaster(out,filename = filename,...)
            
            
            return(out)
          }
)


#-----------
setMethod('categorize', signature(x='numeric'), 
          function(x,nc,probs)  {
            if (missing(nc)) {
              stop("number of classes or a verctor including the break values should be specified...!")
            }
            
            if (missing(probs)) probs <- FALSE
            else if (is.null(probs) || (is.logical(probs) && !probs)) probs <- FALSE
            else {
              if (is.numeric(probs) && length(probs) == 2 && all(probs <= 1) && all(probs >= 0) && probs[2] > probs[1]) {
                probs <- probs
              } else {
                warning('probs is not appropriately specified, e.g. c(0.025,0.975); NULL is considered')
                probs <- FALSE
              }
            }
            #-----
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- range(x,na.rm=TRUE)
              if (is.numeric(probs)) {
                # the quantile is used to avoid the effect of outlier on binning!
                .rq <- quantile(x,probs=probs,na.rm = TRUE)
                n <- (.rq[2] - .rq[1])/ nc
                nc <- seq(.rq[1],.rq[2],n)
                nc[1] <- r[1]
              } else {
                n <- (r[2] - r[1])/ nc
                nc <- seq(r[1],r[2],n)
              }
              
              nc[1] <- nc[1] - n
              if (nc[length(nc)] < r[2]) nc[length(nc)] <- r[2]
            }
            
            .Call('categorize', as.vector(x), as.vector(nc), PACKAGE='elsa')
          }
)



setMethod('categorize', signature(x='list'), 
          function(x,nc,probs)  {
            if (missing(nc)) {
              stop("number of classes or a verctor including the break values should be specified...!")
            }
            
            if (missing(probs)) probs <- FALSE
            else if (is.null(probs) || (is.logical(probs) && !probs)) probs <- FALSE
            else {
              if (is.numeric(probs) && length(probs) == 2 && all(probs <= 1) && all(probs >= 0) && probs[2] > probs[1]) {
                probs <- probs
              } else {
                warning('probs is not appropriately specified, e.g. c(0.025,0.975); NULL is considered')
                probs <- FALSE
              }
            }
            #-----
            
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- lapply(x,range,na.rm=TRUE)
              r <- c(min(sapply(r,function(x) x[1]),na.rm=TRUE),max(sapply(r,function(x) x[2]),na.rm=TRUE))
              if (is.numeric(probs)) {
                # the quantile is used to avoid the effect of outlier on binning!
                .rq <- lapply(x,quantile,probs=probs,na.rm=TRUE)
                .rq <- c(min(sapply(.rq,function(x) x[1]),na.rm=TRUE),max(sapply(.rq,function(x) x[2]),na.rm=TRUE))
                n <- (.rq[2] - .rq[1])/ nc
                nc <- seq(.rq[1],.rq[2],n)
                nc[1] <- r[1]
              } else {
                n <- (r[2] - r[1]) / nc
                nc <- seq(r[1],r[2],n)
              }
              
              nc[1] <- nc[1] - n
              if (nc[length(nc)] < r[2]) nc[length(nc)] <- r[2]
            }
            
            o <- list()
            for (i in 1:length(x)) {
              o[[i]] <- .Call('categorize', as.vector(x[[i]]), as.vector(nc), PACKAGE='elsa')
            }
            o
          }
)
