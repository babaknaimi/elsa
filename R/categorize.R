# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2016
# Last Update : Feb. 2020
# Version 1.6
# Licence GPL v3 



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
              r <- c(apply(r,1,min)[1],apply(r,1,max)[2])
              
              if (is.numeric(probs)) {
                # the quantile is used to avoid the effect of outlier on binning!
                .rq <- t(quantile(x,probs=probs))
                .rq <- c(apply(.rq,1,min)[1],apply(.rq,1,max)[2])
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
            if (filename != '') {
              out <- brick(out,nl=nlayers(x))
              writeRaster(out,filename=filename,...)
              out <- brick(filename)
              
              for (i in 1:nlayers(x)) {
                xx <- .Call('categorize', as.vector(x[[i]][]), as.vector(nc), PACKAGE='elsa')
                out <- update(out,xx,cell=1,band=i)
              }
            } else {
              for (i in 1:nlayers(x)) {
                xx <- raster(out)
                xx[] <- .Call('categorize', as.vector(x[[i]][]), as.vector(nc), PACKAGE='elsa')
                out <- addLayer(out,xx)
              }
            }
            names(out) <- names(x)
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
