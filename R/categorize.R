# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2016
# Version 1.2
# Licence GPL v3 



if (!isGeneric("categorize")) {
  setGeneric("categorize", function(x,nc,...)
    standardGeneric("categorize"))
}	

setMethod('categorize', signature(x='RasterLayer'), 
          function(x,nc,filename='',...)  {
            if (missing(nc)) {
              nc <- nclass(x)
              cat(paste("the optimum number of class has been identified as ",nc,"!\n"))
            }
            
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- cellStats(x,'range')
              n <- (r[2] - r[1])/ nc
              nc <- seq(r[1],r[2],n)
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
            



setMethod('categorize', signature(x='numeric'), 
          function(x,nc)  {
            if (missing(nc)) {
              stop("number of classes or a verctor including the break values should be specified...!")
            }
            
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- range(x,na.rm=TRUE)
              n <- (r[2] - r[1])/ nc
              nc <- seq(r[1],r[2],n)
              nc[1] <- nc[1] - n
              if (nc[length(nc)] < r[2]) nc[length(nc)] <- r[2]
            }
            .Call('categorize', as.vector(x), as.vector(nc), PACKAGE='elsa')
          }
)

