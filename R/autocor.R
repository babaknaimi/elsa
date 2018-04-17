# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2016
# Version 1.1
# Licence GPL v3 

if (!isGeneric("moran")) {
  setGeneric("moran", function(x, d1, d2,...)
    standardGeneric("moran"))
}


setMethod('moran', signature(x='RasterLayer'), 
          function(x, d1, d2,...) {
            
            if (missing(d1)) d1 <- 0
            if (missing(d2)) {
              d2 <- res(x)[1]
              cat("Moran's I is calculated based on d1=",d1," & d2 =",d2,"(eual to ONE cell)\n")
            }
            w <-.Filter(r=res(x)[1],d1=d1,d2=d2)[[2]]
            
            .Call('moran',x[],as.integer(ncol(x)),as.integer(nrow(x)),as.integer(w[,1]),as.integer(w[,2]), PACKAGE='elsa')
          }
)

setMethod('moran', signature(x='Spatial'), 
          function(x, d1, d2,zcol,longlat,...) {
            
            if (missing(d1)) d1 <- 0
            if (missing(d2) && !inherits(d1,'neighbours')) stop('d2 should be specified, or alternatively, put an object in d1 created by dneigh')
            if (missing(longlat)) longlat <- NULL
            
            if (!inherits(d1,'neighbours')) d <- dneigh(x, d1, d2,longlat = longlat)@neighbours
            else d <- d1@neighbours
            
            if (missing(zcol)) {
              if (ncol(x@data) > 1) stop("zcol should be specified!")
              else zcol <- 1
            } else if (is.character(zcol)) {
              w <- which(colnames(x@data) == zcol[1])
              if (w == 0) stop('the specified variable in zcol does not exist in the data')
              zcol <- w
            } else if (is.numeric(zcol)) {
              zcol <- zcol[1]
              if (zcol > ncol(x@data)) stop('the zcol number is greater than the number of columns in data!')
            } else stop("zcol should be a character or a number!")
            
            x <- x@data[,zcol]
            
            if (!is.numeric(x) && !is.integer(x)) stop('the variable specified through zcol is not a numeric variable')
            
            .Call('moran_vector',x,d, PACKAGE='elsa')
          }
)

#---------

if (!isGeneric("geary")) {
  setGeneric("geary", function(x, d1, d2,...)
    standardGeneric("geary"))
}


setMethod('geary', signature(x='RasterLayer'), 
          function(x, d1, d2,...) {
            
            if (missing(d1)) d1 <- 0
            if (missing(d2)) {
              d2 <- res(x)[1]
              cat("Geary's c is calculated based on d1=",d1," & d2 =",d2,"(eual to ONE cell)\n")
            }
            w <-.Filter(r=res(x)[1],d1=d1,d2=d2)[[2]]
            
            .Call('geary',x[],as.integer(ncol(x)),as.integer(nrow(x)),as.integer(w[,1]),as.integer(w[,2]), PACKAGE='elsa')
          }
)

setMethod('geary', signature(x='Spatial'), 
          function(x, d1, d2,zcol,longlat,...) {
            
            if (missing(d1)) d1 <- 0
            if (missing(d2) && !inherits(d1,'neighbours')) stop('d2 should be specified, or alternatively, put an object in d1 created by dneigh')
            if (missing(longlat)) longlat <- NULL
            
            if (!inherits(d1,'neighbours')) d <- dneigh(x, d1, d2,longlat = longlat)@neighbours
            else d <- d1@neighbours
            
            if (missing(zcol)) {
              if (ncol(x@data) > 1) stop("zcol should be specified!")
              else zcol <- 1
            } else if (is.character(zcol)) {
              w <- which(colnames(x@data) == zcol[1])
              if (w == 0) stop('the specified variable in zcol does not exist in the data')
              zcol <- w
            } else if (is.numeric(zcol)) {
              zcol <- zcol[1]
              if (zcol > ncol(x@data)) stop('the specified number in zcol is greater than the number of columns in data!')
            } else stop("zcol should be a character or a number!")
            
            x <- x@data[,zcol,drop=TRUE]
            if (!is.numeric(x) && !is.integer(x)) stop('the variable specified through zcol is not a numeric variable')
            
            .Call('geary_vector',x,d, PACKAGE='elsa')
          }
)