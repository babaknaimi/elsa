# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2016
# Latest Update: May 2019
# Version 1.5
# Licence GPL v3 

if (!isGeneric("Variogram")) {
  setGeneric("Variogram", function(x,width,cutoff,...)
    standardGeneric("Variogram"))
}


setMethod('Variogram', signature(x='RasterLayer'), 
          function(x, width, cutoff, cloud=FALSE,s=NULL,...) {
            re <- res(x)[1]
            if (missing(cutoff)) cutoff<- sqrt((xmin(x)-xmax(x))^2+(ymin(x)-ymax(x))^2) / 3
            if (missing(width)) width <- re
            else if (width < re) width <- re
            
            if (missing(s)) s <- NULL
            
            if (cutoff < width) stop("cutoff should be greater than width size")
            nlag <- ceiling(cutoff / width)
            
            n <- ncell(x) - cellStats(x,'countNA')
            #---
            if (is.null(s)) {
              if (!.checkrasterMemory(n,nlag)) {
                s <- c()
                for (i in (nlag-1):1) s <- c(s,.checkrasterMemory(n,i))
                s <- which(s)
                if (length(s) > 0) {
                  s <- (nlag - s[1]) / (2*nlag)
                  s <- ceiling(n * s)
                  #s <- sampleRandom(x,s,cells=TRUE)[,1]
                  s <- sampleRandom(x,s,sp=TRUE)
                } else {
                  s <- 1 / (2 * nlag)
                  s <- ceiling(n * s)
                  while (!.checkrasterMemory(s,1)) s <- ceiling(s / 2)
                  #s <- sampleRandom(x,s,cells=TRUE)[,1]
                  s <- sampleRandom(x,s,sp=TRUE)
                }
              } 
            } else {
              if (!is.numeric(s)) stop("s argument should be an integer number or NULL!")
              while (!.checkrasterMemory(s[1],1)) s <- ceiling(s[1] * 0.8)
              if (s > n) s <- n
              #s <- sampleRandom(x,s,cells=TRUE)[,1]
              s <- sampleRandom(x,s,sp=TRUE)
            }
            
            #######---------------------
            if (is.null(s)) {
              out <- new("Variogram")
              out@width <- width
              out@cutoff <- cutoff
              d <- seq(width,width*nlag,width) - (width/2)
              out@variogram <- data.frame(distance=d,gamma=rep(NA,length(d)))
              if (cloud) out@variogramCloud <- matrix(NA,nrow=if (is.null(s)) length(x) else length(s),ncol=nlag)
              for (i in 1:nlag) {
                d1 <- (i -1) * width
                d2 <- d1 + width
                w <-.Filter(r=res(x)[1],d1=d1,d2=d2)[[2]]
                w <- .Call('semivar',as.vector(x[]),as.integer(ncol(x)),as.integer(nrow(x)),as.integer(w[,1]),as.integer(w[,2]), PACKAGE='elsa')
                if (cloud) out@variogramCloud[,i] <- w
                w <- w[!is.infinite(w)]
                w <- w[!is.na(w)]
                out@variogram [i,2] <- mean(w)
              }
            } else {
              if (!is.na(projection(x))) {
                longlat <- strsplit(trim(strsplit(projection(x),'\\+')[[1]][2]),'=')[[1]][2] == 'longlat'
                if (is.na(longlat) || !is.logical(longlat)) longlat <- NULL
              } else longlat <- NULL
              out <- Variogram(s,width,cutoff,zcol=names(x),cloud=cloud,s=NULL,longlat=longlat,...)
            }
            out
          }
)
##########


setMethod('Variogram', signature(x='Spatial'), 
          function(x, width, cutoff, zcol,  cloud=FALSE, s=NULL,longlat,...) {
            if (!class(x) %in% c('SpatialPolygonsDataFrame','SpatialPointsDataFrame')) stop('x can only be either of RasterLayer, SpatialPointsDataFrame, SpatialPolygonsDataFrame')
            
            n <- nrow(x)
            
            if (missing(s)) s <- NULL
            
            if (missing(longlat)) longlat <- NULL
            
            if (missing(cutoff)) cutoff<- sqrt((xmin(x)-xmax(x))^2+(ymin(x)-ymax(x))^2) / 3
            if (missing(width)) width <- cutoff / 15
            
            if (cutoff < width) stop("cutoff should be greater than width size")
            
            nlag <- ceiling(cutoff / width)
            
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
            
            xy <- coordinates(x)
            x <- x@data[,zcol]
            #---
            if (!is.null(s) && is.numeric(s) && s < n) {
              s <- sample(n,s)
              x <- x[s]
              n <- length(x)
              xy <- xy[s,]
            }
            #######---------------
            out <- new("Variogram")
            out@width <- width
            out@cutoff <- cutoff
            d <- seq(width,width*nlag,width) - (width/2)
            out@variogram <- data.frame(distance=d,gamma=rep(NA,length(d)))
            if (cloud) {
              out@variogramCloud <- matrix(NA,nrow=n,ncol=nlag)
              for (i in 1:nlag) {
                d1 <- (i -1) * width
                d2 <- d1 + width
                d <- dneigh(xy,d1=d1, d2=d2,longlat = longlat)@neighbours
                w <- .Call('semivar_vector', x, d, PACKAGE='elsa')
                out@variogramCloud[,i] <-w
                w <- w[!is.infinite(w)]
                w <- w[!is.na(w)]
                out@variogram [i,2] <- mean(w,na.rm=TRUE)
              }
            } else {
              for (i in 1:nlag) {
                d1 <- (i -1) * width
                d2 <- d1 + width
                d <- dneigh(xy,d1=d1, d2=d2,longlat = longlat)@neighbours
                w <- .Call('semivar_vector', x, d, PACKAGE='elsa')
                w <- w[!is.infinite(w)]
                w <- w[!is.na(w)]
                out@variogram [i,2] <- mean(w,na.rm=TRUE)
              }
            }
            out
          }
)

