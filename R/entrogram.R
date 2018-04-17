# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2016
# Version 2.1
# Licence GPL v3 

.checkrasterMemory <- function(cells,n=1) {
  cells <- ceiling(sqrt(cells))
  canProcessInMemory(raster(nrows=cells, ncols=cells, xmn=0, xmx=cells,vals=NULL),n)
}

if (!isGeneric("entrogram")) {
  setGeneric("entrogram", function(x, width, cutoff,...)
    standardGeneric("entrogram"))
}


setMethod('entrogram', signature(x='RasterLayer'), 
          function(x, width, cutoff, categorical, nc, dif, cloud=FALSE, s=NULL,...) {
            re <- res(x)[1]
            if (missing(cutoff)) cutoff<- sqrt((xmin(x)-xmax(x))^2+(ymin(x)-ymax(x))^2) / 3
            if (missing(width)) width <- re
            else if (width < re) width <- re
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
                  s <- sampleRandom(x,s,cells=TRUE)[,1]
                } else {
                  s <- 1 / (2 * nlag)
                  s <- ceiling(n * s)
                  while (!.checkrasterMemory(s,1)) s <- ceiling(s / 2)
                  s <- sampleRandom(x,s,cells=TRUE)[,1]
                }
              } else {
                s <- (1:ncell(x))[which(!is.na(x[]))]
              }
            } else {
              if (!is.numeric(s)) stop("s argument should be an integer number or NULL!")
              while (!.checkrasterMemory(s[1],1)) s <- ceiling(s[1] * 0.8)
              if (s > n) s <- n
              s <- sampleRandom(x,s,cells=TRUE)[,1]
            }
            #######---------------
            #----
            if (!missing(nc)) {
              if (missing(categorical)) {
                if (missing(dif)) categorical <- FALSE
                else {
                  categorical <- TRUE
                  cat("input data is considered categorical, and nc is ignored!\n")
                }
              }
            } else {
              if (missing(categorical) && !missing(dif)) categorical <- TRUE
            }
            #----
            if (missing(categorical) || !is.logical(categorical)) {
              # guessing whether the layer is categorical:
              if (.is.categorical(x)) {
                categorical <- TRUE
                cat("the input is considered as a categorical variable...\n")
              } else {
                categorical <- FALSE
                cat("the input is considered as a continuous variable...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
            } else if (categorical) {
              classes <- unique(x)
              nc <- length(classes)
            }
            #-----
            
            if (categorical) {
              if (missing(dif)) {
                dif <- rep(1,nc*nc)
                for (i in 1:nc) dif[(i-1)*nc+i] <-0
              } else {
                dif <- .checkDif(dif,classes)
              }
            }
            #-----
            #######---------------------
            
            if (!categorical) x <- categorize(x,nc)
            
            ncl <- ncol(x)
            nrw <- nrow(x)
            
            out <- new("Entrogram")
            out@width <- width
            out@cutoff <- cutoff
            if (cloud) {
              out@entrogramCloud <- matrix(NA,nrow=length(s),ncol=nlag)
              for (i in 1:nlag) {
                w <-.Filter(r=res(x)[1],d1=0,d2=i*width)
                w <- w[[2]]
                if (categorical) {
                  out@entrogramCloud[,i] <- .Call('v_elsac_cell', as.integer(x[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(s), PACKAGE='elsa')
                } else {
                  out@entrogramCloud[,i] <- .Call('v_elsa_cell', as.integer(x[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(s), PACKAGE='elsa')
                }
                #out@entrogramCloud[,i] <- elsa(x,d=i*width,nc=nc,categorical=categorical,dif=dif,cells=s)
              }
              out@entrogram <- data.frame(distance=seq(width,width*nlag,width) - (width/2),E=apply(out@entrogramCloud,2,mean,na.rm=TRUE))
            } else {
              d <- seq(width,width*nlag,width) - (width/2)
              out@entrogram <- data.frame(distance=d,E=rep(NA,length(d)))
              for (i in 1:nlag) {
                w <-.Filter(r=res(x)[1],d1=0,d2=i*width)[[2]]
                if (categorical) {
                  out@entrogram [i,2] <- mean(.Call('v_elsac_cell', as.integer(x[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(s), PACKAGE='elsa'),na.rm=TRUE)
                } else {
                  out@entrogram [i,2] <- mean(.Call('v_elsa_cell', as.integer(x[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(s), PACKAGE='elsa'),na.rm=TRUE)
                }
                #out@entrogram [i,2] <- mean(elsa(x,d=i*width,nc=nc,categorical=categorical,dif=dif,cells=s),na.rm=TRUE)
              }
            }
            out
          }
)
##########

#-------
setMethod('entrogram', signature(x='SpatialPolygonsDataFrame'), 
          function(x, width, cutoff, categorical, nc, dif, zcol,  cloud=FALSE, s=NULL,method,longlat,...) {
            n <- nrow(x)
            
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
            
            if (missing(method)) method <- 'centroid'
            else {
              if (tolower(method)[1] %in% c('bnd','bound','boundary','bond','b')) method <- 'bound'
              else method <- 'centroid'
            }
            
            if (method == 'centroid') xy <- coordinates(x)
            else xy <- x
            
            x <- x@data[,zcol]
            #---
            if (!is.null(s) && is.numeric(s) && s < n) {
              x <- x[sample(n,s)]
              n <- length(n)
            } 
            #######---------------
            #----
            if (!missing(nc)) {
              if (missing(categorical)) {
                if (missing(dif)) categorical <- FALSE
                else {
                  categorical <- TRUE
                  cat("input data is considered categorical, and nc is ignored!\n")
                }
              }
            } else {
              if (missing(categorical) && !missing(dif)) categorical <- TRUE
            }
            #----
            if (missing(categorical) || !is.logical(categorical)) {
              # guessing whether the layer is categorical:
              if (.is.categorical(x)) {
                categorical <- TRUE
                cat("the input is considered as a categorical variable...\n")
              } else {
                categorical <- FALSE
                cat("the input is considered as a continuous variable...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
              classes <- 1:nc
            } else if (categorical) {
              classes <- unique(x)
              nc <- length(classes)
            }
            #-----
            
            if (categorical) {
              if (missing(dif)) {
                dif <- rep(1,nc*nc)
                for (i in 1:nc) dif[(i-1)*nc+i] <-0
              } else {
                dif <- .checkDif(dif,classes)
              }
            }
            #-----
            if (!categorical) x <- categorize(x,nc)
            #######---------------------
            out <- new("Entrogram")
            out@width <- width
            out@cutoff <- cutoff
            if (cloud) {
              out@entrogramCloud <- matrix(NA,nrow=n,ncol=nlag)
              for (i in 1:nlag) {
                d <- dneigh(xy,d1=0,d2=i*width,method = method,longlat = longlat)@neighbours
                if (categorical) {
                  out@entrogramCloud[,i] <- .Call('v_elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif, PACKAGE='elsa')
                } else {
                  out@entrogramCloud[,i] <-.Call('v_elsa_vector', as.integer(x), d, as.integer(nc), PACKAGE='elsa')
                }
              }
              out@entrogram <- data.frame(distance=seq(width,width*nlag,width) - (width/2),E=apply(out@entrogramCloud,2,mean,na.rm=TRUE))
            } else {
              d <- seq(width,width*nlag,width) - (width/2)
              out@entrogram <- data.frame(distance=d,E=rep(NA,length(d)))
              for (i in 1:nlag) {
                d <- dneigh(xy,d1=0,d2=i*width,method = method,longlat = longlat)@neighbours
                if (categorical) {
                  out@entrogram [i,2] <- mean(.Call('v_elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif),na.rm=TRUE, PACKAGE='elsa')
                } else {
                  out@entrogram [i,2] <- mean(.Call('v_elsa_vector', as.integer(x), d, as.integer(nc)),na.rm=TRUE, PACKAGE='elsa')
                }
              }
            }
            out
          }
)



setMethod('entrogram', signature(x='SpatialPointsDataFrame'), 
          function(x, width, cutoff, categorical, nc, dif, zcol,  cloud=FALSE, s=NULL,longlat,...) {
            n <- nrow(x)
            
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
              x <- x[sample(n,s)]
              n <- length(n)
            } 
            #######---------------
            #----
            if (!missing(nc)) {
              if (missing(categorical)) {
                if (missing(dif)) categorical <- FALSE
                else {
                  categorical <- TRUE
                  cat("input data is considered categorical, and nc is ignored!\n")
                }
              }
            } else {
              if (missing(categorical) && !missing(dif)) categorical <- TRUE
            }
            #----
            if (missing(categorical) || !is.logical(categorical)) {
              # guessing whether the layer is categorical:
              if (.is.categorical(x)) {
                categorical <- TRUE
                cat("the input is considered as a categorical variable...\n")
              } else {
                categorical <- FALSE
                cat("the input is considered as a continuous variable...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
              classes <- 1:nc
            } else if (categorical) {
              classes <- unique(x)
              nc <- length(classes)
            }
            #-----
            
            if (categorical) {
              if (missing(dif)) {
                dif <- rep(1,nc*nc)
                for (i in 1:nc) dif[(i-1)*nc+i] <-0
              } else {
                dif <- .checkDif(dif,classes)
              }
            }
            #-----
            if (!categorical) x <- categorize(x,nc)
            #######---------------------
            out <- new("Entrogram")
            out@width <- width
            out@cutoff <- cutoff
            if (cloud) {
              out@entrogramCloud <- matrix(NA,nrow=n,ncol=nlag)
              for (i in 1:nlag) {
                d <- dneigh(xy,d1=0,d2=i*width,longlat = longlat)@neighbours
                if (categorical) {
                  out@entrogramCloud[,i] <- .Call('v_elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif)
                } else {
                  out@entrogramCloud[,i] <-.Call('v_elsa_vector', as.integer(x), d, as.integer(nc))
                }
              }
              out@entrogram <- data.frame(distance=seq(width,width*nlag,width) - (width/2),E=apply(out@entrogramCloud,2,mean,na.rm=TRUE))
            } else {
              d <- seq(width,width*nlag,width) - (width/2)
              out@entrogram <- data.frame(distance=d,E=rep(NA,length(d)))
              for (i in 1:nlag) {
                d <- dneigh(xy,d1=0,d2=i*width,longlat = longlat)@neighbours
                if (categorical) {
                  out@entrogram [i,2] <- mean(.Call('v_elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif),na.rm=TRUE)
                } else {
                  out@entrogram [i,2] <- mean(.Call('v_elsa_vector', as.integer(x), d, as.integer(nc)),na.rm=TRUE)
                }
              }
            }
            out
          }
)
