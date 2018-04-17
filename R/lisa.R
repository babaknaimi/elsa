# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2016
# Version 1.1
# Licence GPL v3 

#----------
if (!isGeneric("lisa")) {
  setGeneric("lisa", function(x,d1,d2,statistic,...)
    standardGeneric("lisa"))
}


setMethod('lisa', signature(x='RasterLayer'), 
          function(x,d1=0,d2,statistic,mi,filename,...) {
            if (missing(statistic)) stop('statistic should be specified')
            if (length(statistic) > 1) {
              statistic <- statistic[1]
              warning('only first item in statistic is considered!')
            }
            if (missing(mi) || !tolower(mi) %in% c('i','ii','i.i','z','zi','z.i')) mi <- NULL
            else {
              if (tolower(mi) %in% c('i','ii','i.i')) mi <- 1
              else mi <- 2
            }
            
            if (!tolower(statistic) %in% c("i","c","g","g*",'localmoran','moran','localgeary','geary','localg','localg*')) stop("statistic should be one of: I (or localmoran), G, G*, and C (or LocalGeary)")
            if (missing(d2)) d2 <- res(x)[1]
            if (missing(d1)) d1 <- 0
            if (missing(filename)) filename=''
            #-----
            w <-.Filter(r=res(x)[1],d1=d1,d2=d2)
            fdim <- w[[1]]
            w <- w[[2]]
            
            if (fdim < 3) stop("d must be at least equal to the input raster resolution!")
            
            out <- raster(x)
            ncl <- ncol(out)
            nrw <- nrow(out)
            filename=trim(filename)
            
            if (canProcessInMemory(out)) {
              if (tolower(statistic) %in% c("i",'localmoran','moran')) {
                v <- .Call('localmoran', x[], as.integer(ncl), as.integer(nrw),as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                if (is.null(mi)) out[] <- v[[2]]
                else out[] <- v[[mi]]
              } else if (tolower(statistic) %in% c("c",'localgeary')) out[] <- .Call('localgeary', x[], as.integer(ncl), as.integer(nrw), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
              else if (tolower(statistic) %in% c("g",'localg',"g*",'localg*')) {
                v <- .Call('GG', x[], as.integer(ncl), as.integer(nrw), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                if (tolower(statistic) %in% c("g",'localg')) out[] <- v[[1]]
                else out[] <- v[[2]]
              }
              if (filename != '') out <- writeRaster(out, filename, ...)
            } else {
              tr <- blockSize(out, minblocks=3, minrows=fdim)
              
              if (tolower(statistic) %in% c("i",'localmoran','moran')) pb <- pbCreate(tr$n, label='LocalMoran',...)
              else if (tolower(statistic) %in% c("c",'localgeary')) pb <- pbCreate(tr$n, label='LocalGeary',...)
              else if (tolower(statistic) %in% c("g",'localg')) pb <- pbCreate(tr$n, label='LocalG',...)
              else if (tolower(statistic) %in% c("g*",'localg*')) pb <- pbCreate(tr$n, label='LocalG*',...)
              addr <- floor(fdim / 2)
              
              out <- writeStart(out, filename)
              v <- getValues(x, row=1, nrows=tr$nrows[1]+addr)
              
              if (tolower(statistic) %in% c("i",'localmoran','moran')) {
                v <- .Call('localmoran', x[], as.integer(ncl), as.integer(nrw),as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                if (is.null(mi)) v <- v[[2]]
                else v <- v[[mi]]
              } else if (tolower(statistic) %in% c("c",'localgeary')) v <- .Call('localgeary', v, as.integer(ncl), as.integer(nrw), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
              else if (tolower(statistic) %in% c("g",'localg',"g*",'localg*')) {
                v <- .Call('GG', v, as.integer(ncl), as.integer(nrw), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                if (tolower(statistic) %in% c("g",'localg')) v <- v[[1]]
                else v <- v[[2]]
              }
              ex <- length(v) - (addr * ncl)
              out <- writeValues(out, v[1:ex], 1)
              
              for (i in 2:(tr$n-1)) {
                v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr))
                
                if (tolower(statistic) %in% c("i",'localmoran','moran')) {
                  v <- .Call('localmoran', x[], as.integer(ncl), as.integer(nrw),as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                  if (is.null(mi)) v <- v[[2]]
                  else v <- v[[mi]]
                } else if (tolower(statistic) %in% c("c",'localgeary')) v <- .Call('localgeary', v, as.integer(ncl), as.integer(nrw), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                else if (tolower(statistic) %in% c("g",'localg',"g*",'localg*')) {
                  v <- .Call('GG', v, as.integer(ncl), as.integer(nrw), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                  if (tolower(statistic) %in% c("g",'localg')) v <- v[[1]]
                  else v <- v[[2]]
                }
                
                st <- (addr * ncl)+1
                ex <- length(v) - (addr * ncl)
                out <- writeValues(out, v[st:ex], tr$row[i])
                pbStep(pb)
              }
              i <- tr$n
              v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr)
              
              if (tolower(statistic) %in% c("i",'localmoran','moran')) {
                v <- .Call('localmoran', x[], as.integer(ncl), as.integer(nrw),as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                if (is.null(mi)) v <- v[[2]]
                else v <- v[[mi]]
              } else if (tolower(statistic) %in% c("c",'localgeary')) v <- .Call('localgeary', v, as.integer(ncl), as.integer(nrw), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
              else if (tolower(statistic) %in% c("g",'localg',"g*",'localg*')) {
                v <- .Call('GG', v, as.integer(ncl), as.integer(nrw), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                if (tolower(statistic) %in% c("g",'localg')) v <- v[[1]]
                else v <- v[[2]]
              }
              
              st <- (addr * ncl)+1
              ex <- length(v)
              out <- writeValues(out, v[st:ex], tr$row[i])
              pbStep(pb)
              out <- writeStop(out)
              pbClose(pb)
            }
            return(out)
          }
)
#---------

setMethod('lisa', signature(x='Spatial'), 
          function(x,d1,d2,statistic,zcol,longlat,drop,...) {
            
            if (!class(x) %in% c('SpatialPolygonsDataFrame','SpatialPointsDataFrame')) stop('x can only be either of RasterLayer, SpatialPointsDataFrame, SpatialPolygonsDataFrame')
            
            if (missing(statistic)) stop('statistic should be specified')
            if (length(statistic) > 1) {
              statistic <- statistic[1]
              warning('only first item in statistic is considered!')
            }
            
            if (missing(drop) || !is.logical(drop[1])) drop <- FALSE
            else drop <- drop[1]
            
            if (!tolower(statistic) %in% c("i","c","g","g*",'localmoran','moran','localgeary','geary','localg','localg*')) stop("statistic should be one of: I (or localmoran), G, G*, and C (or LocalGeary)")
            if (missing(d1)) d1 <- 0
            if (missing(d2) && !inherits(d1,'neighbours')) stop('d2 should be specified, or put an object created by dneigh in d1')
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
            
            xx <- x
            x <- x@data[,zcol]
            if (is.factor(x) || is.character(x)) stop('lisa statistics only apply on numerical variables')
            
            if (tolower(statistic) %in% c("i",'localmoran','moran')) {
              x <- .Call('localmoran_vector', x, d, PACKAGE='elsa')
              n <- 'LocalMoran'
            } else if (tolower(statistic) %in% c("c",'localgeary')) {
              x <- .Call('localgeary_vector', x, d, PACKAGE='elsa')
              n <- 'LocalGeray'
            }
            else if (tolower(statistic) %in% c("g",'localg',"g*",'localg*')) {
              x <- .Call('GG_vector', x, d, PACKAGE='elsa')
              if (tolower(statistic) %in% c("g",'localg')) {
                x <- x[[1]]
                n <- 'LocalG'
              }
              else {
                x <- x[[2]]
                n <- 'LocalG*'
              }
            }
            
            if (!drop) {
              if (tolower(statistic) %in% c("i",'localmoran','moran')) {
                xx@data$Ii <- x[[1]]
                xx@data$Z.Ii <- x[[2]]
                xx@data <- xx@data[,c('Ii','Z.Ii')]
              } else {
                xx@data$lisa <- x
                xx@data <- xx@data[,'lisa',drop=FALSE]
                colnames(xx@data) <- n
              }
              xx
            } else {
              if (tolower(statistic) %in% c("i",'localmoran','moran')) x[[2]]
              else x
            }
            
          }
)  
