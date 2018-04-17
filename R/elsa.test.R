# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2016
# Version 2.3
# Licence GPL v3 

if (!isGeneric("elsa.test")) {
  setGeneric("elsa.test", function(x, d, n, method, null, nc, categorical, dif,...)
    standardGeneric("elsa.test"))
}


setMethod('elsa.test', signature(x='RasterLayer'), 
          function(x, d, n=99, method, null, nc, categorical, dif,cells,filename,...) {
            
            if (missing(filename)) filename <- ''
            
            if (missing(n)) {
              if (ncell(x) > 20000) n <- 99
              else n <- 999
              
              cat(paste("n (number of runs in Monte Carlo simulations) is set to",n,"...\n"))
            }
            #----------
            if (missing(method)) method <- 2
            else {
              method <- method[1]
              if (method %in% c('boot','bootstrap','b','bo')) method <- 2
              else if (method %in% c('perm','permutation','p','pe')) method <- 1
              else {
                if (!is.numeric(method) || !method %in% 1:2) {
                  warning('method is not identified; default ("boot") is considered')
                  method <- 2
                }
              }
            }
            #------
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
            if (!categorical) {
              if (missing(nc)) nc <- nclass(x)
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
            
            if (missing(null)) {
              null <- calc(x,function(x) { x[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE); x})
            } else if (inherits(null,'numeric') && length(null) == ncell(x)) {
              nullx <- null
              null <- raster(x)
              null <- setValues(null,nullx)
              rm(nullx)
            } else if ((inherits(null,'RasterLayer') && !compareRaster(x,null,crs=FALSE,stopiffalse=FALSE)) || !inherits(null,'RasterLayer')) {
              warning('null is not a numeric vector, or a raster, or is a raster with a different extent, resolution, etc.; so, the null is generated given the default settings!')
              null <- calc(x,function(x) { x[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE); x})
            }
              
            #----------------
            #-----
            w <-.Filter(r=res(x)[1],d1=0,d2=d)
            fdim <- w[[1]]
            w <- w[[2]]
            
            if (fdim < 3) stop("d must be at least equal to the input raster resolution!")
            
            if (!categorical) x <- categorize(x,nc)
            
            out <- raster(x)
            ncl <- ncol(out)
            nrw <- nrow(out)
            filename=trim(filename)
            
            if (canProcessInMemory(out)) {
              if (categorical) {
                if (missing(cells)) {
                  out[] <- .Call('elsac_test', x[],as.vector(null[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif,as.integer(method),as.integer(n), PACKAGE='elsa')
                  if (filename != '') out <- writeRaster(out, filename, ...)
                } else {
                  out <- .Call('elsac_cell_test', x[],as.vector(null[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif,as.integer(cells),as.integer(method),as.integer(n), PACKAGE='elsa')
                }
              } else {
                if (missing(cells)) {
                  out[] <- .Call('elsa_test', x[],as.vector(null[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(method),as.integer(n), PACKAGE='elsa')
                  if (filename != '') out <- writeRaster(out, filename, ...)
                } else {
                  out <- .Call('elsa_cell_test', x[],as.vector(null[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells),as.integer(method),as.integer(n), PACKAGE='elsa')
                }
              }
            } else {
              tr <- blockSize(out, minblocks=3, minrows=fdim)
              pb <- pbCreate(tr$n, label='ELSA',...)
              addr <- floor(fdim / 2)
              
              if (missing(cells)) {
                out <- writeStart(out, filename)
                v <- getValues(x, row=1, nrows=tr$nrows[1]+addr)
                vn <- getValues(null, row=1, nrows=tr$nrows[1]+addr)
                if (!categorical) {
                  v <- .Call('elsa_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(method),as.integer(n), PACKAGE='elsa')
                } else {
                  v <- .Call('elsac_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif,as.integer(method),as.integer(n), PACKAGE='elsa')
                }
                ex <- length(v) - (addr * ncl)
                out <- writeValues(out, v[1:ex], 1)
                
                for (i in 2:(tr$n-1)) {
                  v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr))
                  vn <- getValues(null, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr))
                  if (!categorical) {
                    v <- .Call('elsa_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(method),as.integer(n), PACKAGE='elsa')
                  } else {
                    v <- .Call('elsac_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif,as.integer(method),as.integer(n), PACKAGE='elsa')
                  }
                  st <- (addr * ncl) + 1
                  ex <- length(v) - (addr * ncl)
                  out <- writeValues(out, v[st:ex], tr$row[i])
                  pbStep(pb)
                }
                
                i <- tr$n
                v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr)
                vn <- getValues(null, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr)
                if (!categorical) {
                  v <- .Call('elsa_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(method),as.integer(n), PACKAGE='elsa')
                } else {
                  v <- .Call('elsac_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif,as.integer(method),as.integer(n), PACKAGE='elsa')
                }
                st <- (addr * ncl)+1
                ex <- length(v)
                out <- writeValues(out, v[st:ex], tr$row[i])
                pbStep(pb)
                out <- writeStop(out)      
                pbClose(pb)  
              } else {
                v <- getValues(x, row=1, nrows=tr$nrows[1]+addr)
                vn <- getValues(null, row=1, nrows=tr$nrows[1]+addr)
                cls <- cells[which(cells <= (tr$nrows[1]) * ncl)]
                if (length(cls) > 0) {
                  if (!categorical) {
                    v <- .Call('elsa_cell_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cls),as.integer(method),as.integer(n), PACKAGE='elsa')
                  } else {
                    v <- .Call('elsac_cell_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif,as.integer(cls), as.integer(method),as.integer(n), PACKAGE='elsa')
                  }
                  out <- c(out, v)
                }
                
                for (i in 2:(tr$n-1)) {
                  cls <- cells[which(cells > ((tr$nrow[i] - 1) * ncl) & cells <= ((tr$nrows[i]+ tr$nrows[i] - 1) * ncl))]
                  if (length(cls) > 0) {
                    cls <- cls - ((tr$row[i]-addr-1)*ncl)
                    v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr))
                    vn <- getValues(null, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr))
                    
                    if (!categorical) {
                      v <- .Call('elsa_cell_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cls),as.integer(method),as.integer(n), PACKAGE='elsa')
                    } else {
                      v <- .Call('elsac_cell_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif,as.integer(cls), as.integer(method),as.integer(n), PACKAGE='elsa')
                    }
                    
                    out <- c(out, v)
                    pbStep(pb)
                  }
                }
                
                i <- tr$n
                cls <- cells[which(cells > ((tr$nrow[i] - 1) * ncl) & cells <= ((tr$nrows[i]+ tr$nrows[i] - 1) * ncl))]
                cls <- cls - ((tr$row[i]-addr-1)*ncl)
                if (length(cls) > 0) {
                  v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr)
                  vn <- getValues(null, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr)
                  if (!categorical) {
                    v <- .Call('elsa_cell_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cls),as.integer(method),as.integer(n), PACKAGE='elsa')
                  } else {
                    v <- .Call('elsac_cell_test', v , vn, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif,as.integer(cls), as.integer(method),as.integer(n), PACKAGE='elsa')
                  }
                  out <- c(out, v)
                  pbStep(pb)
                  pbClose(pb)  
                }
              }
            }
            #----------
            return(out)
          }
)

#--------------
# 
# 
# setMethod('elsa.test', signature(x='SpatialPointsDataFrame'), 
#           function(x, d, n=99, method, null, nc, categorical, dif,zcol,...) {
#             
#             if (missing(d)) stop('d is missed!')
#             else if (!class(d) %in% c('numeric','integer','neighbours')) stop('d should be either a number (distance) or an object of class neighbours (created by dneigh function')
#             
#             if (!inherits(d,'neighbours')) d <- dneigh(x, d[1])
#             d <- d@neighbours
#             #-------
#             if (missing(zcol)) {
#               if (ncol(x@data) > 1) stop("zcol should be specified!")
#               else zcol <- 1
#             } else if (is.character(zcol)) {
#               w <- which(colnames(x@data) == zcol[1])
#               if (w == 0) stop('the specified variable in zcol does not exist in the data')
#               zcol <- w
#             } else if (is.numeric(zcol)) {
#               zcol <- zcol[1]
#               if (zcol > ncol(x@data)) stop('the zcol number is greater than the number of columns in data!')
#             } else stop("zcol should be a character or a number!")
#             #-------------
#             if (missing(n)) {
#               if (ncell(x) > 20000) n <- 99
#               else n <- 999
#               
#               cat(paste("n (number of runs in Monte Carlo simulations) is set to",n,"...\n"))
#             }
#             #----------
#             if (missing(method)) method <- 2
#             else {
#               method <- method[1]
#               if (method %in% c('boot','bootstrap','b','bo')) method <- 2
#               else if (method %in% c('perm','permutation','p','pe')) method <- 1
#               else {
#                 if (!is.numeric(method) || !method %in% 1:2) {
#                   warning('method is not identified; default ("boot") is considered')
#                   method <- 2
#                 }
#               }
#             }
#             #------
#             xx <- x
#             xx@data$elsa <- rep(NA,nrow(xx))
#             xx@data$p_value <- rep(NA,nrow(xx))
#             xx@data <- xx@data[,c('elsa','p_value')]
#             
#             x <- x@data[,zcol]
#             
#             if (is.character(x) || is.factor(x)) {
#               x <- as.character(x)
#               if (!missing(categorical) && !categorical) warning("you specified a categorical variable, so categorical changed to TRUE!")
#               categorical <- TRUE
#             }
#             
#             if (!missing(nc)) {
#               if (missing(categorical)) {
#                 if (missing(dif)) categorical <- FALSE
#                 else {
#                   categorical <- TRUE
#                   cat("input data is considered categorical, and nc is ignored!\n")
#                 }
#               } 
#             } else {
#               if (missing(categorical) && !missing(dif)) categorical <- TRUE
#             }
#             #----
#             if (missing(categorical) || !is.logical(categorical)) {
#               # guessing whether the layer is categorical:
#               if (.is.categorical(x)) {
#                 categorical <- TRUE
#                 cat("the specified variable is considered as categorical...\n")
#               } else {
#                 categorical <- FALSE
#                 cat("the specified variable is considered continuous...\n")
#               }
#             }
#             #----
#             if (!categorical && missing(nc)) {
#               nc <- nclass(x)
#               classes <- 1:nc
#             } else if (categorical) {
#               classes <- unique(x)
#               nc <- length(classes)
#             }
#             #-----
#             if (categorical) {
#               if (missing(dif)) {
#                 dif <- rep(1,nc*nc)
#                 for (i in 1:nc) dif[(i-1)*nc+i] <-0
#               } else {
#                 dif <- .checkDif(dif,classes)
#               }
#             }
#             #-----
#             
#             if (!categorical) x <- categorize(x,nc)
#             
#             
#             
#             if (missing(null)) {
#               null <- rep(NA,length(x))
#               null[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE)
#             } else if (class(null) %in% c('numeric','integer')) {
#               if (!length(null) == length(x)) {
#                 warning('null is a numeric vector with a different length; so, the null is generated given the default settings!')
#                 null <- rep(NA,length(x))
#                 null[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE)
#               }
#             } else {
#               warning('null is not a numeric vector, so, the null is generated given the default settings!')
#               null <- rep(NA,length(x))
#               null[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE)
#             }
#             .Call('elsac_vector_test',as.integer(lc[]), as.vector(rc.n[]), z,  as.integer(nc), as.integer(classes), dif, as.integer(2),as.integer(99))
#             if (categorical) {
#               xx@data$elsa <- .Call('elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif)
#               xx@data$p_value <- .Call('elsac_vector_test', as.integer(x),as.integer(null), d,  as.integer(nc), as.integer(classes), dif, as.integer(method),as.integer(n))
#             } else {
#               xx@data$elsa <- .Call('elsa_vector', as.integer(x), d, as.integer(nc))
#               xx@data$p_value <- .Call('elsa_vector_test', as.integer(x),as.integer(null), d,  as.integer(nc), as.integer(method),as.integer(n))
#             }
#             
#           }
# )
# #---------------

setMethod('elsa.test', signature(x='Spatial'), 
          function(x, d, n, method, null, nc, categorical, dif,zcol,longlat,...) {
            
            if (missing(d)) stop('d is missed!')
            else if (!class(d) %in% c('numeric','integer','neighbours')) stop('d should be either a number (distance) or an object of class neighbours (created by dneigh function')
            
            if (missing(longlat) || !is.logical(longlat)) longlat <- NULL
            
            if (class(x) == 'SpatialPolygonsDataFrame') {
              if (!inherits(d,'neighbours')) d <- dneigh(x, 0, d[1],longlat=longlat,method = 'centroid')
            } else if (class(x) == 'SpatialPointsDataFrame') {
              if (!inherits(d,'neighbours')) d <- dneigh(x, 0, d[1],longlat=longlat)
            } else stop('x should be a SpatialPointsDataFrame or SpatialPolygonsDataFrame!')
            
            d <- d@neighbours
            #-------
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
            #-------------
            if (missing(n)) {
              if (nrow(x) > 10000) n <- 99
              else n <- 999
              
              cat(paste("n (number of runs in Monte Carlo simulations) is set to",n,"...\n"))
            }
            #----------
            if (missing(method)) method <- 2
            else {
              method <- method[1]
              if (method %in% c('boot','bootstrap','b','bo')) method <- 2
              else if (method %in% c('perm','permutation','p','pe')) method <- 1
              else {
                if (!is.numeric(method) || !method %in% 1:2) {
                  warning('method is not identified; default ("boot") is considered')
                  method <- 2
                }
              }
            }
            #------
            xx <- x
            xx@data$elsa <- rep(NA,nrow(xx))
            xx@data$p_value <- rep(NA,nrow(xx))
            xx@data <- xx@data[,c('elsa','p_value')]
            
            x <- x@data[,zcol]
            
            if (is.character(x) || is.factor(x)) {
              x <- as.character(x)
              if (!missing(categorical) && !categorical) warning("you specified a categorical variable, so categorical changed to TRUE!")
              categorical <- TRUE
            }
            
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
                cat("the specified variable is considered as categorical...\n")
              } else {
                categorical <- FALSE
                cat("the specified variable is considered continuous...\n")
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
            if (missing(null)) {
              null <- rep(NA,length(x))
              null[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE)
            } else if (class(null) %in% c('numeric','integer')) {
              if (!length(null) == length(x)) {
                warning('null is a numeric vector with a different length; so, the null is generated given the default settings!')
                null <- rep(NA,length(x))
                null[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE)
              }
            } else {
              warning('null is not a numeric vector, so, the null is generated given the default settings!')
              null <- rep(NA,length(x))
              null[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE)
            }
            
            if (categorical) {
              xx@data$elsa <- .Call('v_elsac_vector', x, d, as.integer(nc), as.integer(classes),dif, PACKAGE='elsa')
              xx@data$p_value <- .Call('elsac_vector_test', x,null, d,  as.integer(nc), as.integer(classes), dif, as.integer(method),as.integer(n), PACKAGE='elsa')
            } else {
              xx@data$elsa <- .Call('v_elsa_vector', x, d, as.integer(nc), PACKAGE='elsa')
              xx@data$p_value <- .Call('elsa_vector_test', x,null, d,  as.integer(nc), as.integer(method),as.integer(n), PACKAGE='elsa')
            }
            
          }
)


#--------------
._elsa.testR <- function(x, d, n=99, nc, categorical, dif,cells,filename,...) {
  if (missing(filename)) filename <- ''
  if (missing(n)) n <- 99
  
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
  if (!categorical) {
    if (missing(nc)) nc <- nclass(x)
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
  
  nNA <- which(!is.na(x[]))
  null <- calc(x,function(x) { x[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE); x})
  null <- null[nNA]
  if (missing(cells)) {
    e1 <- elsa(x,d=d,nc=nc,categorical=categorical,dif=dif)
    o2 <- x
    #o2 <- calc(o2,function(x) { x[!is.na(x)] <- 0; x})
    o2[nNA] <- 0
    o1 <- raster(x)
    for (i in 1:n) {
      #o1 <- calc(x,function(x) { x[!is.na(x)] <- sample(null,length(x[!is.na(x)]),replace=TRUE); x})
      o1[nNA] <- sample(null,length(null),replace=TRUE)
      e2 <- elsa(o1,d=d,nc=nc,categorical=categorical,dif=dif)
      ee <- e1 - e2
      ee <- calc(ee,function(x) {x[x > 0] <- 1; x[x <= 0] = 0; x})
      o2 <- o2 + ee
    }
    rm(e1,e2,ee,o1)
    o2 <- (o2+1) / (n+1)
    
    filename <- trim(filename)
    if (filename != '') writeRaster(o2,filename,...)
  } else {
    e1 <- elsa(x,d=d,nc=nc,categorical=categorical,dif=dif,cells=cells)
    o2 <- rep(0,length(cells))
    for (i in 1:n) {
      o1[nNA] <- sample(null,length(null),replace=TRUE)
      #o1 <- calc(x,function(x) { x[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE); x})
      e2 <- elsa(o1,d=d,nc=nc,categorical=categorical,dif=dif,cells=cells)
      ee <- e1 - e2
      ee <- ifelse(ee > 0,1,0)
      o2 <- o2 + ee
    }
    rm(e1,e2,ee,o1)
    o2 <- (o2+1) / (n+1)
  }
  o2
}
