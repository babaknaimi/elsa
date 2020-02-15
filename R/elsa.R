# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2016
# last update: November 2022
# Version 3.3
# Licence GPL v3 


.Filter<-function(r,d1=0,d2) {
  c<- d2%/%r
  x<- y<- seq(-c,c,1)
  eg<- expand.grid(x,y)
  eg[1]<- -eg[1]
  eg[3]<- sqrt(eg[1]^2+eg[2]^2)
  ndim<- c*2+1
  m<-matrix(eg[,3],ncol=ndim,nrow=ndim)*r
  mw<- which(m > d1 & m <= d2)
  m[mw] = 1 
  m[-mw] = NA
  mw <- trunc(length(m)/2)+1
  m[mw]<- 1
  
  w <- ncol(m)
  mr <- matrix(ncol=w,nrow=w) 
  tr <- trunc(w/2)
  
  for (i in 1:w) {
    mr[i,] <-tr 
    tr <- tr - 1
  }
  
  rc = cbind(r=as.vector(m * mr),c=as.vector(m * t(mr)))
  rc=rc[apply(rc,1,function(x) (all(!is.na(x)))),]
  
  return (list(w,rc))
}

#----------
.is.categorical <- function(x) {
  (length(unique(x - round(x))) == 1) && length(unique(x)) < 50
}
#----------
.is.categoricalRaster <- function(x) {
  if (ncell(x) > 1.3e6) {
    x <- sampleRandom(x,1e6)
  }
  (length(unique(x - round(x))) == 1) && length(unique(x)) < 50
}
#----------
.is.categoricalNumeric <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 1e6) {
    x <- x[sample(length(x),1e6)]
  }
  (length(unique(x - round(x))) == 1) && length(unique(x)) < 50
}
#----------
.is.categoricalSpatRaster <- function(x) {
  if (nlyr(x) > 1) {
    x <- x[[1]]
  }
  if (ncell(x) > 1.3e6) {
    x <- spatSample(x,1e6,na.rm=TRUE)[,1]
  }
  (length(unique(x - round(x))[,1]) == 1) && length(unique(x)[,1]) < 50
}


.checkDif <- function(dif,classes) {
  classes <- classes[!is.na(classes)]
  nc <- length(classes)
  classes <- as.character(classes)
  if(is.list(dif)) {
    if(!is.null(names(dif))) {
      if (!all(classes %in% names(dif))) stop("categories in the dif argument do not match with the categories in the input layer!")
      dif <- dif[classes]
    } else {
      if (length(dif) == nc) names(dif) <- classes
      else stop("dif content does not match with the categories in the input layer!")
    }
    
    if (!all(unlist(lapply(classes,function(x) {all(classes %in% names(dif[[x]]))})))) stop("Each item in the dif list should be a named vector with all the classes to which the dissimilarities are specified. It seems that something is wrong with the structure (either some classes are missed in the items, or the vectors are not named appropriately)!")
    
    for (i in 1:length(dif)) dif[[i]] <- dif[[i]][classes]
    
    #if (!all(unlist(lapply(classes,function(x) {length(dif[[x]]) == nc})))) stop("the length of categories' differences in the dif argument does not match with the number of categories!")
    dif <- as.vector(as.matrix(as.data.frame(dif)))
  } else if (is.matrix(dif) || is.data.frame(dif)) {
    #if (ncol(dif) != nrow(dif) || ncol != nc) stop("number of rows and columns in dif should be the same, equal to the number of categories!")
    if (ncol(dif) != nrow(dif)) stop("number of rows and columns in dif should be the same!")
    if (colnames(dif) %in% row.names(dif)) stop("names of columns and rows should be the same, corresponding to classes!")
    if (!all(classes %in% colnames(dif))) stop("Some classes do not exist in dif data.frame/matrix!")
    dif <- dif[classes,classes]
    dif <- as.vector(as.matrix(as.data.frame(dif)))
    
  } else if (is.vector(dif)) {
    if (length(dif) != nc*nc) stop("dif argument does not have an appropriate structure!")
  } else stop("dif argument does not have an appropriate structure!")
  return(dif)
}
########################################
########################################
########################################
if (!isGeneric("elsa")) {
  setGeneric("elsa", function(x,d,nc,categorical,dif,classes,stat,...)
    standardGeneric("elsa"))
}


setMethod('elsa', signature(x='RasterLayer'), 
          function(x,d,nc,categorical,dif,classes,stat,cells,filename,verbose=TRUE,...) {
            if (missing(classes)) classes <- NULL
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(stat) || is.null(stat)) stat <- 'elsa'
            else {
              stat <- tolower(stat)
              if (length(stat) == 1) {
                if (!stat %in% c('elsa','ec','ea')) {
                  stat <- 'elsa'
                  warning('stat should be either of "ELSA", "Ec", "Ea"; the default "ELSA" is considered!')
                }
              } else {
                if (!all(tolower(stat) %in% c('elsa','ec','ea'))) stop('stat should be selected from "ELSA", "Ea", "Ec"')
              }
            }
            #----
            if (missing(d)) d <- res(x)[1] * sqrt(2)
            
            if (missing(filename)) filename <- ''
            
            if (!missing(nc) && !is.null(nc) && !is.na(nc)) {
              if (missing(categorical)) {
                if (missing(dif) && is.null(classes)) categorical <- FALSE
                else {
                  if (!missing(dif) && !is.null(dif) && !is.na(dif) && !is.null(classes) && !is.na(classes) && .is.categoricalRaster(x)) categorical <- TRUE
                  else {
                    if (verbose) cat("the input data seems continues (if not, use categorical=TRUE)!.... dif/classes is ignored!\n")
                  } 
                }
              } 
            } else {
              if (missing(categorical) && !missing(dif) && !is.null(dif) && !is.na(dif) && !is.null(classes) && !is.na(classes)) categorical <- TRUE
            }
            #----
            if (missing(categorical) || !is.logical(categorical)) {
              # guessing whether the layer is categorical:
              if (.is.categoricalRaster(x)) {
                categorical <- TRUE
                if (verbose) cat("the input is considered as a categorical variable...\n")
              } else {
                categorical <- FALSE
                if (verbose) cat("the input is considered as a continuous variable...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
            } else if (categorical) {
              if (is.null(classes) || is.na(classes)) {
                if (missing(dif) || is.null(classes) || is.na(classes) ) {
                  classes <- unique(x)
                } else {
                  if (length(names(dif)) > 1) {
                    classes <- names(dif)
                    .ux <- as.character(unique(x))
                    if (!all(.ux %in% classes)) classes <- .ux
                  } else classes <- unique(x)
                }
              } else {
                .ux <- unique(x)
                if (is.character(classes)) .ux <- as.character(.ux)
                if (!all(.ux %in% classes)) stop('the specified "classes" does not cover all or some of values in the input raster!')
              }
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
                  if (length(stat) == 1 && stat == 'elsa') {
                    out[] <- .Call('v_elsac', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                    names(out) <- 'ELSA'
                  } else {
                  xx <- .Call('elsac', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                  if (length(stat) > 1) {
                    nnn <- c()
                    if ('ea' %in% stat) {
                      outx <- raster(out)
                      outx[] <- xx[[2]]
                      out <- addLayer(out,outx)
                      nnn <- c(nnn,'Ea')
                    }
                    if ('ec' %in% stat) {
                      outx <- raster(out)
                      outx[] <- xx[[1]]
                      out <- addLayer(out,outx)
                      nnn <- c(nnn,'Ec')
                    }
                    if ('elsa' %in% stat) {
                      outx <- raster(out)
                      outx[] <- xx[[2]] * xx[[1]]
                      out <- addLayer(out,outx)
                      nnn <- c(nnn,'ELSA')
                    }
                    names(out) <- nnn
                    
                  } else {
                    if (stat == 'ea') {
                      out[] <- xx[[2]]
                      names(out) <- 'Ea'
                    } else {
                      out[] <- xx[[1]]
                      names(out) <- 'Ec'
                    }
                   }
                  }
                  if (filename != '') out <- writeRaster(out, filename, ...)
                } else {
                  if (length(stat) == 1) {
                    if (stat == 'elsa') out <- .Call('v_elsac_cell', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                    else if (stat == 'ec') out <- .Call('v_elsac_cell_Ec', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                    else if (stat == 'ea') out <- .Call('v_elsac_cell_Ea', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                  } else {
                    xx <- .Call('elsac_cell', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                    out <- list()
                    if ('ea' %in% stat) {
                      out[['Ea']] <- xx[[2]]
                    }
                    if ('ec' %in% stat) {
                      out[['Ec']] <- xx[[1]]
                    }
                    if ('elsa' %in% stat) {
                      out[['ELSA']] <-  xx[[2]] * xx[[1]]
                    }
                  }
                }
              } else {
                if (missing(cells)) {
                  
                  if (length(stat) == 1 && stat == 'elsa') {
                    out[] <- .Call('v_elsa', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                    names(out) <- 'ELSA'
                  } else {
                    xx <- .Call('elsa', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                    if (length(stat) > 1) {
                      nnn <- c()
                      if ('ea' %in% stat) {
                        outx <- raster(out)
                        outx[] <- xx[[2]]
                        out <- addLayer(out,outx)
                        nnn <- c(nnn,'Ea')
                      }
                      if ('ec' %in% stat) {
                        outx <- raster(out)
                        outx[] <- xx[[1]]
                        out <- addLayer(out,outx)
                        nnn <- c(nnn,'Ec')
                      }
                      if ('elsa' %in% stat) {
                        outx <- raster(out)
                        outx[] <- xx[[2]] * xx[[1]]
                        out <- addLayer(out,outx)
                        nnn <- c(nnn,'ELSA')
                      }
                      names(out) <- nnn
                      
                    } else {
                      if (stat == 'ea') {
                        out[] <- xx[[2]]
                        names(out) <- 'Ea'
                      } else {
                        out[] <- xx[[1]]
                        names(out) <- 'Ec'
                      }
                    }
                  }
                  if (filename != '') out <- writeRaster(out, filename, ...)
                  
                } else {
                  if (length(stat) == 1) {
                    if (stat == 'elsa') out <- .Call('v_elsa_cell', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
                    else if (stat == 'ec') out <- .Call('v_elsa_cell_Ec', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
                    else if (stat == 'ea') out <- .Call('v_elsa_cell_Ea', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
                  } else {
                    xx <- .Call('elsa_cell', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
                    out <- list()
                    if ('ea' %in% stat) {
                      out[['Ea']] <- xx[[2]]
                    }
                    if ('ec' %in% stat) {
                      out[['Ec']] <- xx[[1]]
                    }
                    if ('elsa' %in% stat) {
                      out[['ELSA']] <-  xx[[2]] * xx[[1]]
                    }
                  }
                }
              }
            } else {
              if (verbose) cat("\nThe input dataset is considered as a big raster dataset that will be handled out of memory (on the disk), but if you have enough memory on your machine, you can change the settings for maxmemory, and chuncksize, in the rasterOptions function, then the process may be handled in memory that would be much faster...")
              
              tr <- blockSize(out, minblocks=3, minrows=fdim)
              pb <- pbCreate(tr$n, label='ELSA',...)
              addr <- floor(fdim / 2)
              
              if (length(stat) > 1) warning(paste('for big rasters, stat can only have one value, so stat = "',toupper(stat[1]),'", is considered!\n',sep=''))
              stat <- stat[1]
              
              if (missing(cells)) {
                out <- writeStart(out, filename)
                v <- getValues(x, row=1, nrows=tr$nrows[1]+addr)
                if (!categorical) {
                  v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                } else {
                  v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                }
                
                if (stat == 'elsa') v <- v[[1]] * v[[2]]
                else if (stat == 'ea') v <- v[[2]]
                else v <- v[[1]]
                
                ex <- length(v) - (addr * ncl)
                out <- writeValues(out, v[1:ex], 1)
                
                for (i in 2:(tr$n-1)) {
                  v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr))
                  if (!categorical) {
                    v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                  } else {
                    v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                  }
                  
                  if (stat == 'elsa') v <- v[[1]] * v[[2]]
                  else if (stat == 'ea') v <- v[[2]]
                  else v <- v[[1]]
                  
                  st <- (addr * ncl)+1
                  ex <- length(v) - (addr * ncl)
                  out <- writeValues(out, v[st:ex], tr$row[i])
                  pbStep(pb)
                }
                
                i <- tr$n
                v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr)
                if (!categorical) {
                  v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                } else {
                  v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                }
                
                if (stat == 'elsa') v <- v[[1]] * v[[2]]
                else if (stat == 'ea') v <- v[[2]]
                else v <- v[[1]]
                
                st <- (addr * ncl)+1
                ex <- length(v)
                out <- writeValues(out, v[st:ex], tr$row[i])
                pbStep(pb)
                out <- writeStop(out)      
                pbClose(pb)  
              } else {
                v <- getValues(x, row=1, nrows=tr$nrows[1]+addr)
                cls <- cells[which(cells <= (tr$nrows[1]) * ncl)]
                if (length(cls) > 0) {
                  if (!categorical) {
                    v <- .Call('elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls), PACKAGE='elsa')
                  } else {
                    v <- .Call('elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls), PACKAGE='elsa')
                  }
                  
                  if (length(stat) > 1) {
                    out <- list()
                    if ('ea' %in% stat) {
                      out[['Ea']] <- c(out[['L']],v[[2]])
                    }
                    if ('ec' %in% stat) {
                      out[['Ec']] <- c(out[['R']],v[[1]])
                    }
                    if ('elsa' %in% stat) {
                      out[['ELSA']] <-  c(out[['ELSA']],v[[2]] * v[[1]])
                    }
                  } else {
                    out <- c()
                    if (stat == 'ea') {
                      out <- c(out, v[[2]])
                    } else if (stat == 'ec') {
                      out <- c(out, v[[1]])
                    } else out <- c(out, v[[1]]*v[[2]])
                  }
                }
                
                for (i in 2:(tr$n-1)) {
                  cls <- cells[which(cells > ((tr$nrow[i] - 1) * ncl) & cells <= ((tr$nrows[i]+ tr$nrows[i] - 1) * ncl))]
                  if (length(cls) > 0) {
                    cls <- cls - ((tr$row[i]-addr-1)*ncl)
                    v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr))
                    if (!categorical) {
                      v <- .Call('elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls), PACKAGE='elsa')
                    } else {
                      v <- .Call('elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls), PACKAGE='elsa')
                    }
                    if (length(stat) > 1) {
                      if ('ea' %in% stat) {
                        out[['Ea']] <- c(out[['L']],v[[2]])
                      }
                      if ('ec' %in% stat) {
                        out[['Ec']] <- c(out[['R']],v[[1]])
                      }
                      if ('elsa' %in% stat) {
                        out[['ELSA']] <-  c(out[['ELSA']],v[[2]] * v[[1]])
                      }
                    } else {
                      if (stat == 'ea') {
                        out <- c(out, v[[2]])
                      } else if (stat == 'ec') {
                        out <- c(out, v[[1]])
                      } else out <- c(out, v[[1]]*v[[2]])
                    }
                    pbStep(pb)
                  }
                }
                
                i <- tr$n
                cls <- cells[which(cells > ((tr$nrow[i] - 1) * ncl) & cells <= ((tr$nrows[i]+ tr$nrows[i] - 1) * ncl))]
                cls <- cls - ((tr$row[i]-addr-1)*ncl)
                if (length(cls) > 0) {
                  v <- getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr)
                  if (!categorical) {
                    v <- .Call('v_elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls), PACKAGE='elsa')
                  } else {
                    v <- .Call('v_elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls), PACKAGE='elsa')
                  }
                  if (length(stat) > 1) {
                    if ('ea' %in% stat) {
                      out[['Ea']] <- c(out[['Ea']],v[[2]])
                    }
                    if ('ec' %in% stat) {
                      out[['Ec']] <- c(out[['Ec']],v[[1]])
                    }
                    if ('elsa' %in% stat) {
                      out[['ELSA']] <-  c(out[['ELSA']],v[[2]] * v[[1]])
                    }
                  } else {
                    if (stat == 'ea') {
                      out <- c(out, v[[2]])
                    } else if (stat == 'ec') {
                      out <- c(out, v[[1]])
                    } else out <- c(out, v[[1]]*v[[2]])
                  }
                  pbStep(pb)
                  pbClose(pb)  
                }
              }
            }
            return(out)
          }
)  

#---------------

setMethod('elsa', signature(x='SpatialPointsDataFrame'), 
          function(x,d,nc,categorical,dif,classes,stat,zcol,drop,verbose=TRUE,...) {
            if (missing(classes)) classes <- NULL
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(d)) stop('d is missed!')
            else if (!class(d) %in% c('numeric','integer','neighbours')) stop('d should be either a number (distance) or an object of class neighbours (created by dneigh function')
            
            if (!inherits(d,'neighbours')) d <- dneigh(x, 0, d[1])
            d <- d@neighbours
            
            if (missing(drop) || !is.logical(drop[1])) drop <- FALSE
            else drop <- drop[1]
            
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
                  if (verbose) cat("input data is considered categorical, and nc is ignored!\n")
                }
              } 
            } else {
              if (missing(categorical) && !missing(dif) && !is.null(classes)) categorical <- TRUE
            }
            #----
            if (missing(categorical) || !is.logical(categorical)) {
              # guessing whether the layer is categorical:
              if (.is.categorical(x)) {
                categorical <- TRUE
                if (verbose) cat("the specified variable is considered as categorical...\n")
              } else {
                categorical <- FALSE
                if (verbose) cat("the specified variable is considered continuous...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
            } else if (categorical) {
              if (is.null(classes) || is.na(classes)) {
                if (missing(dif) || is.null(classes) || is.na(classes) ) {
                  classes <- unique(x)
                } else {
                  if (length(names(dif)) > 1) {
                    classes <- names(dif)
                    .ux <- as.character(unique(x))
                    if (!all(.ux %in% classes)) classes <- .ux
                  } else classes <- unique(x)
                }
              } else {
                .ux <- unique(x)
                if (is.character(classes)) .ux <- as.character(.ux)
                if (!all(.ux %in% classes)) stop('the specified "classes" does not cover all or some of values in the input raster!')
              }
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
            
            
            if (categorical) {
              x <- .Call('elsac_vector', x, d, as.integer(nc), as.integer(classes),dif, PACKAGE='elsa')
            } else {
              x <- .Call('elsa_vector', x, d, as.integer(nc), PACKAGE='elsa')
            }
            
            xx@data$Ea <- x[[2]]
            xx@data$Ec <- x[[1]]
            xx@data$ELSA <- x[[1]] * x[[2]]
            xx@data <- xx@data[,c('Ea','Ec','ELSA')]
            
            if (!drop) xx
            else xx@data
            
          }
)  


setMethod('elsa', signature(x='SpatialPolygonsDataFrame'), 
          function(x,d,nc,categorical,dif,classes,stat,zcol,drop,method,verbose=TRUE,...) {
            if (missing(classes)) classes <- NULL
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(d)) stop('d is missed!')
            else if (!class(d) %in% c('numeric','integer','neighbours')) stop('d should be either a number (distance) or an object of class neighbours (created by dneigh function')
            
            if (missing(method)) method <- 'centroid'
            
            if (!inherits(d,'neighbours')) d <- dneigh(x, 0, d[1],method = method)
            d <- d@neighbours
            if (missing(drop) || !is.logical(drop[1])) drop <- FALSE
            else drop <- drop[1]
            
            
            
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
                  if (verbose) cat("input data is considered categorical, and nc is ignored!\n")
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
                if (verbose) cat("the specified variable is considered as categorical...\n")
              } else {
                categorical <- FALSE
                if (verbose) cat("the specified variable is considered continuous...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
            } else if (categorical) {
              if (is.null(classes) || is.na(classes)) {
                if (missing(dif) || is.null(classes) || is.na(classes) ) {
                  classes <- unique(x)
                } else {
                  if (length(names(dif)) > 1) {
                    classes <- names(dif)
                    .ux <- as.character(unique(x))
                    if (!all(.ux %in% classes)) classes <- .ux
                  } else classes <- unique(x)
                }
              } else {
                .ux <- unique(x)
                if (is.character(classes)) .ux <- as.character(.ux)
                if (!all(.ux %in% classes)) stop('the specified "classes" does not cover all or some of values in the input raster!')
              }
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
            
            if (categorical) {
              x <- .Call('elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif, PACKAGE='elsa')
            } else {
              x <-.Call('elsa_vector', as.integer(x), d, as.integer(nc), PACKAGE='elsa')
            }
            
            xx@data$Ea <- x[[2]]
            xx@data$Ec <- x[[1]]
            xx@data$ELSA <- x[[1]] * x[[2]]
            xx@data <- xx@data[,c('Ea','Ec','ELSA')]
            
            if (!drop) xx
            else xx@data
            
          }
)  
###########################
###########################
#------------ Adding functions to support input from the terra package (SpatRaster and SpatVector) ----
###########################
###########################


setMethod('elsa', signature(x='SpatRaster'), 
          function(x,d,nc,categorical,dif,classes,stat,cells,filename,verbose=TRUE,...) {
            if (missing(classes)) classes <- NULL
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(stat) || is.null(stat)) stat <- 'elsa'
            else {
              stat <- tolower(stat)
              if (length(stat) == 1) {
                if (!stat %in% c('elsa','ec','ea')) {
                  stat <- 'elsa'
                  warning('stat should be either of "ELSA", "Ec", "Ea"; the default "ELSA" is considered!')
                }
              } else {
                if (!all(tolower(stat) %in% c('elsa','ec','ea'))) stop('stat should be selected from "ELSA", "Ea", "Ec"')
              }
            }
            #----
            if (missing(d)) d <- res(x)[1] * sqrt(2)
            
            if (missing(filename)) filename <- ''
            
            if (!missing(nc) && !is.null(nc) && !is.na(nc)) {
              if (missing(categorical)) {
                if (missing(dif) && is.null(classes)) categorical <- FALSE
                else {
                  if (!missing(dif) && !is.null(dif) && !is.na(dif) && !is.null(classes) && !is.na(classes) && .is.categoricalSpatRaster(x)) categorical <- TRUE
                  else {
                    if (verbose) cat("the input data seems continues (if not, use categorical=TRUE)!.... dif/classes is ignored!\n")
                  } 
                }
              } 
            } else {
              if (missing(categorical) && !missing(dif) && !is.null(dif) && !is.na(dif) && !is.null(classes) && !is.na(classes)) categorical <- TRUE
            }
            #----
            if (missing(categorical) || !is.logical(categorical)) {
              # guessing whether the layer is categorical:
              if (.is.categoricalSpatRaster(x)) {
                categorical <- TRUE
                if (verbose) cat("the input is considered as a categorical variable...\n")
              } else {
                categorical <- FALSE
                if (verbose) cat("the input is considered as a continuous variable...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x[[1]])
            } else if (categorical) {
              if (is.null(classes) || is.na(classes)) {
                if (missing(dif) || is.null(classes) || is.na(classes) ) {
                  classes <- unique(x[[1]],incomparables = TRUE)[[1]]
                  if (nlyr(x) > 1) warning('since multiple categorical layers are in the SpatRaster object and "classes" is not specified, the classes are extracted from the first layer!')
                } else {
                  if (length(names(dif)) > 1) {
                    classes <- names(dif)
                    .ux <- as.character(unique(x[[1]],incomparables = TRUE)[[1]])
                    #.ux <- lapply(.ux,as.character)
                    # if (!all(sapply(.ux,function(x) all(x %in% classes)))) {
                    #   if (ncol(.ux) > 1) classes <- .ux
                    #   else classes <- .ux[[1]]
                    #   #if (any(sapply(.ux,function(x) length(unique(x))) / length(unique(unlist(.ux))) > 0.5)) stop('It seems that the classes in different layers ')
                    # }
                    if (!all(.ux %in% classes)) classes <- .ux
                  } else {
                    #.ux <- unique(x[[1]],incomparables = TRUE)[[1]]
                    classes <- unique(x[[1]],incomparables = TRUE)[[1]]
                    if (nlyr(x) > 1) warning('since multiple categorical layers are in the SpatRaster and "classes" is not specified, the classes are extracted from the first layer!')
                  }
                }
              } else {
                .ux <- unique(x,incomparables = TRUE)
                if (is.character(classes)) .ux <- lapply(.ux,as.character)
                # if (is.list(classes) && length(classes) == length(.ux)) {
                #   for (i in 1:length(classes)) {
                #     if (!all(sapply(.ux,function(x) all(x %in% classes[[i]])))) stop('the specified "classes" does not cover all or some of values in the input raster!')
                #   }
                # } else {
                #   if (!all(sapply(.ux,function(x) all(x %in% classes)))) stop('the specified "classes" does not cover all or some of values in the input raster!')
                # }
                if (!all(sapply(.ux,function(x) all(x %in% classes)))) stop('the specified "classes" does not cover all or some of values in the input raster!')
              }
              # if (is.list(classes)) {
              #   nc <- sapply(classes,length)
              # } else nc <- length(classes)
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
            w <-.Filter(r=res(x)[1],d1=0,d2=d)
            fdim <- w[[1]]
            w <- w[[2]]
            
            if (fdim < 3) stop("d must be at least equal to the input raster resolution!")
            
            if (!categorical) x <- categorize(x,nc)
            
            out <- rast(x)
            ncl <- ncol(out)
            nrw <- nrow(out)
            filename=trim(filename)
            gc()
            #---------------
            if (.canProcessInMemory(out,n=3)) {
              
              if (categorical) {
                if (missing(cells)) {
                  if (length(stat) == 1 && stat == 'elsa') {
                    for (i in 1:nlyr(x)) {
                      out[[i]][] <- .Call('v_elsac', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                    }
                    names(out) <- paste0(names(x),'_ELSA')
                  } else {
                    
                    if (nlyr(x) > 1) {
                      if (length(stat) > 1) out <- rast(x[[1]])
                      for (i in 1:nlyr(x)) {
                        xx <- .Call('elsac', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                        if (length(stat) > 1) {
                          .lyrnames <- names(x)
                          nnn <- c()
                          if ('ea' %in% stat) {
                            outx <- rast(out[[1]])
                            outx[] <- xx[[2]]
                            out <- c(out,outx)
                            nnn <- c(nnn,paste0(.lyrnames[i],'_Ea'))
                          }
                          if ('ec' %in% stat) {
                            outx <- rast(out[[1]])
                            outx[] <- xx[[1]]
                            out <- c(out,outx)
                            nnn <- c(nnn,paste0(.lyrnames[i],'_Ec'))
                          }
                          if ('elsa' %in% stat) {
                            outx <- rast(out[[1]])
                            outx[] <- xx[[2]] * xx[[1]]
                            out <- c(out,outx)
                            nnn <- c(nnn,paste0(.lyrnames[i],'_ELSA'))
                          }
                          names(out) <- nnn
                          
                        } else {
                          if (stat == 'ea') {
                            out[[i]][] <- xx[[2]]
                            names(out) <- 'Ea'
                          } else {
                            out[[i]][] <- xx[[1]]
                            names(out) <- 'Ec'
                          }
                        }
                      }
                      
                    } else {
                      xx <- .Call('elsac', x[][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                      if (length(stat) > 1) {
                        nnn <- c()
                        if ('ea' %in% stat) {
                          outx <- rast(out)
                          outx[] <- xx[[2]]
                          out <- c(out,outx)
                          nnn <- c(nnn,'Ea')
                        }
                        if ('ec' %in% stat) {
                          outx <- rast(out)
                          outx[] <- xx[[1]]
                          out <- c(out,outx)
                          nnn <- c(nnn,'Ec')
                        }
                        if ('elsa' %in% stat) {
                          outx <- rast(out)
                          outx[] <- xx[[2]] * xx[[1]]
                          out <- c(out,outx)
                          nnn <- c(nnn,'ELSA')
                        }
                        names(out) <- nnn
                        
                      } else {
                        if (stat == 'ea') {
                          out[] <- xx[[2]]
                          names(out) <- 'Ea'
                        } else {
                          out[] <- xx[[1]]
                          names(out) <- 'Ec'
                        }
                      }
                    }
                    
                  }
                  if (filename != '') out <- writeRaster(out, filename, ...)
                } else {
                  if (nlyr(x) > 1) {
                    out <- matrix(NA,nrow=length(cells),ncol=nlyr(x)*length(stat))
                    .layernames <- names(x)
                    colnames(out) <- paste0(.layernames,'_',stat)
                    for (i in 1:nlyr(x)) {
                      if (length(stat) == 1) {
                        if (stat == 'elsa') out[,i] <- .Call('v_elsac_cell', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                        else if (stat == 'ec') out[,i] <- .Call('v_elsac_cell_Ec', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                        else if (stat == 'ea') out[,i] <- .Call('v_elsac_cell_Ea', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                      } else {
                        xx <- .Call('elsac_cell', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                        
                        if ('ea' %in% stat) {
                          out[,paste0(.layernames[i],'_','ea')] <- xx[[2]]
                        }
                        if ('ec' %in% stat) {
                          out[,paste0(.layernames[i],'_','ec')] <- xx[[1]]
                        }
                        if ('elsa' %in% stat) {
                          out[,paste0(.layernames[i],'_','elsa')] <- xx[[2]] * xx[[1]]
                        }
                      }
                    }
                  } else {
                    if (length(stat) == 1) {
                      if (stat == 'elsa') out <- .Call('v_elsac_cell', x[][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                      else if (stat == 'ec') out <- .Call('v_elsac_cell_Ec', x[][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                      else if (stat == 'ea') out <- .Call('v_elsac_cell_Ea', x[][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                    } else {
                      xx <- .Call('elsac_cell', x[][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
                      out <- list()
                      if ('ea' %in% stat) {
                        out[['Ea']] <- xx[[2]]
                      }
                      if ('ec' %in% stat) {
                        out[['Ec']] <- xx[[1]]
                      }
                      if ('elsa' %in% stat) {
                        out[['ELSA']] <-  xx[[2]] * xx[[1]]
                      }
                    }
                  }
                  
                }
              } else {
                if (missing(cells)) {
                  
                  if (length(stat) == 1 && stat == 'elsa') {
                    for (i in 1:nlyr(x)) {
                      out[[i]][] <- .Call('v_elsa', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                    }
                    
                    names(out) <- paste0(names(x),'_ELSA')
                  } else {
                    
                    out <- rast(x[[1]])
                    for (i in 1:nlyr(x)) {
                      xx <- .Call('elsa', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                      if ('ea' %in% stat) {
                        outx <- rast(x[[1]])
                        outx[] <- xx[[2]]
                        names(outx) <- paste0(names(x[[i]]),'_Ea')
                        out <- c(out,outx)
                      }
                      
                      if ('ec' %in% stat) {
                        outx <- rast(x[[1]])
                        outx[] <- xx[[1]]
                        names(outx) <- paste0(names(x[[i]]),'_Ec')
                        out <- c(out,outx)
                      }
                      
                      if ('elsa' %in% stat) {
                        outx <- rast(x[[1]])
                        outx[] <- xx[[2]] * xx[[1]]
                        names(outx) <- paste0(names(x[[i]]),'_ELSA')
                        out <- c(out,outx)
                      }
                    }
                  }
                  if (filename != '') out <- writeRaster(out, filename, ...)
                  
                } else {
                  out <- list()
                  for (i in 1:nlyr(x)) {
                    if (length(stat) == 1) {
                      if (stat == 'elsa') out[[paste0(names(x[[i]]),'_ELSA')]] <- .Call('v_elsa_cell', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
                      else if (stat == 'ec') out[[paste0(names(x[[i]]),'Ec')]] <- .Call('v_elsa_cell_Ec', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
                      else if (stat == 'ea') out[[paste0(names(x[[i]]),'_Ea')]] <- .Call('v_elsa_cell_Ea', x[[i]][][,1], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
                    } else {
                      xx <- .Call('elsa_cell', x[[i]][][,i], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
                      if ('ea' %in% stat) {
                        out[[paste0(names(x[[i]]),'_Ea')]] <- xx[[2]]
                      }
                      if ('ec' %in% stat) {
                        out[[paste0(names(x[[i]]),'_Ec')]] <- xx[[1]]
                      }
                      if ('elsa' %in% stat) {
                        out[[paste0(names(x[[i]]),'_ELSA')]] <-  xx[[2]] * xx[[1]]
                      }
                    }
                  }
                }
              }
            } else {
              if (verbose) cat("\nThe input dataset is considered as a big raster dataset that will be handled out of memory (on the disk)...")
              
              if (nlyr(x) > 1) {
                warning("Since the raster dataset cannot handled in memory, the function is applied only to the first layer!")
                x <- x[[1]]
                out <- rast(x[[1]])
              }
              tr <- blocks(out,n=3)
              
              addr <- floor(fdim / 2)
              
              if (missing(cells)) {
                
                if (length(stat) > 1) warning(paste('for big rasters, stat can only have one value, so stat = "',toupper(stat[1]),'", is considered!\n',sep=''))
                stat <- stat[1]
                
                
                readStart(x)
                b <- writeStart(out, filename=filename,...)
                v <- readValues(x, row=1, nrows=b$nrows[1]+addr)
                if (!categorical) {
                  v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                } else {
                  v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                }
                
                if (stat == 'elsa') v <- v[[1]] * v[[2]]
                else if (stat == 'ea') v <- v[[2]]
                else v <- v[[1]]
                
                ex <- length(v) - (addr * ncl)
                writeValues(out, v[1:ex], 1, nrows=b$nrows[1])
                
                for (i in 2:(b$n-1)) {
                  v <- readValues(x, row=tr$row[i]-addr, nrows=b$nrows[i]+(2*addr))
                  if (!categorical) {
                    v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                  } else {
                    v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                  }
                  
                  if (stat == 'elsa') v <- v[[1]] * v[[2]]
                  else if (stat == 'ea') v <- v[[2]]
                  else v <- v[[1]]
                  
                  st <- (addr * ncl)+1
                  ex <- length(v) - (addr * ncl)
                  writeValues(out, v[st:ex], b$row[i],nrows=b$nrows[i])
                }
                
                i <- b$n
                v <- readValues(x, row=b$row[i]-addr, nrows=b$nrows[i])
                if (!categorical) {
                  v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                } else {
                  v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
                }
                
                if (stat == 'elsa') v <- v[[1]] * v[[2]]
                else if (stat == 'ea') v <- v[[2]]
                else v <- v[[1]]
                
                st <- (addr * ncl)+1
                ex <- length(v)
                writeValues(out, v[st:ex], tr$row[i])
                
                writeStop(out)      
                readStop(x)
              } else {
                readStart(x)
                v <- readValues(x, row=1, nrows=tr$nrows[1]+addr)
                cls <- cells[which(cells <= (tr$nrows[1]) * ncl)]
                if (length(cls) > 0) {
                  if (!categorical) {
                    v <- .Call('elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls), PACKAGE='elsa')
                  } else {
                    v <- .Call('elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls), PACKAGE='elsa')
                  }
                  
                  if (length(stat) > 1) {
                    out <- list()
                    if ('ea' %in% stat) {
                      out[['Ea']] <- c(out[['L']],v[[2]])
                    }
                    if ('ec' %in% stat) {
                      out[['Ec']] <- c(out[['R']],v[[1]])
                    }
                    if ('elsa' %in% stat) {
                      out[['ELSA']] <-  c(out[['ELSA']],v[[2]] * v[[1]])
                    }
                  } else {
                    out <- c()
                    if (stat == 'ea') {
                      out <- c(out, v[[2]])
                    } else if (stat == 'ec') {
                      out <- c(out, v[[1]])
                    } else out <- c(out, v[[1]]*v[[2]])
                  }
                }
                
                for (i in 2:(tr$n-1)) {
                  
                  cls <- cells[which((cells > ((tr$row[i] - 1) * ncl)) & (cells <= ((tr$row[i]+ tr$nrows[i] - 1) * ncl)))]
                  if (length(cls) > 0) {
                    cls <- cls - ((tr$row[i]-addr-1)*ncl)
                    v <- readValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr))
                    if (!categorical) {
                      v <- .Call('elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls), PACKAGE='elsa')
                    } else {
                      v <- .Call('elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls), PACKAGE='elsa')
                    }
                    if (length(stat) > 1) {
                      if ('ea' %in% stat) {
                        out[['Ea']] <- c(out[['L']],v[[2]])
                      }
                      if ('ec' %in% stat) {
                        out[['Ec']] <- c(out[['R']],v[[1]])
                      }
                      if ('elsa' %in% stat) {
                        out[['ELSA']] <-  c(out[['ELSA']],v[[2]] * v[[1]])
                      }
                    } else {
                      if (stat == 'ea') {
                        out <- c(out, v[[2]])
                      } else if (stat == 'ec') {
                        out <- c(out, v[[1]])
                      } else out <- c(out, v[[1]]*v[[2]])
                    }
                  }
                }
                
                i <- tr$n
                cls <- cells[which(cells > ((tr$row[i] - 1) * ncl) & cells <= ((tr$row[i]+ tr$nrows[i] - 1) * ncl))]
                cls <- cls - ((tr$row[i]-addr-1)*ncl)
                if (length(cls) > 0) {
                  v <- readValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i])
                  if (!categorical) {
                    v <- .Call('v_elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls), PACKAGE='elsa')
                  } else {
                    v <- .Call('v_elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls), PACKAGE='elsa')
                  }
                  if (length(stat) > 1) {
                    if ('ea' %in% stat) {
                      out[['Ea']] <- c(out[['Ea']],v[[2]])
                    }
                    if ('ec' %in% stat) {
                      out[['Ec']] <- c(out[['Ec']],v[[1]])
                    }
                    if ('elsa' %in% stat) {
                      out[['ELSA']] <-  c(out[['ELSA']],v[[2]] * v[[1]])
                    }
                  } else {
                    if (stat == 'ea') {
                      out <- c(out, v[[2]])
                    } else if (stat == 'ec') {
                      out <- c(out, v[[1]])
                    } else out <- c(out, v[[1]]*v[[2]])
                  }
                  
                }
                readStop(x)
              }
            }
            return(out)
          }
)  


#---------------

setMethod('elsa', signature(x='SpatVector'), 
          function(x,d,nc,categorical,dif,classes,stat,zcol,drop,verbose=TRUE,...) {
            if (missing(classes)) classes <- NULL
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(d)) stop('d is missed!')
            else if (!class(d) %in% c('numeric','integer','neighbours')) stop('d should be either a number (distance) or an object of class neighbours (created by dneigh function')
            
            if (!inherits(d,'neighbours')) d <- dneigh(x, 0, d[1])
            d <- d@neighbours
            
            if (missing(drop) || !is.logical(drop[1])) drop <- FALSE
            else drop <- drop[1]
            
            if (missing(zcol)) {
              if (ncol(x) > 1) stop("zcol should be specified!")
              else zcol <- 1
            } else if (is.character(zcol)) {
              w <- which(names(x) == zcol[1])
              if (w == 0) stop('the specified variable in zcol does not exist in the data')
              zcol <- w
            } else if (is.numeric(zcol)) {
              zcol <- zcol[1]
              if (zcol > ncol(x)) stop('the zcol number is greater than the number of columns in data!')
            } else stop("zcol should be a character or a number!")
            
            xx <- x
            
            x <- data.frame(x)[,zcol]
            
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
                  if (verbose) cat("input data is considered categorical, and nc is ignored!\n")
                }
              } 
            } else {
              if (missing(categorical) && !missing(dif) && !is.null(classes)) categorical <- TRUE
            }
            #----
            if (missing(categorical) || !is.logical(categorical)) {
              # guessing whether the layer is categorical:
              if (.is.categorical(x)) {
                categorical <- TRUE
                if (verbose) cat("the specified variable is considered as categorical...\n")
              } else {
                categorical <- FALSE
                if (verbose) cat("the specified variable is considered continuous...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
            } else if (categorical) {
              if (is.null(classes) || is.na(classes)) {
                if (missing(dif) || is.null(classes) || is.na(classes) ) {
                  classes <- unique(x)
                } else {
                  if (length(names(dif)) > 1) {
                    classes <- names(dif)
                    .ux <- as.character(unique(x))
                    if (!all(.ux %in% classes)) classes <- .ux
                  } else classes <- unique(x)
                }
              } else {
                .ux <- unique(x)
                if (is.character(classes)) .ux <- as.character(.ux)
                if (!all(.ux %in% classes)) stop('the specified "classes" does not cover all or some of values in the input raster!')
              }
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
            
            
            if (categorical) {
              x <- .Call('elsac_vector', x, d, as.integer(nc), as.integer(classes),dif, PACKAGE='elsa')
            } else {
              x <- .Call('elsa_vector', x, d, as.integer(nc), PACKAGE='elsa')
            }
            
            xx$Ea <- x[[2]]
            xx$Ec <- x[[1]]
            xx$ELSA <- x[[1]] * x[[2]]
            
            
            if (!drop) xx
            else data.frame(xx)
            
          }
)  


setMethod('elsa', signature(x='SpatialPolygonsDataFrame'), 
          function(x,d,nc,categorical,dif,classes,stat,zcol,drop,method,verbose=TRUE,...) {
            if (missing(classes)) classes <- NULL
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(d)) stop('d is missed!')
            else if (!class(d) %in% c('numeric','integer','neighbours')) stop('d should be either a number (distance) or an object of class neighbours (created by dneigh function')
            
            if (missing(method)) method <- 'centroid'
            
            if (!inherits(d,'neighbours')) d <- dneigh(x, 0, d[1],method = method)
            d <- d@neighbours
            if (missing(drop) || !is.logical(drop[1])) drop <- FALSE
            else drop <- drop[1]
            
            
            
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
                  if (verbose) cat("input data is considered categorical, and nc is ignored!\n")
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
                if (verbose) cat("the specified variable is considered as categorical...\n")
              } else {
                categorical <- FALSE
                if (verbose) cat("the specified variable is considered continuous...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
            } else if (categorical) {
              if (is.null(classes) || is.na(classes)) {
                if (missing(dif) || is.null(classes) || is.na(classes) ) {
                  classes <- unique(x)
                } else {
                  if (length(names(dif)) > 1) {
                    classes <- names(dif)
                    .ux <- as.character(unique(x))
                    if (!all(.ux %in% classes)) classes <- .ux
                  } else classes <- unique(x)
                }
              } else {
                .ux <- unique(x)
                if (is.character(classes)) .ux <- as.character(.ux)
                if (!all(.ux %in% classes)) stop('the specified "classes" does not cover all or some of values in the input raster!')
              }
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
            
            if (categorical) {
              x <- .Call('elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif, PACKAGE='elsa')
            } else {
              x <-.Call('elsa_vector', as.integer(x), d, as.integer(nc), PACKAGE='elsa')
            }
            
            xx@data$Ea <- x[[2]]
            xx@data$Ec <- x[[1]]
            xx@data$ELSA <- x[[1]] * x[[2]]
            xx@data <- xx@data[,c('Ea','Ec','ELSA')]
            
            if (!drop) xx
            else xx@data
            
          }
)  


