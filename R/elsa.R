# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2016
# last update: February 2020
# Version 3.1
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
  length(unique(x - round(x))) <= 2
}
#----------
.is.categoricalRaster <- function(x) {
  if (ncell(x) > 1.3e6) {
    x <- sampleRandom(x,1e6)
  }
  length(unique(x - round(x))) == 1
}
#----------
.is.categoricalNumeric <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 1e6) {
    x <- x[sample(length(x),1e6)]
  }
  length(unique(x - round(x))) == 1
}
#----------
# .checkDif <- function(dif,classes) {
#   nc <- length(classes)
#   classes <- as.character(classes)
#   if(is.list(dif)) {
#     if(!is.null(names(dif))) {
#       if (!all(unlist(lapply(classes,function(x) {x %in% names(dif)})))) stop("categories in the dif argument do not match with the categories in the input layer!")
#       dif <- dif[classes]
#     } else {
#       if (length(dif) == nc) names(dif) <- classes
#       else stop("dif content does not match with the categories in the input layer!")
#     }
#     if (!all(unlist(lapply(classes,function(x) {length(dif[[x]]) == nc})))) stop("the length of categories' differences in the dif argument does not match with the number of categories!")
#     dif <- as.vector(as.matrix(as.data.frame(dif)))
#   } else if (is.matrix(dif) || is.data.frame(dif)) {
#     if (ncol(dif) != nrow(dif) || ncol != nc) stop("number of rows and columns in dif should be the same, equal to the number of categories!")
#     if (!all(unlist(lapply(colnames(dif),function(x) {x %in% row.names(dif)})))) {
#       if (!all(unlist(lapply(colnames(dif),function(x) {x %in% classes})))) {
#         dif <- as.vector(as.matrix(as.data.frame(dif)))
#       } else {
#         row.names(dif) <- colnames(dif)
#         dif <- dif[classes,classes]
#         dif <- as.vector(as.matrix(as.data.frame(dif)))
#       }
#     } else {
#       if (!all(unlist(lapply(colnames(dif),function(x) {x %in% classes})))) stop("colnames of dif do not match with the name of categories!")
#       else {
#         dif <- dif[classes,classes]
#         dif <- as.vector(as.matrix(as.data.frame(dif)))
#       }
#     }
#     
#   } else if (is.vector(dif)) {
#     if (length(dif) != nc*nc) stop("dif argument does not have an appropriate structure!")
#   } else stop("dif argument does not have an appropriate structure!")
#   return(dif)
# }
#----------
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
            if (missing(d)) d <- res(x)[1]
            
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
