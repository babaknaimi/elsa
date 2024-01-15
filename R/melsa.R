# Author: Babak Naimi, naimi.b@gmail.com
# Date :  May 2018
# last update: March 2019
# Version 1.1
# Licence GPL v3 

#These statistics are under development!


if (!isGeneric("melsa")) {
  setGeneric("melsa", function(x,d,nc,categorical,dif,classes,stat,...)
    standardGeneric("melsa"))
}



setMethod('melsa', signature(x='SpatRaster'), 
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
            
            out <- rast(x[[1]])
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
                  
                  
                  
                  rr <- lapply(1:nlyr(x),function(i) x[[i]][][,1])
                  xx <- .Call('Melsa', rr, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
                  if (length(stat) > 1) {
                    if ('ea' %in% stat) {
                      outx <- rast(x[[1]])
                      outx[] <- xx[[2]]
                      names(outx) <- 'Ea'
                      out <- c(out,outx)
                    }
                    if ('ec' %in% stat) {
                      outx <- rast(x[[1]])
                      outx[] <- xx[[1]]
                      names(outx) <- 'Ec'
                      out <- c(out,outx)
                    }
                    if ('elsa' %in% stat) {
                      outx <- rast(x[[1]])
                      outx[] <- xx[[2]] * xx[[1]]
                      names(outx) <- 'ELSA'
                      out <- c(out,outx)
                    }
                  } else {
                    if (stat == 'ea') {
                      out[] <- xx[[2]]
                      names(out) <- 'Ea'
                    } else if (stat == 'ec') {
                      out[] <- xx[[1]]
                      names(out) <- 'Ec'
                    } else {
                      out[] <- xx[[1]] * xx[[2]]
                      names(out) <- 'ELSA'
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



