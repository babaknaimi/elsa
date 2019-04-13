# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2019
# last update: April 2019
# Version 1.3
# Licence GPL v3 
#------------------

# low-level functions to be used within the other functions

.elsaContinousVector <- function(x,d,nc,ncl,nrw,res) {
  w <-.Filter(r=res[1],d1=0,d2=d)
  fdim <- w[[1]]
  w <- w[[2]]
  
  .Call('v_elsa', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa') 
}
#------
.elsaContinousVectorCell <- function(x,d,nc,ncl,nrw,res,cells) {
  w <-.Filter(r=res[1],d1=0,d2=d)
  fdim <- w[[1]]
  w <- w[[2]]
  .Call('v_elsa_cell', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
}

###########################

.elsaEcEaContinousVector <- function(x,d,nc,ncl,nrw,res) {
  w <-.Filter(r=res[1],d1=0,d2=d)
  fdim <- w[[1]]
  w <- w[[2]]
  
  xx <- .Call('elsa', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), PACKAGE='elsa')
  names(xx) <- c('Ec','Ea')
  return(xx)
}
#---------
.elsaEcEaContinousVectorCells <- function(x,d,nc,ncl,nrw,res,cells) {
  w <-.Filter(r=res[1],d1=0,d2=d)
  fdim <- w[[1]]
  w <- w[[2]]
  
  xx <- .Call('elsa_cell', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells), PACKAGE='elsa')
  names(xx) <- c('Ec','Ea')
  return(xx)
}
#########################

.elsaCategoricalVector <- function(x,d,dif,classes,ncl,nrw,res) {
  #classes <- unique(x)
  #classes <- classes[!is.na(classes)]
  
  nc <- length(classes)
  
  if (missing(dif)) {
    dif <- rep(1,nc*nc)
    for (i in 1:nc) dif[(i-1)*nc+i] <-0
  } else {
    dif <- .checkDif(dif,classes)
  }
  #-----
  w <-.Filter(r=res[1],d1=0,d2=d)
  
  fdim <- w[[1]]
  w <- w[[2]]
  
  .Call('v_elsac', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
}
#-----

.elsaCategoricalVectorCell <- function(x,d,dif,classes,ncl,nrw,res,cells) {
  #if (missing(d)) d <- res[1]
  #----
  # classes <- unique(x)
  # classes <- classes[!is.na(classes)]
   
  nc <- length(classes)
  
  if (missing(dif)) {
    dif <- rep(1,nc*nc)
    for (i in 1:nc) dif[(i-1)*nc+i] <-0
  } else {
    dif <- .checkDif(dif,classes)
  }
  #-----
  w <-.Filter(r=res[1],d1=0,d2=d)
  
  fdim <- w[[1]]
  w <- w[[2]]
  
  .Call('v_elsac_cell', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
}
#---------
.elsaEcEaCategoricalVector <- function(x,d,dif,classes,ncl,nrw,res) {
  # classes <- unique(x)
  # classes <- classes[!is.na(classes)]
  
  nc <- length(classes)
  
  if (missing(dif)) {
    dif <- rep(1,nc*nc)
    for (i in 1:nc) dif[(i-1)*nc+i] <-0
  } else {
    dif <- .checkDif(dif,classes)
  }
  #-----
  w <-.Filter(r=res[1],d1=0,d2=d)
  
  fdim <- w[[1]]
  w <- w[[2]]
  xx <- .Call('elsac', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, PACKAGE='elsa')
  names(xx) <- c('Ec','Ea')
  return(xx)
}
#-----

.elsaEcEaCategoricalVectorCell <- function(x,d,dif,classes,ncl,nrw,res,cells) {
  #if (missing(d)) d <- res[1]
  #----
  # classes <- unique(x)
  # classes <- classes[!is.na(classes)]
  
  nc <- length(classes)
  
  if (missing(dif)) {
    dif <- rep(1,nc*nc)
    for (i in 1:nc) dif[(i-1)*nc+i] <-0
  } else {
    dif <- .checkDif(dif,classes)
  }
  #-----
  w <-.Filter(r=res[1],d1=0,d2=d)
  
  fdim <- w[[1]]
  w <- w[[2]]
  
  xx <- .Call('elsac_cell', x[], as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells), PACKAGE='elsa')
  names(xx) <- c('Ec','Ea')
  return(xx)
}