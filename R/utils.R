# Author: Babak Naimi, naimi.b@gmail.com
# Date :  November 2022
# last update: July 2025
# Version 1.3
# Licence GPL v3 
#-----------------------------

.eval <- function(x,env) {
  eval(parse(text=x),envir=env)
}
#----

.canProcessInMemory <- function(x,n=1) {
  # copied partially from mem_info in the terra package!
  opt <- .eval("terra:::spatOptions()",env=environment())
  opt$ncopies = n
  v <- x@pntr$mem_needs(opt)
  return(round(v[5]) != 0)
}

#------
.is.projected <- function(x) {
  if (inherits(x,'Spatial')) {
    if (!is.na(is.projected(x))) {
      is.projected(x)
    } else {
      all(bbox(x)[1,] <= 180) & all(bbox(x)[1,] >= -180) & all(bbox(x)[2,] <= 90) & all(bbox(x)[2,] >= -90)
    }
  } else if (inherits(x,'matrix') || inherits(x,'data.frame')) {
    all(range(x[,1],na.rm=TRUE) <= 180) & all(range(x[,1],na.rm=TRUE) >= -180) & all(range(x[,2],na.rm=TRUE) <= 90) & all(range(x[,2],na.rm=TRUE) >= -90)
  } else if (inherits(x,'SpatRaster') || inherits(x,'SpatVector')) {
    e <- as.vector(ext(x))
    !all(e >= -180 & e <= 180)
  }
  
  
}