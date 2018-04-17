# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2016
# Version 1.1
# Licence GPL v3 

setClassUnion('data.frameORmatrix',c("data.frame","matrix"))
setClassUnion('matrixORnull',c("NULL","matrix"))

setClass("Entrogram",
         representation(width="numeric",
                        cutoff="numeric",
                        entrogramCloud='matrixORnull',
                        entrogram="data.frame"
                        ),
         prototype(
           entrogramCloud=NULL
           )
)
#-------

setClass("Variogram",
         representation(width="numeric",
                        cutoff="numeric",
                        variogramCloud="matrix",
                        variogram="data.frame")
)
#-------

setClass("Correlogram",
         representation(width="numeric",
                        cutoff="numeric",
                        correlogram="data.frame")
)
#-------
setClass("neighbours",
         representation(distance1="numeric",
                        distance2="numeric",
                        neighbours='list'
         )
)

