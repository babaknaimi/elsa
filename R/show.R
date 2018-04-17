# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2014
# Version 1.0
# Licence GPL v3 


setMethod ('show' , 'Entrogram', 
           function(object) {
             cat('class               :' , class(object), '\n')
             cat('---------------------------------------\n')
             cat('Width               :' , object@width, '\n')
             cat('Cutoff              :' , object@cutoff, '\n')
             cat('Number of lags      : ' , ceiling(object@cutoff / object@width), '\n')
             cat ('\n')
             cat('------ Entrogram data ------','\n')
             if (nrow(object@entrogram) > 10) {
               print(object@entrogram[1:10,])
               cat ('--- and ',nrow(object@entrogram) - 10,' more...!\n')
             } else print(object@entrogram)
           }
)