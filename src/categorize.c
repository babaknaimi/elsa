/* Babak Naimi, August 2014 
   naimi.b@gmail.com
*/
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Rdynload.h"

SEXP categorize(SEXP v, SEXP b) {
  R_len_t i, j;
  int n, nb;
  
  SEXP ans;
  double *xans, *xv, *xb;
  
  n=length(v);
  nb = length(b);
  
  PROTECT(v = coerceVector(v, REALSXP));
  
  PROTECT(ans = allocVector(REALSXP, n));
  
  PROTECT(b = coerceVector(b, REALSXP));
  
  xans=REAL(ans);
  xv=REAL(v);
  xb=REAL(b);
  
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      for (j=0;j < (nb-1);j++) {
        if ((xv[i] > xb[j]) & (xv[i] <= xb[j+1])) {
          xans[i]=j+1;
          break;
        }
      }
    } else {
      xans[i]=R_NaReal;
    }
  }
  UNPROTECT(3);
  return(ans);
}