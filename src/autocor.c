/* Babak Naimi, July 2016
   naimi.b@gmail.com
   July 2016
   v 1.0
*/


#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Random.h"
#include "R_ext/Rdynload.h"


SEXP moran(SEXP v, SEXP nc, SEXP nr, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, nnr, nnc, nrow, ncol, cellnr, n, q;
  double xi;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv, *xans;
  int *xrr, *xcc ;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  
  PROTECT(ans = allocVector(REALSXP, 1));
  ++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  
  int nv=0;
  double sigmaX=0, yi2=0, ybar, yi, yj=0, wij=0;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      nv++;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  ybar=sigmaX/nv;
  for (i=0;i < nv;i++) {
    j=xcells[i];
    yi2=yi2 + pow((xv[j] - ybar),2);
  }
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    row = (j / ncol) + 1;
    col = (j + 1) - ((row - 1) * ncol);
    
    double xn[ngb];
    q=-1;
    for (i=0; i < ngb; i++) {
      nnr= row + xrr[i];
      nnc = col + xcc[i];
      
      if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
        cellnr = ((nnr - 1) * ncol) + nnc;
        if (!R_IsNA(xv[(cellnr-1)])) {
          q+=1;
          xn[q]=xv[(cellnr-1)];
        }
      }
    }
    //----
    yi=(xi - ybar);
    ////----
    
    for (i=0;i <= q;i++) {
      yj=yj+((xn[i] - ybar) * yi);
    }
    
    yj=yj - (yi * yi);
    /////----
    wij=wij+q;
  }
  xans[0]= (nv * yj) / (yi2 * wij);
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}



SEXP geary(SEXP v, SEXP nc, SEXP nr, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, nnr, nnc, nrow, ncol, cellnr, n, q;
  double xi;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv, *xans;
  int *xrr, *xcc ;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  
  PROTECT(ans = allocVector(REALSXP, 1));
  ++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  
  int nv=0;
  double sigmaX=0, yi2=0, ybar, yj=0, wij=0;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      nv++;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  ybar=sigmaX/nv;
  for (i=0;i < nv;i++) {
    j=xcells[i];
    yi2=yi2 + pow((xv[j] - ybar),2);
  }
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    row = (j / ncol) + 1;
    col = (j + 1) - ((row - 1) * ncol);
    
    double xn[ngb];
    q=-1;
    for (i=0; i < ngb; i++) {
      nnr= row + xrr[i];
      nnc = col + xcc[i];
      
      if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
        cellnr = ((nnr - 1) * ncol) + nnc;
        if (!R_IsNA(xv[(cellnr-1)])) {
          q+=1;
          xn[q]=xv[(cellnr-1)];
        }
      }
    }
    //----
    
    for (i=0;i <= q;i++) {
      yj=yj+pow((xi - xn[i]),2);
    }
    wij=wij+q;
  }
  xans[0]= ((nv - 1) * yj) / (2* yi2 * wij);
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
//////

SEXP geary_vector(SEXP v, SEXP nb) {
  int nProtected=0;
  int c, ngb, n, q;
  double xi, a;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv, *xans;
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  
  PROTECT(ans = allocVector(REALSXP, 1));
  ++nProtected;
  
  
  xans=REAL(ans);
  xv=REAL(v);
  
  int nv=0;
  double sigmaX=0, yi2=0, ybar, yj=0, wij=0;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      nv++;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  ybar=sigmaX/nv;
  for (i=0;i < nv;i++) {
    j=xcells[i];
    yi2=yi2 + pow((xv[j] - ybar),2);
  }
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    ngb = length(VECTOR_ELT(nb,c));
    
    double xn[ngb];
    q=-1;
    for (i=0;i < ngb;i++) {
      a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
      if (!R_IsNA(a)) {
        q+=1;
        xn[q]=a;
      }
    }
    //----
    for (i=0;i <= q;i++) {
      yj=yj+pow((xi - xn[i]),2);
    }
    wij=wij+q+1;
  }
  xans[0]= ((nv - 1) * yj) / (2* yi2 * wij);
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
//////

SEXP moran_vector(SEXP v, SEXP nb) {
  int nProtected=0;
  int c, ngb, n, q;
  double xi, a;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv, *xans;
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  
  PROTECT(ans = allocVector(REALSXP, 1));
  ++nProtected;
  
  xans=REAL(ans);
  xv=REAL(v);
  
  int nv=0;
  double sigmaX=0, yi2=0, ybar, yi, yj=0, wij=0;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      nv++;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  ybar=sigmaX/nv;
  for (i=0;i < nv;i++) {
    j=xcells[i];
    yi2=yi2 + pow((xv[j] - ybar),2);
  }
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    
    ngb = length(VECTOR_ELT(nb,c));
    
    double xn[ngb];
    q=-1;
    for (i=0;i < ngb;i++) {
      a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
      if (!R_IsNA(a)) {
        q+=1;
        xn[q]=a;
      }
    }
    //----
    yi=(xi - ybar);
    ////----
    
    for (i=0;i <= q;i++) {
      yj=yj+((xn[i] - ybar) * yi);
    }
    
    /////----
    wij=wij+q+1;
  }
  xans[0]= (nv * yj) / (yi2 * wij);
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
/////
