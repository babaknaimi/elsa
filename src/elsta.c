/* Babak Naimi, May 2018
 * naimi.b@gmail.com
 * last update: November 2022
 * v 1.1
 */

// Entropy-based Local indicator of Spatio-Temporal Associations (ELSTA):
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Rdynload.h"





SEXP elsta(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n, nv; 
  double e, w, s, xi, qq, count, a, xj;
  
  R_len_t i, j;
  
  SEXP ans;
  
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  int *xrr, *xcc;
  //double *xv;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0];
  
  n = length(VECTOR_ELT(v,0));
  
  nv=length(v);
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  //PROTECT(v = coerceVector(v, REALSXP));
  //++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  //xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  
  for (c=0;c < n;c++)  {
    R_CheckUserInterrupt();
    xi=NUMERIC_POINTER(VECTOR_ELT(v,0))[c];
    
    if (R_finite(xi)) {
      row = (c / ncol) + 1;
      col = (c + 1) - ((row - 1) * ncol);
      
      double xn[ngb * nv];
      q=-1;
      for (i=0; i < ngb; i++) {
        nnr= row + xrr[i];
        nnc = col + xcc[i];
        
        if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
          cellnr = ((nnr - 1) * ncol) + nnc;
          for (j=0;j < nv; j++) {
            //xi=NUMERIC_POINTER(VECTOR_ELT(v,j))[c];
            if (R_finite(xi)) {
              xj = NUMERIC_POINTER(VECTOR_ELT(v,j))[(cellnr-1)];
              if (R_finite(xj)) {
                q+=1;
                xn[q]=xj;
              }
            }
          }
        }
      }
      
      // sort
      for (i=0;i <= (q-1);i++) {
        for (j=i+1;j <= q;j++) {
          if (xn[i] > xn[j]) {
            a=xn[i];
            xn[i]=xn[j];
            xn[j]=a;
          }
        }
      }
      //------
      
      a=xn[0];
      count=1;
      e=0;
      qq=q+1;
      
      for (i=1;i <= q;i++) {
        if (xn[i] != a) {
          e = e + ((count / qq) * log2(count / qq));
          a=xn[i];
          count=1;
        } else {
          count+=1;
        }
      }
      e = e + ((count / qq) * log2(count / qq));
      w=0;
      for (i=0; i <= q;i++) {
        w = w + fabs(xn[i] - xi);
      }
      w = w / ((qq - 1) * (ncl - 1));
      
      if (qq > ncl) {
        s = log2(ncl);
      } else {
        s = log2(qq);
      }
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = -e / s;
      //xans[c] = (-e * w) / s;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = w;
    } else {
      //xans[c]=R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
}

//-----------------
