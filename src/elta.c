/* Babak Naimi, May 2018
 * naimi.b@gmail.com
 * last update: November 2022
 * v 1.1
 */
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Rdynload.h"


SEXP elta(SEXP v, SEXP nclass, SEXP lag, SEXP th) {
  // th: threshold specifying the minimum length of the local window within which ELTA can be calculated otherwise NA is returned!
  int nProtected=0;
  int c, ci, ngb, q, ncl, n, nlag, thv;
  double e, w, s, xi, qq, count, a;
  
  R_len_t i, j, k;
  
  SEXP ans;
  
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  double *xv;
  
  nlag=INTEGER(lag)[0];
  thv=INTEGER(th)[0];
  ncl=INTEGER(nclass)[0];
  
  ngb = (nlag * 2) + 1;
  
  thv = thv - 1;
  
  n=length(v);
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  xv=REAL(v);
  
  for (c=nlag;c < (n - nlag);c++)  {
    R_CheckUserInterrupt();
    xi=xv[c];
    if (R_finite(xi)) {
      ci=c-nlag;
      double xn[ngb];
      q=-1;
      for (i=0; i < ngb; i++) {
        if (R_finite(xv[(ci+i)])) {
          q+=1;
          xn[q]=xv[(ci+i)];
        }
      }
      
      if (q > thv) {
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
        NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = R_NaReal;
        NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = R_NaReal;
      }
    } else  {
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = R_NaReal;
    }
  }
  
  //////////////////////////////
  // repeat the procedure for the begining and the end of time series:
  int cr[(nlag * 2)];
  
  q=0;
  for (i=0;i < nlag;i++) {
    cr[i] = i;
    q+=1;
  }
  
  for (i=(n-nlag);i < n;i++) {
    cr[q] = i;
    q+=1;
  }
  /////////////
  
  for (k = 0;k < (nlag * 2); k++)  {
    
    c=cr[k];
    xi=xv[c];
    
    if (R_finite(xi)) {
      if (c < nlag) {
        ngb=(c+nlag+1);
      } else {
        ngb = (n - c + nlag);
      }
      //---
      double xn[ngb];
      //---
      if (c < nlag) {
        q=-1;
        for (i=0; i < ngb; i++) {
          if (R_finite(xv[i])) {
            q+=1;
            xn[q]=xv[i];
          }
        }
      } else {
        q=-1;
        for (i=0; i < ngb; i++) {
          ci = c - nlag + i;
          if (R_finite(xv[ci])) {
            q+=1;
            xn[q]=xv[ci];
          }
        }
      }
      //---
      
      if (q > thv) {
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
