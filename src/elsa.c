/* Babak Naimi, August 2014 
   naimi.b@gmail.com
   August 2016
   v 3.0
*/
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Rdynload.h"


SEXP elsa(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n;
  double e, w, s, xi, qq, count, a;
  
  R_len_t i, j;
  
  SEXP ans;
  
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  int *xrr, *xcc;
  double *xv;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0];
  
  n=length(v);
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  
  for (c=0;c < n;c++)  {
    R_CheckUserInterrupt();
    xi=xv[c];
    if (!R_IsNA(xi)) {
      row = (c / ncol) + 1;
      col = (c + 1) - ((row - 1) * ncol);
      
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

SEXP elsa_cell(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP cells) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n, cn;
  double e, w, s, xi, qq, count, a;
  
  R_len_t i, j;
  
  SEXP ans;
  
  int *xrr, *xcc, *xcells;
  double *xv;
  
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0];
  
  n=length(cells);
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  PROTECT(cells = coerceVector(cells, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  xcells=INTEGER(cells);
  
  for (c=0;c < n;c++)  {
    R_CheckUserInterrupt();
    cn=xcells[c]-1;
    xi=xv[cn];
    if (!R_IsNA(xi)) {
      row = (cn / ncol) + 1;
      col = (cn + 1) - ((row - 1) * ncol);
      
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
  }
  UNPROTECT(nProtected);
  return(ans);
  
}
//----------


SEXP elsa_vector(SEXP v, SEXP nb, SEXP nclass) {
  int nProtected=0;
  int  ncl, n, q, ngb, a;
  double e, w, s,  qq, count,xi;
  R_len_t i, j, c;
  SEXP ans;
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  double *xv;
  ncl=INTEGER(nclass)[0];
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  xv=REAL(v);
  
  for (c=0;c < n;c++)  {
    R_CheckUserInterrupt();
    xi=xv[c];
    if (!R_IsNA(xi)) {
      
      ngb = length(VECTOR_ELT(nb,c));
      
      double xn[ngb+1];
      q=-1;
      for (i=0;i < ngb;i++) {
        a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
        if (!R_IsNA(a)) {
          q+=1;
          xn[i]=a;
        }
      }
      q+=1;
      xn[q]=xi;
      
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
  }
  UNPROTECT(nProtected);
  return(ans);
}


///###########################################
///###########################################
///###########################################
SEXP v_elsa(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n;
  double e, w, s, xi, qq, count, a;
  
  R_len_t i, j;
  
  SEXP ans;
  double *xans, *xv;
  int *xrr, *xcc;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0];
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
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
  
  for (c=0;c < n;c++)  {
    R_CheckUserInterrupt();
    xi=xv[c];
    if (!R_IsNA(xi)) {
      row = (c / ncol) + 1;
      col = (c + 1) - ((row - 1) * ncol);
      
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
      
      xans[c] = (-e * w) / s;
      
    } else {
      xans[c]=R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
  
}

//-----------------

SEXP v_elsa_cell(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP cells) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n, cn;
  double e, w, s, xi, qq, a, count;
  
  R_len_t i, j;
  
  SEXP ans;
  double *xans, *xv;
  int *xrr, *xcc, *xcells;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0];
  
  n=length(cells);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  PROTECT(cells = coerceVector(cells, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  xcells=INTEGER(cells);
  
  for (c=0;c < n;c++)  {
    R_CheckUserInterrupt();
    cn=xcells[c]-1;
    xi=xv[cn];
    if (!R_IsNA(xi)) {
      row = (cn / ncol) + 1;
      col = (cn + 1) - ((row - 1) * ncol);
      
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
      
      xans[c] = (-e * w) / s;
      
    } else {
      xans[c]=R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
  
}
//----------


SEXP v_elsa_vector(SEXP v, SEXP nb, SEXP nclass) {
  int nProtected=0;
  int  ncl, n, a, q, ngb;
  double e, w, s,  qq, count,xi;
  R_len_t i, j, c;
  
  SEXP ans;
  double *xans, *xv;
  
  ncl=INTEGER(nclass)[0];
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  xans=REAL(ans);
  xv=REAL(v);
  
  for (c=0;c < n;c++)  {
    R_CheckUserInterrupt();
    xi=xv[c];
    if (!R_IsNA(xi)) {
      
      ngb = length(VECTOR_ELT(nb,c));
      
      double xn[ngb+1];
      q=-1;
      for (i=0;i < ngb;i++) {
        a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
        if (!R_IsNA(a)) {
          q+=1;
          xn[i]=a;
        }
      }
      q+=1;
      xn[q]=xi;
      
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
      
      xans[c] = (-e * w) / s;
      
    } else {
      xans[c]=R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
}
//////
