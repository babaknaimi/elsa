/* Babak Naimi, July 2016
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



SEXP elsac(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP classes, SEXP dif) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n,  nw;
  double e, w, s, xi, qq, count, maxW, a;
  
  R_len_t i, j;
  
  SEXP ans;
  
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  double *xv, *xdif;
  int  *xrr, *xcc, *xcls;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(v);
  
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  PROTECT(classes = coerceVector(classes, INTSXP));
  ++nProtected;
  
  PROTECT(dif = coerceVector(dif, REALSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  xcls=INTEGER(classes);
  xdif=REAL(dif);
  
  maxW=0;
  for (i=0; i < length(dif);i++) {
    if (xdif[i] > maxW) maxW=xdif[i];
  }
  
  for (c=0;c < n;c++)  {
    xi=xv[c];
    if (!R_IsNA(xi)) {
      row = (c / ncol) + 1;
      col = (c + 1) - ((row - 1) * ncol);
      
      double xw[ngb],xn[ngb];
      
      //------
      for (i=0; i < ncl;i++) {
        if (xcls[i] == (int)xi) {
          nw=i;
          break;
        }
      }
      //-------
      
      q=-1;
      for (i=0; i < ngb; i++) {
        nnr= row + xrr[i];
        nnc = col + xcc[i];
        
        if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
          cellnr = ((nnr - 1) * ncol) + nnc;
          if (!R_IsNA(xv[(cellnr-1)])) {
            q+=1;
            xn[q]=xv[(cellnr-1)];
            for (j=0;j < ncl;j++) {
              if (xcls[j] == xn[q]) {
                xw[q]=xdif[(nw*ncl)+j];
                break;
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
        w = w + xw[i];
      }
      w = w / ((qq - 1) * maxW);
      
      if (qq > ncl) {
        s = log2(ncl);
      } else {
        s = log2(qq);
      }
      
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = -e / s;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = w;
    } else {
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
  
}
//------

SEXP elsac_cell(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP classes, SEXP dif, SEXP cells) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n, nw, cn;
  double e, w, s, xi, qq, count, maxW, a;
  
  R_len_t i, j;
  
  SEXP ans;
  
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  double *xdif, *xv;
  int *xrr, *xcc, *xcls, *xcells;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(cells);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  PROTECT(classes = coerceVector(classes, INTSXP));
  ++nProtected;
  
  PROTECT(dif = coerceVector(dif, REALSXP));
  ++nProtected;
  
  PROTECT(cells = coerceVector(cells, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  xcls=INTEGER(classes);
  xdif=REAL(dif);
  xcells=INTEGER(cells);
  
  maxW=0;
  for (i=0; i < length(dif);i++) {
    if (xdif[i] > maxW) maxW=xdif[i];
  }
  
  for (c=0;c < n;c++)  {
    cn=xcells[c]-1;
    xi=xv[cn];
    if (!R_IsNA(xi)) {
      row = (cn / ncol) + 1;
      col = (cn + 1) - ((row - 1) * ncol);
      
      double xw[ngb], xn[ngb];
      //------
      for (i=0; i < ncl;i++) {
        if (xcls[i] == xi) {
          nw=i;
          break;
        }
      }
      //-------
      
      q=-1;
      for (i=0; i < ngb; i++) {
        nnr= row + xrr[i];
        nnc = col + xcc[i];
        
        if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
          cellnr = ((nnr - 1) * ncol) + nnc;
          if (!R_IsNA(xv[(cellnr-1)])) {
            q+=1;
            xn[q]=xv[(cellnr-1)];
            for (j=0;j < ncl;j++) {
              if (xcls[j] == xn[q]) {
                xw[q]=xdif[(nw*ncl)+j];
                break;
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
      //
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
        w = w + xw[i];
      }
      w = w / ((qq - 1) * maxW);
      
      if (qq > ncl) {
        s = log2(ncl);
      } else {
        s = log2(qq);
      }
      
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = -e / s;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = w;
    } else {
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
  
}

////------

SEXP elsac_vector(SEXP v, SEXP nb,  SEXP nclass, SEXP classes, SEXP dif) {
  int nProtected=0;
  int c,  ngb, q, ncl, n, a, nw;
  double e, w, s, xi, qq, count, maxW;
  
  R_len_t i, j;
  
  SEXP ans;
  
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  double *xv, *xdif;
  int  *xcls;
  
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(v);
  
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  PROTECT(classes = coerceVector(classes, INTSXP));
  ++nProtected;
  
  PROTECT(dif = coerceVector(dif, REALSXP));
  ++nProtected;
  
  xv=REAL(v);
  xcls=INTEGER(classes);
  xdif=REAL(dif);
  
  maxW=0;
  for (i=0; i < length(dif);i++) {
    if (xdif[i] > maxW) maxW=xdif[i];
  }
  
  for (c=0;c < n;c++)  {
    xi=xv[c];
    if (!R_IsNA(xi)) {
      ngb = length(VECTOR_ELT(nb,c));
      double xn[ngb+1], xw[ngb+1];
      //------
      for (i=0; i < ncl;i++) {
        if (xcls[i] == xi) {
          nw=i;
          break;
        }
      }
      //-------
      q=-1;
      for (i=0;i < ngb;i++) {
        a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
        if (!R_IsNA(a)) {
          q+=1;
          xn[q]=a;
          for (j=0;j < ncl;j++) {
            if (xcls[j] == xn[q]) {
              xw[q]=xdif[(nw*ncl)+j];
              break;
            }
          }
        } 
      }
      q+=1;
      xn[q]=xi; //adding also xi to xn array
      for (j=0;j < ncl;j++) {
        if (xcls[j] == xn[q]) {
          xw[q]=xdif[(nw*ncl)+j];
          break;
        }
      }
      ////////
      
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
        w = w + xw[i];
      }
      w = w / ((qq - 1) * maxW);
      
      if (qq > ncl) {
        s = log2(ncl);
      } else {
        s = log2(qq);
      }
      
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = -e / s;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = w;
    } else {
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[c] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[c] = R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
}
//------

//#############################
//#############################
SEXP v_elsac(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP classes, SEXP dif) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n, nw;
  double e, w, s, xi, qq, count, maxW, a;
  
  R_len_t i, j;
  
  SEXP ans;
  
  double *xv, *xans, *xdif;
  int  *xrr, *xcc, *xcls;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(v);
  
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;

  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  PROTECT(classes = coerceVector(classes, INTSXP));
  ++nProtected;
  
  PROTECT(dif = coerceVector(dif, REALSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  xcls=INTEGER(classes);
  xdif=REAL(dif);
  
  maxW=0;
  for (i=0; i < length(dif);i++) {
    if (xdif[i] > maxW) maxW=xdif[i];
  }
  
  for (c=0;c < n;c++)  {
    xi=xv[c];
    if (!R_IsNA(xi)) {
      row = (c / ncol) + 1;
      col = (c + 1) - ((row - 1) * ncol);
      
      double xw[ngb],xn[ngb];
      //------
      for (i=0; i < ncl;i++) {
        if (xcls[i] == xi) {
          nw=i;
          break;
        }
      }
      //-------
      
      q=-1;
      for (i=0; i < ngb; i++) {
        nnr= row + xrr[i];
        nnc = col + xcc[i];
        
        if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
          cellnr = ((nnr - 1) * ncol) + nnc;
          if (!R_IsNA(xv[(cellnr-1)])) {
            q+=1;
            xn[q]=xv[(cellnr-1)];
            for (j=0;j < ncl;j++) {
              if (xcls[j] == xn[q]) {
                xw[q]=xdif[(nw*ncl)+j];
                break;
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
        w = w + xw[i];
      }
      w = w / ((qq - 1) * maxW);
      
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
//------

SEXP v_elsac_cell(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP classes, SEXP dif, SEXP cells) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n, nw, cn, a;
  double e, w, s, xi, qq, count, maxW;
  
  R_len_t i, j;
  
  SEXP ans;
  
  double *xans, *xdif, *xv;
  int *xrr, *xcc, *xcls, *xcells;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(cells);
  
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;

  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  PROTECT(classes = coerceVector(classes, INTSXP));
  ++nProtected;
  
  PROTECT(dif = coerceVector(dif, REALSXP));
  ++nProtected;
  
  PROTECT(cells = coerceVector(cells, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  xcls=INTEGER(classes);
  xdif=REAL(dif);
  xcells=INTEGER(cells);
  
  maxW=0;
  for (i=0; i < length(dif);i++) {
    if (xdif[i] > maxW) maxW=xdif[i];
  }
  
  for (c=0;c < n;c++)  {
    cn=xcells[c]-1;
    xi=xv[cn];
    if (!R_IsNA(xi)) {
      row = (cn / ncol) + 1;
      col = (cn + 1) - ((row - 1) * ncol);
      
      double xw[ngb], xn[ngb];
      //------
      for (i=0; i < ncl;i++) {
        if (xcls[i] == xi) {
          nw=i;
          break;
        }
      }
      //-------
      
      q=-1;
      for (i=0; i < ngb; i++) {
        nnr= row + xrr[i];
        nnc = col + xcc[i];
        
        if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
          cellnr = ((nnr - 1) * ncol) + nnc;
          if (!R_IsNA(xv[(cellnr-1)])) {
            q+=1;
            xn[q]=xv[(cellnr-1)];
            for (j=0;j < ncl;j++) {
              if (xcls[j] == xn[q]) {
                xw[q]=xdif[(nw*ncl)+j];
                break;
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
      //
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
        w = w + xw[i];
      }
      w = w / ((qq - 1) * maxW);
      
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

////------

SEXP v_elsac_vector(SEXP v, SEXP nb,  SEXP nclass, SEXP classes, SEXP dif) {
  int nProtected=0;
  int c,  ngb, q, ncl, n, nw;
  double e, w, s, xi, qq, count, maxW, a;
  
  R_len_t i, j;
  
  SEXP ans;
  
  double *xv, *xans, *xdif;
  int  *xcls;
  
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(v);
  
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  
  PROTECT(classes = coerceVector(classes, INTSXP));
  ++nProtected;
  
  PROTECT(dif = coerceVector(dif, REALSXP));
  ++nProtected;
  
  xans=REAL(ans);
  xv=REAL(v);
  xcls=INTEGER(classes);
  xdif=REAL(dif);
  
  maxW=0;
  for (i=0; i < length(dif);i++) {
    if (xdif[i] > maxW) maxW=xdif[i];
  }
  
  for (c=0;c < n;c++)  {
    xi=xv[c];
    if (!R_IsNA(xi)) {
      ngb = length(VECTOR_ELT(nb,c));
      double xn[ngb+1], xw[ngb+1];
      //------
      for (i=0; i < ncl;i++) {
        if (xcls[i] == xi) {
          nw=i;
          break;
        }
      }
      //-------
      q=-1;
      for (i=0;i < ngb;i++) {
        a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
        if (!R_IsNA(a)) {
          q+=1;
          xn[q]=a;
          for (j=0;j < ncl;j++) {
            if (xcls[j] == xn[q]) {
              xw[q]=xdif[(nw*ncl)+j];
              break;
            }
          }
        } 
      }
      q+=1;
      xn[q]=xi; //adding also xi to xn array
      for (j=0;j < ncl;j++) {
        if (xcls[j] == xn[q]) {
          xw[q]=xdif[(nw*ncl)+j];
          break;
        }
      }
      ////////
      
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
        w = w + xw[i];
      }
      w = w / ((qq - 1) * maxW);
      
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
//------