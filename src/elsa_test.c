/* Babak Naimi, August 2014 
   naimi.b@gmail.com
   August 2016
   v 2.0
*/


#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Random.h"
#include "R_ext/Rdynload.h"
#include <math.h>

//////

static void SampleReplace(int k, int n, int *y)
{
  int i;
  for (i = 0; i < k; i++)
    y[i] = n * unif_rand() + 1;
}

/* Equal probability sampling; without-replacement case */

static void SampleNoReplace(int k, int n, int *y, int *x)
{
  int i, j;
  for (i = 0; i < n; i++)
    x[i] = i;
  for (i = 0; i < k; i++) {
    j = n * unif_rand();
    y[i] = x[j] + 1;
    x[j] = x[--n];
  }
}
//----------------
static void perm(double *xv, int *xNA, int noNA, int type) {
  int i;
  
  // generate random permutation
  SEXP x,y;
  int *xx,*xy;
  PROTECT(y = allocVector(INTSXP, noNA));
  PROTECT(x = allocVector(INTSXP, noNA));
  xx=INTEGER(x);
  xy=INTEGER(y);
  
  GetRNGstate();
  if (type == 1) {
    SampleNoReplace(noNA, noNA, xy, xx);
  } else {
    SampleReplace(noNA, noNA, xy);
  }
  PutRNGstate();
  
  for (i=0;i < noNA;i++) xx[i]=xv[xNA[(xy[i]-1)]];
  for (i=0;i < noNA;i++) xv[xNA[i]]=xx[i];
  UNPROTECT(2);
  
}
//----------------

static void elsaCalc(double *xv, double *xans, int ncol, int nrow, int ncl, int *xrr, int *xcc,int ngb, int n) {
  int c, row, col, q, nnr, nnc, cellnr, i, j;
  double e, w, s, xi, qq, count, a;
  
  for (c=0;c < n;c++)  {
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
}
//-------------
static void elsaCalc_cell(double *xv,double *xans, int ncol, int nrow, int ncl, int *xrr, int *xcc,int ngb, int n, int *xcells) {
  int c, row, col, q, nnr, nnc, cellnr, cn, i, j;
  double e, w, s, xi, qq, count, a;
  
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
}


//---------------
static void elsa_vectorCalc(double *xv, double *xans, int ncl, int n, SEXP nb) {
  int  q, ngb, i, j, c;
  double e, w, s,  qq, count, xi, a;
  
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
}


///////

static void elsacCalc(double *xv, double *xans, int ncol, int nrow, int ncl, int *xrr, int *xcc,int ngb, int n, int *xcls, double *xdif, double maxW) {
  int c, nw, q, nnr, nnc, cellnr, row, col; 
  double e, w, s, xi, qq, count, a;
  
  R_len_t i, j;
  for (c=0;c < n;c++)  {
    xi=xv[c];
    if (!R_IsNA(xi)) {
      row = (c / ncol) + 1;
      col = (c + 1) - ((row - 1) * ncol);
      
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
      //------
      
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
}


static void elsac_cellCalc(double *xv, double *xans, int ncol, int nrow, int ncl, int *xrr, int *xcc,int ngb, int n, int *xcls, double *xdif, double maxW, int *xcells) {
  int c, row, col,  q, nnr, nnc, cellnr, nw, cn, i, j;
  double e, w, s, xi, qq, count, a;
  
  for (c=0;c < n;c++)  {
    cn=xcells[c]-1;
    xi=xv[cn];
    if (!R_IsNA(xi)) {
      row = (cn / ncol) + 1;
      col = (cn + 1) - ((row - 1) * ncol);
      
      double xn[ngb], xw[ngb];
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
}


////
static void elsac_vectorCalc(double *xv, double *xans, int ncl,int n, int *xcls, double *xdif, double maxW, SEXP nb) {
  int c, nw, q, i, j, ngb; 
  double e, w, s, xi, qq, count, a;
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
  
}
//////////////////##############################
/////////////////###############################

// type=1 : shuffle (sampling without replacement)
// type=2 : bootstrap (sampling with replacement)

SEXP elsa_test(SEXP v, SEXP null, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP type, SEXP nperm) {
  int nProtected=0;
  int nrow, ncol, ncl, n, ngb,tp, np,k;
  
  SEXP ans,tmp,ansnull;
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
  
  elsaCalc(xv, xans,ncol,nrow,ncl,xrr, xcc,ngb,n);
  
  //--------------
  double *xtmp, *xnull, *xansnull;
  
  tp=INTEGER(type)[0];
  np=INTEGER(nperm)[0];
  
  PROTECT(null = coerceVector(null, REALSXP));
  ++nProtected;
  PROTECT(tmp = allocVector(REALSXP, n));
  ++nProtected;
  PROTECT(ansnull = allocVector(REALSXP, n));
  ++nProtected;
  
  xnull=REAL(null);
  xtmp=REAL(tmp);
  xansnull=REAL(ansnull);
  //--- Random Permutation:
  int i,j;
  // make a vector of which is not NA
  int noNA=0;
  int *xNA=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xNA[noNA]=i;
      noNA++;
    }
  }
  if (noNA < n) xNA=realloc(xNA,noNA*sizeof(int));
  
  float *xPV=malloc(noNA*sizeof(float));
  for (i=0;i < noNA;i++) xPV[i]=0;
  //------
  for (k=0;k < np;k++) {
    for (i=0;i<n;i++) xtmp[i]=xnull[i];
    
    perm(xtmp,xNA,noNA,tp);
    
    elsaCalc(xtmp, xansnull,ncol,nrow,ncl,xrr, xcc,ngb,n);
    
    for (i=0;i < noNA;i++) {
      j=xNA[i];
      if(xansnull[j] <= xans[j]) xPV[i]=(xPV[i]+1);
    }
  }
  
  for (i=0;i < noNA;i++) xans[xNA[i]]=(xPV[i]+1)/(np+1);
  
  free(xNA);
  free(xPV);
  
  //--------
  UNPROTECT(nProtected);
  return(ans);
  
}
////########### ELSA CELL
SEXP elsa_cell_test(SEXP v, SEXP null, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP cells, SEXP type, SEXP nperm) {
  int nProtected=0;
  int ngb, nrow, ncol, ncl, n, tp, np, k,nv, i, j;
  
  SEXP ans, tmp, ansnull;
  double *xans, *xv;
  int *xrr, *xcc, *xcells;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0];
  
  nv=length(v);
  
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
  
  elsaCalc_cell(xv,xans,ncol,nrow,ncl, xrr, xcc,ngb,n,xcells);
  
  //--------------
  
  double *xtmp, *xnull,*xansnull;
  
  tp=INTEGER(type)[0];
  np=INTEGER(nperm)[0];
  
  PROTECT(null = coerceVector(null, REALSXP));
  ++nProtected;
  PROTECT(tmp = allocVector(REALSXP, nv));
  ++nProtected;
  PROTECT(ansnull = allocVector(REALSXP, n));
  ++nProtected;
  
  xnull=REAL(null);
  xtmp=REAL(tmp);
  xansnull=REAL(ansnull);
  //--- Random Permutation:
  
  // make a vector of which is not NA
  int noNA=0;
  int *xNA=malloc(nv*sizeof(int));
  
  for (i=0;i < nv;i++) {
    if (!R_IsNA(xv[i])) {
      xNA[noNA]=i;
      noNA++;
    }
  }
  if (noNA < nv) xNA=realloc(xNA,noNA*sizeof(int));
  //// noNA and xNA for cells i.e. noNAc and xNAc
  int noNAc=0;
  int *xNAc=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    j=xcells[i]-1;
    if (!R_IsNA(xv[j])) {
      xNAc[noNAc]=i;
      noNAc++;
    }
  }
  if (noNAc < n) xNAc=realloc(xNAc,noNAc*sizeof(int));
  
  float *xPV=malloc(noNAc*sizeof(float));
  for (i=0;i < noNAc;i++) xPV[i]=0;
  //------
  for (k=0;k < np;k++) {
    for (i=0;i < nv;i++) xtmp[i]=xnull[i];
    
    perm(xtmp,xNA,noNA,tp);
    
    elsaCalc_cell(xtmp, xansnull,ncol,nrow,ncl, xrr, xcc,ngb,n,xcells);
    
    
    for (i=0;i < noNAc;i++) {
      j=xNAc[i];
      if(xansnull[j] <= xans[j]) xPV[i]=(xPV[i]+1);
    }
  }
  
  for (i=0;i < noNAc;i++) xans[xNAc[i]]=(xPV[i]+1)/(np+1);
  
  free(xNA);
  free(xPV);
  free(xNAc);
  
  UNPROTECT(nProtected);
  return(ans);
  
}
//----------

///////////////////////

SEXP elsac_vector_test(SEXP v, SEXP null, SEXP nb,  SEXP nclass, SEXP classes, SEXP dif, SEXP type, SEXP nperm) {
  int nProtected=0;
  int ncl, n, k, np, tp;
  double maxW;
  
  R_len_t i, j;
  
  SEXP ans,tmp,ansnull;
  
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
  
  elsac_vectorCalc(xv, xans, ncl,n,xcls,xdif,maxW,nb);
  
  //--------------
  
  double *xtmp, *xnull,*xansnull;
  
  tp=INTEGER(type)[0];
  np=INTEGER(nperm)[0];
  
  
  PROTECT(null = coerceVector(null, REALSXP));
  ++nProtected;
  PROTECT(tmp = allocVector(REALSXP, n));
  ++nProtected;
  PROTECT(ansnull = allocVector(REALSXP, n));
  ++nProtected;
  
  xnull=REAL(null);
  xtmp=REAL(tmp);
  xansnull=REAL(ansnull);
  //--- Random Permutation:
  
  // make a vector of which is not NA
  int noNA=0;
  int *xNA=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xNA[noNA]=i;
      noNA++;
    }
  }
  if (noNA < n) xNA=realloc(xNA,noNA*sizeof(int));
  
  float *xPV=malloc(noNA*sizeof(float));
  for (i=0;i < noNA;i++) xPV[i]=0;
  //------
  for (k=0;k < np;k++) {
    for (i=0;i<n;i++) xtmp[i]=xnull[i];
    
    perm(xtmp,xNA,noNA,tp);
    elsac_vectorCalc(xtmp, xansnull, ncl,n,xcls,xdif,maxW,nb);
    
    for (i=0;i < noNA;i++) {
      j=xNA[i];
      if(xansnull[j] <= xans[j]) xPV[i]=(xPV[i]+1);
    }
  }
  
  for (i=0;i < noNA;i++) xans[xNA[i]]=(xPV[i]+1)/(np+1);
  
  free(xNA);
  free(xPV);
  
  //--------
  UNPROTECT(nProtected);
  return(ans);
  
}
//------


SEXP elsac_test(SEXP v, SEXP null, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP classes, SEXP dif, SEXP type, SEXP nperm) {
  int nProtected=0;
  int ngb,  nrow, ncol,  ncl, n, k, np, tp;
  double maxW;
  
  R_len_t i, j;
  
  SEXP ans,tmp,ansnull;
  
  double *xv, *xans, *xdif;
  int *xrr, *xcc, *xcls;
  
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
  
  
  elsacCalc(xv, xans, ncol, nrow, ncl, xrr, xcc,ngb,n,xcls,xdif,maxW);  
  
  //--------------
  double *xtmp, *xnull,*xansnull;
  
  tp=INTEGER(type)[0];
  np=INTEGER(nperm)[0];
  
  
  PROTECT(null = coerceVector(null, REALSXP));
  ++nProtected;
  PROTECT(tmp = allocVector(REALSXP, n));
  ++nProtected;
  PROTECT(ansnull = allocVector(REALSXP, n));
  ++nProtected;
  
  xnull=REAL(null);
  xtmp=REAL(tmp);
  xansnull=REAL(ansnull);
  //--- Random Permutation:
  // make a vector of which is not NA
  int noNA=0;
  int *xNA=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xNA[noNA]=i;
      noNA++;
    }
  }
  if (noNA < n) xNA=realloc(xNA,noNA*sizeof(int));
  
  float *xPV=malloc(noNA*sizeof(float));
  for (i=0;i < noNA;i++) xPV[i]=0;
  //------
  for (k=0;k < np;k++) {
    for (i=0;i<n;i++) xtmp[i]=xnull[i];
    
    perm(xtmp,xNA,noNA,tp);
    
    elsacCalc(xtmp, xansnull, ncol, nrow, ncl, xrr, xcc,ngb,n,xcls,xdif,maxW);  
    
    for (i=0;i < noNA;i++) {
      j=xNA[i];
      if(xansnull[j] <= xans[j]) xPV[i]=(xPV[i]+1);
    }
  }
  
  for (i=0;i < noNA;i++) xans[xNA[i]]=(xPV[i]+1)/(np+1);
  
  free(xNA);
  free(xPV);
  
  //--------
  UNPROTECT(nProtected);
  return(ans);
}
//----------
SEXP elsac_cell_test(SEXP v, SEXP null, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP classes, SEXP dif, SEXP cells, SEXP type, SEXP nperm) {
  int nProtected=0;
  int ngb, nrow, ncol, ncl, n, np, tp, k, nv, i,j;
  double maxW;
  
  SEXP ans, tmp, ansnull;
  
  double *xans, *xdif, *xv;
  int *xrr, *xcc, *xcls, *xcells;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(cells);
  nv=length(v);
  
  
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
  
  elsac_cellCalc(xv,xans,ncol,nrow,ncl,xrr, xcc, ngb, n, xcls,xdif,maxW,xcells);
  
  //--------------
  double *xtmp, *xnull,*xansnull;
  
  tp=INTEGER(type)[0];
  np=INTEGER(nperm)[0];
  
  
  PROTECT(null = coerceVector(null, REALSXP));
  ++nProtected;
  PROTECT(tmp = allocVector(REALSXP, nv));
  ++nProtected;
  PROTECT(ansnull = allocVector(REALSXP, n));
  ++nProtected;
  
  xnull=REAL(null);
  xtmp=REAL(tmp);
  xansnull=REAL(ansnull);
  
  /////
  
  // make a vector of which is not NA
  int noNA=0;
  int *xNA=malloc(nv*sizeof(int));
  
  for (i=0;i < nv;i++) {
    if (!R_IsNA(xv[i])) {
      xNA[noNA]=i;
      noNA++;
    }
  }
  if (noNA < nv) xNA=realloc(xNA,noNA*sizeof(int));
  //// noNA and xNA for cells i.e. noNAc and xNAc
  int noNAc=0;
  int *xNAc=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    j=xcells[i]-1;
    if (!R_IsNA(xv[j])) {
      xNAc[noNAc]=i;
      noNAc++;
    }
  }
  if (noNAc < n) xNAc=realloc(xNAc,noNAc*sizeof(int));
  
  float *xPV=malloc(noNAc*sizeof(float));
  for (i=0;i < noNAc;i++) xPV[i]=0;
  //------
  for (k=0;k < np;k++) {
    for (i=0;i < nv;i++) xtmp[i]=xnull[i];
    
    perm(xtmp,xNA,noNA,tp);
    elsac_cellCalc(xtmp, xansnull,ncol,nrow,ncl,xrr, xcc, ngb, n, xcls,xdif,maxW,xcells);
    
    for (i=0;i < noNAc;i++) {
      j=xNAc[i];
      if(xansnull[j] <= xans[j]) xPV[i]=(xPV[i]+1);
    }
  }
  
  for (i=0;i < noNAc;i++) xans[xNAc[i]]=(xPV[i]+1)/(np+1);
  
  free(xNA);
  free(xPV);
  free(xNAc);
  //--------
  UNPROTECT(nProtected);
  return(ans);
}

////------

SEXP elsa_vector_test(SEXP v, SEXP null, SEXP nb, SEXP nclass, SEXP type, SEXP nperm) {
  int nProtected=0;
  int  ncl, n, k, tp, np;
  
  SEXP ans, tmp, ansnull;
  double *xans, *xv;
  
  ncl=INTEGER(nclass)[0];
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  xans=REAL(ans);
  xv=REAL(v);
  
  elsa_vectorCalc(xv,xans,ncl,n,nb);
  
  //--------------
  double *xtmp, *xnull,*xansnull;
  
  tp=INTEGER(type)[0];
  np=INTEGER(nperm)[0];
  
  
  PROTECT(null = coerceVector(null, REALSXP));
  ++nProtected;
  PROTECT(tmp = allocVector(REALSXP, n));
  ++nProtected;
  PROTECT(ansnull = allocVector(REALSXP, n));
  ++nProtected;
  
  xnull=REAL(null);
  xtmp=REAL(tmp);
  xansnull=REAL(ansnull);
  //--- Random Permutation:
  int i,j;
  // make a vector of which is not NA
  int noNA=0;
  int *xNA=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xNA[noNA]=i;
      noNA++;
    }
  }
  if (noNA < n) xNA=realloc(xNA,noNA*sizeof(int));
  
  float *xPV=malloc(noNA*sizeof(float));
  for (i=0;i < noNA;i++) xPV[i]=0;
  //------
  for (k=0;k < np;k++) {
    for (i=0;i<n;i++) xtmp[i]=xnull[i];
    
    perm(xtmp,xNA,noNA,tp);
    
    elsa_vectorCalc(xtmp, xansnull,ncl,n,nb);
    
    for (i=0;i < noNA;i++) {
      j=xNA[i];
      if(xansnull[j] <= xans[j]) xPV[i]=(xPV[i]+1);
    }
  }
  
  for (i=0;i < noNA;i++) xans[xNA[i]]=(xPV[i]+1)/(np+1);
  
  free(xNA);
  free(xPV);
  
  //--------
  
  UNPROTECT(nProtected);
  return(ans);
}
//---------
