/* Babak Naimi, July 2016
   naimi.b@gmail.com
   July 2016
   v 1.1
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


//////

//-- Local Moran's I
SEXP localmoran(SEXP v, SEXP nc, SEXP nr, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, nnr, nnc, nrow, ncol, cellnr, n, q;
  double xi;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv;
  int *xrr, *xcc ;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = NEW_LIST(2));
  //PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  //xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  //SET_VECTOR_ELT(ans, 2, NEW_NUMERIC(n));
  
  int nv=0;
  double sigmaX=0, sigmaX2=0, meanXV, s2=0, s4=0, b2, xbar, sx, s,Ii, Wi, EI, VarI, zj, qq,wikh, wi2;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      sigmaX2=sigmaX2 + xv[i]*xv[i];
      nv++;
    } else {
      //xans[i]=R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[i] = R_NaReal;
      //NUMERIC_POINTER(VECTOR_ELT(ans, 2))[i] = R_NaReal;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  meanXV=sigmaX/nv;
  for (i=0;i < nv;i++) {
    j=xcells[i];
    s2=s2 + pow((xv[j] - meanXV),2);
    s4=s4 + pow((xv[j] - meanXV),4);
  }
  s2=s2/(double)nv;
  s4=s4/(double)nv;
  b2=s4/pow(s2,2);
  xbar=meanXV;
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
    sx=(sigmaX-xi);
    xbar= sx / (nv-1);
    s = ((sigmaX2 - (xi * xi)) - (sx * sx / (nv-1))) / (nv-1);
    ////----
    qq=q;
    zj=0;
    for (i=0;i <= q;i++) {
      zj=zj+((1/qq) * (xn[i] - xbar));
    }
    
    zj=zj - ((1/qq) * (xi - xbar));
    
    Ii=((xi - xbar)/s) * zj;
    /////----
    Wi=1/qq;
    EI = -1/(nv-1);
    wi2=0;
    for (i=0;i < q;i++) {
      wi2=wi2+pow(1/qq,2);
    }
    wikh=0;
    for (i=0;i < q;i++) {
      wikh=wikh+pow(1/qq,2)+pow(1/qq,2);
    }
    VarI = (wi2 * (nv - b2) / (nv-1)) + (wikh * ((2*b2 - nv) / ((nv-1)*(nv-2)))) - (pow(Wi,2) / pow((nv-1),2));
    NUMERIC_POINTER(VECTOR_ELT(ans, 0))[j] = Ii;
    //NUMERIC_POINTER(VECTOR_ELT(ans, 1))[j] = VarI;
    NUMERIC_POINTER(VECTOR_ELT(ans, 1))[j] = (Ii - EI) / sqrt(VarI);
    //xans[j] = (Ii - EI) / sqrt(VarI);
  } 
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
/////---------
//-- Local Geary's c
SEXP localgeary(SEXP v, SEXP nc, SEXP nr, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, nnr, nnc, nrow, ncol, cellnr, n;
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
  
  int nv=0;
  double sigmaX=0,sigmaX2=0,s2=0, xbar, Eij,nn;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      sigmaX2=sigmaX2 + xv[i]*xv[i];
      nv++;
    } else {
      xans[i]=R_NaReal;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  nn=nv;
  xbar=sigmaX/nn;
  s2 = (sigmaX2 / nn) - (xbar*xbar);
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    row = (j / ncol) + 1;
    col = (j + 1) - ((row - 1) * ncol);
    
    Eij=0;
    for (i=0; i < ngb; i++) {
      nnr= row + xrr[i];
      nnc = col + xcc[i];
      
      if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
        cellnr = ((nnr - 1) * ncol) + nnc;
        if (!R_IsNA(xv[(cellnr-1)])) {
          Eij=Eij+pow(xi - xv[(cellnr-1)],2);
        }
      }
    }
    //----
    
    xans[j] = Eij / s2;
  } 
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
////----------


//-- LocalG and G*
SEXP GG(SEXP v, SEXP nc, SEXP nr, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, nnr, nnc, nrow, ncol, cellnr, n, q;
  double xi;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv;
  int *xrr, *xcc ;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = NEW_LIST(2));
  //PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  //xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  int nv=0;
  double sigmaX=0, sigmaX2=0, s2=0,xbar, sx, s, Wi, sxbar,ssx, sWi, G, G2, nn ;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      sigmaX2=sigmaX2 + xv[i]*xv[i];
      nv++;
    } else {
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[i] = R_NaReal;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  nn=nv;
  sxbar=sigmaX / nn;
  s2 = (sigmaX2 / nn) - (sxbar*sxbar);
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    row = (j / ncol) + 1;
    col = (j + 1) - ((row - 1) * ncol);
    
    double xn[ngb];
    q=-1;
    sx=0;
    for (i=0; i < ngb; i++) {
      nnr= row + xrr[i];
      nnc = col + xcc[i];
      
      if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
        cellnr = ((nnr - 1) * ncol) + nnc;
        if (!R_IsNA(xv[(cellnr-1)])) {
          q+=1;
          xn[q]=xv[(cellnr-1)];
          sx=sx+xn[q];
        }
      }
    }
    //----
    ssx=sx;
    sx=sx-xi;
    xbar=(sigmaX-xi) / (nn-1) ;
    Wi=q;
    sWi = q+1;
    s = sqrt(((sigmaX2 - (xi*xi)) / (nn-1)) - (xbar*xbar));
    G = (sx - (Wi * xbar)) / (s * sqrt((((nn-1)*Wi) - (Wi*Wi)) / (nn-2)));
    G2= (ssx - (sWi * sxbar)) / (sqrt(s2) * sqrt(((nn*sWi) - (sWi*sWi)) / (nn-1)));
    
    NUMERIC_POINTER(VECTOR_ELT(ans, 0))[j] = G;
    NUMERIC_POINTER(VECTOR_ELT(ans, 1))[j] = G2;
    
  } 
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
////---------####################----------------

SEXP localmoran_vector(SEXP v, SEXP nb) {
  int nProtected=0;
  int c, ngb, n, q;
  double xi, a;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv;
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = NEW_LIST(3));
  ++nProtected;
  
  xv=REAL(v);
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 2, NEW_NUMERIC(n));
  
  int nv=0;
  double sigmaX=0, sigmaX2=0, s2=0, s4=0, b2, xbar, sx, s,Ii, Wi, EI, VarI, zj, qq,wikh, wi2;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      sigmaX2=sigmaX2 + xv[i]*xv[i];
      nv++;
    } else {
      //xans[i]=R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[i] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 2))[i] = R_NaReal;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  xbar=sigmaX/nv;
  for (i=0;i < nv;i++) {
    j=xcells[i];
    s2=s2 + pow((xv[j] - xbar),2);
    s4=s4 + pow((xv[j] - xbar),4);
  }
  s2=s2/(double)nv;
  s4=s4/(double)nv;
  b2=s4/pow(s2,2);
  
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
    sx=(sigmaX-xi);
    xbar= sx / (nv-1);
    s = ((sigmaX2 - (xi * xi)) - (sx * sx / (nv-1))) / (nv-1);
    ////----
    qq=q+1;
    zj=0;
    Wi=1/qq;
    for (i=0;i <= q;i++) {
      zj=zj+(Wi * (xn[i] - xbar));
    }
    
    Ii=((xi - xbar)/s) * zj;
    /////----
    
    EI = -1/(nv-1);
    wi2=0;
    for (i=0;i <= q;i++) {
      wi2=wi2+pow(1/qq,2);
    }
    wikh=0;
    for (i=0;i <= q;i++) {
      wikh=wikh+pow(1/qq,2)+pow(1/qq,2);
    }
    VarI = (wi2 * (nv - b2) / (nv-1)) + (wikh * ((2*b2 - nv) / ((nv-1)*(nv-2)))) - (pow(Wi,2) / pow((nv-1),2));
    NUMERIC_POINTER(VECTOR_ELT(ans, 0))[j] = Ii;
    NUMERIC_POINTER(VECTOR_ELT(ans, 1))[j] = VarI;
    NUMERIC_POINTER(VECTOR_ELT(ans, 2))[j] = (Ii - EI) / sqrt(VarI);
  } 
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
/////---------
SEXP localgeary_vector(SEXP v, SEXP nb) {
  int nProtected=0;
  int c, ngb, n;
  double xi, a;
  R_len_t i, j;
  
  SEXP ans;
  double  *xv, *xans;
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  xans=REAL(ans);
  xv=REAL(v);
  
  int nv=0;
  double sigmaX=0,sigmaX2=0,s2=0, xbar, Eij,nn;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      sigmaX2=sigmaX2 + xv[i]*xv[i];
      nv++;
    } else {
      xans[i]=R_NaReal;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  nn=nv;
  xbar=sigmaX/nn;
  s2 = (sigmaX2 / nn) - (xbar*xbar);
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    
    ngb = length(VECTOR_ELT(nb,c));
    Eij=0;
    for (i=0;i < ngb;i++) {
      a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
      if (!R_IsNA(a)) {
        Eij=Eij+pow(xi - a,2);
      }
    }
    //----
    
    xans[j] = Eij / s2;
  }
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
////----------


//-- LocalG and G* (vector data)
SEXP GG_vector(SEXP v, SEXP nb) {
  int nProtected=0;
  int c, ngb, n, q;
  double xi, a;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv;
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = NEW_LIST(2));
  ++nProtected;
  
  xv=REAL(v);
  
  SET_VECTOR_ELT(ans, 0, NEW_NUMERIC(n));
  SET_VECTOR_ELT(ans, 1, NEW_NUMERIC(n));
  
  int nv=0;
  double sigmaX=0, sigmaX2=0, s2=0,xbar, sx, s, Wi, sxbar,ssx, sWi, G, G2, nn ;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      sigmaX=sigmaX + xv[i];
      sigmaX2=sigmaX2 + xv[i]*xv[i];
      nv++;
    } else {
      NUMERIC_POINTER(VECTOR_ELT(ans, 0))[i] = R_NaReal;
      NUMERIC_POINTER(VECTOR_ELT(ans, 1))[i] = R_NaReal;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  nn=nv;
  sxbar=sigmaX / nn;
  s2 = (sigmaX2 / nn) - (sxbar*sxbar);
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    
    ngb = length(VECTOR_ELT(nb,c));
    
    double xn[ngb];
    q=-1;
    sx=0;
    for (i=0;i < ngb;i++) {
      a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
      if (!R_IsNA(a)) {
        q+=1;
        xn[q]=a;
        sx=sx+a;
      }
    }
    //----
    ssx=sx+xi;
    xbar=(sigmaX-xi) / (nn-1) ;
    Wi=q+1;
    sWi = q+2;
    s = sqrt(((sigmaX2 - (xi*xi)) / (nn-1)) - (xbar*xbar));
    G = (sx - (Wi * xbar)) / (s * sqrt((((nn-1)*Wi) - (Wi*Wi)) / (nn-2)));
    G2= (ssx - (sWi * sxbar)) / (sqrt(s2) * sqrt(((nn*sWi) - (sWi*sWi)) / (nn-1)));
    
    NUMERIC_POINTER(VECTOR_ELT(ans, 0))[j] = G;
    NUMERIC_POINTER(VECTOR_ELT(ans, 1))[j] = G2;
  } 
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
