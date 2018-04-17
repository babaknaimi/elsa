
/*
 * Babak Naimi, July 2016
 * naimi.b@gmail.com
 * 
 * Most parts of this code is copied from spdep package (written by Roger Bivand) 
 * with the following information:
 * 
 *  based on code taken from:
 *  class/class.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-9
 *  and written by Roger Bivand (C) 2001-2014
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
#define ROFFSET 1


void gcdist(double *lon1, double *lon2, double *lat1, double *lat2, 
            double *dist) {
  
  double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
  double w, R, a, f, D, H1, H2;
  double lat1R, lat2R, lon1R, lon2R, DE2RA;
  
  DE2RA = M_PI/180;
  a = 6378.137;              /* WGS-84 equatorial radius in km */
f = 1.0/298.257223563;     /* WGS-84 ellipsoid flattening factor */

lat1R = lat1[0]*DE2RA;
lat2R = lat2[0]*DE2RA;
lon1R = lon1[0]*DE2RA;
lon2R = lon2[0]*DE2RA;

F = ( lat1R + lat2R )/2.0;
G = ( lat1R - lat2R )/2.0;
L = ( lon1R - lon2R )/2.0;

sinG2 = R_pow_di( sin( G ), 2 );
cosG2 = R_pow_di( cos( G ), 2 );
sinF2 = R_pow_di( sin( F ), 2 );
cosF2 = R_pow_di( cos( F ), 2 );
sinL2 = R_pow_di( sin( L ), 2 );
cosL2 = R_pow_di( cos( L ), 2 );

S = sinG2*cosL2 + cosF2*sinL2;
C = cosG2*cosL2 + sinF2*sinL2;

w = atan( sqrt( S/C ) );
R = sqrt( S*C )/w;

D = 2*w*a;
H1 = ( 3*R - 1 )/( 2*C );
H2 = ( 3*R + 2 )/( 2*S );

dist[0] = D*( 1 + f*H1*sinF2*cosG2 - f*H2*cosF2*sinG2 ); 

}



int spOverlapC(double bbi1, double bbi2, double bbi3, double bbi4, double bbj1, double bbj2, double bbj3, double bbj4) {
  
  int overlap=1;
  
  if ((bbi1>bbj3) || (bbi2>bbj4) || 
      (bbi3<bbj1) || (bbi4<bbj2) ) {
    overlap = 0;
  }
  
  return(overlap);
}

int polypolyC(double *px1, double *py1, int n1, double *px2, double *py2,
              int n2, double sn, int crit) {
  int i, j, k=0;
  double dist;
  double x1, x2, y1, y2, xd, yd;
  
  for (i=0; (i < n1) && (k < crit); i++) {
    x1 = px1[i];
    y1 = py1[i];
    for (j=0; (j < n2) && (k < crit); j++) {
      x2 = px2[j];
      y2 = py2[j];
      xd = x1-x2;
      if (fabs(xd)>sn) { continue; }
      yd = y1-y2;
      if (fabs(yd)>sn) { continue; }
      dist = hypot(xd, yd);
      if (dist <= sn) {
        k++;
      }
    }
  }
  
  return(k);
}




SEXP poly_loop2(SEXP n, SEXP i_findInBox, SEXP bb, SEXP pl, SEXP nrs,
                SEXP dsnap, SEXP criterion, SEXP nfIBB) {
  
  int nn = INTEGER_POINTER(n)[0];
  int crit = INTEGER_POINTER(criterion)[0];
  /*    int Scale = INTEGER_POINTER(scale)[0];*/
  int uBound = (int) INTEGER_POINTER(nfIBB)[0]*2;
  int i, j, jj, li, pc = 0;
  int ii = 0;
  int *card, *icard, *is, *jjs, *NRS, *cNRS;
  double *bb1, *bb2, *bb3, *bb4, *plx, *ply;
  double Dsnap = NUMERIC_POINTER(dsnap)[0];
  
  //    struct bbcontainer *bbs;
  
  SEXP ans;
  
  int jhit, khit, nrsi, nrsj;
  
  int xx, yy, zz, ww;
  
  card = (int *) R_alloc((size_t) nn, sizeof(int));
  icard = (int *) R_alloc((size_t) nn, sizeof(int));
  is = (int *) R_alloc((size_t) uBound, sizeof(int));
  jjs = (int *) R_alloc((size_t) uBound, sizeof(int));
  bb1 = (double *) R_alloc((size_t) nn, sizeof(double));
  bb2 = (double *) R_alloc((size_t) nn, sizeof(double));
  bb3 = (double *) R_alloc((size_t) nn, sizeof(double));
  bb4 = (double *) R_alloc((size_t) nn, sizeof(double));
  NRS = (int *) R_alloc((size_t) nn, sizeof(int));
  cNRS = (int *) R_alloc((size_t) nn, sizeof(int));
  
  for (i=0, li=0; i<nn; i++) {
    card[i] = 0;
    icard[i] = 0;
    bb1[i] = NUMERIC_POINTER(bb)[i];
    bb2[i] = NUMERIC_POINTER(bb)[i+(1*nn)];
    bb3[i] = NUMERIC_POINTER(bb)[i+(2*nn)];
    bb4[i] = NUMERIC_POINTER(bb)[i+(3*nn)];
    NRS[i] = INTEGER_POINTER(nrs)[i];
    li += NRS[i];
  }
  
  for (i=0; i<nn; i++) {
    if (i == 0) cNRS[i] = 0;
    else cNRS[i] = NRS[i-1] + cNRS[i-1];
  }
  
  for (i=0; i<uBound; i++) {
    is[i] = 0;
    jjs[i] = 0;
  }
  
  plx = (double *) R_alloc((size_t) li, sizeof(double));
  ply = (double *) R_alloc((size_t) li, sizeof(double));
  
  for (i=0, jj=0; i<nn; i++) {
    nrsi = NRS[i];
    for (j=0; j<nrsi; j++) {
      plx[jj] = NUMERIC_POINTER(VECTOR_ELT(pl, i))[j];
      ply[jj] = NUMERIC_POINTER(VECTOR_ELT(pl, i))[j+nrsi];
      jj++;
      /*            if (i < (nn-1) && jj == li) error("polygon memory overflow");*/
    }
  }
  
  for (i=0; i<(nn-1); i++) {
    li = length(VECTOR_ELT(i_findInBox, i));
    nrsi = NRS[i];
    for (j=0; j<li; j++) {
      jj = INTEGER_POINTER(VECTOR_ELT(i_findInBox, i))[j] - ROFFSET;
      jhit = spOverlapC(bb1[i], bb2[i], bb3[i], bb4[i], bb1[jj],
                        bb2[jj], bb3[jj], bb4[jj]);
      if (jhit > 0) {
        khit = 0;
        nrsj = NRS[jj];
        if (nrsi > 0 && nrsj > 0){
          khit = polypolyC(&plx[cNRS[i]], &ply[cNRS[i]], nrsi,
                           &plx[cNRS[jj]], &ply[cNRS[jj]], nrsj, Dsnap, crit+1L);
        }
        if (khit > crit) {
          card[i]++;
          card[jj]++;
          is[ii] = i;
          jjs[ii] = jj;
          ii++;
          /*                    if (ii == uBound) error("memory error, scale problem");*/
        }
      }
    }
  }
  
  PROTECT(ans = NEW_LIST(nn)); pc++;
  
  for (i=0; i<nn; i++) {
    if (card[i] == 0) {
      SET_VECTOR_ELT(ans, i, NEW_INTEGER(1));
      INTEGER_POINTER(VECTOR_ELT(ans, i))[0] = 0;
    } else {
      SET_VECTOR_ELT(ans, i, NEW_INTEGER(card[i]));
    }
  }
  
  for (i=0; i<ii; i++) {
    xx = is[i];
    yy = jjs[i];
    zz = icard[yy];
    ww = icard[xx];
    /*        if (zz == card[yy]) error("memory error, overflow");
    if (ww == card[xx]) error("memory error, overflow");*/
    INTEGER_POINTER(VECTOR_ELT(ans, yy))[zz] = xx + ROFFSET;
    INTEGER_POINTER(VECTOR_ELT(ans, xx))[ww] = yy + ROFFSET;
    icard[yy]++;
    icard[xx]++;
  }
  
  for (i=0; i<nn; i++) {
    if ((li = length(VECTOR_ELT(ans, i))) > 1) {
      for (j=0; j<li; j++)
        icard[j] = INTEGER_POINTER(VECTOR_ELT(ans, i))[j];
      R_isort(icard, li);
      for (j=0; j<li; j++)
        INTEGER_POINTER(VECTOR_ELT(ans, i))[j] = icard[j];
    }
  }
  
  UNPROTECT(pc);
  return(ans);
}

/////////////////------------

SEXP dnn(SEXP x, SEXP y, SEXP d1, SEXP d2, SEXP lonlat) {
  int nProtected=0;
  int q, n, ll;
  double x1[1], y1[1], x2[1], y2[1], gc[1];
  double dist, xd1, xd2;
  ll = INTEGER_POINTER(lonlat)[0];
  
  R_len_t i, j, k;
  
  SEXP ans;
  double *xx, *yy;
  
  n=length(x);
  
  PROTECT(ans = NEW_LIST(n));
  ++nProtected;
  
  PROTECT(x = coerceVector(x, REALSXP));
  ++nProtected;
  PROTECT(y = coerceVector(y, REALSXP));
  ++nProtected;
  
  xd1 = NUMERIC_POINTER(d1)[0];
  xd2 = NUMERIC_POINTER(d2)[0];
  
  xx=REAL(x);
  yy=REAL(y);
  
  for (i=0;i < n;i++) {
    x1[0]=xx[i];
    y1[0]=yy[i];
    
    int xn[n];
    q=0;
    for (j=0;j < n;j++) {
      if (j == i) continue;
      x2[0]=xx[j];
      y2[0]=yy[j];
      if (ll == 0) dist = hypot((x1[0]-x2[0]), (y1[0]-y2[0]));
      else {
        gcdist(x1, x2, y1, y2, gc);
        dist = gc[0];
      }
      
      if (dist <= xd2 && dist >= xd1) {
        xn[q] = j;
        q=q + 1;
      }
    }
    
    if (q > 0) {
      SET_VECTOR_ELT(ans, i, NEW_INTEGER(q));
      for (k = 0; k < q; k++) {
        INTEGER_POINTER(VECTOR_ELT(ans, i))[k] = xn[k] + 1;
      }
    } else {
      SET_VECTOR_ELT(ans,i,R_NilValue);
    }
  }
  UNPROTECT(nProtected);
  return(ans);
  
}