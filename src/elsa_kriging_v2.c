/*
 * Babak Naimi, May 2023 (last update: Oct. 2025)
 * version: 1.3
 * naimi.b@gmail.com
 * 
 */



#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* ---------- kernels (stationary similarity) ---------- */
static inline double S_exp(double h, double r){ return exp(-h / fmax(r, DBL_EPSILON)); }
static inline double S_mat32(double h, double r){ double a = sqrt(3.0)*h/r; return (1.0+a)*exp(-a); }
static inline double S_mat52(double h, double r){ double a = sqrt(5.0)*h/r; return (1.0+a+(5.0*h*h)/(3.0*r*r))*exp(-a); }
static inline double S_sph(double h, double r){
  if (h >= r) return 0.0; double x = h/r; return 1.0 - 1.5*x + 0.5*x*x*x;
}

/* ---------- Wendland C^2 taper (2D PD): T(r)=(1-r)^4_+(1+4r) ---------- */
static inline double wendland_c2(double d_over_theta){
  if (d_over_theta >= 1.0) return 0.0;
  double r = 1.0 - d_over_theta;
  return r*r*r*r * (1.0 + 4.0*(1.0 - r));
}

/* ---------- PS (Paciorek–Schervish) isotropic nonstationary Gaussian ----------
 Sigma(x)=l(x)^2 I; d-dim prefactor reduces to ((2 l_i l_j)/(l_i^2 + l_j^2))^(d/2)
 exponent uses (x-x')^T ((Sigma_i + Sigma_j)/2)^(-1) (x-x') / 2
 --------------------------------------------------------------------------- */
static inline double S_ps_iso(double dx, double dy, double li, double lj, int dim){
  double li2 = li*li, lj2 = lj*lj;
  double s2 = 0.5 * (li2 + lj2);
  double q = (dx*dx + dy*dy) / fmax(s2, DBL_EPSILON);
  double pref = pow( (2.0*li*lj)/fmax(li2 + lj2, DBL_EPSILON), 0.5*dim );
  return pref * exp(-0.5*q);
}

/* ---------- small utilities ---------- */
static inline double euclid(double x1,double y1,double x2,double y2){
  double dx = x1-x2, dy = y1-y2; return sqrt(dx*dx + dy*dy);
}

/* sort a small array of doubles ascending (for trimmed mean) */
static void sort_dbl(double *a, int n){
  for(int i=1;i<n;i++){
    double key=a[i]; int j=i-1;
    while(j>=0 && a[j]>key){ a[j+1]=a[j]; j--; }
    a[j+1]=key;
  }
}

/* WLS fit of range with penalty; weights w[i] (counts), lags h[i], targets S[i] */
static double fit_r_grid_pen(const double *h, const double *S, const double *w, int L,
                             int type, double r_min, double r_max, double pen_lambda){
  const int G = 40;
  double best_r = r_min, best_obj = HUGE_VAL;
  for(int g=0; g<G; ++g){
    double r = r_min + (r_max - r_min) * ( (double)g / (double)(G-1) );
    double sse = 0.0;
    for(int i=0;i<L;++i){
      double Sh;
      switch(type){
      case 0: Sh = S_exp(h[i], r); break;
      case 1: Sh = S_mat32(h[i], r); break;
      case 2: Sh = S_mat52(h[i], r); break;
      default: Sh = S_sph(h[i], r); break;
      }
      double e = S[i] - Sh;
      sse += (w? w[i] : 1.0) * e*e;
    }
    double obj = sse + pen_lambda * (r / fmax(r_max, DBL_EPSILON)) * (r / fmax(r_max, DBL_EPSILON));
    if (obj < best_obj){ best_obj = obj; best_r = r; }
  }
  return best_r;
}

/* ordinary kriging solve using dgesv on (k+1)x(k+1) system */
static int ok_solve(const double *C, const double *c, int k,
                    double sigma2_plus_nugget,
                    const double *y,
                    double *pred, double *var){
  const int nsys = k + 1;
  double *A = (double*) R_alloc(nsys*nsys, sizeof(double));
  double *rhs = (double*) R_alloc(nsys, sizeof(double));
  for(int i=0;i<k;i++){
    for(int j=0;j<k;j++) A[i + j*nsys] = C[i + j*k];
    A[i + k*nsys] = 1.0; A[k + i*nsys] = 1.0;
  }
  A[k + k*nsys] = 0.0;
  for(int i=0;i<k;i++) rhs[i] = c[i];
  rhs[k] = 1.0;
  
  int *ipiv = (int*) R_alloc(nsys, sizeof(int));
  int info = 0, one = 1;
  F77_CALL(dgesv)(&nsys, &one, A, &nsys, ipiv, rhs, &nsys, &info);
  if (info != 0) return info;
  
  double zhat = 0.0, wc = 0.0;
  for(int i=0;i<k;i++){ zhat += rhs[i]*y[i]; wc += rhs[i]*c[i]; }
  *pred = zhat;
  *var  = fmax(0.0, sigma2_plus_nugget - wc - rhs[k]);
  return 0;
}

/* ---------- .Call entry ---------- */
SEXP C_elsa_local_krige_v2(SEXP Rx, SEXP Ry, SEXP Rz,
                           SEXP Rpx, SEXP Rpy,
                           SEXP Rlags,
                           SEXP RELSA_train,      /* n x L, col-major from R */
SEXP Rls_train,        /* length-scale per training pt or NULL */
SEXP RKentro, SEXP RKkrige,
SEXP Rsigma2,
SEXP Rrmin, SEXP Rrmax,
SEXP Rtau_min, SEXP Rtau_max,
SEXP Rq,
SEXP Rkernel_code,     /* -1 auto; 0..3 fixed; 99 PS */
SEXP Rkernel_menu,     /* int codes 0..3 to consider in auto */
SEXP Rmodel_avg,       /* 0/1 */
SEXP Rtrim_prop,       /* 0..0.4 typical */
SEXP Rlambda_pen,
SEXP Rtaper_theta,     /* 0 off; else >0 */
SEXP Rtaper_mode,      /* 1 xrange, 2 fixed */
SEXP Rthreads){
  /* --- unpack --- */
  const double *x = REAL(Rx), *y = REAL(Ry), *z = REAL(Rz);
  const double *px = REAL(Rpx), *py = REAL(Rpy);
  const double *lags = REAL(Rlags);
  const double *ELSA = REAL(RELSA_train);
  const double *ls_train = (Rls_train == R_NilValue ? NULL : REAL(Rls_train));
  const double sigma2 = REAL(Rsigma2)[0];
  const double rmin = REAL(Rrmin)[0], rmax = REAL(Rrmax)[0];
  const double tmin = REAL(Rtau_min)[0], tmax = REAL(Rtau_max)[0];
  const double q = REAL(Rq)[0];
  const int kernel_code = asInteger(Rkernel_code);
  const int model_avg = asInteger(Rmodel_avg);
  const double trim_prop = REAL(Rtrim_prop)[0];
  const double pen_lambda = REAL(Rlambda_pen)[0];
  const double taper_theta_arg = REAL(Rtaper_theta)[0];
  const int taper_mode = asInteger(Rtaper_mode); /* 1 xrange, 2 fixed */
  
  const int n = LENGTH(Rx);
  const int m = LENGTH(Rpx);
  const int L = LENGTH(Rlags);
  const int K_entro = asInteger(RKentro);
  const int K_krige = asInteger(RKkrige);
  const int nmenu = LENGTH(Rkernel_menu);
  int menu[4] = {0,1,2,3};
  for(int i=0;i<nmenu && i<4; ++i) menu[i] = INTEGER(Rkernel_menu)[i]-1;
  
  /* outputs */
  SEXP Rpred = PROTECT(allocVector(REALSXP, m));
  SEXP Rvar  = PROTECT(allocVector(REALSXP, m));
  double *pred = REAL(Rpred), *vout = REAL(Rvar);
  
  /* scratch */
  double *Evals = (double*) R_alloc(K_entro, sizeof(double)); /* per-lag neighbor ELSA values */
  double *Sbar  = (double*) R_alloc(L, sizeof(double));       /* local similarity 1-E */
  int    *ring_counts = (int*) R_alloc(L, sizeof(int));
  double *wlag = (double*) R_alloc(L, sizeof(double));
  
  int *idx = (int*) R_alloc(n, sizeof(int));
  double *di = (double*) R_alloc(n, sizeof(double));
  
  double *C = (double*) R_alloc(K_krige*K_krige, sizeof(double));
  double *c = (double*) R_alloc(K_krige, sizeof(double));
  double *Xi = (double*) R_alloc(K_krige*2, sizeof(double));
  double *yi = (double*) R_alloc(K_krige, sizeof(double));
  int threads = asInteger(Rthreads);
#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#endif
  
#pragma omp parallel for schedule(static)
  for(int t=0; t<m; ++t){
    /* distances & partial sort (selection) for first max(K_entro, K_krige) */
    for(int j=0;j<n;++j){ di[j] = euclid(px[t], py[t], x[j], y[j]); idx[j] = j; }
    /* simple selection for small n; for large n replace with partial quickselect */
    for(int a=0;a<n-1;++a){
      int min = a;
      for(int b=a+1;b<n;++b) if (di[b] < di[min]) min=b;
      double td=di[a]; di[a]=di[min]; di[min]=td;
      int ti=idx[a]; idx[a]=idx[min]; idx[min]=ti;
    }
    
    const int kE = (K_entro < n? K_entro : n);
    const int kK = (K_krige < n? K_krige : n);
    
    /* Build ring counts from the K_entro neighbors */
    for(int i=0;i<L;++i) ring_counts[i] = 0;
    for(int j=0;j<kE;++j){
      double d = di[j];
      int ring = 0;
      while(ring < L && d > lags[ring]) ring++;
      if (ring >= L) ring = L-1;
      ring_counts[ring]++;
    }
    for(int i=0;i<L;++i) wlag[i] = (double) (ring_counts[i] > 0 ? ring_counts[i] : 1); /* avoid zero weight */
    
    /* Local entrogram Sbar(h) = 1 - trimmed_mean(ELSA over kE neighbors at each lag) */
    for(int i=0;i<L;++i){
      int nvals = 0;
      for(int j=0;j<kE;++j){
        double v = ELSA[idx[j] + i*n];
        if (!ISNAN(v)) Evals[nvals++] = v;
      }
      if (nvals == 0){ Sbar[i] = 0.5; continue; }
      /* trimmed mean */
      sort_dbl(Evals, nvals);
      int a = (int) floor(trim_prop * nvals);
      int b = (int) ceil( (1.0 - trim_prop) * nvals );
      if (b <= a) { a = 0; b = nvals; }
      double sum=0.0; int cnt=0;
      for(int k=a;k<b;k++){ sum += Evals[k]; cnt++; }
      double Ebar = (cnt? sum/cnt : 0.5);
      if (Ebar < 0.0) Ebar = 0.0; if (Ebar > 1.0) Ebar = 1.0;
      Sbar[i] = 1.0 - Ebar;
    }
    
    /* Nugget from inner lag */
    double tau2 = tmin + (tmax - tmin) * pow(1.0 - Sbar[0], q);
    
    /* stationaries or PS? */
    double zhat=NA_REAL, vhat=NA_REAL;
    
    /* --- PS kernel path (nonstationary, isotropic) --- */
    if (kernel_code == 99){
      /* 0) Guard: if very few neighbors in the two nearest rings, fallback to stationary exp */
      int near_ct = (ring_counts[0] + (L>1 ? ring_counts[1] : 0));
      if (near_ct < 5) { /* use exp stationary with WLS-fitted range */
      int tp = 0; /* exp */
      double rhat = fit_r_grid_pen(lags, Sbar, wlag, L, tp, rmin, rmax, pen_lambda);
      /* neighborhood */
      for(int j=0;j<kK;++j){ Xi[j]=x[idx[j]]; Xi[j+kK]=y[idx[j]]; yi[j]=z[idx[j]]; }
      for(int i=0;i<kK;++i){
        for(int j=0;j<kK;++j){
          double dij = euclid(Xi[i], Xi[i+kK], Xi[j], Xi[j+kK]);
          C[i + j*kK] = sigma2 * S_exp(dij, rhat);
        }
        C[i + i*kK] += tau2;
        double di0 = euclid(Xi[i], Xi[i+kK], px[t], py[t]);
        c[i] = sigma2 * S_exp(di0, rhat);
      }
      if (ok_solve(C, c, kK, sigma2 + tau2, yi, &zhat, &vhat) != 0){ zhat=NA_REAL; vhat=NA_REAL; }
      pred[t]=zhat; vout[t]=vhat; continue;
      }
      
      /* 1) Build l0 by trimmed, distance-weighted mean of *smoothed* ls_train */
      double l0;
      {
        int nvals = kE;
        for(int j=0;j<kE;++j) Evals[j] = (ls_train ? ls_train[idx[j]] : rmin); /* reuse scratch */
      sort_dbl(Evals, nvals);
      int a = (int)floor(0.10*nvals), b = (int)ceil(0.90*nvals);
      /* recompute weights on the fly using distances */
      double num=0.0, den=0.0;
      for(int j=0;j<kE;++j){
        if (j<a || j>=b) continue; /* trimmed */
      double w = 1.0 / fmax(di[j], 1e-6); /* inverse distance */
      num += w * Evals[j]; den += w;
      }
      l0 = (den>0? num/den : rmin);
      if (l0 < rmin) l0 = rmin;
      if (l0 > rmax) l0 = rmax;
      }
      const double cap = 3.0; /* cap ratio li/l0 to [1/3, 3] */
      
      /* 2) Neighborhood */
      for(int j=0;j<kK;++j){ Xi[j]=x[idx[j]]; Xi[j+kK]=y[idx[j]]; yi[j]=z[idx[j]]; }
      
      /* 3) Fill C and c using PS with Gaussian base; taper OFF in PS */
      for(int i=0;i<kK;++i){
        double li = (ls_train ? ls_train[idx[i]] : l0);
        /* cap ratio */
        if (li < l0/cap) li = l0/cap;
        if (li > l0*cap) li = l0*cap;
        
        for(int j=0;j<kK;++j){
          double lj = (ls_train ? ls_train[idx[j]] : l0);
          if (lj < l0/cap) lj = l0/cap;
          if (lj > l0*cap) lj = l0*cap;
          
          double sim = S_ps_iso(Xi[i]-Xi[j], Xi[i+kK]-Xi[j+kK], li, lj, 2);
          C[i + j*kK] = sigma2 * sim; /* no taper here */
        }
        C[i + i*kK] += tau2;
        
        double li0 = li; /* reuse last li for cross-cov; benign approximation */
        double sim0 = S_ps_iso(Xi[i]-px[t], Xi[i+kK]-py[t], li0, l0, 2);
        c[i] = sigma2 * sim0;
      }
      if (ok_solve(C, c, kK, sigma2 + tau2, yi, &zhat, &vhat) != 0){ zhat=NA_REAL; vhat=NA_REAL; }
      pred[t]=zhat; vout[t]=vhat; continue;
    }
    
    
    /* stationary kernels: fixed type or auto (+ optional model averaging) */
    int candidates[4]; int M = 0;
    if (kernel_code >= 0 && kernel_code <= 3){ candidates[0]=kernel_code; M=1; }
    else { for(int mm=0; mm<nmenu && mm<4; ++mm) candidates[mm]=menu[mm], M=nmenu; }
    
    double best_aic = HUGE_VAL; int best_type=candidates[0]; double best_r=rmin;
    double zbar = 0.0, vbar = 0.0, wsum = 0.0;
    
    for(int mtype_idx=0; mtype_idx<M; ++mtype_idx){
      int tp = candidates[mtype_idx];
      /* fit local range with WLS (counts) + penalty */
      double rhat = fit_r_grid_pen(lags, Sbar, wlag, L, tp, rmin, rmax, pen_lambda);
      double theta = 0.0;
      if (taper_theta_arg > 0){
        theta = (taper_mode==1 ? taper_theta_arg * rhat : taper_theta_arg);
      }
      /* SSE and AIC for weights */
      double sse = 0.0;
      for(int i=0;i<L;++i){
        double Sh;
        switch(tp){ case 0: Sh=S_exp(lags[i],rhat); break; case 1: Sh=S_mat32(lags[i],rhat); break;
        case 2: Sh=S_mat52(lags[i],rhat); break; default: Sh=S_sph(lags[i],rhat); }
        double e=Sbar[i]-Sh; sse += wlag[i]*e*e;
      }
      double aic = L*log(sse/fmax(L,1)) + 2.0;
      
      if (!model_avg){
        if (aic < best_aic){ best_aic = aic; best_type = tp; best_r = rhat; }
        if (mtype_idx < M-1) continue; /* delay solve to once for best */
      }
      
      /* neighborhood */
      for(int j=0;j<kK;++j){
        Xi[j] = x[idx[j]];
        Xi[j + kK] = y[idx[j]];
        yi[j] = z[idx[j]];
      }
      /* fill C & c */
      for(int i=0;i<kK;++i){
        for(int j=0;j<kK;++j){
          double dij = euclid(Xi[i], Xi[i+kK], Xi[j], Xi[j+kK]);
          double sim;
          switch(tp){ case 0: sim=S_exp(dij,rhat); break; case 1: sim=S_mat32(dij,rhat); break;
          case 2: sim=S_mat52(dij,rhat); break; default: sim=S_sph(dij,rhat); }
          if (taper_theta_arg > 0) sim *= wendland_c2(dij / fmax(theta, DBL_EPSILON));
          C[i + j*kK] = sigma2 * sim;
        }
        C[i + i*kK] += tau2;
        double di0 = euclid(Xi[i], Xi[i+kK], px[t], py[t]);
        double sim0;
        switch(tp){ case 0: sim0=S_exp(di0,rhat); break; case 1: sim0=S_mat32(di0,rhat); break;
        case 2: sim0=S_mat52(di0,rhat); break; default: sim0=S_sph(di0,rhat); }
        if (taper_theta_arg > 0) sim0 *= wendland_c2(di0 / fmax(theta, DBL_EPSILON));
        c[i] = sigma2 * sim0;
      }
      double zhat_m, vhat_m;
      if (ok_solve(C, c, kK, sigma2 + tau2, yi, &zhat_m, &vhat_m) != 0){ zhat_m=NA_REAL; vhat_m=NA_REAL; }
      
      if (!model_avg){
        pred[t] = zhat_m; vout[t] = vhat_m;
      } else {
        /* AIC weight */
        double w = exp(-0.5*(aic)); /* we'll normalize across types */
        zbar += w * zhat_m;
        vbar += w * (vhat_m + zhat_m*zhat_m);
        wsum += w;
        if (mtype_idx == M-1){
          if (wsum > 0){
            zbar /= wsum; vbar = vbar/wsum - zbar*zbar;
          } else { zbar = NA_REAL; vbar = NA_REAL; }
          pred[t] = zbar; vout[t] = fmax(0.0, vbar);
        }
      }
    } /* types loop */
        
        if (!model_avg && kernel_code == -1){
          /* after choosing best, if we didn't compute yet: redo best type solve */
          /* already done inside branch above to avoid duplication */
          (void)best_type; (void)best_r;
        }
  }
  
  SEXP ans = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ans, 0, Rpred);
  SET_VECTOR_ELT(ans, 1, Rvar);
  SEXP nm = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(nm, 0, mkChar("pred"));
  SET_STRING_ELT(nm, 1, mkChar("var"));
  setAttrib(ans, R_NamesSymbol, nm);
  
  UNPROTECT(4);
  return ans;
}
