# Author: Babak Naimi, naimi.b@gmail.com
# Date :  May 2023
# Last Update :  Oct. 2025
# Version 1.3
# Licence GPL v3 
#--------------------------

.fit_lengthscale_row <- function(S_row, lags, w = NULL,
                                 r_lo_min = 1e-6, r_hi_fac = 10) {
  ok <- is.finite(S_row)
  h  <- lags[ok]; S <- S_row[ok]
  if (length(h) < 2L) return(NA_real_)   # let caller fill fallback
  
  # clamp similarity into (0,1) to avoid log/exp edge issues
  S <- pmin(pmax(S, 1e-6), 1 - 1e-6)
  if (is.null(w)) w <- rep(1, length(h)) else w <- w[ok]
  
  # strictly positive bounds
  hpos <- h[h > 0]
  r_lo <- if (length(hpos)) min(hpos)/10 else r_lo_min
  r_hi <- max(h) * r_hi_fac
  
  obj <- function(r) {
    r <- max(r, r_lo_min)
    Sh <- exp(-h / r)         # exponential similarity
    sum(w * (S - Sh)^2)
  }
  stats::optimize(obj, c(r_lo, r_hi))$minimum
}

####
# Gaussian kernel smoother for log-lengthscales
.smooth_log_ls <- function(XY, ls, bw = NULL, min_n = 10L) {
  XY <- as.matrix(XY); stopifnot(ncol(XY) == 2L)
  n  <- nrow(XY)
  ls <- as.numeric(ls)
  if (length(ls) != n) stop(sprintf("length(ls) [%d] != nrow(XY) [%d]", length(ls), n))
  logls <- log(pmax(ls, .Machine$double.eps))
  
  if (is.null(bw)) {
    d <- as.matrix(dist(XY)); diag(d) <- Inf
    bw <- 2 * stats::median(apply(d, 1, min), na.rm = TRUE)
    if (!is.finite(bw) || bw <= 0) bw <- 1
  }
  
  out <- numeric(n)
  for (i in seq_len(n)) {
    dx <- XY[,1] - XY[i,1]; dy <- XY[,2] - XY[i,2]
    r2 <- dx*dx + dy*dy
    w  <- exp(-r2 / (2*bw*bw))        # length n
    if (sum(w > 1e-6) < min_n) {     # ensure minimum effective n
      ord <- order(r2); w[ord[seq_len(min_n)]] <- 1
    }
    w  <- w / sum(w)
    out[i] <- as.numeric(crossprod(w, logls))  # same length dot product
  }
  pmax(exp(out), .Machine$double.eps)
}


#------


if (!isGeneric("eKrige")) {
  setGeneric("eKrige", function(x, zcol, width=NULL, cutoff=NULL, longlat,Lmax=20,
                                K_entro = 60,
                                ps_lengthscale = TRUE,
                                kernel_menu = c("exp","mat32","mat52","sph"))
    standardGeneric("eKrige"))
}



setMethod('eKrige', signature(x='SpatVector'), 
          function(x, zcol, width=NULL, cutoff=NULL, longlat,Lmax=20,
                   K_entro = 60,
                   ps_lengthscale = TRUE,
                   kernel_menu = c("exp","mat32","mat52","sph")) {
            
            if (missing(zcol)) {
              if (ncol(x) > 1) stop("zcol should be specified!")
              else zcol <- 1
            } else if (is.character(zcol)) {
              w <- which(names(x) == zcol[1])
              if (w == 0) stop('the specified variable in zcol does not exist in the data')
              zcol <- w
            } else if (is.numeric(zcol)) {
              zcol <- zcol[1]
              if (zcol > ncol(x)) stop('the zcol number is greater than the number of columns in data!')
            } else stop("zcol should be a character or a number!")
            #-----------
            
            XY <- as.matrix(crds(x))
            z  <- as.numeric(terra::values(x)[, zcol])
            #-----
            dd <- neighd(x, 0, 2000000, longlat = longlat)
            r  <- range(unlist(dd), na.rm = TRUE)
            d     <- if (is.null(width)) max(r[1]*2, .Machine$double.eps) else width
            dmax  <- if (is.null(cutoff)) r[2]/3 else cutoff
            lags  <- seq(d, dmax, by = d)
            if (length(lags) > Lmax) lags <- lags[seq_len(Lmax)]
            L <- length(lags); n <- nrow(XY)
            
            # training-point ELSA at each lag (neighbor-averaged):
            ELSA_train <- matrix(NA_real_, n, L, dimnames = list(NULL, paste0("lag_", seq_len(L))))
            for (i in seq_len(L)) {
              dn <- dneigh(x, 0, lags[i], longlat = longlat)
              ei <- elsa(x, d = dn, zcol = zcol, categorical = is.factor(z))
              evec <- as.data.frame(ei)[, "ELSA"]
              nb <- dn@neighbours
              for (j in seq_along(nb)) if (!is.null(nb[[j]])) {
                ELSA_train[j, i] <- mean(evec[nb[[j]]], na.rm = TRUE)
              }
            }
            #-------
            
            ls_train <- NULL
            
            if (isTRUE(ps_lengthscale)) {
              S_train <- 1 - pmin(pmax(ELSA_train, 0), 1)
              wlag    <- c(lags[1], diff(lags))
              .fit_r <- function(S_row) {
                ok <- is.finite(S_row); h <- lags[ok]; S <- S_row[ok]
                if (length(h) < 2) return(NA_real_)
                S <- pmin(pmax(S, 1e-6), 1 - 1e-6)
                obj <- function(r) { Sh <- exp(-h / max(r, 1e-6)); sum(wlag[ok] * (S - Sh)^2) }
                optimize(obj, c(max(min(h[h>0])/10, 1e-6), max(h)*10))$minimum
              }
              ls_raw <- apply(S_train, 1, .fit_r)
              if (anyNA(ls_raw)) {
                med <- if (any(is.finite(ls_raw))) stats::median(ls_raw, na.rm = TRUE) else 2 * width
                ls_raw[!is.finite(ls_raw)] <- med
              }
              
              
              ls_raw <- as.numeric(ls_raw)
              
              ls_smooth <- .smooth_log_ls(XY, ls_raw, bw = 2.5 * d, min_n = 12L)
              
              # store the SMOOTHED ls in the object (PS branch will use it)
              object_ls <- ls_smooth
            } else {
              object_ls <- NULL
            }
            
            structure(list(
              coords = XY, z = z,
              lags = lags, width = d, cutoff = dmax,
              ELSA_train = ELSA_train,
              ls_train   = object_ls,           # <-- smoothed if available
              K_entro = as.integer(K_entro),
              meta = list(longlat = longlat, kernel_menu = match(kernel_menu, c("exp","mat32","mat52","sph")))
            ), class = "elsakrig")
          }
)

elsa_krige_fit <- function(train, zcol,
                           width = NULL, cutoff = NULL,
                           longlat = FALSE, Lmax = 20,
                           K_entro = 60,
                           ps_lengthscale = TRUE,
                           kernel_menu = c("exp","mat32","mat52","sph")) {
  
  if (inherits(train, "SpatVector")) {
    XY <- as.matrix(terra::geom(train)[, c("x","y")])
    z  <- as.numeric(terra::values(train)[, zcol])
  } else {
    XY <- sp::coordinates(train)
    z  <- as.numeric(train[[zcol]])
  }
  
  # auto lag setup (your logic)
  dd <- neighd(train, 0, 2000000, longlat = longlat)
  r  <- range(unlist(dd), na.rm = TRUE)
  d     <- if (is.null(width)) max(r[1]*2, .Machine$double.eps) else width
  dmax  <- if (is.null(cutoff)) r[2]/3 else cutoff
  lags  <- seq(d, dmax, by = d)
  if (length(lags) > Lmax) lags <- lags[seq_len(Lmax)]
  L <- length(lags); n <- nrow(XY)
  
  # training-point ELSA at each lag (neighbor-averaged like your snippet)
  ELSA_train <- matrix(NA_real_, n, L, dimnames = list(NULL, paste0("lag_", seq_len(L))))
  for (i in seq_len(L)) {
    dn <- dneigh(train, 0, lags[i], longlat = longlat)
    ei <- elsa(train, d = dn, zcol = zcol, categorical = is.factor(z))
    evec <- as.data.frame(ei)[, "ELSA"]
    nb <- dn@neighbours
    for (j in seq_along(nb)) if (!is.null(nb[[j]])) {
      ELSA_train[j, i] <- mean(evec[nb[[j]]], na.rm = TRUE)
    }
  }
  
  
  ls_train <- NULL
  
  if (isTRUE(ps_lengthscale)) {
    S_train <- 1 - pmin(pmax(ELSA_train, 0), 1)
    wlag    <- c(lags[1], diff(lags))
    .fit_r <- function(S_row) {
      ok <- is.finite(S_row); h <- lags[ok]; S <- S_row[ok]
      if (length(h) < 2) return(NA_real_)
      S <- pmin(pmax(S, 1e-6), 1 - 1e-6)
      obj <- function(r) { Sh <- exp(-h / max(r, 1e-6)); sum(wlag[ok] * (S - Sh)^2) }
      optimize(obj, c(max(min(h[h>0])/10, 1e-6), max(h)*10))$minimum
    }
    ls_raw <- apply(S_train, 1, .fit_r)
    if (anyNA(ls_raw)) {
      med <- if (any(is.finite(ls_raw))) stats::median(ls_raw, na.rm = TRUE) else 2 * width
      ls_raw[!is.finite(ls_raw)] <- med
    }
    # NEW: smooth on the log-scale
    ls_smooth <- .smooth_log_ls(XY, ls_raw, bw = 2.5 * d, min_n = 12L)
    
    # store the SMOOTHED ls in the object (PS branch will use it)
    object_ls <- ls_smooth
  } else {
    object_ls <- NULL
  }
  
  structure(list(
    coords = XY, z = z,
    lags = lags, width = d, cutoff = dmax,
    ELSA_train = ELSA_train,
    ls_train   = object_ls,           # <-- smoothed if available
    K_entro = as.integer(K_entro),
    meta = list(longlat = longlat, kernel_menu = match(kernel_menu, c("exp","mat32","mat52","sph")))
  ), class = "elsakrig")
}



# elsa_krige_fit <- function(train, zcol,
#                            width = NULL, cutoff = NULL,
#                            longlat = FALSE, Lmax = 20,
#                            K_entro = 60) {
#   # Coords + values (terra or sp)
#   
#   if (inherits(train, "SpatVector")) {
#     stopifnot(requireNamespace("terra", quietly = TRUE))
#     XY <- terra::geom(train)[, c("x","y")]
#     z  <- terra::values(train)[, zcol]
#   } else {
#     XY <- sp::coordinates(train)
#     z  <- train[[zcol]]
#   }
#   
#   
#   dd <- neighd(train, 0, 20000000, longlat = longlat)
#   r  <- range(unlist(dd), na.rm = TRUE)
#   d  <- if (is.null(width)) max(r[1]*2, .Machine$double.eps) else width
#   dmax <- if (is.null(cutoff)) r[2]/3 else cutoff
#   lags <- seq(d, dmax, by = d); if (length(lags) > Lmax) lags <- lags[seq_len(Lmax)]
#   
#   # 2) ELSA at training points for each lag (neighbor-averaged like your snippet)
#   n <- nrow(XY); L <- length(lags)
#   ELSA_train <- matrix(NA_real_, n, L)
#   colnames(ELSA_train) <- paste0("lag_", seq_len(L))
#   for (i in seq_len(L)) {
#     dn <- dneigh(train, 0, lags[i], longlat = longlat)
#     ei <- elsa(train, d = dn, zcol = zcol, categorical = is.factor(z))
#     evec <- as.data.frame(ei)[, "ELSA"]
#     nb <- dn@neighbours
#     for (j in seq_along(nb)) {
#       if (!is.null(nb[[j]])) ELSA_train[j, i] <- mean(evec[nb[[j]]], na.rm = TRUE)
#     }
#   }
#   
#   structure(list(
#     coords = XY, z = as.numeric(z),
#     lags = as.numeric(lags),
#     width = d, cutoff = dmax,
#     ELSA_train = ELSA_train,
#     K_entro = as.integer(K_entro),
#     meta = list(longlat = longlat)
#   ),
#   class = "elsakrig"
#   )
# }
# #------------
# 
# predict.elsakrig <- function(object, newdata,
#                              K_krige = 40,
#                              sigma2 = stats::var(object$z, na.rm = TRUE),
#                              range_min = object$width,
#                              range_max = 8 * object$width,
#                              nugget_min = 0.05 * sigma2,
#                              nugget_max = 0.5  * sigma2,
#                              q = 1,
#                              kernel_menu = c("exp","mat32","mat52","sph"),
#                              threads = 1,
#                              out_var = TRUE) {
#   
#   stopifnot(inherits(object, "elsakrig"))
#   # Collect prediction coordinates
#   if (inherits(newdata, "SpatRaster")) {
#     stopifnot(requireNamespace("terra", quietly = TRUE))
#     XYp <- terra::xyFromCell(newdata, 1:terra::ncell(newdata))
#   } else if (inherits(newdata, "SpatVector")) {
#     XYp <- as.matrix(terra::geom(newdata)[, c("x","y")])
#   } else {
#     XYp <- as.matrix(newdata) # matrix/data.frame of x,y
#   }
#   
#   # Call the C engine via .Call
#   res <- .Call(
#     "C_elsa_local_krige",
#     PACKAGE = "elsa",
#     as.numeric(object$coords[,1]), as.numeric(object$coords[,2]),
#     as.numeric(object$z),
#     as.numeric(XYp[,1]), as.numeric(XYp[,2]),
#     as.numeric(object$lags),
#     object$ELSA_train,
#     as.integer(object$K_entro),
#     as.integer(K_krige),
#     as.numeric(sigma2),
#     as.numeric(range_min), as.numeric(range_max),
#     as.numeric(nugget_min), as.numeric(nugget_max),
#     as.numeric(q),
#     as.integer(match(kernel_menu, c("exp","mat32","mat52","sph"))),
#     as.integer(threads)
#   )
#   
#   # res is a list: pred, var, range, nugget, kernel_type
#   pred <- res$pred; s2 <- res$var
#   if (inherits(newdata, "SpatRaster")) {
#     r_pred <- newdata[[1]]; terra::values(r_pred) <- pred
#     if (out_var) {
#       r_var <- newdata[[1]]; terra::values(r_var) <- s2
#       return(list(pred = r_pred, var = r_var))
#     } else return(r_pred)
#   }
#   if (out_var) cbind(XYp, pred = pred, var = s2) else cbind(XYp, pred = pred)
# }
# #--------

predict.elsakrig <- function(object, newdata,
                             K_krige = 40,
                             sigma2 = stats::var(object$z, na.rm = TRUE),
                             range_min = object$width,
                             range_max = 4 * object$width,   # tighter default than before
                             nugget_min = 0.05 * sigma2,
                             nugget_max = 0.5  * sigma2,
                             q = 1.5,                         # stronger nugget response
                             kernel = c("auto","exp","mat32","mat52","sph","ps"),
                             model_average = TRUE,            # AIC-weighted combine for auto
                             trim_prop = 0.1,                 # trimmed mean over neighbors
                             range_penalty = 0.01,            # SSE + lambda*(r/rmax)^2
                             taper_theta = 0,                 # 0=off; else range multiplier or absolute (see mode)
                             taper_mode = c("xrange","fixed"),# "xrange": theta = fact * local range
                             threads = 1,
                             out_var = TRUE) {
  
  kernel <- match.arg(kernel)
  taper_mode <- match.arg(taper_mode)
  
  # coords for prediction
  XYp <- if (inherits(newdata, "SpatRaster")) {
    terra::xyFromCell(newdata, 1:terra::ncell(newdata))
  } else if (inherits(newdata, "SpatVector")) {
    as.matrix(terra::geom(newdata)[, c("x","y")])
  } else {
    as.matrix(newdata)
  }
  
  ls_arg <- if (is.null(object$ls_train)) NULL else as.numeric(object$ls_train)
  
  res <- .Call("C_elsa_local_krige_v2", PACKAGE = "elsa",
               as.numeric(object$coords[,1]), as.numeric(object$coords[,2]),
               as.numeric(object$z),
               as.numeric(XYp[,1]), as.numeric(XYp[,2]),
               as.numeric(object$lags),
               object$ELSA_train,
               ls_arg,                     # <— SMOOTHED training ls for PS kernel
               as.integer(object$K_entro), as.integer(K_krige),
               as.numeric(sigma2),
               as.numeric(range_min), as.numeric(range_max),
               as.numeric(nugget_min), as.numeric(nugget_max),
               as.numeric(q),
               as.integer(switch(kernel, auto=-1L, exp=0L, mat32=1L, mat52=2L, sph=3L, ps=99L)),
               as.integer(object$meta$kernel_menu),
               as.integer(model_average),
               as.numeric(trim_prop),
               as.numeric(range_penalty),
               as.numeric(taper_theta),
               as.integer(match(taper_mode, c("xrange","fixed"))),
               as.integer(threads)
  )
  
  pred <- res$pred; s2 <- res$var
  
  if (inherits(newdata, "SpatRaster")) {
    r_pred <- newdata[[1]]; terra::values(r_pred) <- pred
    if (!out_var) return(r_pred)
    r_var <- newdata[[1]]; terra::values(r_var) <- s2
    return(list(pred = r_pred, var = r_var))
  }
  if (out_var) cbind(XYp, pred = pred, var = s2) else cbind(XYp, pred = pred)
}
