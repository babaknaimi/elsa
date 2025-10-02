\name{predict}

\alias{predict}
\alias{predict,SpatVector-method}


\title{Predict/Interpolate with ELSA-driven local kriging (stationary kernels or Paciorek–Schervish)}

\description{
Based on the fitted ELSA-driven local kriging model (using \code{eKrige}), this function predict (interpolate) values for new coordinates (provided by a SpatRaster),
}

\usage{
predict(
  object, newdata,
  K_krige = 40,
  sigma2 = stats::var(object$z, na.rm = TRUE),
  range_min = object$width,
  range_max = 4 * object$width,
  nugget_min = 0.05 * sigma2,
  nugget_max = 0.5  * sigma2,
  q = 1.5,
  kernel = c("auto","exp","mat32","mat52","sph","ps"),
  model_average = TRUE,
  trim_prop = 0.10,
  range_penalty = 0.01,
  taper_theta = 0,
  taper_mode = c("xrange","fixed"),
  threads = 1,
  out_var = TRUE
)
}

\arguments{

  \item{object}{ An "elsakrig" model from \code{eKrige()}.}
  
  \item{newdata}{ Where to predict. A SpatRaster (returns rasters), a SpatVector of points, or a numeric matrix/data.frame with columns (x, y)}
  
  \item{K_krige}{ Integer. Number of nearest training points used in each local kriging system (moving neighborhood).}
  
  \item{sigma2}{ Numeric. Sill (process variance). Default: variance of training z.}
  
  \item{range_min}{Numerics. Lower bound for the local range when fitting the kernel to the local entrogram. Tighter bounds → less over-smoothing.}
  
  \item{range_max}{Numerics. upper bound for the local range when fitting the kernel to the local entrogram. Tighter bounds → less over-smoothing.}
  
  \item{nugget_min}{Numerics. Lower bound for the local nugget. The nugget at a pixel is mapped from the innermost-lag ELSA via nugget_min + (nugget_max - nugget_min) * ELSA^q.}
  
  \item{nugget_max}{Numerics. Upper Bound for the local nugget. The nugget at a pixel is mapped from the innermost-lag ELSA via nugget_min + (nugget_max - nugget_min) * ELSA^q.}
  
  \item{q}{Numeric (≥1). Exponent controlling how strongly high-ELSA (patchy) areas inflate the nugget. Larger q → more edge preservation.}
  
  \item{kernel}{"auto" (default), a fixed stationary type ("exp", "mat32", "mat52", "sph"), or "ps" for the Paciorek–Schervish nonstationary Gaussian kernel.}
  
  \item{model_average}{Logical. If TRUE and kernel="auto", combine candidate kernels using AIC weights; if FALSE, use the single best by AIC.}
  
  \item{trim_prop}{Proportion (0–0.4). Trim used when averaging ELSA across neighbor curves at each lag (robust entrogram; higher trims reduce over-smoothing).}
  
  \item{range_penalty}{Numeric ≥0. Penalty term added to the WLS objective during range fitting (SSE + lambda*(r/r_max)^2) to prevent long-range creep.}
  
  \item{taper_theta}{Numeric. Compact support taper radius. 0 disables tapering. For stationary kernels, a small taper (≈ 2–3× local range) sharpens maps. For kernel="ps", tapering is disabled internally (nonstationary + taper may disconnect neighborhoods).}
  
  \item{taper_mode}{"xrange" or "fixed". If "xrange", the effective taper is theta = taper_theta × local_range; if "fixed", theta is taken as an absolute distance.}
  
  \item{threads}{Integer. Number of OpenMP threads (if compiled with OpenMP). Use >1 for multi-core predictions.}
  
  \item{out_var}{Logical. If TRUE, also return the kriging variance.}
  
  
  
  
}

\details{
Using this function, around each prediction location, the local entrogram is built by averaging nearby training-point entrograms (ELSA vs distance) across K_entro neighbors, using a trimmed mean per lag.

Then, it is converted to a similarity curve S(h) = 1 − ELSA(h) and fit a valid kernel locally (range only if type fixed). With kernel="auto", several types are tried and combined or chosen by AIC.

It also maps the inner-lag ELSA to a local nugget using exponent q.

It then solves an ordinary kriging system in a K_krige moving window to get the prediction and variance.

For kernel="ps", a smoothed length-scale field drives the globally PD Paciorek–Schervish kernel; local ratio capping avoids “islands”.

The kernel parameter defines how the local kernel should be fit:

- auto: fit range for each candidate in kernel_menu (from fit) and combine predictions (see model_average).

- fixed: fit only the local range; the type is held constant.

- ps: uses the smoothed length-scale field (ls_train) and a stabilized PS construction (smoothed ℓ0, ratio capping).



}

\value{
If newdata is a SpatRaster:

list(pred = SpatRaster, var = SpatRaster) when out_var=TRUE; otherwise a single SpatRaster of predictions.
}

\references{

Naimi, B., Hamm, N. A., Groen, T. A., Skidmore, A. K., Toxopeus, A. G., & Alibakhshi, S. (2019). ELSA: Entropy-based local indicator of spatial association. Spatial statistics, 29, 66-88.

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{

library(terra); data(meuse, package="sp"); data(meuse.grid, package="sp")
meuse$logZn <- log(meuse$zinc)

v <- vect(meuse, geom=c("x","y"))

# fit the elsa_based kriging (generates object/model)
fit <- eKrige(v, zcol="logZn", K_entro=45, ps_lengthscale=TRUE)

# Stationary (auto types, AIC model averaging), with mild taper tied to local range

g <- rast(meuse.grid)[[1]] # a base raster providing locations for interpolation

pr1 <- predict(fit, g, kernel="auto", model_average=TRUE,
               K_krige=35, q=1.5, taper_theta=2.5, taper_mode="xrange")

# Fixed exponential (often best RMSE)
pr2 <- predict(fit, g, kernel="exp", K_krige=35, q=1.5, range_max=3*fit$width)

# Nonstationary PS kernel (uses smoothed length-scales from fit)
pr3 <- predict(fit, g, kernel="ps", K_krige=35)

pr1 <- mask(rast(pr1),g) # let's mask it based on g
pr2 <- mask(rast(pr2),g) 
pr3 <- mask(rast(pr3),g) 

plot(pr1)
# plot(pr2)
# plot(pr3)

}

\keyword{spatial}