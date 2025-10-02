\name{eKrige}

\alias{eKrige}
\alias{eKrige,SpatVector-method}


\title{ELSA-based Local Kriging}

\description{
Fit ELSA-driven local kriging model (stores training entrograms and optional length-scale field)
}

\usage{
eKrige(x, zcol,width = NULL, cutoff = NULL,longlat = FALSE, Lmax = 20,
  K_entro = 60,ps_lengthscale = TRUE,kernel_menu = c("exp","mat32","mat52","sph"))
}

\arguments{

  \item{x}{ Training locations. A \code{SpatVector}. Must contain the target variable column specified in \code{zcol}.}
  
  \item{zcol}{ Character. Name of the target variable in \code{x}. (Continuous or categorical; ELSA handles both.)}
  
  \item{width}{ Numeric. Lag width used to build the training entrograms (ELSA vs distance). If NULL, chosen automatically ≈ 2× median nearest-neighbor spacing.}
  
  \item{cutoff}{ Numeric. Maximum lag distance. If NULL, chosen automatically ≈ 1/3 of the maximum neighbor distance.}
  
  \item{longlat}{Logical. If TRUE, neighbor distances are treated as geographic (lon/lat). Default FALSE (projected/planar)}
  
  \item{Lmax}{Integer. Upper bound on number of lags to store (guard for very fine width).}
  
  \item{K_entro}{Integer. Number of nearest training points used later (at prediction time) to form a local entrogram around each prediction cell.}
  
  \item{ps_lengthscale}{Logical. If TRUE, fit a per-site exponential length-scale from each training entrogram, then smooth it spatially (log-scale kernel smoother). Used by the Paciorek–Schervish (PS) kernel.}
  
  \item{kernel_menu}{Character vector. Stationary kernel types considered when \code{predict(..., kernel="auto")}. Allowed values: "exp", "mat32", "mat52", "sph".}
  
  
  
  
}

\details{
The function computes per-lag ELSA at each training point (using the \code{elsa} function + \code{dneigh} function), building an n × L matrix of ELSA values (the training entrograms).

If ps_lengthscale=TRUE, an exponential similarity is fitted to each entrogram row to get a raw length-scale, then smoothed over space to stabilize the PS kernel.

}

\value{
An object of class "elsakrig" with elements:

- coords (n×2), z (n),

- lags, width, cutoff,

- ELSA_train (n×L matrix),

- ls_train (length n; smoothed length-scale if computed, otherwise NULL),

- K_entro (integer),

- meta (list; includes longlat and numeric codes for kernel_menu).

}

\references{

Naimi, B., Hamm, N. A., Groen, T. A., Skidmore, A. K., Toxopeus, A. G., & Alibakhshi, S. (2019). ELSA: Entropy-based local indicator of spatial association. Spatial statistics, 29, 66-88.

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{

# meuse example (sp -> SpatVector)

library(terra); data(meuse, package="sp")
meuse$logZn <- log(meuse$zinc)

v <- vect(meuse, geom=c("x","y"))

fit <- elsa_krige_fit(v, zcol="logZn", K_entro=60, ps_lengthscale=TRUE)

}

\keyword{spatial}