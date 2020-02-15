\name{categorize}

\alias{categorize}
\alias{categorize,RasterLayer-method}
\alias{categorize,RasterStackBrick-method}
\alias{categorize,numeric-method}
\alias{categorize,list-method}

\title{Categorizing continious data}

\description{
A function to categorize (discretize or binning) numerical data in the form of a raster layer, or a vector.
}

\usage{
categorize(x,nc,probs,...)
}

\arguments{
  \item{x}{ A \code{RasterLayer} or a \code{numerical vector}}
  \item{nc}{ Number of classes, if missing, it will be automatically detected}
  \item{probs}{ a two-item numeric vector of the lower and upper probabilities within the range of [0,1], e.g., c(0.025, 0.975) (default: NULL) to consider a quantile range of values in categorization (to avoid the effects of outliers); if \code{NULL} or \code{FALSE}, then the whole range is used.}
  \item{...}{ Additional parameters for \code{\link[raster]{writeRaster}} function (\code{filename} should be specified as additional argument)}
  
}
\details{
If nc is not specified, the function calls the function \code{nclass} to find the best number of classes. For the details of how nclass works, see the reference.

When outliers do exist in \code{x}, it is likely that they affect the quality of categorization. To avoid their affects, a quantile with the specified probabilities can be used. default is recommended: c(0.025,0.975), but to ignore it, you can turn the probs argument to \code{NULL}.
}

\value{
\item{RasterLayer}{if x is a \code{RasterLayer}}

\item{numeric vector}{if x is a numeric \code{vector}}
}

\references{

Naimi, B., Hamm, N. A., Groen, T. A., Skidmore, A. K., Toxopeus, A. G., & Alibakhshi, S. (2019). ELSA: Entropy-based local indicator of spatial association. Spatial statistics, 29, 66-88.
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{

categorize(1:10,3) # categorizeing a numeric vector into 3 categories

file <- system.file('external/dem_example.grd',package='elsa')
r <- raster(file)
plot(r,main='a continuous raster map')
rc <- categorize(r,nc=4)
plot(rc, main='categorized map')

}

\keyword{spatial}
