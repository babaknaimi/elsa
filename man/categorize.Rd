\name{categorize}

\alias{categorize}
\alias{categorize,RasterLayer-method}
\alias{categorize,numeric-method}

\title{Categorizing continious data}

\description{
A function to categorize (discretize) numerical data in the form of a raster layer, or a vector.
}

\usage{
categorize(x,nc,...)
}

\arguments{
  \item{x}{ A \code{RasterLayer} or a \code{numerical vector}}
  \item{nc}{ Number of classes, if missing, it will be automatically detected}
  \item{...}{ Additional parameters for \code{\link[raster]{writeRaster}} function (\code{filename} should be specified as additional argument)}
  
}
\details{
If nc is not specified, the function calls the function \code{nclass} to find the best number of classes. For the details of how nclass works, see the reference.
}

\value{
\item{RasterLayer}{if x is a \code{RasterLayer}}

\item{numeric vector}{if x is a numeric \code{vector}}
}

\references{
Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, Geographical Analysis;
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
