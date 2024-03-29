\name{moran}

\alias{moran}
\alias{geary}
\alias{moran,Spatial-method}
\alias{moran,RasterLayer-method}
\alias{moran,SpatRaster-method}
\alias{moran,SpatVector-method}
\alias{geary,Spatial-method}
\alias{geary,RasterLayer-method}
\alias{geary,SpatRaster-method}
\alias{geary,SpatVector-method}



\title{Global Spatial Autocorrelation Statistics}

\description{
Functions to calculate Moran's I and Geary's c statistics.
}

\usage{
moran(x,d1,d2,...) 
geary(x,d1,d2,...) 
}

\arguments{
  \item{x}{ a raster object (\code{RasterLayer} or \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}}
  \item{d1}{ lower bound local distance, or an object of class neighbours created by dneigh when x is SpatialPoints or SpatialPolygons}
  \item{d2}{ upper bound local distance}
  \item{...}{ additional arguments including \code{zcol} (when \code{x} is Spatial* object, specifies the name of the variable in the dataset; \code{longlat} (when \code{x} is Spatial* object, spacifies whether the dataset has a geographic coordinate system}
}
\details{
moran and geary are two functions to measure global spatial autocorrelation within the range of distance specified through d1 and d2. It returns a single numeric value than can show the degree of spatial autocorrelation in the whole dataset.
}

\value{
A numeric value.

}

\references{

Naimi, B., Hamm, N. A., Groen, T. A., Skidmore, A. K., Toxopeus, A. G., & Alibakhshi, S. (2019). ELSA: Entropy-based local indicator of spatial association. Spatial statistics, 29, 66-88.

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{
file <- system.file('external/dem_example.grd',package='elsa')
r <- raster(file)

moran(r, d1=0, d2=2000)

geary(r, d1=0, d2=2000)
}

\keyword{spatial}