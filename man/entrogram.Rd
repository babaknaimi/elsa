\name{entrogram}

\alias{entrogram}
\alias{entrogram,RasterLayer-method}
\alias{entrogram,SpatialPointsDataFrame-method}
\alias{entrogram,SpatialPolygonsDataFrame-method}

\title{Entrogram for Spatial Data}

\description{
Compute sample (empirical) entrogram from spatial data. The function returns a binned entrogram and an entrogram cloud.
}

\usage{
entrogram(x, width, cutoff,...) 
}

\arguments{
  \item{x}{ a spatial object (\code{RasterLayer} or \code{SpatialPoints} or \code{SpatialPolygons}) }
  \item{width}{ the lag size (width of subsequent distance intervals) into which cell pairs are grouped for ELSA estimates. If missing, the cell size (raster resolution) is assigned. }
  \item{cutoff}{ spatial separation distance up to which cell pairs are included in ELSA estimates; as a default, the length of the diagonal of the box spanning the data is divided by three. }
  \item{...}{Additional arguments including \code{zcol} (when \code{x} is Spatial* object, specifies the name of the variable in the dataset; \code{longlat} (when \code{x} is Spatial* object, spacifies whether the dataset has a geographic coordinate system); \code{s} (only when x is a Raster object, it would be useful when the dataset is big, so then by specifying \code{s}, the calculation would be based on a sample with size \code{s} drawn from the dataset, default is \code{NULL} means all cells should be contributed in the calculations)}
  
}
\details{
Entrogram is a variogram-like graph to explore spatial structure in a single variable. An entrogram summarizes the spatial relations in the data, and can be used to understand within what range (distance) the data is spatially autocorrelated.
}

\value{
\item{Entrogram}{an object containing entrogram cloud and the entrogram within each distance interval}
}

\references{

Naimi, B., Hamm, N. A., Groen, T. A., Skidmore, A. K., Toxopeus, A. G., & Alibakhshi, S. (2019). ELSA: Entropy-based local indicator of spatial association. Spatial statistics, 29, 66-88.
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{
\donttest{
file <- system.file('external/dem_example.grd',package='elsa')
r <- raster(file)
plot(r,main='a continuous raster map')

en <- entrogram(r, width=2000)

plot(en)
}


}

\keyword{spatial}
