\name{correlogram}

\alias{correlogram}
\alias{correlogram,RasterLayer-method}
\alias{correlogram,Spatial-method}


\title{Correlogram from Spatial Data}

\description{
Compute Correlogram from spatial data. The function returns a binned correlogram by calculating Moran's I (spatial autocorrelation) in different range of distances.
}

\usage{
correlogram(x, width, cutoff,...) 
}

\arguments{
  \item{x}{ a spatial object (\code{RasterLayer} or \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}) }
  \item{width}{ the lag size (width of subsequent distance intervals) into which cell pairs are grouped for semivariance estimates. If missing, the cell size (raster resolution) is assigned. }
  \item{cutoff}{ spatial separation distance up to which cell pairs are included in semivariance estimates; as a default, the length of the diagonal of the box spanning the data is divided by three. }
  \item{...}{Additional arguments including \code{zcol} (when \code{x} is Spatial* object, specifies the name of the variable in the dataset; \code{longlat} (when \code{x} is Spatial* object, spacifies whether the dataset has a geographic coordinate system); \code{s} (only when x is a Raster object, it would be useful when the dataset is big, so then by specifying \code{s}, the calculation would be based on a sample with size \code{s} drawn from the dataset, default is \code{NULL} means all cells should be contributed in the calculations)}
  
}
\details{
Correlogram is a graph to explore spatial structure in a single variable. A correlogram summarizes the spatial relations in the data, and can be used to understand within what range (distance) the data is spatially autocorrelated.
}

\value{
\item{Correlogram}{an object containing Moran's I values within each distance interval}
}

\references{
Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, Geographical Analysis;
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{
\dontrun{
file <- system.file('external/dem_example.grd',package='elsa')
r <- raster(file)
plot(r,main='a continuous raster map')

co <- correlogram(r, width=2000,cutoff=30000)

plot(co)
}


}

\keyword{spatial}
