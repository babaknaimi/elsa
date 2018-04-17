\name{lisa}

\alias{lisa}
\alias{lisa,Spatial-method}
\alias{lisa,RasterLayer-method}


\title{Local indicators of Spatial Associations}

\description{
Calculate local indicators of spatial association (LISA) for a continuous (numeric) variable at each location in a Raster layer or a SpatialPointsDataFrame or a SpatialPolygonsDataFrame.
}

\usage{
lisa(x,d1,d2,statistic,...) 
}

\arguments{
  \item{x}{ a raster object (\code{RasterLayer} or \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}}
  \item{d1}{ numeric lower bound of local distance (default=0), or an object of class neighbours created by dneigh when x is SpatialPoints or SpatialPolygons}
  \item{d2}{ numeric upper bound of local distance, not needed if d1 is a neighbours object,}
  \item{statistic}{ a character string specifying the LISA statistic that should be calculated. This can be one of "I" (or "localmoran" or "moran"), "c" (or "localgeary" or "geary"), "G" (or "localG"), "G*" (or "localG*")}
  \item{...}{ additional arguments including \code{filename} (only when x is Raster, specifies the name of the raster file when the output should be written; additional arguments for \code{\link[raster]{writeRaster}} function can also be specified);
  \code{mi} (only when x is Raster and \code{statistic='I'}, specifies whether raw Local Moran's I statistic (Ii) should be returned, or standardized value (Z.Ii). e.g., mi="I", mi='Z' (default)); 
  \code{zcol} (only when x is a Spatial* object specifies the name of the variable column in the data); 
  \code{longlat} (logical, only when x is a Spatial* object specifies whether the coordinate system is geographic); 
  \code{drop} (logical, only when x is a Spatial* object, if TRUE, the original data structure (Spatial* object) is returned, otherwise a numeric vector is returned) }
}
\details{
This function can calculate different LISA statistics at each location in the input dataset. The statistics, implemented in this function, include local Moran's I ("I"), local Geary's c ("c"), local G and G* ("G" and "G*"). This function returns standardized value (Z) for Moran, G and G*.
}

\value{
\item{RasterLayer}{if \code{x} is a \code{RasterLayer}}

\item{Spatial*}{if \code{x} is a \code{Spatial*} and \code{drop=FALSE}}

\item{numeric vector}{if \code{x} is a \code{Spatial*} and \code{drop=TRUE}}
}

\references{

Anselin, L. 1995. Local indicators of spatial association, Geographical Analysis, 27, 93--115;

Getis, A. and Ord, J. K. 1996 Local spatial statistics: an overview. In P. Longley and M. Batty (eds) \emph{Spatial analysis: modelling in a GIS environment} (Cambridge: Geoinformation International), 261--277.

Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, \emph{Geographical Analysis}

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{
\dontrun{

file <- system.file('external/dem_example.grd',package='elsa')
r <- raster(file)

plot(r,main='a continuous raster map')

mo <- lisa(r,d2=2000,statistic='i') # local moran's I (Z.Ii value)

plot(mo, main="local Moran's I (Z.Ii)")

mo <- lisa(r,d2=2000,statistic='i',mi='I') # local moran's I (Ii value  (non-standardized))

plot(mo, main="local Moran's I (Ii))")

gc <- lisa(r,d2=2000,statistic='c') # local Geary's c

plot(gc, main="local Geary's c")

g <- lisa(r,d2=2000,statistic='g') # local G

plot(g, main="local G")

}
}

\keyword{spatial}