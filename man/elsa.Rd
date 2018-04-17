\name{elsa}

\alias{elsa}
\alias{elsa,SpatialPointsDataFrame-method}
\alias{elsa,SpatialPolygonsDataFrame-method}
\alias{elsa,RasterLayer-method}


\title{Entropy-based Local indicators of Spatial Association}

\description{
Calculate ELSA indicator for either a categorical or continuous spatial dataset.
}

\usage{
elsa(x,d,nc,categorical,dif,...) 
}

\arguments{
  \item{x}{ a raster object (\code{RasterLayer} or \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}}
  \item{d}{ numeric local distance, or an object of class neighbours created by dneigh when x is SpatialPoints or SpatialPolygons}
  \item{nc}{ optional, for continuous data it specifies the number of classes through categorizing the variable. If missed, it is automatically calculated (recommended)  }
  \item{categorical}{ logical, specified whether x is a continuous or categorical. If missed the function tries to detect it}
  \item{dif}{the difference between categories, only for categorical }
  \item{...}{ additional arguments including: 
  \code{stat} - specifies which components of ELSA should be returned and can be "E" (default), "L", and "R" (see details);
  \code{cells} - when the elsa is needed for only some cells in the raster data;
  \code{filename} - the name of the raster file when the output should be written
  }
}
\details{
{dif} can be used when categorical values are sorted into hierarchical system (e.g., CORINE land cover). This make it possible to difine different weights of similarity between each pairs of categories when the level of similarity is not the same between different classes in the variable. For example, two categories belong to two forest types are more similar than two categories, one a forest type and the other one an agriculture type. So, it can take this differences into account when the spatial autocorrelation for categorical variables is quantified.

the ELSA statistics has two terms that is called the left term ("L") and the right term ("R") in the reference. It can be specified through the \code{stat} argument if either of these terms should be returned from the function or ELSA ("E"), which is the multiplication of \code{L * R}. All three terms can also be selected.
}

\value{
\item{RasterLayer}{if \code{x} is a \code{RasterLayer}}

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

e <- elsa(r,d=2000,categorical=FALSE)

plot(e)

}
}

\keyword{spatial}