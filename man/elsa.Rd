\name{elsa}

\alias{elsa}
\alias{elsa,SpatialPointsDataFrame-method}
\alias{elsa,SpatialPolygonsDataFrame-method}
\alias{elsa,RasterLayer-method}


\title{Entropy-based Local indicator of Spatial Association}

\description{
Calculate ELSA statistic for a categorical or continuous spatial dataset.
}

\usage{
elsa(x,d,nc,categorical,dif,classes,stat,...) 
}

\arguments{
  \item{x}{ a raster object (\code{RasterLayer} or \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}}
  \item{d}{ numeric local distance, or an object of class neighbours created by dneigh when x is SpatialPoints or SpatialPolygons}
  \item{nc}{ optional, for continuous data it specifies the number of classes through categorizing the variable. If missing, it is automatically calculated (recommended)}
  \item{categorical}{ logical, specified whether x is a continuous or categorical. If missed the function tries to detect it}
  \item{dif}{the difference between categories, only for categorical }
  \item{classes}{Optional, only when \code{x} is categorical is a character vector contains classes; would be useful when the dataset is part of a bigger dataset or when it does not contain all the categories, then by specifying the full set of categories, they will be taken into account to calculate ELSA, and therefore, it would be comparable with the other dataset with the same list of classes (these classes may alternatively introduce by \code{dif}, as the classes considered to specify dissimilarities in \code{dif} list, would be used as \code{classes})}
  \item{stat}{specifies which statistic should be calculated by the function; it can be "elsa" (default), or either of the two components of the statistic, "Ea", or "Ec"; ELSA is the product of Ea and Ec. (it is possible to select more than one statistic); this argument is ignored if \code{x} is \code{Spatial*} object as all the three statistics are returned (see details)}
  
  \item{...}{ additional arguments including: 
  
  \code{cells} - a numeric vector to specify for which raster cells the ELSA statistic should be calculated; it works when \code{x} is \code{RasterLayer}, and if it is specified, ELSA is calculated only for the specified cells;
  
  \code{filename} - only if \code{x} is \code{RasterLayer}, specifies the name of the raster file to be written in the working directory.
  
  \code{zcol} - only if \code{x} is \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}, specifies the name of the column (variable) in the attribute table for which the elsa statistics are calcualted (i.e., ELSA, Ea, Ec).
  
 \code{drop} - logical; only if \code{x} is \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}, specifies whether the output should be a data.frame or a \code{Spatial*} object.
  
  \code{method} - only if \code{SpatialPolygonsDataFrame}, specifies the method for identifying the neighbourhood polygons; default: 'centroid' (see \code{dneigh}).
  
  
  }
}
\details{
{dif} can be used when categorical values are sorted into hierarchical system (e.g., CORINE land cover). This make it possible to difine different weights of similarity between each pairs of categories when the level of similarity is not the same between different classes in the variable. For example, two categories belong to two forest types are more similar than two categories, one a forest type and the other one an agriculture type. So, it can take this differences into account when the spatial autocorrelation for categorical variables is quantified.

the ELSA statistics has two terms, "Ea" and "Ec", in the reference. It can be specified in the \code{stat} argument if either of these terms should be returned from the function or ELSA ("E"), which is the product of these two terms, \code{Ea * Ec}. All three terms can also be selected.
}

\value{
\item{Raster*}{if \code{x} is a \code{RasterLayer}}
\item{Spatial* or data.frame}{if \code{x} is a \code{Spatial*}}

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

plot(r, main='a continuous raster map')

e <- elsa(r,d=2000,categorical=FALSE)

plot(e)
}

\keyword{spatial}