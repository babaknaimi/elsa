\name{dneigh}

\alias{dneigh}
\alias{dneigh,SpatialPoints-method}
\alias{dneigh,SpatialPolygons-method}
\alias{dneigh,data.frameORmatrix-method}


\title{Construct neighbours list}

\description{
This function identifies the neighbours features (points or polygons) given the specified distance and builds a list of neighbours.
}

\usage{
dneigh(x,d1,d2,longlat,method,...)
}

\arguments{
  \item{x}{ a SpatialPoints, or SpatialPolygons or a matrix (or data.frame) of point coordinates or a SpatialPoints object}
  \item{d1}{ lower local distance bound }
  \item{d2}{ upper local distance bound }
  \item{longlat}{ TRUE if point coordinates are longitude-latitude}
  \item{method}{if x is SpatialPolygons, specifies the method to identify the neighbour polygons; see details}
  \item{...}{ additional arguments; see details}
}
\details{
The function is mostly based on \code{dnearneigh} (for points), and \code{poly2nb} (for polygons), implemented in the spdep package by Roger Bivand.

When \code{x} is SpatialPolygons, there is two methods (can be specified through \code{method}) to identify the neighbour polygons. The default method ('bound') seeks the polygons that has one or more points in their boundaries within the specified distance (\code{d}), while the method 'centroid' considers any polygon with a centriod within the given distance.

One additional argument is queen (default is TRUE), can beused only when x is SpatialPolygons, and method='bound', if TRUE, a single shared boundary point meets the contiguity condition, if FALSE, more than one shared point is required

}

\value{
An object of class \code{neighbours}

}

\references{
Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, Geographical Analysis;

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{
#
}

\keyword{spatial}
