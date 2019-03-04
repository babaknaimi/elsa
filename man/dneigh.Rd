\name{dneigh}

\alias{dneigh}
\alias{neighd}
\alias{dneigh,SpatialPoints-method}
\alias{dneigh,SpatialPolygons-method}
\alias{dneigh,data.frameORmatrix-method}
\alias{neighd,SpatialPoints-method}
\alias{neighd,SpatialPolygons-method}
\alias{neighd,data.frameORmatrix-method}



\title{Construct neighbours list}

\description{
This function identifies the neighbours features (points or polygons) given the specified distance (in kilometer for geographic coordinates, i.e., if longlat=TRUE; and in the map unit for projected datasets, i.e., if longlat = FALSE) and builds a list of neighbours.

The neighd function returns a list including distance of each feature to neighbourhood locations.
}

\usage{
dneigh(x,d1,d2,longlat,method,...)

neighd(x,d1,d2,longlat,...)
}

\arguments{
  \item{x}{ a SpatialPoints, or SpatialPolygons or a matrix (or data.frame) of point coordinates or a SpatialPoints object}
  \item{d1}{ lower local distance bound (if longlat = TRUE, in kilometer; otherwise in the spatial unit of the dataset, e.g., meter) }
  \item{d2}{ upper local distance bound (if longlat = TRUE, in kilometer; otherwise in the spatial unit of the dataset, e.g., meter)}
  \item{longlat}{ TRUE if point coordinates are longitude-latitude}
  \item{method}{if x is SpatialPolygons, specifies the method to identify the neighbour polygons; see details}
  \item{...}{ additional arguments; see details}
}
\details{
The function is mostly based on \code{dnearneigh} (for points), and \code{poly2nb} (for polygons), implemented in the spdep package by Roger Bivand.

When \code{x} is SpatialPolygons, there is two methods (can be specified through \code{method}) to identify the neighbour polygons. The default method ('bound') seeks the polygons that has one or more points in their boundaries within the specified distance (\code{d}), while the method 'centroid' considers any polygon with a centriod within the given distance.

One additional argument is queen (default is TRUE), can beused only when x is SpatialPolygons, and method='bound', if TRUE, a single shared boundary point meets the contiguity condition, if FALSE, more than one shared point is required.

\code{neighd} for SpatialPolygons returns distances of each polygon to centroids of neighbor polygons.

}

\value{
An object of class \code{neighbours}

}

\references{

Naimi, B., Hamm, N. A., Groen, T. A., Skidmore, A. K., Toxopeus, A. G., & Alibakhshi, S. (2019). ELSA: Entropy-based local indicator of spatial association. Spatial statistics, 29, 66-88.

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{
#
}

\keyword{spatial}
