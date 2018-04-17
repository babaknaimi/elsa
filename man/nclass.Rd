\name{nclass}

\alias{nclass}
\alias{nclass,RasterLayer-method}
\alias{nclass,numeric-method}


\title{Best number of classes for categorizing a continuous variable}

\description{
This function explores the best number of classes to categorize (discretize) a continuous variable.
}

\usage{
nclass(x,th)
}

\arguments{
  \item{x}{ a RasterLayer or a numeric vector}
  \item{th}{A threshold (default = 0.005) used to find the best number of classes}
}
\details{
The function uses an approach introduced in Naimi et al. (under review), to find the best number of classes (categories) when a continuous variable is discretizing. The threhold is corresponding to the acceptable level of information loose through discretizing procedure. For the details, see the reference.

}

\value{
An object with the same class as the input \code{x}

}

\references{
Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, Geographical Analysis;

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{

file <- system.file('external/dem_example.grd',package='elsa')
r <- raster(file)
plot(r,main='a continuous raster map')

nclass(r)

nclass(r, th=0.01)

nclass(r, th=0.1)
}

\keyword{spatial}
