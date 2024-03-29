\name{elsa.test}

\alias{elsa.test}
\alias{elsa.test,RasterLayer-method}
\alias{elsa.test,SpatRaster-method}
\alias{elsa.test,Spatial-method}



\title{Elsa test for local spatial autocorrelation}

\description{
This function uses a non-parametric approach to test whether local spatial autocorrelation (characterised by ELSA) is significant. It generates a p-value at each spatial location (a raster cell or spatial point/polygon) that can be used to infer the significancy of local spatial autocorrelation.
}

\usage{
elsa.test(x, d, n, method, null, nc, categorical, dif,classes,...)
}

\arguments{
  \item{x}{ A Raster or Spatial* dataset}
  \item{d}{the local distance, or an object of class neighbours created by dneigh function}
  \item{n}{number of simulation, default is 999 for small datasets, and 99 for large datasets}
  \item{null}{ Optional, a null distribution of data (a Raster if x is Raster or a numerical vector if x is either Raster or Spatial dataset ); if not provided, a null distribution is generated by the function}
  \item{method}{resampling method for nonparametric simulation, can be either 'boot' (bootstraping; default) or 'perm' (permutation) }
  \item{nc}{number of classes (only if x is a continuous variable); if not specified, it is estimated using nclass function}
  \item{categorical}{logical, specifies whether x is a categorical; if not specified, it is guessed by the function}
  \item{dif}{the level of dissimilarities between different categories (only if x is a categorical variable); see \code{\link[elsa]{dif2list}} for more details}
  \item{classes}{Optional, only when \code{x} is categorical is a character vector contains classes; would be useful when the dataset is part of a bigger dataset or when it does not contain all the categories, then by specifying the full set of categories, they will be taken into account to calculate ELSA, and therefore, it would be comparable with the other dataset with the same list of classes (these classes may alternatively introduce by \code{dif}, as the classes considered to specify dissimilarities in \code{dif} list, would be used as \code{classes})}
  \item{...}{Aditional arguments passed to writeRaster function (applied only when x is Raster)}
}
\details{
This function test how significant the local spatial autocorrelation is at each location, so it generates a p-value at each location through a Monte Carlo simulation and a non-parametric approach. See the reference for the details about the method.

If \code{null} distribution is not provided, the function generates a null distribution by randomly shuffling the values in the dataset. 

}

\value{
An object same as the input (\code{x})

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

et <- elsa.test(r,d=2000,n=99, categorical=FALSE)

plot(et)

}

}

\keyword{spatial}
