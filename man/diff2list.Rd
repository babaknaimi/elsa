\name{dif2list}

\alias{dif2list}
\alias{dif2list,data.frameORmatrix,ANY-method}
\alias{dif2list,data.frameORmatrix-method}



\title{Convert differences in the level of categorical map to a list}

\description{
This function can be used to convert a data.frame including the values specifying the differences (digree of dissimilarities) between the classes (categories) in a categoricl map, to a list.
}

\usage{
dif2list(x,pattern)
}

\arguments{
  \item{x}{ a data.frame containing digrees of dissimilarity between different categories}
  \item{pattern}{a numeric vector, specifies the pattern the data is organised in x; (is needed if the number of classes or subclasses is greater than 9; see examples)}
}
\details{
A user may only need this function when intends to apply elsa for a categorical map with NOT the same level of dissimilarities between different categories. For example, a landuse map may contain several classes of forest, and several classes of agriculture. The level of dissimilarity between a class of forest and a class of agriculture is not the same as between two classes of forest. Elsa can incorporate different levels of dissimilarity by introducing them through a list or a data.frame. Using \code{diff2list} is not necessary to introduce the differences between different categories, as a user can either specify them through a data.frame or a list. However, defining them in a data.frame would be easier specially when the categorical maps are in the form of hierarchical.
For the details about the approach used by ELSA, and some examples, see the reference.

}

\value{
A list

}

\references{
Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, Geographical Analysis;

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{
# imagine we have a categorical map including 4 classes (values 1:4), and the first two classes
# (i.e., 1 and 2) belong to the major class 1 (so can have legends of 11, 12, respectively), and
# the second two classes (i,e, 3 and 4) belong to the major class 2 (so can have legends of 21,
# and 22 respectively). Then we can construct the data.frame as:

d <- data.frame(g=c(1,2,3,4),leg=c(11,12,21,22))

d

# dif2list generates a list including 4 values each corresponding to each value (class in the map
#, i,e, 1:4). Each item then has a numeric vector containing a relative dissimilarity between the
# main class (the name of the item in the list) and the other classes. If one wants to change 
# the relative dissimilarity between two specific classes, then the list can easily be edited and
# used in the elsa function

dif2list(d)


# As you see in the legend, each value contains a sequence of numbers specifying the class, 
#subclass, sub-subclass, .... and so on in a hierarchical manner (for example, 12 means class 1
# and subclass 2). In case if there is more than 9 classes or subclasses (for example, 112 should
# refer to class 1, and subclass 12, not class 1 , subclass 1, and sub-subclass 2), then the
# pattern should be specified as a vector like c(1,2) means that the length of the major class in
# the hierarchy is 1, while the length of the sub class is 2.

d <- data.frame(g=c(1,2,3,4),leg=c(101,102,201,202))

d

dif2list(d,pattern=c(1,2))

}

\keyword{spatial}
