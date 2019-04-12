\name{dif2list}

\alias{dif2list}
\alias{dif2list,data.frameORmatrix,ANY-method}
\alias{dif2list,data.frameORmatrix-method}



\title{Convert differences in the level of categorical map to a list}

\description{
This function converts a data.frame including the values specifying the differences (contrast or degrees of dissimilarities) between the classes (categories) in a categoricl map, to a list.
}

\usage{
dif2list(x,pattern,fact=1)
}

\arguments{
  \item{x}{ a data.frame containing degrees of dissimilarities between different categories}
  \item{pattern}{a numeric vector, specifies the pattern the data is organised in x; (only when the number of classes or subclasses is greater than 9; see the examples)}
  \item{fact}{ a numeric value (default = 1), specifies the factor that multiplied to the estimated dissimilarities}
}
\details{
When ELSA is calculated for a categorical map, by default it is assumed that the level of dissimilarities (or level of contrast) between different classes are the same. For example, if a categorical map has four classes of "A","B","C", and "D", dissimilarity or contrast between "A" & "B" is the same as between "A" & "C", or between "C" & "D", etc. Sometimes, it is not a valid assumption as some classes might be more similar than the others. For example, a landuse map may contain several classes of forest, and several classes of agriculture. The level of dissimilarity between a class of forest and a class of agriculture is not the same as between two different types of forest.

ELSA is flexible enough to incorporate different levels of dissimilarity. To do that, a list can be specified in which you specify the contrast or dissimilarity between each class (name of an item in the list) to the other classes (named vector with items corresponsing to the other classes). For instance, list(A=c(A=0,B=1,C=2,D=2), B=c(A=1,B=0,C=2,D=2),...) simply specifies the dissimilarities between class A and class B (the first two items in the list) to the four classes of A, B, C, and D (the vectors assigned to the first and second item in the list). 

Alternatively, a simple coding approach can be used when the classes can be organised in a hierarchical way (see Naimi et al. (2019) for more details), and if that is the case, a code can be assigned to each class within a data.frame. Then, \code{dif2list} function can be used to convert the data.frame to a list with the structure required by the function (like the above example). Using \code{diff2list} is not necessary to introduce the differences between different categories, as a user can either specify them directly in a list. However, defining them in a data.frame would be easier specially when the classes are in related in a hierarchical form. 

In the data.frame required by \code{dif2list}, a code is assigned to each class. The code is specified according to a hierarchical structure and the number of levels in it. For instance, a land cover may have four classes including forest broadleaves, forest needle-leaves, cropland rainfed, and cropland irrigated. These four classes can be organised hierarchically at two levels, the first level would be forst and cropland, and the second level would be the name of the classes. The codes assigned to these four classes can be 11, 12, 21, and 22, respectively. The dissimilarities between, for example, classes of 11 & 12 would be 1 and between 11 & 21 would be 2. The estimated dissimilarities can be adjusted by changing the \code{fact} that is multiplied to the dissimilarity values.

}

\value{
A list

}

\references{

Naimi, B., Hamm, N. A., Groen, T. A., Skidmore, A. K., Toxopeus, A. G., & Alibakhshi, S. (2019). ELSA: Entropy-based local indicator of spatial association. Spatial statistics, 29, 66-88.

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

dif2list(d,pattern=c(1,2), fact=2) # dissimilarities are multiplied by 2 (fact=2).

}

\keyword{spatial}
