\name{evalContour}

\alias{evalContour}

\title{Evaluating Convex Polytopes}

\description{
Given the system of inequalities \code{AAMat\%*\%ZVec <= BBVec} 
describing a convex polytope/contour with an interior point 
\code{IPVec} in the Euclidean space of dimension two to six, this 
function identifies all nonredundant constraints and computes 
some characteristics of the resulting convex polytope such as 
its vertices, facets, volume and surface area.
}

\usage{
evalContour(AAMat, BBVec = NULL, IPVec = NULL)
}

\arguments{
\item{AAMat}{the constraints matrix from the system 
             of inequalities defining the convex polytope.  
             It should be a numeric matrix with two to six 
             columns.}

\item{BBVec}{the right-hand side from the system of 
             inequalities defining the convex polytope.  
             It should be a numeric column vector 
             of the same length as the first column 
             of \code{AAMat}.}

\item{IPVec}{an interior point of the investigated convex 
             polytope.  This argument can be omitted or set 
             equal to a numeric column vector of the same 
             length as the first row of \code{AAMat}.  
             If \code{IPVec} is \code{NULL} or if given 
             \code{IPVec} does not lie well 
             inside the convex contour, then a 
             well-positioned interior point is searched for 
             internally, which may slow down the computation 
             and make it less reliable.}
}

\details{
This function is included to be used for evaluating 
(regression) quantile contours or their cuts.  

In fact, the function analyzes not the polytope 
itself, but its regularized intersection with 
the zero-centered hypercube of the edge length 
2 000 that is employed as an artificial bounding 
box to avoid the problems with unbounded contours.  
The regularization consists of rounding the vertices 
(i.e., all of their coordinates) of such 
an intersection to the seventh decimal digit and of 
considering only the polytope determined by all 
the distinct rounded vertices for the final analysis.  
}

\value{
evalContour returns a list with the following components 
describing the resulting convex polytope:  
\item{Status}{
0 - OK.

2 - the contour seems virtually empty.

3 - the search for a well-positioned interior 
    point \code{IPVec} failed.

4 - the number of input parameters is too low.

5 - \code{AAMat} is not a numeric matrix with two to six 
    columns.

6 - \code{BBVec} is not a numeric column vector of 
    the right length.

7 - \code{IPVec} is not a numeric column vector of 
    the right length.
}

\item{TVVMat}{the matrix with clearly distinct contour 
              vertices (in rows).}
  
\item{TKKMat}{the matrix with clearly distinct elementary 
              facets (in rows).  
              Each row contains the indices of the rows 
              of \code{TVVMat} where the facet vertices are stored.  
              Each facet has the same number of vertices 
              equal to the number of columns of \code{AAMat}.  
              See also \code{help(convhulln)} for 
              the meaning of \code{TKKMat}.}

\item{NumF}{the number of clearly distinct contour facets.}

\item{NumV}{the number of clearly distinct contour vertices.}

\item{Vol}{the volume of the contour (the area in 2D).}

\item{Area}{the surface area of the contour (the 
              circumference in 2D and the surface in 3D).}

}

\examples{
##- a simple example using a tilted zero-centered square 
AAMat <- rbind(c(-1,-1), c(1,-1), c(1,1), c(-1,1))
BBVec <- c(1, 1, 1, 1)
IPVec <- c(0, 0)
CST <- evalContour(AAMat, BBVec, IPVec)
print(CST)

##- computing and evaluating the 0.15-quantile contour of 199 
##random points uniformly distributed in the unit square 
##centered at zero
Tau  <- 0.15
YMat <- matrix(runif(2*199, -0.5, 0.5), 199, 2)
C <- compContourM1u(Tau, YMat)
CST <- evalContour(-C$CharST$HypMat[,1:2], -C$CharST$HypMat[,3])
print(CST)

##See also the examples ExampleA to ExampleE for some 
##more elaborate ways of computing, evaluating and 
##plotting the (regression) quantile contours.
}
