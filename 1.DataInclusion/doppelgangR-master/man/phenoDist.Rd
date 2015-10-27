\name{phenoDist}
\alias{phenoDist}
\title{Calculate distance between two vectors, rows of one matrix/dataframe, or rows of two matrices/dataframes.}
\description{This function does some simple looping to allow x and y to be
various combinations of vectors and matrices/dataframes.}
\usage{phenoDist(x, y = NULL, bins = 10, vectorDistFun = vectorWeightedDist, 
    ...)}
\arguments{
  \item{x}{A vector, matrix or dataframe}
  \item{y}{NULL, a vector, matrix, or dataframe.  If x is a vector, y must also be specified.}
  \item{bins}{discretize continuous fields in the specified number of bins}
  \item{vectorDistFun}{A function of two vectors that returns the distance between those vectors.}
  \item{\dots}{Extra arguments passed on to vectorDistFun}
}

\value{a matrix of distances between pairs of rows of x (if y is
unspecified), or between all pairs of rows between x and y (if
both are provided).}

\author{Levi Waldron and Markus Riester}





