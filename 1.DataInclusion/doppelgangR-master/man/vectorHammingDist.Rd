\name{vectorHammingDist}
\alias{vectorHammingDist}
\title{Calculate Hamming Distance between two vectors, using pairwise complete observations.}
\description{Simple function to count the fraction of different elements (in
the same position) between two vectors of the same length, after
removing elements from both vectors corresponding to positions
that are NA in either vector.}
\usage{vectorHammingDist(x, y, k, l)}
\arguments{
  \item{x}{a matrix }
  \item{y}{a matrix with the same number of columns as x}
  \item{k}{row in x to test for differences}
  \item{l}{row in y to test for differences}
}

\value{Returns a numeric value, the Hamming Distance (the number of
non-equal values between x and y).}

\author{Levi Waldron and Markus Riester}





