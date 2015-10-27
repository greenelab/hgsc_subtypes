\name{corFinder}
\alias{corFinder}
\title{Calculate pair-wise correlations between samples using the expr() slots of a list of two ExpressionSets.}
\description{This function acts as a wrapper around ComBat (sva package) and
cor(), to calculate pairwise correlations within one or between
two ExpressionSets.}
\usage{corFinder(eset.pair, separator = ":", use.ComBat = TRUE, ...)}
\arguments{
  \item{eset.pair}{a list of ExpressionSets, with two elements.  If the two
elements are identical, return the correlation matrix for pairs
of samples in the first element.  If not identical, return pairs
between the two elements.}
  \item{separator}{Separator between dataset name and sample name.  Dataset names are
added to sample names to keep track of dataset of origin.}
  \item{use.ComBat}{Use the sva::ComBat function for batch correction of the expr()
data between the two datasets.}
  \item{\dots}{Extra arguments passed to the cor() function.}
}

\value{Returns a matrix of sample-wise Pearson Correlations.}

\author{Levi Waldron and Markus Riester}





