\name{smokingGunFinder}
\alias{smokingGunFinder}
\title{Find doppelgangers based on "smoking gun" phenotypes - those that should be unique to each patient.}
\description{Checks all pairwise combinations of samples for values of the
"smoking" gun phenotypes that are identical.}
\usage{smokingGunFinder(eset.pair, smokingguns, transFun = I, separator = ":")}
\arguments{
  \item{eset.pair}{a list of ExpressionSets, with two elements.  If the two elements
are identical, the function will check for duplicate IDs within
one element. If not identical, it will check for duplicate IDs
between elements.}
  \item{smokingguns}{phenoData column names found in multiple elements of eset.pair
that may contain "smoking guns" such as identifiers that should be
unique to each sample.}
  \item{transFun}{a function to apply to IDs before comparing.  By default apply no transformation.}
  \item{separator}{Separator between dataset name and sample name.  Dataset names are
added to sample names to keep track of dataset of origin.}
}

\value{Returns an adjacency matrix for samples where matches have value
1, non-matches have value zero.  Value for a sample against itself
is NA.}

\author{Levi Waldron and Markus Riester}





