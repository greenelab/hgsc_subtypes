\name{outlierFinder}
\alias{outlierFinder}
\title{##Identifies outliers in a similarity matrix.}
\description{By default uses the
Fisher z-transform for Pearson correlation (atanh), and
identifies outliers as those above the quantile of a skew-t
distribution with mean and standard deviation estimated from the
z-transformed matrix.  The quantile is calculated from the
Bonferroni-corrected cumulative probability of the upper tail.}
\usage{outlierFinder(similarity.mat, bonf.prob = 0.05, transFun = atanh, 
    normal.upper.thresh = NULL, tail = "upper")}
\arguments{
  \item{similarity.mat}{A matrix of similarities - larger values mean more similar.}
  \item{bonf.prob}{Bonferroni-corrected probability.  A raw.prob is calculated by
dividing this by the number of non-missing values in
similarity.mat, and the rejection threshold is qnorm(1-raw.prob,
mean, sd) where mean and sd are estimated from the
transFun-transformed similarity.mat.}
  \item{transFun}{A function applied to the numeric values of similarity.mat, that
should result in normally-distributed values.}
  \item{normal.upper.thresh}{Instead of specifying bonf.prob and transFun, an upper similarity
threshold can be set, and values above this will be considered
likely duplicates.  If specified, this over-rides bonf.prob.}
  \item{tail}{"upper" to look for samples with very high similarity values,
"lower" to look for very low values, or "both" to look for both.}
}

\value{Returns either NULL or a dataframe with three columns: sample1, sample2, and similarity.}

\author{Levi Waldron and Markus Riester}





