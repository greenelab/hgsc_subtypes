\name{doppelgangR}
\alias{doppelgangR}
\title{doppelgangR}
\description{Identify samples with suspiciously high correlations and phenotype similarities}
\usage{doppelgangR(esets, separator = ":", corFinder.args = list(separator = separator,
    use.ComBat = TRUE, method = "pearson"), phenoFinder.args = list(separator = separator,
    vectorDistFun = vectorWeightedDist), outlierFinder.expr.args = list(bonf.prob = 0.5,
    transFun = atanh, tail = "upper"), outlierFinder.pheno.args = list(normal.upper.thresh = 0.99,
    bonf.prob = NULL, tail = "upper"), smokingGunFinder.args = list(transFun = I),
    impute.knn.args = list(k = 10, rowmax = 0.5, colmax = 0.8,
        maxp = 1500, rng.seed = 362436069), manual.smokingguns = NULL,
    automatic.smokingguns = FALSE, within.datasets.only = FALSE,
    intermediate.pruning = FALSE, cache.dir = "cache", verbose = TRUE)}
\arguments{
  \item{esets}{a list of ExpressionSets, containing the numeric and phenotypic data to be analyzed.}
  \item{separator}{a delimitor to use between dataset names and sample names}
  \item{corFinder.args}{a list of arguments to be passed to the corFinder function.}
  \item{phenoFinder.args}{a list of arguments to be passed to the phenoFinder function.  If
NULL, samples with similar phenotypes will not be searched for.}
  \item{outlierFinder.expr.args}{a list of arguments to be passed to outlierFinder when called for expression data}
  \item{outlierFinder.pheno.args}{a list of arguments to be passed to outlierFinder when called for phenotype data}
  \item{smokingGunFinder.args}{a list of arguments to be passed to smokingGunFinder}
  \item{impute.knn.args}{a list of arguments to be passed to impute::impute.knn.  Set to
NULL to do no knn imputation.}
  \item{manual.smokingguns}{a character vector of phenoData columns that, if identical, will
be considered evidence of duplication}
  \item{automatic.smokingguns}{automatically look for "smoking guns."  If TRUE, look for
phenotype variables that are unique to each patient in dataset 1,
also unique to each patient in dataset 2, but contain exact
matches between datasets 1 and 2.}
  \item{within.datasets.only}{If TRUE, only search within each dataset for doppelgangers.}
  \item{intermediate.pruning}{The default setting FALSE will result in output with no missing
values, but uses extra memory because all results from the
expression, phenotype, and smoking gun doppelganger searches must
be saved until the end.  Setting this to TRUE will save memory for
very large searches, but distance metrics will only be available
if that value was identified as a doppelganger (for example,
phenotype doppelgangers will have missing values for the
expression and smoking gun similarity).}
  \item{cache.dir}{The name of a directory in which to cache or look up results to save
re-calculating correlations.  Set to NULL for no caching.}
  \item{verbose}{Print progress information}
}

\value{Returns an object of S4-class "DoppelGang".  See ?DoppelGang-class.}

\author{Levi Waldron and Markus Riester}




\examples{
        library(curatedOvarianData)
        data(GSE32062.GPL6480_eset)
        data(GSE32063_eset)
        data(GSE12470_eset)
        data(GSE17260_eset)
        testesets <- list(JapaneseA=GSE32062.GPL6480_eset,
                          JapaneseB=GSE32063_eset,
                          Yoshihara2009=GSE12470_eset,
                          Yoshihara2010=GSE17260_eset)
        testesets <- lapply(testesets, function(X){
            ## standardize the sample ids to improve matching based on clinical annotation
            X$alt_sample_name <- paste(X$sample_type, gsub("[^0-9]", "", X$alt_sample_name), sep="_")
            pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
            X <- X[, 1:20]  ##speed computations
            return(X) })
        results1 <- doppelgangR(testesets, cache.dir=NULL)
        results1
        plot(results1)
        summary(results1)
    ## Set phenoFinder.args=NULL to ignore similar phenotypes, and
    ## turn off ComBat batch correction:
##    results2 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=FALSE), phenoFinder.args=NULL, cache.dir=NULL)
##    summary(results2)
}
