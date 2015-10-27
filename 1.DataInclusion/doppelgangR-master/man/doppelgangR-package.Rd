\name{doppelgangR-package}
\alias{doppelgangR-package}
\docType{package}
\title{Identify likely duplicate samples from genomic or meta-data.}
\description{The main function is doppelgangR(), which takes as minimal
input a list of ExpressionSet object, and searches all list pairs for
duplicated samples.  The search is based on the genomic data (exprs(eset)),
phenotype/clinical data (pData(eset)), and "smoking guns" - supposedly unique
identifiers found in pData(eset).}
\details{
\tabular{ll}{Package: \tab doppelgangR\cr
Maintainer: \tab Levi Waldron <levi.waldron@hunter.cuny.edu>\cr
Depends: \tab R (>= 3.0), affy, methods, Biobase, BiocParallel\cr
Suggests: \tab curatedOvarianData, knitr, ROCR, pROC\cr
Imports: \tab sva, impute, digest, mnormt\cr
Author: \tab Levi Waldron and Markus Riester\cr
Version: \tab 0.10.0\cr
Date: \tab 2015-01-27\cr
License: \tab GPL (>=2.0)\cr
Title: \tab Identify likely duplicate samples from genomic or meta-data.\cr
URL: \tab https://github.com/lwaldron/doppelgangR\cr
BugReports: \tab https://github.com/lwaldron/doppelgangR/issues\cr
VignetteBuilder: \tab knitr\cr}
}
\author{Levi Waldron and Markus Riester}

\keyword{ package }


