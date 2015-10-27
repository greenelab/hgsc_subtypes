\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,DoppelGang-method}
\alias{plot.DoppelGang}
\alias{plot.doppelgangR}
\title{ Histograms of all pairwise sample correlations, showing
  identified doppelgangers.}
\description{Identified doppelgangers are shown with a red vertical line
  overlaid on a histogram of pairwise sample correlations.  One plot is
  made per pair of datasets.}
\arguments{
  \item{x}{An object of class \code{\link{DoppelGang}}}
  \item{skip.no.doppels}{(default FALSE) If TRUE, do not plot histograms where no
    doppelgangers were identified.}
  \item{plot.pair}{An optional character vector of length two, providing
    the names of two datasets.  If provided, only the comparison of
    these two datasets will be plotted.}
  \item{\dots}{ Additional arguments passed on to \code{\link{hist}}. }
}
\value{None}
\examples{
library(curatedOvarianData)
data(TCGA_eset)
data(GSE26712_eset)
## Remove some TCGA samples to speed computation:
keep.tcga <-
c("TCGA.13.2060", "TCGA.24.2290", "TCGA.25.2392", "TCGA.25.2404",
"TCGA.59.2349", "TCGA.09.2044", "TCGA.24.2262", "TCGA.24.2293",
"TCGA.25.2393", "TCGA.25.2408", "TCGA.59.2350", "TCGA.09.2045",
"TCGA.24.2267", "TCGA.59.2351", "TCGA.09.2048", "TCGA.24.2271",
"TCGA.24.2298", "TCGA.25.2398", "TCGA.59.2354", "TCGA.09.2050",
"TCGA.24.2281", "TCGA.09.2051", "TCGA.29.2428", "TCGA.09.2055",
"TCGA.24.2289", "TCGA.29.2414", "TCGA.59.2352", "TCGA.36.2532",
"TCGA.36.2529", "TCGA.36.2551", "TCGA.42.2590", "TCGA.13.2071",
"TCGA.29.2432", "TCGA.36.2537", "TCGA.36.2547", "TCGA.04.1369",
"TCGA.42.2591", "TCGA.23.2641", "TCGA.29.2434", "TCGA.36.2538",
"TCGA.36.2548", "TCGA.04.1516", "TCGA.42.2593", "TCGA.36.2549",
"TCGA.04.1644", "TCGA.13.2057", "TCGA.23.2647", "TCGA.36.2530",
"TCGA.36.2552", "TCGA.42.2587", "TCGA.13.2061", "TCGA.42.2588",
"TCGA.36.2544", "TCGA.42.2589", "TCGA.13.2066", "TCGA.61.2613",
"TCGA.61.2614", "TCGA.24.1852", "TCGA.29.1704", "TCGA.13.1819"
)
keep.tcga <- unique(c(keep.tcga, sampleNames(TCGA_eset)[1:200]))
testesets <- list(Bonome08=GSE26712_eset, TCGA=TCGA_eset[, keep.tcga])
results1 <- doppelgangR(testesets,
                        corFinder.args=list(use.ComBat=FALSE), phenoFinder.args=NULL, cache.dir=NULL)
plot(results1)
}
\author{Levi Waldron}
\section{Methods}{
\describe{
\item{\code{signature(x = "DoppelGang")}}{
  Histograms of all pairwise sample correlations, showing
  identified doppelgangers.
}
}}
\keyword{methods}

