############################################
# Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes
#
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will install all CRAN and Bioconductor packages that the workflow requires

####################################
# Install Libraries
####################################
cranLib <- c("reshape2", "outliers", "grid", "gplots", "RColorBrewer", "ggplot2", "cluster", "NMF", 
             "GGally", "sfsmisc", "survival", "plyr", "Rcpp", "phylobase", "devtools")
bioCLib <- c("curatedOvarianData", "siggenes", "BiocParallel", "sva", "impute", "limma")
customLib <- "1.DataInclusion/doppelgangR-master/"

# Install CRAN libraries
install.packages(cranLib, repos = "http://cran.us.r-project.org")

# Install Bioconductor libraries
source("http://bioconductor.org/biocLite.R")
biocLite(bioCLib)

# Install custom package doppelgangR from source file included in the repository
install.packages(customLib, repos = NULL, type = "source")

####################################
# Load Libraries
####################################
library(curatedOvarianData)
library(doppelgangR)
library(reshape2)
library(outliers)
library(limma)
library(grid)
library(curatedOvarianData)
library(gplots)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(cluster)
library(NMF)
library(siggenes)
library(GGally)
library(sfsmisc)
library(survival)
library(mnormt)
library(phylobase)
library(devtools)

sessionInfo()
