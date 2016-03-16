############################################
# Cross-population analysis of high-grade serous 
# ovarian cancer does not support four subtypes
#
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will ensure reproducibility of workflow
############################################

mirror <- "http://cran.us.r-project.org"

######################
# Install CRAN packages 
######################
cran_pkgs <- c(
  'plyr',
  'gridExtra',
  'reshape',
  'reshape2',
  'scales',
  'proto', 
  'stringr',
  'mnormt',
  'sfsmisc',
  'ggplot2',
  'cluster',
  'devtools',
  'GGally',
  'gplots',
  'NMF',
  'outliers',
  'survival'
)

install.packages(cran_pkgs, repos = mirror)

######################
# Install bioconductor packages
######################
source("https://bioconductor.org/biocLite.R")
bioc_pkgs <- c(
  'affy',
  'Biobase',
  'BiocParallel',
  'sva',
  'impute',
  'siggenes',
  'limma'
)

biocLite(bioc_pkgs, suppressUpdates = TRUE)

######################
# Install source packages
######################
# doppelgangR
customLib <- "hgsc_subtypes/1.DataInclusion/doppelgangR-master/"
install.packages(customLib, repos = NULL, type = "source")

# curatedOvarianData
customLib <- "curatedOvarianData"
install.packages(customLib, repos = NULL, type = "source")
