############################################
# Cross-population analysis of high-grade serous
# ovarian cancer does not support four subtypes
#
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will ensure reproducibility of workflow
############################################
# USAGE: This script is to be run upon docker image build

library('methods')
library('checkpoint')

######################
# Install CRAN packages
######################
dir.create('.checkpoint')
checkpoint('2016-03-01', checkpointLocation = '.')

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
  'outliers',
  'NMF',
  'survival',
  'VennDiagram',
  'optparse'
)

install.packages(cran_pkgs)

######################
# Install ESTIMATE
######################
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos = rforge, dependencies = TRUE)

######################
# Install MCPcounter
######################
library(devtools)
install_github("ebecht/MCPcounter",ref="master", subdir="Source")
