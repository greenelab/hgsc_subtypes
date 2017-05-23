#!/bin/bash

#exec &>hgsc_analysis.out 

############################################
# Based on shell script from:
# Cross-population analysis of high-grade serous ovarian cancer does
# not support four subtypes
#
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script stores instructions to reproduce the HGSC subtyping analysis
# across populations. All scripts and relevant files are included in the 
# repository and the workflow depends on the running sequential scripts within 
# the larger folder structure. See the README for general information and 
# INSTALL.R for package dependencies.
# This script is run using a Docker image 
# (see <https://hub.docker.com/r/gregway/hgsc_subtypes/>)
# The actual 
# ~~~~~~~~~~~~~~~~~~~~~
############################################
# INSTALL DEPENDENCIES  
Rscript INSTALL.R

#################
# PART ZERO: 
# Download Mayo data
#################
# COMBAT adjust Mayo data
Rscript 1.DataInclusion/Scripts/processMayoEset/\
Agilent1and2and3_COMBAT_datamerge.R

# Create an eset from the Mayo data
# NOTE: This requires the Normalizer function from the Sleipnir library 
# (http://libsleipnir.bitbucket.org/)
Rscript 1.DataInclusion/Scripts/processMayoEset/createMayoEset.R

# Define Constants
DATASETS="TCGA_eset mayo.eset GSE32062.GPL6480_eset GSE9891_eset aaces.eset"
KMIN=2
KMAX=4
SEED=123
NSTARTS=20
NO_SHUFFLE=FALSE
SHUFFLE=TRUE
SAM_SUBSET='commongenes'
AACES_PATH='aaces_expression.tsv'

#################
# PART ONE: 
# Dataset Selection and Inclusion
#################
# ~~~~~~~~~~~~~~~~~~~~~
# This section will determine which samples meet a specific inclusion criteria 
# for use in downstream analyses
# ~~~~~~~~~~~~~~~~~~~~~

# We are using data from curatedOvarianData version 1.8.0
# NOTE: The Mayo Clinic Data is not currently in curatedOvarianData.

# Output the samples for each dataset that pass the inclusion criteria
#
R --no-save --args $AACES_PATH 1.DataInclusion/Scripts/A.getInclusion.R  # (Table 1)

# Output the common genes and the MAD (Median Absolute Deviation) genes to be 
# used in developing moderated t score vectors and in clustering, respectively. 
# This script will also output Venn diagrams for visualizing overlapping genes 
# (Sup. Fig. S1) #NOTE: Bonome (GSE12672) is removed following the across 
# dataset correlations analysis. Add it here.
R --no-save --args $DATASETS "GSE26712_eset" < 1.DataInclusion/Scripts/\
B.getGenes.R

#################
# PART TWO: 
# Run k means and SAM
#################
# ~~~~~~~~~~~~~~~~~~~~~
# The scripts will take as inputs the samples and genes from the previous 
# section that passed the inclusion criteria. It will also run k means for 
# k min - k max, output several figures (moderated t score heatmaps, kmeans bar 
# chart distributions, correlation matrices) and tables (cluster membership by 
# dataset, within dataset cluster correlations)
# ~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~
# SAM with MAD genes
# ~~~~~~~~~~~~~
# Output across dataset correlations for MAD genes
# NOTE: common genes used in downstream analyses
R --no-save --args $KMIN $KMAX $NSTARTS $SEED FALSE $NO_SHUFFLE "madgenes" \
$DATASETS "GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R

# ~~~~~~~~~~~~~
# k means & SAM (with common genes)
# ~~~~~~~~~~~~~
# Perform k means and SAM (Figure 1)
R --no-save --args $KMIN $KMAX $NSTARTS $SEED FALSE $NO_SHUFFLE $SAM_SUBSET \
$DATASETS "GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R

# Output correlation matrices (Sup. Fig. S2)	
R --no-save --args $KMIN $KMAX $SEED Figures/CorrelationMatrix/ $DATASETS \
"GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/B.CorrelationMatrix.R

