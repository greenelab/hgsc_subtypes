#!/bin/bash
#exec &>hgsc_analysis.out

# Amy Campbell, 2017
# This shell script performs differential clustering analyses
# from Way et al's 'Cross-population analysis of high-grade serous ovarian 
# cancer does not support four subtypes'
# This script assumes you have activated the conda hgsc_subtypes environment
# encoded in environment.yml

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
# 
# We are using data from curatedOvarianData version 1.8.0
# NOTE: The Mayo Clinic Data is not currently in curatedOvarianData.

# Output the samples for each dataset that pass the inclusion criteria
Rscript 1.DataInclusion/Scripts/A.getInclusion.R --aaces $AACES_PATH

# Output the common genes and the MAD (Median Absolute Deviation) genes to be
# used in developing moderated t score vectors and in clustering, respectively.

Rscript 1.DataInclusion/Scripts/B.getGenes.R $DATASETS "GSE26712_eset"

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
Rscript 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R $KMIN $KMAX $NSTARTS $SEED FALSE $NO_SHUFFLE \
"madgenes" $DATASETS "GSE26712_eset" 

# ~~~~~~~~~~~~~
# k means & SAM (with common genes)
# ~~~~~~~~~~~~~
# Perform k means and SAM
Rscript 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R $KMIN $KMAX $NSTARTS $SEED FALSE $NO_SHUFFLE $SAM_SUBSET \
$DATASETS "GSE26712_eset" 

# Output correlation matrices
Rscript 2.Clustering_DiffExprs/Scripts/B.CorrelationMatrix.R $KMIN $KMAX $SEED Figures/CorrelationMatrix/ $DATASETS \
"GSE26712_eset"

# Output k-means barcharts
Rscript 2.Clustering_DiffExprs/Scripts/C.KMeansBarCharts.R $KMIN $KMAX $DATASETS 

# Shuffle genes to compare across population correlations in real data
Rscript 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R $KMIN $KMAX $NSTARTS $SEED FALSE $SHUFFLE $SAM_SUBSET \
$DATASETS "GSE26712_eset" 

# ~~~~~~~~~~~~~
# NMF
# ~~~~~~~~~~~~~
# Output consensus matrices, NMF cluster membership files  and
# cophenetic coefficients
Rscript 2.Clustering_DiffExprs/Scripts/D.NMF.R $KMIN $KMAX $NSTARTS $SEED $DATASETS "GSE26712_eset"

# Run SAM on NMF clusters (TRUE argument forces NMF analysis)
Rscript 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R $KMIN $KMAX $NSTARTS $SEED TRUE $NO_SHUFFLE $SAM_SUBSET \
$DATASETS "GSE26712_eset"

# ~~~~~~~~~~~~~
# k means vs. NMF
# ~~~~~~~~~~~~~
# Compare k-means defined clusters with NMF defined clusters 
Rscript 2.Clustering_DiffExprs/Scripts/E.kmeans_v_nmf.R $DATASETS  

# Compile table with all cluster membership information
Rscript 2.Clustering_DiffExprs/Scripts/F.clusterMembership.R $DATASETS 
