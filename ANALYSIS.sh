#!/bin/bash

exec &>hgsc_analysis.out 

############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support
# four subtypes
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
# ~~~~~~~~~~~~~~~~~~~~~
############################################

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
DATASETS="TCGA_eset Mayo GSE32062.GPL6480_eset GSE9891_eset"
KMIN=2
KMAX=4
SEED=123
NSTARTS=20
NO_SHUFFLE=FALSE
SHUFFLE=TRUE

#################
# PART ONE: 
# Dataset Selection and Inclusion
#################
# ~~~~~~~~~~~~~~~~~~~~~
# This section will determine which samples meet a specific inclusion criteria 
# for use in downstream analyses
# ~~~~~~~~~~~~~~~~~~~~~

# We are using data from curatedOvarianData version 1.8
# NOTE: The Mayo Clinic Data is not currently in curatedOvarianData.

# Output the samples for each dataset that pass the inclusion criteria
R --no-save --args < 1.DataInclusion/Scripts/A.getInclusion.R  # (Table 1)

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
# k means & SAM
# ~~~~~~~~~~~~~
# Perform k means and SAM (args: kmin, kmax, nstarts, seed, bNMF) (Figure 1)
R --no-save --args $KMIN $KMAX $NSTARTS $SEED FALSE $NO_SHUFFLE $DATASETS \
"GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R

# Output correlation matrices (Sup. Fig. S2)	
R --no-save --args $KMIN $KMAX $SEED Figures/CorrelationMatrix/ $DATASETS \
"GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/B.CorrelationMatrix.R

# Output k-means barcharts (Figure S8)
R --no-save --args $KMIN $KMAX $DATASETS < 2.Clustering_DiffExprs/Scripts/\
C.KMeansBarCharts.R

# Shuffle genes to compare across population correlations in real data
R --no-save --args $KMIN $KMAX $NSTARTS $SEED FALSE $SHUFFLE $DATASETS \
"GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R

# ~~~~~~~~~~~~~
# NMF
# ~~~~~~~~~~~~~
# Output consensus matrices, NMF cluster membership files (Sup. Figure S2) and 
# cophenetic coefficients (Sup. Figures S3-S7)
R --no-save --args $KMIN $KMAX $NSTARTS $SEED $DATASETS "GSE26712_eset" \
< 2.Clustering_DiffExprs/Scripts/D.NMF.R

# Run SAM on NMF clusters (TRUE argument forces NMF analysis)
R --no-save --args $KMIN $KMAX $NSTARTS $SEED TRUE $NO_SHUFFLE $DATASETS \
"GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R

# ~~~~~~~~~~~~~
# k means vs. NMF
# ~~~~~~~~~~~~~
# Compare k-means defined clusters with NMF defined clusters (Figure 2)
R --no-save --args $DATASETS < 2.Clustering_DiffExprs/Scripts/E.kmeans_v_nmf.R

# Compile table with all cluster membership information (Sup. Table S2)
R --no-save --args $DATASETS < 2.Clustering_DiffExprs/Scripts/\
F.clusterMembership.R

# ~~~~~~~~~~~~~
# Dataset Concordance
# ~~~~~~~~~~~~~
# Investigate the similarities in cluster membership in original TCGA 2011 
# paper, the Konecny 2014 paper, and the Tothill 2008 paper (Table 4)
R --no-save < 2.Clustering_DiffExprs/Scripts/G.Dataset_concordance.R

# ~~~~~~~~~~~~~
#  Tothill LMP  #
# ~~~~~~~~~~~~~
# Observe consensus matrices and cophenetic coefficients for Tothill dataset if 
# LMP samples are not removed. This is similar to the results presented by TCGA 
# supplementary figure S6.2 (Figure 3)
R --no-save < 2.Clustering_DiffExprs/Scripts/H.TCGA_LMP_TothillPrediction.R

#################
# PART THREE: 
# Goodness of Fit
#################
# ~~~~~~~~~~~~~~~~~~~~~
# The section will output several figures, all different goodness of fit 
# metrics, for each dataset. The metrics include AIC, BIC, Gap Statistic, and 
# Silhouette Widths (Note: Cophenetic obtained in step 2)
# ~~~~~~~~~~~~~~~~~~~~~

# Perform AIC, BIC, and silhouette width analyses
R --no-save --args 2 8 20 $NSTARTS $SEED $DATASETS < \
3.Fit/Scripts/A.GoodnessFit.R

# Determine gap statistic
# NOTE: This step was performed on the Discovery Cluster at Dartmouth College
R --no-save --args 8 250 20 50 $SEED $DATSETS < 3.Fit/Scripts/\
B.GAP_GoodnessFit.R

#################
# PART FOUR: 
# Survival
#################
# ~~~~~~~~~~~~~~~~~~~~~
# This section will perform all survival analyses and output summary tables
# ~~~~~~~~~~~~~~~~~~~~~

# Output Kaplan-Meier survival curves, and perform a cox proportional hazards 
# regression model (Sup. Fig. S11)
R --no-save --args $DATASETS < 4.Survival/Scripts/A.Survival.R

# Summarize the results of the survival analysis (Sup. Table S3)
R --no-save --args $DATASETS < 4.Survival/Scripts/B.Summarize_Survival.R

#################
# PART FIVE: 
# Gene and pathway Analyses
#################
# ~~~~~~~~~~~~~~~~~~~~~
# Describe the cluster driving genes in common across populations
# ~~~~~~~~~~~~~~~~~~~~~

# Output tables of cluster specific genes based on significant SAM FDR values
R --no-save --args $DATASETS < 5.Pathway/Scripts/A.GeneEnrichment.R

# The output of this script is input into a PANTHER pathways analysis
# (http://pantherdb.org/)
