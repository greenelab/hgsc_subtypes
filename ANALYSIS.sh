#!/bin/bash

############################################
# Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes 
#
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script stores instructions to reproduce the HGSC subtyping analysis
# across populations. All scripts and relevant files are included in the repository
# and the workflow depends on the running sequential scripts within the larger folder
# structure sequence. See the README for general information and package dependencies.
# ~~~~~~~~~~~~~~~~~~~~~
############################################

#################
# PART ZERO: 
# Install Packages and Download Mayo data
#################
Rscript INSTALL.R  # Note: refer to the README for version control

# COMBAT adjust Mayo data
Rscript 1.DataInclusion/Scripts/processMayoEset/Agilent1and2and3_COMBAT_datamerge.R

# Create an eset from the Mayo data
# NOTE: This requires the Normalizer function from the Sleipnir library (http://libsleipnir.bitbucket.org/)
Rscript 1.DataInclusion/Scripts/processMayoEset/createMayoEset.R

# Define Constants
DATASETS="TCGA_eset Mayo GSE32062.GPL6480_eset GSE9891_eset"
KMIN=3
KMAX=4
SEED=123
NSTARTS=20

#################
# PART ONE: 
# Dataset Selection and Inclusion
#################
# ~~~~~~~~~~~~~~~~~~~~~
# This section will determine which samples meet a specific inclusion criteria for use in downstream analyses
# ~~~~~~~~~~~~~~~~~~~~~

# We are using data from curatedOvarianData version 1.3.4 
# NOTE: The Mayo Clinic Data is not currently in curatedOvarianData.

# Output the samples for each dataset that pass the inclusion criteria
R --no-save --args < 1.DataInclusion/Scripts/A.getInclusion.R  # (Table 1)

# Output the common genes and the MAD (Median Absolute Deviation) genes to be used in developing
# moderated t score vectors and in clustering, respectively. This script will also output Venn diagrams
# for visualizing overlapping genes (Sup. Fig. S1) #NOTE: Bonome (GSE12672) is removed following the
# across dataset correlations analysis. Add it here.
R --no-save --args $DATASETS "GSE26712_eset" < 1.DataInclusion/Scripts/B.getGenes.R

#################
# PART TWO: 
# Run k means and SAM
#################
# ~~~~~~~~~~~~~~~~~~~~~
# The scripts will take as inputs the samples and genes from the previous section that passed the inclusion 
# criteria. It will also run k means for k min - k max, output several figures (moderated t score heatmaps, ggpairs 
# scatterplots, kmeans bar chart distributions, correlation matrices) and tables (cluster membership by dataset,
# within dataset cluster correlations)
# ~~~~~~~~~~~~~~~~~~~~~

	#################
	#~k means & SAM~#
	#################
	# Perform k means and SAM (args: kmin, kmax, nstarts, seed, bNMF) (Figure 3)
	R --no-save --args $KMIN $KMAX $NSTARTS $SEED FALSE $DATASETS "GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R

	# Output correlation matrices (Figure 1)	
	R --no-save --args $KMIN $KMAX $SEED Figures/CorrelationMatrix/ $DATASETS "GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/B.CorrelationMatrix.R

	# Output k-means barcharts (Figure 2)
	R --no-save --args $KMIN $KMAX $DATASETS < 2.Clustering_DiffExprs/Scripts/C.KMeansBarCharts.R

	#################
	#~~~~~~NMF~~~~~~#
	#################
	# Output consensus matrices, NMF cluster membership files (Sup. Figure S2) and cophenetic coefficients (Sup. Figures S3-S7)
	R --no-save --args $KMIN $KMAX $NSTARTS $SEED $DATASETS "GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/D.NMF.R

	# Run SAM on NMF clusters (TRUE argument forces NMF analysis)
	R --no-save --args $KMIN $KMAX $NSTARTS $SEED TRUE $DATASETS "GSE26712_eset" < 2.Clustering_DiffExprs/Scripts/A.run_kmeans_SAM.R

	#################
	#k means vs. NMF#
	#################
	# Compare k-means defined clusters with NMF defined clusters(Figure 4)
	R --no-save --args $DATASETS < 2.Clustering_DiffExprs/Scripts/E.kmeans_v_nmf.R

	# Compile table with all cluster membership information (Sup. Table S2)
	R --no-save --args $DATASETS < 2.Clustering_DiffExprs/Scripts/F.clusterMembership.R

	#################
	#Other Clusters~#
	#################
	# Investigate the similarities in cluster membership in original TCGA 2011 paper, the Konecny 2014 paper, 
	# and the Tothill 2008 paper (Table 4)
	R --no-save < 2.Clustering_DiffExprs/Scripts/G.Dataset_concordance.R
        
        #################
        #  Tothill LMP  #
        #################
        # Observe consensus matrices and cophenetic coefficients for Tothill dataset if LMP samples are not removed.
        # This is similar to the results presented by TCGA supplementary figure S6.2
        R --no-save < 2.Clustering_DiffExprs/Scripts/H.TCGA_SupFig6.2_TothillPrediction.R

#################
# PART THREE: 
# Goodness of Fit
#################
# ~~~~~~~~~~~~~~~~~~~~~
# The section will output several figures, all different goodness of fit metrics, for each dataset
# The metrics include AIC, BIC, Gap Statistic, and Silhouette Widths (Note: Cophenetic obtained in step 2)
# ~~~~~~~~~~~~~~~~~~~~~

# Perform AIC, BIC, and silhouette width analyses
R --no-save --args 2 8 20 $NSTARTS $SEED $DATASETS < 3.Fit/Scripts/A.GoodnessFit.R

# Determine gap statistic
# NOTE: This step was performed on the Discovery Cluster at Dartmouth College
R --no-save --args 8 250 20 50 $SEED $DATSETS < 3.Fit/Scripts/B.GAP_GoodnessFit.R

#################
# PART FOUR: 
# Survival
#################
# ~~~~~~~~~~~~~~~~~~~~~
# This section will perform all survival analyses and output summary tables
# ~~~~~~~~~~~~~~~~~~~~~

# Output Kaplan-Meier survival curves, and perform a cox proportional hazards regression model
R --no-save --args $DATASETS < 4.Survival/Scripts/A.Survival.R

# Summarize the results of the survival analysis
R --no-save --args $DATASETS < 4.Survival/Scripts/B.Summarize_Survival.R

#################
# PART FIVE: 
# Gene and pathway Analyses
#################
# ~~~~~~~~~~~~~~~~~~~~~
# This section will describe the cluster driving genes in common across populations
# ~~~~~~~~~~~~~~~~~~~~~

# Output tables of cluster specific genes based on significant SAM FDR values (Sup. Table S4)
R --no-save --args $DATASETS < 5.Pathway/Scripts/A.GeneEnrichment.R

# The output of this script is input into a PANTHER pathways analysis (http://pantherdb.org/)
