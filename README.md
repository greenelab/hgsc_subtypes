
############################################
# Analytical Code for "Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes"

#### Way, G., Rudd, J., Wang, C., Hamidi, H., Fridley, B., Konecny, G., Goode, E., Greene, C., Doherty, J. 

[![DOI](https://zenodo.org/badge/18957/gwaygenomics/hgsc_subtypes.svg)](https://zenodo.org/badge/latestdoi/18957/gwaygenomics/hgsc_subtypes)
############################################

#######################
# SUMMARY
#######################
The repository contains instructions to replicate the analysis of identifying high-grade serous ovarian cancer subtypes across Australian, American, and Japanese populations. We leverage data from a variety of studies extracted from the bioconductor package curatedOvarianData (Ganzfried et al. 2013) as well as a dataset we uploaded to GEO (GSE74357). We apply a unified, unsupervised bioinformatics pipeline to compare subtypes across these populations and determine that specific subtypes are reliably identified. The most replicable subtypes were mesenchymal-like and proliferative-like and their sample representation was highly concordant with other independent clustering studies performed on single populations.

#######################
# CONTACT
#######################
Please report all bugs and direct coding questions to:
GregWay@mail.med.upenn.edu

Please direct questions regarding the analysis or other correspondence to:
Jennifer.A.Doherty@dartmouth.edu

#######################
# ANALYSIS
#######################
All scripts are intended to be run in sequential order as defined by the ANALYSIS.sh file. The bash script can be run directly in the terminal to replicate the findings, but this is not recommended since there are elements that require heavy computational resources that are not suitable to be performed on a standard personal computer. We do recommend to follow the bash script line by line to replicate the findings and to reproduce the figures and tables included in the manuscript. 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e.g.
# Define Constants
DATASETS="TCGA_eset Mayo GSE32062.GPL6480_eset GSE9891_eset"
KMIN=3
KMAX=4
SEED=123
NSTARTS=100

# Perform k-means and SAM (args: kmin, kmax, nstarts, seed, bNMF) (Also, this outputs Figure 3)
R --no-save --args $KMIN $KMAX $NSTARTS $SEED $DATASETS FALSE < 2.Clustering_DiffExprs/Scripts/A.runKMeans_SAM.R
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note: many more figures and tables are included in this repository than in the manuscript.

#######################
# Data
#######################
All data was retrieved from curatedOvarianData except for the Mayo data (GSE74357).

#######################
# DEPENDENCIES
#######################
To install all required CRAN and bioconductor packages, please run the INSTALL.R script

~~~~~~~~~~
Rscript INSTALL.R
~~~~~~~~~~

#######################
# VERSIONS TESTED
#######################

# R Version
* R 3.1.2

# R Packages
* affy_1.44.0
* Biobase_2.26.0
* BiocGenerics_0.12.0 
* BiocParallel_1.0.3
* cluster_1.15.3
* doppelgangR_0.10.3
* GGally_0.5.0
* ggplot2_1.0.1
* gplots_2.17.0
* limma_3.22.1
* multtest_2.220
* NMF_0.20.6
* RColorBrewer_1.1-2
* registry_0.2
* reshape_0.8.5
* reshape2_1.4.1
* rngtools_1.2.4
* sfsmisc_1.0-27
* survival_2.38-1

# Python Version
* Python 2.7.6

# Tested on Operating System
* Ubuntu 14.04.2 LTS
* Windows 8

# Other
* Sleipner (Normalizer version 1.0)

#######################
# ACKNOWLEDGEMENTS
#######################

This work was supported by the Institute for Quantitative Biomedical Sciences; the Norris Cotton Cancer Center Developmental Funds; the National Cancer Institute at the National Institutes of Health (R01 CA168758 to J.A.D., F31 CA186625 to J.R., R01 CA122443 to E.L.G.); the Mayo Clinic Ovarian Cancer SPORE (P50 CA136393 to E.L.G.); the Mayo Clinic Comprehensive Cancer Center-Gene Analysis Shared Resource (P30 CA15083); the Gordon and Betty Moore Foundation’s Data-Driven Discovery Initiative (grant number GBMF 4552 to C.S.G.); and the American Cancer Society (grant number IRG 8200327 to C.S.G.).
