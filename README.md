
############################################
# Analytical Code for "Comprehensive cross-population analysis of high-grade serous ovarian cancer supports no more than three subtypes"

#### Way, G., Rudd, J., Wang, C., Hamidi, H., Fridley, B., Konecny, G., Goode, E., Greene, C., Doherty, J. 

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32906.svg)](http://dx.doi.org/10.5281/zenodo.32906)
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
Jennifer.A.Doherty@dartmouth.edu and/or CSGreene@mail.med.upenn.edu

#######################
# ANALYSIS
#######################
For ease of use and to ensure reproducibility, all analyses should be performed in our Docker image <https://hub.docker.com/r/gregway/hgsc_subtypes/>
To install docker, please follow the friendly instructions provided here: <https://docs.docker.com/linux/>

After installing Docker, install our docker image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$docker pull gregway/hgsc_subtypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone this github repository and run docker to perform analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$git clone https://github.com/greenelab/hgsc_subtypes.git
$docker run -d -v ~/<PATH>/hgsc_subtypes/:/hgsc_subtypes/ -p 5000:80 -i gregway/hgsc_subtypes sh hgsc_subtypes/docker/docker_command.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Where PATH is the directory structure where the github repo is cloned. Note, runtime is on the order of several hours.

####  FULL BUILD

To retrieve the full build without installing docker download from <https://zenodo.org/record/53990>
# DATA
#######################
All data was retrieved from curatedOvarianData except for the Mayo data (GSE74357).

#######################
# DEPENDENCIES
#######################
All dependencies are pre-installed in the Docker image. 

* For specific R package installations, view INSTALL.R. 
* The analysis also requires the Sleipnir Normalizer function

#######################
# ACKNOWLEDGEMENTS
#######################

This work was supported by the Institute for Quantitative Biomedical Sciences; the Norris Cotton Cancer Center Developmental Funds; the National Cancer Institute at the National Institutes of Health (R01 CA168758 to J.A.D., F31 CA186625 to J.R., R01 CA122443 to E.L.G.); the Mayo Clinic Ovarian Cancer SPORE (P50 CA136393 to E.L.G.); the Mayo Clinic Comprehensive Cancer Center-Gene Analysis Shared Resource (P30 CA15083); the Gordon and Betty Moore Foundation’s Data-Driven Discovery Initiative (grant number GBMF 4552 to C.S.G.); and the American Cancer Society (grant number IRG 8200327 to C.S.G.).
