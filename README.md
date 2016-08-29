# High-Grade Serous Ovarian Cancer Subtypes - Why has the field settled on four?

**(C) Trustees of the University of Pennsylvania**

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32906.svg)](http://dx.doi.org/10.5281/zenodo.32906)

## Summary

The repository contains instructions to replicate the analysis of identifying
high-grade serous ovarian cancer subtypes across Australian, American, and
Japanese populations. We leverage data from a variety of studies extracted from
the bioconductor package curatedOvarianData
([Ganzfried et al. 2013](http://doi.org/10.1093/database/bat013)) as well as a
 dataset we uploaded to GEO (GSE74357). We apply a unified, unsupervised
bioinformatics pipeline to compare subtypes across these populations and
determine that specific subtypes are reliably identified. The most replicable
subtypes are mesenchymal-like and proliferative-like and their sample
representation was highly concordant with other independent clustering studies
performed on single populations.

## Contact

For all analysis or coding related questions please file a
[GitHub issue](https://github.com/greenelab/hgsc_subtypes/issues)

Please direct all other correspondence to:
Jennifer.A.Doherty@dartmouth.edu and/or CSGreene@mail.med.upenn.edu


## Analysis

For ease of use and to ensure reproducibility, all analyses should be performed
in our [Docker image](https://hub.docker.com/r/gregway/hgsc_subtypes/)
To install docker, please follow these
[friendly instructions](https://docs.docker.com/linux/)

After installing Docker, please follow these steps:

```sh
docker pull gregway/hgsc_subtypes
git clone https://github.com/greenelab/hgsc_subtypes.git

# Run the analysis and save results to local volume
docker run \
--detach \
--volume ~/<PATH>/hgsc_subtypes/:hgsc_subtypes/ \
--publish 5000:80
--interactive gregway/hgsc_subtypes bash hgsc_subtypes/docker/docker_command.sh

```

Where PATH is the local directory structure where the github repo is cloned. The
`docker run` command will output the container ID of the analysis. Runtime is on
the order of several hours.

## Full Build

To retrieve the full build without installing docker download our 
[Zenodo record](https://zenodo.org/record/53990)

## Data
All data was retrieved from curatedOvarianData except for the Mayo data
(GSE74357).

## Dependencies

All dependencies are pre-installed in the Docker image. 

* For specific R package installations, view INSTALL.R. 
* The analysis also requires the Sleipnir Normalizer function

## Acknowledgements

This work was supported by the Institute for Quantitative Biomedical Sciences
(Dartmouth); The graduate program in Genomics and Computational Biology (Penn);
The Norris Cotton Cancer Center Developmental Funds;
the National Cancer Institute at the National Institutes of Health (R01 CA168758
to J.A.D., F31 CA186625 to J.R., R01 CA122443 to E.L.G.); The Mayo Clinic
Ovarian Cancer SPORE (P50 CA136393 to E.L.G.); The Mayo Clinic Comprehensive
Cancer Center-Gene Analysis Shared Resource (P30 CA15083); The Gordon and Betty
Moore Foundation’s Data-Driven Discovery Initiative (grant number GBMF 4552 to
C.S.G.); and The American Cancer Society (grant number IRG 8200327 to C.S.G.).
