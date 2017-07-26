# High-Grade Serous Ovarian Cancer Subtypes - Why has the field settled on four?

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32906.svg)](http://dx.doi.org/10.5281/zenodo.32906)

## Summary

In this repository, we compare high-grade serous ovarian cancer (HGSC) subtypes
across Australian, American, and Japanese populations. We determine that two or
three subtypes are most consistent across different datasets. A full report of
this analysis is published in _G3: Genes, Genomes, Genetics_
([Way et al. 2016](https://doi.org/10.1534/g3.116.033514)).
Instructions are provided in
[release version 1.3](https://github.com/greenelab/hgsc_subtypes/tree/1.3)
to reproduce the analysis.

We leverage data extracted from the bioconductor package `curatedOvarianData`
([Ganzfried et al. 2013](http://doi.org/10.1093/database/bat013)) as well as a
 dataset we uploaded to GEO 
([GSE74357](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74357)). We
apply a unified, unsupervised bioinformatics pipeline to compare subtypes across
these populations and determine that specific subtypes are reliably identified.
The most replicable subtypes are mesenchymal-like and proliferative-like and
their sample representation was highly concordant with other independent
clustering studies performed on single populations.

We are currently working on adding African American HGSC samples to this
pipeline to determine the representation of HGSC subtypes in an additional
population. This project is in development and will be associated with a
future release.

## Contact

For all analysis or coding related questions please file a
[GitHub issue](https://github.com/greenelab/hgsc_subtypes/issues)

## Environment

To ensure analysis reproducibility, all R/Bioconductor/Python packages are
versioned using a combination of
[checkpoint](https://cran.r-project.org/web/packages/checkpoint/index.html)
and [bioconda](https://github.com/bioconda/bioconda-recipes).

To create an instance of this environment run the following:

```sh
conda env create --force --file environment.yml
source activate hgsc_subtypes

Rscript INSTALL.R
```

## Analyses

There are currently two pipelines in place to analyze hgsc subtypes, to
reproduce the results of either pipeline, run:

```sh
# Cross-population HGSC subtypes analysis 
bash hgsc_subtypes_pipeline.sh

# African American HGSC subtypes analysis
bash aaces_subtypes_pipeline.sh
```

## Data

All data was retrieved from curatedOvarianData except for the
[Mayo data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74357)
and AACES data.

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
