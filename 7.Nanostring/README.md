# Gene Expression Correlations with Nanostring Subtype Classifier

**Gregory Way, Casey Greene, and Jennifer Doherty 2018**

Previous research identified a sparse classifier of few genes that can classify High Grade Serous Ovarian Cancer (HGSC) subtypes.
This classifier (Random Forest) was used to find a sparse set of genes to be placed on a [nanostring](https://www.nanostring.com/) panel.

Here, we test the correlation of these genes against all other genes measured across four HGSC gene expression datasets.
The datasets include TCGA, Tothill, Yoshihara, and Mayo.
We use [curatedOvarianData](https://bioconductor.org/packages/release/data/experiment/html/curatedOvarianData.html) to access datasets.

## Procedure

### Part I

We take the 59 nanostring genes and consider their gene expression vectors across the four datasets independently.
We then take pairwise Pearson correlations of each of these 59 genes against **all** other genes for each dataset.
Lastly, we consider the genes in the top 5% and top 1% of these correlations for all 59 genes for each dataset.
A useful output is the long data frame (`results/all_threshold_classifier_gene_correlations.tsv`), which stores the highest correlated genes against the nanostring classifier genes.

### Part II

In addition, we output geneset files ([`.gmt` format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)) under 4 different tiers of confidence (see below).
The file is located in `results/correlated_hgsc_classifier_genes.gmt`.
The format of each line in this tab separated file is: [classifier gene and confidence tier, correlated geneset (1, 2, ..., n)].
The first gene in the correlated geneset is the classifier gene itself.

### Part III

We then run a geneset overrepresentation analysis (ORA) using `Tier 1B` genes using all unique genes measured in the 4 datasets as background.
We use [WebGestalt](https://doi.org/10.1093/nar/gkx356 "WebGestalt 2017: a more comprehensive, powerful, flexible and interactive gene set enrichment analysis toolkit") for ORA.
The results of the pathway analysis are stored in `results/overrepresented_pathways_classifier_genes.tsv`.

## Computational Environment

All processing and analysis scripts were performed using the conda environment specified in `environment.yml`.
To build and activate this environment run:

```bash
# conda version 4.5.1
conda env create --force --file environment.yml

conda activate hgsc-nanostring
```

## Reproduce Pipeline

Each of the scripts are designed to be executed sequentially.
The scripts are located in the `scripts/` folder and include:

| Script | Output |
| :----- | :----- |
| `A.get_correlation_output.R` | (A) Correlations between HGSC genes and subtype classifier genes (B) Missing classifier genes by dataset |
| `B.explore_correlations.py` | (A) Distribution of correlations across classifier genes and datasets (B) 95% and 99% threshold correlation dataframes in long format |
| `C.threshold_venns.R` | Venn diagrams describing correlated genes across datasets for both thresholds and against all 29 classifier genes |
| `D.get_overlap_genes.R` | A summary dataframe sorted by highest correlated genes with the highest support (99% threshold across all 4 datasets) |
| `E.gmt_genesets_pathway_analysis.R` | A single `.gmt` file for use in downstream gene set enrichment-like analyses |
| `F.summarize_pathways.R` | A single table describing top 15 enriched gene ontology terms for each classifier gene |

## Output

### All thresholded gene correlations

The file `results/all_threshold_classifier_gene_correlations.tsv` contains 57,551 rows describing correlations of HGSC gene expression values against 59 nanostring classifier genes.
A description of the columns are:

1. `classifier_gene` - The gene used in the nanostring subtype classifier
2. `gene` - The gene compared against in the HGSC gene expression datasets
3. `num_datasets` - How many datasets the nanostring classifier gene is measured in
4. `percent_datasets` - Of `num_datasets` how many times is the gene correlation pair observed
5. `min_gene_cor` - The minimum correlation observed for the gene by gene (for those not filtered by threshold)
6. `max_gene_cor` - The maximum correlation observed for the gene by gene (for those not filtered by threshold)
7. `95_threshold` - Indicator variable if the gene by gene correlation was observed in the top 5% of correlations
8. `99_threshold` - Indicator variable if the gene by gene correlation was observed in the top 1% of correlations

The file is sorted by percentage of datasets encountered, if the correlation exists in the 99% threshold, and then by maximum gene correlation.
Note that the scenario `(0, 1)` is possible for the 95% threshold and 99% threshold.
This happens when a significant gene by gene correlation exists in a different number of datasets per threshold.

### Gene set file

This file is compiled into a `.gmt` file with confidence tiers.
The confidence tiers include

| Confidence Tier | Description |
| :-------------: | :---------- |
| Tier 1A | Genes in 99% threshold correlations for all 4 datasets |
| Tier 1B | Genes in 99% threshold correlations for all 3 datasets measured in (allows for a single dataset with missing measurements) |
| Tier 2A | Genes in 95% threshold correlations for all 4 datasets |
| Tier 2B | Genes in 95% threshold correlations for all 3 datasets measured in (allows for a single dataset with missing measurements) |
| Tier 3 | Genes in 95% threshold correlations for 3/4 datasets (allows for a single dataset to not have correlation) |
| Tier 4 | Genes in 95% threshold correlations for 2/4 datasets (allows for two datasets to not have correlations) |

The tier system is built into the single `.gmt`.
For high confident genesets use Tier 1A but at the cost of smaller sets.

### Overrepresented pathways

The table stores the top 15 enriched gene ontology (GO) biological process terms against correlated genesets for all classifier genes.
We use `Tier 1B` genesets as input.
A description of the columns:

1. `classifier_gene` - The gene used in the nanostring subtype classifier
2. `tier` - which geneset tier was used (in this case `tier-1B` was used)
3. `geneset` - GO identifier
4. `description` - GO term description
5. `link` - a direct link to the GO term online information
6. `PValue` - unadjusted P value for enrichment
7. `FDR` - Benjamini-Hochberg adjusted P value
8. `OverlapGene_UserID` - Which input genes were found in the given pathway

The analyses is designed so that others may make different decisions on which genesets, tiers, or pathways to use.
