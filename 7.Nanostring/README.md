# Gene Expression Correlations with Nanostring Subtype Classifier

**Gregory Way and Casey Greene, 2018**

Previous research identified a sparse classifier of few genes that can classify High Grade Serous Ovarian Cancer (HGSC) subtypes.
This classifier (Lasso) was selected to induce sparsity in the solutions so that the identified genes can be placed on a [nanostring](https://www.nanostring.com/) panel.

Here, we test the correlation of these genes against all other genes measured across four HGSC gene expression datasets.
The datasets include TCGA, Tothill, Yoshihara, and Mayo.
We use [curatedOvarianData](https://bioconductor.org/packages/release/data/experiment/html/curatedOvarianData.html) to access datasets.

## Procedure

We take the 29 nanostring genes and consider their gene expression vectors across the four datasets independently.
We then take pairwise correlations of each of these 29 genes against **all** other genes for each dataset.
Lastly, we consider the genes in the top 5% and top 1% of these correlations for all 29 genes for each dataset.
The final output is a long data frame (`results/all_threshold_classifier_gene_correlations.tsv`) that stores the highest correlated genes against the nanostring classifier genes.

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
| `D.get_overlap_genes.R` | A final summary dataframe sorted by highest correlated genes with the highest support (99% threshold across all 4 datasets) |

## Output DataFrame

The file `results/all_threshold_classifier_gene_correlations.tsv` contains 57,551 rows describing correlations of HGSC gene expression values against 29 nanostring classifier genes.
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

