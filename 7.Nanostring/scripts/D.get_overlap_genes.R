# Nanostring Classifier Genes
# Gregory Way 2018
# 7.Nanostring/scripts/D.get_overlap_genes.R
#
# Determine gene candidates that are consistently strongly correlated with
# classifier genes' expression across different HGSC datasets
#
# Output:
# A single long dataframe consisting of classifier gene, correlated gene,
# range of correlations across datasets, whether it exists in 99% or 95%
# threshold, and how many and what percentage of datasets are included in the
# overlap.

library(dplyr)

base_file <- file.path("7.Nanostring", "results")

# Load thresholded dataframes
file <- file.path(base_file, "rf_high_99_quantile_correlated_genes.tsv")
high_df <- readr::read_tsv(file)
high_df$threshold <- "high_thresh"

file <- file.path(base_file, "rf_relaxed_95_quantile_correlated_genes.tsv")
relaxed_df <- readr::read_tsv(file)
relaxed_df$threshold <- "relaxed_thresh"

# Not all classifier genes were measured by each platform - find out how many
classifier_gene_per_dataset <- high_df %>% 
  dplyr::group_by(classifier_gene, dataset) %>% 
  dplyr::summarise(num_genes = n()) %>% 
  dplyr::mutate(num_datasets = n()) %>% 
  dplyr::select(classifier_gene, num_datasets) %>% 
  dplyr::distinct()

# Combine datasets
all_cor_df <- dplyr::bind_rows(high_df, relaxed_df) %>% 
  dplyr::left_join(classifier_gene_per_dataset, by = 'classifier_gene')

# Wrangle classifier data into interpretable output dataframe and output
f <- file.path(base_file, "rf_all_thresholded_classifier_gene_correlations.tsv")

all_cor_df %>%
  # Collect groups of genes and identifiers across datasets
  dplyr::group_by(classifier_gene, gene, threshold, num_datasets) %>%

  # Determine count and the extent the correlation was observed across datasets
  dplyr::summarise(num_datasets_gene = n(),
                   min_gene_cor = min(gene_cor),
                   max_gene_cor = max(gene_cor)) %>%

  # Create a new variable determining the percentage of dataset observations
  dplyr::mutate(percent_datasets = (num_datasets_gene / num_datasets)) %>%

  # Cast the dataframe stratified by threshold and display membership (length)
  reshape2::dcast(classifier_gene + gene + num_datasets + percent_datasets +
                    min_gene_cor + max_gene_cor ~
                    threshold,
                  fun.aggregate = function(x) length(x),
                  value.var = "classifier_gene") %>%

  # Convert back to a tibble object
  dplyr::as.tbl() %>%
  
  # Sort tibble by confidence in correlation
  dplyr::arrange(desc(percent_datasets),
                 desc(high_thresh),
                 desc(max_gene_cor)) %>%
  
  # Write the output file
  readr::write_tsv(f)
