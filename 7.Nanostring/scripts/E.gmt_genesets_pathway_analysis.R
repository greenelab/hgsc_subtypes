# Nanostring Classifier Genes
# Gregory Way 2018
# 7.Nanostring/scripts/E.gmt_genesets_pathway_analysis.R
#
# Define a tier system based on confidence of correlation between classifier
# genes and HGSC gene expression. The tier system is as follows:
#
# Tier 1A - Genes in 99% threshold correlations for all 4 datasets
# Tier 1B - Genes in 99% threshold correlations for all 3 datasets measured in
# Tier 2A - Genes in 95% threshold correlations for all 4 datasets
# Tier 2B - Genes in 95% threshold correlations for all 3 datasets measured in
# Tier 3 - Genes in 95% threshold correlations for 3/4 datasets
# Tier 4 - Genes in 95% threshold correlations for 2/4 datasets
#
# Output:
# A Gene Matrix Transposed (.gmt) file describing gene sets of correlated genes
# against random forest classifier genes. And a series of overrepresentation
# pathway analyses of these correlated genes (tier 1b) in the 
# 7.Nanostring/results/gestalt/ folder.

library(dplyr)
library(WebGestaltR)

compile_gmt <- function(tier_matrix, tier_name) {
  # Extract geneset in gmt format based on input matrix
  #
  # Arguments
  # tier_matrix - binary matrix
  #               classifier genes (rows) by correlated genes (columns)
  # tier_name - the name of the tier to store in gmt gene set names
  #
  # Output:
  # a list of genesets
  #    the first element of the geneset is the classifier gene followed by the
  #    correlated geneset
  
  geneset_list <- list()
  for (geneset_idx in 1:nrow(tier_matrix)) {
    class_gene <- rownames(tier_matrix)[geneset_idx]
    geneset_name <- paste(class_gene, tier_name, sep = "_")
    geneset <- tier_matrix[geneset_idx, ]
    geneset <- names(geneset[geneset >= 1])
    geneset_list[[geneset_name]] <- c(geneset_name, class_gene, geneset)
  }
  return(geneset_list)
}

results_dir <- file.path("7.Nanostring", "results")
f <- file.path(results_dir,
               "rf_all_thresholded_classifier_gene_correlations.tsv")
df <- readr::read_tsv(f,
                      col_types = list(.default = "c",
                                       num_datasets = readr::col_integer(),
                                       percent_datasets = readr::col_double(),
                                       min_gene_cor = readr::col_double(),
                                       max_gene_cor = readr::col_double(),
                                       high_thresh = readr::col_integer(),
                                       relaxed_thresh = readr::col_integer()))

# Collect Gene Set Tiers
# Tier 1a - is the strictest cutoff for correlated genesets -
# Genes were correlated in 4/4 datasets in the 99% percentile of correlations
tier_1a  <- df %>%
  dplyr::filter(high_thresh == 1, num_datasets == 4, percent_datasets == 1) %>%
  dplyr::group_by(classifier_gene) %>%
  reshape2::acast(classifier_gene ~ gene, value.var = "high_thresh", fill = 0)

# Tier 1b - is the second strictest cutoff
# Genes were correlated in at least 3 datasets but 100% of the time in the
# 99% percentile of correlations
# This includes the situation where a gene was not measured in one dataset
tier_1b  <- df %>%
  dplyr::filter(high_thresh == 1, num_datasets >= 3, percent_datasets == 1) %>%
  dplyr::group_by(classifier_gene) %>%
  reshape2::acast(classifier_gene ~ gene, value.var = "high_thresh", fill = 0)

# Tier 2a - Exactly as tier 1a but includes genes in the 95% percentile
tier_2a  <- df %>%
  dplyr::filter(relaxed_thresh == 1, num_datasets == 4,
                percent_datasets == 1) %>%
  dplyr::group_by(classifier_gene) %>%
  reshape2::acast(classifier_gene ~ gene, value.var = "relaxed_thresh",
                  fill = 0)
# Tier 2b - Exactly as tier 1b but includes genes in the 95% percentile
tier_2b  <- df %>%
  dplyr::filter(relaxed_thresh == 1, num_datasets >= 3,
                percent_datasets == 1) %>%
  dplyr::group_by(classifier_gene) %>%
  reshape2::acast(classifier_gene ~ gene, value.var = "relaxed_thresh",
                  fill = 0)

# Tier 3 - 95% percentile genes that is correlated in 3/4 datasets
tier_3  <- df %>%
  dplyr::filter(relaxed_thresh == 1, num_datasets >= 3,
                percent_datasets >= 0.75) %>%
  dplyr::group_by(classifier_gene) %>%
  reshape2::acast(classifier_gene ~ gene, value.var = "relaxed_thresh",
                  fill = 0)

# Tier 4 - 95% percentile genes that is correlated in 2/4 datasets
tier_4  <- df %>%
  dplyr::filter(relaxed_thresh == 1, num_datasets >= 2,
                percent_datasets >= 0.50) %>%
  dplyr::group_by(classifier_gene) %>%
  reshape2::acast(classifier_gene ~ gene, value.var = "relaxed_thresh",
                  fill = 0)

# Create gmt geneset lists
tier_1a_gmt <- compile_gmt(tier_matrix = tier_1a, tier_name = "tier-1A")
tier_1b_gmt <- compile_gmt(tier_matrix = tier_1b, tier_name = "tier-1B")
tier_2a_gmt <- compile_gmt(tier_matrix = tier_2a, tier_name = "tier-2A")
tier_2b_gmt <- compile_gmt(tier_matrix = tier_2b, tier_name = "tier-2B")
tier_3_gmt <- compile_gmt(tier_matrix = tier_3, tier_name = "tier-3")
tier_4_gmt <- compile_gmt(tier_matrix = tier_4, tier_name = "tier-4")

big_gmt_list <- c(tier_1a_gmt, tier_1b_gmt, tier_2a_gmt, tier_2b_gmt,
                   tier_3_gmt, tier_4_gmt)

# Write contents of gmt list into file line by line
gmt_file <- file.path(results_dir, "correlated_hgsc_classifier_genes.gmt")
for (gmt in big_gmt_list) {
  write.table(t(gmt), file = gmt_file, sep = "\t", append = TRUE,
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}

# Perform an overrepresentation analysis on tier 1B genesets
ref_file <- file.path("7.Nanostring", "results", "background_genes.txt")
ref_genes <- readr::read_tsv(ref_file, col_names = FALSE)

geneset_dir <- file.path("7.Nanostring", "results", "gestalt")
if (!dir.exists(geneset_dir)) {
  dir.create(geneset_dir)
}

for (gmt_idx in 1:length(tier_1b_gmt)) {
  geneset_name <- names(tier_1b_gmt)[gmt_idx]
  geneset <- tier_1b_gmt[[geneset_name]]
  webgestalt_output <-
    WebGestaltR::WebGestaltR(enrichMethod = "ORA",
                             organism = "hsapiens",
                             interestGene = geneset,
                             interestGeneType = "genesymbol",
                             minNum = 4,
                             sigMethod = "top",
                             topThr = 15,
                             fdrMethod = "BH",
                             is.output = TRUE,
                             outputDirectory = geneset_dir,
                             referenceGene = as.vector(ref_genes$X1),
                             referenceGeneType = "genesymbol",
                             projectName = geneset_name)
}
