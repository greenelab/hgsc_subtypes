# Nanostring Classifier Genes
# Gregory Way 2018
# 7.Nanostring/scripts/A.get_correlation_output.R
#
# Load in classifier and gene expression matrices and obtain correlation
# between classifier genes and HGSC expressed genes across populations.
# The datasets we use include TCGA, Mayo Clinic, Yoshihara, and Tothill
#
# Output:
# 1) Correlation dataframe (long format)
# 2) Gene by Dataset missing info table

library(dplyr)

file <- file.path("data", "nanostring_classifier.tsv")
classifier_df <- readr::read_tsv(file)

data <- c("TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset", "GSE9891_eset")

# The script loads the ovarian cancer datasets
load.ovca.path <- file.path("..", "1.DataInclusion", "Scripts", "Functions",
                            "LoadOVCA_Data.R")
source(load.ovca.path)

# Use the LoadOVCA_Data function to read in the datasets subset by commongenes
ExpData <- LoadOVCA_Data(datasets = data, genelist_subset = "None")

all_cor <- list()
missing_info <- c()
comp_idx <- 1
# Loop through each unique gene in the classifier
for (gene in unique(classifier_df$gene)) {

  # Loop through each of the four datasets
  for (dataset in names(ExpData)) {
    exprs_df <- ExpData[[dataset]]

    # Make sure the gene exists in the dataset
    if (gene %in% rownames(exprs_df)) {
      # Subset gene expression vector of classifier gene
      gene_vector <- exprs_df[gene, ]

      # Obtain Pearson correlation of the gene vector and all other genes
      gene_cor <- apply(exprs_df, 1, function(x) {
        round(cor(x, gene_vector, method = "pearson"), 3)
      })

      # Convert the output to a dataframe that can be saved downstream
      gene_cor <- dplyr::as.tbl(data.frame(gene_cor)) %>%
        tibble::rownames_to_column(var = "gene")
      gene_cor$classifier_gene <- gene
      gene_cor$dataset <- dataset
      
      # Build results in iterative list
      all_cor[[comp_idx]] <- gene_cor
      comp_idx <- comp_idx + 1
    }

    # If the gene is not in the data frame store the missingness
    else {
      missing_gene_dataset <- c(gene, dataset)
      missing_info <- rbind(missing_info, missing_gene_dataset)
    }
  }
}

# Output results
# 1) Long correlation dataframe
output_df <- dplyr::bind_rows(all_cor)
file <- file.path("results", "classifier_correlation.tsv")
readr::write_tsv(output_df, file)

# 2) Missing gene by dataset information
missing_info <- as.data.frame(missing_info)
colnames(missing_info) <- c("gene", "dataset")
file <- file.path("results", "missing_gene_by_dataset.tsv")
readr::write_tsv(missing_info, file)
