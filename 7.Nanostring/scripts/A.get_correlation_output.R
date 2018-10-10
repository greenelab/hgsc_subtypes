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

file <- file.path("7.Nanostring", "data", "overallFreqs.csv")

# OPTION 1: Classifier genes
# Column rf stores the classifier genes - sort and take top 59
# top_n_genes <- 59
# classifier_df <- readr::read_csv(file) %>%
#   dplyr::arrange(desc(rfFreq)) %>%
#   dplyr::top_n(n = top_n_genes) %>%
#   dplyr::mutate(genes = toupper(genes))

# OPTION 2: 10 genes for remaining variability
# Select 10 genes that are thought to capture most remaining variability 
# that is not captured by the classifier genes
classifier_df <- readr::read_csv(file) %>%
  dplyr::filter(genes == "BOP1" 
                  |genes == "DNAI1"
                  |genes == "HSF1"
                  |genes == "LRRC50"
                  |genes == "MS4A3"
                  |genes == "NTN2L"
                  |genes == "SHARPIN"
                  |genes == "SLC12A3"
                  |genes == "SOX10"
                  |genes == "TSNAXIP1")
  
# OPTION 3: All 513 genes
#classifier_df <- readr::read_csv(file) %>%
#   dplyr::mutate(genes = toupper(genes))

# OPTION 4: 454 genes that are not part of the classifier
  # Column rf stores the classifier genes - sort and take bottom 454
  # bottom_n_genes <- 454
  # classifier_df <- readr::read_csv(file) %>%
  #   dplyr::arrange(desc(rfFreq)) %>%
  #   dplyr::bottom_n(n = bottom_n_genes) %>%
  #   dplyr::mutate(genes = toupper(genes))

data <- c("TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset", "GSE9891_eset")

# The script loads the ovarian cancer datasets
load.ovca.path <- file.path("1.DataInclusion", "Scripts", "Functions",
                            "LoadOVCA_Data.R")
source(load.ovca.path)

# Use the LoadOVCA_Data function to read in the datasets subset by commongenes
ExpData <- LoadOVCA_Data(datasets = data, genelist_subset = "None")

all_cor <- list()
missing_info <- c()
comp_idx <- 1
background_genes <- c()
# Loop through each unique gene in the classifier
for (gene in unique(classifier_df$genes)) {

  # Loop through each of the four datasets
  for (dataset in names(ExpData)) {
    exprs_df <- ExpData[[dataset]]
    
    # Track background genes for pathway analysis comparisons
    background_genes <- unique(c(background_genes, rownames(exprs_df)))

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
base_dir <- file.path("7.Nanostring", "results")
dir.create(base_dir)

fig_dir <- file.path("7.Nanostring", "figures")
dir.create(fig_dir)

# 1) Long correlation dataframe
output_df <- dplyr::bind_rows(all_cor)
file <- file.path(base_dir, "rf_classifier_correlation.tsv")
readr::write_tsv(output_df, file)

# 2) Missing gene by dataset information
rownames(missing_info) <- 1:nrow(missing_info)
missing_info <- as.data.frame(missing_info)
colnames(missing_info) <- c("gene", "dataset")
file <- file.path(base_dir, "rf_missing_gene_by_dataset.tsv")
readr::write_tsv(missing_info, file)

# 3) Background genes
background_df <- as.data.frame(background_genes[order(background_genes)])
file <- file.path(base_dir, "background_genes.txt")
readr::write_tsv(background_df, file, col_names = FALSE)
