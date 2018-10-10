# Nanostring Classifier Genes
# Gregory Way 2018
# 7.Nanostring/scripts/C.threshold_venns.R
#
# Obtain highly correlated genes (by 95% or 99% quantile) for each dataset
# across all 29 genes. The four datasets are TCGA, Mayo, Yoshihara, and Tothill
#
# Output:
# Two Venn Diagrams for Every Classifier Gene
# The two venn diagrams will display overlap across the four datasets for the
# two threshold levels across all genes

library(dplyr)
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

output_all_venns <- function(df, output_dir, threshold) {
  # Get the overlap of genes correlating with classifier genes across datasets
  #
  # Arguments:
  # df - Thresholded list of dataset gene correlations to classifier genes
  # output_dir - the folder to save the results
  # threshold - the percentile of filtered gene correlations
  #
  # Output:
  # Venn Diagram figure of gene overlap

  # Obtain list of lists that stores sets of correlated genes by classifier
  # genes by dataset - this will be used for the Venn Diagram plotting
  set_list <- list()
  range_list <- list()

  # Loop over each classifier gene
  for (class_gene in unique(df$classifier_gene)) {
    
    # Subset the dataframe to each unique classifier gene
    gene_subset <- df %>% dplyr::filter(classifier_gene == class_gene)
    
    # Initialize an empty list within the list named by the classifier gene
    set_list[[class_gene]] <- list()
    range_list[[class_gene]] <- list()

    # Loop over each dataset
    for (data in unique(gene_subset$dataset)) {
      
      # Subset the gene subset data frame to each unique dataset
      dataset_gene_subset <- gene_subset %>% dplyr::filter(dataset == data)

      # Store a set of correlated genes per dataset per classifier gene
      set_list[[class_gene]][[data]] <- unique(dataset_gene_subset$gene)

      # Store the range of the correlations
      cor_range <- paste(round(range(dataset_gene_subset$gene_cor_abs), 2),
                         collapse = " - ")
      range_list[[class_gene]][[data]] <- cor_range
    }
  }

  # Create a Venn Diagram for each gene - the circles will represent datasets
  for (gene_name in names(set_list)) {

    # Obtain the list for the current classifier gene
    current_set_list <- set_list[[gene_name]]
    current_range_list <- range_list[[gene_name]]

    # Set the file name to output
    file <- file.path(output_dir, paste0("venn_", gene_name, ".png"))
    
    # Which datasets are represented?
    # Not all datasets have all classifier genes measured
    datasets <- names(current_set_list)

    # Set labels in Venn Diagram
    yos <- paste0("Yoshihara\n", current_range_list[["GSE32062.GPL6480_eset"]])
    tot <- paste0("Tothill\n", current_range_list[["GSE9891_eset"]])
    tcga <- paste0("TCGA\n", current_range_list[["TCGA_eset"]])
    mayo <- paste0("Mayo\n", current_range_list[["mayo.eset"]])

    # Recode dataset eset nomenclature and fill color
    category_names <- dplyr::recode(datasets,
                                    "GSE32062.GPL6480_eset" = yos,
                                    "GSE9891_eset" = tot,
                                    "TCGA_eset" = tcga,
                                    "mayo.eset" = mayo)
    fill_colors <- dplyr::recode(datasets,
                                 "GSE32062.GPL6480_eset" = "blue",
                                 "GSE9891_eset" = "purple",
                                 "TCGA_eset" = "red",
                                 "mayo.eset" = "green")
    
    # Output the Venn Diagrams
    VennDiagram::venn.diagram(
      current_set_list,
      main = paste0(gene_name, "\n", threshold, "% Threshold"),
      main.cex = 1.5,
      main.pos = c(0.5, 1.03),
      category.names = category_names,
      filename = file,
      output = FALSE,
      imagetype = "png",
      fill = fill_colors,
      lwd = 2)
  }
}

# Load thresholded dataframes
base_file <- file.path("7.Nanostring", "results")
fig_base <- file.path("7.Nanostring", "figures")

file <- file.path(base_file, "rf_high_99_quantile_correlated_genes.tsv")
high_df <- readr::read_tsv(file)

file <- file.path(base_file, "rf_relaxed_95_quantile_correlated_genes.tsv")
relaxed_df <- readr::read_tsv(file)

# Output all Venn Diagrams
high_out <- file.path(fig_base, "rf_venns", "high_threshold")
dir.create(high_out, recursive = TRUE)
output_all_venns(high_df, output_dir = high_out, threshold = "99")

relaxed_out <- file.path(fig_base, "rf_venns", "relaxed_threshold")
dir.create(relaxed_out, recursive = TRUE)
output_all_venns(relaxed_df, output_dir = relaxed_out, threshold = "95")
