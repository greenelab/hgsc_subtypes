# Timothy Chang 2017
# ssGSEA Analysis
# 6.Immune_Infiltrate/Scripts/B.ssGSEA.R
#
# Analyzes gene set enrichment of 22 leukocyte signatures (LM22) with 
# ssGSEA and plots enrichment by cancer subtype
#
# LM22 obtained from CIBERSORT (Newman et al., 2015)
# [https://doi.org/10.1038/nmeth.3337]
#
# Usage: Run by selecting ovarian cancer datasets to be analyzed
#
#     R --no-save --args $DATASETS < 6.Immune_Infiltrate/Scripts/B.ssGSEA.R
#
#     with the argument:
#         $DATASETS     a comma separated string of dataset names
#
# Output:
# Tables with ssGSEA scores for each cell class
# Box plots by subtype for ssGSEA scores

suppressMessages(library(checkpoint))
suppressMessages(checkpoint('2016-03-01', checkpointLocation = "."))

args <- commandArgs(trailingOnly = TRUE)

# args <- c("TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset", "GSE9891_eset")

############################################
# Load Libraries
############################################
library(curatedOvarianData)
library(GSVA)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)

# Loads curatedOvarianData
loadovcaPath <- file.path("1.DataInclusion", "Scripts", "Functions",
                          "LoadOVCA_Data.R")
source(loadovcaPath)

############################################
# Data in curatedOvarianData
############################################
# Get data from the curatedOvarianData package  
detailedData <- data(package = "curatedOvarianData")[3]
# Get the eset IDs
detailedData.names <- detailedData$results[ ,3]
# Determine the datasets in curatedOvarianData
inCuratedOvarianData <- args %in% detailedData.names

# Here are the args that where curated
argsCurated <- args[inCuratedOvarianData]

# Add Mayo to argsCurated
if ("mayo.eset" %in% args) {
  argsCurated = c(argsCurated[1], "mayo.eset", 
                  argsCurated[2:length(argsCurated)])
}

############################################
# Load Data
############################################
# Use the LoadOVCA_Data function to read in the datasets
ExpData <- LoadOVCA_Data(datasets = argsCurated, 
                         genelist_subset = "commongenes",
                         zscore = TRUE)

############################################
# Get Cluster Membership Info
############################################
# This loop will pull all cluster membership files and load them into a list 
# according to the order of the datasets passed as arguments
ClusterMembershipList <- list()
ClusterMembershipList_NMF <- list()

kmeansPath <- file.path("2.Clustering_DiffExprs", "Tables", 
                        "ClusterMembership", "kmeans")
datasetMembers <- list.files(path = kmeansPath)

nmfPath <- file.path("2.Clustering_DiffExprs", "Tables", 
                     "ClusterMembership", "nmf")
datasetMembers_NMF <- list.files(path = nmfPath) 

for (dataset in args) {  
  # kmeans cluster membership
  ClusMember <- datasetMembers[grep(dataset, datasetMembers)]
  kmeansFile <- file.path(kmeansPath, ClusMember)
  ClusData <- read.csv(kmeansFile, stringsAsFactors = FALSE)
  names(ClusData)[1] <- "Sample"
  ClusterMembershipList[[dataset]] <- ClusData
  
  # NMF cluster membership
  ClusMember <- datasetMembers_NMF[grep(dataset, datasetMembers_NMF)]
  ClusMember <- ClusMember[grepl("mapped", ClusMember)]
  nmfFile <- file.path(nmfPath, ClusMember)
  ClusData <- read.csv(nmfFile, stringsAsFactors = FALSE)
  names(ClusData)[1] <- "Sample"
  ClusterMembershipList_NMF[[dataset]] <- ClusData
}

############################################
# Prepare LM22 Data
############################################
# Load LM22 gene signatures by cell type
lm22Path <- file.path("6.Immune_Infiltrate", "Gene_signatures", 
                      "lm22_genes.tsv")
lm22_genes <- read.table(lm22Path, sep = "\t", header = TRUE, row.names = 1,
                         stringsAsFactors = FALSE)
lm22_geneset <- lapply(lm22_genes, function(x) rownames(lm22_genes)[x == 1])

# Load LM22 cell class list
lm22_agg <- list()
lm22ClassPath <- file.path("6.Immune_Infiltrate", "Gene_signatures", 
                           "lm22_classes.csv")
lm22_agg_read <- read.csv(lm22ClassPath, header = TRUE, 
                          stringsAsFactors = FALSE, fill = TRUE)

# Create list of cell classes and included cell types
data_iter <- 1
for (cell_class in lm22_agg_read) {
  types <- cell_class[!is.na(cell_class)]
	lm22_agg[[data_iter]] <- names(lm22_geneset)[types]

  data_iter <- data_iter + 1
}
names(lm22_agg) <- colnames(lm22_agg_read)

# Create gene signature list for each cell class
lm22_agg_geneset <- list()
for (cell_class in 1:length(lm22_agg)) {
  types <- lm22_agg[[cell_class]]
  unique_genes <- unique(unlist(lm22_geneset[types], use.names = FALSE))
  lm22_agg_geneset[[cell_class]] <- unique_genes
}
names(lm22_agg_geneset) <- names(lm22_agg)

############################################
# Compute ssGSEA Scores
############################################
ssGSEA_Result <- list()
ssGSEA_Result_class <- list()

for (dataset in args) {
  # Load gene expression data
  gene_exp <- ExpData[[dataset]]
  gene_input <- rownames(gene_exp)

  ## Run ssGSEA function by cell type
  ssgsea_type <- GSVA::gsva(expr = gene_exp,
                            gset.idx.list = lm22_geneset,
                            method = "ssgsea",
                            verbose = FALSE)
  ssgsea_type <- t(ssgsea_type)
  ssgsea_type <- tibble::rownames_to_column(as.data.frame(ssgsea_type), 
                                             "Sample")
  ssGSEA_Result[[dataset]] <- ssgsea_type

  # Write results to disk
  fName <- paste0(dataset, "_ssGSEA.csv")
  fPath <- file.path("6.Immune_Infiltrate", "Tables", "ssGSEA", fName)
  write.csv(ssgsea_type, fPath, quote = FALSE)

  ## Run ssGSEA function by cell class
  ssgsea_class <- GSVA::gsva(expr = gene_exp, 
                             gset.idx.list = lm22_agg_geneset,
                             method = "ssgsea",
                             verbose = FALSE)
  ssgsea_class <- t(ssgsea_class)
  ssgsea_class <- tibble::rownames_to_column(as.data.frame(ssgsea_class), 
                                             "Sample")
  ssGSEA_Result_class[[dataset]] <- ssgsea_class

  # Write results to disk
  fName <- paste0(dataset, "_ssGSEA_Class.csv")
  fPath <- file.path("6.Immune_Infiltrate", "Tables", "ssGSEA", fName)
  write.csv(ssgsea_class, fPath, quote = FALSE)
}

############################################
# Merge ssGSEA Scores and Clustering Data
############################################
for (dataset in args) {
  ssGSEA_Result[[dataset]] <-
    inner_join(ssGSEA_Result[[dataset]],
               ClusterMembershipList[[dataset]],
               by = "Sample")

  ssGSEA_Result_class[[dataset]] <-
    inner_join(ssGSEA_Result_class[[dataset]],
               ClusterMembershipList[[dataset]],
               by = "Sample")
}

############################################
# Analysis of ssGSEA Scores
############################################
k_list <- c("K2", "K3", "K4")
boxplot_theme <- theme(panel.background = element_blank(),
                       panel.border = element_rect(color = "black", fill = NA),
                       axis.text.y = element_text(size = 12, angle = 90,
                                                  hjust = 0.4, 
                                                  colour = "black"),
                       axis.text.x = element_text(size = 12, colour = "black"),
                       axis.title.y = element_text(size = 12, colour = "black",
                                                   face = "bold"))

# Analysis with k-means clustering
for (dataset in args) {
  # Make a glob for each clustering type
  for (k in k_list) {
    ssGSEAPlots <- list()
    ssGSEAPlotsClass <- list()

    data_iter <- 1
    data_iter_class <- 1

    num_clus <- gsub("K", "", k)

    ## Analysis for cell types
    for (type in names(lm22_geneset)) {
      plot_cols <- c("Sample", type, paste0("Cluster", k))
      plot_df <- ssGSEA_Result[[dataset]][plot_cols]
      colnames(plot_df)[2] <- "Score"
      colnames(plot_df)[3] <- "Cluster"

      # Plot ssGSEA score by subtype
      c <- ggplot(data = plot_df,
                  mapping = aes(x = Cluster, y = Score, group = Cluster)) +
        geom_boxplot() +
        scale_x_discrete(limits = paste(1:num_clus)) +
        labs(title = "", x = "", y = type) +
        boxplot_theme
      ssGSEAPlots[[data_iter]] <- c

      data_iter <- data_iter + 1
    }

    # Build string to combine ssGSEA plots
    plot_eval <- "full_grobs <- grid.arrange("
    for (i in 1:length(ssGSEAPlots)) {
      plot_eval <- paste0(plot_eval, "ssGSEAPlots[[", i, "]]", ", ")
    }
    plot_eval <- paste0(plot_eval, "ncol = 4", 
                        ", nrow = ceiling(length(ssGSEAPlots) / 4))")

    # Write plot to PNG
    fName <- paste0(dataset, "_", k, "_ssGSEA.png")
    fPath <- file.path("6.Immune_Infiltrate", "Figures", "ssGSEA", 
                       "Cell_type", fName)
    png(filename = fPath,
        width = 300 * 4,
        height = 300 * ceiling(length(ssGSEAPlots) / 4))
    eval(parse(text = plot_eval))
    dev.off()

    ## Analysis for cell classes
    for (class in names(lm22_agg_geneset)) {
      plot_cols <- c("Sample", class, paste0("Cluster", k))
      plot_df <- ssGSEA_Result_class[[dataset]][plot_cols]
      colnames(plot_df)[2] <- "Score"
      colnames(plot_df)[3] <- "Cluster"

      # Plot ssGSEA score by subtype
      g <- ggplot(data = plot_df,
                  mapping = aes(x = Cluster, y = Score, group = Cluster)) +
        geom_boxplot() +
        scale_x_discrete(limits = paste(1:num_clus)) +
        labs(title = "", x = "", y = class) +
        boxplot_theme
      ssGSEAPlotsClass[[data_iter_class]] <- g

      data_iter_class <- data_iter_class + 1
    }

    # Build string to combine ssGSEA plots
    plot_eval <- "full_grobs <- grid.arrange("
    for (i in 1:length(ssGSEAPlotsClass)) {
      plot_eval <- paste0(plot_eval, "ssGSEAPlotsClass[[", i, "]]", ", ")
    }
    plot_eval <- paste0(plot_eval, "ncol = 3", 
                        ", nrow = ceiling(length(ssGSEAPlotsClass) / 3))")

    # Write plot to PNG
    fName <- paste0(dataset, "_", k, "_ssGSEA_Class.png")
    fPath <- file.path("6.Immune_Infiltrate", "Figures", "ssGSEA", 
                       "Cell_class", fName)
    png(filename = fPath,
        width = 300 * 3,
        height = 300 * ceiling(length(ssGSEAPlotsClass) / 3))
    eval(parse(text = plot_eval))
    dev.off()
  }
}
