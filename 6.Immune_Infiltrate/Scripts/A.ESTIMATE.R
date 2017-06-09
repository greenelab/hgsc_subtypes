# Timothy Chang 2017
# ESTIMATE Tumor Purity Analysis
# 6.Immune_Infiltrate/Scripts/A.ESTIMATE.R
#
# Analyzes tumor purity with the ESTIMATE algorithm and plots 
# Stromal/Immune infiltrate and tumor purity by subtype
#
# Usage: Run by selecting ovarian cancer datasets to be analyzed
#
#     Rscript 6.Immune_Infiltrate/Scripts/A.ESTIMATE.R
#
#     with the required argument:
#         --args        Datasets to be analyzed
#
# Output:
# Tables with Stromal, Immune and ESTIMATE score along with tumor purity
# Box plots by subtype for ESTIMATE score and tumor purity
# Scatter plots by subtype for comparing Stromal and Immune score

suppressMessages(library(checkpoint))
suppressMessages(checkpoint('2016-03-01', checkpointLocation = "."))

args <- commandArgs(trailingOnly = TRUE)

# args <- c("TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset", "GSE9891_eset")

############################################
# Load Libraries
############################################
library(curatedOvarianData)
library(estimate)
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

# Write data to disk for ESTIMATE functions
for (dataset in names(ExpData)) {
  fName <- paste0(dataset, "_processed.gct")
  fPath <- file.path("6.Immune_Infiltrate", "Tables", "ESTIMATE", fName)
  write.table(ExpData[[dataset]], fPath, sep = '\t', col.names = NA, 
              quote = FALSE)
}

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
# Compute ESTIMATE Scores
############################################
datasetPlatform <- c()

for (dataset in args) {
  processedFile <- paste0(dataset, "_processed.gct")
  filteredFile <- paste0(dataset, "_filtered.gct")
  estimateFile <- paste0(dataset, "_estimate.gct")

  # Filter for common genes in the ESTIMATE function
  inputFile <- file.path("6.Immune_Infiltrate", "Tables", "ESTIMATE", 
                         processedFile)
  outputFile <- file.path("6.Immune_Infiltrate", "Tables", "ESTIMATE", 
                          filteredFile)
  filterCommonGenes(input.f = inputFile,
                    output.f = outputFile,
                    id = "GeneSymbol")

  # Check if Mayo (Agilent)
  if (dataset == "mayo.eset") {
    datasetPlatform[dataset] <- "agilent"
  # Check if dataset is based on the Agilent platform
  } else if (grepl("Agilent", 
             experimentData(eval(as.name(dataset)))@other$platform_manufacturer,
             ignore.case = TRUE)) {
    datasetPlatform[dataset] <- "agilent"
  # Run remaining sets as Affymetrix
  } else {
    datasetPlatform[dataset] <- "affymetrix"
  }

  # Run ESTIMATE algorithm
  inputFile <- file.path("6.Immune_Infiltrate", "Tables", "ESTIMATE", 
                         filteredFile)
  outputFile <- file.path("6.Immune_Infiltrate", "Tables", "ESTIMATE", 
                          estimateFile)
  estimateScore(input.ds = inputFile,
                output.ds = outputFile,
                platform = datasetPlatform[dataset])
}

############################################
# Read ESTIMATE Scores
############################################
# This loop will pull all ESTIMATE files and load them into a list 
# according to the order of the datasets passed as arguments
ESTIMATEResultsList <- list()

estimatePath <- file.path("6.Immune_Infiltrate", "Tables", "ESTIMATE")
datasetMembers <- list.files(path = estimatePath)

for (dataset in args) {
  EstMembers <- datasetMembers[grep(dataset, datasetMembers)]
  EstMember <- EstMembers[grep("estimate", EstMembers)]

  # Read and clean up input table
  estimateFile <- file.path(estimatePath, EstMember)
  EstData <- read.table(file = estimateFile,
                        sep = "\t",
                        skip = 2,
                        header = TRUE,
                        stringsAsFactors = FALSE)
  rownames(EstData) <- EstData[, 2]
  EstData$NAME <- NULL
  EstData$Description <- NULL
  EstData <- t(EstData)
  EstData <- tibble::rownames_to_column(as.data.frame(EstData), "Sample")

  ESTIMATEResultsList[[dataset]] <- EstData
}

# Duplicate table for NMF
ESTIMATEResultsList_NMF <- ESTIMATEResultsList

############################################
# Merge ESTIMATE scores and clustering data
############################################
for (dataset in args) {
  ESTIMATEResultsList[[dataset]] <- 
    inner_join(ESTIMATEResultsList[[dataset]], 
               ClusterMembershipList[[dataset]],
               by = "Sample")
  ESTIMATEResultsList_NMF[[dataset]] <- 
    inner_join(ESTIMATEResultsList[[dataset]],
               ClusterMembershipList_NMF[[dataset]],
               by = "Sample")
}

############################################
# Analaysis of ESTIMATE and Tumor Purity
############################################
k_list <- c("K2", "K3", "K4")
boxplot_theme <- theme(panel.background = element_blank(),
                       panel.border = element_rect(color = "black", fill = NA),
                       axis.text.y = element_text(size = 14, angle = 90,
                                                  hjust = 0.4, 
                                                  colour = "black"),
                       axis.text.x = element_text(size = 14, colour = "black"))

# Analysis with k-means clustering
for (dataset in args) {
  ESTIMATEPlots <- list()
  PurityPlots <- list()

  data_iter <- 1

  # Make a glob for each clustering type
  for (k in k_list) {
    num_clus <- gsub("K", "", k)

    # Extract dataframe with scores and clustering
    est_rows <- c("Sample", "ESTIMATEScore", paste0("Cluster", k))
    est_df <- ESTIMATEResultsList[[dataset]][est_rows]
    colnames(est_df)[3] <- "Cluster"

    # Plot ESTIMATE Score by subtype
    e <- ggplot(data = est_df, 
                mapping = aes(x = Cluster, y = ESTIMATEScore, 
                              group = Cluster)) +
      geom_boxplot() +
      scale_x_discrete(limits = paste(1:num_clus)) + 
      labs(title = "", x = "", y = "") +
      boxplot_theme
    ESTIMATEPlots[[data_iter]] <- e

    # Analyze tumor purity for Affymetrix datasets
    if (datasetPlatform[dataset] == "affymetrix") {
      # Extract dataframe with purity and clustering
      pur_rows <- c("Sample", "TumorPurity", paste0("Cluster", k))
      pur_df <- ESTIMATEResultsList[[dataset]][pur_rows]
      colnames(pur_df)[3] <- "Cluster"

      # Plot Tumor Purity by subtype
      p <- ggplot(data = pur_df, 
                  mapping = aes(x = Cluster, y = TumorPurity, 
                                group = Cluster)) +
        geom_boxplot() +
        scale_x_discrete(limits = paste(1:num_clus)) + 
        labs(title = "", x = "", y = "") +
        boxplot_theme
      PurityPlots[[data_iter]] <- p
    }

    data_iter <- data_iter + 1
  }

  # Build a string to combine ESTIMATE plots
  est_plot_eval <- "full_grobs <- grid.arrange("
  for (i in 1:length(ESTIMATEPlots)) {
    est_plot_eval <- paste0(est_plot_eval, "ESTIMATEPlots[[", i,"]]", ", ")
  }
  est_plot_eval <- paste0(est_plot_eval, "ncol = length(ESTIMATEPlots)", 
                          ", nrow = 1)")

  # Write plot to PNG
  fName <- paste0(dataset, "_ESTIMATEScore_kmeans.png")
  fPath <- file.path("6.Immune_Infiltrate", "Figures", "ESTIMATE", fName)
  png(filename = fPath,
      width = 400 * length(ESTIMATEPlots),
      height = 400)
  eval(parse(text = est_plot_eval))
  dev.off()

  if (datasetPlatform[dataset] == "affymetrix") {
    # Build a string to evaluate Tumor Purity plots
    pur_plot_eval <- "full_grobs <- grid.arrange("
    for (i in 1:length(PurityPlots)) {
      pur_plot_eval <- paste0(pur_plot_eval, "PurityPlots[[", i, "]]", ", ")
    }
    pur_plot_eval <- paste0(pur_plot_eval, "ncol = length(PurityPlots)", 
                            ", nrow = 1)")

    # Write plot to PNG
    fName <- paste0(dataset, "_TumorPurity_kmeans.png")
    fPath <- file.path("6.Immune_Infiltrate", "Figures", "ESTIMATE", fName)
    png(filename = fPath,
        width = 400 * length(PurityPlots),
        height = 400)
    eval(parse(text = pur_plot_eval))
    dev.off()
  }
}

############################################
# Analaysis of Stromal and Immune Scores
############################################
for (dataset in args) {
  StromalImmunePlots <- list()

  data_iter <- 1

  for (k in k_list) {
    num_clus <- gsub("K", "", k)

    # Extract dataframe with scores and clustering
    svi_rows <- c("Sample", "StromalScore", "ImmuneScore", paste0("Cluster", k))
    svi_df <- ESTIMATEResultsList[[dataset]][svi_rows]
    colnames(svi_df)[4] <- "Cluster"

    # Plot scatter plot of Stromal vs. Immune Score by subtype
    s <- ggplot(data = svi_df,
                mapping = aes(x = StromalScore, y = ImmuneScore, 
                              color = factor(Cluster))) +
      geom_point(cex = 2.5) +
      scale_colour_manual(values = c("skyblue1", "tomato", "springgreen", 
                                     "violet")) +
      labs(title = "",
           color = "Subtype",
           x = "Stromal Score",
           y = "Immune Score") +
      theme(panel.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            axis.text.y = element_text(size = 16, angle = 90,
                                       hjust = 0.4, colour = "black"),
            axis.text.x = element_text(size = 16, colour = "black"),
            axis.title.x = element_text(size = 18, colour = "black"),
            axis.title.y = element_text(size = 18, angle = 90,
                                       hjust = 0.4, colour = "black"),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14))
    StromalImmunePlots[[data_iter]] <- s

    data_iter <- data_iter + 1
  }

  # Build a string to combine StromalImmune plots
  svi_plot_eval <- "full_grobs <- grid.arrange("
  for (i in 1:length(StromalImmunePlots)) {
    svi_plot_eval <- paste0(svi_plot_eval, "StromalImmunePlots[[", i,"]]", ", ")
  }
  svi_plot_eval <- paste0(svi_plot_eval, "ncol = length(StromalImmunePlots)", 
                          ", nrow = 1)")

  # Write plot to PNG
  fName <- paste0(dataset, "_StromalvImmune_kmeans.png")
  fPath <- file.path("6.Immune_Infiltrate", "Figures", "ESTIMATE", fName)
  png(filename = fPath,
      width = 800 * length(StromalImmunePlots),
      height = 650)
  eval(parse(text = svi_plot_eval))
  dev.off()
}
