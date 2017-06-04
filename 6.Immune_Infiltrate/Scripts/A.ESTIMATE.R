############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will perform ESTIMATE analysis on the datasets and output infiltration
# results for all input datasets

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

# The script loads the ovarian cancer datasets
source("1.DataInclusion/Scripts/Functions/LoadOVCA_Data.R")

############################################
# Data in curatedOvarianData
############################################
# Get data from the curatedOvarianData package  
detailedData <- data(package="curatedOvarianData")[3]
# Get the eset IDs
detailedData.names <- detailedData$results[ ,3]
# Determine the datasets in curatedOvarianData
inCuratedOvarianData <- args %in% detailedData.names

# Here are the args that where curated
argsCurated <- args[inCuratedOvarianData]

# Add Mayo to argsCurated
if("mayo.eset" %in% args) {
  argsCurated = c(argsCurated[1], "mayo.eset", argsCurated[2:length(argsCurated)])
}

############################################
# Load Data
############################################
# Use the LoadOVCA_Data function to read in the datasets subset by commongenes
ExpData <- LoadOVCA_Data(datasets = argsCurated, genelist_subset = "commongenes")

############################################
# Get Cluster Membership Info
############################################
# This loop will pull all cluster membership files and load them into a list 
# according to the order of the datasets passed as arguments
ClusterMembershipList <- list()
ClusterMembershipList_NMF <- list()
datasetMembers <- list.files(path = "2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/")
datasetMembers_NMF <- list.files(path = "2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/") 
for (dataset in 1:length(args)) {
  holder <- args[dataset]
  ClusMember <- datasetMembers[grep(holder, datasetMembers)]
  
  ClusData <- read.csv(paste("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/", 
                             ClusMember, sep = ""), row.names = 1)
  
  ClusterMembershipList[[dataset]] <- ClusData
  names(ClusterMembershipList)[dataset] <- args[dataset]
  
  # NMF cluster membership
  ClusMember <- datasetMembers_NMF[grep(holder, datasetMembers_NMF)]
  ClusMember <- ClusMember[grepl("mapped", ClusMember)]
  ClusData <- read.csv(paste("2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/", 
                             ClusMember, sep = ""), row.names = 1)
  
  ClusterMembershipList_NMF[[dataset]] <- ClusData
  names(ClusterMembershipList_NMF)[dataset] <- args[dataset]
}

############################################
# Compute z-Score
############################################
# Function to compute z scores
z_score_func <- function(x) {
  z <- (x - mean(x)) / sd(x)
  return(z)
}

# Compute z-score by gene
for (dataset in 1:length(ExpData)) {
  for (gene in 1:dim(ExpData[[dataset]])[2]) {
    ExpData[[dataset]][, gene] <- z_score_func(ExpData[[dataset]][, gene])
  }

  # Write dataset to disk
  fName <- paste0(names(ExpData)[dataset], '_processed.gct')
  write.table(ExpData[[dataset]], 
              file = file.path("6.Immune_Infiltrate", "Tables", "ESTIMATE", fName),
              sep = "\t",
              col.names = NA,
              quote = FALSE)
}

############################################
# Compute ESTIMATE Scores
############################################
isAffymetrix <- c()

for (dataset in 1:length(ExpData)) {
  fName <- paste0(names(ExpData)[dataset], '_processed.gct')
  fName1 <- paste0(names(ExpData)[dataset], '_filtered.gct')
  fName2 <- paste0(names(ExpData)[dataset], '_estimate.gct')

  # Filter for common genes in the ESTIMATE function
  filterCommonGenes(input.f = file.path("6.Immune_Infiltrate", "Tables",
                                        "ESTIMATE", fName),
                    output.f = file.path("6.Immune_Infiltrate", "Tables",
                                         "ESTIMATE", fName1),
                    id = "GeneSymbol")

  # Check if Mayo (Agilent)
  if (grepl("mayo", names(ExpData)[dataset], ignore.case = TRUE)) {
    estimateScore(input.ds = file.path("6.Immune_Infiltrate", "Tables",
                                       "ESTIMATE", fName1),
                  output.ds = file.path("6.Immune_Infiltrate", "Tables",
                                        "ESTIMATE", fName2),
                  platform = "agilent")
    isAffymetrix[dataset] <- FALSE
  # Run Agilent-based datasets
  } else if (grepl("Agilent", 
             experimentData(eval(as.name(names(ExpData)[dataset])))@other$platform_manufacturer,
             ignore.case = TRUE)) {
    estimateScore(input.ds = file.path("6.Immune_Infiltrate", "Tables",
                                       "ESTIMATE", fName1),
                  output.ds = file.path("6.Immune_Infiltrate", "Tables",
                                        "ESTIMATE", fName2),
                  platform = "agilent")
    isAffymetrix[dataset] <- FALSE
  # Run Affymetrix-based datasets
  } else {
    estimateScore(input.ds = file.path("6.Immune_Infiltrate", "Tables",
                                       "ESTIMATE", fName1 ),
                  output.ds = file.path("6.Immune_Infiltrate", "Tables",
                                        "ESTIMATE", fName2),
                  platform = "affymetrix")
    isAffymetrix[dataset] <- TRUE
  }
}

############################################
# Read ESTIMATE Scores
############################################
# This loop will pull all ESTIMATE files and load them into a list 
# according to the order of the datasets passed as arguments
ESTIMATEResultsList <- list()
datasetMembers <- list.files(path = "6.Immune_Infiltrate/Tables/ESTIMATE/")
for (dataset in 1:length(args)) {
  holder <- args[dataset]
  EstMembers <- datasetMembers[grep(holder, datasetMembers)]
  EstMember <- EstMembers[grep("estimate", EstMembers)]

  # Read and clean up input table
  EstData <- read.table(file = paste0("6.Immune_Infiltrate/Tables/ESTIMATE/", 
                                      EstMember),
                        sep = "\t",
                        skip = 2,
                        header = T)
  rownames(EstData) <- EstData[, 2]
  EstData$NAME <- NULL
  EstData$Description <- NULL
  EstData <- t(EstData)

  ESTIMATEResultsList[[dataset]] <- EstData
  names(ESTIMATEResultsList)[dataset] <- args[dataset]
}

# Duplicate table for NMF
ESTIMATEResultsList_NMF <- ESTIMATEResultsList

############################################
# Merge ESTIMATE scores and clustering data
############################################
for (dataset in 1:length(args)) {
  ESTIMATEResultsList[[dataset]] <- merge(x = ESTIMATEResultsList[[dataset]], 
                                          y = ClusterMembershipList[[dataset]], 
                                          by = "row.names")
  ESTIMATEResultsList_NMF[[dataset]] <- merge(x = ESTIMATEResultsList[[dataset]],
                                              y = ClusterMembershipList_NMF[[dataset]],
                                              by = "row.names")
}

############################################
# Analaysis of ESTIMATE and Tumor Purity
############################################
k_list <- c("K2", "K3", "K4")

# Analysis with k-means clustering
for (dataset in 1:length(args)) {
  ESTIMATEPlots <- list()
  PurityPlots <- list()

  data_iter <- 1

  # Make a glob for each clustering type
  for (k in k_list) {
    num_clus <- gsub("K", "", k)

    # Extract dataframe with scores and clustering
    est_df <- ESTIMATEResultsList[[dataset]][c("Row.names", "ESTIMATEScore",
                                               paste0("Cluster", k))]
    colnames(est_df)[3] <- "Cluster"

    # Plot ESTIMATE Score by subtype
    e <- ggplot(data = est_df, 
                mapping = aes(x = Cluster, y = ESTIMATEScore, group = Cluster)) +
      geom_boxplot() +
      scale_x_discrete(limits = paste(1:num_clus)) + 
      labs(title = "", x = "", y = "") +
      theme(panel.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            axis.text.y = element_text(size = 14, angle = 90,
                                       hjust = 0.4, colour = "black"),
            axis.text.x = element_text(size = 14, colour = "black"))
    ESTIMATEPlots[[data_iter]] <- e

    # Analyze tumor purity for Affymetrix datasets
    if (isAffymetrix[dataset]) {
      # Extract dataframe with purity and clustering
      pur_df <- ESTIMATEResultsList[[dataset]][c("Row.names", "TumorPurity",
                                                 paste0("Cluster", k))]
      colnames(pur_df)[3] <- "Cluster"

      # Plot Tumor Purity by subtype
      p <- ggplot(data = pur_df, 
                  mapping = aes(x = Cluster, y = TumorPurity, group = Cluster)) +
        geom_boxplot() +
        scale_x_discrete(limits = paste(1:num_clus)) + 
        labs(title = "", x = "", y = "") +
        theme(panel.background = element_blank(),
              panel.border = element_rect(color = "black", fill = NA),
              axis.text.y = element_text(size = 14, angle = 90,
                                        hjust = 0.4, colour = "black"),
              axis.text.x = element_text(size = 14, colour = "black"))
      PurityPlots[[data_iter]] <- p
    }

    data_iter <- data_iter + 1
  }

  # Build a string to combine ESTIMATE plots
  est_plot_eval <- "full_grobs <- grid.arrange("
  for (i in 1:length(ESTIMATEPlots)) {
    est_plot_eval <- paste0(est_plot_eval, "ESTIMATEPlots[[", i,"]]", ", ")
  }
  est_plot_eval <- paste0(est_plot_eval, "ncol = length(ESTIMATEPlots)", ", nrow = 1)")

  # Write plot to PNG
  fName <- paste0(args[dataset], "_ESTIMATEScore_kmeans.png")
  png(filename = file.path("6.Immune_Infiltrate", "Figures", "ESTIMATE", fName),
      width = 400 * length(ESTIMATEPlots),
      height = 400)
  eval(parse(text = est_plot_eval))
  dev.off()

  if (isAffymetrix[dataset]) {
    # Build a string to evaluate Tumor Purity plots
    pur_plot_eval <- "full_grobs <- grid.arrange("
    for (i in 1:length(PurityPlots)) {
      pur_plot_eval <- paste0(pur_plot_eval, "PurityPlots[[", i, "]]", ", ")
    }
    pur_plot_eval <- paste0(pur_plot_eval, "ncol = length(PurityPlots)", ", nrow = 1)")

    # Write plot to PNG
    fName <- paste0(args[dataset], "_TumorPurity_kmeans.png")
    png(filename = file.path("6.Immune_Infiltrate", "Figures", "ESTIMATE", fName),
        width = 400 * length(PurityPlots),
        height = 400)
    eval(parse(text = pur_plot_eval))
    dev.off()
  }
}

############################################
# Analaysis of Stromal and Immune Scores
############################################
for (dataset in 1:length(args)) {
  StromalImmunePlots <- list()

  data_iter <- 1

  for (k in k_list) {
    num_clus <- gsub("K", "", k)

    # Extract dataframe with scores and clustering
    svi_df <- ESTIMATEResultsList[[dataset]][c("Row.names", "StromalScore",
                                               "ImmuneScore", paste0("Cluster", k))]
    colnames(svi_df)[4] <- "Cluster"

    # Plot scatter plot of Stromal vs. Immune Score by subtype
    s <- ggplot(data = svi_df,
                mapping = aes(x = StromalScore, y = ImmuneScore, color = factor(Cluster))) +
      geom_point(cex = 2.5) +
      scale_color_brewer(palette = "Set1") +
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
  svi_plot_eval <- paste0(svi_plot_eval, "ncol = length(StromalImmunePlots)", ", nrow = 1)")

  # Write plot to PNG
  fName <- paste0(args[dataset], "_StromalvImmune_kmeans.png")
  png(filename = file.path("6.Immune_Infiltrate", "Figures", "ESTIMATE", fName),
      width = 800 * length(StromalImmunePlots),
      height = 650)
  eval(parse(text = svi_plot_eval))
  dev.off()
}
