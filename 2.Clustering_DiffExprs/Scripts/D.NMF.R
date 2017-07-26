############################################
# Cross-population analysis of high-grade serous ovarian cancer does
# not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will input a series of datasets and perform Non-negative Matrix
# Factorization (NMF)

args <- commandArgs(trailingOnly = TRUE)
# args <- c(2, 4, 100, 123, "TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset",
#          "GSE9891_eset", "GSE26712_eset", "aaces.eset")
################################  
# Load Libraries
################################
library(curatedOvarianData)
library(cluster)
library(NMF)

# Load important kmeans and SAM functions
kmeansfxnpath <- file.path("2.Clustering_DiffExprs", "Scripts", "Functions", 
                           "kmeans_SAM_functions.R")
source(kmeansfxnpath)

# Loads curatedOvarianData
loadovcapath <- file.path("1.DataInclusion", "Scripts", "Functions",
                          "LoadOVCA_Data.R")
source(loadovcapath)

################################
# Constants
################################
rank1 = as.numeric(paste(args[1]))
rank2 = as.numeric(paste(args[2]))
nruns = as.numeric(paste(args[3]))
set.seed(args[4])

############################################
# Load Data
############################################
# Separate the eset arguments from the rest of the commandArgs
argsCurated <- args[grep("eset", args)]

# Add Mayo to argsCurated
if("Mayo" %in% args) {
  argsCurated = c(argsCurated[1], "Mayo", argsCurated[2:length(argsCurated)])
}

# Use the LoadOVCA_Data function to read in the datasets subset by MAD genes
ExpData <- LoadOVCA_Data(datasets = argsCurated, genelist_subset = "madgenes")

############################################
# Run NMF
############################################
nmfClusters <- list()
for (dataset in 1:length(ExpData)) {
  
  # Load corresponding k-means cluster membership file
  membfile <- file.path("2.Clustering_DiffExprs",
                        "Tables", "ClusterMembership", "kmeans",
                        paste0("KMembership_", names(ExpData)[[dataset]],
                               ".csv"))
  ClusterAssign <- read.csv(file = membfile, row.names = 1)
  
  # Initialize a blank matrix with the same number of rows as samples and 
  # kmin to kmax columns
  newMatrix <- matrix(NA, nrow = nrow(ClusterAssign),
                      ncol = length(rank1:rank2))
  rownames(newMatrix) <- rownames(ClusterAssign)
  
  # For each centroid assignment, perform NMF and assign cluster 
  # membership files to newMatrix
  for (centroid in 1:length(rank1:rank2)) {
    fname <- paste0("NMF_Clusters_", names(ExpData)[dataset], "_K",
                   centroid + 1)
    
    # Run custom NMF function with 100 runs for k = 2, 3, and 4
    newMatrix[ , centroid] <- runNMF(ExpData[[dataset]],
                                     k = centroid + 1, nrun = nruns, 
                                     KClusterAssign = ClusterAssign,
                                     fname = fname)
  }
  
  # Assign results to list
  nmfClusters[[dataset]] <- newMatrix
  
  # Run custom NMF function with 10 runs to output cophenetic coefficient
  coph_coeff <- runNMF(ExpData[[dataset]], coph = T)
  
  # Base name to save figures
  sup_fname <- paste("NMF_Clusters_", names(ExpData)[dataset], sep = "")
  
  # Write the consensus mapping as a figure
  fname <- file.path("2.Clustering_DiffExprs", "Figures", "nmf",
                     "CopheneticMaps", paste0(sup_fname, "_k2-6.png"))
  png(fname, width = 950, height = 550)
  
  # Plot consensus maps
  consensusmap(coph_coeff$consensus[1:5], labCol = NA, labRow = NA,
               main = "", fontsize = 12)
  
  # Close the connection
  dev.off()
  
  # Write the Cophenetic Coefficient as a figure
  fname <- file.path("2.Clustering_DiffExprs", "Figures", "nmf",
                     "CopheneticMaps", paste0(sup_fname, "coph_coeff", ".png"))
  
  png(fname, width = 270, height = 230)
  
  # Plot cophenetic coefficient
  par(mar = c(4.5, 4.5, 1.5, 1))
  plot(coph_coeff$measures$cophenetic, xaxt = "n", cex.axis = 1.5,
       cex.lab = 1.5, xlab = 'k', ylab = 'Cophenetic Correlation')
  axis(1, at=1:7, labels=2:8, cex.axis = 1.5)
  lines(coph_coeff$measures$cophenetic, lwd = 2)
  points(coph_coeff$measures$cophenetic, col = 'black', pch = 19, cex = 1.2)
  
  # Close the connection
  dev.off()
}

############################################
# Write NMF clusters
############################################
for (dataset in 1:length(nmfClusters)) {
  filepath <- file.path("2.Clustering_DiffExprs", "Tables", "ClusterMembership",
                        "nmf", paste0(names(ExpData[dataset]), "_nmf.csv"))
  write.table(as.data.frame(nmfClusters[[dataset]]), 
              file = filepath, sep = ",", row.names = T, col.names = NA)
}
