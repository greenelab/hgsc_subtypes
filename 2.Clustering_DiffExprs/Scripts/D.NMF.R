############################################
# High-grade serous ovarian cancer subtypes are similar across populations
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will input a series of datasets and perform Non-negative Matrix Factorization (NMF)

args <- commandArgs(trailingOnly=TRUE)
# args <- c(3, 4, 100, 123, "TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset", "GSE26712_eset")
################################
# Load Libraries
################################
library(curatedOvarianData)
library(cluster)
library(NMF)

# Load important kmeans and SAM functions
source("2.Clustering_DiffExprs/Scripts/Functions/kmeans_SAM_functions.R")

# Loads curatedOvarianData
source("1.DataInclusion/Scripts/Functions/LoadOVCA_Data.R")

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
  membfile <- paste("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/KMembership_", 
                    names(ExpData)[[dataset]], ".csv", sep = "")
  ClusterAssign <- read.csv(file = membfile, row.names = 1)
  
  # Initialize a blank matrix with the same number of rows as samples and kmin to kmax columns
  newMatrix <- matrix(NA, nrow = nrow(ClusterAssign), ncol = length(rank1:rank2))
  rownames(newMatrix) <- rownames(ClusterAssign)
  
  # For each centroid assignment, perform NMF and assign cluster membership files to newMatrix
  for (centroid in rank1:rank2) {
    fname <- paste("NMF_Clusters_", names(ExpData)[dataset], "_K", centroid, sep = "")
    newMatrix[ ,(centroid-2)] <- runNMF(ExpData[[dataset]], k = centroid, nrun = nruns, 
                                        KClusterAssign = ClusterAssign, fname = fname)
  } 
  
  # Assign results to list
  nmfClusters[[dataset]] <- newMatrix
}

############################################
# Write NMF clusters
############################################
for (dataset in 1:length(nmfClusters)) {
  write.table(as.data.frame(nmfClusters[[dataset]]), 
              file = paste("2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/", names(ExpData)[dataset], 
                           "_nmf.csv", sep = ""), sep = ",", row.names = T, col.names = NA)
}
