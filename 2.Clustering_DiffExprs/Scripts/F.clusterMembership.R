############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will output a data frame of cluster membership for all clustering events

suppressMessages(library(checkpoint))
suppressMessages(checkpoint('2016-03-01', checkpointLocation = "."))

args <- commandArgs(trailingOnly=TRUE)
# args <- c("TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset", "GSE9891_eset", "GSE26712_eset", "aaces.eset")

# Get dataset name of the esets that we're loading
datasets <- c("TCGA", "mayo.eset", "Yoshihara", "Tothill", "Bonome", "aaces.eset")

# List all kmeans files
allkmeansFiles <- list.files("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/")

# List all NMF files
allnmfFiles <- list.files("2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/")
############################################
# Load data, process data and output table
############################################
# Get kmeans and nmf cluster membership files
clusterMembership <- c()
for (dataset in 1:length(args)) {
  
  # Load kmeans file
  kmeansFile <- allkmeansFiles[grepl(args[dataset], allkmeansFiles)]
  kmeansMemb <- read.table(paste("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/", kmeansFile, sep = ""), 
                           row.names = 1, header = T, sep = ",")
  colnames(kmeansMemb) <- c("ClusterK3_kmeans", "ClusterK4_kmeans")
  			
  # Load nmf file
  nmfFile <- allnmfFiles[grepl(args[dataset], allnmfFiles)]
  nmfFile <- nmfFile[grepl("mapped", nmfFile)]
  nmfMemb <- read.table(paste("2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/", nmfFile, sep = ""), 
                        row.names = 1, header = T, sep = ",")
  colnames(nmfMemb) <- c("ClusterK3_NMF", "ClusterK4_NMF")
  
  # get a vector of the dataset name of the same length as the number of samples
  Dataset <- rep(datasets[dataset], nrow(kmeansMemb))
  names(Dataset) <- rownames(kmeansMemb)
  
  # Make sure the nmf cluster membership is in the same order
  nmfMemb <- nmfMemb[match(rownames(nmfMemb), rownames(kmeansMemb)),]
  
  # Combine all the information
  datasetMemb <- cbind(Dataset, kmeansMemb, nmfMemb) 
  
  #Compile all datasets clustermembership information
  if(dataset == 1) {
    clusterMembership <- datasetMemb
  } else {
    clusterMembership <- rbind(clusterMembership, datasetMemb)
  }
}

colnames(clusterMembership) <- c('Dataset', 'ClusterK2_kmeans', 'ClusterK3_kmeans', 'ClusterK4_kmeans',
                                 'ClusterK2_NMF', 'ClusterK3_NMF', 'ClusterK4_NMF')

# This is a supplemental table, write it to file
write.table(clusterMembership, "2.Clustering_DiffExprs/Tables/ClusterMembership/FullClusterMembership.csv", 
            row.names = T, col.names = NA, sep = ",")

