############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will input a series of datasets and output sample by sample 
# correlation matrix heatmaps

suppressMessages(library(checkpoint))
suppressMessages(checkpoint('2016-03-01', checkpointLocation = "."))

args <- commandArgs(trailingOnly=TRUE)
# args <- c(2, 4, 123, "Figures/CorrelationMatrix/", "TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset", "GSE9891_eset", "GSE26712_eset", "aaces.eset")
###########################################
# Load Libraries
############################################
library(curatedOvarianData)
library(gplots)
library(RColorBrewer)
library(reshape)

# The script loads the ovarian cancer datasets
source("1.DataInclusion/Scripts/Functions/LoadOVCA_Data.R")

# This script loads the heatmap.3 function
source("2.Clustering_DiffExprs/Scripts/Functions/heatmap3.R")

############################################
#Constants
############################################
k = as.numeric(paste(args[1])) 
k2 = as.numeric(paste(args[2]))  #Number of clusters 
SEED = as.numeric(paste(args[3])) 
filepath = as.character(paste(args[4]))
set.seed(SEED)
krange = k:k2

############################################
# Functions
############################################
CorMatrixOrder <- function (DataSet, ClusterMemb, ClusterColumn) {  
  # ~~~~~~~~~~~~~~
  # Outputs a dataframe of cluster membership for k = 3 and k = 4 using Global MAD genes
  #
  # Args: 
  # DataSet: A gene expression matrix with genes as rows and samples as columns
  # ClusterMemb: a dataframe with cluster assignments
  # ClusterColumn: which column the cluster memberships are stored
  # 
  # Returns:
  # The order by which the samples in each correlation matrix heatmap are plotted
  # ~~~~~~~~~~~~~~

  # Initialize an order for the correlation matrix
  DataSet_order <- c()
  
  # Determine the number of unique clusters in the given cluster column
  uniqueClusters <- length(unique(ClusterMemb[ ,ClusterColumn]))
  
  # For each cluster, perform hierarchical clustering to determine the order, 
  # within each cluster assignment, for presentation in the correlation matrix
  for (clus in 1:uniqueClusters) {
    # Select the gene expression of the samples with the cluster assignment
    cluster.subset <- DataSet[ ,ClusterMemb[ ,ClusterColumn] == clus]
    
    # Observe the distance matrix for these samples
    cluster.dist <- dist(t(cluster.subset))
    
    # Perform heirarchical clustering on this distance matrix
    cluster.hclust <- hclust(cluster.dist)
    
    # Order the subset of samples based on the order resulting from the heirarchical clustering
    cluster.order <- cluster.subset[ ,cluster.hclust$order]
    cluster.order <- as.matrix(cluster.order)
    
    # Reorder the gene expression matrices based on first the k means clusters and next based on 
    # the heirarchical clusters
    DataSet_order <- cbind(DataSet_order, cluster.order)
  }
   
  # Obtain the correlation matrix and return
  corReady <- cor(DataSet_order)
  return(corReady)
}
                          
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

# Load the GlobalMAD genes
GlobalMAD <- read.csv(file = "1.DataInclusion/Data/Genes/GlobalMAD_genelist.csv", header = T, 
                      stringsAsFactors = F)

############################################
# Build Sample by Sample Correlation Matrices 
############################################
for (dataset in 1:length(ExpData)) {
  # Store the dataset name
  name <- names(ExpData)[[dataset]]
  
  # Obtain cluster membership file
  membfile <- paste("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/KMembership_", 
                    name, ".csv", sep = "")
  
  # Read in cluster membershipfile
  ClusterAssign <- read.csv(file = membfile, row.names = 1)

  # Loop over the range of k
  for (centroid in 1:length(krange)) {
    
    # Reorder the gene expression dataset
    dataclustertmp <- CorMatrixOrder(ExpData[[dataset]], ClusterAssign, as.numeric(paste(centroid)))   
    
    # Store side color information depending on cluster num
    colorPool <- c('skyblue1', 'tomato')
    if (centroid == 2) {
      colorPool <- c(colorPool, "springgreen")
    } else if (centroid == 3) {
      colorPool <- c(colorPool, "springgreen", "violet")
    }
    
    # The side colors are different in Bonome since the clusters are not similar
    if (name == "GSE26712_eset" & centroid == 1) {
      colorPool <- c("gray75", "gray60")
    } else if (name == "GSE26712_eset" & centroid == 2) {
      colorPool <- c("gray75", "gray60", "gray45")
    } else if (name == "GSE26712_eset" & centroid == 3) {
      colorPool <- c("gray75", "gray60", "gray45", "gray30")
    }
    
    # Prepare side colors from the color pool
    sideColors <- as.matrix(rep(colorPool, as.vector(table(ClusterAssign[ ,centroid]))))
    
    # Store plot names
    heatFile <- paste(name, krange[centroid], sep = "_")
    heatTitle <- c("Correlation Matrix\n", paste(unlist(strsplit(name, "_"))[1], "k =", 
                                                 krange[centroid], sep = " "))
    
    # Store the color range for the correlation colors
    color_range <- colorRampPalette(c("blue", "White", "red"))(n = 2999)
    
    # Store color breaks
    mid_cut = mean(dataclustertmp) - (3 * sd(dataclustertmp))
    mid_cut_high = quantile(dataclustertmp, probs = .65)
    col_breaks = c(seq(min(dataclustertmp), mid_cut, length = 100), seq(mid_cut,mid_cut_high, length = 1600), 
                   seq(mid_cut_high, 1, length = 1300))

    # Compile a heatmap for each dataset for each k
    png(paste("2.Clustering_DiffExprs/Figures/cormatrix/", heatFile, ".png", sep = ""), 
        width = 1000, height = 900)
    
    # Plot the heatmap
    heatmap.3(dataclustertmp, symm = T, trace = 'none', Rowv = NULL, Colv = "Rowv", dendrogram = 'none', 
              key = F, labRow = F, labCol = F, col = color_range, breaks = col_breaks, main = heatTitle, 
              ColSideColors = sideColors, ColSideColorsSize = 2, RowSideColors = t(sideColors), 
              RowSideColorsSize = 2)
    
    dev.off()
  }
}

