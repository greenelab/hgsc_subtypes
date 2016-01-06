############################################
# Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will compare clusters identified by k-means and NMF

args <- commandArgs(trailingOnly=TRUE)
# args <- c("TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset")

############################################
#Load Libraries
############################################
library(reshape2)
library(ggplot2)
library(gridExtra)

############################################
# Load Data
############################################
# Get k-means and nmf SAM Scores (moderated t scores; D scores)
DScoreList_kmeans <- list()
DScoreList_nmf <- list()
datasetMembers_kmeans <- list.files(path = "2.Clustering_DiffExprs/Tables/DScores/")
datasetMembers_nmf <- list.files(path = "2.Clustering_DiffExprs/Figures/nmf/DscoreVectors/" )

# Loop over each dataset given as arguments
for (dataset in 1:length(args)) {
  holder <- args[dataset]
  
  # Get the given membership files for NMF and kmeans
  DMembership <- datasetMembers_kmeans[grep(holder, datasetMembers_kmeans)]
  DMembershipNMF <- datasetMembers_nmf[grep(holder, datasetMembers_nmf)]
  
  # Load k-means File
  DScoreList_kmeans[[args[dataset]]] <- read.table(paste("2.Clustering_DiffExprs/Tables/DScores/", 
                                                   DMembership, sep = ""), sep = ",", row.names = 1, header = T)
  
  # Load and process NMF file
  NMFfile <- read.table(paste("2.Clustering_DiffExprs/Figures/nmf/DscoreVectors/", DMembershipNMF, sep = ""), 
                     sep = ",", row.names = 1, header = T)
  colnames(NMFfile) <- c(paste(args[dataset], "_ClusterK3_1_NMF", sep = ""), paste(args[dataset], "_ClusterK3_2_NMF", sep = ""), 
                         paste(args[dataset], "_ClusterK3_3_NMF", sep = ""), paste(args[dataset], "_ClusterK4_1_NMF", sep = ""), 
                         paste(args[dataset], "_ClusterK4_2_NMF", sep = ""), paste(args[dataset], "_ClusterK4_3_NMF", sep = ""), 
                         paste(args[dataset], "_ClusterK4_4_NMF", sep = ""))
  
  # Store NMF file
  DScoreList_nmf[[args[dataset]]] <- NMFfile
}

############################################
# Across Method Correlations within Datasets
############################################
Method_Cor <- list()
for (dataset in 1:length(DScoreList_kmeans)) {
  
  # Get within dataset correlation matrix between kmeans and NMF methods
  Method_Cor[[args[dataset]]] <- cor(DScoreList_kmeans[[dataset]], DScoreList_nmf[[dataset]])
  
  # Name the output plot
  fName <- paste(args[dataset],"_kmeansVnmf.tiff")
  
  # Plot a correlation heatmap
  ggplot(data = melt(Method_Cor[[dataset]]), aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = "white") + 
    scale_fill_gradient2(high = "red", low = "blue", mid = "white", midpoint = 0, limit = c(-1, 1), 
                         name = "Moderated t Score\nPearson\nCorrelation") +
    theme_minimal() + theme(legend.position = "none", axis.text.x = element_blank(), 
                            axis.text.y = element_blank()) + coord_fixed() +
    labs(x = "", y = "", title = "") + theme(plot.margin = unit(rep(0.1, 4), "cm"))
  
  # Save plot to file (Figure 4)
  ggsave(file.path("2.Clustering_DiffExprs/", "Figures", "kmeans_v_nmf", "datasets", fName), width=8.25, height=8.25)
}
