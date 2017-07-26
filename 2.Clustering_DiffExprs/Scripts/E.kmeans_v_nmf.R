############################################
# Cross-population analysis of high-grade serous ovarian cancer reveals
# only two robust subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will compare clusters identified by k-means and NMF

args <- commandArgs(trailingOnly=TRUE)
# args <- c("TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset",
# "GSE9891_eset", "aaces.eset")

############################################
#Load Libraries
############################################
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

############################################
# Load Data
############################################
# Get k-means and nmf SAM Scores (moderated t scores; D scores)
DScoreList_kmeans <- list()
DScoreList_nmf <- list()

kmeansclusters.path <- file.path("2.Clustering_DiffExprs", "Tables", "DScores")
datasetMembers_kmeans <- list.files(path = kmeansclusters.path)

nmfclusters.path <- file.path("2.Clustering_DiffExprs",
                              "Figures", "nmf", "DscoreVectors")
datasetMembers_nmf <- list.files(path = nmfclusters.path)

# Loop over each dataset given as arguments
for (dataset in args) {

  # Get the given membership files for NMF and kmeans
  DMembership <- datasetMembers_kmeans[grep(dataset, datasetMembers_kmeans)]
  DMembershipNMF <- datasetMembers_nmf[grep(dataset, datasetMembers_nmf)]

  # Load k-means File
  kmeans.path <- file.path("2.Clustering_DiffExprs", "Tables",
                           "DScores", DMembership)
  DScoreList_kmeans[[dataset]] <- read.table(kmeans.path,
                                             sep = ",", row.names = 1,
                                             header = T)

  # Load and process NMF file
  NMFpath <- file.path("2.Clustering_DiffExprs", "Figures", "nmf",
                       "DscoreVectors", DMembershipNMF)
  NMFfile <- read.table(NMFpath, sep = ",", row.names = 1, header = T)

    colnames(NMFfile) <- c(paste(dataset, "_ClusterK2_1_NMF", sep = ""),
                           paste(dataset, "_ClusterK2_2_NMF", sep = ""),
                           paste(dataset, "_ClusterK3_1_NMF", sep = ""),
                           paste(dataset, "_ClusterK3_2_NMF", sep = ""), 
                           paste(dataset, "_ClusterK3_3_NMF", sep = ""),
                           paste(dataset, "_ClusterK4_1_NMF", sep = ""), 
                           paste(dataset, "_ClusterK4_2_NMF", sep = ""),
                           paste(dataset, "_ClusterK4_3_NMF", sep = ""), 
                           paste(dataset, "_ClusterK4_4_NMF", sep = ""))
  
  # Store NMF file
  DScoreList_nmf[[dataset]] <- NMFfile
}

############################################
# Across Method Correlations within Datasets
############################################
Method_Cor <- list()
k_list <- c('K2', 'K3', 'K4')
all_centroid_plots <- list()

for (dataset in 1:length(DScoreList_kmeans)) {
  # Get within dataset correlation matrix between kmeans and NMF methods
  Method_Cor[[args[dataset]]] <- cor(DScoreList_kmeans[[dataset]],
                                     DScoreList_nmf[[dataset]])
  
  # Name the output plot
  fName <- paste(args[dataset],"_kmeansVnmf.png")
  
  # Melt the dataset and extract centroid specific data
  melted_cor <- melt(Method_Cor[[dataset]])
  
  # Make a glob for each centroid assignment
  data_iter <- 1
  for (k in k_list) {
    melted_cor_sub <- melted_cor[grepl(k, melted_cor[, 1]), ]
    melted_cor_sub <- melted_cor_sub[grepl(k, melted_cor_sub[, 2]), ]
    num_clus <- gsub('K', '', k)
    
    # Plot a correlation heatmap
    g <- ggplot(data = melted_cor_sub, aes(x = Var1, y = Var2, fill = value)) + 
      geom_tile(color = "white") + 
      scale_fill_gradient2(high = "red", low = "blue", mid = "white",
                           midpoint = 0, limit = c(-1, 1)) +
      scale_x_discrete(labels = paste(1:num_clus)) + 
      scale_y_discrete(labels = paste(1:num_clus)) +
      xlab("") + 
      ylab("") +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(), 
            panel.background = element_blank(),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            plot.background = element_blank(), 
            legend.position = "none",
            axis.text.y = element_text(face = 'bold', size = 18, angle = 90,
                                       hjust = 0.4, colour = 'black'),
            axis.text.x = element_text(face = 'bold', size = 18,
                                       colour = 'black'),
            plot.margin = unit(c(0, -0.3, -0.3, -0.3), 'cm'))
    
    all_centroid_plots[[data_iter]] <- g
    data_iter <- data_iter + 1
  }
  
  # Build a string to evaluate
  plot_eval <- 'full_grobs <- grid.arrange('
  for (p in 1:length(all_centroid_plots)) {
    plot_eval <- paste0(plot_eval, 'all_centroid_plots[[', p, ']]', ', ')
  }
  plot_eval <- paste0(plot_eval, 'ncol = 3', ',nrow = 1)')
  
  png(file.path("2.Clustering_DiffExprs", "Figures", "kmeans_v_nmf",
                "datasets", fName), width = 800, height = 280)
  eval(parse(text = plot_eval))
  dev.off()
  all_centroid_plots <- list()
}

