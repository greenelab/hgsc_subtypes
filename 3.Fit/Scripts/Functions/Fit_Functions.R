############################################
# Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script stores a series of functions required for Goodness of Fit analyses

############################################
# Functions
############################################
KMeansElbow <- function (DataSet, Name, kmax = 8, itermax = 20, starts = 20) {
  # ~~~~~~~~~~~~~~
  # Output elbow quantification plots
  #
  # Args: 
  # DataSet: A gene expression matrix with genes as rows and samples as columns
  # Name: The name of the eset 
  # kmax: The maximum number of clusters to find (from 0 to kmax)
  # itermax: The number of iterations to determine optimal kmeans clusters
  # starts: The number of random starts in the kmeans algorithm
  #
  # Returns:
  # A graph depicting the AIC and BIC from kmin to kmax and will determine the optimal k
  # ~~~~~~~~~~~~~~
  
  # Ensure the dataset is a matrix
  DataSet <- as.matrix(DataSet)
  
  tmp_totwss <- rep(0, kmax)
  tmp_kmfit <- list()
  # Perform kmeans for each k from 1 to kmax and store the total within sum of squares for each k
  for (centroid in 1:kmax) {
    
    # Perform kmeans finding centroids from 1 to kmax
    kclus <- kmeans(t(DataSet), centers = centroid, iter.max = itermax, nstart = starts)
    
    # Store the total within sum of squares
    tmp_totwss[centroid] <- kclus$tot.withinss
    tmp_kmfit[[centroid]] <- kclus
  }
  
  # Calculate AIC and BIC and store both in a dataframe
  tmp_AICBIC <- sapply(tmp_kmfit, function(fit) {
    m = ncol(fit$centers)
    n = length(fit$cluster)
    k = nrow(fit$centers)
    D = fit$tot.withinss
    
    # Determine and output the AIC and BIC for all clusters
    return(data.frame(AIC = D + 2*m*k, BIC = D + log(n)*m*k))
    })
  
  # Prepare the plots and plotting margins
  tiff(file.path(paste("3.Fit/Figures/AICBIC/", Name, "_AICBIC.tiff", sep = "")), 
       width = 2400, height = 1200)
  par(mfrow = c(1,2))
  par(mar = rep(2.4, 4))
  par(cex = 3)
  
  # Plot AIC
  plot(seq(1, kmax), tmp_AICBIC[1, ], xlab = "", main = "", ylab = "", pch = 20, cex = 2, cex.axis = 1.2)
  
  # Identify AIC elbow
  v = -diff(as.numeric(tmp_AICBIC[1, ]))
  nv = length(v)
  tmp_fom_A = (v[1:(nv - 1)])/abs(v[2:nv])
  
  # The max difference
  nclus = which.max(tmp_fom_A) + 1
  points(nclus, tmp_AICBIC[1, ][nclus], col = 2, pch = 20, cex = 2)
  
  # Plot BIC
  plot(seq(1, kmax), tmp_AICBIC[2, ], xlab = "", main = "", ylab = "", pch = 20, cex = 2, cex.axis = 1.2)
  
  # Identify BIC elbow
  v1 = -diff(as.numeric(tmp_AICBIC[2, ]))
  nv1 = length(v1)
  tmp_fom_B = (v1[1:(nv1-1)])/abs(v1[2:nv1])
  
  # The max BIC difference
  nclus2 = which.max(tmp_fom_B) + 1
  points(nclus2, tmp_AICBIC[2, ][nclus2], col = 2,pch = 20, cex = 2)
  dev.off()
}

#This function outputs cluster specific silhouette widths for each dataset for each k
RunSilhouette <- function (dataset, Name, MADgenes) {
  # ~~~~~~~~~~~~~~
  # Output silhouette width plots
  #
  # Args: 
  # dataset: A gene expression matrix with genes as rows and samples as columns
  # Name: The name of the eset 
  # MADgenes: A string of most variably expressed genes
  # kmax: The maximum number of clusters to find (from 0 to kmax)
  # itermax: The number of iterations to determine optimal kmeans clusters
  # starts: The number of random starts in the kmeans algorithm
  #
  # Returns:
  # Silhouette width plots
  # ~~~~~~~~~~~~~~
  
  # First check to see if the global mad genes are able to intersect with the given dataset
  if (length(intersect(rownames(dataset), MADgenes)) != length(MADgenes)) {
    subGlobalMAD <- intersect(rownames(dataset), MADgenes)
    dataUse <- t(dataset[subGlobalMAD, ])
  } else {
    dataUse <- t(dataset[MADgenes, ])
  }
  
  # Get data in appropriate format
  dataExprs <- data.frame(dataUse)
  
  # Obtain distance matrix of expression values
  dataDist <- dist(dataExprs)
  
  # Obtain cluster membership file
  membfile <- list.files(path = "2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/", 
                         pattern = paste(Name, "*", sep = ""))

  # Read in cluster membershipfile
  ClusterAssign <- read.csv(file = file.path("2.Clustering_DiffExprs", "Tables", 
                                             "ClusterMembership", "kmeans", membfile), row.names = 1)
  
  # Perform silhouette width analyses for k = 3 and k = 4
  K2_sil <- silhouette(as.numeric(paste(ClusterAssign$ClusterK2)), dataDist)
  K3_sil <- silhouette(as.numeric(paste(ClusterAssign$ClusterK3)), dataDist)
  K4_sil <- silhouette(as.numeric(paste(ClusterAssign$ClusterK4)), dataDist)
  
  # Output figures
  png(file.path("3.Fit/", "Figures", "Silhouette", paste(Name, "_Silhouette.png", sep="")), 
       width = 1100, height = 800)
  
  # Get appropriate margins
  par(mfrow = c(1, 3))
  par(mar = c(2, 1.5, 1, 1.5))
  par(cex.axis = 0.7, cex = 2, cex.main = 2)
  plot(K2_sil, main = "", xlab = "", sub = "", col = c('skyblue1', 'tomato'), do.col.sort = T, do.n.k = F)
  plot(K3_sil, main = "", xlab = "", sub = "", col = c('skyblue1', 'tomato', 'springgreen'), do.col.sort = T, do.n.k = F)
  plot(K4_sil, main = "",  xlab = "", sub = "", col = c('skyblue1', 'tomato', 'springgreen', 'violet'), do.col.sort = T, do.n.k = F)
  dev.off()
}
