############################################
# High-grade serous ovarian cancer subtypes are similar across populations
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script stores several functions that are required for several analyses 
# in the pipeline.

############################################
#Functions
############################################
MADgenes <- function (Dataset, numGenes = 1500) {
  # ~~~~~~~~~~~~~~
  # Outputs a given number of most variably expressed genes in a given dataset.
  #
  # Args: 
  # Dataset: A gene expression matrix with genes as rows and samples as columns
  # numGenes: The number of genes to return, ordered by variability
  #
  # Returns:
  # A string of the top n most variably expressed genes
  # ~~~~~~~~~~~~~~
  
  # Get the median expression of all genes in a dataset
  Dataset.MAD <- data.frame(apply(Dataset, 1, 
                                  function(x) {
                                    thisMedian <- median(x, na.rm = TRUE)
                                    result <- median(abs(x - thisMedian), na.rm = TRUE)
                                    return(result)
                                    }
                                  ))
  
  #Ensure numGenes is numeric
  numGenes <- as.numeric(paste(numGenes))
  
  # Order the genes and select the top MAD genes
  Dataset.MAD.numGenes <- rownames(Dataset)[order(Dataset.MAD, decreasing = T)][1:numGenes]
  return(Dataset.MAD.numGenes)
}


KmeansGlobal <- function (Dataset, FilteredGenes, kmin = 3, kmax = 4, starts = 20) {
  # ~~~~~~~~~~~~~~
  # Performs k-means clustering for a given Dataset over a given k-range
  #
  # Args: 
  # Dataset: A gene expression matrix with genes as rows and samples as columns
  # FilteredGenes: A string of genes to subset the given dataset
  # kmin: An integer of the first k-means centroid assignment
  # kmax: An integer of the last k-means centroid assignment
  # starts: The number of k-means iterations to consider
  #
  # Returns:
  # A dataframe of cluster membership for each sample in the dataset
  # ~~~~~~~~~~~~~~
  
  #transpose the dataset after filtering a defined set of genes
  FilteredGenes <- intersect(rownames(Dataset), FilteredGenes)
  DataSet_Global <- t(Dataset[FilteredGenes,])
  
  #run kmeans over each number of clusters inputted from kmin to kmax
  krange <- seq(kmin, kmax)
  DataSetClusterMembership <- data.frame()
  for (clusterK in 1:length(krange)) {
    #Print status to screen
    cat("Clustering...", paste("k = ", krange[clusterK], sep = ""), "\n")
    
    #run the kmeans function
    dta_k <- kmeans(DataSet_Global, krange[clusterK], nstart = starts)
    
    #Save the cluster membership info
    if (clusterK == 1) {
      DataSetClusterMembership <- as.data.frame(dta_k$cluster)
    } else {
      DataSetClusterMembership <- cbind(DataSetClusterMembership, as.data.frame(dta_k$cluster))
    }
     
  }
  DataSetClusterMembership <- data.frame(DataSetClusterMembership)
  colnames(DataSetClusterMembership) <- paste("ClusterK", seq(kmin, kmax), sep = "")
  return(DataSetClusterMembership)
}


getSamNames <- function (krange)
{
  # ~~~~~~~~~~~~~~
  # Retrieves column names from the SAM results
  #
  # Args: 
  # krange: range(kmin, kmax)
  #
  # Returns:
  # A string of column names
  # ~~~~~~~~~~~~~~
  
  # Get SAM colnames
  krange_grid <- expand.grid(c("SAMD", "SAMP"), min(krange):max(krange))
  
  # Ensure both columns are characters
  krange_grid[ ,1] <- as.character(krange_grid[,1])
  krange_grid[ ,2] <- as.character(krange_grid[,2])
  
  # Output custom column names based on krange
  samColNames <- c()
  for (gridrow in 1:nrow(krange_grid)) {
    samColNames <- c(samColNames, paste(krange_grid[gridrow,1], krange_grid[gridrow,2], sep="_"))
  }
  
  return(samColNames) 
}


RunSam <- function (DataList, ClusterList, Method = "d.stat", type = "normal", kFDR = 0.01) {
  # ~~~~~~~~~~~~~~
  # Perform SAM on all the input datasets for each k given
  # Using the siggenes package
  #
  # Args: 
  # DataList: a list holding the expression values for all datasets being considered
  # ClusterList: a list holding the cluster membership files for the datasets
  # Method: internal to sam() from siggenes package
  # type: controls whether or not a Delta is computed; defaults to no Delta
  # kFDR: the False Discovery Rate to determine the appropriate Delta in SAM
  # 
  # Returns:
  # A list of SAM results for krange in ClusterList
  # Each element of the SAM list is a dataframe of d statistics and p values
  # If type = "Delta" then the second element of the SAM list holds delta values
  # ~~~~~~~~~~~~~~
  
  #Determine the number of clusters were considered from kmeans
  numClusts <- ncol(ClusterList[[1]])
  
  SamList <- list()
  Deltas <- c()
  for (dataset in 1:length(DataList)) {
    cat("Running SAM on:", names(DataList)[dataset], "\n")
    
    for (centroid in 1:ncol(ClusterList[[dataset]]))
    {
      # Determine how many different cluster assignments there are
      uniqueIDs <- sort(unique(ClusterList[[dataset]][ ,centroid]))
      thisCol <- colnames(ClusterList[[dataset]])[centroid]
      
      # For each column of cluster assignments, loop over the unique elements, 
      # performing a 1 versus all SAM.
      SamD <- c()
      for (kclus in 1:length(uniqueIDs)) {
        
        # Build a comparison vector where the cases are 1's and the controls are 0's
        comparisonVector <- rep(0, nrow(ClusterList[[dataset]]))
        comparisonVector[ClusterList[[dataset]][ ,centroid] == uniqueIDs[kclus]] <- 1
        
        # Run the sam function
        cat("      For k =", length(uniqueIDs), ", Control: Cluster", kclus, "\n")
        tmpsam <- sam(DataList[[dataset]], comparisonVector, method=Method)
        
        # Extract the delta required to observe a given fdr threshhold
        if (type == "Delta") {
          # Find and report the delta for the given number of genes
          delta <- findDelta(tmpsam, fdr = kFDR, prec = 15)
          
          if (is.matrix(delta)) {
            d <- delta[1,1]
          } else {
            d <- delta[1]
          }
          
          names(d) <- paste(names(DataList)[dataset], thisCol, uniqueIDs[kclus], sep = "_")
          Deltas <- c(Deltas, d)
        }
        
        # Extract the d statistic and p value from the sam
        SamD <- cbind(SamD, tmpsam@d, tmpsam@p.value)
      }
      
      SamD <- data.frame(SamD)
      colnames(SamD) <- getSamNames(uniqueIDs)
      
      # Store the SAM results in a list separated by dataset
      listIndex <- ((dataset - 1) * numClusts) + centroid
      SamList[[listIndex]] <- SamD
      
      # Name the list
      thisName <- paste(names(DataList)[dataset], thisCol, sep="_")
      names(SamList)[length(SamList)] <- thisName
    }
  }  
  
  if(type == "Delta") {
    SamList <- list(SamList, Deltas)
  }
  
  return(SamList)
}


AssignReference <- function (Reference, Cluster, Cor, ClusterList = Clusters) {
  # ~~~~~~~~~~~~~~
  # This piece of code below will assign the appropriate clusters to a reference 
  # cluster. (for us, "TCGA" is the reference)
  #
  # Args: 
  # Reference: a string identifier for the dataset of interest to map clusters
  # Cluster: Either "ClusterK3" or "ClusterK4" to begin mapping
  # Cor: A list of correlation matrices of within cluster correlations
  # ClusterList: A list of cluster membership dataframes
  # 
  # Returns:
  # New cluster assignments for the reference dataset
  # ~~~~~~~~~~~~~~
  
  # Test if input for "Cluster" is given in the correct format
  if (Cluster != "ClusterK3" & Cluster != "ClusterK4") {
    stop("Please input for Cluster either 'ClusterK3' or 'ClusterK4'")
  }
  
  # Get the within dataset correlation matrix for the given reference 
  corMat <- Cor[[grep(Reference, names(Cor))]]
  
  require(reshape2)
  # Melt the corerlation matrix and assign column names
  corMelted <- melt(corMat)
  colnames(corMelted) <- c("Var1", "Var2", "value")
  
  # Order the melted correlation matrix by decreasing correlations and reassign colnames
  corMeltOrder <- corMelted[order(corMelted$value, decreasing = T),]
  colnames(corMeltOrder) <- c("Var1", "Var2", "value")
  
  # Remove all self correlations comparisons
  newCorframe <- c()
  for (sample in 1:nrow(corMeltOrder)) {
    if (corMeltOrder$Var1[sample] != corMeltOrder$Var2[sample]) {
      newCorframe <- rbind(newCorframe, corMeltOrder[sample, ])
    }
  }
  
  # Observe the top six correlating clusters
  topCorClusters <- newCorframe[1:6, 1:3]
  
  # Split the cluster names and store in list
  splitClusterA <- strsplit(as.character(topCorClusters$Var1), "_")
  splitClusterB <- strsplit(as.character(topCorClusters$Var2), "_")
  
  # Initialize the new mapping based on decreasing correlation
  clusMap <- c()
  for (ID in 1:length(splitClusterA)) {
    clusA <- splitClusterA[[ID]][3] #ClusterK3 or ClusterK4
    clusB <- splitClusterB[[ID]][3]
    numA <- splitClusterA[[ID]][5] #1, 2, 3, or 4
    numB <- splitClusterB[[ID]][5]
    newAssgn <- cbind(clusA, numA, clusB, numB)
    clusMap <- rbind(clusMap, newAssgn)
  }
  
  # Focus on mapping the given Cluster (Either ClusterK3 or ClusterK4)
  NewMapping <- clusMap[clusMap[,1] == Cluster,]
  
  # Obtain the original cluster membership dataframe for the given reference dataset
  ref <- ClusterList[[grep(Reference, names(ClusterList))]]
  for (map in 1:nrow(NewMapping)) {
    # Temporarily replace original cluster assignments
    ref[ ,1][ref[ ,1] == as.integer(NewMapping[map,2])] <- paste(map, "tmp", sep = "")
    ref[ ,2][ref[ ,2] == as.integer(NewMapping[map,4])] <- paste(map, "tmp", sep = "")
  }
  
  # If there is no "tmp" found in the column, make it cluster 4
  ref[ ,2][!grepl("tmp", ref[ ,2])] <- 4
  
  # Return all cluster assignments to numeric values
  ref[ ,1] <- as.integer(unlist(strsplit(ref[ ,1], "tmp")))
  ref[ ,2] <- as.integer(unlist(strsplit(ref[ ,2], "tmp")))
  
  # Reassign the ClusterList and return new assignments
  ClusterList[[grep(Reference, names(ClusterList))]] <- ref
  return(ClusterList)
}


AssignReference_NMF <- function (kmeans_dscore_dir, nmf_Dlist, nmf_cluster_list, Reference = "TCGA") {
  # ~~~~~~~~~~~~~~
  # This function will use a the dscores extracted in the kmeans mapping to map nmf clusters
  #
  # Args: 
  # kmeans_dscore_location: the directory where the cluster specific d scores are stored
  # nmf_Dlist: a list holding d score vectors for the nmf algorithm
  # nmf_cluster_list: a list of cluster assignments that has already been mapped to TCGA
  # Reference: defaults to TCGA; which population to use for mapping clusters across populations
  #
  # Returns:
  # New cluster assignments for the NMF reference dataset
  # ~~~~~~~~~~~~~~
  
  # Upload the dscore vector for the reference
  d_score_file <- paste(kmeans_dscore_dir, 
                        list.files(kmeans_dscore_dir)[grepl(Reference, list.files(kmeans_dscore_dir))],
                        sep = "/")
  kmeans_dscore <- read.table(d_score_file, sep = ",", header = T, row.names = 1, stringsAsFactors = F)
  
  # Subset the Dlist to only the reference dataset
  Dlist_subset <- Dlist[[grep(Reference, names(Dlist))]]
  colnames(Dlist_subset) <- paste("nmf_", colnames(Dlist_subset), sep = "")
  
  # Subset the list of cluster memberships
  clus_memb <- nmf_cluster_list[[grep(Reference, names(nmf_cluster_list))]]
  
  # Get correlation matrix of nmf and kmeans results
  cor_matrix <- cor(Dlist_subset, kmeans_dscore)
  
  # Melt the correlation matrix
  cor_matrix_melted <- melt(cor_matrix)
  cor_matrix_melted <- cor_matrix_melted[order(cor_matrix_melted$value, decreasing = T),]
  
  # Get cluster specific matrices
  cor_matrix_k3 <- cor_matrix_melted[grep("K3", cor_matrix_melted$Var1), ]
  cor_matrix_k3 <- cor_matrix_k3[grep("K3", cor_matrix_k3$Var2), ]
  
  cor_matrix_k4 <- cor_matrix_melted[grep("K4", cor_matrix_melted$Var1), ]
  cor_matrix_k4 <- cor_matrix_k4[grep("K4", cor_matrix_k4$Var2), ]
  
  # Make Maps
  k3_maps <- c()
  for (row in 1:3) {
    nmf_mapping <- unlist(strsplit(as.character(cor_matrix_k3$Var1)[row], "_"))
    nmf_mapping <- nmf_mapping[length(nmf_mapping)]
    
    k_mapping <- unlist(strsplit(as.character(cor_matrix_k3$Var2)[row], "_"))
    k_mapping <- k_mapping[length(k_mapping)]
    k3_maps <- rbind(k3_maps, c(as.numeric(paste(nmf_mapping)), as.numeric(paste(k_mapping))))
  }
  
  # We only care about the top three highest correlations
  k4_maps <- c()
  for (row in 1:3) {
    nmf_mapping <- unlist(strsplit(as.character(cor_matrix_k4$Var1)[row], "_"))
    nmf_mapping <- nmf_mapping[length(nmf_mapping)]
    
    k_mapping <- unlist(strsplit(as.character(cor_matrix_k4$Var2)[row], "_"))
    k_mapping <- k_mapping[length(k_mapping)]
    k4_maps <- rbind(k4_maps, c(as.numeric(paste(nmf_mapping)), as.numeric(paste(k_mapping))))
  }
  
  for(assgn in 1:nrow(k3_maps)) {
    map <- k3_maps[assgn, 1]
    mapk4 <- k4_maps[assgn, 1]
    clus_memb$ClusterK3[grep(paste("^", map, "$", sep = ""), clus_memb$ClusterK3)] <- paste(k3_maps[assgn, 2], "tmp", sep = "_")
    clus_memb$ClusterK4[grep(paste("^", mapk4, "$", sep = ""), clus_memb$ClusterK4)] <- paste(k4_maps[assgn, 2], "tmp", sep = "_")
  }

  # If there is no "tmp" found in the column, make it cluster 4
  clus_memb[ ,2][!grepl("tmp", clus_memb[ ,2])] <- 4
  
  # Return all cluster assignments to numeric values
  clus_memb[ ,1] <- as.integer(unlist(strsplit(clus_memb[ ,1], "_tmp")))
  clus_memb[ ,2] <- as.integer(unlist(strsplit(clus_memb[ ,2], "_tmp")))
  
  # Reassign the clusters and return
  nmf_cluster_list[[grep(Reference, names(nmf_cluster_list))]] <- clus_memb
  return(nmf_cluster_list)
}

MapClusters <- function (DistMatrixList, Reference = "TCGA") {
  # ~~~~~~~~~~~~~~
  # This function will use a distance matrix and hierarchical clustering to map clusters 
  # across populations.
  #
  # Args: 
  # DistMatrixList: a list of distance matrices for each value of k
  # Reference: defaults to TCGA; which population to use for mapping clusters across populations
  #
  # Returns:
  # New cluster assignments for the reference dataset
  # ~~~~~~~~~~~~~~
  
  # This is initializing the list that we will be returning
  resultList <- list()
        
  # Loop over each number of centroids
  for (centroid in substr(names(DistMatrixList), 2, 2)) {
    
    # Convert the correlation matrix to a dist object
    thisDistance <- DistMatrixList[[paste("K", centroid, sep = "")]]

    # Perform hierarchical clusterhing
    thisHC <- hclust(thisDistance)
    
    # Group the subtypes into k groups
    thisGroups <- as.data.frame(cutree(thisHC, k = centroid))
    
    result <- list()
    # Determine the group ID for each of the reference subtypes
    for (hclusts in 1:as.numeric(centroid)) {
      pattern <- paste(Reference, ".*_", hclusts, sep = "")
      rowID <- grep(pattern, rownames(thisGroups))
      oldGroupLabel <- thisGroups[rowID, 1]
      
      result[[hclusts]] <- rownames(thisGroups)[thisGroups[ ,1] == oldGroupLabel]
    }
    resultList[[paste("K", centroid, sep = "")]] <- result
  }
  return(resultList)
}



runNMF <- function (Data, k, fname, KClusterAssign, nruns = 10) {
  # ~~~~~~~~~~~~~~
  # This function will perform NMF clustering on gene expression data
  #
  # Args: 
  # Data: Sample by gene expression dataframe
  # k: The number of clusters NMF should find
  # fname: The output file name to store result
  # KClusterAssign: The k-means cluster assignments
  # nruns: How many iterations NMF should perform 
  #
  # Returns:
  # 1) Consensus matrices with silhouette and k-means tracks
  # 2) Cluster membership files for each sample
  # ~~~~~~~~~~~~~~
  
  # Determine the minimum of the data
  minimum <- min(Data)
  
  # If the minimum is less than zero, renormalize because NMF cannot handle negative values
  if (minimum < 0) {
    # Renormalize removing negative values
    const <- -1 * minimum
    Data <- apply(Data, 2, function(x){x + const})
  }
  
  # Perform nmf clustering
  clus.nmf <- nmf(Data, k, nrun = nruns, .options = "t")
  
  # Get Clustermembership predictions
  clusterMemb <- predict(clus.nmf)
  
  # Write the consensus mapping as a figure
  png(paste("2.Clustering_DiffExprs/Figures/nmf/ConsensusMaps/",fname, ".png", sep = ""), 
       width = 700, height = 570)
  
  # The plot colors are different according to the number of clusters
  if (k == 3) {
    consensusmap(clus.nmf, labCol = NA, labRow = NA, tracks = c("silhouette:"), 
                 annCol = list("kmeans" = as.character(paste(KClusterAssign$ClusterK3))), 
                 annColors = list("kmeans" = c("blue", "red", "green")), main = "", 
                 fontsize = 20, treeheight = 60)
  } else {
    consensusmap(clus.nmf, labCol = NA, labRow = NA, tracks = c("silhouette:"), 
                 annCol = list("kmeans" = as.character(paste(KClusterAssign$ClusterK4))), 
                 annColors = list("kmeans" = c("blue", "red", "green", "purple")), main = "", 
                 fontsize = 20, treeheight = 60)
  }
  dev.off()
  
  return(as.character(clusterMemb))
}

