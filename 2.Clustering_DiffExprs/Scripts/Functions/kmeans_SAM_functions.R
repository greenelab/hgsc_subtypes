############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support
# four subtypes
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
                                    result <- median(abs(x - thisMedian),
                                                     na.rm = TRUE)
                                    return(result)
                                    }
                                  ))
  
  #Ensure numGenes is numeric
  numGenes <- as.numeric(paste(numGenes))
  
  # Order the genes and select the top MAD genes
  Dataset.MAD.numGenes <- rownames(Dataset)[order(Dataset.MAD,
                                                  decreasing = TRUE)][1:numGenes]
  return(Dataset.MAD.numGenes)
}


KmeansGlobal <- function (Dataset, FilteredGenes, kmin = 2,
                          kmax = 4, starts = 20) {
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
  
  # Transpose the dataset after filtering a defined set of genes
  FilteredGenes <- intersect(rownames(Dataset), FilteredGenes)
  DataSet_Global <- t(Dataset[FilteredGenes, ])
  
  # Run kmeans over each number of clusters inputted from kmin to kmax
  krange <- seq(kmin, kmax)
  DataSetClusterMembership <- data.frame()
  for (clusterK in 1:length(krange)) {
    # Print status to screen
    cat("Clustering...", paste0("k = ", krange[clusterK]), "\n")
    
    # Run the kmeans function
    dta_k <- kmeans(DataSet_Global, krange[clusterK], nstart = starts)
    
    # Save the cluster membership info
    if (clusterK == 1) {
      DataSetClusterMembership <- as.data.frame(dta_k$cluster)
    } else {
      DataSetClusterMembership <- cbind(DataSetClusterMembership,
                                        as.data.frame(dta_k$cluster))
    }
  }
  DataSetClusterMembership <- data.frame(DataSetClusterMembership)
  colnames(DataSetClusterMembership) <- paste0("ClusterK", seq(kmin, kmax))
  return(DataSetClusterMembership)
}


getSamNames <- function (krange) {
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
  krange_grid[, 1] <- as.character(krange_grid[, 1])
  krange_grid[, 2] <- as.character(krange_grid[, 2])
  
  # Output custom column names based on krange
  samColNames <- c()
  for (gridrow in 1:nrow(krange_grid)) {
    samColNames <- c(samColNames, paste(krange_grid[gridrow, 1],
                                        krange_grid[gridrow, 2], sep = "_"))
  }
  
  return(samColNames) 
}


RunSam <- function (DataList, ClusterList, Method = "d.stat",
            type = "normal", kFDR = 0.01) {
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
  
  # Determine the number of clusters were considered from kmeans
  numClusts <- ncol(ClusterList[[1]])
  
  SamList <- list()
  Deltas <- c()
  for (dataset in 1:length(DataList)) {
    cat("Running SAM on:", names(DataList)[dataset], "\n")
    
    for (centroid in 1:ncol(ClusterList[[dataset]])) {
      # Determine how many different cluster assignments there are
      uniqueIDs <- sort(unique(ClusterList[[dataset]][, centroid]))
      thisCol <- colnames(ClusterList[[dataset]])[centroid]
      
      # For each column of cluster assignments, loop over the unique elements, 
      # performing a 1 versus all SAM.
      SamD <- c()
      for (kclus in 1:length(uniqueIDs)) {
        
        # Build a comparison vector where the cases are 1's and the controls are 0's
        comparisonVector <- rep(0, nrow(ClusterList[[dataset]]))
        comparisonVector[ClusterList[[dataset]][ , centroid] == uniqueIDs[kclus]] <- 1
        
        # Run the sam function
        cat("      For k =", length(uniqueIDs), ", Control: Cluster", kclus, "\n")
        tmpsam <- sam(DataList[[dataset]], comparisonVector, method = Method)
        
        # Extract the delta required to observe a given fdr threshhold
        if (type == "Delta") {
          # Find and report the delta for the given number of genes
          delta <- findDelta(tmpsam, fdr = kFDR, prec = 15)
          
          if (is.matrix(delta)) {
            d <- delta[1, 1]
          } else {
            d <- delta[1]
          }
          
          names(d) <- paste(names(DataList)[dataset], thisCol,
                            uniqueIDs[kclus], sep = "_")
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
      thisName <- paste(names(DataList)[dataset], thisCol, sep = "_")
      names(SamList)[length(SamList)] <- thisName
    }
  }  
  
  if (type == "Delta") {
    SamList <- list(SamList, Deltas)
  }
  
  return(SamList)
}


AssignReference <- function (Reference, Cluster, Cor, ClusterList) {
  # ~~~~~~~~~~~~~~
  # This piece of code below will assign the appropriate clusters to a
  # reference cluster. (for us, "TCGA" is the reference)
  #
  # Args: 
  # Reference: a string identifier for the dataset of interest to map clusters
  # Cluster: Either "ClusterK3" or "ClusterK4" to begin mapping from
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
  
  # Order the melted correlation matrix by decreasing correlations and
  # reassign colnames
  corMeltOrder <- corMelted[order(corMelted$value, decreasing = TRUE),]
  colnames(corMeltOrder) <- c("Var1", "Var2", "value")
  
  # Remove all self correlations comparisons
  newCorframe <- c()
  for (sample in 1:nrow(corMeltOrder)) {
    if (corMeltOrder$Var1[sample] != corMeltOrder$Var2[sample]) {
      newCorframe <- rbind(newCorframe, corMeltOrder[sample, ])
    }
  }
  
  # Observe the top correlating clusters
  unique_clusters <- as.character(unique(newCorframe$Var1))
  
  # Find top correlating cluster for each unique cluster
  topCorClusters <- c()
  for (uniq_clus in unique_clusters) {
    subset_cor <- newCorframe[newCorframe$Var1 == uniq_clus, ]
    topCorClusters <- rbind(topCorClusters, subset_cor[1, ])
  }
  
  # Split the cluster names and store in list
  splitClusterA <- strsplit(as.character(topCorClusters$Var1), "_")
  splitClusterB <- strsplit(as.character(topCorClusters$Var2), "_")
  
  # Initialize the new mapping based on decreasing correlation
  clusMap <- c()
  for (ID in 1:length(splitClusterA)) {
    clusA <- splitClusterA[[ID]][3]  # Cluster K
    clusB <- splitClusterB[[ID]][3]
    numA <- splitClusterA[[ID]][5] # 1, 2, 3, or 4
    numB <- splitClusterB[[ID]][5]
    newAssgn <- cbind(clusA, numA, clusB, numB)
    clusMap <- rbind(clusMap, newAssgn)
  }
  
  # Focus on mapping the given Cluster (given as an argument)
  NewMapping <- clusMap[clusMap[, 3] == Cluster, ]
  
  # Split the multiple comparisons with the input cluster
  Comparisons <- c()
  for (i in unique(NewMapping[, 1])) {
    Comparisons[[i]] <- NewMapping[NewMapping[, 1] == i, ]
  }
  
  # Obtain the original cluster membership dataframe for the given reference
  # dataset
  ref <- ClusterList[[grep(Reference, names(ClusterList))]]
  
  # Focus on K = 3 vs. K = 4 first
  for (map in 1:length(unique(Comparisons[[1]][, 4]))) {
    # Temporarily replace original cluster assignments
    ref[, 2][ref[, 2] == as.integer(Comparisons[[1]][map, 4])] <- paste0(map,
                                                                         "tmp")
    ref[, 3][ref[, 3] == as.integer(Comparisons[[1]][map, 2])] <- paste0(map,
                                                                         "tmp")
  }
  
  # If there is no "tmp" found in the column, make it cluster 4
  ref[, 3][!grepl("tmp", ref[, 3])] <- 4
  
  # The Clusters assigned to K = 2 are mapped automatically
  
  # Return all cluster assignments to numeric values
  ref[, 2] <- as.integer(unlist(strsplit(ref[, 2], "tmp")))
  ref[, 3] <- as.integer(unlist(strsplit(ref[, 3], "tmp")))
  
  # Reassign the ClusterList and return new assignments
  ClusterList[[grep(Reference, names(ClusterList))]] <- ref
  return(ClusterList)
}


AssignReference_NMF <- function (kmeans_dscore_dir, nmf_Dlist,
                                 nmf_cluster_list, Reference = "TCGA") {
  # ~~~~~~~~~~~~~~
  # This function will use a the dscores extracted in the kmeans mapping to map nmf
  # clusters
  #
  # Args: 
  # kmeans_dscore_location: the directory where the cluster specific d scores
  # are stored
  # nmf_Dlist: a list holding d score vectors for the nmf algorithm
  # nmf_cluster_list: a list of cluster assignments that has already been
  # mapped to TCGA
  # Reference: defaults to TCGA; which population to use for mapping clusters
  # across populations
  #
  # Returns:
  # New cluster assignments for the NMF reference dataset
  # ~~~~~~~~~~~~~~
  
  # Upload the dscore vector for the reference
  kmeans_files <- list.files(kmeans_dscore_dir)
  kmeans_ref <- kmeans_files[grepl(Reference, kmeans_files)]
  d_score_file <- paste(kmeans_dscore_dir, kmeans_ref, sep = "/")
  
  kmeans_dscore <- read.table(d_score_file, sep = ",", header = TRUE,
                              row.names = 1, stringsAsFactors = FALSE)
  
  # Subset the Dlist to only the reference dataset
  Dlist_subset <- Dlist[[grep(Reference, names(Dlist))]]
  colnames(Dlist_subset) <- paste0("nmf_", colnames(Dlist_subset))
  
  # Subset the list of cluster memberships
  clus_memb <- nmf_cluster_list[[grep(Reference, names(nmf_cluster_list))]]
  
  # Get correlation matrix of nmf and kmeans results
  cor_matrix <- cor(Dlist_subset, kmeans_dscore)
  
  # Melt the correlation matrix
  cor_matrix_melted <- melt(cor_matrix)
  cor_matrix_melted <- cor_matrix_melted[order(cor_matrix_melted$value,
                                               decreasing = TRUE), ]
  
  # Get cluster specific matrices
  cor_matrix_k2 <- cor_matrix_melted[grep("K2", cor_matrix_melted$Var1), ]
  cor_matrix_k2 <- cor_matrix_k2[grep("K2", cor_matrix_k2$Var2), ]
  
  cor_matrix_k3 <- cor_matrix_melted[grep("K3", cor_matrix_melted$Var1), ]
  cor_matrix_k3 <- cor_matrix_k3[grep("K3", cor_matrix_k3$Var2), ]
  
  cor_matrix_k4 <- cor_matrix_melted[grep("K4", cor_matrix_melted$Var1), ]
  cor_matrix_k4 <- cor_matrix_k4[grep("K4", cor_matrix_k4$Var2), ]
  
  # Make Maps
  k2_maps <- c()
  for (row in 1:2) {
    nmf_mapping <- unlist(strsplit(as.character(cor_matrix_k2$Var1)[row], "_"))
    nmf_mapping <- nmf_mapping[length(nmf_mapping)]
    
    k_mapping <- unlist(strsplit(as.character(cor_matrix_k2$Var2)[row], "_"))
    k_mapping <- k_mapping[length(k_mapping)]
    k2_maps <- rbind(k2_maps, c(as.numeric(paste(nmf_mapping)),
                                as.numeric(paste(k_mapping))))
  }
  
  k3_maps <- c()
  for (row in 1:3) {
    nmf_mapping <- unlist(strsplit(as.character(cor_matrix_k3$Var1)[row], "_"))
    nmf_mapping <- nmf_mapping[length(nmf_mapping)]
    
    k_mapping <- unlist(strsplit(as.character(cor_matrix_k3$Var2)[row], "_"))
    k_mapping <- k_mapping[length(k_mapping)]
    k3_maps <- rbind(k3_maps, c(as.numeric(paste(nmf_mapping)),
                                as.numeric(paste(k_mapping))))
  }
  
  # We only care about the top three highest correlations
  k4_maps <- c()
  for (row in 1:4) {
    nmf_mapping <- unlist(strsplit(as.character(cor_matrix_k4$Var1)[row], "_"))
    nmf_mapping <- nmf_mapping[length(nmf_mapping)]
    
    k_mapping <- unlist(strsplit(as.character(cor_matrix_k4$Var2)[row], "_"))
    k_mapping <- k_mapping[length(k_mapping)]
    k4_maps <- rbind(k4_maps, c(as.numeric(paste(nmf_mapping)),
                                as.numeric(paste(k_mapping))))
  }
  
  # Map labels
  for (assgn in 1:nrow(k2_maps)) {
    map <- k2_maps[assgn, 1]
    place_ref <- grep(paste0("^", map, "$"),
                      clus_memb$ClusterK2)
    return_tmp <- paste(k2_maps[assgn, 2], "tmp", sep = "_")
    clus_memb$ClusterK2[place_ref] <- return_tmp
  }
  
  for (assgn in 1:nrow(k3_maps)) {
    map <- k3_maps[assgn, 1]
    place_ref <- grep(paste0("^", map, "$"),
                      clus_memb$ClusterK3)
    return_tmp <- paste(k3_maps[assgn, 2], "tmp", sep = "_")
    clus_memb$ClusterK3[place_ref] <- return_tmp
  
  }
  
  for (assgn in 1:nrow(k4_maps)) {
    map <- k4_maps[assgn, 1]
    place_ref <- grep(paste0("^", map, "$"),
                      clus_memb$ClusterK4)
    return_tmp <- paste(k4_maps[assgn, 2], "tmp", sep = "_")
    clus_memb$ClusterK4[place_ref] <- return_tmp
    
  }

  # Return all cluster assignments to numeric values
  clus_memb[, 1] <- as.integer(unlist(strsplit(clus_memb[, 1], "_tmp")))
  clus_memb[, 2] <- as.integer(unlist(strsplit(clus_memb[, 2], "_tmp")))
  clus_memb[, 3] <- as.integer(unlist(strsplit(clus_memb[, 3], "_tmp")))
  
  # Reassign the clusters and return
  nmf_cluster_list[[grep(Reference, names(nmf_cluster_list))]] <- clus_memb
  return(nmf_cluster_list)
}


MapClusters <- function (DistMatrixList, dataset_names, Reference = "TCGA") {
  # ~~~~~~~~~~~~~~
  # This function will use a distance matrix and hierarchical
  # clustering to map clusters across populations.
  #
  # Args: 
  # DistMatrixList: a list of distance matrices for each value of k
  # dataset_names: a character string of all names in the dataset
  # Reference: defaults to TCGA; which population to use for mapping clusters
  # across populations
  #
  # Returns:
  # New cluster assignments for the reference dataset
  # ~~~~~~~~~~~~~~
  
  # Curate the dataset names
  dataset_names <- gsub("_eset", "", dataset_names)
  dataset_names <- dataset_names[!grepl(Reference, dataset_names)]
  
  # This is initializing the list that we will be returning
  resultList <- list()
        
  # Loop over each number of centroids
  for (centroid in substr(names(DistMatrixList), 2, 2)) {
    
    # Add to the result list
    resultList[[paste0("K", centroid)]] <- list()
    
    # Convert the correlation matrix to a dist object
    thisDistance <- DistMatrixList[[paste0("K", centroid)]]
    
    # Subset to only consider correlations with the reference dataset
    correlation_subset <- thisDistance[grep(Reference, thisDistance[, 1]), ]
    
    # Remove Reference compared to Reference
    correlation_subset <- correlation_subset[!grepl(Reference,
                                                    correlation_subset[, 2]), ]
    
    # Now, compare the already mapped reference clusters one by
    # one to each dataset
    reference_clus <- as.character(unique(correlation_subset[, 1]))
    
    # Make sure the reference clusters are in order
    reference_clus <- reference_clus[order(reference_clus)]
    
    for (clus in 1:length(reference_clus)) {
      # Keep adding to resultList
      resultList[[paste0("K", centroid)]][[clus]] <- c(reference_clus[clus])
      
      clus_subset <- correlation_subset[
        correlation_subset[, 1] == reference_clus[clus], ]
      for (dataset in dataset_names) {
        clus_cor <- clus_subset[grep(dataset, clus_subset[, 2]), ]
        other_clus <- as.character(unique(clus_cor[, 2]))
        max_cor <- -1
        for (cur_clus in other_clus) {
          cur_val <- clus_cor$value[clus_cor[, 2] == cur_clus]
          if (cur_val > max_cor) {
            max_cor <- cur_val
            max_clus <- cur_clus
          }
        }
        
        # After this logic you have the cluster with the highest correlation
        # Add it to the result list
        resultList[[paste0("K", centroid)]][[clus]] <- c(resultList[[
          paste0("K", centroid)]][[clus]], max_clus)
        # Remove it from correlation_subset to never be considered again
        correlation_subset <- correlation_subset[
          !grepl(max_clus, correlation_subset[, 2]), ]
      }
    }
  }
  return(resultList)
}


runNMF <- function (Data, k, fname, KClusterAssign, nruns = 10, coph = FALSE,
                    coph_range = 2:8) {
  # ~~~~~~~~~~~~~~
  # This function will perform NMF clustering on gene expression data
  #
  # Args: 
  # Data: Sample by gene expression dataframe
  # k: The number of clusters NMF should find
  # fname: The output file name to store result
  # KClusterAssign: The k-means cluster assignments
  # nruns: How many iterations NMF should perform
  # coph: Decide whether or not to output cophentic coefficient plots
  # coph_range: for what values of k is the cophenetic coefficient determined?
  #
  # Returns:
  # 1) Consensus matrices with silhouette and k-means tracks
  # 2) Cluster membership files for each sample
  # 3) Cophenetic coefficient plots
  # ~~~~~~~~~~~~~~
  
  # Determine the minimum of the data
  minimum <- min(Data)
  
  # If the minimum is less than zero, renormalize because NMF cannot handle
  # negative values
  if (minimum < 0) {
    # Renormalize removing negative values
    const <- -1 * minimum
    Data <- apply(Data, 2, function(x) {x + const})
  }
  
  if (!coph) {
    # Perform nmf clustering
    clus.nmf <- nmf(Data, k, nrun = nruns, .options = "t")
    
    # Get Clustermembership predictions
    clusterMemb <- predict(clus.nmf)
    
    # Write the consensus mapping as a figure
    fpath <- file.path("2.Clustering_DiffExprs", "Figures", "nmf",
                       "ConsensusMaps", paste0(fname, ".png"))
    png(fpath, width = 700, height = 570)
    
    # Show a plot of k = 4 k means clusters overlayed with the NMF clusters
    consensusmap(clus.nmf, labCol = NA, labRow = NA, tracks = c("silhouette:"), 
                 annCol = list("kmeans" =
                                 as.character(paste(KClusterAssign$ClusterK4))),
                 annColors = list("kmeans" =
                                    c("skyblue1", "tomato",
                                      "springgreen", "violet")),
                 main = "",
                 fontsize = 20, treeheight = 60)
    dev.off()
    
    return(as.character(clusterMemb))
    
  } else {
    estim.coeff <- nmf(Data, coph_range, nrun = nruns, seed = 123456)
    return(estim.coeff)
  }
}


# Organize confidence limits table for pearson correlations
organize_confidence <- function(confidence_row) {
  # ~~~~~~~~~~~~~~
  # This function will return a character string with information about confidence
  # correlations
  #
  # Args: 
  # confidence_row - a row of the confidence dataframe that stores correlation info
  #
  # Returns: a list with two elements
  # 1) The character string of the entry
  # 2) The row and column of the matrix belonging to the entry
  # ~~~~~~~~~~~~~~
  
  return_list <- list()
  name <- confidence_row[1]
  r <- sprintf("%.2f", round(as.numeric(paste(confidence_row[2])), 2))
  up <- sprintf("%.2f", round(as.numeric(paste(confidence_row[3])), 2))
  down <- sprintf("%.2f", round(as.numeric(paste(confidence_row[4])), 2))
  entry_output <- paste0(r, " (", up, " to ", down, ")")
  
  # Split the name
  first_comp <- unlist(strsplit(name, ":"))[1]
  second_comp <- unlist(strsplit(name, ":"))[2]
  
  # Split each of the clusters
  first_comp_info <- unlist(strsplit(first_comp, "_"))
  second_comp_info <- unlist(strsplit(second_comp, "_"))
  
  # Logic to determine where each entry should go
  first_cluster <- first_comp_info[grepl("Cluster", first_comp_info)]
  second_cluster <- second_comp_info[grepl("Cluster", second_comp_info)]
  
  if (first_cluster != second_cluster) {
    entry_result <- NA
    return(return_list)
  } else {
    k <- substr(first_cluster, 9, 10)
  }
  
  # Get the dataset
  first_dataset <- first_comp_info[1]
  second_dataset <- second_comp_info[1]
  
  # Get the subtype
  first_subtype <- first_comp_info[length(first_comp_info)]
  second_subtype <- second_comp_info[length(second_comp_info)]
  
  # Get dataset entry
  if (first_dataset == second_dataset) {
    entry_result <- NA
  } else {
    entry_result <- dataset_entry(first_dataset, first_subtype, second_dataset,
                                  second_subtype, k)
  }
  
  return_list <- list("conf" = entry_output, "entry" = entry_result)
  
  return(return_list)
}


# Determine Dataset entry 
dataset_entry <- function(dataset_one, one_subtype,
                          dataset_two, two_subtype, clus) {
  # ~~~~~~~~~~~~~~
  # This function will return a number depending on the dataset_name input
  #
  # Args: 
  # dataset_one - the name of the first dataset to consider
  # one_subtype - the subtype of the first dataset
  # dataset_two - the name of the second dataset to consider
  # two_subtype - the subtype of the second dataset
  # clus - the k that is currently being compared
  #
  # Returns:
  # two numbers indicating the row and column of interest
  # ~~~~~~~~~~~~~~
 
  if (dataset_one == "TCGA" |
      dataset_two == "GSE9891" |
      dataset_one == "GSE26712" |
      dataset_two == "GSE26712") {
    return(NA)
  } else if (dataset_one == "Mayo" && dataset_two == "GSE32062.GPL6480") {
    return(NA)
  }
  
  
  if (clus == "2") {
    if (dataset_one == "Mayo") {
      if (one_subtype == "1") {
        row <- 1
      } else {
        row <- 2
      }
    }
    
    if (dataset_one == "GSE32062.GPL6480") {
      if (one_subtype == "1") {
        row <- 3
      } else {
        row <- 4
      }
    }
    
    if (dataset_one == "GSE9891") {
      if (one_subtype == "1") {
        row <- 5
      } else {
        row <- 6
      }
    }
    
    if (dataset_two == "TCGA") {
      if (two_subtype == "1") {
        column <- 1
      } else {
        column <- 2
      }
    }
      
    if (dataset_two == "Mayo") {
      if (two_subtype == "1") {
        column <- 3
      } else {
        column <- 4
      }
    }
      
    if (dataset_two == "GSE32062.GPL6480") {
      if (two_subtype == "1") {
        column <- 5
      } else {
        column <- 6
      }
    }
  } else if (clus == "3") {
    
    if (dataset_one == "Mayo") {
      if (one_subtype == "1") {
        row <- 7
      } else if (one_subtype == "2") {
        row <- 8
      } else {
        row <- 9
      }
    }
    
    if (dataset_one == "GSE32062.GPL6480") {
      if (one_subtype == "1") {
        row <- 10
      } else if (one_subtype == "2") {
        row <- 11
      } else {
        row <- 12
      }
    }
    
    if (dataset_one == "GSE9891") {
      if (one_subtype == "1") {
        row <- 13
      } else if (one_subtype == "2") {
        row <- 14
      } else {
        row <- 15
      }
    }
    
    if (dataset_two == "TCGA") {
      if (two_subtype == "1") {
        column <- 7
      } else if (two_subtype == "2"){
        column <- 8
      } else {
        column <- 9
      }
    }
    
    if (dataset_two == "Mayo") {
      if (two_subtype == "1") {
        column <- 10
      } else if (two_subtype == "2") {
        column <- 11
      } else {
        column <- 12
      }
    }
    
    if (dataset_two == "GSE32062.GPL6480") {
      if (two_subtype == "1") {
        column <- 13
      } else if (two_subtype == "2") {
        column <- 14
      } else {
        column <- 15
      }
    }
  } else {
    
    if (dataset_one == "Mayo") {
      if (one_subtype == "1") {
        row <- 16
      } else if (one_subtype == "2") {
        row <- 17
      } else if (one_subtype == "3") {
        row <- 18
      } else {
        row <- 19
      }
    }
    
    if (dataset_one == "GSE32062.GPL6480") {
      if (one_subtype == "1") {
        row <- 20
      } else if (one_subtype == "2") {
        row <- 21
      } else if (one_subtype == "3") {
        row <- 22
      } else {
        row <- 23
      }
    }
    
    if (dataset_one == "GSE9891") {
      if (one_subtype == "1") {
        row <- 24
      } else if (one_subtype == "2") {
        row <- 25
      } else if (one_subtype == "3") {
        row <- 26
      } else {
        row <- 27
      }
    }
    
    if (dataset_two == "TCGA") {
      if (two_subtype == "1") {
        column <- 16
      } else if (two_subtype == "2"){
        column <- 17
      } else if (two_subtype == "3") {
        column <- 18
      } else {
        column <- 19
      }
    }
    
    if (dataset_two == "Mayo") {
      if (two_subtype == "1") {
        column <- 20
      } else if (two_subtype == "2") {
        column <- 21
      } else if (two_subtype == "3") {
        column <- 22
      } else {
        column <- 23
      }
    }
    
    if (dataset_two == "GSE32062.GPL6480") {
      if (two_subtype == "1") {
        column <- 24
      } else if (two_subtype == "2") {
        column <- 25
      } else if (two_subtype == "3") {
        column <- 26
      } else {
        column <- 27
      }
    }
  }

  return_list <- list("row" = row, "column" = column)
  return(return_list)
}

plot_reassigned_heatmaps <- function(shuffle, bNMF, Dlist.mapped.cor) {
    if (shuffle) {
      shuffle.string <- "shuffle"
    } else {
      shuffle.string <- ""
    }
    for (plot in 1:length(Dlist.mapped.cor)) {
      
      num_clus <- krange[plot]
      # Build indices for plot comparisions
      all_centroid_plots <- list()
      data_iter <- 1
      for (i in datasets[1:(length(datasets) - 1)]) {
        other_comparisons <- Dlist.mapped.cor[[plot]][
          grepl(i, Dlist.mapped.cor[[plot]][, 1]), ]
        for (j in datasets[2:length(datasets)]) {
          builder <- paste(i, j, sep = "-")
          print(builder)
          
          final_comparison <- other_comparisons[
            grepl(j, other_comparisons[, 2]), ]
          final_comparison <- final_comparison[1:(nrow(final_comparison) / 2), ]
          
          if (data_iter == 1) {
            blank_glob <- ggplot(data = data.frame(final_comparison),
                                 aes(x = Var1, y = Var2,
                                     fill = as.numeric(paste(value)))) +
              geom_blank() + xlab("") + ylab("") + 
              theme(axis.line = element_blank(),
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    plot.background = element_blank(), 
                    legend.position = "none")
          }
          
          g <- ggplot(data = data.frame(final_comparison),
                      aes(x = Var1, y = Var2, fill = as.numeric(paste(value))))

          x.label <- strsplit(toString(final_comparison$Var1[1]),
                              split = "_Cluster")[[1]][1]
          y.label <- strsplit(toString(final_comparison$Var2[1]),
                              split = "_Cluster")[[1]][1]
          
          g <- g + 
            geom_tile(color = "white") + 
            scale_fill_gradient2(high = "red", low = "blue", mid = "white",
                                 midpoint = 0, limit = c(-1, 1)) + 
            xlab(x.label) + 
            ylab(y.label) + 
            scale_x_discrete(labels = paste(1:num_clus)) + 
            scale_y_discrete(labels = paste(1:num_clus)) +
            theme(axis.line = element_blank(),
                  axis.ticks = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  plot.background = element_blank(), 
                  legend.position = "none")
        
          all_centroid_plots[[data_iter]] <- g
          data_iter <- data_iter + 1
          
        }
      }
      fname <- paste0(shuffle.string, "k", num_clus, ".pdf")
      pdf(file.path("2.Clustering_DiffExprs", "Figures", fname))
      for (p in 1:length(all_centroid_plots)) {
        plot(all_centroid_plots[[p]])
      }
      dev.off()
    }
  return(NULL)
}
