############################################
# Cross-population analysis of high-grade serous ovarian cancer does not
# support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will run k-means and SAM on the datasets and map the clusters 
# across k and across datasets, it will output cluster membership files and 
# correlations

suppressMessages(library(checkpoint))
suppressMessages(checkpoint("2016-03-01", checkpointLocation = "."))

args <- commandArgs(trailingOnly = TRUE)

args <- c(2, 4, 20, 123, FALSE, FALSE, "madgenes", "TCGA_eset", "mayo.eset",
          "GSE32062.GPL6480_eset", "GSE9891_eset", "GSE26712_eset",
          "aaces.eset")
############################################
# Load Libraries
############################################
library(cluster)
library(curatedOvarianData)
library(ggplot2)
library(RColorBrewer)
library(siggenes)
library(GGally)
library(reshape2)
library(gridExtra)
library(grid)
# The script holds custom functions to run kmeans and SAM
source("2.Clustering_DiffExprs/Scripts/Functions/kmeans_SAM_functions.R")
source("2.Clustering_DiffExprs/Scripts/Functions/heatmap3.R")

# The script loads the ovarian cancer datasets
source("1.DataInclusion/Scripts/Functions/LoadOVCA_Data.R")

############################################
# Constants
############################################
k <- as.numeric(paste(args[1]))
k2 <- as.numeric(paste(args[2]))
kStarts <- as.numeric(paste(args[3]))
kSEED <- as.numeric(paste(args[4]))
bNMF <- as.logical(args[5])
shuffle <- as.logical(args[6])
SAM_subset <- args[7]

# Separate the eset arguments from the rest of the commandArgs
argsCurated <- args[grep("eset", args)]

# Set the seed
set.seed(kSEED)

############################################
# Load Data
############################################
# Add Mayo to argsCurated
if ("Mayo" %in% args) {
  argsCurated <- c(argsCurated[1], "Mayo",
                  argsCurated[2:length(argsCurated)])
}

# Use the LoadOVCA_Data function to read in the datasets subset by commongenes
ExpData <- LoadOVCA_Data(datasets = argsCurated, genelist_subset = SAM_subset,
                         shuffle = shuffle)

# Read in common genes csv file. The csv was generated as an intersection of
# the genes in TCGA, Yoshihara, Mayo, Tothill, and Bonome
CommonGenes <- read.csv("1.DataInclusion/Data/Genes/CommonGenes_genelist.csv",
                        header = T,
                        stringsAsFactors = F)

# Read in mad genes csv file. The file was generated as an intersection
#of the top 1500 most variably expressed genes
GlobalMAD <-
  read.csv(file = "1.DataInclusion/Data/Genes/GlobalMAD_genelist.csv",
           header = T,
           stringsAsFactors = F)

############################################
# Perform k-means clustering using Global MAD
############################################
# Initialize a Clusters list, which will hold the cluster memberships
# for all datasets. The cluster membership will either be derived from
# the k-means algorithm, or the NMF algorithm.

Clusters <- list()
# If bNMF = TRUE in the arguments, then k-means will not be performed
if (!bNMF) {
  # Loop over each dataset with expression data
  for (dataset in names(ExpData)) {
    # Store the results of the k-means clustering to the Clusters list
    Clusters[[dataset]] <-
      KmeansGlobal(ExpData[[dataset]], GlobalMAD[ ,1], k, k2, starts = kStarts)
  }
  # Extract cluster membership files for NMF results
  } else {
  for (dataset in names(ExpData)) {
    # Load NMF cluster membership files
    NMF_files <-
      list.files("2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/")
    NMF_files <-
      NMF_files[grepl("nmf.csv", NMF_files)]
    file <-
      NMF_files[grepl(dataset, NMF_files)]
    NMFclusterMemb <-
      read.csv(paste("2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/",
                     file, sep = ""), sep = ",", row.names = 1, header = T)
    colnames(NMFclusterMemb) <- c("ClusterK2", "ClusterK3", "ClusterK4")
    Clusters[[dataset]] <- NMFclusterMemb
  }
}

############################################
# Run preliminary SAM function to get moderated t score vectors 
############################################
SamList <- RunSam(ExpData, Clusters, d.stat)

############################################
# Combine d statistics
############################################
# The following section will map clusters to each other using Pearson 
# correlationsof moderated t score vectors generated from each gene's moderated
# t statistic. A moderated t statistic is a measure of mean gene expression
# difference weighted by variance. Note - t statistic is equal to d statistic
# from sam function.

# Each element of the Dlist is a dataframe with genes as rows and each
# subtype (3+4) specific d-score as the columns. The naming convention for the
# columns is DatasetName_ClusterK?_SAMf_?? where ? is the k value and ?? is the
# Cluster ID.

Dlist <- list()
krange <- seq(k, k2)
for (dataset in argsCurated) {
  ptrn <- unlist(strsplit(dataset, "_"))[1]

  # SamList holds lists of lists depending on the krange
  whichList <- grep(ptrn, names(SamList))

  # Loop over all the elements in each SAM list to get a list of dataset
  # specific dataframes
  tmpSAMdata <- c()
  for (elem in 1:length(whichList)) {
    listIndex <- whichList[elem]
    # Find the d score vector indices
    indeces <- grep("D", colnames(SamList[[listIndex]]))
    tmpSam <- SamList[[listIndex]][indeces]
    colnames(tmpSam) <-
      paste(names(SamList)[listIndex], colnames(tmpSam), sep = "_")
    if (elem == 1) {
      tmpSAMdata <- tmpSam
    } else {
      tmpSAMdata <- cbind(tmpSAMdata, tmpSam)
    }
  }
  Dlist[[dataset]] <- tmpSAMdata
}

############################################
# Correlate moderated t scores
############################################
# Obtain within dataset, across kmeans correlations
WithinDatasetCor <- list()
for (dataset in 1:length(Dlist)) {
  # Pearson Correlations
  tmpCor <- cor(Dlist[[dataset]], use = "p")
  WithinDatasetCor[[names(Dlist)[dataset]]] <- tmpCor
}

# Each element in Dlist.cor is a data.frame containing all pairs of clusters 
# (k=3 and k=4) and the correlation of the d statistic. The "value" column
# contains the correlation coefficient and the "Var1" and "Var2" columns
# contain cluster IDs
Dlist.cor <- list()
for (centroids in krange) {
  dataUse <- c()
  ptrn <- paste("K", centroids, sep = "")

  # Loop over the entire Dlist holding sam results for all datasets
  for (samResults in 1:length(Dlist)) {
    # Identify which centroid to compare and subset
    whichCentroid <- grep(ptrn, names(Dlist[[samResults]]))
    data_to_cor <- Dlist[[samResults]][whichCentroid]

    tmpdataUse <- c()
    for (samResults_compare in 1:length(Dlist)) {
      # Prepare a similar subset to compare and perform all correlations
      whichCentroid <- grep(ptrn, names(Dlist[[samResults_compare]]))
      data_to_cor_comp <- Dlist[[samResults_compare]][whichCentroid]

      # Ensure that there are overlapping genes and subsetting will execute
      genes <- intersect(rownames(data_to_cor), rownames(data_to_cor_comp))
      data_to_cor <- data_to_cor[genes, ]
      data_to_cor_comp <- data_to_cor_comp[genes, ]

      # Obtain correlations across cluster t scores
      corData <- cbind(data_to_cor, data_to_cor_comp)
      tmpcor <- cor(corData, use = "p")

      # Remove upper triangle of correlation matrix
      tmpcor[upper.tri(tmpcor, diag = TRUE)] <- NA
      # Melt the correlation matrix for a vector by vector correlation
      tmpcor <- melt(tmpcor)
      colnames(tmpcor) <- c("Var1", "Var2", "value")

      # Remove the results in the upper triangle
      tmpcor <- tmpcor[!is.na(tmpcor$value), ]
      tmpdataUse <- rbind(tmpdataUse, tmpcor)
    }
    dataUse <- rbind(dataUse, tmpdataUse)
  }
  Dlist.cor[[paste("K", centroids, sep = "")]] <- dataUse
}

############################################
# Map Clusters Accross Datasets
############################################
# Assign a reference category for all other datasets to map to
# Assign nmf clusters to kmeans clusters
if (bNMF) {
  NewClusters <-
    AssignReference_NMF(kmeans_dscore_dir =
                          "2.Clustering_DiffExprs/Tables/DScores/",
                        nmf_Dlist = Dlist, nmf_cluster_list = Clusters,
                        Reference = "TCGA")
} else {
  NewClusters <- AssignReference("TCGA", Cluster = "ClusterK3",
                                 Cor = WithinDatasetCor,
                                 ClusterList = Clusters)
}

# Run SAM again but with TCGA clusters assigned as a reference
ReferenceSamList <- RunSam(ExpData, NewClusters, d.stat)

# NewDlist will house dataframes of SAM moderated t score vectors
# for each dataset passed as an argument after remapping to reference clusters.
NewDlist <- list()
for (dataset in argsCurated) {
  ptrn <- unlist(strsplit(dataset, "_"))[1]

  # ReferenceSamList holds lists of lists depending on the krange
  whichList <- grep(ptrn, names(ReferenceSamList))

  # Loop over all the elements in each SAM list to get a list of dataset
  # specific dataframes
  tmpSAMdata <- c()
  for (elem in 1:length(whichList)) {
    listIndex <- whichList[elem]

    #Find the d score vector indices
    indeces <- grep("D", colnames(ReferenceSamList[[listIndex]]))
    tmpSam <- ReferenceSamList[[listIndex]][indeces]
    colnames(tmpSam) <- paste(names(ReferenceSamList)[listIndex],
                              colnames(tmpSam), sep = "_")
    if (elem == 1) {
      tmpSAMdata <- tmpSam
    } else {
      tmpSAMdata <- cbind(tmpSAMdata, tmpSam)
    }
  }
  NewDlist[[dataset]] <- tmpSAMdata
}

# Each element in newDlist.cor is a data.frame containing all pairs of clusters
# (k=3 and k=4) and the correlation of the d statistic. The "value" column
# contains the correlation coefficient and the "Var1" and "Var2" columns
# contain cluster IDs
NewDlist.cor <- list()
for (centroids in krange) {
  dataUse <- c()
  ptrn <- paste("K", centroids, sep = "")

  # Loop over the entire newDlist holding sam results for all datasets
  for (samResults in 1:length(NewDlist)) {
    # Identify which centroid to compare and subset
    whichCentroid <- grep(ptrn, names(NewDlist[[samResults]]))
    data_to_cor <- NewDlist[[samResults]][whichCentroid]

    tmpdataUse <- c()
    for (samResults_compare in 1:length(NewDlist)) {
      # Prepare a similar subset to compare and perform all correlations
      whichCentroid <- grep(ptrn, names(NewDlist[[samResults_compare]]))
      data_to_cor_comp <- NewDlist[[samResults_compare]][whichCentroid]

      # Ensure that there are overlapping genes and subsetting will execute
      genes <- intersect(rownames(data_to_cor), rownames(data_to_cor_comp))
      data_to_cor <- data_to_cor[genes, ]
      data_to_cor_comp <- data_to_cor_comp[genes, ]

      # Obtain correlations across cluster t scores
      corData <- cbind(data_to_cor, data_to_cor_comp)
      tmpcor <- cor(corData, use = "p")

      # Remove upper triangle of correlation matrix
      tmpcor[upper.tri(tmpcor, diag = TRUE)] <- NA

      # Melt the correlation matrix for a vector by vector correlation
      tmpcor <- melt(tmpcor)
      colnames(tmpcor) <- c("Var1", "Var2", "value")

      # Remove the results in the upper triangle
      tmpcor <- tmpcor[!is.na(tmpcor$value), ]
      tmpdataUse <- rbind(tmpdataUse, tmpcor)
    }
    dataUse <- rbind(dataUse, tmpdataUse)
  }
  NewDlist.cor[[paste("K", centroids, sep = "")]] <- dataUse
}
  
# Map the rest of the clusters, after the reference clusters are in place
newClus <-
  MapClusters(NewDlist.cor, dataset_names = names(ExpData), Reference = "TCGA")

# Output a final cluster mappings list
MapList <- list()
for (i in 1:length(newClus)) {
  dframe <- c()
  for (j in 1:length(newClus[[i]])) {
    tmp <- newClus[[i]][[j]]
    dframe <- rbind(dframe, tmp)
  }
  MapList[[i]] <- dframe
}

############################################
# Re-assign cluster labels using the maps
############################################
Clusters.mapped <- list()
for (centroid in 1:length(NewClusters)) {

  # Identify the original cluster assignments for the given datasets
  thisDF <- NewClusters[[centroid]]
  thisDataset <- names(ExpData)[centroid]

  # Duplicate the dataframe of cluster assignments
  outputDF <- thisDF
  for (assgn in 1:length(MapList)) {
    for (clus in 1:nrow(MapList[[assgn]])) {

      # Get the dataset name and orignal cluster ID
      thisColName <- grep(thisDataset, MapList[[assgn]][clus, ])
      thisColName <- MapList[[assgn]][clus, thisColName]
      tmpVec <- unlist(strsplit(thisColName, "_"))
      origClusterID <- tmpVec[length(tmpVec)]
      newClusterID <- clus

      # Reassign the clusters according to the new cluster ID
      outputDF[thisDF[, assgn] == origClusterID, assgn] <- newClusterID

      message <-
        paste("K", krange[assgn], ": Mapping Cluster ", origClusterID, " to ",
                       newClusterID, " in ", thisDataset, "\n", sep = "")
      cat(message)
    }
  }
  Clusters.mapped[[centroid]] <- outputDF
}

names(Clusters.mapped) <- names(NewClusters)

# Write the cluster re-assignments
if (!shuffle) {
  if (!bNMF) {
    for (i in 1:length(Clusters.mapped)) {
      fName <- paste("KMembership_", names(ExpData)[i], ".csv", sep = "")
      write.csv(Clusters.mapped[[i]],
                file = file.path("2.Clustering_DiffExprs", "Tables",
                                 "ClusterMembership", "kmeans", fName),
                row.names = TRUE)
    }
  } else {
    for (i in 1:length(Clusters.mapped)) {
      fName <-
        paste("KMembership_", names(ExpData)[i], "_mapped.csv", sep = "")
      write.csv(Clusters.mapped[[i]],
                file = file.path("2.Clustering_DiffExprs/",
                                 "Tables", "ClusterMembership",
                                 "nmf", fName), row.names = TRUE)
    }
  }
}

############################################
# Re-Run SAM Function to get new moderated t score vectors 
############################################
# Output Delta measurements for the given FDR
FDR <- 0.05 / nrow(CommonGenes)

# Run SAM one last time with all the clusters mapped
if (!shuffle) {
  SamList.mapped <- RunSam(ExpData, Clusters.mapped, d.stat,
                           type = "Delta", kFDR = FDR)
  Deltas <- SamList.mapped[[2]]
  SamList.mapped <- SamList.mapped[[1]]
} else {
  SamList.mapped <- RunSam(ExpData, Clusters.mapped, d.stat)
}

# We are only interested in the SAM results using k-means
if (!shuffle) {
  if (!bNMF) {
    write.table(Deltas, "2.Clustering_DiffExprs/Tables/SAM_Deltas.csv",
                sep = ",")

    # Write the results of the SAM to a folder
    for (i in 1:length(SamList.mapped)) {
      file.name <-
        paste("2.Clustering_DiffExprs/Tables/SAM_results/SAM_pVal-Stat_",
                         names(SamList.mapped)[i], ".csv", sep = "")
      write.table(SamList.mapped[[i]], file = file.name,
                  sep = ",", row.names = T,
                  col.names = NA)
    }
  }
}

############################################
# Combine moderated t score vectors
############################################
# This will map the moderated t score vectors for each dataset with the
# reassigned cluster labels
Dlist.mapped <- list()
for (dtaset in 1:length(argsCurated)) {
  ptrn <- unlist(strsplit(argsCurated[dtaset], "_"))[1]

  # Identify the index of the dataset
  whichList <- grep(ptrn, names(SamList.mapped))
  
  tmpData <- c()
  for (j in 1:length(whichList)) {
    listIndex <- whichList[j]
    
    #Find the colnames that store the D score vectors
    indeces <- grep("D", colnames(SamList.mapped[[listIndex]]))
    SAMtmp <- SamList.mapped[[listIndex]][indeces]
    colnames(SAMtmp) <- paste(names(SamList.mapped)[listIndex],
                              colnames(SAMtmp), sep="_")
    
    if (j == 1) {
      tmpData <- SAMtmp
    } else {
      tmpData <- cbind(tmpData, SAMtmp)
    }
  }
  
  # Write to the disk the dataframe of dataset specific moderated t score
  # vectors and specify if nmf or not
  if (!shuffle) {
    if (!bNMF) {
      write.table(tmpData,
                  file = paste("2.Clustering_DiffExprs/Tables/DScores/", 
                               argsCurated[dtaset], "_kmeans_DScoreVectors.csv",
                               sep = ""),
                  sep = ",", row.names = T, col.names = NA)
    } else {
      write.table(tmpData,
                  file =
                    paste("2.Clustering_DiffExprs/Figures/nmf/DscoreVectors/",
                    argsCurated[dtaset], "_nmf_DScoreVectors.csv", sep = ""), 
                  sep = ",", row.names = T, col.names = NA)
    }
  }

  Dlist.mapped[[argsCurated[dtaset]]] <- tmpData
}

############################################
# Re-Correlate D Scores
############################################
# Obtain the correlated list for the remapped clusters
Dlist.mapped.cor <- list()
confidence_data <- c()
for (centroid in 1:length(krange)) {
  dataUse <- c()
  ptrn <- paste("K", krange[centroid], sep = "")
  for (clus in 1:length(Dlist.mapped)) {
    which <- grep(ptrn, names(Dlist.mapped[[clus]]))
    
    # Extract the data from the mapped list according to the given centroid
    data <- Dlist.mapped[[clus]][which]
    tmpdataUse <- c()
    
    # Loop over the same mapped Dlist to get correlations
    for (clus2 in 1:length(Dlist.mapped)) {
      which <- grep(ptrn, names(Dlist.mapped[[clus2]]))
      data1 <- Dlist.mapped[[clus2]][which]
      
      # To be safe, ensure that the genes exist in both datasets
      genes <- intersect(rownames(data), rownames(data1))
      dataSub <- data[genes, ]
      data1Sub <- data1[genes, ]
      
      # Get ready to observe Pearson correlations
      corData <- cbind(dataSub, data1Sub)
      tmpcor <- cor(corData, use = "p")
      
      # Obtain Confidence Limits for each comparison
      for (vect_A in 1:ncol(corData)) {
        for (vect_B in 1:ncol(corData)) {
          comp_name <- paste0(colnames(corData)[vect_A], ":",
                              colnames(corData)[vect_B])
          clus_A <- corData[, vect_A]
          clus_B <- corData[, vect_B]
          cor_result <- cor.test(clus_A, clus_B)
          confidence_int <- cor_result$conf.int
          output_result <- c(comp_name, cor_result$estimate[[1]],
                             cor_result$conf.int[[1]],
                             cor_result$conf.int[[2]])
          confidence_data <- rbind(confidence_data, output_result)
        }
      }
      
      # Melt the correlation matrix for a vector by vector correlation
      tmpcor <- melt(tmpcor)
      colnames(tmpcor) <- c("Var1", "Var2", "value")
      tmpdataUse <- rbind(tmpdataUse, tmpcor)
    }
    dataUse <- rbind(dataUse, tmpdataUse)
  }
  Dlist.mapped.cor[[paste("K", krange[centroid], sep = "")]] <- dataUse
}

############################################
# Write Correlations to File
############################################
# Within dataset correlations file with mapped clusters
WithinDatasetCor <- list()
for (centroid in 1:length(Dlist.mapped)) {
  tmpCor <- cor(Dlist.mapped[[centroid]], use = "p")
  WithinDatasetCor[[centroid]] <- tmpCor
  names(WithinDatasetCor)[centroid] <- names(Dlist.mapped)[centroid]
  if (!shuffle) {
    if (!bNMF) {
      write.table(WithinDatasetCor[[centroid]], file = 
                    paste("2.Clustering_DiffExprs/Tables/WithinCor/",
                          names(WithinDatasetCor)[centroid],
                          "_", SAM_subset,
                          "_WithinDatasetCorrelations.csv", sep = ""),
                  row.names = T, col.names = NA, sep = ",")
    } else {
      write.table(WithinDatasetCor[[centroid]],
                  file = paste("2.Clustering_DiffExprs/Figures/nmf/WithinCor/",
                               names(WithinDatasetCor)[centroid],
                               "_", SAM_subset,
                               "nmf_WithinDatasetCorrelations.csv", sep = ""),
                  row.names = T, col.names = NA, sep = ",")
    }
  }
}

############################################
# Write Confidence Limits to File
############################################
colnames(confidence_data) <- c("comparison", "value", "lower", "upper")

t1 <-
  unlist(apply(confidence_data, 1,
               function(x) return((strsplit(x["comparison"], ":"))[[1]][1])))
t2 <-
  unlist(apply(confidence_data, 1,
               function(x) return((strsplit(x["comparison"], ":"))[[1]][2])))

confidence.frame <- data.frame(confidence_data) 
confidence.frame$test1 <- t1
confidence.frame$test2 <- t2


############################################
# Recast melted correlation dataframes for across dataset comparisons
############################################

for (centroid in 1:length(Dlist.mapped.cor)) {
  tmpCor <-
    dcast(Dlist.mapped.cor[[centroid]], Var1~Var2, mean, value = value)
  rownames(tmpCor) <- tmpCor[, 1]
  tmpCor <- tmpCor[, -1]
  if (!shuffle) {
    if (!bNMF) {
      write.table(tmpCor, 
                  file =
                    paste(
                      "2.Clustering_DiffExprs/Tables/AcrossCor/AcrossDatasetCor_K",
                      krange[centroid], "_", SAM_subset, ".csv", sep = ""),
                  sep = ",", row.names = T, col.names = NA)
    } else {
      write.table(tmpCor,
                  file =
                    paste(
                      "2.Clustering_DiffExprs/Figures/nmf/AcrossCor/AcrossDatasetCor_nmf_K",
                      krange[centroid], "_", SAM_subset, ".csv", sep = ""),
                  sep = ",", row.names = T, col.names = NA)
    }
  }
}

############################################
# Before plotting figures, remove Bonome dataset
############################################
for (cor.list in 1:length(Dlist.mapped.cor)) {
  for (centroid in 1:(ncol(Dlist.mapped.cor[[cor.list]]) - 1)) {
    Dlist.mapped.cor[[cor.list]] <-
      Dlist.mapped.cor[[cor.list]][-grep(
        "GSE26712_eset", paste(Dlist.mapped.cor[[cor.list]][, centroid])), ]
  }
}

datasets <- argsCurated[-grep("GSE26712_eset", argsCurated)]

############################################
# Plot Re-Assigned Heatmaps
############################################
plot_reassigned_heatmaps(shuffle, bNMF, Dlist.mapped.cor)
