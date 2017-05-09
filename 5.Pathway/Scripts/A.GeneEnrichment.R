############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# The script will output a table of the differentially expressed genes for each cluster, 
# both positively and negatively regulated

suppressMessages(library(checkpoint))
suppressMessages(checkpoint('2016-03-01', checkpointLocation = "."))

args <- commandArgs(trailingOnly=TRUE)
#args <- c("TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset")

sig_alpha <- 0.05 / 10930 # 0.05 / number of common genes

############################################
# Load Libraries
############################################
library(limma)
library(plyr)

############################################
# Load Data
############################################
# Get moderated t scores and p values
datasetMembers <- list.files(path = "2.Clustering_DiffExprs/Tables/SAM_results/")

# SAMlist is a list of lists holding SAM results for each dataset for each cluster 
# and for each centroid assignment
SAMList <- list()
for (dataset in 1:length(args)) {
  holder <- args[dataset]
  DMember <- datasetMembers[grep(holder, datasetMembers)]
  tmpList <- list()
  
  for(cluster in 1:length(DMember)) {
    tmpList[[DMember[cluster]]] <- read.table(paste("2.Clustering_DiffExprs/Tables/SAM_results/", 
                                                    DMember[cluster], sep = ""), sep = ",", 
                                              row.names = 1, header = T)
  }
  
  SAMList[[args[dataset]]] <- tmpList
}

############################################
# Identify cluster genes
############################################
# Select the genes that meet the alpha FDR threshhold
ExplanatoryGenesList <- list()
for (dataset in 1:length(SAMList)) {
  DatasetGenesList <- list()
  for (centroid in 1:length(SAMList[[dataset]])) {
    
    # Extract Data
    useData <- SAMList[[dataset]][[centroid]]
    
    # Get Data name and cluster name
    dataname <- unlist(strsplit(names(SAMList[[dataset]])[centroid], "Stat_"))[2]
    dataname_cluster <- unlist(strsplit(paste("Cluster", unlist(strsplit(dataname, 
                                                                         "Cluster"))[2], 
                                              sep = ""), "[.]"))[1] 
    tmp_ExpGene_list <- list()
    for(cluster in 1:(ncol(useData) / 2)) { # Divide by two because each cluster has two columns

      # We are interested in cluster specific moderated t scores and p values
      useCols <- ((cluster * 2) - 1):(cluster * 2)
      SAMdataSub <- useData[ ,useCols]
      
      # Get the genes that meet the cutoff
      explanatoryGenes <- SAMdataSub[SAMdataSub[ , 2] < sig_alpha, ]
      
      # Store these genes into the temporary explanatory genes list
      tmp_ExpGene_list[[paste(dataname_cluster, cluster, sep = "_")]] <- explanatoryGenes
    }
    DatasetGenesList[[dataname_cluster]] <- tmp_ExpGene_list
  }
  ExplanatoryGenesList[[args[dataset]]] <- DatasetGenesList
}

############################################
# Compare cluster genes
############################################
# The second filtering step is to observe cluster specific overlaps 
# across the 4 main datasets separately for k = 2, k = 3, and k = 4
subtypeOverlapingList <- centroidOverlapingList <- list()

#Loop over the different numbers of centroids in the "TCGA_eset" gene list to begin
for (centroid in 1:length(ExplanatoryGenesList[["TCGA_eset"]])) {
  
  # Loop over the different numbers of subtypes
  for(subtype in 1:length(ExplanatoryGenesList[["TCGA_eset"]][[centroid]])) {
    
    # Loop over the different datasets to observe centroid and subtype specific overlapping genes 
    for(dataset in 1:length(ExplanatoryGenesList)) {
      
      # Get the dataframe that holds the significantly differentially expressed genes with p values and d stats
      ClusterGenes <- ExplanatoryGenesList[[dataset]][[centroid]][[subtype]]
      
      # Get the genes with d statistics greater than 0, indicating up regulated genes
      ClusterGenesPos <- rownames(ClusterGenes[ClusterGenes[ ,1] > 0, ])
      
      # Get the genes with d statistics less than 0, indicating down regulated genes
      ClusterGenesNeg <- rownames(ClusterGenes[ClusterGenes[ ,1] < 0, ])
      
      # initialize an intersecting character string of genes
      if (dataset == 1) {
        intersectingPos <- ClusterGenesPos
        intersectingNeg <- ClusterGenesNeg
        
        # Observe the cluster and centroid specific overlaps
      } else {
        intersectingPos <- intersect(intersectingPos, ClusterGenesPos)
        intersectingNeg <- intersect(intersectingNeg, ClusterGenesNeg)
      }
    }
    
    # Store subtpye specific overlapping genes
    listname <- paste('Cluster', subtype, sep = "")
    subtypeOverlapingList[[listname]] <- list("Pos" = intersectingPos, "Neg" = intersectingNeg)
  }
  
  # Store centroid specific overlapping genes
  centroidname <- paste("K", centroid + 2, sep = "")
  centroidOverlapingList[[centroidname]] <- subtypeOverlapingList
}

############################################
# Get data in format to write to file
############################################
SigGeneMatrix <- matrix("", nrow = 800, ncol = 18)
set_columns <- list(c(1,2), c(3, 4), c(5, 6), c(7, 8))
for (centroid in 1:length(centroidOverlapingList)) {
  
  # Get Centroid Specific Data
  CentroidSpecific <- centroidOverlapingList[[centroid]]
  
  data_matrix <- matrix("", nrow = 800, ncol = length(CentroidSpecific) * 2)
  for (clus in 1:length(CentroidSpecific)) {
    
    # Get Cluster Specific Data
    ClusterSpecific <- CentroidSpecific[[clus]]
    
    # Initialize an empty matrix to store cluster specific pos and neg genes
    tmpMatrix <- matrix("", nrow = 800, ncol = 2) #800 is bigger than any gene list
    
    # Combine cluster specific pos and neg genes into single matrix
    for (regdir in 1:length(ClusterSpecific)) {
      exprsSpecific <- ClusterSpecific[[regdir]]
      tmpMatrix[1:length(exprsSpecific), regdir] <- exprsSpecific
    }
    data_matrix[, c(set_columns[[clus]])] <- tmpMatrix
  }
  
  # Append these matrices into a full summary of pos and neg genes across clusters
  if (length(CentroidSpecific) == 2) {
    SigGeneMatrix[ ,1:4] <- data_matrix
  } else if(length(CentroidSpecific) == 3)  {
    SigGeneMatrix[ ,5:10] <- data_matrix
  } else {
    SigGeneMatrix[ ,11:18] <- data_matrix
  }
}

colnames(SigGeneMatrix) <- c("K2Cluster1_Pos", "K2Cluster1_Neg", "K2Cluster2_Pos", "K2Cluster2_Neg",
                             "K3Cluster1_Pos", "K3Cluster1_Neg", "K3Cluster2_Pos", "K3Cluster2_Neg", 
                             "K3Cluster3_Pos", "K3Cluster3_Neg", "K4Cluster1_Pos", "K4Cluster1_Neg", 
                             "K4Cluster2_Pos", "K4Cluster2_Neg", "K4Cluster3_Pos", "K4Cluster3_Neg", 
                             "K4Cluster4_Pos", "K4Cluster4_Neg")

# Write supplementary table to file
write.table(SigGeneMatrix, file = "5.Pathway/Tables/Overlappling_Signinficant_GeneLists.csv", 
            sep = ",", row.names = F)

############################################
# Write Cluster Specific Genes to File
############################################
k2_one_genes <- c(SigGeneMatrix[ , 1], SigGeneMatrix[ , 2])
k2_one_genes <- k2_one_genes[k2_one_genes != ""]
k2_one_genes <- as.matrix(k2_one_genes)
colnames(k2_one_genes) <- "Gene symbol"
write.table(k2_one_genes, "5.Pathway/Gene_lists/K2_Cluster1_genes.csv", sep = ",", row.names = F, quote = F)

k2_two_genes <- c(SigGeneMatrix[ , 3], SigGeneMatrix[ , 4])
k2_two_genes <- k2_two_genes[k2_two_genes != ""]
k2_two_genes <- as.matrix(k2_two_genes)
colnames(k2_two_genes) <- "Gene symbol"
write.table(k2_two_genes, "5.Pathway/Gene_lists/K2_Cluster2_genes.csv", sep = ",", row.names = F, quote = F)

k3_one_genes <- c(SigGeneMatrix[ , 5], SigGeneMatrix[ , 6])
k3_one_genes <- k3_one_genes[k3_one_genes != ""]
k3_one_genes <- as.matrix(k3_one_genes)
colnames(k3_one_genes) <- "Gene symbol"
write.table(k3_one_genes, "5.Pathway/Gene_lists/K3_Cluster1_genes.csv", sep = ",", row.names = F, quote = F)

k3_two_genes <- c(SigGeneMatrix[ , 7], SigGeneMatrix[ , 8])
k3_two_genes <- k3_two_genes[k3_two_genes != ""]
k3_two_genes <- as.matrix(k3_two_genes)
colnames(k3_two_genes) <- "Gene symbol"
write.table(k3_two_genes, "5.Pathway/Gene_lists/K3_Cluster2_genes.csv", sep = ",", row.names = F, quote = F)

k3_three_genes <- c(SigGeneMatrix[ , 9], SigGeneMatrix[ , 10])
k3_three_genes <- k3_three_genes[k3_three_genes != ""]
k3_three_genes <- as.matrix(k3_three_genes)
colnames(k3_three_genes) <- "Gene symbol"
write.table(k3_three_genes, "5.Pathway/Gene_lists/K3_Cluster3_genes.csv", sep = ",", row.names = F, quote = F)

k4_one_genes <- c(SigGeneMatrix[ , 11], SigGeneMatrix[ , 12])
k4_one_genes <- k4_one_genes[k4_one_genes != ""]
k4_one_genes <- as.matrix(k4_one_genes)
colnames(k4_one_genes) <- "Gene symbol"
write.table(k4_one_genes, "5.Pathway/Gene_lists/K4_Cluster1_genes.csv", sep = ",", row.names = F, quote = F)

k4_two_genes <- c(SigGeneMatrix[ , 13], SigGeneMatrix[ , 14])
k4_two_genes <- k4_two_genes[k4_two_genes != ""]
k4_two_genes <- as.matrix(k4_two_genes)
colnames(k4_two_genes) <- "Gene symbol"
write.table(k4_two_genes, "5.Pathway/Gene_lists/K4_Cluster2_genes.csv", sep = ",", row.names = F, quote = F)

k4_three_genes <- c(SigGeneMatrix[ , 15], SigGeneMatrix[ , 16])
k4_three_genes <- k4_three_genes[k4_three_genes != ""]
k4_three_genes <- as.matrix(k4_three_genes)
colnames(k4_three_genes) <- "Gene symbol"
write.table(k4_three_genes, "5.Pathway/Gene_lists/K4_Cluster3_genes.csv", sep = ",", row.names = F, quote = F)

k4_four_genes <- c(SigGeneMatrix[ , 17], SigGeneMatrix[ , 18])
k4_four_genes <- k4_four_genes[k4_four_genes != ""]
k4_four_genes <- as.matrix(k4_four_genes)
colnames(k4_four_genes) <- "Gene symbol"
write.table(k4_four_genes, "5.Pathway/Gene_lists/K4_Cluster4_genes.csv", sep = ",", row.names = F, quote = F)

commongenes <- read.csv("1.DataInclusion/Data/Genes/CommonGenes_genelist.csv", sep = ",")
commongenes <- as.matrix(commongenes)
colnames(commongenes) <- "Gene symbol"
write.table(commongenes, "5.Pathway/Gene_lists/Reference.csv", sep = ",", row.names = F, quote = F)
