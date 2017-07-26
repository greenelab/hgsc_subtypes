############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will output kaplan-meier curves and survival assessments for all input datasets

args <- commandArgs(trailingOnly=TRUE)

# args <- c("TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset")
############################################
# Load Libraries
############################################
library(curatedOvarianData)
library(survival)

# The script that holds custom functions
source("4.Survival/Scripts/Functions/Survival_Functions.R")

############################################
# Data in curatedOvarianData
############################################
detailedData <- data(package="curatedOvarianData")[3]
detailedData.names <- detailedData$results[,3]
inCuratedOvarianData <- args %in% detailedData.names

argsCurated <- args[inCuratedOvarianData]

############################################
# Load Phenotype Data
############################################
phenoData <- list()

for (phenoset in 1:length(args)) {
  # Load the expression data
  dta <- c()
  if (args[phenoset] %in% detailedData.names) {
    # Get the list of samples to remove
    fName <- paste(args[phenoset], "_samplesRemoved.csv", sep = "")
    good <- read.csv(file = file.path("1.DataInclusion", "Data", "GoodSamples", fName))
    accept <- c(paste(good$x))
    
    #load the phenotype data
    data(list = args[phenoset], package = "curatedOvarianData")
    tmp <- get(args[phenoset])
    dta <- pData(tmp)
    phenoData[[phenoset]] <- dta[accept, ]
    names(phenoData)[phenoset] <- args[phenoset]
    
    # The phenotype data is held separately for the Mayo data
  } else {
    if (args[phenoset] == "Mayo") {
      dta <- pData(get(load("1.DataInclusion/Data/Mayo/MayoEset.Rda")))
      Mayogood <- read.csv(file = file.path("1.DataInclusion", "Data", "GoodSamples", 
                                            "Mayo_samplesRemoved.csv"))
      Mayoaccept <- c(paste(Mayogood$x))
      rownames(dta) <- dta[ ,1]
      phenoData[[phenoset]] <- dta[Mayoaccept, ]
      names(phenoData)[phenoset] <- args[phenoset]
    }
  }
}

############################################
# Get Cluster Membership Info
############################################
# This loop will pull all cluster membership files and load them into a list 
# according to the order of the datasets passed as arguments
ClusterMembershipList <- list()
ClusterMembershipList_NMF <- list()
datasetMembers <- list.files(path = "2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/")
datasetMembers_NMF <- list.files(path = "2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/") 
for (dataset in 1:length(args)) {
  holder <- args[dataset]
  ClusMember <- datasetMembers[grep(holder, datasetMembers)]
  
  ClusData <- read.csv(paste("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/", 
                             ClusMember, sep = ""), row.names = 1)
  
  ClusterMembershipList[[dataset]] <- ClusData
  names(ClusterMembershipList)[dataset] <- args[dataset]
  
  # NMF cluster membership
  ClusMember <- datasetMembers_NMF[grep(holder, datasetMembers_NMF)]
  ClusMember <- ClusMember[grepl("mapped", ClusMember)]
  ClusData <- read.csv(paste("2.Clustering_DiffExprs/Tables/ClusterMembership/nmf/", 
                             ClusMember, sep = ""), row.names = 1)
  
  ClusterMembershipList_NMF[[dataset]] <- ClusData
  names(ClusterMembershipList_NMF)[dataset] <- args[dataset]
}

############################################
# Get data frames ready for CoxPh
############################################
CoxPHdata <- GetCoxPHready(phenoData, ClusterMembershipList)
CoxPHdata_NMF <- GetCoxPHready(phenoData, ClusterMembershipList_NMF)

############################################
# Kaplan Meier curves
############################################
for (coxmodel in 1:length(CoxPHdata)) {
  doCoxPH_KM(CoxPHdata[[coxmodel]], as.character(paste(names(CoxPHdata)[coxmodel])))
  doCoxPH_KM(CoxPHdata_NMF[[coxmodel]], paste0(as.character(paste(names(CoxPHdata_NMF)[coxmodel])), '_NMF_'))
}

############################################
# CoxPH tables and KM curves
############################################
# Perform Cox Proportional Hazards, each dataset has different covariate data to select
# Run custom cox function on all available variables ("multi") and all available variables minus
# age ("removeAge")
type = c("multi", "removeAge")

PropHaz <- list()
for (dataset in 1:length(CoxPHdata)) {
  typeList <- list()
  for (model in 1:length(type)) {
    typeList[[type[model]]] <- customCoxPH(CoxPHdata[[dataset]], type[model])
  }
  PropHaz[[names(CoxPHdata)[dataset]]] <- typeList
}

############################################
# Write out individual survival files
############################################
for(dataset in 1:length(PropHaz))
{
  for(type in 1:length(PropHaz[[dataset]]))
  {
    for(centroid in 1:length(PropHaz[[dataset]][[type]]))
    {
      # Get the name of the file
      name <- paste(names(PropHaz)[dataset], names(PropHaz[[dataset]])[type], 
                    names(PropHaz[[dataset]][[type]])[centroid], "survival", sep = "_")
      
      # Run the function to output cox proportional hazards analysis
      coxSum(PropHaz[[dataset]][[type]][[centroid]], name)
    }
  }
}
