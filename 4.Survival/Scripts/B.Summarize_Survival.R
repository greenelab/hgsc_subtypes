############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will summarize all of the multi-variate and univariate survival analyses

############################################
args <- commandArgs(trailingOnly=TRUE)
# args <- c("TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset")

library(curatedOvarianData)

# How many centroids assigned to k means algorithm
clus <- c(2, 3, 4)

############################################
# Load Results
############################################
files <- list.files("4.Survival/Tables/")

# Upload dataset specific survival tables
survivalData <- list()
for (dataset in 1:length(args)) {
  
  # All the files associated with the given dataset
  survs <- files[grepl(args[dataset], files)]

  datasetSurv <- list()
  for (type in 1:length(survs)) {
    
    # Upload all types of cox models and store in dataset indexed list
    datasetSurv[[survs[type]]] <- read.csv(paste("4.Survival/Tables/",survs[type], sep = ""), 
                                           header = T, row.names = 1, sep = ",")
  }
  survivalData[[args[dataset]]] <- datasetSurv
}

############################################
# Compile Months to Death Data
############################################
detailedData <- data(package = "curatedOvarianData")[3]
detailedData.names <- detailedData$results[ ,3]
inCuratedOvarianData <- args %in% detailedData.names

argsCurated <- args[inCuratedOvarianData]

# First, load all pheno data, and convert while loading
phenoData <- list()
monthstodeath <- list()
for (dataset in 1:length(args)) {
  
  # Load the phenotype data
  dta <- c()
  if (args[dataset] %in% detailedData.names) {
    # Get the list of samples to remove
    fName <- paste(args[dataset], "_samplesRemoved.csv",sep="")
    good <- read.csv(file = file.path("1.DataInclusion", "Data", "GoodSamples", fName))
    accept <- c(paste(good$x))

    # Load the phenotype data
    data(list = args[dataset], package = "curatedOvarianData")
    tmp <- get(args[dataset])
    dta <- pData(tmp)
    phenoData[[dataset]] <- dta[accept, ]
    
    # Get months to death median value and 95% range
    dtdmedian <- round(median(as.numeric(paste(dta$days_to_death)), na.rm = T) / 30.4375, digits = 0)
    dtdmin <- round(quantile(as.numeric(paste(dta$days_to_death)), probs = 0.05, na.rm = T) / 30.4375, digits = 0)
    dtdmax <- round(quantile(as.numeric(paste(dta$days_to_death)), probs = 0.95, na.rm = T) / 30.4375, digits = 0)
    
    # Combine into one vector and store in list
    dtdvector <- c(dtdmedian, paste(dtdmin, dtdmax, sep = " - "))
    names(dtdvector) <- c("Median", "95% CI" )
    monthstodeath[[dataset]] <- dtdvector
    
    # Name the list
    names(phenoData)[dataset] <- names(monthstodeath)[dataset] <- args[dataset]
   
  } else {
    if (args[dataset] == "Mayo") {
      
      # Read in file and subset according to good samples
      dta <- pData(get(load("1.DataInclusion/Data/Mayo/MayoEset.Rda")))
      Mayogood <- read.csv(file = file.path("1.DataInclusion", "Data", "GoodSamples", "Mayo_samplesRemoved.csv"))
      Mayoaccept <- c(paste(Mayogood$x))
      rownames(dta) <- dta[ ,1]
      phenoData[[dataset]] <- dta[Mayoaccept,]
      
      # Get months to death median value and 95% range
      dtdmedian <- round(median(as.numeric(paste(dta$days_to_death)), na.rm = T)  / 30.4375, digits = 0)
      dtdmin <- round(quantile(as.numeric(paste(dta$days_to_death)), probs = 0.05, na.rm = T) / 30.4375, digits = 0)
      dtdmax <- round(quantile(as.numeric(paste(dta$days_to_death)), probs = 0.95, na.rm = T) / 30.4375, digits = 0)
      
      # Combine into one vector and store in list
      dtdvector <- c(dtdmedian, paste(dtdmin, dtdmax, sep = " - "))
      names(dtdvector) <- c("Median", "95% CI" )
      monthstodeath[[dataset]] <- dtdvector
      
      # Name the list
      names(phenoData)[dataset] <- names(monthstodeath)[dataset] <- args[dataset]    
    }
  }
}

############################################
# Get Cluster Membership Info
############################################
# This loop will pull all cluster membership files and load them into a list 
# according to the order of the datasets passed as arguments
ClusterMembershipList <- list()
datasetMembers <- list.files(path = "2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/")
for (dataset in 1:length(args)) {
  holder <- args[dataset]
  ClusMember <- datasetMembers[grep(holder, datasetMembers)]
  
  ClusData <- read.csv(paste("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/", 
                             ClusMember, sep = ""), row.names = 1)
  
  ClusterMembershipList[[dataset]] <- ClusData
  names(ClusterMembershipList)[dataset] <- args[dataset]
}

############################################
# Compile and Write Results
############################################
dataSurvivalCuration <- list()
for (dataset in 1:length(survivalData)) {
  
  # Load all the files associated with the particular dataset
  survtmp <- survivalData[[dataset]]
  ptmp <- phenoData[[dataset]]
  clusterM <- ClusterMembershipList[[dataset]]
  
  # Initialize a matrix that will hold the survival information
  survMatrix <- matrix(NA, nrow = 9, ncol = 7)
  
  for (centroid in 1:length(clus)) {
    subset <- survtmp[grepl(paste("K", clus[centroid], sep = ""), names(survtmp))]
    
    # Make a new matrix that will hold the results of the full and partially adjusted models
    modMatrix <- matrix(NA, nrow = clus[centroid], ncol = 4)
    
    # Loop over each of the models
    for (type in 1:length(subset)) {
      clusterInfo <- c()
      
      # We are only interested in showing the results for each of the cluster variables in the model
      for(clustering in 1:length(grep("Cluster", rownames(subset[[type]])))) {
        
        # Round the hazard ratios to 1 decimal places
        hazard <- format(round(as.numeric(paste(subset[[type]][clustering, 3])), digits = 1), nsmall = 1)
        
        # Write a string that holds the confidence intervals in parenthesis
        CI <- paste("(", format(round(as.numeric(paste(subset[[type]][clustering, 1])), digits = 1), nsmall = 1), 
                    " to ", format(round(as.numeric(paste(subset[[type]][clustering, 2])), digits = 1), nsmall = 1),
                    ")", sep = "")
        
        # Round p value to 3 digits if it is greater than 0.0005
        if (as.numeric(paste(subset[[type]][clustering, 4])) > 0.0005) {
          p <- format(round(as.numeric(paste(subset[[type]][clustering, 4])), digits = 3), nsmall = 3)
        } else {
          p <- as.numeric(paste(subset[[type]][clustering, 4]))
        }
        
        # Get the info for the specific cluster
        info <- c(paste(hazard, CI), p)
        
        if(clustering == 1) {
          clusterInfo <- info
        } else {
          clusterInfo <- rbind(clusterInfo, info)
        }
      }
      
      # Begin building the "modMatrix" matrix (modMatrix will have different dimensions 
      # depending on clustering)
      if (type == 1 & centroid == 1) {
        # the first row will be cluster 1 (which is the reference category in each model)
        modMatrix[1, 1:2] <- c("1.0 (ref.)", "")
        modMatrix[2, 3:4] <- clusterInfo
        
      } else if (type == 2 & centroid == 1) {
        modMatrix[1, 3:4] <- c("1.0 (ref.)", "")
        modMatrix[2, 1:2] <- clusterInfo
        survMatrix[1:2, 4:7] <- modMatrix
        
      } else if (type == 1 & centroid == 2) {
        # the first row will be cluster 1 (which is the reference category in each model)
        modMatrix[1, 1:2] <- c("1.0 (ref.)", "")
        modMatrix[2:3, 3:4] <- clusterInfo
        
      } else if (type == 2 & centroid == 2) {
        modMatrix[1, 3:4] <- c("1.0 (ref.)", "")
        modMatrix[2:3, 1:2] <- clusterInfo
        survMatrix[3:5, 4:7] <- modMatrix
        
      } else if (type == 1 & centroid == 3) {
        modMatrix[1, 1:2] <- c("1.0 (ref.)", "")
        modMatrix[2:4, 3:4] <- clusterInfo  
        
      } else if (type == 2 & centroid == 3) {
        modMatrix[1, 3:4] <- c("1.0 (ref.)", "")
        modMatrix[2:4, 1:2] <- clusterInfo
        survMatrix[6:9, 4:7] <- modMatrix
      }
    }
  }
  
  # Obtain survival information
  for (memb in 1:ncol(clusterM)) {
    clusterAssign <- clusterM[ ,memb]
    for (subtype in 1:length(unique(clusterAssign))) {
      samples <- rownames(clusterM[clusterM[ ,memb] == subtype, ])
      
      # Count number of samples and perfom the appropriate subset
      numT <- length(samples)
      dta <- ptmp[samples, ]
        
      # Get months to death
      mtdmedian <- round(median(as.numeric(paste(dta$days_to_death)), na.rm = T) / 30.4375, 
                         digits = 0)
      mtdmin <- round(quantile(as.numeric(paste(dta$days_to_death)), 
                               probs = 0.05, na.rm = T) / 30.4375, digits = 0)
      mtdmax <- round(quantile(as.numeric(paste(dta$days_to_death)), probs = 0.95, 
                               na.rm = T) / 30.4375, digits = 0)
      
      num <- length(dta$days_to_death[complete.cases(dta$days_to_death)])
         
      # Get months to death vector
      mtdvector <- paste(mtdmedian, paste("(", mtdmin, " - ", mtdmax, ")",sep = ""))
      
      # Store survival information to extract later
      if (memb == 1) {
        survMatrix[subtype, 2] <- paste(num, "/", numT, sep = "")
        survMatrix[subtype, 3] <- mtdvector
      } else if (memb == 2) {
        survMatrix[(subtype + 2), 2] <- paste(num, "/", numT, sep = "")
        survMatrix[(subtype + 2), 3] <- mtdvector
      } else {
        survMatrix[(subtype + 5), 2] <- paste(num, "/", numT, sep = "")
        survMatrix[(subtype + 5), 3] <- mtdvector
      }
    }
  }
  
  # Store specific info in the survMatrix
  survMatrix[1, 1] <- paste(monthstodeath[[dataset]][1], " (",
                            monthstodeath[[dataset]][2], ")", sep = "")
  
  
  colnames(survMatrix) <- c("Total Median Survival", "Events", "Median Survival (95% CI)", 
                            "Partial: HR (95% CI)", "Partial: Wald's P", "Full: HR (95% CI)", "Full: Wald's P")
  rownames(survMatrix) <- c("k=2: Cluster 1", "k=2: Cluster 2", "k=3: Cluster 1",
                            "k=3: Cluster 2", "k=3: Cluster 3", "k=4: Cluster 1", 
                            "k=4: Cluster 2", "k=4: Cluster 3", "k=4: Cluster 4")
  
  write.table(survMatrix, file = paste("4.Survival/Summary/", names(survivalData)[dataset], "_", 
                                       "survivalSummary.csv", sep = ""), sep = ",", row.names = T, 
              col.names = NA)
  dataSurvivalCuration[[names(survivalData)[dataset]]] <- survMatrix
}
