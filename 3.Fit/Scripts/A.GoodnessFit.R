############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will output the AIC, BIC, and Silhouette Width plots

args <- commandArgs(trailingOnly=TRUE)
# args <- c(2, 8, 20, 20, 123, "TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset")
############################################
# Load Libraries
############################################
library(curatedOvarianData)
library(sfsmisc)
library(cluster)
library(ggplot2)

############################################
# Constants
############################################
kmin = as.numeric(paste(args[1]))
kmax = as.numeric(paste(args[2]))  
itermax = as.numeric(paste(args[3])) 
starts = as.numeric(paste(args[4])) 
SEED = as.numeric(paste(args[5])) 

# Load custom functions
source("3.Fit/Scripts/Functions/Fit_Functions.R")

# The script loads the ovarian cancer datasets
source("1.DataInclusion/Scripts/Functions/LoadOVCA_Data.R")

############################################
#Data in curatedOvarianData
############################################
# Get data from the curatedOvarianData package  
detailedData <- data(package="curatedOvarianData")[3]
# Get the eset IDs
detailedData.names <- detailedData$results[ ,3]
# Determine the datasets in curatedOvarianData
inCuratedOvarianData <- args %in% detailedData.names

# Here are the args that where curated
argsCurated <- args[inCuratedOvarianData]

# Add Mayo to argsCurated
if("Mayo" %in% args) {
  argsCurated = c(argsCurated[1], "Mayo", argsCurated[2:length(argsCurated)])
}

############################################
# Load Data
############################################
# Use the LoadOVCA_Data function to read in the datasets subset by commongenes
ExpData <- LoadOVCA_Data(datasets = argsCurated, genelist_subset = "commongenes")

# Get gene lists
commongenes <- read.table("1.DataInclusion/Data/Genes/CommonGenes_genelist.csv", sep = ",", header = T)
GlobalMAD <- read.table("1.DataInclusion/Data/Genes/GlobalMAD_genelist.csv", sep = ",", header =T)

############################################
#Output AIC, BIC, and Silhouette Plots
############################################
set.seed(SEED)

for (dataset in 1:length(ExpData)) {
  Name <- names(ExpData)[dataset]
  
  # Run the two functions to output AIC, BIC and Silhouette plots for each dataset
  #this <- KMeansElbow(ExpData[[dataset]], Name)
  that <- RunSilhouette(ExpData[[dataset]], Name, MADgenes = GlobalMAD$x)
}
