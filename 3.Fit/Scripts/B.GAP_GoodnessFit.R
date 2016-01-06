############################################
# Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will output the GAP statistic for a predefined number of bootstraps. 
# The script is computationally expensive so it would be useful to run on a cluster.

args <- commandArgs(trailingOnly=TRUE)
args <- c(8, 200, 20, 20, 123, "TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset")
############################################
# Load Libraries
############################################
library(curatedOvarianData)
library(cluster)

# The script loads the ovarian cancer datasets
source("1.DataInclusion/Scripts/Functions/LoadOVCA_Data.R")

############################################
# Constants
############################################
KMAX <- as.numeric(paste(args[1]))
Boot <- as.numeric(paste(args[2])) #Number of bootstrapping iterations
K.iter <- as.numeric(paste(args[3])) #number of kmeans iterations
nstart <- as.numeric(paste(args[4])) #nstarts for kmeans
SEED <- as.numeric(paste(args[5]))

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
GlobalMAD <- read.table("1.DataInclusion/Data/Genes/GlobalMAD_genelist.csv", sep = ",", header =T)

############################################
# Run GAP
############################################
set.seed(SEED)
for (dataset in 1:length(ExpData)) {
  # The genes in the global MAD variable do not overlap with the datasets completely
  SubsetMAD <- intersect(rownames(ExpData[[dataset]]), GlobalMAD$x)
  
  # Run the clusGap function
  GAP <- clusGap(t(ExpData[[dataset]][SubsetMAD, ]), FUN = kmeans, nstart = nstart, 
                 iter.max = K.iter, K.max = KMAX, B = Boot)

  # Extract the tibs criterion from the results
  tibs_criterion <- maxSE(GAP$Tab[ ,"gap"], GAP$Tab[ ,"SE.sim"])
  cat("200 bootstraps completed:", paste(names(ExpData)[dataset], "\n"))
  
  # Get the plot ready
  tiff(file.path("3.Fit", "Figures", "Gap", paste(names(ExpData)[dataset], "GAP.tiff",sep = "_")), 
       width = 1200, height=1200)
  
  # Get proper margins
  par(mar = rep(2.4, 4))
  par(cex = 3)
  
  plot(GAP, main = "", xlab = "", ylab = "", cex = 2, cex.axis = 1.2, pch = 20)
  abline(v = tibs_criterion, col = 4, lty = 3, lwd = 2)
  dev.off()
  
  # Since the script takes a long time to run, save the outputs to tables
  write.table(GAP$Tab, file = paste("3.Fit/Tables/", paste(names(ExpData)[dataset], 
                                                           "_GAP_data.csv", sep = ""), sep = "/"), 
              row.names = T, col.names = NA, sep = ",")
}
