############################################
# High-grade serous ovarian cancer subtypes are similar across populations
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will input a series of datasets and output bar charts

args <- commandArgs(trailingOnly=TRUE)
# args <- c("3", "4", "TCGA_eset", "Mayo", "GSE32062.GPL6480_eset", "GSE9891_eset", "GSE26712_eset")

############################################
# Load Libraries
############################################
library(ggplot2)
library(RColorBrewer)
library(reshape2)

############################################
# Constants
############################################
k = as.numeric(paste(args[1])) 
k2 = as.numeric(paste(args[2]))  # Number of clusters 
files = as.character(paste(args[3:length(args)]))

for (dataset in 1:length(files)) {
  # Get membership file name
  membfile = as.character(paste(files[dataset]))
  
  # Load in cluster assignment dataframe
  ClusterAssign <- read.csv(file = file.path("2.Clustering_DiffExprs", "Tables", "ClusterMembership", "kmeans", 
                                             paste("KMembership_", membfile, ".csv", sep = "")), row.names = 1)

  # Output the k-means colored bars
  fName <- paste(paste("2.Clustering_DiffExprs/Figures/barcharts/", membfile, sep = ""), 
                 "KMeansBarChart.tiff", sep ="_")
  
  ggplot(ClusterAssign, aes(factor(ClusterAssign$ClusterK4), fill = factor(ClusterAssign$ClusterK3))) + 
    geom_bar() + xlab("") + ylab("") + scale_fill_manual(values = c("blue", "red", "green"), 
                                                                          name = "") +
    theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = "none", panel.background = element_blank(), panel.border = element_blank(), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              plot.background = element_blank(), axis.text.y = element_text(size = 18))
  
  ggsave(file.path(fName), width = 8, height = 8)
}