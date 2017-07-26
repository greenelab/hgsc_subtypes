###########################################
# Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# Reproduce Supplementary Figure S6.2 using our genes in TCGA OvCa paper 
# except without removing LMP samples to contrast our biological filtering
# findings

set.seed(123)

################################
# Load Libraries
################################
library(curatedOvarianData)
library(cluster)
library(NMF)

# Load important kmeans and SAM functions
source("2.Clustering_DiffExprs/Scripts/Functions/kmeans_SAM_functions.R")

# Loads curatedOvarianData
source("1.DataInclusion/Scripts/Functions/LoadOVCA_Data.R")

################################
# Load and Subset Tothill Data
################################
# Load eset
data("GSE9891_eset")

# Load phenotype data
Tothill_pdata <- pData(GSE9891_eset)

# Load samples used in clustering
fileName <- "GSE9891_eset_samplesRemoved.csv"
Tothill_goodsamples <- read.csv(file = paste("1.DataInclusion/Data/GoodSamples/", fileName, sep = ""))

# Get the list of samples removed
Tothill_badsamples <- setdiff(rownames(Tothill_pdata), Tothill_goodsamples[ ,2])

# Get the list of LMP samples removed
Tothill_lmpsamples <- Tothill_pdata[Tothill_badsamples, ]
Tothill_lmpsamples <- rownames(Tothill_lmpsamples[Tothill_lmpsamples$sample_type == 'borderline', ])

# Add these samples to the list of good samples 
Tothill_goodsampleswlmp <- c(as.character(paste(Tothill_goodsamples[ ,2])), Tothill_lmpsamples)

# Load and subset expression data
madgenes <- read.csv("1.DataInclusion/Data/Genes/GlobalMAD_genelist.csv", header = T, stringsAsFactors = F)
Tothill_exp <- exprs(GSE9891_eset)[madgenes$x, Tothill_goodsampleswlmp]

################################
# Perform NMF consensus analysis
################################
# Run custom NMF function with 10 runs to output cophenetic coefficient
coph_coeff <- runNMF(Tothill_exp, coph = T)

# Base name to save figures
sup_fname <- "NMF_ConsensusMatrices_Tothill_LMP"

# Write the consensus mapping as a figure
png(paste("2.Clustering_DiffExprs/Figures/nmf/CopheneticMaps/", 
          paste(sup_fname, '_k2-6', sep = ""), ".png", sep = ""), 
    width = 950, height = 550)

# Plot consensus maps
consensusmap(coph_coeff$consensus[1:5], labCol = NA, labRow = NA, main = "", fontsize = 12)

# Close the connection
dev.off()

# Write the Cophenetic Coefficient as a figure
png(paste("2.Clustering_DiffExprs/Figures/nmf/CopheneticMaps/", 
          paste(sup_fname, 'coph_coeff', sep = ""), ".png", sep = ""), 
    width = 270, height = 230)

# Plot cophenetic coefficient
par(mar = c(4.5, 4.5, 1.5, 1))
plot(coph_coeff$measures$cophenetic, xaxt = "n", cex.axis = 1.5, cex.lab = 1.5,
     xlab = 'k', ylab = 'Cophenetic Correlation')
axis(1, at=1:7, labels=2:8, cex.axis = 1.5)
lines(coph_coeff$measures$cophenetic, lwd = 2)
points(coph_coeff$measures$cophenetic, col = 'black', pch = 19, cex = 1.2)

# Close the connection
dev.off()

