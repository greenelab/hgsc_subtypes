############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will identify k-means clustering concordance with the original TCGA 
# subtypes identified in the 2011 Nature Paper, the original Tothill subtypes 
# identified in the 2008 Clin Cancer Research Paper, and the original Konecny 
# subtypes identified in the 2014 JNCI paper

library(curatedOvarianData)

############################################
# Load our TCGA cluster assignments
############################################
TCGAClusterAssign <- read.csv("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/KMembership_TCGA_eset.csv", 
                              row.names = 1)

# Data from Verhaak et al. 2013 Supplemental Information
# (Where we are extracting TCGA cluster membership information from)
Verhaak <- read.csv("2.Clustering_DiffExprs/Data/Verhaak2013_Supplemental.csv")

############################################
#Subset Data
############################################
# Subset the Verhaak file to only those samples labeled by the original TCGA discovery
TCGA_Nature <- subset(Verhaak, Verhaak$X.1 == "TCGA-discovery")[ ,c(1, 3)]

# The names are separated by a "-" instead of a "."
rownames(TCGA_Nature) <- gsub("-", ".", TCGA_Nature$X)

# Get the samples in common
commonSamples <- intersect(rownames(TCGA_Nature), rownames(TCGAClusterAssign))

# Subset our cluster assignment to only accept commonsamples
KmeansSub <- TCGAClusterAssign[commonSamples, ]

# Combine the two assignments
withTCGAAssignments <- cbind(KmeansSub, TCGA_Nature[commonSamples, ])
withTCGAAssignments <- withTCGAAssignments[ ,c(1,2,3,5)]

# Look at what was not mapped to the cluster assignments
whatsleft <- TCGAClusterAssign[setdiff(rownames(TCGAClusterAssign), commonSamples), ]

############################################
# Analyze concordance
############################################
# Look at which of our clusters map to the TCGA assigned clusters
K2 <- table(withTCGAAssignments$ClusterK2, withTCGAAssignments$Patient.characteristics)
K3 <- table(withTCGAAssignments$ClusterK3, withTCGAAssignments$Patient.characteristics)
K4 <- table(withTCGAAssignments$ClusterK4, withTCGAAssignments$Patient.characteristics)

# What samples did not map?
K2NA <- table(whatsleft$ClusterK2)
K3NA <- table(whatsleft$ClusterK3)
K4NA <- table(whatsleft$ClusterK4)

# Reorder the table
reorderK2 <-cbind(K2[ ,3], K2[ ,4], K2[ ,2], K2[ ,1], K2NA)
colnames(reorderK2) <- c("Mesenchymal", "Proliferative", "Immunoreactive", "Differentiated", "NotMapped")
reorderK3 <- cbind(K3[ ,3], K3[ ,4], K3[ ,2], K3[ ,1], K3NA)
colnames(reorderK3) <- c("Mesenchymal", "Proliferative", "Immunoreactive", "Differentiated", "NotMapped")
reorderK4 <- cbind(K4[ ,3], K4[ ,4], K4[ ,2], K4[ ,1], K4NA)
colnames(reorderK4) <- c("Mesenchymal", "Proliferative", "Immunoreactive", "Differentiated", "NotMapped")

############################################
# Write to File
############################################
write.table(reorderK2, "2.Clustering_DiffExprs/Tables/Data_Concordance/TCGA_KmeansClusterK2.csv", 
            row.names = T, col.names = NA, sep = ",")
write.table(reorderK3, "2.Clustering_DiffExprs/Tables/Data_Concordance/TCGA_KmeansClusterK3.csv", 
            row.names = T, col.names = NA, sep = ",")
write.table(reorderK4, "2.Clustering_DiffExprs/Tables/Data_Concordance/TCGA_KmeansClusterK4.csv", 
            row.names = T, col.names = NA, sep = ",")

############################################
# Analyze Konecny cluster membership
############################################
# Load the file that stores our cluster assignments
KonecnyClusterAssign <- read.csv("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/KMembership_mayo.eset.csv", row.names = 1)

# This is a supplementary table in the Konecney supplemental material, it is indexed by "ov number"
Kon_Sup <- read.csv("2.Clustering_DiffExprs/Data/Konecny_supplemental.csv", row.names = 1)

# Load the covariate file to obtain a dictionary describing overlapping sample IDs
Kon_eset <- get(load("1.DataInclusion/Data/Mayo/MayoEset.Rda"))
Kon_Cov <- pData(Kon_eset)

# Subeset covariate file to only consider samples we used in our inclusion set
Kon_Cov <- Kon_Cov[intersect(rownames(KonecnyClusterAssign), rownames(Kon_Cov)), ]

# Ov number samples in common between covariate file and supplemental table
commonSamplesOV <- intersect(as.character(Kon_Cov$ovnum), rownames(Kon_Sup))

# Subset the covariate file supplemental table
Kon_Cov <- Kon_Cov[match(commonSamplesOV, Kon_Cov$ovnum), ]
Kon_Sup <- Kon_Sup[commonSamplesOV, ]

# Subset to only samples that were used in the original publication
OvNumberDict <- cbind(rownames(Kon_Cov[match(commonSamplesOV, as.character(Kon_Cov$ovnum)), ]), 
                      Kon_Cov[match(commonSamplesOV, as.character(Kon_Cov$ovnum)), ][ ,2])

# Finally, get the total intersection appropriately ordered
Kon_Comp <- cbind(KonecnyClusterAssign[match(OvNumberDict[ ,1],rownames(KonecnyClusterAssign)),], 
                  Kon_Sup[match(OvNumberDict[ ,2], rownames(Kon_Sup)), ][, 8:11])

# Investigate the samples that were not mapped by the original publication
whatsleft <- KonecnyClusterAssign[setdiff(rownames(KonecnyClusterAssign), rownames(Kon_Comp)), ]
K2NA <- table(whatsleft$ClusterK2)
K3NA <- table(whatsleft$ClusterK3)
K4NA <- table(whatsleft$ClusterK4)

# Observe concordance
K2_kon <- table(Kon_Comp$ClusterK2, Kon_Comp$MAYO.C4)
K3_kon <- table(Kon_Comp$ClusterK3, Kon_Comp$MAYO.C4)
K4_kon <- table(Kon_Comp$ClusterK4, Kon_Comp$MAYO.C4)

# Prepare files for output
k2Ready <- cbind(K2_kon, K2NA)
colnames(k2Ready)[ncol(k2Ready)] <- "NotMapped"

k3Ready <- cbind(K3_kon, K3NA)
colnames(k3Ready)[ncol(k3Ready)] <- "NotMapped"

k4Ready <- cbind(K4_kon, K4NA)
colnames(k4Ready)[ncol(k4Ready)] <- "NotMapped"

############################################
# Write to file
############################################
write.table(k2Ready, "2.Clustering_DiffExprs/Tables/Data_Concordance/Konecny_KmeansClusterK2.csv", 
            row.names = T, col.names = NA, sep = ",")
write.table(k3Ready, "2.Clustering_DiffExprs/Tables/Data_Concordance/Konecny_KmeansClusterK3.csv", 
            row.names = T, col.names = NA, sep = ",")
write.table(k4Ready, "2.Clustering_DiffExprs/Tables/Data_Concordance/Konecny_KmeansClusterK4.csv", 
            row.names = T, col.names = NA, sep = ",")

############################################
# Analyze Tothill cluster membership
############################################
TothillClusterAssign <- read.csv("2.Clustering_DiffExprs/Tables/ClusterMembership/kmeans/KMembership_GSE9891_eset.csv", 
                                 row.names = 1)

Tot_Sup <- read.csv("2.Clustering_DiffExprs/Data/Tothill_supplemental.csv")
Tot_Sup$ID <- paste("X", as.character(paste(Tot_Sup$ID)), sep = "")

# Load Tothill pData
library(curatedOvarianData)
data(GSE9891_eset)
Tothill <- pData(GSE9891_eset)

# Match sample names to appropriate alternate sample names stored in the pData
index <- cbind(rownames(Tothill)[match(Tot_Sup$ID, Tothill$alt_sample_name)], Tot_Sup$ID, Tot_Sup$kmeans)

# Combine our kmeans cluster assignments versus the assignments presented by Tothill et al. 2008
Comparison <- cbind(index, TothillClusterAssign[match(index[,1], rownames(TothillClusterAssign)), ])
colnames(Comparison) <- c("Sample", "Alt", "Tothill_C", "ClusterK2", "ClusterK3", "ClusterK4")

K2_tot <- table(Comparison$Tothill_C, Comparison$ClusterK2)
rownames(K2_tot) <- c(paste("C", 1:6, sep = ""), "NotMapped")

K3_tot <- table(Comparison$Tothill_C, Comparison$ClusterK3)
rownames(K3_tot) <- c(paste("C", 1:6, sep = ""), "NotMapped")

K4_tot <- table(Comparison$Tothill_C, Comparison$ClusterK4)
rownames(K4_tot) <- c(paste("C", 1:6, sep = ""), "NotMapped")

############################################
#Write to File
############################################
write.table(t(K2_tot), "2.Clustering_DiffExprs/Tables/Data_Concordance/Tothill_KmeansClusterK2.csv", 
            row.names = T, col.names = NA, sep = ",")
write.table(t(K3_tot), "2.Clustering_DiffExprs/Tables/Data_Concordance/Tothill_KmeansClusterK3.csv", 
            row.names = T, col.names = NA, sep = ",")
write.table(t(K4_tot), "2.Clustering_DiffExprs/Tables/Data_Concordance/Tothill_KmeansClusterK4.csv", 
            row.names = T, col.names = NA, sep = ",")

############################################
# Combine all tables together to create table 4
############################################
k2 <- cbind(reorderK2, t(K2_tot), k2Ready)
k3 <- cbind(reorderK3, t(K3_tot), k3Ready)
k4 <- cbind(reorderK4, t(K4_tot), k4Ready)

write.table(k2, "2.Clustering_DiffExprs/Tables/Data_Concordance/TotalConcordance_KmeansClusterK2.csv", 
            row.names = T, col.names = NA, sep = ",")
write.table(k3, "2.Clustering_DiffExprs/Tables/Data_Concordance/TotalConcordance_KmeansClusterK3.csv", 
            row.names = T, col.names = NA, sep = ",")
write.table(k4, "2.Clustering_DiffExprs/Tables/Data_Concordance/TotalConcordance_KmeansClusterK4.csv", 
            row.names = T, col.names = NA, sep = ",")
