############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
#  
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# Output gene lists (common genes and MAD genes) and
# Overlapping Genes Venn Diagram: Supplemental Figure 1

suppressMessages(library(checkpoint))
suppressMessages(checkpoint('2016-03-01', checkpointLocation = "."))

args <- commandArgs(trailingOnly = T)
# args <- c("TCGA_eset", "mayo.eset", "GSE32062.GPL6480_eset", "GSE9891_eset", "GSE26712_eset", "aaces.eset")
# args <- c("TCGA_eset", "aaces.eset")
################################
# Load Libraries
################################
library(limma)
library(grid)
library(curatedOvarianData)
library(VennDiagram)

# The script loads the ovarian cancer datasets
source("1.DataInclusion/Scripts/Functions/LoadOVCA_Data.R")

#Library has some important custom functions
source("2.Clustering_DiffExprs/Scripts/Functions/kmeans_SAM_functions.R")

# Load venn diagram function
getvenn <- function(venngenes, data_set_column) {
  venngenes <- as.data.frame(venngenes)
  return_venn_count <- c()
  getcolumn <- venngenes[ ,data_set_column]
  for (g in 1:length(getcolumn)) {
    if (getcolumn[g] == 1) {
      return_venn_count <- c(return_venn_count, g)
    }
  }
  return(return_venn_count)
}

################################
# PART I: Load Data
################################
# Use the LoadOVCA_Data function to read in the datasets subset by commongenes
ExpData <- LoadOVCA_Data(datasets = args, genelist_subset = "None")

 # Get all the common genes and universe genes
commongeneSet <- c()
universe <- c()
for (dataset in 1:length(ExpData)) {
  datasetExprs <- ExpData[[dataset]]
  
  # Extract the gene names for each dataset
  datasetGenes <- rownames(datasetExprs)
  #datasetGenes.aliases <- datasetGenes[grepl("///", datasetGenes)]
  datasetgenes.alone <- datasetGenes[!grepl("///", datasetGenes)]
  
  # If treating multiple mappings as unique genes, which for now we aren't:
  # if (length(datasetGenes.aliases) > 0){
  #    datasetgenes.alone <- c(datasetgenes.alone,
  #                            lapply( datasetgenes.alone,
  #                            function(x) unlist(strsplit(x, "///"))[[1]]))
  #
  #     }
  #   }
  # } 
  
  # We are also interested in all the genes that are measured
  universe <- unique(c(universe, datasetgenes.alone))
  if (dataset == 1) {
    commongeneSet <- datasetgenes.alone
  } else {
    commongeneSet <- intersect(commongeneSet, datasetgenes.alone)
  }
}
# Get MAD genes
MADgeneSet <- c()
for (dataset in 1:length(ExpData)) {
  # Subset all datasets to common genes
  datasetExprs <- ExpData[[dataset]][commongeneSet, ]
  
  # Use the MADgenes function to extract the top 1500 most variably expressed genes
  madSet <- MADgenes(datasetExprs, numGenes = 1500)
  
  MADgeneSet <- unique(c(MADgeneSet, madSet))
}

# Write to file
write.csv(commongeneSet, file = "1.DataInclusion/Data/Genes/CommonGenes_genelist.csv", 
          row.names = F)
write.csv(MADgeneSet, file = "1.DataInclusion/Data/Genes/GlobalMAD_genelist.csv", 
          row.names = F)

################################
# PART II: Venn Diagram (Supplementary Figure S1a)
################################
# Initialize the venn matrix, which will serve as a template for the observing
# overlaps between the datasets (Supplementary Figure S1)
# Note: This step is only performed when AACES is excluded
if (!("aaces.eset" %in% args)){
  venn <- matrix(NA, nrow = length(universe), ncol = length(args))
  rownames(venn) <- universe
  colnames(venn) <- c("TCGA", "Mayo", "Yoshihara", "Tothill", "Bonome")
  
  # Create the Venn Matrix
  for (gene in 1:length(universe)) {
    vector <- c()
    for (gene_sp in 1:length(args)) {
      if (universe[gene] %in% rownames(ExpData[[gene_sp]])) {
        tmp <- 1
      } else {
        tmp <- 0
      }
      vector <- c(vector, tmp)
    }
    venn[gene, ] <- vector
  }
  
  # Generate and output Venn Diagram
  TCGA_venn <- getvenn(venn, 1)
  Mayo_venn <- getvenn(venn, 2)
  Yoshihara_venn <- getvenn(venn, 3)
  Tothill_venn <- getvenn(venn, 4)
  Bonome_venn <- getvenn(venn, 5)
  venn.file.path <- file.path("1.DataInclusion", "Figures",
                              "OverlappingGenesVenn")
  venn.plot <- venn.diagram(x = list('TCGA' = TCGA_venn,
                                     'Mayo' = Mayo_venn,
                                     'Yoshihara' = Yoshihara_venn,
                                     'Tothill' = Tothill_venn,
                                     'Bonome' = Bonome_venn),
                            filename = venn.file.path,
                            height = 3000, width = 3000,
                            fill = c("red", "green", "yellow", "blue", "purple"),
                            cat.cex = rep(1.6, 5),
                            cat.pos = c(0, -20, 20, -20, 20),
                            cat.dist = c(.19, .21, -.18, -.2, .21),
                            margin = .1, cex = 0.8,
                            main.cex = 2)
  
}
