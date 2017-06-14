# Timothy Chang 2017
# MCPcounter Analysis
# 6.Immune_Infiltrate/Scripts/C.MCPcounter.R
#
# Estimate the abundance of immune and stromal infiltration in Tothill
# using the MCPcounter algorithm
#
# Usage: Run the script
#
#     Rscript 6.Immune_Infiltrate/Scripts/C.MCPcounter.R
#
# Output:
# Tables with MCPcounter results for each sample
# Box plots by subtype for results

suppressMessages(library(checkpoint))
suppressMessages(checkpoint('2016-03-01', checkpointLocation = "."))

############################################
# Load Libraries
############################################
library(MCPcounter)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(GEOquery)

############################################
# Load Data
############################################
# Download expression data
Tothill_GEO <- getGEO(GEO = "GSE9891", getGPL = FALSE)
Tothill_eset <- exprs(Tothill_GEO[[1]])

############################################
# Get Cluster Membership Info
############################################
# Load cluster membership info for Tothill
kmeansPath <- file.path("2.Clustering_DiffExprs", "Tables", 
                        "ClusterMembership", "kmeans", 
                        "KMembership_GSE9891_eset.csv")

ClusterMembership <- read.csv(kmeansPath, stringsAsFactors = FALSE)
names(ClusterMembership)[1] <- "Sample"

############################################
# Run MCPcounter Analysis
############################################
# Run MCPcounter
mcp_result <- MCPcounter.estimate(expression = Tothill_eset, 
                                  featuresType = "affy133P2_probesets")
infiltration <- test_for_infiltration(mcp_result, platform = "133P2")

# Format result data
mcp_result <- t(mcp_result)
mcp_result <- tibble::rownames_to_column(as.data.frame(mcp_result),
                                         "Sample")

# Write results to disk
fPath <- file.path("6.Immune_Infiltrate", "Tables", "MCPcounter", 
                   "GSE9891_mcp.csv")
write.csv(mcp_result, fPath, quote = FALSE)

# Add cluster membership
mcp_result_cluster <- inner_join(mcp_result, ClusterMembership,
                                 by = "Sample")

############################################
# Output plots
############################################
k_list <- c("K2", "K3", "K4")
boxplot_theme <- theme(panel.background = element_blank(),
                       panel.border = element_rect(color = "black", fill = NA),
                       axis.text.y = element_text(size = 12, angle = 90,
                                                  hjust = 0.4, 
                                                  colour = "black"),
                       axis.text.x = element_text(size = 12, colour = "black"),
                       axis.title.y = element_text(size = 12, colour = "black",
                                                   face = "bold"))
cell_types <- colnames(mcp_result)[-1]

# Analysis with k-means clustering
for (k in k_list) {
  MCPPlots <- list()
  data_iter <- 1
  num_clus <- gsub("K", "", k)

  for (type in cell_types) {
    plot_cols <- c("Sample", type, paste0("Cluster", k))
    plot_df <- mcp_result_cluster[plot_cols]
    colnames(plot_df)[2] <- "Abundance"
    colnames(plot_df)[3] <- "Cluster"

    # Plot by subtype
    g <- ggplot (data = plot_df,
                 mapping = aes (x = Cluster, y = Abundance, group = Cluster)) +
      geom_boxplot() +
      scale_x_discrete(limits = paste(1:num_clus)) +
      labs(title = "", x = "", y = type) +
      boxplot_theme
    MCPPlots[[data_iter]] <- g

    data_iter <- data_iter + 1
  }

  # Build string to combine plots
  plot_eval <- "full_grobs <- grid.arrange("
    for (i in 1:length(MCPPlots)) {
      plot_eval <- paste0(plot_eval, "MCPPlots[[", i, "]]", ", ")
    }
    plot_eval <- paste0(plot_eval, "ncol = 3", 
                        ", nrow = ceiling(length(MCPPlots) / 3))")

  # Write plot to PNG
  fName <- paste0("GSE9891", "_", k, "_MCP.png")
  fPath <- file.path("6.Immune_Infiltrate", "Figures", "MCPcounter", 
                     fName)
  png(filename = fPath,
        width = 400 * 3,
        height = 400 * ceiling(length(MCPPlots) / 3))
    eval(parse(text = plot_eval))
    dev.off()
}
