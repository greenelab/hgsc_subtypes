# Install rforge and custom github packages
#
# Usage:
#
#    After activating conda environment (`source activate hgsc_subtypes`):
#
#          R --no-save < install_custom.R

library(devtools)

# Install ESTIMATE
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos = rforge, dependencies = TRUE)

# Install MCPcounter
install_github("ebecht/MCPcounter", ref = "master", subdir = "Source",
               force = TRUE)
