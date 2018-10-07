#!/bin/bash
#
# The following pipeline defines how the scripts should be run in sequential order.
# By running `bash 7.Nanostring/run_analysis.sh`, the entire nanostring module will be reproduced.
#
# NOTE: all scripts should be run from the top `hgsc_subtypes` directory

set -o errexit
exec &>7.Nanostring/nanostring_analysis.out

# Step 0 - Install package inside conda environment
Rscript 7.Nanostring/scripts/0.install_custom.R

# Step 1 - Obtain the correlation output files
Rscript 7.Nanostring/scripts/A.get_correlation_output.R

# Step 2 - Explore these correlations, output figures, and several descriptive results
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=7.Nanostring/scripts \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000 \
        --execute 7.Nanostring/scripts/B.explore_correlations.ipynb

# Step 3 - Visualize differences with a series of Venn diagrams
Rscript 7.Nanostring/scripts/C.threshold_venns.R

# Step 4 - Obtain overlapping genes across datasets
Rscript 7.Nanostring/scripts/D.get_overlap_genes.R

# Step 5 - Perform pathway analysis on genesets
Rscript 7.Nanostring/scripts/E.gmt_genesets_pathway_analysis.R

# Step 6 - Summarize the results of the pathway analysis
Rscript 7.Nanostring/scripts/F.summarize_pathways.R
