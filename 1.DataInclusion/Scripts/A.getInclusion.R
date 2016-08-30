############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
#
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will perform our inclusion criteria on all datasets included 
# in curatedOvarianData and a dataset from the Mayo Clinic

suppressMessages(library(checkpoint))
suppressMessages(checkpoint('2016-03-01', checkpointLocation = "."))

####################################
# Load Libraries
####################################
library(curatedOvarianData)
library(doppelgangR)
library(reshape2)
library(outliers)

# This R script holds custom inclusion functions
source("1.DataInclusion/Scripts/Functions/Inclusion_functions.R")

####################################
# Load Constants
####################################
vars <- c("sample_type", "histological_type", "grade", "primarysite", "arrayedsite", "summarystage", 
          "tumorstage", "substage", "pltx", "tax", "neo", "recurrence_status", "vital_status", "os_binary",
          "relapse_binary", "site_of_tumor_first_recurrence", "primary_therapy_outcome_success", "debulking")

minimumSamples <- 100

# Esets that should never be included: Dressman and Bentink
excludeEsets <- c("PMID17290060_eset", "E.MTAB.386_eset")

####################################
# Analysis
####################################
# All the datasets within the curatedOvarianData package
esets <- getAllDataSets("curatedOvarianData")

# load the Konecny data (If an eset is not given here, please skip this step)
load("1.DataInclusion/Data/Mayo/MayoEset.Rda")

# Use the inclusion/exclusion decision tree to filter samples in all curatedOvarainData datasets
inclusionTable <- exclusionTable(esets)

# Use the inclusion/exclusion decision tree on the Mayo data
inclusionTable.mayo <- simpleExclusion(mayo.eset)

# Combine the inclusion results from curatedOvarianData and Mayo
inclusionTable[[1]] <- cbind(inclusionTable[[1]], inclusionTable.mayo[[1]])

# The first list element in the inclusionTable is the data.frame which details the 
# creation of the analytic set and how many samples were excluded and why
colnames(inclusionTable[[1]])[ncol(inclusionTable[[1]])] <- "Mayo.eset"

# Save the data.frame to the harddrive
write.csv(inclusionTable[[1]], "1.DataInclusion/Data/Inclusions.csv")

# The second list element in inclusionTable is a list of dataset specific 'good' sample IDs
# i.e. samples which should be included before applying doppelgangR
goodSamples <- inclusionTable[[2]]

# Remove the extra TCGA data (rnaseq and mirna)
goodSamples <- goodSamples[-1 * grep("rna|RNA", names(goodSamples))]

# Only consider the esets with the minimum number of samples
esetList.chosen <- list()
goodSamples.chosen <- list()
for (i in 1:(length(goodSamples))) {
  if (length(goodSamples[[i]]) >= minimumSamples & !(names(goodSamples)[i] %in% excludeEsets)) {
    goodSamples.chosen[[names(goodSamples)[i]]] <- goodSamples[[i]]  # get the samples
    exprs <- paste("data(",names(goodSamples)[i],")",sep="")  # load the eset
    eval(parse(text=exprs))
    
    rm(exprs)
    
    # Limit it to only the good samples
    exprsString <- paste(names(goodSamples)[i], " <- ",names(goodSamples)[i],"[,goodSamples[[i]]]",sep="" )
    eval(parse(text=exprsString))
    
    # add the eset to the eset List
    exprsString <- paste("esetList.chosen[[",length(goodSamples.chosen),"]] <- ",names(goodSamples)[i], sep="")
    eval(parse(text=exprsString))
    
    # delete the eset
    exprsString <- paste("rm(",names(goodSamples)[i],")",sep="")
    eval(parse(text=exprsString))
  }
}
names(esetList.chosen) <- names(goodSamples.chosen)

# Add the Mayo data to the esetList.chosen and goodSamples.chosen
esetList.chosen[[length(esetList.chosen)+1]] <- mayo.eset[,inclusionTable.mayo[[2]]]
goodSamples.chosen[[length(esetList.chosen)]] <- inclusionTable.mayo[[2]]
names(esetList.chosen)[length(esetList.chosen)] <- names(goodSamples.chosen)[length(esetList.chosen)] <- "mayo.eset"

# Pre-process the esets to improve matching
testesets <- lapply(esetList.chosen, function(X){  
  X$alt_sample_name <- paste(X$sample_type, gsub("[^0-9]", "", X$alt_sample_name), sep="_")
  pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
  return(X) })

# Use doppelgangR to find similar sample pairs accross datasets except for mayo
# this code assumes that the Mayo data is the last eset in the list
doppelResult <- doppelgangR(testesets[-1 * length(testesets)], corFinder.args = list(use.ComBat=TRUE), 
                            cache.dir = NULL)

# Process the doppelgangR results into data.frames and write to the harddrive
doppelResult_sum <- summary(doppelResult)
doppelResult.full <- doppelResult@fullresults

# Write the doppelgangR results to the harddrive
doppelResult.full_out <- c()
for (i in 1:length(doppelResult.full)) {
  tmp <- merge(doppelResult.full[[i]][["expr.doppels"]][["outlierFinder.res"]], 
               doppelResult.full[[i]][["pheno.doppels"]][["outlierFinder.res"]], 
               by = c("sample1", "sample2"))
  
  colnames(tmp) <- c("sample1", "sample2", "expr.similarity", "expr.doppel", 
                     "pheno.similarity", "pheno.doppel")
  
  doppelResult.full_out <- rbind(doppelResult.full_out, tmp)
}

write.table(doppelResult.full_out, file = "1.DataInclusion/Data/doppelgangR/pairwiseSampleComparisons.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Write the good samples to file
for (i in 1:length(goodSamples.chosen)) {
  if (names(goodSamples.chosen)[i] != "mayo.eset") {
    # Investigate within Dataset Low Correlating Samples
    sub <- doppelResult.full_out[lapply(strsplit(doppelResult.full_out$sample1, ":"), function(x){x[1]}) %in% 
                             names(goodSamples.chosen)[i] & lapply(strsplit(doppelResult.full_out$sample2, ":"), 
                                                                   function(x){x[1]}) %in% names(goodSamples.chosen)[i], ]
    cut <- mean(sub$expr.similarity) - (2 * sd(sub$expr.similarity))
    inc <- sub[sub$expr.similarity > cut,]
    samp <- setdiff(sub$sample1, inc$sample1)
    lowcorSamples <- c()
    if (length(samp) != 0) {
      lowcorSamples <- unlist(strsplit(samp,":"))[seq(2, length(samp)*2, 2)]
      cat("Remove", length(lowcorSamples), "low-correlating samples from", names(goodSamples.chosen)[i], "\n")
    }
    
    doppelSamples <- c(sub$sample1[sub$expr.doppel & sub$expr.similarity > 0.95], 
                       sub$sample2[sub$expr.doppel & sub$expr.similarity > 0.95])
    
    if (length(doppelSamples) != 0) {
      doppelSamples <- unlist(strsplit(doppelSamples, ":"))[seq(2, length(doppelSamples)*2, 2)]
      doppelSamples <- unique(doppelSamples)
      doppelSamples <- c(doppelSamples, lowcorSamples)
      
      sampleList <- setdiff(goodSamples.chosen[[i]], doppelSamples)
      outFName <- paste("1.DataInclusion/Data/GoodSamples/", names(goodSamples.chosen)[i], 
                        "_samplesRemoved.csv", sep="")
      write.csv(sampleList, outFName)
    } else {
      sampleList <- setdiff(goodSamples.chosen[[i]], lowcorSamples)
      outFName <- paste("1.DataInclusion/Data/GoodSamples/", names(goodSamples.chosen)[i], 
                        "_samplesRemoved.csv", sep="")
      write.csv(sampleList, outFName)
    }
  } else {
    sampleList <- goodSamples.chosen[[i]]
    outFName <- paste("1.DataInclusion/Data/GoodSamples/Mayo_samplesRemoved.csv", sep="")
    write.csv(sampleList, outFName)
  }
}
