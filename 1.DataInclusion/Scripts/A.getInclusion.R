############################################
# Amy Campbell
# Updated from Greg Way/ James Rudd pipeline
# Cross-population analysis of high-grade serous ovarian cancer
# does not support four subtypes
#
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will perform our inclusion criteria on all datasets included
# in curatedOvarianData and a dataset from the Mayo Clinic

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
vars <- c("sample_type", "histological_type", "grade", "primarysite",
          "arrayedsite", "summarystage", "tumorstage", "substage",
          "pltx", "tax", "neo", "recurrence_status", "vital_status",
          "os_binary", "relapse_binary", "site_of_tumor_first_recurrence",
          "primary_therapy_outcome_success", "debulking")

minimumSamples <- 100

# Esets that should never be included: Dressman and Bentink
excludeEsets <- c("PMID17290060_eset", "E.MTAB.386_eset")

####################################
# Load Data
####################################
# All the datasets within the curatedOvarianData package
esets <- getAllDataSets("curatedOvarianData")

# Load the Konecny data from GEO
mayo.eset <- GEOquery::getGEO("GSE74357", getGPL = FALSE)
mayo.eset <- mayo.eset[[1]]
 
# Load the AACES expression data
aaces.exprs <- readr::read_tsv("expression.tsv")
aaces.eset <- ExpressionSet(assayData = as.matrix(aaces.exprs))


##################################
# ANALYSIS
##################################
# Process mayo.eset to include only high-grade serous ovarian tumor samples
mayo.pheno <- readr::read_csv(file.path("1.DataInclusion", "Data", "Mayo",
                                        "Mayo_Pheno_Data.csv"))
p <- pData(mayo.eset)

mapper <- data.frame(p$geo_accession)
mapper["unique_patient_ID"] <-
  unlist(lapply(p$title, function(x) strsplit(toString(x), split = "  ")[[1]][2]))

pheno.map <- dplyr::inner_join(mayo.pheno, mapper, by = "unique_patient_ID")
pheno.map.inclusion <- dplyr::filter(pheno.map, histological_type == "ser",
                                     grade == 3)

included.geo <- unlist(pheno.map.inclusion["p.geo_accession"])
mayo.eset <- mayo.eset[, sampleNames(mayo.eset) %in% included.geo]


# Use the inclusion/exclusion decision tree to filter samples in all
# curatedOvarainData datasets
inclusionTable <- exclusionTable(esets)

# Use the inclusion/exclusion decision tree on the Mayo data
inclusionTable.mayo <- sampleNames(mayo.eset)

# Use the inclusion/exclusion decision tree on aaces data
inclusionTable.aaces <- sampleNames(aaces.eset)

# Combine the inclusion results from curatedOvarianData and Mayo

inclusionTable[[1]] <- cbind(inclusionTable[[1]],
                             inclusionTable.mayo[[1]],
                             inclusionTable.aaces[[1]])

# The first list element in the inclusionTable is the data.frame which details
# the creation of the analytic set and how many samples were excluded and why

colnames(inclusionTable[[1]])[(ncol(inclusionTable[[1]]) - 1)] <- "Mayo.eset"
colnames(inclusionTable[[1]])[(ncol(inclusionTable[[1]]))] <- "aaces.eset"

# Save the data.frame to the harddrive
write.csv(inclusionTable[[1]], "1.DataInclusion/Data/Inclusions.csv")

# The second list element in inclusionTable is a list of dataset specific
# 'good' sample IDs i.e. samples which should be included before applying
# doppelgangR
goodSamples <- inclusionTable[[2]]

# Remove the extra TCGA data (rnaseq and mirna)
goodSamples <- goodSamples[-1 * grep("rna|RNA", names(goodSamples))]

# Only consider the esets with the minimum number of samples
esetList.chosen <- list()
goodSamples.chosen <- list()

for (i in 1:(length(goodSamples))) {
    if (length(goodSamples[[i]]) >= minimumSamples &
        !(names(goodSamples)[i] %in% excludeEsets)) {
    # get the samples
    goodSamples.chosen[[names(goodSamples)[i]]] <- goodSamples[[i]]
    # load the eset
    exprs <- paste("data(", names(goodSamples)[i] , ")", sep = "")

    eval(parse(text = exprs))

    rm(exprs)

    # Limit it to only the good samples
    exprsString <- paste(names(goodSamples)[i], " <- ",
                         names(goodSamples)[i], "[,goodSamples[[i]]]", sep = "")

    eval(parse(text = exprsString))

    # add the eset to the eset List
    exprsString <- 
      paste("esetList.chosen[[", length(goodSamples.chosen), "]] <- ",
            names(goodSamples)[i], sep = "")
    eval(parse(text = exprsString))

    # delete the eset
    exprsString <- paste("rm(", names(goodSamples)[i], ")", sep = "")

    eval(parse(text = exprsString))
  }
}

names(esetList.chosen) <- names(goodSamples.chosen)


esetList.chosen[[length(esetList.chosen) + 1]] <- mayo.eset
esetList.chosen[[length(esetList.chosen) + 1]] <- aaces.eset

goodSamples.chosen[[(length(esetList.chosen) - 1)]] <- sampleNames(mayo.eset)
goodSamples.chosen[[length(esetList.chosen)]] <- sampleNames(aaces.eset)


names(esetList.chosen)[(length(esetList.chosen) - 1)] <- 
  names(goodSamples.chosen)[(length(esetList.chosen) - 1)] <- "mayo.eset"
names(esetList.chosen)[length(esetList.chosen)] <-
  names(goodSamples.chosen)[length(esetList.chosen)] <- "aaces.eset"

testesets <- esetList.chosen
testesets[1:(length(esetList.chosen) - 2)] <- 
  lapply(testesets[1:(length(testesets) - 2)],
         function(X) {X$alt_sample_name <-
                      paste(X$sample_type,
                            gsub("[^0-9]", "", X$alt_sample_name), 
                            sep = "_")
                    pData(X) <- pData(X)[, !grepl("uncurated_author_metadata",
                                                  colnames(pData(X)))]
                    return(X)
                    })


doppel.result <-
  doppelgangR::doppelgangR(testesets, corFinder.args = list(use.ComBat = TRUE),
                           cache.dir = NULL)


# Process the doppelgangR results into data.frames and write to the harddrive

doppelResult.full <- summary(doppel.result)
doppelResult.full_out <-
  doppelResult.full[c("sample1", "sample2",
                      "expr.similarity", "expr.doppel",
                      "pheno.similarity", "pheno.doppel")]
doppel.fname <-
  "1.DataInclusion/Data/doppelgangR/pairwiseSampleComparisons.tsv"
write.table(doppelResult.full_out,
            file = doppel.fname,
            sep = "\t", quote = FALSE, row.names = FALSE)

for (i in 1:length(goodSamples.chosen)) {
    sub <-
      (doppelResult.full_out)[
        lapply(strsplit(doppelResult.full_out$sample1, ":"),
               function(x) {x[1]}) %in% names(goodSamples.chosen)[i] &
          lapply(strsplit(doppelResult.full_out$sample2, ":"),
                 function(x) {x[1]}) %in% names(goodSamples.chosen)[i], ]
    
    cut <- mean(sub$expr.similarity) - (2 * sd(sub$expr.similarity))
    inc <- sub[sub$expr.similarity > cut, ]
    samp <- setdiff(sub$sample1, inc$sample1)
    lowcorSamples <- c()
    if (length(samp) != 0) {
      lowcorSamples <- unlist(strsplit(samp, ":"))[seq(2, length(samp) * 2, 2)]
      cat("Remove", length(lowcorSamples), "low-correlating samples from",
          names(goodSamples.chosen)[i], "\n")
    }

    doppelSamples <-
      c(sub$sample1[sub$expr.doppel & sub$expr.similarity > 0.95],
                    sub$sample2[sub$expr.doppel & sub$expr.similarity > 0.95])

    if (length(doppelSamples) != 0) {
      doppelSamples <- 
        unlist(strsplit(doppelSamples, ":"))[seq(2, length(doppelSamples) * 2, 2)]
      doppelSamples <- unique(doppelSamples)
      doppelSamples <- c(doppelSamples, lowcorSamples)

      sampleList <- setdiff(goodSamples.chosen[[i]], doppelSamples)
      outFName <- paste("1.DataInclusion/Data/GoodSamples/",
                        names(goodSamples.chosen)[i],
                        "_samplesRemoved.csv", sep = "")
      write.csv(sampleList, outFName)
    } else {
      sampleList <- setdiff(goodSamples.chosen[[i]], lowcorSamples)
      outFName <- paste("1.DataInclusion/Data/GoodSamples/",
                        names(goodSamples.chosen)[i],
                        "_samplesRemoved.csv", sep = "")
      write.csv(sampleList, outFName)
    }
    }
