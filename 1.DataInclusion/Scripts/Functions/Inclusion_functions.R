############################################
# Cross-population analysis of high-grade serous ovarian cancer does not support four subtypes
#
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# The following script contains all functions required for the inclusion and doppelgangR script

#####################################
# FUNCTIONS
#####################################
# This function will load all of the data in curatedOvarianData
getAllDataSets <- function (packageName) {
  tmp <- data(package=packageName)
  datasets <- tmp[3][[1]][ ,3]
  return(datasets)
}

# This function reformats/renames the input 
# sample's alt_sample_name, sample_type
rename.esets <- function(sample.info) {
  sample.info$alt_sample_name <-
    paste(sample.info$sample_type,
          gsub("[^0-9]", "", sample.info$alt_sample_name), 
          sep = "_")
  pData(sample.info) <- pData(sample.info)[, !grepl("uncurated_author_metadata",
                                                    colnames(pData(sample.info)))]
  return(sample.info)
}

# This function will output an exclusion table of all the good samples
exclusionTable <- function(esetVec) {
  outTable <- data.frame()
  sampleList <- list()
  for (i in 1:length(esetVec)) {
    # Load the expression data for each eset vector
    exprs <- paste("data(", esetVec[i], ")", sep="")
    eval(parse(text=exprs))
    
    # Load the phenotype data for each eset vector
    exprs <- paste("pheno <- pData(", esetVec[i], ")", sep="")
    eval(parse(text=exprs))
    
    # Here are the names of the samples
    goodSamples <- rownames(pheno)
    nTumor <- nSerous <- nEndo <- nHighSerous <- nLowSerous <- nHighEndo <- 0
    
    # how many tumors>
    nTumor <- sum(pheno$sample_type == "tumor", na.rm=T)
    
    # reassign good samples to only tumors
    goodSamples <- rownames(pheno)[pheno$sample_type == "tumor" %in% TRUE]
    
    if (sum(is.na(pheno$histological_type), na.rm=T) < nrow(pheno)) {
      # histological type must be assigned
      nSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser", na.rm=T)
      nEndo <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "endo", na.rm=T)
      
      # Good samples are serous or endo tumors
      goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser") | (pheno$sample_type == "tumor" & pheno$histological_type == "endo"))  %in% TRUE]
      
      # How many samples are missing grade assignments
      nMissingGrade <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & is.na(pheno$grade), na.rm=T)
      
      if(sum(is.na(pheno$grade), na.rm=T) < nrow(pheno)) {
        # grade must be assigned
        nHighSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1, na.rm=T)
        nLowSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade == 1, na.rm=T)
        nHighEndo <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2, na.rm=T)
        
        # Good samples are high grade serous or high grade endo
        goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1) | (pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2)) %in% TRUE]
      }
    }
    totSamples <- nrow(pheno)
    
    sampleList[[i]] <- NA
    if (length(goodSamples) > 1) {
      sampleList[[i]] <- goodSamples
    }
    
    # Build the column for each eset
    thisInclusion <- c(totSamples, nTumor, nSerous, nMissingGrade, nHighSerous, nLowSerous, 
                       nEndo, nHighEndo, length(goodSamples))
    if(ncol(outTable) < 1) {
      outTable <- data.frame(X = thisInclusion)
    } else {
      outTable <- cbind(outTable, thisInclusion)
    }
    
    exprs <- paste("rm(",esetVec[i], ", envir=.GlobalEnv)", sep = "")
    eval(parse(text = exprs))
  }
  # name the outTable
  colnames(outTable) <- esetVec
  rownames(outTable) <- c("TotalSamples", "Tumor", "Serous", "Missing Grade","HighGradeSerous", 
                          "LowGradeSerous", "Endo", "HighGradeEndo", "AnalyticSet")
  
  # Name the list of good samples
  names(sampleList) <- esetVec
  output <- list(outTable, sampleList)
  return(output)
}

# This function performs the exclusion step exclusively on the Konecny eset data
simpleExclusion <- function (eset) {
  
  pheno <- pData(eset)
  
  goodSamples <- rownames(pheno)
  nTumor <- nSerous <- nEndo <- nHighSerous <- nLowSerous <- nHighEndo <- 0
  
  
  nTumor <- sum(pheno$sample_type == "tumor", na.rm=T)
  goodSamples <- rownames(pheno)[pheno$sample_type == "tumor" %in% TRUE]
  
  if (sum(is.na(pheno$histological_type), na.rm=T) < nrow(pheno)) {
    # histological type must be assigned
    nSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser", na.rm=T)
    nEndo <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "endo", na.rm=T)
    
    goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser") | (pheno$sample_type == "tumor" & pheno$histological_type == "endo"))  %in% TRUE]
    nMissingGrade <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & is.na(pheno$grade), na.rm=T)
    
    if (sum(is.na(pheno$grade), na.rm=T) < nrow(pheno)) {
      # grade must be assigned
      nHighSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1, na.rm=T)
      nLowSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade == 1, na.rm=T)
      nHighEndo <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2, na.rm=T)
      goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1) | (pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2)) %in% TRUE]
    }
  }
  totSamples <- nrow(pheno)

  sampleList <- NA
  if (length(goodSamples) > 1) {
    sampleList <- goodSamples
  }
  
  thisInclusion <- c(totSamples, nTumor, nSerous, nMissingGrade, nHighSerous, nLowSerous, nEndo, nHighEndo, 
                     length(goodSamples))
  outTable <- data.frame(X=thisInclusion)
  
  result <- list()
  result[[1]] <- outTable
  result[[2]] <- sampleList
  return(result)
}
