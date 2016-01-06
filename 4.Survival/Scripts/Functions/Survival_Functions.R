############################################
# Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script stores several functions that are required for survival analyses

############################################
#Load Functions
############################################
#pData is a list of full phenodata for each dataset
#clusterMemb is a list of cluster memberships for each sample in each dataset

GetCoxPHready <- function (pData, clusterMemb, recodeAge = F) {
  # ~~~~~~~~~~~~~~
  # Retrieves and subsets pertinent covariates and samples to be used in building
  # the survival models. This function will output dataframes to be used in coxph.
  # Note: the pData list and clusterMemb list are required to be in the same order.
  #
  # Args: 
  # pData: phenotype data for a list of datasets
  # clusterMemb: cluster membership information for a list of datasets
  # recodeAge: a boolean expression indicating whether or not age should be binned 
  #
  # Returns:
  # A list object holding the important covariate data for all samples with full data
  # ~~~~~~~~~~~~~~
  
  # Load PhenoData it is in a list; first step is to subset to only the covariates we need
  PhenoDataSubset <- list()
  for (pheno in 1:length(pData)) {
    # Covariates need to be recoded for the Mayo clinic data
    if(names(pData)[pheno] == "Mayo") {
      phenoData <- pData[[pheno]]
      
      # Compile columns of use, combine them, and then rename them
      # Age column
      age <- phenoData$age_at_initial_pathologic_diagnosis
      
      # Bin the ages into 10 bins if recodAge is True
      if (recodeAge == T) {
        age <- recode(age, "0:14=1; 15:39=2; 40:44=3; 45:49=4; 50:54=5; 55:59=6; 60:64=7; 
                      65:69=8; 70:74=9; 75:150=10")
      }
      
      # Rename other important variables
      stage <- phenoData$tumorstage
      grade <- phenoData$grade
      debulking <- phenoData$debulking
      vital <- phenoData$vital_status
      days <- phenoData$days_to_death
      
      # Convert days into months 
      days <- as.numeric(paste(days)) / 30.4375
      
      # Combine important covariate data
      usePheno <- cbind(age, stage, grade, debulking, vital, days)
      colnames(usePheno) <- c("age_at_initial_pathologic_diagnosis", "tumorstage", "grade", 
                              "debulking", "vital_status", "days_to_death")
      rownames(usePheno) <- phenoData$ID
      
      # Prepare the vital status codes in advance (0 = alive; 1 = deceased)    
      usePheno <- as.data.frame(usePheno)
      usePheno$vital_status <- gsub(1, 0, usePheno$vital_status)
      usePheno$vital_status <- gsub(2, 1, usePheno$vital_status)
      usePheno$vital_status <- as.numeric(paste(usePheno$vital_status))
      PhenoDataSubset[[names(pData)[pheno]]] <- usePheno
 
    } else {
      phenoData <- pData[[pheno]]
      usePheno <- phenoData[ ,c("age_at_initial_pathologic_diagnosis", "tumorstage", "grade", 
                               "debulking", "vital_status", "days_to_death")]
      colnames(usePheno) <- c("age_at_initial_pathologic_diagnosis", "tumorstage", "grade", 
                              "debulking", "vital_status", "days_to_death")
      if(recodeAge == T) {
        recoding <- recode(usePheno[ ,"age_at_initial_pathologic_diagnosis"], "0:14=1; 15:39=2; 
                           40:44=3; 45:49=4; 50:54=5; 55:59=6; 60:64=7; 65:69=8; 70:74=9; 75:150=10")
        usePheno[ ,"age_at_initial_pathologic_diagnosis"] <- recoding
      }
      
      # Recode vital status
      usePheno$vital_status <- gsub("living", 0, usePheno$vital_status)
      usePheno$vital_status <- gsub("deceased", 1, usePheno$vital_status)
      usePheno$vital_status <- as.numeric(paste(usePheno$vital_status))
      rownames(usePheno) <- rownames(phenoData)
      
      # Convert days into months 
      usePheno$days_to_death <- as.numeric(paste(usePheno$days_to_death)) / 30.4375
      
      # Store to list
      PhenoDataSubset[[names(pData)[pheno]]] <- usePheno
    }
  }
  
  # Create a list of dataframes that includes phenodata and cluster membership data
  CoxDataFrameList <- list()
  for (clus in 1:length(clusterMemb)) {
    
    # Get cluster membership and PhenoData for each dataset
    Clusters <- clusterMemb[[clus]]
    PhenoData <- PhenoDataSubset[[clus]]
    
    # Combine and store in new dataframe
    coxUse <- cbind(PhenoData, Clusters)
    CoxDataFrameList[[names(PhenoDataSubset)[clus]]] <- coxUse
  }
  
  # Run the cox proportional hazards model correcting for covariates, when they are available
  CoxPH <- list()
  for (coxmodel in 1:length(CoxDataFrameList)) {
    
    # compile a list of useable covariates for each dataset in CoxDataFrameList
    covariates <- c()
    for(covar in 1:ncol(CoxDataFrameList[[coxmodel]]))
    {
      
      # Covariate information for the particular cox model
      tmp <- CoxDataFrameList[[coxmodel]][ ,covar]
      if (sum(is.na(tmp)) == length(tmp)) {
        
        # If there is no information for all samples for the given factor, print info to screen
        cat(names(CoxDataFrameList)[coxmodel], "do not use this factor:", 
            colnames(CoxDataFrameList[[coxmodel]])[covar], "\n")
      } else {
        tmpfactor <- as.character(paste(colnames(CoxDataFrameList[[coxmodel]])[covar]))
        covariates <- c(covariates, tmpfactor)
      }
    }
    
    # Subset the given model to useable covariates and print to screen which ones are being used
    CoxDataFrameList[[coxmodel]] <- CoxDataFrameList[[coxmodel]][ ,covariates]
    cat(names(PhenoDataSubset)[coxmodel], "use these covariates:", covariates, "\n")
    
    # Here, you must subset usecox to only the complete.cases (subset samples)
    usecox <- CoxDataFrameList[[coxmodel]]
    usecox <- usecox[complete.cases(usecox), ] 
    
    # Correct for "unknown" coding for debulking status in Mayo
    if (names(CoxDataFrameList)[coxmodel] == "Mayo") {
      usecox <- usecox[usecox$debulking == 1 | usecox$debulking == 2, ]
    }
    
    # Store to cox model list and prepare for output
    CoxPH[[names(CoxDataFrameList)[coxmodel]]] <- usecox
  }
  return(CoxPH)
}


doCoxPH_KM <- function (coxphdata, fname) {
  # ~~~~~~~~~~~~~~
  # This function will write the KM plots to file
  #
  # Args: 
  # coxphdata: a sample by covariate dataframe
  # fname: the name of the dataset
  #
  # Returns:
  # A Kaplan Meier curve representing survival for each cluster
  # ~~~~~~~~~~~~~~
  
  # Get a survfit object and save it to a list
  tmpsurvfit3 <- survfit(Surv(days_to_death, vital_status) ~ ClusterK3, data = coxphdata)
  tmpsurvfit4 <- survfit(Surv(days_to_death, vital_status) ~ ClusterK4, data = coxphdata)
  survlist <- list(K3 = tmpsurvfit3, K4 = tmpsurvfit4)
  
  # Store sample size variable
  n <- nrow(coxphdata)
  
  # Initialize graphics device
  png(paste("4.Survival/Figures/", fname, "KM_survival.png", sep = ""), width = 525, height = 900)
  
  # Obtain appropriate plotting margins
  par(mfrow=c(2,1))
  par(mar = c(5,5,5,2) + 0.1)
  
  # K = 3 plot and legend
  plot(survlist[[1]], main = "", col = c("blue", "red", "green"), lwd = 3, cex = 2, cex.lab = 2, 
       cex.axis = 1.75, lty = c(4, 5), xlab = "", ylab = "") 
  legend("topright", legend = paste("Survival\nn = ", n), bty = "n", cex = 2)
  
  # K = 4 plot and legend
  plot(survlist[[2]], main = "", col = c("blue", "red", "green", "purple"), lwd = 3, cex = 2, 
       cex.lab = 2, cex.axis = 1.75, lty = c(4,5), xlab = "", ylab = "") 
  legend("topright", legend = paste("Survival\nn = ", n), bty = "n", cex = 2)
  dev.off()
}


customCoxPH <- function (data, type = "multi") {
  # ~~~~~~~~~~~~~~
  # This function will prepare the cox models and evaluate them
  #
  # Args: 
  # data: a sample by covariate dataframe
  # type: the completeness by which the cox model is built. Can be either
  # "multi", "uni", or "removeAge"
  #
  # Returns:
  # The results of a Cox proportional hazards model separated by clustering methods
  # ~~~~~~~~~~~~~~
  
  # Get the covariate variables
  variables <- colnames(data)
  
  # Get the important days to death and vital status variables
  dtd <- variables[grepl("days_to_death", variables)]
  vs <- variables[grepl("vital_status", variables)]
  
  # Get cluster information
  clus <- variables[grepl("Cluster", variables)]
  
  # Separate the variables into covariates
  covariates <- setdiff(variables, c(dtd, vs, clus))
  
  # Build the cox models for each cluster method
  coxPH_list <- list()
  for (centroid in 1:length(clus)) {
    if(type == "multi") {
      
      # Initialize the model; this will be different according to the covariates available
      model <- c()
      for (covar in 1:(length(covariates) + 1)) {
        if (covar == 1) {
          model <- paste("coxph(Surv(", dtd, ", ", vs, ")~ as.factor(", clus[centroid], ") + ", 
                         covariates[covar], sep = "")
        } else if (covar <= length(covariates)) {
          model <- paste(model, covariates[covar], sep = " + ")
        } else {
          model <- paste(model, ", data = data)", sep = "")
        }     
      }
      
      # Evaluate the model and store into list
      coxPH_list[[clus[centroid]]] <- eval(parse(text = model))

    # Get the univariate model information
    } else if (type == "uni") {
      model <- paste("coxph(Surv(", dtd, ", ", vs, ")~ as.factor(", clus[centroid], "), data = data)")
      coxPH_list[[clus[centroid]]]  <- eval(parse(text = model))

    # Remove age adjustments for all models
    } else if (type == "removeAge") {
      model <- c()
      if ("age_at_initial_pathologic_diagnosis" %in% covariates) {
        covariates <- covariates[-grep("age_at_initial_pathologic_diagnosis", covariates)]
      }
      for (covar in 1:(length(covariates) + 1)) {
        
        if (covar == 1) {
          model <- paste("coxph(Surv(", dtd, ", ", vs, ")~ as.factor(", clus[centroid], ") + ", 
                         covariates[covar], sep = "")
        } else if (covar <= length(covariates)) {
          model <- paste(model, covariates[covar], sep = " + ")
        } else {
          model <- paste(model, ", data = data)", sep = "")
        }   
      }
      
      # Evaluate the model
      coxPH_list[[clus[centroid]]] <- eval(parse(text = model))
    }
  }
  return(coxPH_list)
}


coxSum <- function (coxData, name) {
  # ~~~~~~~~~~~~~~
  # This function will write out the hazards ratios, confidence intervals, pvalues,
  # and Wald's P for each cox proportional hazards model
  #
  # Args: 
  # coxData: a coxph object
  # name: the name of the dataset
  #
  # Returns:
  # a dataframe of pertinent information regarding the cox model summary
  # ~~~~~~~~~~~~~~
  
  cox.sum <- summary(coxData)
  
  # 95% Confidence Intervals
  conf.int <- cox.sum$conf.int[ ,3:4]
  
  # Hazard ratio
  hazard <- exp(coxData$coefficients)
  
  # p value
  Pvalues <- cox.sum$coefficients[ ,5]
  
  # Wald's P
  WaldsP <- cox.sum$waldtest["pvalue"]
  
  # Combine together and write to file
  CoxPHtmp <- cbind(conf.int, hazard, Pvalues, WaldsP)
  write.csv(CoxPHtmp, paste("4.Survival/Tables/", name, ".csv", sep = ""), row.names = T)
}
