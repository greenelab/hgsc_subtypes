doppelgangR <- structure(function
### Identify samples with suspiciously high correlations and phenotype similarities
(esets,
### a list of ExpressionSets, containing the numeric and phenotypic data to be analyzed.
separator=":",
### a delimitor to use between dataset names and sample names
corFinder.args=list(separator=separator, use.ComBat=TRUE, method="pearson"),
### a list of arguments to be passed to the corFinder function.
phenoFinder.args=list(separator=separator, vectorDistFun=vectorWeightedDist),
### a list of arguments to be passed to the phenoFinder function.  If
### NULL, samples with similar phenotypes will not be searched for.
outlierFinder.expr.args=list(bonf.prob=0.5, transFun=atanh, tail="upper"),
### a list of arguments to be passed to outlierFinder when called for expression data
outlierFinder.pheno.args=list(normal.upper.thresh=0.99, bonf.prob=NULL, tail="upper"),
### a list of arguments to be passed to outlierFinder when called for phenotype data
smokingGunFinder.args=list(transFun=I),
### a list of arguments to be passed to smokingGunFinder
impute.knn.args=list(k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069),
### a list of arguments to be passed to impute::impute.knn.  Set to
### NULL to do no knn imputation.
manual.smokingguns=NULL,
### a character vector of phenoData columns that, if identical, will
### be considered evidence of duplication
automatic.smokingguns=FALSE,
### automatically look for "smoking guns."  If TRUE, look for
### phenotype variables that are unique to each patient in dataset 1,
### also unique to each patient in dataset 2, but contain exact
### matches between datasets 1 and 2.
within.datasets.only=FALSE,
### If TRUE, only search within each dataset for doppelgangers.
intermediate.pruning=FALSE,
### The default setting FALSE will result in output with no missing
### values, but uses extra memory because all results from the
### expression, phenotype, and smoking gun doppelganger searches must
### be saved until the end.  Setting this to TRUE will save memory for
### very large searches, but distance metrics will only be available
### if that value was identified as a doppelganger (for example,
### phenotype doppelgangers will have missing values for the
### expression and smoking gun similarity).
cache.dir="cache",
### The name of a directory in which to cache or look up results to save
### re-calculating correlations.  Set to NULL for no caching.
verbose=TRUE
### Print progress information
 ){
    ##Save input args except for esets:
    input.argnames <- ls()[-match("esets", ls())]
    input.args <- lapply(input.argnames, function(x) get(x))
    names(input.args) <- input.argnames
    input.args$esets.names <- names(esets)
    if(!is.null(cache.dir))
        dir.create(cache.dir, showWarnings=FALSE)
    if (is.null(names(esets)))
        names(esets) <- paste("dataset", 1:length(esets), sep="_")
    if(length(esets) > length(unique(names(esets))))
        names(esets) <- make.unique(names(esets))
    for (i in 1:length(esets)){
        sampleNames(esets[[i]]) <- paste(names(esets)[i], sampleNames(esets[[i]]), sep=separator)
        if(min(exprs(esets[[i]]), na.rm=TRUE) == -Inf | max(exprs(esets[[i]]), na.rm=TRUE) == Inf){
            warning(paste("Replacing -+Inf with min/max expression values for dataset", names(esets)[i]))
            exprs(esets[[i]])[exprs(esets[[i]])==-Inf] <- min(exprs(esets[[i]])[is.finite(exprs(esets[[i]]))], na.rm=TRUE)
            exprs(esets[[i]])[exprs(esets[[i]])==Inf] <- max(exprs(esets[[i]])[is.finite(exprs(esets[[i]]))], na.rm=TRUE)
        }
        if(!is.null(impute.knn.args) & any(!complete.cases(exprs(esets[[i]])))){
            ##KNN imputation
            message(paste("KNN imputation for", names(esets)[i]))
            impute.knn.args$data <- exprs(esets[[i]])
            impute.knn.output  <- do.call(impute::impute.knn, args=impute.knn.args)
            exprs(esets[[i]]) <- impute.knn.output$data
            .Random.seed <- impute.knn.output$rng.state  ##restore original RNG state
        }
    }
    ds.combns <- lapply(1:length(esets), function(i) c(i, i))
    if(!within.datasets.only)
        ds.combns <- c(ds.combns, combn(1:length(esets), 2, simplify=FALSE))
    output.full <- bplapply(ds.combns, function(ij){
        i <- ij[1]
        j <- ij[2]
        if (verbose) message(paste("Working on datasets", names(esets)[i], "and", names(esets)[j]))
        ## Create a negative result dataframe that will be used when a similarity-finder
        ## (corFinder, phenoFinder, smokingGunFinder) is not called.
        na.output <- matrix(0, nrow=ncol(esets[[i]]), ncol=ncol(esets[[j]]))
        rownames(na.output) <- sampleNames(esets[[i]])
        colnames(na.output) <- sampleNames(esets[[j]])
        if(i == j){
            na.output[!upper.tri(na.output)] <- NA
        }
        na.output <- outlierFinder(na.output, normal.upper.thresh=1)
        na.output$outlierFinder.res$similarity <- NA
        ## output object:
        output3 <- list()
        ## calculate correlation matrix
        if( sum(featureNames(esets[[i]]) %in% featureNames(esets[[j]])) < 2){
            warning(paste(names(esets)[i], "and", names(esets)[j], "have no featureNames in common, skipping corFinder for this dataset pair."))
            corFinder.args <- NULL
        }
        if(!is.null(cache.dir) & !is.null(corFinder.args)){
            cache.file <- paste(cache.dir, "/", digest::digest(list(corFinder, corFinder.args, esets[c(i, j)])), ".rda", sep="")
            if(file.exists(cache.file)) {
		if (verbose) message("\tSkipping corFinder, loading cached results.")
                load(cache.file)
            }
        }
        if(is.null(corFinder.args)){
            cor.sim <- na.output
        }else{
            corFinder.args$eset.pair <- esets[c(i, j)]
            if(!exists("cor.sim")){
                if(verbose) message("Calculating correlations...")
                cor.sim <- do.call(corFinder, corFinder.args)
            }
        }
        if(!is.null(cache.dir) && !file.exists(cache.file))
            save(cor.sim, file=cache.file)
        ## find numeric (expression) doppelgangers
        outlierFinder.expr.args$similarity.mat <- cor.sim
        if(is.null(corFinder.args)){
            output3[["expr.doppels"]] <- na.output
        }else{
            if(verbose) message("Identifying correlation doppelgangers...")
            output3[["correlations"]] <- cor.sim
            output3[["expr.doppels"]] <-
                do.call(outlierFinder, outlierFinder.expr.args)
        }
        ## If there is no phenoData in one of the esets, do not
        ## search for phenotype doppelgangers:
        if(min(c(dim(pData(esets[[1]])), dim(pData(esets[[2]])))) == 0)
            phenoFinder.args <- NULL
        ## If column names don't match, do not search for phenotype doppelgangers:
        if(!identical(colnames(pData(esets[[i]])), colnames(pData(esets[[j]])))){
            warning(paste(names(esets)[i], "and",names(esets)[j], "have different column names in phenoData.  Skipping phenotype checking for this pair.  Set phenoFinger.args=NULL to disable phenotype checking altogether."))
            phenoFinder.args <- NULL
        }
        ## automatically find potential "smoking gun" phenotypes
        if(automatic.smokingguns & !is.null(phenoFinder.args)){
            new.smokinggun.phenotypes <- unlist(sapply(colnames(pData(esets[[i]])), function(cname){
                if(cname %in% colnames(pData(esets[[j]]))){
                    if((sum(!is.na(pData(esets[[i]])[, cname])) > 2 &
                        sum(!is.na(pData(esets[[j]])[, cname])) > 2) &
                       (identical(length(pData(esets[[i]])[, cname]), length(unique(pData(esets[[i]])[, cname]))) |
                        identical(length(pData(esets[[j]])[, cname]), length(unique(pData(esets[[j]])[, cname]))))
                       )
                        return(cname)}}))
            manual.smokingguns <- unique(c(manual.smokingguns, new.smokinggun.phenotypes))
        }
##        output3[["smokingguns"]] <- manual.smokingguns  ##FIXME: write full arguments somewhere else
        if(is.null(manual.smokingguns)){
            output3[["smokinggun.doppels"]] <- na.output
        }else{
            ## find smokinggun doppelgangers
            if(verbose) message("Identifying smoking-gun doppelgangers...")
            smokingGunFinder.args$eset.pair <- esets[c(i, j)]
            smokingGunFinder.args$smokingguns <- manual.smokingguns
            outlierFinder.smokinggun.args <- list()
            outlierFinder.smokinggun.args$similarity.mat <- do.call(smokingGunFinder, smokingGunFinder.args)
            outlierFinder.smokinggun.args$normal.upper.thresh <- 0.5
            output3[["smokinggun.doppels"]] <- do.call(outlierFinder, outlierFinder.smokinggun.args)
        }
        ## calculate phenotype similarity matrix
        if(is.null(phenoFinder.args)){
            output3[["pheno.doppels"]] <- na.output
        }else{
            phenoFinder.args$eset.pair <- esets[c(i, j)]
            ##keep no rows because we only need pData(x):
            for (k in 1:2)
                phenoFinder.args$eset.pair[[k]] <- phenoFinder.args$eset.pair[[k]][0, ]
            ## do not use "smoking guns" when calculating phenotype distances
            for (k in 1:2)
                if(any(manual.smokingguns %in% colnames(pData(phenoFinder.args$eset.pair[[k]]))))
                    pData(phenoFinder.args$eset.pair[[k]]) <-
                        pData(phenoFinder.args$eset.pair[[k]])[, (!colnames(pData(phenoFinder.args$eset.pair[[k]])) %in% manual.smokingguns)]
            if(!is.null(cache.dir)){
                cache.file <- paste(cache.dir, "/", digest::digest(list(phenoFinder, phenoFinder.args)), ".rda", sep="")
                if(file.exists(cache.file)) {
		    if (verbose) message("\tSkipping phenoFinder, loading cached results.")
                    load(cache.file)
		}
            }
            if(!exists("pheno.sim")){
                if(verbose) message("Calculating phenotype similarities...")
                pheno.sim <- do.call(phenoFinder, phenoFinder.args)
            }
            if(!is.null(cache.dir) && !file.exists(cache.file))
                save(pheno.sim, file=cache.file)
            ## find phenotype doppelgangers
            outlierFinder.pheno.args$similarity.mat <- pheno.sim
            if(verbose) message("Identifying phenotype doppelgangers...")
            output3[["pheno.doppels"]] <- do.call(outlierFinder, outlierFinder.pheno.args)
        }
        ##If a sample is identified as a doppelganger by any method,
        ##then keep that sample for all methods, so we don't lose the
        ##data when wrapping up:
        if(intermediate.pruning){
            keep.rows <- output3[["pheno.doppels"]]$outlierFinder.res$doppel |
                output3[["expr.doppels"]]$outlierFinder.res$doppel |
                    output3[["smokinggun.doppels"]]$outlierFinder.res$doppel
            output3[["smokinggun.doppels"]]$outlierFinder.res <-
                output3[["smokinggun.doppels"]]$outlierFinder.res[keep.rows, ]
            output3[["pheno.doppels"]]$outlierFinder.res <-
                output3[["pheno.doppels"]]$outlierFinder.res[keep.rows, ]
            output3[["expr.doppels"]]$outlierFinder.res <-
                output3[["expr.doppels"]]$outlierFinder.res[keep.rows, ]
        }
        return(output3)
    })
    names(output.full) <- sapply(ds.combns, function(ij){
        paste(names(esets)[c(ij[1], ij[2])], collapse=separator)})
    if (verbose) message("Finalizing...")
    wrapUp <- function(object, element){
        tmp <- lapply(object, function(x) x[[element]]$outlierFinder.res)
        tmp <- tmp[!sapply(tmp, is.null)]
        do.call(rbind, tmp)
    }
    pheno.doppels <- wrapUp(output.full, "pheno.doppels")
    expr.doppels <- wrapUp(output.full, "expr.doppels")
    smokinggun.doppels <- wrapUp(output.full, "smokinggun.doppels")
    addCols <- function(orig, add){
        newrows1 <- paste(add[, 1], add[, 2])
        newrows2 <- paste(add[, 2], add[, 1])
        oldrows <- paste(orig[, 2], orig[, 1])
        if(any(!(newrows1 %in% oldrows | newrows2 %in% oldrows))){
            tmp <- add[!(newrows1 %in% oldrows | newrows2 %in% oldrows), 1:2]
            tmp <- cbind(tmp, matrix(NA, nrow=nrow(tmp), ncol=(ncol(orig) - 2)))
            colnames(tmp) <- colnames(orig)
            orig <- rbind(orig, tmp)
        }
        match.rows <- match(paste(add[, 1], add[, 2]), paste(orig[, 1], orig[, 2]))
        match.rows2 <- match(paste(add[, 2], add[, 1]), paste(orig[, 1], orig[, 2]))
        match.rows[is.na(match.rows)] <- match.rows2[is.na(match.rows)]
        if(nrow(add) > 0 & nrow(orig) == 0){
            ##adding to a zero-row dataframe
            newdf <- add[, 1:2]
            for (i in 3:ncol(orig))
                newdf[[colnames(orig)[i]]] <- NA
            orig <- newdf
        }
        if(nrow(add) == 0 & nrow(orig) == 0){
            return(cbind(orig, add[, -1:-2]))
        }
        orig[[colnames(add)[3]]] <- NA
        orig[[colnames(add)[3]]][match.rows] <- add[, 3]
        orig[[colnames(add)[4]]] <- FALSE
        orig[[colnames(add)[4]]][match.rows] <- add[, 4]
        orig
    }
    ## merge all doppelganger types
    all.doppels <- expr.doppels
    if(intermediate.pruning){
        colnames(all.doppels)[3:4] <- paste("expr.", colnames(all.doppels)[3:4], sep="")
        all.doppels <- addCols(all.doppels, pheno.doppels)
        colnames(all.doppels)[5:6] <- paste("pheno.", colnames(all.doppels)[5:6], sep="")
        all.doppels <- addCols(all.doppels, smokinggun.doppels)
        colnames(all.doppels)[7:8] <- sub("pheno", "smokinggun", colnames(all.doppels)[5:6])
    }else{
        if(identical(expr.doppels[, 1:2], pheno.doppels[, 1:2]) &
           identical(expr.doppels[, 1:2], smokinggun.doppels[, 1:2])){
            all.doppels <- cbind(all.doppels, pheno.doppels[, 3:4], smokinggun.doppels[, 3:4])
            colnames(all.doppels)[3:8] <- paste(c("expr", "expr", "pheno", "pheno",
                                                  "smokinggun", "smokinggun"),
                                                colnames(all.doppels)[3:8], sep=".")
        }else{
            stop("Intermediate pruning off but no addCols shortcut available.")
        }
    }
    ##Following 2 lines needed if pruning is done above on output3:
    for (k in c(4, 6, 8))
        all.doppels[, k][is.na(all.doppels[, k])] <- FALSE
    all.doppels <- all.doppels[all.doppels[, 4] | all.doppels[, 6] | all.doppels[, 8], ]
    rownames(all.doppels) <- NULL
    ##If all esets have the same colnames of pdata, add merged
    ##pairwise sample pdata to all.doppels:
    if( identical(nrow(all.doppels) > 0, TRUE) &&
       identical(colnames(pData(esets[[1]])),
                 unique(unlist(lapply(esets, function(eset) colnames(pData(eset)))))) ){
        all.pdat <- lapply(esets, pData)
        for (k in 1:length(all.pdat))
            all.pdat[[k]]$sampleid <- rownames(all.pdat[[k]])
        all.pdat <- do.call(rbind, all.pdat)
        pdat.sample1 <- all.pdat[match(all.doppels[, 1], all.pdat$sampleid), ]
        pdat.sample2 <- all.pdat[match(all.doppels[, 2], all.pdat$sampleid), ]
        merged.pdat <- lapply(1:ncol(pdat.sample1), function(k){
            paste(pdat.sample1[, k],
                  pdat.sample2[, k],
                  sep=separator)
        })
        merged.pdat <- do.call(cbind, merged.pdat)
        merged.pdat <- subset(merged.pdat, select=-ncol(merged.pdat))
        colnames(merged.pdat) <- colnames(pData(esets[[1]]))
        all.doppels <- cbind(all.doppels, merged.pdat)
    }
    new("DoppelGang", fullresults=output.full, summaryresults=all.doppels, inputargs=input.args)
### Returns an object of S4-class "DoppelGang".  See ?DoppelGang-class.
}, ex=function(){
        library(curatedOvarianData)
        data(GSE32062.GPL6480_eset)
        data(GSE32063_eset)
        data(GSE12470_eset)
        data(GSE17260_eset)
        testesets <- list(JapaneseA=GSE32062.GPL6480_eset,
                          JapaneseB=GSE32063_eset,
                          Yoshihara2009=GSE12470_eset,
                          Yoshihara2010=GSE17260_eset)
        testesets <- lapply(testesets, function(X){
            ## standardize the sample ids to improve matching based on clinical annotation
            X$alt_sample_name <- paste(X$sample_type, gsub("[^0-9]", "", X$alt_sample_name), sep="_")
            pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
            X <- X[, 1:20]  ##speed computations
            return(X) })
        results1 <- doppelgangR(testesets, cache.dir=NULL)
        results1
        plot(results1)
        summary(results1)
    ## Set phenoFinder.args=NULL to ignore similar phenotypes, and
    ## turn off ComBat batch correction:
##    results2 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=FALSE), phenoFinder.args=NULL, cache.dir=NULL)
##    summary(results2)
})
