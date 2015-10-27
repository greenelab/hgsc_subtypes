library(BiocParallel)
multicoreParam <- MulticoreParam()

breast.packages <- c("breastCancerMAINZ", "breastCancerNKI", "breastCancerTRANSBIG", "breastCancerUNT", "breastCancerUPP", "breastCancerVDX")

other.packages <- "WGCNA"

if (!require(BiocInstaller))
    stop("You need to install Bioconductor, which includes BiocInstaller.")


for (pkg in breast.packages){
    if(!require(package=pkg, character.only=TRUE)){
        print(paste("Need to install", pkg))
        biocLite(pkg, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
    }
}


for (pkg in other.packages){
    if(!require(package=pkg, character.only=TRUE)){
        print(paste("Need to install", pkg))
        biocLite(pkg, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
    }
}



esets <- bplapply(breast.packages, function(pkg){
    print(pkg)
    library(affy)
    esetname <- tolower(sub("breastCancer", "", pkg))
    data(list=esetname)
    output <- get(esetname)
    output <- output[!is.na(featureData(output)$EntrezGene.ID), ]
    merge.probeset <- WGCNA::collapseRows(datET=exprs(output),
                                          rowGroup=featureData(output)$EntrezGene.ID,
                                          rowID=featureNames(output))
    output <- output[merge.probeset$selectedRow, ]
    featureNames(output) <- featureData(output)$EntrezGene.ID
    return(output)
})
names(esets) <- sub("breastCancer", "", breast.packages)

save(esets, file="breast_esets.rda", compress="bzip2")
##load("esets.rda")

#set.seed(1)
#eset2 <- lapply(esets, function(x) x[sample(1:nrow(x), 300), ])

library(doppelgangR)

dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL,
                   outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
save(dop, file="breast_dop_1.0.rda")
##load("breast_dop.rda")

write.csv(dop@summaryresults, file="breast_dop_1.0.csv")
pdf("breast_dop_1.0.pdf")
plot(dop)
dev.off()

