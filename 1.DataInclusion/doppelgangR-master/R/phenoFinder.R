phenoFinder <- function # Calculate pairwise similarities of phenoData between samples for a list containing two ExpressionSets
### This function acts as a wrapper to phenoDist to handle cases of
### one ExpressionSet, a list of two identical ExpressionSets, or a
### list of two different ExpressionSets.
(eset.pair,
### input: a list of ExpressionSets with two elements, or an
### ExpressionSet.  If the two elements are identical, return the
### correlation matrix for pairs of samples in the first element.  If
### not identical, return pairs between the two elements.
separator=":",
### a separator between dataset name (taken from the list names) and
### sample name (taken from sampleNames(eset), to keep track of which
### samples come from which dataset.
...
### Extra arguments passed on to phenoDist
){
    if((!is(eset.pair, "list") | length(eset.pair) != 2))
        stop("eset.pair should be a list of length 2")
    if(!identical(colnames(pData(eset.pair[[1]])), colnames(pData(eset.pair[[2]]))))
        stop("pData slots of esets must have identical column names, e.g.
identical(colnames(pData(eset[[1]])), colnames(pData(eset[[2]])))")
    matrix.one <- as.matrix(pData(eset.pair[[1]]))
    if(is.null(rownames(matrix.one)))
        rownames(matrix.one) <- make.names(1:nrow(matrix.one))
    ##This part removes columns that are all NA.  The complication
    ##with keep.col is necessary because Surv objects get turned
    ##into two elements in keep.col.
    keep.col <- apply(matrix.one, 2, function(x) !all(is.na(x)))
    keep.col <- keep.col[names(keep.col) %in% colnames(matrix.one)]
    matrix.one <- subset(matrix.one, select=match(names(keep.col), colnames(matrix.one)))
    matrix.one <- subset(matrix.one, select=keep.col)
    if( identical(all.equal(pData(eset.pair[[1]]), pData(eset.pair[[2]]) ), TRUE)){
        ##Calculate similarity matrix for a single ExpressionSet:
        similarity.mat <- 1 - phenoDist(matrix.one, ...)
        similarity.mat[!upper.tri(similarity.mat)] <- NA  ##NA for all but upper triangle.
    }else{
        ##Calculate similarity matrix for two distinct ExpressionSets:
        matrix.two <- as.matrix(pData(eset.pair[[2]]))
        matrix.two <- matrix.two[, match(colnames(matrix.one), colnames(matrix.two))]
        if(is.null(rownames(matrix.two)))
            rownames(matrix.two) <- make.names(1:nrow(matrix.two))
        similarity.mat <- 1 - phenoDist(matrix.one, matrix.two, ...)
    }
    return(similarity.mat)
### A matrix of similarities between the phenotypes of pairs of samples.
}
