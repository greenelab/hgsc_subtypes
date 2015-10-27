
phenoDist <- function #Calculate distance between two vectors, rows of one matrix/dataframe, or rows of two matrices/dataframes.
### This function does some simple looping to allow x and y to be
### various combinations of vectors and matrices/dataframes.
(x,
### A vector, matrix or dataframe
y=NULL,
### NULL, a vector, matrix, or dataframe.  If x is a vector, y must also be specified.
bins=10,
### discretize continuous fields in the specified number of bins
vectorDistFun=vectorWeightedDist,
### A function of two vectors that returns the distance between those vectors.
...
### Extra arguments passed on to vectorDistFun
){
    if (is.vector(x) && is.vector(y)) {
        z <- vectorDistFun(matrix(x, nrow=1), matrix(y, nrow=1), 1, 1, ...)
    }
    else {
        x <- .discretizeDataFrame(x, bins)
        if(is.null(y)){
            z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
            for (k in 1:(nrow(x) - 1)) {
                for (l in (k + 1):nrow(x)) {
                    z[k, l] <- vectorDistFun(x, x, k, l, ...)
                    z[l, k] <- z[k, l]
                }
            }
            dimnames(z) <- list(rownames(x), rownames(x))
        }else{
            y <- .discretizeDataFrame(y, bins)
            z <- matrix(0, nrow = nrow(x), ncol = nrow(y))
            for (k in 1:(nrow(x))) {
                for (l in 1:nrow(y)) {
                    z[k, l] <- vectorDistFun(x, y, k, l, ...)
                }
            }
            dimnames(z) <- list(rownames(x), rownames(y))
        }
    }
    z
### a matrix of distances between pairs of rows of x (if y is
### unspecified), or between all pairs of rows between x and y (if
### both are provided).
}


.discretizeDataFrame <- function(X, bins=10) {
    
    .discretizeRow <- function(x) {
        if (length(levels(as.factor(x))) > bins)
            return(cut(x, breaks=bins))
        as.factor(x)    
    }
    idx <- sapply(X, is.numeric)
    if (sum(idx)==0) return(X)
    X[,idx] <- as.data.frame(apply(X[,idx, drop=FALSE], 2, .discretizeRow))
    X
}

