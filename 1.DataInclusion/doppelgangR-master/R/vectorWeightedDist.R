vectorWeightedDist <- function #Calculate a weighted distance between two vectors, using pairwise complete observations.
### Simple function to count the fraction of different elements (in
### the same position) between two vectors of the same length, after
### removing elements from both vectors corresponding to positions
### that are NA in either vector. Distance is the probability for observing
### the matches and mismatches in two random patients.
(x,
### a matrix 
y,
### a matrix with the same number of columns as x
k,
### row in x to test for differences
l
### row in y to test for differences
){
    idx <- !(is.na(x[k,]) | is.na(y[l,]))

    if (sum(idx) < 2) return(1)

    x <- x[,idx, drop=FALSE]
    y <- y[,idx, drop=FALSE]

    p.x <- sapply(1:ncol(x), function(i) sum(x[k,i]==x[,i],
    na.rm=TRUE)/nrow(x[!is.na(x[,i]),]))
    p.y <- sapply(1:ncol(y), function(i) sum(y[l,i]==y[,i], 
    na.rm=TRUE)/nrow(y[!is.na(y[,i]),]))

    idx <- x[k,] != y[l,]
    w <- 1 - p.x * p.y
    if (sum(w) < 2 ) return(1)
    1-(sum(ifelse(idx,-1,1)*w)/sum(w)+1)/2
### Returns a numeric value, the log of the probability of observing the
### matches in x and y
}
