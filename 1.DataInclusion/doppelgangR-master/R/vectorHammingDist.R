vectorHammingDist <- function #Calculate Hamming Distance between two vectors, using pairwise complete observations.
### Simple function to count the fraction of different elements (in
### the same position) between two vectors of the same length, after
### removing elements from both vectors corresponding to positions
### that are NA in either vector.
(x,
### a matrix 
y,
### a matrix with the same number of columns as x
k,
### row in x to test for differences
l
### row in y to test for differences
){
    z <- as.vector(x[k,] != y[l,])
    z <- sum(z, na.rm=TRUE) / length(na.omit(z))
    return(z)
### Returns a numeric value, the Hamming Distance (the number of
### non-equal values between x and y).
}
