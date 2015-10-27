.outer2df <- function  #going from matrices to two columns
### Output rows contain all possible
### pairwise combinations of x and y
(x, y,
### If x is a vector and y is a vector, return all or subsets of
### outer(x, y) depending on the values of bidirectional and diag.  If
### x is a matrix and y is unspecified, return the elements of x in
### the same order as would be returned if x=rownames and y=colnames.
bidirectional=TRUE,
### If TRUE, include separate rows for elements i, j and j, i
diag=TRUE
### If TRUE, include i, i elements.
 ){
    bizarre.and.unlikely.separator=" as3a2s5df5hjnm4qwe2rxcvb "
    if(class(x) == "matrix"){
        output.samples <- x
    }else if(is.vector(x) & is.vector(y)){
        output.samples <- outer(x, y, paste, sep=bizarre.and.unlikely.separator)
    }else{
        stop("Require either x to be a matrix, or both x and y to be vectors.")
    }
    if(!bidirectional & !diag){
        output.samples <- output.samples[upper.tri(output.samples)]
    }else if(!bidirectional & diag){
        output.samples <- output.samples[!lower.tri(output.samples)]
    }else if(bidirectional & !diag){
        output.samples <- output.samples[upper.tri(output.samples) | lower.tri(output.samples)]
    }else if(bidirectional & diag){
        output.samples <- as.vector(output.samples)
    }
    if(any(grepl(bizarre.and.unlikely.separator, output.samples, fixed=TRUE))){
        output <- data.frame(do.call(rbind, strsplit(output.samples, split=bizarre.and.unlikely.separator)),
                             stringsAsFactors=FALSE)
    }else{
        output <- output.samples
    }
    return(output)
### A two-column dataframe if x and y are vectors, or a vector if x is
### a matrix and y is NULL.
}
