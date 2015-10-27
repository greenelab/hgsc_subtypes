library(affy)
efiles <- dir(pattern="^.*_esets\\.rda$")

res <- sapply(rev(efiles), function(efile){
    print(efile)
    load(efile)
    output <- c(length(esets), sum(sapply(esets, ncol)))
    names(output) <- c("n.datasets", "n.samples")
    output
})

write.csv(t(res), file="database_totalpatients.csv")
