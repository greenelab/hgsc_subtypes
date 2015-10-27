library(gdata)

snames <- paste(c("bladder", "CRC", "breast", "ovarian"), "_dop.csv", sep="")

if(file.exists("manual_duplicates_curation.rda")){
    load("manual_duplicates_curation.rda")
}else{
    manuals <- lapply(snames, function(x) read.xls("manual_duplicates_curation.xls", sheet=x, as.is=TRUE))
    names(manuals) <- snames
    save(manuals, file="manual_duplicates_curation.rda")
}

autos <- lapply(snames, function(x) read.csv(x, as.is=TRUE))
names(autos) <- snames

output <- sapply(1:length(autos), function(i){
    auto.ids <- with(autos[[i]], paste(sample1, sample2))
    if(grepl("breast", names(manuals)[i])){
        is.healthy <- rep(FALSE, nrow(manuals[[i]]))
    }else{
        is.healthy <- grepl("healthy|adjacentnormal", manuals[[i]]$sample_type)
    }
    manual.ids.healthy <- with(manuals[[i]][is.healthy & manuals[[i]]$manually.confirmed, ], paste(sample1, sample2))
    manual.ids.tumor <- with(manuals[[i]][!is.healthy & manuals[[i]]$manually.confirmed, ], paste(sample1, sample2))
    ##Dataframe for counting within and between-dataset:
    manual.ds <- strsplit(manual.ids.tumor, split=" ")
    manual.ds <- t(sapply(manual.ds, function(x) sub(":.*$", "", x)))
    auto.ds <- strsplit(auto.ids, split=" ")
    auto.ds <- t(sapply(auto.ds, function(x) sub(":.*$", "", x)))
    df <- data.frame(truepos=sum(auto.ids %in% manual.ids.tumor),
                     falsepos=sum(!auto.ids %in% manual.ids.tumor),
                     falsepos.is.healthy=sum(!auto.ids %in% manual.ids.healthy),
                     totalautopos=length(auto.ids),
                     within.ds.doppels=sum(manual.ds[, 1] == manual.ds[, 2]),
                     between.ds.doppels=sum(manual.ds[, 1] != manual.ds[, 2]),
                     total.dspairs.with.dops=length(unique(paste(manual.ds[, 1], manual.ds[, 2]))),
                     total.dspairs.found=sum(unique(paste(manual.ds[, 1], manual.ds[, 2])) %in% paste(auto.ds[, 1], auto.ds[, 2])),
                     total.manual.pos=length(manual.ids.tumor)
                     )
    rownames(df) <- names(autos)[i]
    df
})
colnames(output) <- names(autos)
output <- t(output)
output

write.csv(output, file="doppelgangers_summary.csv")
