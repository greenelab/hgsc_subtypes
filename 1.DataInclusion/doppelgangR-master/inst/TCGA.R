if( !require(RTCGAToolbox) ){
    library(devtools)
    install_github("mksamur/RTCGAToolbox")
    library(RTCGAToolbox)
}

all.dates <- getFirehoseRunningDates()
all.datasets <- getFirehoseDatasets()

if(file.exists("tcga.res")){
    load("/scratch/lw391/doppelgangR/inst")
}else{
    tcga.res <- list()
    for (i in 1:length(all.datasets)){
        (ds.name <- all.datasets[i])
        if(!ds.name %in% names(tcga.res)){
            res <- try(getFirehoseData(ds.name, runDate=all.dates[1], RNAseq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, mRNA_Array=TRUE))
        }
        if(!is(res, "try-error")){
            tcga.res[[ds.name]] <- res
        }
    }
    save(tcga.res, file="/scratch/lw391/doppelgangR/inst/TCGA.rda")
}else{
    load("/scratch/lw391/doppelgangR/inst")
}

extractRTCGA <- function(object, type){
    typematch <- match.arg(type,
          choices=c("RNAseq_Gene", "Clinic", "miRNASeq_Gene",
              "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
              "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
              "miRNA_Array", "RPPA"))
    if(identical(typematch, "RNAseq_Gene")){
        output <- object@RNASeqGene
    }else if(identical(typematch, "RNAseq2_Gene_Norm")){
        output <- object@RNASeq2GeneNorm
    }else if(identical(typematch, "mRNA_Array")){
        if(is(object@mRNAArray, "FirehosemRNAArray")){
            output <- object@mRNAArray@Datamatrix
        }else if(is(object@mRNAArray, "list")){
            output <- lapply(object@mRNAArray, function(tmp){
                tmp@DataMatrix
            })
            if(length(output) == 0){
                output <- matrix(NA, nrow=0, ncol=0)
            }else if(length(output) == 1){
                output <- output[[1]]
            }else{
                ## just silently take the platform with the greatest
                ## number of samples:
                keeplist <- which.max(sapply(output, ncol))
                output <- output[[keeplist]]
                warning(paste("Taking the mRNA_Array platform with the greatest number of samples:", keeplist))
            }
        }
    }else{
        stop(paste("Type", typematch, "not yet supported."))
    }
    return(output)
}

library(affy)
if(file.exists("/scratch/lw391/doppelgangR/inst/tcga.esets.rda")){
    load("/scratch/lw391/doppelgangR/inst/tcga.esets.rda")
}else{
    tcga.esets <- list()
    for (i in 1:length(tcga.res)){
        print(names(tcga.res)[i])
        tmp <- list()
        tmp[["mrna"]] <- extractRTCGA(tcga.res[[i]], "mRNA_Array")
        tmp[["rnaseq"]] <- extractRTCGA(tcga.res[[i]], "RNAseq_Gene")
        tmp[["rnaseq2"]] <- extractRTCGA(tcga.res[[i]], "RNAseq2_Gene_Norm")
        pickplat <- which.max(sapply(tmp, ncol))
        tcga.esets[[paste(names(tcga.res)[i], names(tmp)[pickplat])]] <- ExpressionSet(tmp[[pickplat]])
    }
    save(tcga.esets, file="/scratch/lw391/doppelgangR/inst/tcga.esets.rda")
}

cor.list <- lapply(tcga.esets, function(eset){
    output <- cor(exprs(eset))
    output[upper.tri(output)]
})
names(cor.list) <- names(tcga.esets)

save(cor.list, file="/scratch/lw391/doppelgangR/inst/cor.list.rda")

load("/scratch/lw391/doppelgangR/inst/cor.list.rda")

par(ask=TRUE)
for (i in 1:length(cor.list)){
    hist(cor.list[[i]], main=names(cor.list)[i], breaks="FD")
    abline(v=0.99, col="red"); abline(v=0.95, col="red")
}

100*sort(sapply(cor.list, function(x) sum(x > 0.95) / length(x)))
sort(sapply(cor.list, median))
sort(sapply(cor.list, quantile, 0.999))
     
##bimodal: KICH (RNAseq2), maybe KIRC and KIRP, LUSC, PAAD, PCPG, THCA, GBM

##very high correlations normally: THCA, HNSC, LIHC, maybe KIRP, LAML, PCPG, PRAD, STAD, THYM, ESCA

##looks good: ACC, BLCA, BRCA, CESC, COAD, COADREAD, DLBC, LGG, 
