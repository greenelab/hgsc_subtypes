library(curatedCRCData)
library(affy)
library(BiocParallel)
multicoreParam <- MulticoreParam()

source(system.file("extdata",
    "patientselection_all.config",package="curatedCRCData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedCRCData"))
rm(list=ls(pattern="_eset"))

library(doppelgangR)
##esets <- esets[c("GSE11237_eset", "GSE14095_eset")]
##esets <- esets[c("GSE11237_eset", "GSE3294_eset")]
esets <- esets[-match(c("GSE3294_eset","TCGA.RNASeqV2_eset", "TCGA.RNASeqV2.READ_eset"), names(esets))]  ##uses Genbank IDs instead of gene symbols.
save(esets, file="CRC_esets.rda", compress="bzip2")
dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL, outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
warnings()
#dop <- doppelgangR(esets)
save(dop, file="crc_dop_1.0.rda")

load("crc_dop_1.0.rda")
write.csv(dop@summaryresults, file="crc_dop_1.0.csv")
pdf("CRC_dop_1.0.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
