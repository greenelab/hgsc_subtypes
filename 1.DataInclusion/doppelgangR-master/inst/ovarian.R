library(curatedOvarianData)
library(affy)
library(logging)
library(BiocParallel)
multicoreParam <- MulticoreParam()

source(system.file("extdata",
    "patientselection_all.config",package="curatedOvarianData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedOvarianData"))
rm(list=ls(pattern="_eset"))

library(doppelgangR)
save(esets, file="ovarian_esets.rda")

dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL,
                   outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
#dop <- doppelgangR(esets)
save(dop, file="ovarian_dop_1.0.rda")

load("ovarian_dop_1.0.rda")
write.csv(dop@summaryresults, file="ovarian_dop_1.0.csv")
pdf("ovarian_dop_1.0.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
