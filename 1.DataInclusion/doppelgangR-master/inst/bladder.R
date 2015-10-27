library(curatedBladderData)
library(affy)
library(BiocParallel)
multicoreParam <- MulticoreParam()

source(system.file("extdata",
    "patientselection_all.config",package="curatedBladderData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedBladderData"))
rm(list=ls(pattern="_eset"))

library(doppelgangR)
save(esets, file="bladder_esets.rda", compress="bzip2")

dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL,
                   outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
#dop <- doppelgangR(esets)
save(dop, file="bladder_dop_1.0.rda")

load("bladder_dop_1.0.rda")
write.csv(dop@summaryresults, file="bladder_dop_1.0.csv")

pdf("bladder_dop_1.0.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
