######################################################## 
##  Program Name  :/projects/epi/ovarian/s103108.goode/projects/MethylationProject/expression/rpgm/Agilent1and2and3_COMBAT_datamerge_forGO.r
##  Study Title   : OVCA Agilent expression tumor data
##  Programmer    : Sebastian Armasu
##  Investigator  : Dr. Ellen Goode
##  Statistician  : Dr. Brooke Fridley, Dr. Chen Wang, Sebastian Armasu
##  Date Created  : Thursday, 17 April 2014 02:13 PM CDT
##  Last Modified : Monday, 14 September 2015 05:02 PM CDT
##  Study Number  : s103108.goode
##  Function      : adjust for batch effect using COMBAT, applied on merged 3 sets of Agilent expression data
##
##  Input Lib(s)  :
##  Members(s)    :
##
##  Output Lib(s) :
##  Member(s)     :
##  Modified by   : Gregory Way 
##  Modified date : 19 September 2015 12:00 PM EST
######################################################## 

######################## RUN under R version 2.15.0 (2012-03-30)
######################## TESTED under R version 3.1.2 (2014-10-31)
###################################################################

#Setting stringsAsFactors to FALSE
#Loading rlocal for Mayo HSR
#
#
#Attaching package: rlocal
#
#The following object(s) are masked frompackage:base:
#
#    table
#
#Loading required package: stats
#Loading required package: utils
#Loading required package: graphics
#Loading required package: splines
#Loaded the rlocal and survival packages
## get the older versions of all the packages in the gentools locations
## addLibGentools(first=TRUE)
### ComBat for batch adjustment
###--------------------------------------------------------------------------
#=== sva is the R package with ComBat
#require(sva)
## load library 'sva' that is needed for ComBat batch adjustment
library(sva)
library(Biobase)
#Loading required package: corpcor
#Loading required package: mgcv
#This is mgcv 1.7-16. For overview type 'help("mgcv-package")'.

## set the path where the data is located
combatoutdir<-'1.DataInclusion/Data/Mayo/'

# Input logic to import raw data from GEO HERE
#load(paste(combatoutdir, 'PreComBatAgilent3BatchesExpr.RData', sep = ''))
mayo.GEO.entry <- GEOquery::getGEO("GSE74357", getGPL = FALSE)
mayo.eset <- mayo.GEO.entry[[1]]
mayo.p <- phenoData(mayo.eset)
mayo.expression <- exprs(mayo.eset)

# Map GEO sample accession IDs to "X" format
mapper <- data.frame(geo.accession = mayo.p$geo_accession)
mapper["unique_patient_ID"] <-
  unlist(lapply(mayo.p$title,
                function(x) strsplit(toString(x), split = "  ")[[1]][2]))

set.order <- data.frame(geo.accession = colnames(mayo.expression))
mapper <- dplyr::inner_join(set.order, mapper, by = "geo.accession")
colnames(mayo.expression) <- unlist(mapper["unique_patient_ID"])

load(paste(combatoutdir, 'PreComBatAgilent3BatchesInfo.RData', sep = ''))

objects()
#[1] "combatoutdir"    "r_batch4"        "mayo.expression"

mch.batch.idx <- match(colnames(mayo.expression), r_batch4$UniqueID)
r_batch_info <- r_batch4[mch.batch.idx, ]

table(r_batch_info$AgilentPart)
#  1   2   3 
#257 216  56 

table(r_batch_info$cy5cy3dates)
# 2009-8or9#2/24/10   2009-8or9#4#6/10  2010-2or3#2/24/10   2010-2or3#4#6/10 
#                23                 11                 76                 55 
#  2010-4or5#4#6/10      2011-8#4#6/10            2013-11 2013-1or2# 1.18.13 
#                81                 11                 56                 90 
# 2013-1or2# 2.5.13 
#               126 


merg_batch_vec <- paste(r_batch_info$cy5cy3dates)
all(colnames(mayo.expression) == r_batch_info$UniqueID)
#[1] TRUE
sum(is.na(mayo.expression))
#[1] 98057

#=== 04/22/2013: Chen decided to exclude both probes with low SD (SD<=0.05) and probes missing in more than 10% of the samples (>=52 samples)
#=== Combat will fail if including probes with very small variance (equivalently very small SD)
sd_vec <- apply(mayo.expression, 1, sd,na.rm = TRUE)

sel.prb.idx <- which(!is.na(sd_vec) & sd_vec >= 0.05)
length(sel.prb.idx)
#[1] 37270

combat.set3 <- mayo.expression[sel.prb.idx, ]
dim(combat.set3)
#[1] 37270   529
sum(is.na(combat.set3))
#[1] 96299

na.vec3 <- apply(is.na(combat.set3), 1, sum)
summary(na.vec3)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   0.000   0.000   2.584   1.000 429.000 

miss.prb.idx <- na.vec3 < 52
table(miss.prb.idx)
#miss.prb.idx
#FALSE  TRUE 
#  395 36875 

combat.set4 <- combat.set3[miss.prb.idx, ]
dim(combat.set4)
#[1] 36875   529

sum(is.na(combat.set4))
#[1] 29815

combat.set5 <- t(combat.set4)
dim(combat.set5)
#[1]   529 36875

imputeMiss <- function(x) {
  mean.x<-colMeans(x,na.rm = T)
  x[is.na(x)]<-mean.x[col(x)[is.na(x)]]
  return(x)
}

### impute the missing values (samples on rows and probes on columns) with the probe mean
combat.set6 <- imputeMiss(combat.set5)
dim(combat.set6)
#[1]   529 36875

table(is.na(combat.set6))
#   FALSE 
#19506875 

combat.set7 <- t(combat.set6)

all(colnames(combat.set7) == r_batch4$UniqueID)
#[1] TRUE
merg_batch_vec <- paste(r_batch4$cy5cy3dates)

comb.mtx <- ComBat(batch = merg_batch_vec, dat = combat.set7, mod = NULL)

dim(comb.mtx)
#[1] 36875   529


########## put the missing values back after ComBat and save the data
dim(combat.set4)
#[1] 36875   529

sum(is.na(combat.set4))
#[1] 29815

table(rownames(combat.set4) == rownames(comb.mtx))
# TRUE 
#36875 

table(colnames(combat.set4) == colnames(comb.mtx))
#TRUE 
# 529 

sum(is.na(comb.mtx))
#[1] 0

comb.mtx2 <- comb.mtx
comb.mtx2[is.na(combat.set4)] <- NA

dim(comb.mtx2)
#[1] 36875   529
sum(is.na(comb.mtx2))
#[1] 29815

# Write comb.mtx99 to file
write.table(comb.mtx2, "1.DataInclusion/Data/Mayo/COMBATadj_withNAcy5cy3.tsv", sep = "\t")
