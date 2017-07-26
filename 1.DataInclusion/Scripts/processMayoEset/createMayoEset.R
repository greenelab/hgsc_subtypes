############################################
# Cross-population analysis of high-grade serous ovarian cancer reveals only two robust subtypes
# 
# Way, G.P., Rudd, J., Wang, C., Hamidi, H., Fridley, L.B,  
# Konecny, G., Goode, E., Greene, C.S., Doherty, J.A.
# ~~~~~~~~~~~~~~~~~~~~~
# This script will take the normalized matrix from the Mayo Data and output 
# an eset to be used in subsequent analyses

##########################
# Load Libraries
##########################
library(Biobase)
library(biomaRt)

##########################
# CONSTANTS
##########################
# Input Files
ExpressionLocation <- "1.DataInclusion/Data/Mayo/COMBATadj_withNAcy5cy3.tsv"
PhenotypeLocation <- "1.DataInclusion/Data/Mayo/Mayo_Pheno_Data.csv"
ProbeToGeneMapLocation <- "1.DataInclusion/Data/Mayo/efg_agilent_wholegenome_4x44k_v1.csv"

# Output Files
GeneLevelExpressionOutputLocation <- "1.DataInclusion/Data/Mayo/MayoDataWithMap_Normalizer.txt"
ExpressionSetOutputLocation <- "1.DataInclusion/Data/Mayo/MayoEset.Rda"

#########################
# PERFORM ANALYSIS
#########################
#Load the expression data
comb.mtx2 <- read.delim(ExpressionLocation)

#Load the phenotype data
mayo.pheno <- read.table(PhenotypeLocation, row.names = 1, header = T, sep = ",")

# The following are taken from the data dictionary for the mayo pheno variables
# casecon: 1 Ovarian Case; 2 Control; 3 Fallopian Tube Case; 4 Primary Peritoneal Case)
# histology: 1 Serous; 2 Mucinous; 3 Endometroid; 4 Clear Cell; 5 Mixed Cell; 6 Carcinosarcoma; 7 Other)

# Restrict expression data to only those samples for which are high grade serousui
mayo.exprs <- comb.mtx2[ ,rownames(mayo.pheno)]

#load the GPL to map probes to genes
datGPL <- read.csv(file = ProbeToGeneMapLocation)

#Reorder the columns in datGPL so that we have gene name then probe ID
map <- datGPL[ ,c(2,1)]

# match the rownames of expression data with the IDs (2nd column) of the map:
map <- map[map[ ,2] %in% rownames(mayo.exprs), ]
map <- map[match(rownames(mayo.exprs), map[ ,2]), ]
map[ ,2] <- rownames(mayo.exprs)

rnames <- rownames(mayo.exprs)
mayo.exprs <- apply(mayo.exprs, 2, function(x) as.numeric(as.character(x)))
rownames(mayo.exprs) <- rnames

if(identical(all.equal(map[,2], rownames(mayo.exprs)), TRUE)) { 

    # Map to gene IDs.  Requires the Normalizer function from the Sleipnir library.
  expr.withmap <- cbind(map, mayo.exprs)
  expr.withmap <- na.omit(expr.withmap)
  rownames(expr.withmap) <- expr.withmap$probeset
  mayo.exprs <- na.omit(mayo.exprs)
  expr.withmap["probeset"] <- NULL
  
  write.table(expr.withmap, row.names = FALSE, file = "1.DataInclusion/Data/Mayo/MayoDataWithMap.txt", sep = "\t")
  expr.withmap <- WGCNA::collapseRows(mayo.exprs,
                                      rowGroup = expr.withmap$hgnc,
                                      rowID = rownames(expr.withmap))
} else {
  stop("Could not match map IDs to expression IDs.")
}

# Convert the phenotype data.frame into an annotated data.frame
mayo.pheno.ADF <- new("AnnotatedDataFrame", data = mayo.pheno)

# Create the expressionset object
mayo.eset <- new("ExpressionSet", exprs = expr.withmap$datETcollapsed, phenoData = mayo.pheno.ADF)

# Save the expressionset to the harddrive
save(mayo.eset, file=ExpressionSetOutputLocation)

