# ==== Short Reads ====
biocLite("ShortRead")
library(ShortRead)
# Read fastq
reads <- readFastq("PATH")
reads
# Fastq has 4 lines per read
# Retrieve different information
sread(reads)
ShortRead::FastqQuality(reads)
ShortRead::id(reads)
# To get fret scale qualities, must convert the characters to matix
as.matrix(FastqQuality(reads))
# Read chunks of fastq files


# ==== RSAM Tools ====
biocLite("Rsamtools")
library(Rsamtools)
bamPath <- system.file("extdata", "ex1.bam", package = "Rsamtools")
bamFile <- BamFile(bamPath)
bamFile

# BAM supports seqlevels and seqinfo
seqinfo(bamFile)
seqlevels(bamFile)

# Read the bam file at once
aln <- scanBam(bamFile)
aln[[1]]

# Reads can align to the genome, or RNA (depending on the data)

# Read BAM in small chunks
# Determine chunk size
yieldSize(bamFile) <- 1
# Must open the bam file for this to work
open(bamFile)
# Now each call to scanBam will read _yieldSize_ elements
scanBam(bamFile)[[1]]$seq

# Close bam file
close(bamFile)
yieldSize(bamFile) <- NA

gr <- GRanges("seq2", ranges = IRanges(start = c(100, 1000), end = c(1500, 2000)))
gr
# Can run a query on the bam file, similar to BSGenome package
params <- ScanBamParam(which = gr, what = scanBamWhat())
# This determines the pieces of the BAM file that will be read; can be modified
scanBamWhat()

aln <- scanBam(bamFile, param = params)
aln
names(aln)
# All are in the specified gr range
head(aln[[1]]$pos)

# One sample, no range
bamView <- BamViews(bamPath)
aln <- scanBam(bamView)
aln
names(aln[[1]])
names(aln[[1]][[1]])

bamRanges(bamView) <- gr
aln <- scanBam(bamView)
aln
names(aln[[1]])

# Quick summary of the bamfile
quickBamFlagSummary(bamFile)


# ==== OLIGO ====
biocLite("oligo")
library(oligo)
library(GEOquery)
getGEOSuppFiles("GSE38792")
# Supplementary Affymetrix files are in a special binary format called
# CEL
list.files("GSE38792")
# Unzip and export to a directory called CEL
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL")
celFiles <- list.files("GSE38792/CEL", full.names = TRUE)
celFiles
rawData <- oligo::read.celfiles(celFiles)

rawData
# GeneFeatureSet Similar to a expression set containter
getClass("GeneFeatureSet")

# These are raw reads
exprs(rawData)[1:4, 1:3]
# Affymetrix 16 bit scanner, 2^16 - 1 levels (65536 - 1)
max(exprs(rawData))

# Clean up the names
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
?ifelse
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
pData(rawData)
class(rawData)
boxplot(rawData)
# Each box is a different sample
# The y axis is on a log scale
# There are different mean, spread
# Controls 5 to 7 have very low measurements

# Normalize the data using rma() function
# Backgroun correction
# Quantile normalization
# Expression calculation
normData <- rma(rawData)

# From 1 million features to 33 thousand
normData
featureNames(normData)[1:10]
# These are Affymetrix IDs, can use biomaRt to translate to gene IDs

boxplot(normData)
# The distributions are now a lot more similar

exprs(normData)[1:3, 1:3]


# ==== limma ====
# Linear models for microarrays
biocLite("limma")
library(limma)
biocLite("leukemiasEset")
library(leukemiasEset)
data("leukemiasEset")
leukemiasEset
table(leukemiasEset$LeukemiaType)
# Subset by ALL and NoL
myData <- leukemiasEset[,leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]
# Note that LeukemiaType is a factor
myData$LeukemiaType
# The levels still remain, must get rid of unneeded levels
?factor
myData$LeukemiaType = factor(myData$LeukemiaType)

design <- model.matrix(~ myData$LeukemiaType)
head(design)
# Intercept: Average gene expression across ALL samples
# Right column: Difference of gene expression between NoL and ALL
# Usual setup for two comparison

fit <- lmFit(myData, design)
# 
fit
fit <- eBayes(fit)
fit
# These are the list of differentially expression genes
topTable(fit)
topTable(fit, number = 1)
# This is similar to a feature name
geneName <- rownames(topTable(fit, number = 1))

typeMean <- tapply(exprs(myData)[geneName, ], myData$LeukemiaType, mean)
typeMean
typeMean["NoL"] - typeMean["ALL"]
# This means that the gene has lower expression in ALL
?"~"

# Get a new matrix where the column names are different, and two paramters are
# expression level ALL, and expression NoL (before, the second column was the 
# _difference_ in expression)
design2 <- model.matrix(~ myData$LeukemiaType - 1)
head(design2)
fit2 <- lmFit(design2)
# A hypothesis test: whether ALL - NoL is zero
?makeContrasts
contrast.matrix <- makeContrasts("ALL-NoL", levels = design2)
fit2C <- constrast.fit(fit2, constrast.matrix)
fit2C <- eBayes(fit2C)

topTable(fit2C)



# ==== minfi ====
# Methylation happens in CpG islands
biocLite("minfi")
library(minfi)
library(GEOquery)

getGEOSuppFiles("GSE68777")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
head(list.files("GSE68777/idat/", pattern = "idat"))

# Unzip IDAT files for reading
idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

# Can read all idat files in the directory using read.450k.exp
rgSet <- read.metharray.exp("GSE68777/idat/")
rgSet

# No information...
pData(rgSet)
head(sampleNames(rgSet))
# We have to load the original dataset, as that dataset has the original dataset
# associated with it
geoMat <- getGEO("GSE68777")
pD.all <- pData(geoMat[[1]])
colnames(pD.all)
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)

names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)

# common sample identifier and we make sure we re-order the phenotype data in
# the same order as the methylation data.
sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))
rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]
pData(rgSet) <- pD
rgSet

# Normalizes the data and maps it to the genome
grSet <- preprocessQuantile(rgSet)
# We get a GenomicRatioSet
grSet
granges(grSet)
# Useful for knowing distance from the CpG islands
head(getIslandStatus(grSet))

getBeta(grSet)[1:3, 1:3]

# Clusters of CpG's that change in the same direction


# ==== Count-Based RNA Seq Analysis
biocLite("DESeq2")
biocLite("edgeR")
library(DESeq2)
library(edgeR)

library(airway)
data("airway")
airway

granges(airway)

?DGEList
dge <- DGEList(counts = assay(airway, "counts"), group = airway$dex)
dge
dge$samples <- merge(dge$samples, as.data.frame(colData(airway)), by = 0)
dge$samples

head(dge$genes)

dge$genes <- data.frame(name = names(rowRanges(airway)), stringsAsFactors = FALSE)
head(dge$genes)

dge <- calcNormFactors(dge)
