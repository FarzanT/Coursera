library(BiocInstaller)
biocLite(c("ALL", "Biobase", "hgu95av2.db"))
library(ALL)
?data
data(ALL)
ALL
?ALL
# Get the experminent related information, and the PMID of related papers!
experimentData(ALL)
abstract(ALL)
exprs(ALL)[1:4, 1:4]
head(sampleNames(ALL))
head(featureNames(ALL))
# pData for phenotype data
pData(ALL)[1,]
# Recall that rows are features, and columns are samples; can subset accordingly
featureData(ALL)
# Can contain information about the genes
ids = featureNames(ALL)[1:5]
ids
# Can use hgu95av2.db to map the Affymatrix IDs to gene IDs
library(hgu95av2.db)
hgu95av2ENTREZID[ids]
# Get the ENTREZ IDs of Affymetrix IDs
as.list(hgu95av2ENTREZID[ids])
# Note the difference between between pData and phenoData
phenoData(ALL)
pData(ALL)[1,]
pData(phenoData(ALL))[1,]
varLabels(ALL)


# ==== airway ====
library(GenomicRanges)
library(Biobase)
biocLite("airway")
library(airway)
data("airway")
airway
?airway
# Use colData in place of pData for SummarizedExperiment sets
colData(airway)
# colnames for sample names
colnames(airway)
# rownames for feature names
rownames(airway)

# To get data associated with the experiment, use the assay() function
airway
# We can see that there is only one assay: "counts"
assayNames(airway)
assay(airway, "counts")[1:4, 1:4]
# Note that each row has a GRanges associated with it!
# Use the rowRanges() function
length(rowRanges(airway))
# There are 64102 features
rowRanges(airway)
# Returns a GRangesList, where each gene has the exons associated with it in a
# GRanges object
# e.g. The first gene has 17 exons
# Use elementNROWS to see the number of exons per element (i.e. gene)
sum(elementNROWS(rowRanges(airway)))
# There are 745,593 exons associated with 64102

# start() returns the starting positions of the exons associated with each gene
start(airway)
# Can use subsetByOverlaps() to retrieve a GRange
gr = GRanges("1", ranges = IRanges(1, 10 ^ 7))
subsetByOverlaps(airway, gr)
# There are 329 genes in the given range on chromosome 1

# ==== GEOquery ====
biocLite("GEOquery")
library(GEOquery)
# Can retrieve a dataset associated with a study
# Must have a GEO identifier to query the database
eList = getGEO("GSE11675")
eList
# The data is returned in a list since some experiments might have multiple sets
# associated with them
length(eList)
names(eList)
# Subset the list 
eData = eList[[1]]
# This is a simple expression set
eData
# This is the data associated with the expression set
names(pData(eData))
# We can also get the RAW DATA for the same experiments by using
# getGEOSuppFiles()

eList2 <- getGEOSuppFiles("GSE11675")
length(eList2)
names(eList2)

# ==== biomaRt ====
biocLite("biomaRt")
library(biomaRt)

# Use listMarts() to list providers
listMarts()
# Determine the database to use by using useMart()
mart <- useMart("ensembl")
# Get a list of available datasets by using listDatasets()
head(listDatasets(mart))
# Find dataset of choice
listDatasets(mart)
# Retrieve the dataset of choice
ensembl_hsapiens <- useDataset("hsapiens_gene_ensembl", mart = mart)
ensembl_hsapiens

# Retrieve gene_ids based on Affymatrix IDs
values = c("202763_at")
# In order to retrieve the one-to-one mapping, must provide both the ensembl as
# well as the affy_hg 
getBM(attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"), 
      filters = "affy_hg_u133_plus_2", values = values, mart = ensembl_hsapiens)

# To get a list of available attributes:
att <- listAttributes(ensembl_hsapiens)
# There are 1477 attributes
nrow(att)
head(att)
# Some of these attributes allow matching a gene against its orthologue, e.g.
# olive baboon genes

# Can also listFilters()
filters <- listFilters(ensembl_hsapiens)
nrow(filters)
head(filters)

# Pages are internal to the database structure
attributePages(ensembl_hsapiens)

# There maybe duplicate rows in the queries, beware


# ==== Quiz 3 ====
# === Q1
# What is the mean expression across all features for sample 5 in the ALL 
# dataset (from the ALL package)?
library(ALL)
mean(exprs(ALL)[,5])

# === Q2
mart <- useMart(host = "feb2014.archive.ensembl.org",
                biomart = "ENSEMBL_MART_ENSEMBL")
mart
# Using this version of Ensembl, annotate each feature of the ALL dataset with
# the Ensembl gene id. How many probesets (features) are annotated with more
# than one Ensembl gene id?
ALL
head(featureNames(ALL))
allData <- listDatasets(mart)
grep("Homo", allData$description)
# 31
allData[31,]
# hsapiens_gene_ensembl
HsapGenes <- useDataset("hsapiens_gene_ensembl", mart)

att <- listAttributes(HsapGenes)
# Get annotation type
annotation(ALL)
# Find Affymetrix related annotations
listFilters(mart = HsapGenes)[120:190,]
listFilters(mart = HsapGenes)[164,]
vals <- featureNames(ALL)


annot <- getBM(attributes = c("ensembl_gene_id", "affy_hg_u95av2"),
      filters = c("affy_hg_u95av2", "chromosome_name"), values = vals, HsapGenes)
nrow(annot)
any(duplicated(annot$affy_hg_u95av2))
sum(duplicated(annot$affy_hg_u95av2)) / 2

# === Q3
# How many probesets (Affymetrix IDs) are annotated with one or more genes on
# the autosomes (chromosomes 1 to 22)?
filt <- listFilters(HsapGenes)
head(filt)
vals <- list(affyid = featureNames(ALL), chromosome = as.character(seq(1, 22)))
annot <- getBM(attributes = c("ensembl_gene_id", "affy_hg_u95av2"),
               filters = c("affy_hg_u95av2", "chromosome_name"), values = vals, HsapGenes)
nrow(annot)
length(unique(annot$affy_hg_u95av2))

# === Q4
# Use the MsetEx dataset from the minfiData package. Part of this question is
# to use the help system to figure out how to address the question.

# Question: What is the mean value of the Methylation channel across the
# features for sample “5723646052_R04C01”?
biocLite("minfiData")
library(minfiData)
minfiData::MsetEx
?MsetEx
names(pData(MsetEx))
varLabels(MsetEx)
pData(MsetEx)
featureNames(MsetEx)
# Get sample methylation mean
mean(getMeth(MsetEx[,"5723646052_R04C01"]), na.rm = TRUE)

# === Q5
# Access the processed data from NCBI GEO Accession number GSE788. What is the
# mean expression level of sample GSM9024?
geo <- GEOquery::getGEO("GSE788")
geo <- geo[[1]]
mean(exprs(geo[,"GSM9024"]))

# === Q6
# We are using the airway dataset from the airway package.

# Question: What is the average of the average length across the samples in the
# expriment?
# Recall that airway is of type RangedSummarizedExperiment
### Nope
# class(airway)
# # Must use colNames and rowNames instead
# colnames(airway)
# nrow(airway)
# ncol(airway)
# # Average row wise then column wise
# colMeans(assay(airway, "counts"))
# mean(colMeans(assay(airway, "counts")))
###
# What the question was asking for
mean(colData(airway)$avgLength)

# === Q7
# We are using the airway dataset from the airway package. The features in this
# dataset are Ensembl genes.
# 
# Question: What is the number of Ensembl genes which have a count of 1 read or
# more in sample SRR1039512?
sum(assay(airway, "counts")[,"SRR1039512"] > 0)

# === Q8
# Question: The airway dataset contains more than 64k features. How many of
# these features overlaps with transcripts on the autosomes (chromosomes 1-22)
# as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?
# 
# Clarification: A feature has to overlap the actual transcript, not the intron
# of a transcript. So you will need to make sure that the transcript
# representation does not contain introns.
library(airway)
data("airway")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# Remember to restore seq levels
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
restoreSeqlevels(txdb)
seqlevels(txdb)

txdb <- keepStandardChromosomes(txdb)
txdb <- dropSeqlevels(txdb, c("chrY", "chrX", "chrM"))
seqlevels(txdb)

# Retain only the exons
exons <- exonsBy(txdb, "tx")
# Retrieve airway exon range
airwayExons <- rowRanges(airway)
airwayExons <- keepStandardChromosomes(airwayExons)
airwayExons <- dropSeqlevels(airwayExons, c("X", "Y", "MT"))
seqlevels(airwayExons)
# Aha! Must rename exons 
airwayExons <- renameSeqlevels(airwayExons, paste0("chr", seq(1, 22)))
# Find overlaps
?findOverlaps
overlaps <- subsetByOverlaps(airwayExons, exons)
length(overlaps)

# === Q9
# For sample SRR1039508, how big a percentage (expressed as a number between 0
# and 1) of the total reads in the airway dataset for that sample, are part of
# a feature which overlaps an autosomal TxDb.Hsapiens.UCSC.hg19.knownGene
# transcript?
seqlevels(txdb)
trs <- transcripts(txdb)
airwayExons <- rowRanges(airway[, "SRR1039508"])
airwayExons <- keepStandardChromosomes(airwayExons)
airwayExons <- dropSeqlevels(airwayExons, c("X", "Y", "MT"))
seqlevels(airwayExons)
# Aha! Must rename exons 
airwayExons <- renameSeqlevels(airwayExons, paste0("chr", seq(1, 22)))

overlaps <- subsetByOverlaps(trs, airwayExons)
sum(overlaps) / sum(elementNROWS(rowRanges(airway[, "SRR1039508"])))


# === Q10
# Obtain the H3K4me3 narrowPeaks from the E096 sample using the AnnotationHub
# package.
# What is the median number of counts per feature (for sample SRR1039508)
# containing a H3K4me narrowPeak in their promoter (only features which overlap
# autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene are considered)?
AH <- AnnotationHub()
q <- query(AH, c("E096", "narrow", "H3K4me3"))
narrow <- q[[1]]
# Retrieve features that contain H3K4me3 narrowPeak in their promoter
narrow <- keepSeqlevels(narrow, paste0("chr", seq(1, 22)))

# Overlap with those from sample SRR1039508
exons <- rowRanges(airway[, "SRR1039508"])
exons <- keepSeqlevels(exons, as.character(seq(1, 22)))
exons <- renameSeqlevels(exons, paste0("chr", seq(1, 22)))
exons <- sapply(exons, '[', 1)
exons <- GRangesList(exons)
proms <- promoters(exons)

overlaps <- subsetByOverlaps(proms, narrow)

# Get overlap ENG IDs
allNames <- names(overlaps)
# Now retieve counts using names
counts <- assay(airway[, "SRR1039508"], "counts")
median(counts[allNames,])
mean(counts[allNames,])

