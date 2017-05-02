# Refer to:
# http://kasperdanielhansen.github.io/genbioconductor/html/Usecase_AnnotationHub_GRanges.html
# For extra bits
library(AnnotationHub)
library(GenomicRanges)

AH <- AnnotationHub()
AH

AH <- subset(AH, species == "Homo sapiens")
methQuery <- query(AH, c("H3K4me3", "Gm12878"))
methQuery

# First, methylation peak
gr1 = methQuery[[1]]
gr2 = methQuery[[4]]

summary(width(gr1))
summary(width(gr2))
table(width(gr2))

peaks <- gr2

# Second, promoter regions (from RefSeq)
(refQuery <- query(AH, "RefSeq"))
refGenes <- refQuery[[1]]
refGenes
# Number of isoforms per gene name
table(table(refGenes$name))
# Blocks here represent exons
refGenes$blocks

promoters <- promoters(refGenes)
# Simple as that...

# Now we ask, where and how much do promoters and methylation peaks overlap?
(OV <- findOverlaps(promoters, peaks))

length(unique(queryHits(OV)))
length(unique(subjectHits(OV)))

# Areas of overlaps
subsetByOverlaps(promoters, peaks, ignore.strand = TRUE)

# Percentage of overlaps
length(subsetByOverlaps(promoters, peaks, ignore.strand = TRUE))
length(subsetByOverlaps(promoters, peaks, ignore.strand = TRUE)) / length(peaks) * 100
# 33 percent overlaps of promoter over peaks
# i.e. 33 percent of tri-methylations are in promoters
length(subsetByOverlaps(peaks, promoters, ignore.strand = TRUE)) / length(promoters) * 100
# 44 percent overlaps of peaks over promoters
# i.e. 44 percent of tri-methylations are in promoters (where are the rest?)

# Lets see how many bases the peaks cover?
sum(width(reduce(peaks, ignore.strand = TRUE))) / 10^6
# 11 mega bases

sum(width(reduce(promoters, ignore.strand = TRUE))) / 10^6
# 62 mega bases

# How much do they intersect in our data? (note that there is no call to reduce)
sum(width(intersect(promoters, peaks, ignore.strand = TRUE))) / 10^6
# 3 mega bases

# Let's see if these are significant
# First create a 2-by-2 matrix
inOut <- matrix(0, nrow = 2, ncol = 2)
rownames(inOut) <- c("peak_in", "peak_out")
colnames(inOut) <- c("prom_in", "prom_out")
inOut

# Now, how many bases of promoters and peaks overlap?
(inOut[1,1] <- sum(width(intersect(promoters, peaks, ignore.strand = TRUE))))

# How many in promoters, but not in peaks?
(inOut[2,1] <- sum(width(setdiff(promoters, peaks, ignore.strand = TRUE))))

# How many in peaks, but not in promoters?
(inOut[1,2] <- sum(width(setdiff(peaks, promoters, ignore.strand = TRUE))))

rowSums(inOut) # 11 mega bases
colSums(inOut) # 62 mega bases

# And finally, how many in neither? (suppose the human genome is 3 billion base pairs)
(inOut[2,2] <- 3*10^9 - sum(inOut))





# ==== Quiz 1 ====

AH <- AnnotationHub()
AH <- subset(AH, species == "Homo sapiens")

cpg <- query(AH, "CpG islands")
cpgGenes <- cpg[[1]] 
(cpgGenes <- keepStandardChromosomes(cpgGenes))
seqinfo(cpgGenes)
(autosomalCpG <- dropSeqlevels(cpgGenes, c("chrX", "chrY", "chrM")))
seqnames(autosomalCpG)
length(autosomalCpG)
length(unique(autosomalCpG))

# How many CpG islands in the autosomes? 26,641
# How many on Chromosome 4? 1031

# Obtain data for H3K4me3
AH <- AnnotationHub()
AH <- subset(AH, species == "Homo sapiens")

(Q = query(AH, c("H3K4me3")))
hist <- display(Q)
H1methylation <- AH[[names(hist)]]

seqinfo(H1methylation)
H1methylation <- keepStandardChromosomes(H1methylation)
seqinfo(H1methylation)
autosomalH1 <- dropSeqlevels(H1methylation, c("chrX", "chrY", "chrM"))
seqinfo(autosomalH1)

# How many bases do these regions cover?
sum(width(reduce(autosomalH1, ignore.strand = TRUE)))

# Obtain data for H3K27me3
Q <- query(AH, "H3K27me3")
hist <- display(Q)
(K27meth <- AH[[names(hist)]])
seqinfo(K27meth)
K27meth <- keepStandardChromosomes(K27meth)
autosomalK27 <- dropSeqlevels(K27meth, c("chrX", "chrY", "chrM"))
seqinfo(autosomalK27)

# What is the mean signalValue accross all regions?
mean(autosomalK27$signalValue)

# How many bases are bivalently methylated?
sum(width(intersect(autosomalH1, autosomalK27, ignore.strand = TRUE)))

# What fraction of bivalent regions overlap CpG islands?
(bivalentRegions <- intersect(autosomalH1, autosomalK27, ignore.strand = TRUE))
OV <- findOverlaps(bivalentRegions, autosomalCpG, ignore.strand = TRUE)
length(unique(queryHits(OV))) / length(bivalentRegions)

# What fraction of _bases_ in CpG islands are bivalently marked?
sum(width(intersect(autosomalCpG, bivalentRegions))) / sum(width(autosomalCpG))

# How many bases are bivalently marked within 10kb of CpG Islands?
?resize

withinTen <- flank(autosomalCpG, width = 10000, both = TRUE)
autosomalCpG
withinTen

OV2 <- findOverlaps(withinTen, bivalentRegions, ignore.strand = TRUE)
sum(width(intersect(bivalentRegions, withinTen)))


# How big a fraction of the human genome is contained in a CpG Island?
sum(width(reduce(autosomalCpG))) / sum(as.numeric(seqlengths(autosomalCpG)))
