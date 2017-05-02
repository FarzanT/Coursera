library(Biostrings)

dna = Biostrings::DNAString("TCGAAA")
dna
reverse(dna)
rev(dna)
reverseComplement(dna)
complement(dna)
translate(dna)
alphabetFrequency(dna)
letterFrequency(dna, "GC")
dinucleotideFrequency(dna)
# Good for PWMs
consensusMatrix(dna)

# BSGenome
library(BSgenome)
?BSgenome
available.genomes()
library(BiocInstaller)
biocLite("BSgenome.Scerevisiae.UCSC.sacCer2")
library("BSgenome.Scerevisiae.UCSC.sacCer2")
Scerevisiae


seqnames(Scerevisiae)
seqlengths(Scerevisiae)

Scerevisiae$chrI

letterFrequency(Scerevisiae$chrI, "GC")
letterFrequency(Scerevisiae$chrI, "GC", as.prob = TRUE)

# Better than apply for genomes, loads and unloads genomes as needed
?bsapply
param = new("BSParams", X = Scerevisiae, FUN = letterFrequency)
# Put additional arguments in the bsapply function call
bsapply(param, "GC")
bsapply(param, "GC", as.prob = TRUE)
unlist(bsapply(param, "GC", as.prob = TRUE))

# Overall GC content
sum(unlist(bsapply(param, "GC"))) / sum(seqlengths(Scerevisiae))

# Biostring - Matching
dnaseq <- DNAString("ACGTACGT")
matchPattern(dnaseq, Scerevisiae$chrI)
# Count the number of matches
countPattern(dnaseq, Scerevisiae$chrI)

# Count against a set of sequences
vmatchPattern(dnaseq, Scerevisiae)
# Hits both on the reverse strands

dnaseq == reverseComplement(dnaseq)

# Position weight matrix matching
matchPWM()

# Global or local
pairwiseAlignment()

# Trim on the flanking regions?
trimLRPatterns()

# View
vi = matchPattern(dnaseq, Scerevisiae$chrI)
vi
ranges(vi)
Scerevisiae$chrI[57932:57939] == dnaseq

gr = vmatchPattern(dnaseq, Scerevisiae)
gr

# Can be useful when the matching sequence must be in the same object, e.g. PWM
vi2 = Views(Scerevisiae, gr)
vi2

library(AnnotationHub)
AH <- AnnotationHub()
Q <- query(AH, c("sacCer2", "genes"))
Q

genes = Q[[1]]

prom = promoters(genes)
prom

# Trim anything outside of sequence lengths of the genome
prom <- trim(prom)

promView = Views(Scerevisiae, prom)
promView

# Get the GC content
gcProm = letterFrequency(promView, "GC", as.prob = TRUE)
gcProm

plot(density(gcProm))
abline(v = 0.38)

# Is there a difference in the GC content of the promoters and non-protomer regions?

withoutProm = overlapsAny(Scerevisiae, prom)


# ==== RLE ====
library(GenomicRanges)
rl = Rle(c(1,1,1,1,1,1,2,2,2,2,2,5,5,5,5,5,6,6,2,2,1,1))
rl
runLength(rl)
runValue(rl)

# Display the original RLE by using as.numeric() 
as.numeric(rl)

ir = IRanges(start = c(2,8), width = 5)
ir

# Find the mean of values at the specified range ir in rl
aggregate(rl, ir, FUN = mean)

ir = IRanges(start = 1:5, width = 3)
ir
# e.g. After RNA seq, how many reads cover a certain base?
coverage(ir)

# Areas where this vector is large
as.numeric(rl)
# Below contains all the positions where the vector is geq to two
slice(rl, 2)

rl = Rle(c('A', 'T', 'C', 'G', 'G', 'A', 'A'))
rl

# Doesn't mean much
slice(rl, 2)

aggregate(rl, ir, print)


# ==== Genomic Features ====
library(BiocInstaller)
biocLite("GenomicFeatures")
library(GenomicFeatures)
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
# Create shortcut
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# txdb contains information regarding genes, trascripts and exons
txdb
# Transcript can both refer to RNA before or after splicing
gr <- GRanges(seqnames = "chr1", strand = "+", ranges = IRanges(start = 11874, end = 14409))
# Note that the gene_id uses ENTREZ Gene ID
GenomicFeatures::genes(txdb)

subsetByOverlaps(genes(txdb), gr, ignore.strand = TRUE)

transcripts(txdb)
# See overlaps with transcripts in the gr range
subsetByOverlaps(transcripts(txdb), gr)
# Each transcript has a separate name
# The transcripts are different (although they have the same start and end
# positions) because their exons are differnt
subsetByOverlaps(exons(txdb), gr)
# How are exons combined together to create transcripts?
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)
# Exons 5 and 6 overlap
# The $ sign associates with the tx_id

# Coding sequences can be hard to figure out
# Not all transcripts have coding sequences
# A transcript might have multiple open reading frames
subsetByOverlaps(cds(txdb), gr)
# These are not open reading frames...
subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)
# There are 3 transcripts in the gr sequence, but only transcript number 2 ($2)
# has coding potential
# Newer versions of UCSC all have cdcBy

subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)["2"]
# Note the 3' and 5' untranslated regions of the transcript


subset(transcriptLengths(txdb, with.cds_len = TRUE), gene_id == "100287102")
# 3 different transcripts, each have 3 exons with different lengths
# cds_len of 0 means there is no coding sequence
# Confirm the coding sequence length
sum(width(subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)[["2"]]))

# Users are advised to create their own TXDB objects, from their preferred
# sources
makeTxDbFromBiomart()
makeTxDbFromUCSC()


# ==== rtracklayer DataImport ====
library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)
# Check the supported formats

# A Wiggle file is a file that holds a signal along the genome
# A BigWig file contains a single vector along the length of the genome; it is
# compressed and is used to extract the values for a given region
AH <- AnnotationHub()
table(AH$rdataclass)
# There are 10247 BigWigFiles

AH.bw = subset(AH, rdataclass == "BigWigFile" & species == "Homo sapiens")
AH.bw
bw = AH.bw[[1]]
bw
# Since bw is a very large vector, it's best to only read the required part
?import
gr.chr22 = import(bw, which = GRanges("chr22", ranges = IRanges(1, 10^8)))
gr.chr22
# Note that the file is a GRanges object


# Lift over is a tool from UCSC that allows conversion between different genome
# versions
AH.chain = subset(AH, rdataclass == "ChainFile")
AH.chain
# Subset further to get different genomes (e.g. Homo Sapiens)
AH.chain = subset(AH.chain, species == "Homo sapiens")
AH.chain

# Now to convert between different genome versions
query(AH.chain, c("hg18", "hg19"))
# Take the first one
chain = query(AH.chain, c("hg18", "hg19"))[[1]]
# Simply convert between the two
gr.hg18 = liftOver(gr.chr22, chain)
# There is no check that the genomes are correct
# The returned object is a GRangesList, and the length of the list is equal to 
# the length of the object initially passed (gr.chr22 in this case)

elementLengths(gr.hg18)
table(elementLengths(gr.hg18))
# We can see how sequences were lifted over, how many could not and the 
# intervals that got broken up into pieces


# ==== QUIZ 2 ====
# ==== Q1
library(BSgenome)
library(AnnotationHub)
library(GenomicRanges)
AH <- AnnotationHub()
# table(AH$rdataclass)
# AHsub <- subset(AH, species == "Homo sapiens" & rdataclass == "BigWigFile")
# AHsub
# query(AHsub, c("hg19"))
# BW = query(AHsub, "hg19")[[1]]
# gr.chr22 = rtracklayer::import(BW, which = GRanges(seqnames = "chr22", ranges = IRanges(1, 10^8)))
# gr.chr22
# class(gr.chr22)
# 
# # Now use bsapply with letterFrequency
# param = new("BSParams", X = BW, FUN = letterFrequency)

library(BiocInstaller)
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
Hsapiens
seqinfo(Hsapiens)
seqlengths(Hsapiens)
# Now for chromosome 22
Hsapiens$chr22
letterFrequency(Hsapiens$chr22, "GC") / letterFrequency(Hsapiens$chr22, "ACGT")

# ==== Q2
library(BSgenome)
library(IRanges)
library(GenomicFeatures)
??GenomicFeatures
# What is mean GC content of H3K27me3 “narrowPeak” regions from Epigenomics
# Roadmap from the H1 stem cell line on chr 22?
AH <- AnnotationHub()
AH <- subset(AH, species == "Homo sapiens")
query(AH, c("H3K27me3", "H1", "narrow"))
# AH29892
narrow <- AH[["AH29892"]]
narrow
class(narrow)
seqinfo(narrow)
seqnames(narrow)
seqlengths(narrow)
# This is a nice trick for choosing a certain sequence (e.g. chromosome)
seqlevels(narrow, force = TRUE) <- "chr22"
seqlevels(narrow)
narrow
# Now import the same regions (narrow peaks) from hg19, and compute GC content
available.genomes()
class(Hsapiens)
# Use Views to overlap ranges with actual sequences
chr22.peaks <- Views(Hsapiens, narrow)
# Compute the GC content for each peak region as a percentage and then average
# those percentages to compute a number between 0 and 1.
peakGC <- mean(letterFrequency(chr22.peaks, "GC", as.prob = TRUE))

# === Q3
# What is the correlation between GC content and “signalValue” of these regions
# (on chr22)?]
narrow
cor(letterFrequency(chr22.peaks, "GC", as.prob = TRUE), narrow$signalValue)


# === Q4
library(GenomicRanges)
# what is the correlation between the “signalValue” of the “narrowPeak” regions
# and the average “fc.signal” across the same regions?
# Obetain the fc.signal data for the same cell line and histone modification
AH <- AnnotationHub()
q <- query(AH, c("fc.signal", "hg19", "Homo sapiens", "H3K27me3"))
q
# E003 is AH32033
fc.sig <- AH[["AH32033"]]
class(fc.sig)
seqinfo(fc.sig)
# Only take in chromosome 22
fc <- import(fc.sig, which = granges(narrow), as = "RleList")
fc.chr22 <- fc$chr22
class(fc.chr22)
# First compute the average “fc.signal” for across each region, for example
# using “Views”
# narrow is a GRanges object, so is fc.chr22
# Use aggregate.......
fcMeans <- aggregate(fc.chr22, narrow, mean)
cor(fcMeans, narrow$signalValue)

# === Q5

fc.chr22 <- import(fc.sig, which = GRanges("chr22", ranges = IRanges(1, 3 * 10 ^ 8)))
sum(width(fc.chr22[score(fc.chr22) >= 1]))

# === Q6
# Identify the regions of the genome where the signal in E003 is 0.5 or lower
# and the signal in E055 is 2 or higher.
AH <- AnnotationHub()
q <- query(AH, c("E055", "H3K27me3", "fc.signal"))
q
E55 <- q[[1]]
seqinfo(E55)
E55.chr22 <- import(E55, which = GRanges("chr22", ranges = IRanges(1, 3 * 10 ^ 8)))
E55.chr22
# Recall that chr22 fc.signals for E003 is in fc.chr22
E3V <- Views(Hsapiens, fc.chr22[(score(fc.chr22) <= 0.5)])
E55V <- Views(Hsapiens, chr22[(score(E55.chr22) >= 2)])
E3Vgranges <- granges(E3V)
E55Vgranges <- granges(E55V)
intersection <- intersect(E3Vgranges, E55Vgranges)
sum(width(intersection))

# === Q7
# What is the average observed-to-expected ratio of CpG dinucleotides for CpG
# Islands on chromosome 22?
AH <- AnnotationHub()
q <- query(AH, c("CpG Islands"))
q
cpg <- q[[1]]
cpg22 <- keepSeqlevels(cpg, "chr22")
cpgRegions <- IRanges::Views(Hsapiens, cpg22)
BSgenome::alphabetFrequency()
cpgRegions
cpgRegions
observed <- sum(dinucleotideFrequency(cpgRegions)[,"CG"])
temp <- sum(letterFrequency(cpgRegions, "G")) / sum(width(cpgRegions))
expected <- temp * sum(letterFrequency(cpgRegions, "C"))
observed / expected

# === Q8
# How many TATA boxes are there on chr 22 of build hg19 of the human genome?
# Both forward and reverse strands
seqinfo(Hsapiens)
tata <- PDict("TATAAA")
as.character(reverseComplement(DNAString("TATAAA")))
atat <- PDict("TTTATA")
countPDict(tata, Hsapiens$chr22) + countPDict(atat, Hsapiens$chr22)

# === Q9
# How many promoters of transcripts on chromosome 22 containing a coding
# sequence, contains a TATA box on the same strand as the transcript?
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Keep only chr22
seqlevels(txdb) <- "chr22"
proms <- promoters(txdb, upstream = 900, downstream = 100)
promsPos <- proms[strand(proms) == "+"]
promsNeg <- proms[strand(proms) == "-"]
sum(vcountPattern("TATAAA", DNAStringSet(posSeq)))
sum(vcountPattern("TATAAA", DNAStringSet(negSeq)))
tata <- PDict("TATAAA", )
posSeq <- Views(Hsapiens, promsPos)
negSeq <- Views(Hsapiens, promsNeg)
tataCountPos <- sum(vcountPDict(tata, DNAStringSet(posSeq)))
tataCountNeg <- sum(vcountPDict(tata, DNAStringSet(negSeq)))
tataCountNeg + tataCountPos

# Get the transcripts that have _coding sequences_ in them
coding <- cdsBy(txdb, by = "tx")
overlaps <- subsetByOverlaps(coding, promoters(txdb, upstream = 900, downstream = 100))
overlapsSeq <- Views(Hsapiens, unlist(overlaps))
pos <- overlapsSeq[strand(overlapsSeq) == "+"]
neg <- overlapsSeq[strand(overlapsSeq) == "-"]

sum(width(overlapsSeq))
sum(vcountPattern("TATAAA", DNAStringSet(pos))) + sum(vcountPattern("TATAAA", DNAStringSet(neg)))
sum(vcountPattern("TATAAA", DNAStringSet(overlapsSeq)))
reverse <- reverseComplement(DNAStringSet(overlapsSeq))
# Not sure
sum(vcountPattern("TATAAA", reverse)) + sum(vcountPattern("TATAAA", DNAStringSet(overlapsSeq)))

# === Q10
# How many bases on chr22 are part of more than one promoter of a coding 
# sequence?
# coding <- cdsBy(txdb, by = "tx")
# prom <- promoters(coding, upstream = 900, downstream = 100)
# cov <- coverage(prom)
# sum(width(slice(cov, 2)))

# Must subset for transcripts that have non-zero lengths
df <- subset(transcriptLengths(txdb, with.cds_len = TRUE))
df <- df[df$cds_len > 0,]
# Take the transcript names (or IDs)
names <- df$tx_name
# Subset all transcripts by the coding transcripts that have non-zero lengths,
# using the IDs obtained above
codingTrans <- transcripts(txdb)[transcripts(txdb)$tx_name %in% names]
# Find the promoters of these transcripts
proms <- promoters(codingTrans, upstream = 900, downstream = 100)
# Coverage returns an RLE for number of times a base is repeated (or an index
# in a range)
# Sum the number of the bases that are in at least 2 promoters
sum(width(slice(coverage(proms), 2)))
