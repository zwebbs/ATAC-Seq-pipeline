# File Name: plot_ATAC_diagnostics.R
# Created On: 2022-08-05
# Created By: ZW
# Purpose: Generate plots for fragment length distrubition
# and TSS enrichment. Requires R 4.1.0

# Setup Command Line Arguments
# -----------------------------------------------------------------------------

library(optparse)
library(Rsamtools)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

# build flag options
option_list <- list(
    make_option(c("-b", "--bam"), type="character", default=NULL, 
                help="ATAC-Seq alignments BAM", metavar="character"),
    make_option(c("-n", "--sample-name"), type="character", default=NULL, 
                help="Analysis name for plotting", metavar="character"),
    make_option(c("-g", "--genome"), type="character", default="hg38",
                help="whether to run TSSe from hg19 or hg38 [default %default]"),
    make_option(c("-o","--outfile"), type="character",
                help="output file name for the plot (should be .pdf)")
)

# parse commandline args
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# throw error if passed genome is not implemented
# otherwise load the appropriate genome database
if (opt$genome == "hg38") {
    # implement hg38
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
    whichp <- as(seqinformation, "GRanges")
} else if (opt$genome == "hg19") {
    # implement hg19
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
    whichp <- as(seqinformation, "GRanges")
} else {
    # other genomes not implemented
    errmsg <- sprintf("Error: Passed Genome: (%s) not implemented.", opt$genome)
    stop(errmsg)
}


# Plot TSS enrichment and Fragment Size Distribution
# -----------------------------------------------------------------------------

bam.path <- opt$bam
bam.name <- opt$`sample-name`

bf.ga <- readBamFile(bam.path, which=whichp, asMates=FALSE, bigFile=TRUE)
tsse <- TSSEscore(bf.ga, txs)

# generate plots
pdf(opt$outfile, width=12, height=6)
par(mfrow=c(1,2))
plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score",
     main=sprintf("TSS enrichment %s",bam.name)
)
fragSizeDist(c(bam.path), bam.name)
dev.off()

