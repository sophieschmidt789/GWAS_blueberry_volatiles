
rm(list = ls())

library(QTLseqr)
library(data.table)

setwd("C:/Users/lidys/UFL Dropbox/Lidysce Mata Cantarero/Lidysce - Heat Tolerance/Lidysce/Classes/HOS6236 - Molecular Marker Assisted Plant Breeding/GWAS Project/filtered")

# Load data
G      <- fread("geno_filtered.tsv",   sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
PH     <- fread("phenotype_filtered.tsv", sep="\t", header=TRUE, data.table=FALSE)
snpmap <- fread("snp_map_filtered.tsv",sep="\t", header=TRUE, data.table=FALSE)
covar  <- fread("covariates_filtered.tsv", sep="\t", header=TRUE, data.table=FALSE)

#######################################################################################
############################# TRAIT #1 ################################################
#######################################################################################

# Define bulks (LOW/HIGH) from phenotype 
trait     <- "V1"     # change to V2/V3/V4 if you want a different trait
tail_prop <- 0.15     # bottom/top 15%

# Ensure sample order matches between G and PH (Taxa column)
stopifnot(identical(G$Taxa, PH$Taxa))

# Order by trait and pick tails
PH_ord     <- PH[order(PH[[trait]]), ]
n_per_bulk <- floor(nrow(PH_ord) * tail_prop)

LOW_ids  <- PH_ord$Taxa[1:n_per_bulk]
HIGH_ids <- PH_ord$Taxa[(nrow(PH_ord) - n_per_bulk + 1):nrow(PH_ord)]

# Sizes used later by runQTLseqAnalysis()
bulkSize <- c(length(LOW_ids), length(HIGH_ids))

# Align SNPs to the map and subset
# Figure out which SNPs are shared and put them in MAP order
snps_in_G   <- setdiff(names(G), "Taxa")
keep_snps   <- intersect(snpmap$SNP, snps_in_G)
snp_ordered <- snpmap$SNP[snpmap$SNP %in% keep_snps]  # preserves map order

# Subset G to those SNP columns (still keeping 'Taxa' first)
G_sub <- cbind(Taxa = G$Taxa, G[, snp_ordered, drop = FALSE])

# Build sample x SNP matrix and split into LOW/HIGH
SNPmat <- as.matrix(G_sub[, -1, drop = FALSE])   # numeric 0/1/2 expected
rownames(SNPmat) <- G_sub$Taxa

LOW_mat  <- SNPmat[LOW_ids,  , drop = FALSE]
HIGH_mat <- SNPmat[HIGH_ids, , drop = FALSE]

# Sum allele counts per tail & build BSA table 
# ALT = sum of genotypes; REF = 2 * (# non-NA samples) - ALT
ALT_LOW   <- colSums(LOW_mat,  na.rm = TRUE)
ALT_HIGH  <- colSums(HIGH_mat, na.rm = TRUE)
nLOW      <- colSums(!is.na(LOW_mat))
nHIGH     <- colSums(!is.na(HIGH_mat))
REF_LOW   <- 2L * nLOW  - ALT_LOW
REF_HIGH  <- 2L * nHIGH - ALT_HIGH

# quick sanity
stopifnot(all(REF_LOW  + ALT_LOW  == 2L * nLOW))
stopifnot(all(REF_HIGH + ALT_HIGH == 2L * nHIGH))

# build table in the shape QTLseqR expects 
BSA <- data.table::data.table(
  CHROM        = snpmap$Chromosome[snpmap$SNP %in% colnames(SNPmat)],
  POS          = as.integer(snpmap$Position[snpmap$SNP %in% colnames(SNPmat)]),
  REF          = NA_character_,
  ALT          = NA_character_,
  AD_REF.LOW   = as.integer(REF_LOW),
  AD_ALT.LOW   = as.integer(ALT_LOW),
  AD_REF.HIGH  = as.integer(REF_HIGH),
  AD_ALT.HIGH  = as.integer(ALT_HIGH)
)

# optional: add pseudo-depths (useful for filtering/plots)
BSA[, DP.LOW  := AD_REF.LOW  + AD_ALT.LOW ]
BSA[, DP.HIGH := AD_REF.HIGH + AD_ALT.HIGH]

# Write file for QTLseqr 
outfile <- "BSA_input_table.tsv"
data.table::fwrite(
  BSA[, .(CHROM, POS, REF, ALT, AD_REF.LOW, AD_ALT.LOW, AD_REF.HIGH, AD_ALT.HIGH, DP.LOW, DP.HIGH)],
  file = outfile, sep = "\t"
)

# Import the summed table into QTLseqR
QTL <- importFromTable(
  file     = "BSA_input_table.tsv",
  highBulk = "HIGH",
  lowBulk  = "LOW",
  sep      = "\t"
)

# SNP quality visualization ----
library(ggplot2)
library(gridExtra)

p1 <- ggplot(data = QTL) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW), bins = 50, fill = "gray30") +
  xlim(0, 1000) +
  ggtitle("read depth")

p2 <- ggplot(data = QTL) +
  geom_histogram(aes(x = SNPindex.HIGH), bins = 50, fill = "gray30") +
  ggtitle("per-bulk SNP-index")

gridExtra::grid.arrange(p1, p2, ncol = 2)

# SNP filtering 
QTLfilt <- filterSNPs(
  SNPset          = QTL,
  refAlleleFreq   = 0.20,
  minTotalDepth   = 100,
  maxTotalDepth   = 400,
  depthDifference = 100,
  minSampleDepth  = 40,
  verbose         = TRUE
)

# QTL-seq analysis 
# make sure POS is integer & sorted (safe even if already done)
QTLfilt$POS <- as.integer(QTLfilt$POS)
QTLfilt <- QTLfilt[order(QTLfilt$CHROM, QTLfilt$POS), ]

QTLseq <- runQTLseqAnalysis(
  SNPset       = QTLfilt,
  windowSize   = 1e6,        # 1 Mb
  popStruc     = "F2",
  bulkSize     = bulkSize,   # c(nLOW, nHIGH) from Step 3
  replications = 10000,
  intervals    = c(95, 99)
)

chr_map <- c(
  "VaccDscaff1"  = "Chr1",
  "VaccDscaff2"  = "Chr2",
  "VaccDscaff4"  = "Chr3",
  "VaccDscaff6"  = "Chr4",
  "VaccDscaff7"  = "Chr5",
  "VaccDscaff11" = "Chr6",
  "VaccDscaff12" = "Chr7",
  "VaccDscaff13" = "Chr8",
  "VaccDscaff17" = "Chr9",
  "VaccDscaff20" = "Chr10",
  "VaccDscaff21" = "Chr11",
  "VaccDscaff22" = "Chr12"
)

# apply mapping to your dataset (QTLseq object)
QTLseq$CHROM <- ifelse(QTLseq$CHROM %in% names(chr_map),
                       chr_map[QTLseq$CHROM],
                       QTLseq$CHROM)

# Plot Δ(SNP-index) with 95%/99% CI bands
plotQTLStats(
  SNPset        = QTLseq,
  var           = "deltaSNP",
  plotIntervals = TRUE
) + 
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  ggtitle("Δ(SNP-index) across genome")

#######################################################################################
############################# TRAIT #2 ################################################
#######################################################################################

# Define bulks (LOW/HIGH) from phenotype 
trait     <- "V2"     # change to V2/V3/V4 if you want a different trait
tail_prop <- 0.15     # bottom/top 15%

# ensure sample order matches between G and PH (Taxa column)
stopifnot(identical(G$Taxa, PH$Taxa))

# order by trait and pick tails
PH_ord     <- PH[order(PH[[trait]]), ]
n_per_bulk <- floor(nrow(PH_ord) * tail_prop)

LOW_ids  <- PH_ord$Taxa[1:n_per_bulk]
HIGH_ids <- PH_ord$Taxa[(nrow(PH_ord) - n_per_bulk + 1):nrow(PH_ord)]

# sizes used later by runQTLseqAnalysis()
bulkSize <- c(length(LOW_ids), length(HIGH_ids))

# Align SNPs to the map and subset
# Figure out which SNPs are shared and put them in MAP order
snps_in_G   <- setdiff(names(G), "Taxa")
keep_snps   <- intersect(snpmap$SNP, snps_in_G)
snp_ordered <- snpmap$SNP[snpmap$SNP %in% keep_snps]  # preserves map order

# Subset G to those SNP columns (still keeping 'Taxa' first)
G_sub <- cbind(Taxa = G$Taxa, G[, snp_ordered, drop = FALSE])

# Build sample x SNP matrix and split into LOW/HIGH
SNPmat <- as.matrix(G_sub[, -1, drop = FALSE])   # numeric 0/1/2 expected
rownames(SNPmat) <- G_sub$Taxa

LOW_mat  <- SNPmat[LOW_ids,  , drop = FALSE]
HIGH_mat <- SNPmat[HIGH_ids, , drop = FALSE]

# Sum allele counts per tail & build BSA table
# ALT = sum of genotypes; REF = 2 * (# non-NA samples) - ALT
ALT_LOW   <- colSums(LOW_mat,  na.rm = TRUE)
ALT_HIGH  <- colSums(HIGH_mat, na.rm = TRUE)
nLOW      <- colSums(!is.na(LOW_mat))
nHIGH     <- colSums(!is.na(HIGH_mat))
REF_LOW   <- 2L * nLOW  - ALT_LOW
REF_HIGH  <- 2L * nHIGH - ALT_HIGH

# quick sanity
stopifnot(all(REF_LOW  + ALT_LOW  == 2L * nLOW))
stopifnot(all(REF_HIGH + ALT_HIGH == 2L * nHIGH))

# build table in the shape QTLseqR expects
BSA <- data.table::data.table(
  CHROM        = snpmap$Chromosome[snpmap$SNP %in% colnames(SNPmat)],
  POS          = as.integer(snpmap$Position[snpmap$SNP %in% colnames(SNPmat)]),
  REF          = NA_character_,
  ALT          = NA_character_,
  AD_REF.LOW   = as.integer(REF_LOW),
  AD_ALT.LOW   = as.integer(ALT_LOW),
  AD_REF.HIGH  = as.integer(REF_HIGH),
  AD_ALT.HIGH  = as.integer(ALT_HIGH)
)

# optional: add pseudo-depths (useful for filtering/plots)
BSA[, DP.LOW  := AD_REF.LOW  + AD_ALT.LOW ]
BSA[, DP.HIGH := AD_REF.HIGH + AD_ALT.HIGH]

# write the file we’ll feed to QTLseqr next
outfile <- "BSA_input_table.tsv"
data.table::fwrite(
  BSA[, .(CHROM, POS, REF, ALT, AD_REF.LOW, AD_ALT.LOW, AD_REF.HIGH, AD_ALT.HIGH, DP.LOW, DP.HIGH)],
  file = outfile, sep = "\t"
)

# Import the summed table into QTLseqR 
QTL <- importFromTable(
  file     = "BSA_input_table.tsv",
  highBulk = "HIGH",
  lowBulk  = "LOW",
  sep      = "\t"
)

# SNP quality visualization 
library(ggplot2)
library(gridExtra)

p1 <- ggplot(data = QTL) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW), bins = 50, fill = "gray30") +
  xlim(0, 1000) +
  ggtitle("read depth")

p2 <- ggplot(data = QTL) +
  geom_histogram(aes(x = SNPindex.HIGH), bins = 50, fill = "gray30") +
  ggtitle("per-bulk SNP-index")

gridExtra::grid.arrange(p1, p2, ncol = 2)

# SNP filtering 
QTLfilt <- filterSNPs(
  SNPset          = QTL,
  refAlleleFreq   = 0.20,
  minTotalDepth   = 100,
  maxTotalDepth   = 400,
  depthDifference = 100,
  minSampleDepth  = 40,
  verbose         = TRUE
)

# Step 8: QTL-seq analysis 
# make sure POS is integer & sorted (safe even if already done)
QTLfilt$POS <- as.integer(QTLfilt$POS)
QTLfilt <- QTLfilt[order(QTLfilt$CHROM, QTLfilt$POS), ]

QTLseq <- runQTLseqAnalysis(
  SNPset       = QTLfilt,
  windowSize   = 1e6,        # 1 Mb
  popStruc     = "F2",
  bulkSize     = bulkSize,   # c(nLOW, nHIGH) from Step 3
  replications = 10000,
  intervals    = c(95, 99)
)

chr_map <- c(
  "VaccDscaff1"  = "Chr1",
  "VaccDscaff2"  = "Chr2",
  "VaccDscaff4"  = "Chr3",
  "VaccDscaff6"  = "Chr4",
  "VaccDscaff7"  = "Chr5",
  "VaccDscaff11" = "Chr6",
  "VaccDscaff12" = "Chr7",
  "VaccDscaff13" = "Chr8",
  "VaccDscaff17" = "Chr9",
  "VaccDscaff20" = "Chr10",
  "VaccDscaff21" = "Chr11",
  "VaccDscaff22" = "Chr12"
)

# apply mapping to your dataset (QTLseq object)
QTLseq$CHROM <- ifelse(QTLseq$CHROM %in% names(chr_map),
                       chr_map[QTLseq$CHROM],
                       QTLseq$CHROM)
# Plot Δ(SNP-index) with 95%/99% CI bands
plotQTLStats(
  SNPset        = QTLseq,
  var           = "deltaSNP",
  plotIntervals = TRUE
) + 
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  ggtitle("Δ(SNP-index) across genome")

#######################################################################################
############################# TRAIT #3 ################################################
#######################################################################################

# Define bulks (LOW/HIGH) from phenotype ----
trait     <- "V3"     # change to V2/V3/V4 if you want a different trait
tail_prop <- 0.15     # bottom/top 15%

# ensure sample order matches between G and PH (Taxa column)
stopifnot(identical(G$Taxa, PH$Taxa))

# order by trait and pick tails
PH_ord     <- PH[order(PH[[trait]]), ]
n_per_bulk <- floor(nrow(PH_ord) * tail_prop)

LOW_ids  <- PH_ord$Taxa[1:n_per_bulk]
HIGH_ids <- PH_ord$Taxa[(nrow(PH_ord) - n_per_bulk + 1):nrow(PH_ord)]

# sizes used later by runQTLseqAnalysis()
bulkSize <- c(length(LOW_ids), length(HIGH_ids))

# Align SNPs to the map and subset ----
# Figure out which SNPs are shared and put them in MAP order
snps_in_G   <- setdiff(names(G), "Taxa")
keep_snps   <- intersect(snpmap$SNP, snps_in_G)
snp_ordered <- snpmap$SNP[snpmap$SNP %in% keep_snps]  # preserves map order

# Subset G to those SNP columns (still keeping 'Taxa' first)
G_sub <- cbind(Taxa = G$Taxa, G[, snp_ordered, drop = FALSE])

# Build sample x SNP matrix and split into LOW/HIGH
SNPmat <- as.matrix(G_sub[, -1, drop = FALSE])   # numeric 0/1/2 expected
rownames(SNPmat) <- G_sub$Taxa

LOW_mat  <- SNPmat[LOW_ids,  , drop = FALSE]
HIGH_mat <- SNPmat[HIGH_ids, , drop = FALSE]

# Sum allele counts per tail & build BSA table 
# ALT = sum of genotypes; REF = 2 * (# non-NA samples) - ALT
ALT_LOW   <- colSums(LOW_mat,  na.rm = TRUE)
ALT_HIGH  <- colSums(HIGH_mat, na.rm = TRUE)
nLOW      <- colSums(!is.na(LOW_mat))
nHIGH     <- colSums(!is.na(HIGH_mat))
REF_LOW   <- 2L * nLOW  - ALT_LOW
REF_HIGH  <- 2L * nHIGH - ALT_HIGH

# quick sanity
stopifnot(all(REF_LOW  + ALT_LOW  == 2L * nLOW))
stopifnot(all(REF_HIGH + ALT_HIGH == 2L * nHIGH))

# build table in the shape QTLseqR expects
BSA <- data.table::data.table(
  CHROM        = snpmap$Chromosome[snpmap$SNP %in% colnames(SNPmat)],
  POS          = as.integer(snpmap$Position[snpmap$SNP %in% colnames(SNPmat)]),
  REF          = NA_character_,
  ALT          = NA_character_,
  AD_REF.LOW   = as.integer(REF_LOW),
  AD_ALT.LOW   = as.integer(ALT_LOW),
  AD_REF.HIGH  = as.integer(REF_HIGH),
  AD_ALT.HIGH  = as.integer(ALT_HIGH)
)

# optional: add pseudo-depths (useful for filtering/plots)
BSA[, DP.LOW  := AD_REF.LOW  + AD_ALT.LOW ]
BSA[, DP.HIGH := AD_REF.HIGH + AD_ALT.HIGH]

# write the file we’ll feed to QTLseqr next
outfile <- "BSA_input_table.tsv"
data.table::fwrite(
  BSA[, .(CHROM, POS, REF, ALT, AD_REF.LOW, AD_ALT.LOW, AD_REF.HIGH, AD_ALT.HIGH, DP.LOW, DP.HIGH)],
  file = outfile, sep = "\t"
)

# Import the summed table into QTLseqR 
QTL <- importFromTable(
  file     = "BSA_input_table.tsv",
  highBulk = "HIGH",
  lowBulk  = "LOW",
  sep      = "\t"
)

# SNP quality visualization 
library(ggplot2)
library(gridExtra)

p1 <- ggplot(data = QTL) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW), bins = 50, fill = "gray30") +
  xlim(0, 1000) +
  ggtitle("read depth")

p2 <- ggplot(data = QTL) +
  geom_histogram(aes(x = SNPindex.HIGH), bins = 50, fill = "gray30") +
  ggtitle("per-bulk SNP-index")

gridExtra::grid.arrange(p1, p2, ncol = 2)

# SNP filtering 
QTLfilt <- filterSNPs(
  SNPset          = QTL,
  refAlleleFreq   = 0.20,
  minTotalDepth   = 100,
  maxTotalDepth   = 400,
  depthDifference = 100,
  minSampleDepth  = 40,
  verbose         = TRUE
)

# QTL-seq analysis 
# make sure POS is integer & sorted (safe even if already done)
QTLfilt$POS <- as.integer(QTLfilt$POS)
QTLfilt <- QTLfilt[order(QTLfilt$CHROM, QTLfilt$POS), ]

QTLseq <- runQTLseqAnalysis(
  SNPset       = QTLfilt,
  windowSize   = 1e6,        # 1 Mb, like Felipe
  popStruc     = "F2",
  bulkSize     = bulkSize,   # c(nLOW, nHIGH) from Step 3
  replications = 10000,
  intervals    = c(95, 99)
)

chr_map <- c(
  "VaccDscaff1"  = "Chr1",
  "VaccDscaff2"  = "Chr2",
  "VaccDscaff4"  = "Chr3",
  "VaccDscaff6"  = "Chr4",
  "VaccDscaff7"  = "Chr5",
  "VaccDscaff11" = "Chr6",
  "VaccDscaff12" = "Chr7",
  "VaccDscaff13" = "Chr8",
  "VaccDscaff17" = "Chr9",
  "VaccDscaff20" = "Chr10",
  "VaccDscaff21" = "Chr11",
  "VaccDscaff22" = "Chr12"
)

# apply mapping to your dataset (QTLseq object)
QTLseq$CHROM <- ifelse(QTLseq$CHROM %in% names(chr_map),
                       chr_map[QTLseq$CHROM],
                       QTLseq$CHROM)

# Plot Δ(SNP-index) with 95%/99% CI bands
plotQTLStats(
  SNPset        = QTLseq,
  var           = "deltaSNP",
  plotIntervals = TRUE
) + 
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  ggtitle("Δ(SNP-index) across genome")

#######################################################################################
############################# TRAIT #4 ################################################
#######################################################################################

# Define bulks (LOW/HIGH) from phenotype 
trait     <- "V4"     # change to V2/V3/V4 if you want a different trait
tail_prop <- 0.15     # bottom/top 15%

# Ensure sample order matches between G and PH (Taxa column)
stopifnot(identical(G$Taxa, PH$Taxa))

# order by trait and pick tails
PH_ord     <- PH[order(PH[[trait]]), ]
n_per_bulk <- floor(nrow(PH_ord) * tail_prop)

LOW_ids  <- PH_ord$Taxa[1:n_per_bulk]
HIGH_ids <- PH_ord$Taxa[(nrow(PH_ord) - n_per_bulk + 1):nrow(PH_ord)]

# sizes used later by runQTLseqAnalysis()
bulkSize <- c(length(LOW_ids), length(HIGH_ids))

# Align SNPs to the map and subset ----
# Figure out which SNPs are shared and put them in MAP order
snps_in_G   <- setdiff(names(G), "Taxa")
keep_snps   <- intersect(snpmap$SNP, snps_in_G)
snp_ordered <- snpmap$SNP[snpmap$SNP %in% keep_snps]  # preserves map order

# Subset G to those SNP columns (still keeping 'Taxa' first)
G_sub <- cbind(Taxa = G$Taxa, G[, snp_ordered, drop = FALSE])

# Build sample x SNP matrix and split into LOW/HIGH
SNPmat <- as.matrix(G_sub[, -1, drop = FALSE])   # numeric 0/1/2 expected
rownames(SNPmat) <- G_sub$Taxa

LOW_mat  <- SNPmat[LOW_ids,  , drop = FALSE]
HIGH_mat <- SNPmat[HIGH_ids, , drop = FALSE]

# Sum allele counts per tail & build BSA table 
# ALT = sum of genotypes; REF = 2 * (# non-NA samples) - ALT
ALT_LOW   <- colSums(LOW_mat,  na.rm = TRUE)
ALT_HIGH  <- colSums(HIGH_mat, na.rm = TRUE)
nLOW      <- colSums(!is.na(LOW_mat))
nHIGH     <- colSums(!is.na(HIGH_mat))
REF_LOW   <- 2L * nLOW  - ALT_LOW
REF_HIGH  <- 2L * nHIGH - ALT_HIGH

# quick sanity
stopifnot(all(REF_LOW  + ALT_LOW  == 2L * nLOW))
stopifnot(all(REF_HIGH + ALT_HIGH == 2L * nHIGH))

# build table in the shape QTLseqR expects 
BSA <- data.table::data.table(
  CHROM        = snpmap$Chromosome[snpmap$SNP %in% colnames(SNPmat)],
  POS          = as.integer(snpmap$Position[snpmap$SNP %in% colnames(SNPmat)]),
  REF          = NA_character_,
  ALT          = NA_character_,
  AD_REF.LOW   = as.integer(REF_LOW),
  AD_ALT.LOW   = as.integer(ALT_LOW),
  AD_REF.HIGH  = as.integer(REF_HIGH),
  AD_ALT.HIGH  = as.integer(ALT_HIGH)
)

# optional: add pseudo-depths (useful for filtering/plots)
BSA[, DP.LOW  := AD_REF.LOW  + AD_ALT.LOW ]
BSA[, DP.HIGH := AD_REF.HIGH + AD_ALT.HIGH]

# write the file we’ll feed to QTLseqr next
outfile <- "BSA_input_table.tsv"
data.table::fwrite(
  BSA[, .(CHROM, POS, REF, ALT, AD_REF.LOW, AD_ALT.LOW, AD_REF.HIGH, AD_ALT.HIGH, DP.LOW, DP.HIGH)],
  file = outfile, sep = "\t"
)

# Import the summed table into QTLseqR 
QTL <- importFromTable(
  file     = "BSA_input_table.tsv",
  highBulk = "HIGH",
  lowBulk  = "LOW",
  sep      = "\t"
)

# SNP quality visualization 
library(ggplot2)
library(gridExtra)

p1 <- ggplot(data = QTL) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW), bins = 50, fill = "gray30") +
  xlim(0, 1000) +
  ggtitle("read depth")

p2 <- ggplot(data = QTL) +
  geom_histogram(aes(x = SNPindex.HIGH), bins = 50, fill = "gray30") +
  ggtitle("per-bulk SNP-index")

gridExtra::grid.arrange(p1, p2, ncol = 2)

# SNP filtering 
QTLfilt <- filterSNPs(
  SNPset          = QTL,
  refAlleleFreq   = 0.20,
  minTotalDepth   = 100,
  maxTotalDepth   = 400,
  depthDifference = 100,
  minSampleDepth  = 40,
  verbose         = TRUE
)

# QTL-seq analysis 
# make sure POS is integer & sorted (safe even if already done)
QTLfilt$POS <- as.integer(QTLfilt$POS)
QTLfilt <- QTLfilt[order(QTLfilt$CHROM, QTLfilt$POS), ]

QTLseq <- runQTLseqAnalysis(
  SNPset       = QTLfilt,
  windowSize   = 1e6,        # 1 Mb, like Felipe
  popStruc     = "F2",
  bulkSize     = bulkSize,   # c(nLOW, nHIGH) from Step 3
  replications = 10000,
  intervals    = c(95, 99)
)

# Step: Rename scaffolds to chromosome numbers 
chr_map <- c(
  "VaccDscaff1"  = "Chr1",
  "VaccDscaff2"  = "Chr2",
  "VaccDscaff4"  = "Chr3",
  "VaccDscaff6"  = "Chr4",
  "VaccDscaff7"  = "Chr5",
  "VaccDscaff11" = "Chr6",
  "VaccDscaff12" = "Chr7",
  "VaccDscaff13" = "Chr8",
  "VaccDscaff17" = "Chr9",
  "VaccDscaff20" = "Chr10",
  "VaccDscaff21" = "Chr11",
  "VaccDscaff22" = "Chr12"
)

# apply mapping to your dataset (QTLseq object)
QTLseq$CHROM <- ifelse(QTLseq$CHROM %in% names(chr_map),
                       chr_map[QTLseq$CHROM],
                       QTLseq$CHROM)

# Plot Δ(SNP-index) with 95%/99% CI bands
plotQTLStats(
  SNPset        = QTLseq,
  var           = "deltaSNP",
  plotIntervals = TRUE
) + 
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  ggtitle("Δ(SNP-index) across genome")