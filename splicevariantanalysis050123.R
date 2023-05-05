if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SGSeq")
BiocManager::install("GenomicRanges")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library(SGSeq)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# bams7709 <- data.frame()
# bams7709 <- (c("7709-LB-1Aligned.sortedByCoord.out.bam","7709-LB-2Aligned.sortedByCoord.out.bam", "7709-LB-3Aligned.sortedByCoord.out.bam", "7709-LB-4Aligned.sortedByCoord.out.bam", "7709-LB-5Aligned.sortedByCoord.out.bam", "7709-LB-6Aligned.sortedByCoord.out.bam", "7709-LB-7Aligned.sortedByCoord.out.bam"))
# bams7709 <- cbind(c("WT", "si5B1", "si5B2", "si1B1", "si1B4", "siScramble", "5BKO"), bams7709)
# 
# colnames(bams7709)[1] = "sample_name"
# colnames(bams7709)[2] = "file_bam"
# 
# bams7709 <- as.data.frame(bams7709)
# 
# path <- "\\\\umms-kjverhey-win.turbo.storage.umich.edu/umms-kjverhey/abonil_RNAseq_Files/7709-LB_RNAseq/7709-LB_aligned/bams"
# bams7709$file_bam <- file.path(path, bams7709$file_bam)
# 
# complete7709 <- getBamInfo(bams7709)
# saveRDS(complete7709, file = "complete7709.Rds")
# 
# 
# bams7804 <- data.frame()
# bams7804 <- (c("7804-LB-1Aligned.sortedByCoord.out.bam","7804-LB-2Aligned.sortedByCoord.out.bam", "7804-LB-3Aligned.sortedByCoord.out.bam", "7804-LB-4Aligned.sortedByCoord.out.bam", "7804-LB-5Aligned.sortedByCoord.out.bam", "7804-LB-6Aligned.sortedByCoord.out.bam", "7804-LB-7Aligned.sortedByCoord.out.bam"))
# bams7804 <- cbind(c("WT", "si5B1", "si5B2", "si1B1", "si1B4", "siScramble", "5BKO"), bams7804)
# 
# colnames(bams7804)[1] = "sample_name"
# colnames(bams7804)[2] = "file_bam"
# 
# bams7804 <- as.data.frame(bams7804)
# 
# path <-"\\\\umms-kjverhey-win.turbo.storage.umich.edu/umms-kjverhey/abonil_RNAseq_Files/7804-LB_RNAseq/7804-LB_aligned/bams"
# bams7804$file_bam <- file.path (path, bams7804$file_bam)
# 
# complete7804 <- getBamInfo(bams7804)
# saveRDS(complete7804, file = "complete7804.Rds")


# Start here 

complete7709 <- readRDS(file = "complete7709.Rds")
complete7804 <- readRDS(file = "complete7804.Rds")

# Reorder rows 

complete7709 <- complete7709[c(1,7,4,5,2,3,6),]
complete7804 <- complete7804[c(1,7,4,5,2,3,6),]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- keepSeqlevels(txdb, "chr1")
seqlevelsStyle(txdb) <- "NCBI"

edb <- EnsDb.Hsapiens.v86
KIF1B <- genes(edb, filter = GeneNameFilter("KIF1B"))


txf_ucsc <- convertToTxFeatures(txdb)
KIF1B@seqinfo@genome <- "GRCh38.p14"
txf_ucsc <- txf_ucsc[txf_ucsc %over% KIF1B]

sgf_ucsc <- convertToSGFeatures(txf_ucsc)

# Changing path of bam file, cannot find with given path 

sgf_ucsc_7709 <- analyzeFeatures(complete7709, features = txf_ucsc)
sgf_ucsc_7804 <- analyzeFeatures(complete7804, features = txf_ucsc)

df <- plotFeatures(sgf_ucsc_7709, geneID = 1)
df.2 <- plotFeatures(sgf_ucsc_7804, geneID = 1)

sgfc_pred <- analyzeFeatures(complete7709, which = KIF1B)
sgfc_pred <- annotate(sgfc_pred, txf_ucsc)
sgfc_pred.2 <- analyzeFeatures(complete7804, which = KIF1B)
sgfc_pred.2 <- annotate(sgfc_pred.2, txf_ucsc)


df7709 <- plotFeatures(sgfc_pred, geneID = 1, color_novel = "red", Rowv = NA, main = "De Novo Prediction 7709 replicate")
df7804 <- plotFeatures(sgfc_pred.2, geneID = 1, color_novel = "red", Rowv = NA, main = "De Novo Prediction 7804 replicate")

sgvc_pred <- analyzeVariants(sgfc_pred)
sgvc_pred.2 <- analyzeVariants(sgfc_pred.2)
replicate7709 <- sgvc_pred
replicate7804 <- sgvc_pred.2
variantFreq(replicate7709)
variantFreq(replicate7804)

plotVariants(sgvc_pred, eventID = 1, color_novel = "red", Rowv = NA, main = "7709 replicate")
plotVariants(sgvc_pred.2, eventID= 1, color_novel = "red", Rowv = NA, main = "7804 replicate")

seqlevelsStyle(Hsapiens) <- "NCBI"

vep <- predictVariantEffects(sgvc_pred, txdb, Hsapiens)
