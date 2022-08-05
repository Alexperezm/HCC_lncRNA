#                 02.Normalization + Derregulation Analysis - NASIR
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Gathering of all RNA counts, its normalization and Gene Expression
#Analysis.



# Construct the SummarizedExperiment object:

setwd("")

library(SummarizedExperiment)
library(tidyverse)

NASIR_P <- readRDS("NASIR_P.rds")
NASIR_T <- readRDS("NASIR_T.rds")

NASIR <- left_join(NASIR_P,NASIR_T,by="Gene",suffix=c("P","T"))
G <- NASIR$Gene
NASIR <- NASIR%>%
  select(-Gene)
rownames(NASIR) <- G
NASIR <- NASIR[1:58721,]
nrows <- dim(NASIR)[1]
ncols <- length(NASIR)
rownames(NASIR) <- G[1:58721]

counts <- (NASIR)
counts <- as.matrix(counts)

rowRanges <- GRanges(rep(c("chr1", "chr2"), c(58720,1)),
                     IRanges(floor(runif(58721, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 58721, TRUE),
                     feature_id=list(NASIR_P[1:58721,1]),
                     row.names(NASIR_P[1:58721,1]))

colData <- DataFrame(Type=c(rep("Peritumor",13),rep("Tumor",30)),
                     Patient=str_replace(colnames(NASIR),"[:alpha:]","")) #row.names=colnames(NASIR)


se <- SummarizedExperiment(assays=list(counts=counts), colData=colData,rowRanges = rowRanges)

colSums(assay(se))
head(assay(se), 3)
rowRanges(se)
str(metadata(rowRanges(se)))
colData(se)
class(se)

#Once developed the Summarized Experiment object, is time to start with the Gene
#Expression Analysis:

library("DESeq2")

# Order of factors before relevel
se$Type

#This experiment will be based on the study between patients and sample type
#(Peritumor vs. Tumor):
dds_raw <- DESeqDataSet(se, design = ~ Patient + Type)
colData(dds_raw)

head(assay(dds_raw))
rowRanges(dds_raw)
class(dds_raw)


#Visualization of raw counts:
head(counts(dds_raw))
rownames(dds_raw)
dds_raw

#pre-filtering of the counts, so as to remove considerable amount of
#insignificant records:
keep <- rowSums(counts(dds_raw)) >= 50
dds <- dds_raw[keep,]


dds <- DESeq(dds)
res <- results(dds)
res

#Normalization:
colSums(counts(dds_raw)) 
dds_raw <- estimateSizeFactors(dds_raw)
sizeFactors(dds_raw)

# counts()  allows to immediately retrieve the **normalized** read counts.
counts_normalized  <- counts(dds_raw, normalized = TRUE)
# Now take a look at the sum of the total depth after normalization 
colSums(counts(dds_raw, normalized=TRUE))

# Observe now the normalized counts
head(counts_normalized)

#Several normalization methods:

# Log^2^ normalization adding 1 pseudocount.
counts_log_normalized <- log2(counts_normalized + 1)

par(mfrow=c(1,2)) # to plot  the  following  two  images side by side  each  other

# first , boxplots  of non -transformed  read  counts (one  per  sample)
boxplot(counts_normalized , notch = TRUE , las=2, cex.axis = 0.7,
        main = "untransformed  read  counts", ylab = "read  counts")

# box  plots  of log^2^ -transformed  read  counts
boxplot(counts_log_normalized , notch = TRUE , las=2, cex.axis = 0.7,
        main = "log2 -transformed  read  counts",
        ylab = "log2(read  counts)")



library(vsn)
library(ggplot2)
library(hexbin)

# When the expected amount of variance is approximately the same across different
#mean values, the data is said to be homoskedastic. For RNA-seq counts, however,
#the expected variance grows with the mean.
SdPlot <- meanSdPlot(counts_normalized, ranks = FALSE, plot = FALSE)  
SdPlot$gg + ggtitle("sequencing depth normalized") + ylab("standard deviation")


# The logarithm with a small pseudocount amplifies differences when the values
#are close to 0. The low count genes with low signal-to-noise ratio will overly
#contribute to sample-sample distances and PCA plots.
SdPlot_log <- meanSdPlot(counts_log_normalized, ranks = FALSE, plot = FALSE)
SdPlot_log$gg + ggtitle("sequencing depth normalized log2(read counts)") + ylab("standard deviation")


# We can see this property of count with simulated data. Poison counts in a range
#of lambda 0.1 to 100
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)


# Standard deviation of each row against the mean on raw counts
sim_SdPlot <- meanSdPlot(cts, ranks = FALSE, plot = FALSE)
sim_SdPlot$gg + ggtitle("simulated data sequencing depth normalized)") + ylab("standard deviation")

# Standard deviation of each row against the mean on log2 pseudocount transformed data
log.cts.one <- log2(cts + 1)
sim_SdPlot_log <- meanSdPlot(log.cts.one, ranks = FALSE, plot = FALSE)
sim_SdPlot_log$gg + ggtitle("simulated data sequencing depth normalized log2(read counts)") + ylab("standard deviation")

#Variance stabilizing transformation
vsd <- vst(dds_raw, blind = FALSE)
head(assay(vsd), 3)

#rlog
rld <- rlog(dds_raw, blind = FALSE)
head(assay(rld), 3)

par(mfrow=c(2,2))
# Raw data normalized by sequencing depth
plot(counts_normalized [,1:2], cex=.1, main = "Normalized by sequencing depth")

# Log^2^ normalization adding 1 pseudocount.
plot(counts_log_normalized [,1:2], cex=.1, main = "Normalized log2(read counts)")

# rlog transformed
rlog_norm_counts <- assay(rld)
plot(rlog_norm_counts[,1:2], cex=.1, main = "rlog transformed", xlim=c(0,18), ylim=c(0,18))

# Variance stabilizing transformation
vsd_norm_counts <- assay(vsd)
plot(vsd_norm_counts[,1:2], cex=.1, main = "Variance stabilizing transformation", xlim=c(0,18), ylim=c(0,18))

#Clustering:
#Variance stabilizing transformation

library(pheatmap)
library(RColorBrewer)

# Sample distances with VST
sampleDists_vsd <- dist(t(assay(vsd)))
sampleDists_vsd

# Transform sample distances to matrix
sampleDistMatrix_vsd <- as.matrix( sampleDists_vsd )
# Modify row names joining treatment and cell type
rownames(sampleDistMatrix_vsd) <- paste( vsd$dex, vsd$cell, sep = " - " )
# Remove col names
colnames(sampleDistMatrix_vsd) <- NULL
# Colors palette
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# Draw heatmap
heatmap <- pheatmap(sampleDistMatrix_vsd,
                    clustering_distance_rows = sampleDists_vsd,
                    clustering_distance_cols = sampleDists_vsd,
                    col = colors)

# Plot PCA on Patient type
plotPCA(vsd, intgroup="Patient")

# Plot PCA on Type: Tumor/Peritumor
plotPCA(vsd, intgroup="Type")

# Plot PCA on Patient and Type treatment
plotPCA(vsd, intgroup=c("Patient", "Type"))


#rlog:

# Sample distances with rld
sampleDists_rld <- dist(t(assay(rld)))
sampleDists_rld

# Transform sample distances to matrix
sampleDistMatrix_rld <- as.matrix( sampleDists_rld )
# Modify row names joining treatment and cell type
rownames(sampleDistMatrix_rld) <- paste( rld$dex, rld$cell, sep = " - " )
# Remove col names
colnames(sampleDistMatrix_rld) <- NULL
# Colors palette
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# Draw heatmap
heatmap <- pheatmap(sampleDistMatrix_rld,
                    clustering_distance_rows = sampleDists_rld,
                    clustering_distance_cols = sampleDists_rld,
                    col = colors)

# Plot PCA on Patient type
plotPCA(rld, intgroup="Patient")

# Plot PCA on Type
plotPCA(rld, intgroup="Type")

# Plot PCA on cell type and dexamethasone treatment
plotPCA(rld, intgroup=c("Patient", "Type"))

#######
#Differential gene expression
#######

# Run DGE DeSeq2 pipeline
dds_DGE <- DESeq(dds_raw)

# Mean-dispersion relationship 
par(mfrow=c(1,1))
plotDispEsts(dds_DGE)

#Results table:
# Calling results() will build the base means across samples, log2 fold variation,
#standard errors, p value and p adjusted. 
dds_DGE_results <- results(dds_DGE)
head(dds_DGE_results)

# Result object can be filter like a data frame. We have 4023 significant genes
#with a p adjusted < 0.05.  
table(dds_DGE_results$padj < 0.05)

# rownames(subset(dds_DGE_results, padj < 0.01))

# We can also access to the metadata contained on the columns and check the
#pairwise contrast studied ()
mcols(dds_DGE_results, use.names = TRUE)

# Summary of DGE
summary(dds_DGE_results)

# Select significant genes with a p adjusted lower than 0.05
resSig <- subset(dds_DGE_results, padj < 0.05)
# Order significant genes by p adjusted 
resSig <- resSig[ order(resSig$padj), ]
# Order significant genes with the strongest down-regulation 
head(resSig[ order(resSig$log2FoldChange), ])

# And with the strongest up-regulation
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

saveRDS(resSig, file="resSig.rds")
saveRDS(rld, file="rld.rds")


#######
#Annotation of significant genes:
#######
library("AnnotationDbi")
library("org.Hs.eg.db")
library("org.Hs.egGENETYPE")
resSig%>%mutate(names (resSig))
resSig$SYMBOL = mapIds(org.Hs.eg.db,
                       keys=rownames(resSig), 
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")

# Add gene name column
resSig$GENENAME = mapIds(org.Hs.eg.db,
                         keys=rownames(resSig), 
                         column="GENENAME",
                         keytype="ENSEMBL",
                         multiVals="first")

# Add functional path column
resSig$PATH = mapIds(org.Hs.eg.db,
                     keys=rownames(resSig), 
                     column="PATH",
                     keytype="ENSEMBL",
                     multiVals="first")

# Add entrex ID column
resSig$ENTREZID = mapIds(org.Hs.eg.db,
                         keys=rownames(resSig), 
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
resSig$Type = mapIds(org.Hs.eg.db,
                       keys=rownames(resSig), 
                       column="GENETYPE",
                       keytype="ENSEMBL",
                       multiVals="first")


annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = NASIR$Gene,  # data to use for retrieval
                                           columns = c("SYMBOL", "ENTREZID","GENENAME"), # information to retreive for given data
                                           keytype = "ENSEMBL") # type of data given in 'keys' argument

# Add ENSEMBL column
resSig$ENSEMBL <- rownames(resSig)

resSig_DF <- as.data.frame(resSig)


# Save significant DGE results to file:
write.table(resSig, file = "DESeq2_Sig_results_tum_Vs_pertum.tab", sep = "\t", quote = FALSE , row.names = TRUE)

#Exploratory Plots:

# Individual gene raw count plot
plotCounts(dds_DGE, gene = "ENSG00000256811.1", intgroup=c("Type"))


par(mfrow=c(1,1))
# Histogram of p-values frequencies  
hist(dds_DGE_results$pvalue , col = "blue",  xlab = "", , border = "white", ylab = "Frequency", breaks =0:40/40, main = "frequencies of p-values")

par(mfrow=c(1,2))
# MA plot without shrinkage
plotMA(dds_DGE_results, alpha = 0.05,  main = "Peritumor vs. Tumor", ylim = c(-5,5))

library("apeglm")
# Check possible contrast to create the MA plot
resultsNames(dds_DGE)

res <- lfcShrink(dds_DGE, coef="Tumor vs Peritumor", type="apeglm")
plotMA(res, alpha = 0.05, main = "Tumor vs Peritumor log2 shrink",ylim = c(-5, 5))


#Volcano Plot:
library("dplyr")

dds_DGE_results <- dds_DGE_results[order(dds_DGE_results$padj),]
results_order <- as.data.frame(dplyr::mutate(as.data.frame(dds_DGE_results), sig=ifelse(dds_DGE_results$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(dds_DGE_results))


p = ggplot2::ggplot(results_order, ggplot2::aes(log2FoldChange, -log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Significant genes Peritumor vs Tumor")

p

#Heatmap:

# Select the 20 most variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 200)

# Transform to matrix and use vsd transformation
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Patient","Type")])
pheatmap(mat, annotation_col = anno,show_rownames = FALSE,show_colnames = FALSE)


save.image("02.NASIR_Differential_gene_expression_analysis.Rdata")