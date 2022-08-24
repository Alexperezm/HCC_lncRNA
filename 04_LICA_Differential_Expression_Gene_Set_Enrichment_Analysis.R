#         14_LICA_Differential_Expression_Gene set enrichment analysis:
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Gene expression processing and classification. The survival rate has
#been split into low and high, based on the first and fourth quantiles. The gene
#expression of those two groups will be studied in order to extract a curated
#list of genes directly related to the difference in survival rate. 

setwd("")

library(SummarizedExperiment)
library(tidyverse)
library(readxl)

#Gene counts and survival data is loaded an pre-processed:

LICA <- read_excel("Jessica_subread_counts.xlsx")
survival <- read_excel("survival.xlsx", na = "NA")
tnorm <- data.frame(t(LICA))
colnames(tnorm) <- tnorm[1,]
tnorm <- slice(tnorm,2:162)
patient_id <- survival[,1]

survival <- na.omit(survival)
surv <- survival[,2]
cuantiles<-quantile(surv, probs = c(.25,.50,.75,1),na.rm=TRUE)

#The survival data is split, and the first and fourth quantiles are extracted,
#only those patients will be studied.

L_survival <- survival%>%
  filter(`Last survival (days)`<= cuantiles["25%"])
L_survival["Group"] <- "L"

H_survival <- survival%>%
  filter(`Last survival (days)`>= cuantiles["75%"])
H_survival["Group"] <- "H"

survival <- bind_rows(L_survival,H_survival)

tnorm["Patient"] <- rownames(tnorm)
LICA_survival <- inner_join(tnorm,survival,by="Patient")

traspuesto <- data.frame(t(LICA_survival))
write.csv(traspuesto,"LICA_58K_quantile_sumarized_experiment_type.csv")

#Loading the final dataset, which already contains all the modifications:

Jesica <- read_excel("LICA_58K_quantile_sumarized_experiment_type.xlsx")
G <- Jesica$GENE
Jesica <- Jesica%>%
  select(-GENE)
rownames(Jesica) <- G
Jesica <- Jesica[1:58288,]
nrows <- dim(Jesica)[1]
ncols <- length(Jesica)
rownames(Jesica) <- G[1:58288]

counts <- (Jesica)
counts <- as.matrix(counts)

rowRanges <- GRanges(rep(c("chr1"), 58288),
                     IRanges(floor(runif(58288, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 58288, TRUE),
                     feature_id=list(Jesica[1:58288,1]),
                     row.names(Jesica[1:58288,1]))

grupos_sumarized_experiment_type <- read_excel("grupos_quantile_sumarized_experiment_type.xlsx")
grupos_sumarized_experiment_type["Patient"] <- str_replace(colnames(Jesica),"[:alpha:]","")
colData <- grupos_sumarized_experiment_type

#So as to carry the Differential Gene Expression Analysis, data is introduced
#into a Summarized Experiment object:

se <- SummarizedExperiment(assays=list(counts=counts), colData=colData,rowRanges = rowRanges)
saveRDS(se,"SummarizedExperiment_quantile_58k.rds")


#DESeq2:
#===============================================================================
#Here starts the analysis itself:
#Loading of the file, pre-processing and visualization, to check that everything
#is correct.


se <- readRDS("SummarizedExperiment_quantile_58k.rds")

colSums(assay(se))
head(assay(se), 3)
rowRanges(se)
str(metadata(rowRanges(se)))
colData(se)
class(se)

library("DESeq2")

# Order of factors before relevel
se$Group
#se$Group <- relevel(se$Group, "L")


dds_raw <- DESeqDataSet(se, design = ~ Group)
colData(dds_raw)

head(assay(dds_raw))
rowRanges(dds_raw)
class(dds_raw)


#Visualization of the raw data:

head(counts(dds_raw))
rownames(dds_raw)
dds_raw


#pre-filter, those rows that do not sum 50, will be removed as that counts will
#be insignificant:

keep <- rowSums(counts(dds_raw)) >= 50
dds <- dds_raw[keep,]


dds <- DESeq(dds)
res <- results(dds)
res


#Normalization, so as to analyze the pros&cons of each normalization procedure,
# will be using the log2 normalization and the vsd and rlog:

colSums(counts(dds_raw)) 
dds_raw <- estimateSizeFactors(dds_raw)
sizeFactors(dds_raw)

# counts()  allows to immediately retrieve the **normalized** read counts.
counts_normalized  <- counts(dds_raw, normalized = TRUE)
# Now take a look at the sum of the total depth after normalization 
colSums(counts(dds_raw, normalized=TRUE))

# Observe now the normalized counts
head(counts_normalized)

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

# When the expected amount of variance is approximately the same across different mean values, the data is said to be homoskedastic. For RNA-seq counts, however, the expected variance grows with the mean.

SdPlot <- meanSdPlot(counts_normalized, ranks = FALSE, plot = FALSE)  
SdPlot$gg + ggtitle("sequencing depth normalized") + ylab("standard deviation")


# The logarithm with a small pseudocount amplifies differences when the values are close to 0. The low count genes with low signal-to-noise ratio will overly contribute to sample-sample distances and PCA plots.

SdPlot_log <- meanSdPlot(counts_log_normalized, ranks = FALSE, plot = FALSE)
SdPlot_log$gg + ggtitle("sequencing depth normalized log2(read counts)") + ylab("standard deviation")


# We can see this property of count with simulated data. Poison counts in a range of lambda 0.1 to 100

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

#rlog: Regularized log transformation:

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
#Variance stabilizing transformation:

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


# Plot PCA on Type: Low / High Survival rate

plotPCA(vsd, intgroup="Group")


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

# Plot PCA on Type

plotPCA(rld, intgroup="Group")


#Differential gene expression:

# Run DGE DeSeq2 pipeline

dds_DGE <- DESeq(dds_raw)

# Mean-dispersion relationship 

par(mfrow=c(1,1))
plotDispEsts(dds_DGE)

#Results table:
# Calling results() will build the base means across samples, log2 fold variation, standard errors, p value and p adjusted. 

dds_DGE_results <- results(dds_DGE)
head(dds_DGE_results)

# Result object can be filter like a data frame. We have 4023 significant genes with a p adjusted < 0.05.  
table(dds_DGE_results$padj < 0.05)


# We can also access to the metadata contained on the columns and check the pairwise contrast studied ()
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

#The results obtained are saved so as to be accessible. So as not to repeat the
#whole procedure each time:
saveRDS(resSig, file="resSig_quantile_58k.rds")
saveRDS(rld, file="rld_quantile_58k.rds")

write.csv(resSig, "resultados_resSig.csv")

#ANNOTATE:
#===============================================================================
#Now, the results will be annotated, so as to know more about those relevant and
#significant genes; their symbol, gene name, gene type, entrezid and path will
#be included.

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")

library("AnnotationDbi")
library("org.Hs.eg.db")
library(stringr)



resSig_names <- rownames(resSig)
resSig_names <- str_sub(resSig_names,end = 15)

#resSig%>%mutate(names (resSig))
resSig$SYMBOL = mapIds(org.Hs.eg.db,
                       keys=resSig_names, #rownames(resSig)
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")

# Add gene name column
resSig$GENENAME = mapIds(org.Hs.eg.db,
                         keys=resSig_names, #rownames(resSig) 
                         column="GENENAME",
                         keytype="ENSEMBL",
                         multiVals="first")

# Add functional path column
resSig$PATH = mapIds(org.Hs.eg.db,
                     keys=resSig_names, #rownames(resSig)
                     column="PATH",
                     keytype="ENSEMBL",
                     multiVals="first")

# Add entrex ID column
resSig$ENTREZID = mapIds(org.Hs.eg.db,
                         keys=resSig_names, #rownames(resSig) 
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
resSig$Type = mapIds(org.Hs.eg.db,
                     keys=resSig_names, #rownames(resSig)
                     column="GENETYPE",
                     keytype="ENSEMBL",
                     multiVals="first")

saveRDS(resSig,"resSig_anotado_cuantiles58k.rds")
resSig$ENSEMBL <- rownames(resSig)
resSig_DF <- as.data.frame(resSig)


# Save significant DGE results to file

write.table(resSig, file = "DESeq2_Sig_results_anotado_cuantiles58k.tab", sep = "\t", quote = FALSE , row.names = TRUE)

#The annotated dataset is exported.

#Exploratory Plots:
#===============================================================================
#Developed in order to visualize and develop further studies
#around those significant genes:

par(mfrow=c(1,1))
# Histogram of p-values frequencies  
hist(dds_DGE_results$pvalue , col = "blue",  xlab = "", , border = "white", ylab = "Frequency", breaks =0:40/40, main = "frequencies of p-values")

library("apeglm")
# Check possible contrast to create the MA plot
resultsNames(dds_DGE)

res <- lfcShrink(dds_DGE, coef="Group_L_vs_H", type="apeglm")
plotMA(res, alpha = 0.05, main = "L vs H log2 shrink",ylim = c(-5, 5))


#Volcano Plot:

library("dplyr")

dds_DGE_results <- dds_DGE_results[order(dds_DGE_results$padj),]
results_order <- as.data.frame(dplyr::mutate(as.data.frame(dds_DGE_results), sig=ifelse(dds_DGE_results$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(dds_DGE_results))

p = ggplot2::ggplot(results_order, ggplot2::aes(log2FoldChange, -log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Significant genes Low vs High Survival")

p

#Heatmap:
# Select the 20 most variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 200)

# Transform to matrix and use vsd transformation
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group")])
rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno,show_rownames = TRUE,show_colnames = TRUE)


#GSEA:
#===============================================================================
#Gene set enrichment analysis:
#Identify classes of genes that might be over-represented in our whole gene set,
#and may have an association with our Low/High survival phenotypes.


library(readxl)
#BiocManager::install("clusterProfiler", version = "3.8")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(readr)
library(tidyverse)

setwd("")
resultados_resSig <- read_csv("resultados_resSig.csv")


#Prepare input:
# we want the log2 fold change 
original_gene_list <- resultados_resSig$log2FoldChange


# name the vector
names(original_gene_list) <- resultados_resSig$...1
names(original_gene_list) <-  str_sub(names(original_gene_list),start = 1,end = 15)


# omit any NA values 
gene_list<-na.omit(original_gene_list)


# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

library(clusterProfiler)


#Carry the gseGO analysis:
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             #maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
resultados <- gse@result[,1:10]

library(writexl)
write_xlsx(gse@result)


#Dotplot:
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


#Enrichment map
x2 <- pairwise_termsim(gse)
emapplot(x2,showCategory = 10,layout = "nicely")

emapplot(x2,
         showCategory = 10,
         color = "p.adjust",
         layout = "nicely",
         split = TRUE,
         pie = "equal",
         pie_scale = TRUE,
         cex_line = 1,
         min_edge = 0.2,
         cex_label_category = 1,
         node_label_size = TRUE,
)


#category network
#categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 6)


#ridgeplot (study the distribution of a numeric variable for several groups):
ridgeplot(gse,showCategory = 10) + labs(x = "enrichment distribution")


#GSEA plot:
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


# KEGG Enrichment Analysis:
#===============================================================================

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)


# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]


# Create a new dataframe df2 which has only the genes which were successfully
#mapped using the bitr function above:
df2 = resultados_resSig[names(original_gene_list) %in% dedup_ids$ENSEMBL,]


# Create a new column in df2 with the corresponding ENTREZ IDs
df2$trans = dedup_ids$ENTREZID


# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange


# Name vector with ENTREZ ids
names(kegg_gene_list) <- dedup_ids$ENTREZID


# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)


# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_gene_list
