#                         11. Normalization - LICA
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: LICA gene expression data set normalization script.

library(readr)
library(tibble)
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)

setwd("")

just.raw.counts <-read.delim(file="LICA-read-counts_no_cero.csv",sep = ";", row.names = 1)
just.raw.counts <- data.matrix(just.raw.counts)
head(just.raw.counts)
dim(just.raw.counts)

meta.data = read.delim(file="meta_data.csv",sep = ";", row.names = 1)

head(meta.data)

count.data.set <- DESeqDataSetFromMatrix(countData=just.raw.counts, 
                                         colData=meta.data, design= ~ 1) 

count.data.set.object <- DESeq(count.data.set)


vsd <- vst(count.data.set.object)

norm.data = assay(vsd)

write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off_no_cero.txt", row.names=TRUE,col.names=NA,quote=FALSE)
