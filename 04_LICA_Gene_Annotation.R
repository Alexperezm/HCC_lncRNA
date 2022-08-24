#                         12. Gene annotation - LICA
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Annotation of the normalized LICA gene expressiona dataset.
#Including the symbol, gene name and gene type; so as to filter the dataset for
#posterior procedures.

setwd("")
library(readxl)
library(tidyverse)


library(readxl)

LICA_norm <- read_excel("normalizado_sin_ceros_LICA.xlsx")

dim(LICA_norm)
colnames(LICA_norm)
LICA_norm$gene_id <- str_sub(LICA_norm$gene_id,end = 15)



#ANNOTATE:
library("AnnotationDbi")
library("org.Hs.eg.db")
LICA_norm$SYMBOL = mapIds(org.Hs.eg.db,
                                            keys=LICA_norm$gene_id, 
                                            column="SYMBOL",
                                            keytype="ENSEMBL",
                                            multiVals="first")

# Add gene name column
LICA_norm$GENENAME = mapIds(org.Hs.eg.db,
                                              keys=LICA_norm$gene_id, 
                                              column="GENENAME",
                                              keytype="ENSEMBL",
                                              multiVals="first")

LICA_norm$GENETYPE = mapIds(org.Hs.eg.db,
                                              keys=LICA_norm$gene_id, 
                                              column="GENETYPE",
                                              keytype="ENSEMBL",
                                              multiVals="first")

write.csv(LICA_norm,"LICA_normalized_genetype.csv",row.names = TRUE)
