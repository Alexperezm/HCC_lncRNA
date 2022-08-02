#                  07. Gene Annotation & Gene Ontology - NASIR
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Datasets will be enriched with key data as: genetype, genename,
#symbol... So, in following studies datasets will be filtered based on genetype,
#and enrichment analysis will be developed.

setwd("")
library(readxl)
library(tidyverse)

R_NICO100Table <- read_excel("R-NICO100Table.xlsx")
HCCDeReg <- read_excel("HCCDeReg.xlsx")
nDereg <- HCCDeReg[,1]

#MAnipulate gene_id, for dataset joining:
NASIR_FINAL_TABLE_COUNTS <- read_xlsx("CUENTAS_TUMORES_PERITUMORES.xlsx")
dim(NASIR_FINAL_TABLE_COUNTS)
NASIR_FINAL_TABLE_COUNTS <- as.data.frame(t(NASIR_FINAL_TABLE_COUNTS))
GENE <- rownames(NASIR_FINAL_TABLE_COUNTS)
NASIR_FINAL_TABLE_COUNTS$GENE <- GENE
num <- str_sub(NASIR_FINAL_TABLE_COUNTS$GENE,start=6,end = 15)
numerico <- as.integer(num)
NASIR_FINAL_TABLE_COUNTS$NUM <- numerico

#Extract gene_id from R100:
num100 <- str_sub(R_NICO100Table$ENSEMBLEID,start=6,end = 15)
numerico100 <- as.integer(num100)
R_NICO100Table$NUM <- numerico100

#Matching both datasets:
match <- inner_join(NASIR_FINAL_TABLE_COUNTS,R_NICO100Table,by=("NUM"))



################################################################################
#Same procedures, in this case using the other dataset:HCCDeReg.

num <- str_sub(HCCDeReg$gene_id,start=6,end = 15)
numerico <- as.integer(num)
HCCDeReg$NUM <- numerico


match_De_Reg <- inner_join(NASIR_FINAL_TABLE_COUNTS,HCCDeReg,by=("NUM"))


match_filt_de_reg <- t(match_De_Reg)
colnames(match_filt_de_reg) <- match_filt_de_reg[1,]
match_filt_de_reg <- match_filt_de_reg[2:31,]
match_filt_de_reg <- as.data.frame(match_filt_de_reg)
summary(match_filt_de_reg)
match_filt_de_reg <- as.double(match_filt_de_reg)
match_filt_de_reg <- as.data.frame(sapply(match_filt_de_reg,as.numeric))
str(match_filt_de_reg)

summary(match_filt_de_reg)      


################################################################################
################################################################################
#Gene annotation:
library(readxl)
CUENTAS_TUMORES_PERITUMORES <- read_excel("CUENTAS_TUMORES_PERITUMORES.xlsx")

dim(CUENTAS_TUMORES_PERITUMORES)
colnames(CUENTAS_TUMORES_PERITUMORES)
CUENTAS_TUMORES_PERITUMORES$GENE <- str_sub(CUENTAS_TUMORES_PERITUMORES$GENE,end = 15)

library("AnnotationDbi")
library("org.Hs.eg.db")

#Add gene symbol column:
CUENTAS_TUMORES_PERITUMORES$SYMBOL = mapIds(org.Hs.eg.db,
                                            keys=CUENTAS_TUMORES_PERITUMORES$GENE, 
                                            column="SYMBOL",
                                            keytype="ENSEMBL",
                                            multiVals="first")

# Add gene name column:
CUENTAS_TUMORES_PERITUMORES$GENENAME = mapIds(org.Hs.eg.db,
                                              keys=CUENTAS_TUMORES_PERITUMORES$GENE, 
                                              column="GENENAME",
                                              keytype="ENSEMBL",
                                              multiVals="first")

# Add gene type column:
CUENTAS_TUMORES_PERITUMORES$GENETYPE = mapIds(org.Hs.eg.db,
                                              keys=CUENTAS_TUMORES_PERITUMORES$GENE, 
                                              column="GENETYPE",
                                              keytype="ENSEMBL",
                                              multiVals="first")

write.csv(CUENTAS_TUMORES_PERITUMORES,"CUENTAS_TUMORES_PERITUMORES_genetype.csv",row.names = TRUE)


#Same procedure for the filtered dataset:
library(readxl)
Nasir_filtered <- read_excel("Nasir filtered.xlsx", sheet = "Hoja1 (2)")

dim(Nasir_filtered)
colnames(Nasir_filtered)
Nasir_filtered$GENE <- str_sub(Nasir_filtered$GENE,end = 15)



#ANNOTATE:
library("AnnotationDbi")
library("org.Hs.eg.db")
Nasir_filtered$SYMBOL = mapIds(org.Hs.eg.db,
                               keys=Nasir_filtered$GENE, 
                               column="SYMBOL",
                               keytype="ENSEMBL",
                               multiVals="first")

# Add gene name column
Nasir_filtered$GENENAME = mapIds(org.Hs.eg.db,
                                 keys=Nasir_filtered$GENE, 
                                 column="GENENAME",
                                 keytype="ENSEMBL",
                                 multiVals="first")

Nasir_filtered$GENETYPE = mapIds(org.Hs.eg.db,
                                 keys=Nasir_filtered$GENE, 
                                 column="GENETYPE",
                                 keytype="ENSEMBL",
                                 multiVals="first")

write.csv(Nasir_filtered,"Nasir_filtered_with_genetype.csv",row.names = TRUE)

################################################################################
################################################################################
#Gene ontology and gene clustering:

library("clusterProfiler")
library("enrichplot")
# BiocManager::install("ggnewscale")
library("ggnewscale")

OrgDb <- org.Hs.eg.db # can also be other organisms
# Get ENTREZID as vector. genes list are gonna be used to call GO terms
CUENTAS_TUMORES_PERITUMORES$GENE <- str_sub(CUENTAS_TUMORES_PERITUMORES$GENE,start=5,end = 15)
numerico <- as.integer(CUENTAS_TUMORES_PERITUMORES$GENE)
CUENTAS_TUMORES_PERITUMORES$GENE <- numerico
genes <- as.character(CUENTAS_TUMORES_PERITUMORES$GENE)


# GeneOntology terms

## Cellular compartment
# In clusterProfiler, groupGO is designed for gene classification based on GO
#distribution at a specific level.
ggo <- clusterProfiler::groupGO(gene     = genes,
                                OrgDb    = OrgDb,
                                ont      = "CC",     #MF , BP CC
                                level    = 3,
                                readable = TRUE)
head(as.data.frame(ggo)[,-5])

barplot(ggo, drop=TRUE, showCategory=30, vertex.label.cex=0.8)

# GO over-representation test:
ego <- clusterProfiler::enrichGO(gene          = genes,
                                 OrgDb         = OrgDb,
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.01, 
                                 readable      = TRUE)


barplot(ego, showCategory=30)


# Gene Concept Network
edox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

## Gene list with Fold Change data
geneList <- resSig$log2FoldChange
names(geneList) <- as.character(unique(resSig$ENTREZID))
geneList <- sort(geneList, decreasing = TRUE)
#de <- names(geneList)[abs(geneList) > 2]
head(geneList)

## categorySize can be scaled by 'pvalue' or 'geneNum'
p1 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)

# Enrichment Map
# Enrichment map organizes enriched terms into a network with edges connecting
#overlapping gene sets. 
# In this way, mutually overlapping gene sets are tend to cluster together,
#making it easy to identify functional module.
library(DOSE)

edox <- pairwise_termsim(edox)
p2 <- emapplot(edox, layout="kk")


# Plot Gene Concept Network and Enrichment Map
pdf(file=paste("Merge_GO_Gene_Set_Enrichment_Analysis", "CC", ".pdf", sep="_"),
    width=25, height=10)
cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

dev.off()