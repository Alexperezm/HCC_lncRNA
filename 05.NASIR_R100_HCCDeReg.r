#                         05. R100/HCCDEREG - NASIR
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: A curated set of 100 genes will be studied, so in this script a new
#file will be created with the corresponding expression data for the required genes.

library(readxl)
R_NICO100Table <- read_excel("R-NICO100Table.xlsx", 
                             range = "c1:c101")

HCCDeReg <- read_excel("HCCDeReg.xlsx")
nDereg <- HCCDeReg[,1]


n <- rownames(NASIR)
n <- str_sub(n,end = -4)
nDereg <- str_sub(nDereg,end = -2)

str_which(R_NICO100Table,nDereg)

str_detect()


library(readr)
NASIR_FINAL_TABLE_COUNTS <- read_csv("NASIR_FINAL_TABLE_COUNTS.csv")
NASIR_FINAL_TABLE_COUNTS <-  t(NASIR_FINAL_TABLE_COUNTS)
dim(NASIR_FINAL_TABLE_COUNTS)

head(NASIR_FINAL_TABLE_COUNTS)
rownames(NASIR_FINAL_TABLE_COUNTS) <- str_sub(rownames(NASIR_FINAL_TABLE_COUNTS),end=15)
NASIR_FINAL_TABLE_COUNTS[,43] <- 0

mutate(NASIR_FINAL_TABLE_COUNTS,rn)


left_join(nDereg,NASIR_FINAL_TABLE_COUNTS,by=c())
