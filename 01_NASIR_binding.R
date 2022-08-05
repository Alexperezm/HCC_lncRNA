################################################################################
#########################01.Binding NASIR raw files#############################
########Author: Alexperezm | Master's End of Degree Project - 2021-2022#########
################################################################################
#Objective: This script has been developed to bind all the NASIR RNA count files
#and gather all data in two main datasets: RNA counts for Tumor (T) and
#Peritumour (P).

library(tidyverse)
library(stringr)

#1st: IS obtained a list with the files available in the direction and based on
#that, consequent joining are done in order to gather all data in one object,
#which will be saved and used all along the project.
#Raw files available on: "NASIR_htseq-20220109T070602Z-001".

setwd("")
archivos <- list.files(path=".")
i <- 1
dataframe <- read_table(file = archivos[i],col_names = FALSE)
names(dataframe) <- c("Gene",i)

i=i+1
while (i<length(archivos)+1) {
  dataframe_1 <- read_delim(archivos[i],col_names = FALSE)
  numero <- str_sub(archivos[i],start=13,end=-24)
  names(dataframe_1) <- c("Gene",numero)
  dataframe <- left_join(x = dataframe,y = dataframe_1,by=c("Gene"="Gene"))
  i <- i+1
  }
write.table(dataframe,"NASIR_P.csv")
saveRDS(dataframe,file = "NASIR_P.rds")


#Same procedure for Tumor:

setwd("")
archivos <- list.files(path=".")
i <- 1
dataframe <- read_delim(archivos[i],col_names = FALSE)
names(dataframe) <- c("Gene",i)
i=i+1

while (i<length(archivos)+1) {
  dataframe_1 <- read_delim(archivos[i],col_names = FALSE)
  numero <- str_sub(archivos[i],start=12,end=-24)
  names(dataframe_1) <- c("Gene",numero)
  dataframe <- left_join(x = dataframe,y = dataframe_1,by=c("Gene"="Gene"))
  i <- i+1
}
write.table(dataframe,"NASIR_T.csv")
saveRDS(dataframe,"NASIR_T.rds")
