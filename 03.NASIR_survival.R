#                 03. Survival - NASIR
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Adding the variable "Overall Survival (OS)" to the normalized dataset.
#This will simplify following procedures, as OS will be key variable in several steps.

#Uploading of the normalized dataset, processing, biding of the clinical variable
#and following writing of the whole dataset: gene expression + outcome variable.

DF <- (assay(rld))

Rnames <- rownames(assay(rld))
Cnames <- colnames(assay(rld))

tDF <- t(DF)

rownames(tDF) <- Cnames

View(tDF)

TDFDF <- data.frame(tDF)

TDFDF[,58722] <- se$Type
TDFDF[,58723] <- se$Patient
colnames(TDFDF)[58722] <- "Type"
colnames(TDFDF)[58723] <- "ID"
TDFDF <- TDFDF[-c(4), ]
saveRDS(TDFDF,file="tabla_lnc_type_OS.rds")
write_csv(TDFDF,file = "tabla_lnc_type_OS.csv")


numbers <- TDFDF[,58723]
numbers <- as.numeric(numbers)
numbers <- sort(numbers)
numbers

unumbers <- unique(numbers)

library(readxl)
survival <- read_xlsx("NASIR_Main_Table_001_survival.xlsx")
OS <- data.frame(survival$NumId_3)
OS[,2] <- survival$OS_days_68
colnames(OS) <- c("ID","OS")


IDs <- data.frame(as.numeric(TDFDF$ID))
#data_frame(IDs)
#IDs[,2] <- TDFDF$ID
colnames(IDs) <- "ID"
View(TDFDF)

  
OS_f <- inner_join(IDs,OS,by="ID")
TDFDF[,58724] <- OS_f$OS
#  mutate(OS_f$OS)
colnames(TDFDF)[58724] <- "OS"


write_csv(TDFDF, "NASIR_FINAL_TABLE_COUNTS.csv")

library(writexl)

writexl::write_xlsx(TDFDF,path = "NASIR_FINAL_TABLE_COUNTS.xlsx")

################################################################################

DF <- (assay(rld))

Rnames <- rownames(assay(rld))
Cnames <- colnames(assay(rld))

tDF <- t(DF)

rownames(tDF) <- Cnames

View(tDF)

TDFDF <- data.frame(tDF)

relevantes <- rownames(resSig)

TDFDF <-  TDFDF %>%
  select(relevantes)

TDFDF[,10305] <- se$Type
TDFDF[,10306] <- se$Patient
colnames(TDFDF)[10305] <- "Type"
colnames(TDFDF)[10306] <- "ID"
saveRDS(TDFDF,file="tabla_lnc_type_OS_selected10000.rds")
write_csv(TDFDF,file = "tabla_lnc_type_OS_selected10000.csv")


numbers <- TDFDF[,10306]
numbers <- as.numeric(numbers)
numbers <- sort(numbers)
numbers

unumbers <- unique(numbers)

library(readxl)
survival <- read_xlsx("NASIR_Main_Table_001_survival.xlsx")

OS <- survival$NumId_3
OS <- data.frame(OS)
OS[,2] <- survival$OS_days_68
colnames(OS) <- c("ID","OS")
#OS[,1] <- as.integer(OS[,1])
#OS[,2] <- as.integer(OS[,2])


IDs <- TDFDF$ID
#IDs <- as.integer(IDs)
data_frame(IDs)
IDs[,2] <- TDFDF$Type

IDs <- as.integer(IDs)
OS <- as.integer(OS)
colnames(IDs) <- "ID"
View(TDFDF)
OS_f <- inner_join(IDs, OS, by="ID")
TDFDF %>%
  mutate(OS_f$OS)
colnames(TDFDF)[58724] <- "OS"


write_csv(TDFDF, "NASIR_FINAL_TABLE_COUNTS.csv")

library(writexl)

writexl::write_xlsx(TDFDF,path = "NASIR_FINAL_TABLE_COUNTS.xlsx")