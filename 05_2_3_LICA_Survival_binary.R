#           13.7.1. Survival Clustering binary - LICA
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Survival variable is split into two different groups, binary,
#so, the model that will be trained in order to predict survival will not predict 
#exact values, but ranges. This is done so as to simplify the model.

library(readxl)

setwd("LICA_genetype")
survival <- read_excel("survival.xlsx", na = "NA")

patient_id <- survival[,1]

survival <- na.omit(survival)
surv <- survival[,2]
cuantiles<-quantile(surv, probs = c(.5,1),na.rm=TRUE)


G1 <- survival%>%
  filter(`Last survival (days)`< cuantiles["50%"])
G1["Group"] <- "G1"
G2 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["100%"]) %>%
  filter(`Last survival (days)`> cuantiles["50%"])
G2["Group"] <- "G2"


survival <- bind_rows(G1,G2)

survival_origin <- read_excel("survival.xlsx", na = "NA")
survival_binary <- full_join(survival_origin,survival)
survival_binary <- survival_binary%>%
  select(-c(`Last survival (days)`))

write.csv(survival_binary,"survival_binary.csv")

