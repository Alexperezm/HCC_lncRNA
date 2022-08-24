#           13.6.2. Survival Clustering into tercils - LICA
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Survival variable is split into three different groups, into tercils,
#so, the model that will be trained in order to predict survival will not predict 
#exact values, but ranges. This is done so as to simplify the model.

library(readxl)

setwd("")
survival <- read_excel("survival.xlsx", na = "NA")

patient_id <- survival[,1]

survival <- na.omit(survival)
surv <- survival[,2]
cuantiles<-quantile(surv, probs = c(.33,.66,1),na.rm=TRUE)


G1 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["33%"])
G1["Group"] <- "G1"
G2 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["66%"]) %>%
  filter(`Last survival (days)`> cuantiles["33%"])
G2["Group"] <- "G2"
G3 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["100%"])%>%
  filter(`Last survival (days)`> cuantiles["66%"])
G3["Group"] <- "G3"


survival <- bind_rows(G1,G2,G3)

survival_origin <- read_excel("survival.xlsx", na = "NA")
survival_tercil <- full_join(survival_origin,survival)
survival_tercil <- survival_tercil%>%
  select(-c(`Last survival (days)`))

write.csv(survival_tercil,"survival_tercil.csv")

