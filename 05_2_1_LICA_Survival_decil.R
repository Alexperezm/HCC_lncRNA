#           13.6.1. Survival Clustering into decils - LICA
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Survival variable is split into ten different groups, into decils,
#so, the model that will be trained in order to predict survival will not predict 
#exact values, but ranges (10 groups). This is done so as to simplify the model.


setwd("")
survival <- read_excel("survival.xlsx", na = "NA")

patient_id <- survival[,1]

survival <- na.omit(survival)
surv <- survival[,2]
cuantiles<-quantile(surv, probs = c(.10,.20,.30,.40,.50,.60,.70,.80,.90,1),na.rm=TRUE)


G1 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["10%"])
G1["Group"] <- "G1"
G2 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["20%"]) %>%
  filter(`Last survival (days)`> cuantiles["10%"])
G2["Group"] <- "G2"
G3 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["30%"])%>%
  filter(`Last survival (days)`> cuantiles["20%"])
G3["Group"] <- "G3"
G4 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["40%"])%>%
  filter(`Last survival (days)`> cuantiles["30%"])
G4["Group"] <- "G4"
G5 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["50%"])%>%
  filter(`Last survival (days)`> cuantiles["40%"])
G5["Group"] <- "G5"
G6 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["60%"])%>%
  filter(`Last survival (days)`> cuantiles["50%"])
G6["Group"] <- "G6"
G7 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["70%"])%>%
  filter(`Last survival (days)`> cuantiles["60%"])
G7["Group"] <- "G7"
G8 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["80%"])%>%
  filter(`Last survival (days)`> cuantiles["70%"])
G8["Group"] <- "G8"
G9 <- survival%>%
  filter(`Last survival (days)`<= cuantiles["90%"])%>%
  filter(`Last survival (days)`> cuantiles["80%"])
G9["Group"] <- "G9"
G10<- survival%>%
  filter(`Last survival (days)`<= cuantiles["100%"])%>%
  filter(`Last survival (days)`> cuantiles["90%"])
G10["Group"] <- "G10"
survival <- bind_rows(G1,G2,G3,G4,G5,G6,G7,G8,G9,G10)

survival_origin <- read_excel("survival.xlsx", na = "NA")
survival_decil <- full_join(survival_origin,survival)
survival_decil <- survival_decil%>%
  select(-c(`Last survival (days)`))

write.csv(survival_decil,"survival_decil.csv")

