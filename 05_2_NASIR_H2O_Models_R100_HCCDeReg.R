#                         06. R100 H2O Modeling - NASIR
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Seen that model development with ~58.000 variables, makes no sense,
#now the dataset will be first filtered and then with a smaller cohort, the same
#procedure as in: "04.NASIR_H2O_Models" will be developed.

library(tidyverse)
library(readxl)

setwd("")

R_NICO100Table <- read_excel("R-NICO100Table.xlsx")
HCCDeReg <- read_excel("HCCDeReg.xlsx")
nDereg <- HCCDeReg[,1]

#Manipulating the gene_id:
NASIR_FINAL_TABLE_COUNTS <- read_excel("CUENTAS_TUMORES.xlsx")
dim(NASIR_FINAL_TABLE_COUNTS)
num <- str_sub(NASIR_FINAL_TABLE_COUNTS$GENE,start=6,end = 15)
numerico <- as.integer(num)
NASIR_FINAL_TABLE_COUNTS$NUM <- numerico

#Extract gene_id for the curated R100 dataset:
num100 <- str_sub(R_NICO100Table$ENSEMBLEID,start=6,end = 15)
numerico100 <- as.integer(num100)
R_NICO100Table$NUM <- numerico100

#Matching both tables:

match <- inner_join(NASIR_FINAL_TABLE_COUNTS,R_NICO100Table,by=("NUM"))



match_filt <- t(match)
colnames(match_filt) <- match_filt[1,]
match_filt <- match_filt[2:31,]
match_filt <- as.data.frame(match_filt)
summary(match_filt)
match_filt <- as.double(match_filt)
match_filt <- as.data.frame(sapply(match_filt,as.numeric))
str(match_filt)

summary(match_filt)                            

# Applying filters to the dataset:
i <- 1
for(i in 1:95){
match_filt[31,i] <- min(match_filt[1:30,i]) #minimo expresion.
match_filt[32,i] <- max(match_filt[1:30,i]) #max expresion.
match_filt[33,i] <- (mean(match_filt[1:30,i])==0) #la expresion es 0 en todos los pacientes.
match_filt[34,i] <- max(match_filt[1:30,i])>50 #el valor maximo es mayor que 50.
match_filt[35,i] <- (sum(match_filt[1:30,i]>10))>5 # expresion mayor que 10 en mas de 5 pacientes.
match_filt[36,i] <- (max(match_filt[1:30,i])/min(match_filt[1:30,i]))>10 #el log_fold es mayor que 10. 
i <- i+1
}

t_resumen <-as.data.frame(c("Exp. =0 en todos los pacientes","Max >50", ">10 en > 5 pacientes", "Fold Change > 10"))
colnames(t_resumen) <- "Filtro"
t_resumen[1,2] <- sum(match_filt[33,])
t_resumen[2,2] <- sum(match_filt[34,])
t_resumen[3,2] <- sum(match_filt[35,])
t_resumen[4,2] <- sum(match_filt[36,],na.rm=TRUE)
colnames(t_resumen) <- c("Filtro","Nº filtrados")


i <- 1
j <- 1
to_remove <- as.data.frame(1)

for(i in 1:95){
  if (match_filt[33,i]==1){
    to_remove[j,1] <- i
    j <- j+1}
}


################################################################################
#NASIR vs. HCCDeReg:
#Same procedure, but now using the other curated dataset: "HCCDeReg".
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

i <- 1
for(i in 1:861){
  match_filt_de_reg[31,i] <- min(match_filt_de_reg[1:30,i]) #minimo expresion.
  match_filt_de_reg[32,i] <- max(match_filt_de_reg[1:30,i]) #max expresion.
  match_filt_de_reg[33,i] <- (mean(match_filt_de_reg[1:30,i])==0) #la expresion es 0 en todos los pacientes.
  match_filt_de_reg[34,i] <- max(match_filt_de_reg[1:30,i])>50 #el valor maximo es mayor que 50.
  match_filt_de_reg[35,i] <- (sum(match_filt_de_reg[1:30,i]>10))>5 # expresion mayor que 10 en mas de 5 pacientes.
  match_filt_de_reg[36,i] <- (max(match_filt_de_reg[1:30,i])/min(match_filt_de_reg[1:30,i]))>10 #el log_fold es mayor que 10. 
  i <- i+1
}

t_resumen_de_reg <-as.data.frame(c("Exp. =0 en todos los pacientes","Max >50", ">10 en > 5 pacientes", "Fold Change > 10"))
colnames(t_resumen_de_reg) <- "Filtro"
t_resumen_de_reg[1,2] <- sum(match_filt_de_reg[33,])
t_resumen_de_reg[2,2] <- sum(match_filt_de_reg[34,])
t_resumen_de_reg[3,2] <- sum(match_filt_de_reg[35,])
t_resumen_de_reg[4,2] <- sum(match_filt_de_reg[36,],na.rm=TRUE)
colnames(t_resumen_de_reg) <- c("Filtro","Nº filtrados")


i <- 1
j <- 1
to_remove_de_reg <- as.data.frame(1)

for(i in 1:861){
  if (match_filt_de_reg[33,i]==1){
    to_remove_de_reg[j,1] <- i
    j <- j+1}
}

################################################################################

numero_paciente <- as.double(colnames(NASIR_FINAL_TABLE_COUNTS)[2:31])

Clinical_data <- read_excel("NASIR_Main_Table_001_survival.xlsx")
OS <- Clinical_data[,c(3,12)]
colnames(OS) <- c("NUM", "OS")

################################################################################
#Model development:H2O

match_filt[1:30,96] <- numero_paciente
match_filt <- as.data.frame(sapply(match_filt,as.numeric))

colnames(match_filt)[96] <- "NUM"
prueba <- inner_join(match_filt,OS,by=("NUM"))
prueba <- prueba[c(1:30),]
prueba <- prueba%>%
  select(-NUM)

Data <- prueba

library(tidyverse)
library(lubridate)
library(tidymodels)
library(skimr)
library(DataExplorer)
library(ggpubr)
library(mosaicData)
library(h2o)
library(readxl)

# Split data into train and test:
# ==============================================================================
set.seed(123)
split_inicial <- initial_split(
  data   = Data,
  prop   = 0.8,
  strata = OS
)

datos_train <- training(split_inicial)
datos_test  <- testing(split_inicial)

summary(datos_train$OS)
summary(datos_test$OS)


# Se almacenan en un objeto `recipe` todos los pasos de preprocesado y, finalmente,
# se aplican a los datos.
transformer <- recipe(
  formula = OS ~ .,
  data =  datos_train
) %>%
  step_naomit(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_center(all_numeric(), -all_outcomes()) %>%
  step_scale(all_numeric(), -all_outcomes()) %>%
  step_dummy(all_nominal(), -all_outcomes())

transformer


# Se entrena el objeto recipe
transformer_fit <- prep(transformer)

# Se aplican las transformaciones al conjunto de entrenamiento y de test
datos_train_prep <- bake(transformer_fit, new_data = datos_train)
datos_test_prep  <- bake(transformer_fit, new_data = datos_test)

glimpse(datos_train_prep)


# Cluster initialization:
# ==============================================================================
h2o.init(
  nthreads = -1,
  max_mem_size = "4g"
)


# Se eliminan los datos del cluster por si ya habÃ­a sido iniciado.
h2o.removeAll()
h2o.no_progress()

datos_train  <- as.h2o(datos_train_prep, key = "datos_train")
datos_test   <- as.h2o(datos_test_prep, key = "datos_test")

# Defining hyperparameters:
# ==============================================================================
hiperparametros <- list(
  epochs = c(50, 100, 500, 1000, 2000),
  hidden = list(5, 10, 25, 50,100, c(10, 10))
)


# Search of best models by cross-validation:
# ==============================================================================
variable_respuesta <- 'OS'
predictores <- setdiff(colnames(datos_train), variable_respuesta)

grid <- h2o.grid(
  algorithm    = "deeplearning",
  activation   = "Rectifier",
  x            = predictores,
  y            = variable_respuesta,
  training_frame  = datos_train,
  nfolds       = 20, #validacion cruzada
  standardize  = FALSE,
  hyper_params = hiperparametros,
  search_criteria = list(strategy = "Cartesian"),
  seed         = 123,
  grid_id      = "grid"
)

resultados_grid <- h2o.getGrid(
  sort_by = 'rmse',
  grid_id = "grid",
  decreasing = FALSE
)
data.frame(resultados_grid@summary_table)

# Best model found:
# ==============================================================================
modelo_final <- h2o.getModel(resultados_grid@model_ids[[1]])

predicciones <- h2o.predict(
  object  = modelo_final,
  newdata = datos_test
)

predicciones <- predicciones %>%
  as_tibble() %>%
  mutate(valor_real = as.vector(datos_test$M))

predicciones %>% head(5)

modelo_final@allparameters

mojo_destination <- h2o.save_mojo(modelo_final, path = getwd(),filename = "modelo_M")

# Before finishing the script, the H2O cluster must be closed:
# ==============================================================================
h2o.shutdown(prompt = FALSE)
