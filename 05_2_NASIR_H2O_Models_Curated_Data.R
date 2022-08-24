#                  08. H2O Modeling - Curated Data - NASIR
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Seen that model development with ~58.000 variables, is not possible,
#now the same procedure will be followed, but with the curated dataset.


setwd("")
library(readxl)
library(tidyverse)

R100_DeReg_Tumor_Peritumor_REVISADO_PURI <- read_excel("R100_DeReg_Tumor&Peritumor_REVISADO_PURI.xlsx", 
                                                       sheet = "r100_Tumor_Peritumor_FILTER")
R100_DeReg_Tumor_Peritumor_REVISADO_PURI <- R100_DeReg_Tumor_Peritumor_REVISADO_PURI[,1:6]
NASIR_FINAL_TABLE_COUNTS <- read_xlsx("CUENTAS_TUMORES.xlsx")

num <- str_sub(NASIR_FINAL_TABLE_COUNTS$GENE,start=6,end = 15)
numerico <- as.integer(num)
NASIR_FINAL_TABLE_COUNTS$NUM <- numerico

target_genes <- R100_DeReg_Tumor_Peritumor_REVISADO_PURI$GENE
num <- str_sub(R100_DeReg_Tumor_Peritumor_REVISADO_PURI$GENE,start=6,end = 15)
R100_DeReg_Tumor_Peritumor_REVISADO_PURI$NUM <- as.integer(R100_DeReg_Tumor_Peritumor_REVISADO_PURI$NUM)



target_exp <- inner_join(R100_DeReg_Tumor_Peritumor_REVISADO_PURI,NASIR_FINAL_TABLE_COUNTS,by="NUM")
target_exp <- t(target_exp)
colnames(target_exp) <- target_exp[1,]

target_exp_f <- target_exp[8:37,]
target_exp_f <- as.data.frame(target_exp_f)
target_exp_f$Paciente <- as.integer(rownames(target_exp_f))


################################################################################
#Survival:
################################################################################

NASIR_Main_Table_survival <- read_excel("NASIR_Main_Table_001_survival.xlsx")
Supervivencia <- as.data.frame(NASIR_Main_Table_survival$NumId_3)
Supervivencia <- mutate(Supervivencia,NASIR_Main_Table_survival$OS_days_68)
colnames(Supervivencia) <- c("Paciente","OS")

Datos <- inner_join(target_exp_f,Supervivencia,by="Paciente")
Datos_sin_code <- subset (Datos, select = -c(Paciente))

#Turning binary the OS variable:
Datos_sin_code$OSb <- Datos_sin_code$OS>mean(Datos_sin_code$OS)
Datos_sin_code <- subset (Datos_sin_code, select = -c(OS))
Datos_sin_code$OSb <- as.integer(Datos_sin_code$OSb)
write.csv(Datos_sin_code,"Datos_clustering.csv")

################################################################################
#H2O:
################################################################################

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
  data   = Datos_sin_code,
  prop   = 0.8,
  strata = OSb
)

datos_train <- training(split_inicial)
datos_test  <- testing(split_inicial)

summary(datos_train$OS)
summary(datos_test$OS)


#Raw data must be processed to remove NA, dummy variables, etc. A recipe has been
#developed.
transformer <- recipe(
  formula = OSb ~ .,
  data =  datos_train
) %>%
  step_naomit(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_center(all_numeric(), -all_outcomes()) %>%
  step_scale(all_numeric(), -all_outcomes()) %>%
  step_dummy(all_nominal(), -all_outcomes())

transformer


# The recipe is trained:
transformer_fit <- prep(transformer)

# Transformations are developed for both: trainand test datasets:
datos_train_prep <- bake(transformer_fit, new_data = datos_train)
datos_test_prep  <- bake(transformer_fit, new_data = datos_test)

glimpse(datos_train_prep)


# Cluster initialization:
# ==============================================================================
h2o.init(
  nthreads = -1,
  max_mem_size = "4g"
)


# Remove all data already available in the cluster:
h2o.removeAll()
h2o.no_progress()

datos_train  <- as.h2o(datos_train_prep, key = "datos_train")
datos_test   <- as.h2o(datos_test_prep, key = "datos_test")

# Defining hyperparameters:
# ==============================================================================
hiperparametros <- list(
  epochs = c(50, 100, 500,1000),
  hidden = list(5, 10, 25, 50,100, c(10, 10))
)


# Search of best models by cross-validation:
# ==============================================================================
variable_respuesta <- 'OSb'
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
