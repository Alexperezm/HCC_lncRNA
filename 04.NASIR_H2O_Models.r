#                 04. H2O - DL Modeling - NASIR
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Key aspect of the project, first approach to model synthesis; based 
#on the expression data, by using H2O algorithms for Deep Learning model developing.
#Predictor variables: all gene expresion. Response variable: Overall Survival.

#Library loading:
library(tidyverse)
library(lubridate)
library(tidymodels)
library(skimr)
library(DataExplorer)
library(ggpubr)
library(mosaicData)
library(h2o)
library(readxl)

setwd("")

# Data processing:
#########################################################
Data <- read_csv("NASIR_FINAL_TABLE_COUNTS.csv")
Datos_mitad1 <- Data[,1:(dim(Data)[2]/5)]
Datos_mitad1[,11743:11745] <- Data[,58722:58724]

Datos_mitad1 <- Datos_mitad1[,-c(dim(Datos_mitad1)-1)]
Data <- Datos_mitad1
#########################################################
Data <- read_csv("NASIR_FINAL_TABLE_COUNTS.csv")
Datos_mitad1 <- Data[,1:(dim(Data)[2]/20)]
dim(Datos_mitad1)
Datos_mitad1[,2937:2939] <- Data[,58722:58724]

Datos_mitad1 <- Datos_mitad1[,-c(dim(Datos_mitad1)-1)]
Data <- Datos_mitad1
##############################################################


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


#Raw data must be processed to remove NA, dummy variables, etc. A recipe has been
#developed.

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
  epochs = c(50, 100, 500, 1000),
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
  nfolds       = 10, #validacion cruzada
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