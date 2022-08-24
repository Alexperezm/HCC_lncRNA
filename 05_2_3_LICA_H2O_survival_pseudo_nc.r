#           13.6. H2O DL Survival prediction model - LICA
#       Author: Alexperezm | Master's End of Degree Project - 2021-2022

#Objective: Development of a Deep Learning model, that based on the expression of
#pesudo and non-coding gene set of LICA, will pretend to the survival of that 
#patient (in days), taking into account that the survival rates have been
#clusterized previously, based on deciles and terciles.


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
LICA_normalized_filtered_genetype <- read_excel("LICA_normalized_filtered_genetype.xlsx")

ncRNA <- filter(LICA_normalized_filtered_genetype, GENETYPE == "ncRNA"| GENETYPE=="pseudo")

cnames <- (ncRNA$gene_id)

ncRNA <- select(ncRNA,-c(SYMBOL,GENENAME,GENETYPE))
ncRNA <- select(ncRNA,-gene_id)
ncRNA <- as.tibble(t(ncRNA))

#survival <- read_excel("survival_decil.xlsx", na = "NA")
survival <- read_excel("survival_tercil.xlsx", na = "NA")
survival <- survival[,2]
survival <- as.tibble(survival)


ncRNA[,(dim(ncRNA)[2]+1)] <- survival
colnames(ncRNA) <- cnames
colnames(ncRNA)[dim(ncRNA)[2]] <- "survival"
Data <- ncRNA

# Split data into train and test:
# ==============================================================================
#set.seed(123)
split_inicial <- initial_split(
  data   = Data,
  prop   = 0.8,
  strata = survival
)

datos_train <- training(split_inicial)
datos_train$survival <-  as.factor(datos_train$survival)
datos_test  <- testing(split_inicial)
datos_test$survival <-  as.factor(datos_test$survival)

summary(datos_train$survival)
summary(datos_test$survival)

#Raw data must be processed to remove NA, dummy variables, etc. A recipe has been
#developed.
transformer <- recipe(
  formula = survival ~ .,
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
  hidden = list(5, 10, 25, 50,100, c(10, 10)) # 200.
)


# Search of best models by cross-validation:
# ==============================================================================
variable_respuesta <- 'survival'
predictores <- setdiff(colnames(datos_train), variable_respuesta)

grid <- h2o.grid(
  algorithm    = "deeplearning",
  activation   = "Rectifier",
  x            = predictores,
  y            = variable_respuesta,
  training_frame  = datos_train,
  nfolds       = 20, #validacion cruzada #30
  standardize  = FALSE,
  hyper_params = hiperparametros,
  search_criteria = list(strategy = "Cartesian"),
  #seed         = 123,
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
  mutate(valor_real = as.vector(datos_test$survival))

predicciones %>% head(5)

modelo_final@allparameters

mojo_destination <- h2o.save_mojo(modelo_final, path = getwd(),filename = "modelo_pseudo+nc_nfold20_survival_terciles_200_30")

# Before finishing the script, the H2O cluster must be shutdown:
# ==============================================================================
h2o.shutdown(prompt = FALSE)
