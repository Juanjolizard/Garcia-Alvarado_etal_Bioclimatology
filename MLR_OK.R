####script for Itc and Tp interpolations##
###using MLR + OK of the residuals###
###Garcia-Alvarado et al. 2023##

####required packages###

if(!require(raster)) install.packages('raster')
if(!require(terra)) install.packages('terra')
if(!require(cli)) install.packages('cli')
if(!require(stringr)) install.packages('stringr')
if(!require(dplyr)) install.packages('dplyr')
if(!require(metrica)) install.packages('metrica')
if(!require(MASS)) install.packages('MASS')
if(!require(ggplot2)) install.packages('ggplot2')
if(!require(rgdal)) install.packages('rgdal')
if(!require(caret)) install.packages('caret')
if(!require(gstat)) install.packages('gstat')
if(!require(automap)) install.packages('automap')
if(!require(sf)) install.packages('sf')
if(!require(sp)) install.packages('sp')
if(!require(leaps)) install.packages('leaps')


###Packages##
library(raster)
library(terra)
library(leaps)
library(dplyr)
library(MASS)
library(metrica)
library(automap)
library(gstat)
library(sf)
library(sp)
library(rgdal)
library(caret)
library(ggplot2)
library(stringr)
library(cli)


setwd('C:/BIOCLIMATOLOGIA_TENERIFE')


###create a destiny folder###

unidad <- 'C:/' ###disk##

directorio <- 'BIOCLIMATOLOGIA_TENERIFE/MLR_OK_VALIDACION/' ###directory to save rasters and validation results###


dst_folder <- paste0(unidad, directorio)


if(!dir.exists(dst_folder)){
  dir.create(dst_folder)
  print("Path OK")
} else{
  print("Path already exists")
}

###dataset###

df <- read.csv('DATOS/estaciones_temperatura.csv', header = T)



###covariables###
YUTM <- raster('RASTERS_ENTRADA/YUTM.tif')
Elevation <- raster('RASTERS_ENTRADA/DEM_TENERIFE.tif')
names(Elevation) <- 'Elev1'
XUTM <- raster('RASTERS_ENTRADA/XUTM.tif')
names(XUTM) <- 'XUTM.'
names(YUTM) <- 'YUTM'

r <- raster::stack(YUTM, Elevation, XUTM)


####Itc model###
train.control <- trainControl(method = "LOOCV")


step.model1 <- train(Itc ~ XUTM + YUTM + Elev1 , data = df,
                     method = "lmStepAIC",
                     trControl = train.control,
                     trace = FALSE
)


# Model accuracy
step.model1$results
# Final model coefficients
step.model1$finalModel
# Summary of the model
summary(step.model1$finalModel)


####Spatial model###
prediction_model1 <- predict(r, step.model1$finalModel)

plot(prediction_model1)


###Ordinary Kriging Residuals interpolation###
####ordinary kringing de los residuos

df_validacion <- extract(prediction_model1, df[2:3])

df_final_val <- cbind(df, df_validacion)

df_final_val$Residuos <- df_final_val$Itc - df_final_val$df_validacion

coordinates(df_final_val) <- ~ XUTM + YUTM

proj4string(df_final_val) <- CRS("+proj=utm +zone=28 +datum=WGS84 +units=m +no_defs")

crs_stations <- CRS(SRS_string = "EPSG:32628")
slot(df_final_val, "proj4string") <- crs_stations
g <- automap::autofitVariogram(formula = Residuos ~ 1, input_data = df_final_val)

plot(g)
g <- gstat::gstat(formula = Residuos ~ 1, data = df_final_val, model = g$var_model)

z1 <- interpolate(prediction_model1, g) ###el prediction_model1 sólo se usa como grid

plot(z1)

####MLR + OK

mlr_ok <- z1 + prediction_model1

superficies <- raster::brick(prediction_model1, mlr_ok)

plot(superficies)

writeRaster(mlr_ok, filename = paste0(dst_folder, 'Itc.tif'), format = 'Gtiff', overwrite = T)



#########VALIDATION####################
#### Manual LOOCV    ###

###test dataset to compare predicted vs observed values###

columns <- c('Réplica','Obs', 'Pred')
df_test_plot <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df_test_plot) <- c('Réplica','Obs','Pred')


tecnica <- "Multiple Linear Regression + Ordinary Kriging of the residuals"


for (i in seq(1:nrow(df))){
  n <- 10 + i
  set.seed(n)
  df_train<- df[-i,]
  df_test <- df[i,]
  
  
  train.control <- trainControl(method = "LOOCV")
  
  
  step.model1 <- train(Itc ~ XUTM + YUTM + Elev1 , data = df_train,
                       method = "lmStepAIC",
                       trControl = train.control,
                       trace = FALSE
  )
  
  prediction_model1 <- predict(r, step.model1$finalModel)
  

  ###predichos y observados###
  pred <- step.model1$pred
  df_train_val <- cbind(df_train, pred)
  
  df_train_val$Residuos <- df_train_val$Itc - df_train_val$pred
  
  coordinates(df_train_val) <- ~ XUTM + YUTM
  
  proj4string(df_train_val) <- CRS("+proj=utm +zone=28 +datum=WGS84 +units=m +no_defs")
  
  crs_stations <- CRS(SRS_string = "EPSG:32628")
  slot(df_train_val, "proj4string") <- crs_stations
  g <- automap::autofitVariogram(formula = Residuos ~ 1, input_data = df_train_val)
  
  plot(g)
  g <- gstat::gstat(formula = Residuos ~ 1, data = df_train_val, model = g$var_model)
  
  z1 <- interpolate(prediction_model1, g) ###grid###
  
  plot(z1)
  
  ####MLR + OK
  
  mlr_ok <- z1 + prediction_model1
  
  superficies <- raster::brick(prediction_model1, mlr_ok)
  
  plot(superficies)
  
  
  ###test###
  itc_interpolado_test <- raster::extract(mlr_ok, df_test[2:3])
  df_test_plot[nrow(df_test_plot) + 1,] = c(i,df_test$Itc, itc_interpolado_test)
  print(paste0('iteration ', i, '/', nrow(df),' : ', tecnica))
  
}

####Validation metrics and plots###
###global metrics###

r2 <- metrica::R2(df_test_plot, df_test_plot$Obs, df_test_plot$Pred)
mae <- metrica::MAE(df_test_plot, df_test_plot$Obs, df_test_plot$Pred)
rmse <- metrica::RMSE(df_test_plot, df_test_plot$Obs, df_test_plot$Pred)

###obs vs pred plot###

tiles.plot <- 
  tiles_plot(data = df_test_plot, 
             obs = Obs, 
             pred = Pred,
             bins = 10, 
             orientation = "PO",
             colors = c(low = "green", high = "yellow"))

tiles.plot















