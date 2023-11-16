####Script interpolation###
###Juan José García Alvarado###
###November 2023###

###set directory##

unidad <- 'C:/' ####set unit###
setwd(paste0(unidad,'BIOCLIMATOLOGIA_TENERIFE/Entrega/stations')) ###set folder with bioclimate data ###


###if required statement###
if(!require(raster)) install.packages('raster')
if(!require(terra)) install.packages('terra')
if(!require(cli)) install.packages('cli')
if(!require(leaps)) install.packages('leaps')
if(!require(stringr)) install.packages('stringr')
if(!require(dplyr)) install.packages('dplyr')
if(!require(fields)) install.packages('fields')
if(!require(metrica)) install.packages('metrica')
if(!require(spm2)) install.packages('spm2')
if(!require(MASS)) install.packages('MASS')
if(!require(ggplot2)) install.packages('ggplot2')
if(!require(rgdal)) install.packages('rgdal')


####packages###
library(raster)
library(terra)
library(leaps)
library(dplyr)
library(MASS)
library(metrica)
library(rgdal)
library(ggplot2)
library(stringr)
library(cli)
library(fields)
library(spm2)

df <- read.csv('ombrothermic_stations.csv', header = T)


####dst folder###
dst_folder <- 'C:/BIOCLIMATOLOGIA_TENERIFE/TPS/'


###auxiliary variable##
ALTITUD <- raster('../../RASTERS_ENTRADA/DEM_TENERIFE.tif')


###Interpolation whole dataset###
###default lambda###

tps1 <- Tps(df[,2:4], df$Io, lambda = 0)


tps_model <- interpolate(ALTITUD, tps1, xyOnly = F)


plot(tps_model)


####100 replicates 10-fold Cross-Validation###

n <- 100
set.seed(1234)
tpsvecv2 <- NULL

for(i in 1:n){
  
  tpscv1 <- tpscv(df[,2:4], df$Io, validation = 'CV', cv.fold = 10, predacc = "ALL")
  tpsvecv2[i] <- tpscv1
  
  
}

cv_df <-as.data.frame(tpscv1) ##data resume###



####Cross-Validation: 70% Train - 30% Test###
set.seed(123)

columns <- c('Réplica','R2', 'RMSE', 'MAE')

length(columns)

df_train_plot <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df_train_plot) <- c('Réplica','R2','RMSE','MAE')

columns_test <- c('Réplica','R2','RMSE','MAE')

df_test_plot <- data.frame(matrix(nrow = 0, ncol = length(columns_test)))
colnames(df_test_plot) <- c('Réplica','R2','RMSE','MAE')

tecnica <- 'Thin Plate Splines using elevation as covariable'

for (i in seq(1:10)){
  df_train<- sample_frac(df, 0.7)
  df_test <- subset(df, !(fid %in% c(df_train$fid))) ###not included in train ###
  
  tps <- Tps(df_train[,c(2:4)], df_train$Io)
  
  
  ###interpolation##
  
  p3 <- interpolate(ALTITUD, tps, xyOnly = F)
  
  plot(p3, main = paste0('Interpolación iteración ', i, '  ', 'Tenerife', ' ', tecnica))
  
  ###export##
  
  writeRaster(p3, filename = paste0(dst_folder,'TPS_Io_TF','_',i,'.tif'), format = 'Gtiff', overwrite = T)
  
  
  df_pred <- predict(tps, df_train[,c(2:4)])
  
  df_pred_final_train <- data.frame(df_train$Io, df_pred)
  
  r2_train <- metrica::R2(data = df_pred_final_train, obs = df_pred_final_train$df_train.Io, pred = df_pred_final_train$df_pred)
  rmse_train <- metrica::RMSE(data = df_pred_final_train, obs = df_pred_final_train$df_train.Io, pred = df_pred_final_train$df_pred)
  mae_train <- metrica::MAE(data = df_pred_final_train, obs = df_pred_final_train$df_train.Io, pred = df_pred_final_train$df_pred)
  df_train_plot[nrow(df_train_plot) + 1,] = c(i,r2_train, rmse_train, mae_train)
  
  ###testing##
  
  df_pred_test <- predict(tps, df_test[,c(2:4)])
  
  
  
  df_pred_final_test <- data.frame(df_test$Io, df_pred_test)
  
  r2_test <- metrica::R2(data = df_pred_final_test, obs = df_pred_final_test$df_test.Io, pred = df_pred_final_test$df_pred_test)
  rmse_test <- metrica::RMSE(data = df_pred_final_test, obs = df_pred_final_test$df_test.Io, pred = df_pred_final_test$df_pred_test)
  mae_test <- metrica::MAE(data = df_pred_final_test, obs = df_pred_final_test$df_test.Io, pred = df_pred_final_test$df_pred_test)
  df_test_plot[nrow(df_test_plot) + 1,] = c(i,r2_test, rmse_test, mae_test)
}


write.csv(df_train_plot, 'C:/BIOCLIMATOLOGIA_TENERIFE/validacion/train.csv')


df_test_promedio <- as.data.frame(colMeans(df_test_plot[,c(2:4)]))
colnames(df_test_promedio) <- 'Averaged Test Metrics'
###bloxplot

df_test_plot$Modo <- 'Test'
df_train_plot$Modo <- 'Train'

df_final <- rbind(df_train_plot, df_test_plot)


df_boxplot_pivot <- tidyr::pivot_longer(df_final, !Modo, names_to = 'Métrica', values_to = 'Valor')


ggplot(df_boxplot_pivot, aes(x=Métrica, y = Valor, fill = Métrica))+
  geom_boxplot()+
  facet_wrap(~ Modo)







