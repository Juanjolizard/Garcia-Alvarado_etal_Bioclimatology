setwd('C:/BIOCLIMATOLOGIA_TENERIFE/SCRIPTS/')

####MLR###

library(raster)
library(terra)
library(leaps)
library(dplyr)
library(MASS)
library(metrica)
library(corrplot)
library(automap)
library(gstat)
library(sf)
library(sp)
library(rgdal)
library(caret)
library(ggplot2)
library(stringr)
library(cli)
library(fields)
library(spm2)

df <- read.csv('../ESTACIONES/estaciones_io.csv', header = T)

# df1 <- read.csv('../DATOS/estaciones_temperatura.csv', header = T)

####dst folder###
dst_folder <- 'C:/BIOCLIMATOLOGIA_TENERIFE/TPS/'


###covariable##
ALTITUD <- raster('../RASTERS_ENTRADA/DEM_TENERIFE.tif')


####tps para Pp

tps<- Tps(df[,3:5], df$estacion10, lambda = 0)



p3 <- interpolate(ALTITUD, tps, xyOnly = F)


plot(p3)


writeRaster(p3, 'C:/BIOCLIMATOLOGIA_TENERIFE/PCA/Pp.tif', format = 'Gtiff', overwrite = T)

###tps1 con optimización de parámetros###

tps.parameters <- fastTpsMLE(df[,3:5], df$estacion10, theta = 3)

lam <- tps.parameters$lambda.fixed

tps1 <- Tps(df[,3:5], df$estacion10, lambda = 0)


p3 <- interpolate(ALTITUD, tps1, xyOnly = F)


plot(p3)


###croosvalidation 100 veces, para poder predictivo##
n <- 100
set.seed(1234)
tpsvecv2 <- NULL

for(i in 1:n){
  
  tpscv1 <- tpscv(df[,3:5], df$estacion10, validation = 'CV', cv.fold = 10, predacc = "ALL")
  tpsvecv2[i] <- tpscv1

  
}


for(i in 1:n){
  
  tpscv1 <- tpscv(df[,3:5], df$estacion10, validation = 'LOO', predacc = "ALL")
  tpsvecv2[i] <- tpscv1
  
  
}

tpsvecv2





mean(tpsvecv2)

mae <- as.data.frame(tpsvecv2)


ve <- as.data.frame(tpsvecv2)
mae <- as.data.frame(tpsvecv2)



df_validacion <- as.data.frame()

mean(tpsvecv2) ###0.83
range(tpsvecv2)

names(ALTITUD) <- 'ALTITUD'


# ###tps 
# tps2 <- Tps(df[,c(3,4,5)], df$estacion10)
# tps2
# 
# p3 <- interpolate(ALTITUD, tps2, xyOnly = F)
# 
# plot(p3)
# 
# df_pred <- predict(tps2, df[,c(3,4,5)])
# 
# df_pred_final <- data.frame(df$estacion10, df_pred)
# 
# metrica::R2(data = df_pred_final, obs = df_pred_final$df.estacion10, pred = df_pred_final$df_pred)
# 
# writeRaster(p3, 'C:/BIOCLIMATOLOGIA_TENERIFE/RASTERS_GENERADOS/Io_TPS.tif', format = 'Gtiff')


#####create test y training datasets###

columns <- c('Réplica','R2', 'RMSE', 'MAE')

length(columns)

df_train_plot <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df_train_plot) <- c('Réplica','R2','RMSE','MAE')

columns_test <- c('Réplica','R2','RMSE','MAE')

df_test_plot <- data.frame(matrix(nrow = 0, ncol = length(columns_test)))
colnames(df_test_plot) <- c('Réplica','R2','RMSE','MAE')

tecnica <- 'Thin Plate Splines using elevation as covariable'

# ####empezar interpolación###
# for (i in seq(1:nrow(df))){
#   n <- 10 + i
#   set.seed(n)
#   df_train<- df[-i,]
#   df_test <- df[i,]
#   
#   tps <- Tps(df_train[,c(3,4,5)], df_train$estacion10)
#   
#   
#   ###hacer la interpolación##
#   
#   p3 <- interpolate(ALTITUD, tps, xyOnly = F)
#   
#   plot(p3, main = paste0('Interpolación iteración ', i, '  ', 'Tenerife', ' ', tecnica))
#   
#   ###exportar el raster de la interpolación###
#   
#   # writeRaster(p3, filename = paste0(dst_folder,'TPS_LOOCV','_',i,'.tif'), format = 'Gtiff', overwrite = T)
#   
#   
#   df_pred <- predict(tps, df_train[,c(3,4,5)])
#   
#   df_pred_final_train <- data.frame(df_train$estacion10, df_pred)
#   
#   r2_train <- metrica::R2(data = df_pred_final_train, obs = df_pred_final_train$df_train.estacion10, pred = df_pred_final_train$df_pred)
#   rmse_train <- metrica::RMSE(data = df_pred_final_train, obs = df_pred_final_train$df_train.estacion10, pred = df_pred_final_train$df_pred)
#   mae_train <- metrica::MAE(data = df_pred_final_train, obs = df_pred_final_train$df_train.estacion10, pred = df_pred_final_train$df_pred)
#   df_train_plot[nrow(df_train_plot) + 1,] = c(i,r2_train, rmse_train, mae_train)
#   
#   ###fase de testing##
#   
#   df_pred_test <- predict(tps, df_test[,c(3,4,5)])
#   
#   # p3 <- interpolate(r_stack, tps, xyOnly = F, fun = pfun)
#   ####también se puede cambiar por un predict con la estación correspondiente###
#   
#   df_pred_final_test <- data.frame(df_test$estacion10, df_pred_test)
#   
#   r2_test <- metrica::R2(data = df_pred_final_test, obs = df_pred_final_test$df_test.estacion10, pred = df_pred_final_test$df_pred_test)
#   rmse_test <- metrica::RMSE(data = df_pred_final_test, obs = df_pred_final_test$df_test.estacion10, pred = df_pred_final_test$df_pred_test)
#   mae_test <- metrica::MAE(data = df_pred_final_test, obs = df_pred_final_test$df_test.estacion10, pred = df_pred_final_test$df_pred_test)
#   df_test_plot[nrow(df_test_plot) + 1,] = c(i,r2_test, rmse_test, mae_test, df_pred_test)
# }
# 
# df_validacion_loocv <- data.frame(df$estacion10, df_test_plot$Pred)
# r2 <- metrica::R2(data = df_validacion_loocv, obs = estacion10, pred = df_test_plot.Pred)
# 
# ###JLME##
# 
# fit<-lm(df_validacion_loocv$df.estacion10 ~ df_validacion_loocv$df_test_plot.Pred, data=df_validacion_loocv)
# R2<-format(summary(fit)$adj.r.squared, digits=2)
# 
# 
# tiles_plot(data = df_validacion_loocv, obs = df.estacion10, pred = df_test_plot.Pred,
#            bins = 10,
#            colors = c(low = 'yellow', high = 'green'))


###bucle pero con particiones 70:30###

for (i in seq(1:10)){
  n <- 10 + i
  set.seed(n)
  df_train<- sample_frac(df, 0.7)
  df_test <- subset(df, !(COD %in% c(df_train$COD)))
  
  tps <- Tps(df_train[,c(3,4,5)], df_train$estacion10)
  
  
  ###hacer la interpolación##
  
  p3 <- interpolate(ALTITUD, tps, xyOnly = F)
  
  plot(p3, main = paste0('Interpolación iteración ', i, '  ', 'Tenerife', ' ', tecnica))
  
  ###exportar el raster de la interpolación###
  
  writeRaster(p3, filename = paste0(dst_folder,'TPS_IO_TF','_',i,'.tif'), format = 'Gtiff', overwrite = T)
  
  
  df_pred <- predict(tps, df_train[,c(3,4,5)])
  
  df_pred_final_train <- data.frame(df_train$estacion10, df_pred)
  
  r2_train <- metrica::R2(data = df_pred_final_train, obs = df_pred_final_train$df_train.estacion10, pred = df_pred_final_train$df_pred)
  rmse_train <- metrica::RMSE(data = df_pred_final_train, obs = df_pred_final_train$df_train.estacion10, pred = df_pred_final_train$df_pred)
  mae_train <- metrica::MAE(data = df_pred_final_train, obs = df_pred_final_train$df_train.estacion10, pred = df_pred_final_train$df_pred)
  df_train_plot[nrow(df_train_plot) + 1,] = c(i,r2_train, rmse_train, mae_train)
  
  ###fase de testing##
  
  df_pred_test <- predict(tps, df_test[,c(3,4,5)])
  
  # p3 <- interpolate(r_stack, tps, xyOnly = F, fun = pfun)
  ####también se puede cambiar por un predict con la estación correspondiente###
  
  df_pred_final_test <- data.frame(df_test$estacion10, df_pred_test)
  
  r2_test <- metrica::R2(data = df_pred_final_test, obs = df_pred_final_test$df_test.estacion10, pred = df_pred_final_test$df_pred_test)
  rmse_test <- metrica::RMSE(data = df_pred_final_test, obs = df_pred_final_test$df_test.estacion10, pred = df_pred_final_test$df_pred_test)
  mae_test <- metrica::MAE(data = df_pred_final_test, obs = df_pred_final_test$df_test.estacion10, pred = df_pred_final_test$df_pred_test)
  df_test_plot[nrow(df_test_plot) + 1,] = c(i,r2_test, rmse_test, mae_test)
}

# ###JLME##
# 
# fit<-lm(df_validacion_loocv$df.estacion10 ~ df_validacion_loocv$df_test_plot.Pred, data=df_validacion_loocv)
# R2<-format(summary(fit)$adj.r.squared, digits=2)
# 
# 
# tiles_plot(data = df_validacion_loocv, obs = df.estacion10, pred = df_test_plot.Pred,
#            bins = 10,
#            colors = c(low = 'yellow', high = 'green'))


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


####promedio de los rasters###
setwd(dst_folder)
lista_rasters <- list.files(pattern = '*.tif')

raster_stack <- raster::stack(lista_rasters)
raster_promedio <- mean(raster_stack)
plot(raster_promedio)
writeRaster(raster_promedio, filename = 'promedio.tif', format = 'Gtiff')