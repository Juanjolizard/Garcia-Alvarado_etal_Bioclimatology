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


df <- read.csv('../ESTACIONES/estaciones_temp_2006.csv', header = T)

# df1 <- read.csv('../DATOS/estaciones_temperatura.csv', header = T)

####dst folder###
dst_folder <- 'C:/BIOCLIMATOLOGIA_TENERIFE/MLR_OK_VALIDACION/'



#####model1###

model1 <- lm(Itc ~ X.UTM. + Y.UTM. + Height, data = df )

summary(model1)

train.control <- trainControl(method = "LOOCV")


step.model1 <- train(Itc ~ X.UTM. + Y.UTM. + Height , data = df,
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

step.model1$pred

###RASTERS###
Y.UTM. <- raster('../RASTERS_ENTRADA/YUTM.tif')
Height <- raster('../RASTERS_ENTRADA/DEM_TENERIFE.tif')
names(Height) <- 'Height'
X.UTM. <- raster('../RASTERS_ENTRADA/XUTM.tif')
names(X.UTM.) <- 'X.UTM.'
names(Y.UTM.) <- 'Y.UTM.'

r <- raster::stack(Y.UTM., Height, X.UTM.)

prediction_model1 <- predict(r, step.model1)

plot(prediction_model1)

####ordinary kringing de los residuos

df_validacion <- extract(prediction_model1, df[4:5])

df_final_val <- cbind(df, df_validacion)

df_final_val$Residuos <- df_final_val$Itc - df_final_val$df_validacion

coordinates(df_final_val) <- ~ X.UTM. + Y.UTM.

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

writeRaster(mlr_ok, filename = '../SUPERFICIES_INTERPOLADAS/itc_interpolado.tif', format = 'Gtiff', overwrite = T)

io_interpolado <- raster::extract(mlr_ok, df[3:4])

df_validacion_2 <- data.frame(df$estacion10, io_interpolado)

metrica::R2(df_validacion_2, obs = df_validacion_2$df.estacion10, pred = df_validacion_2$io_interpolado)
metrica::RMSE(df_validacion_2, obs = df_validacion_2$df.estacion10, pred = df_validacion_2$io_interpolado)
metrica::MAE(df_validacion_2, obs = df_validacion_2$df.estacion10, pred = df_validacion_2$io_interpolado)

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


#####test y training###

columns <- c('Réplica','R2', 'RMSE', 'MAE')

length(columns)

df_train_plot <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df_train_plot) <- c('Réplica','R2','RMSE','MAE')

columns_test <- c('Réplica','R2','RMSE','MAE')

df_test_plot <- data.frame(matrix(nrow = 0, ncol = length(columns_test)))
colnames(df_test_plot) <- c('Réplica','R2','RMSE','MAE')

tecnica <- 'Thin Plate Splines using elevation as covariable'

####empezar interpolación###
for (i in seq(1:nrow(df))){
  n <- 10 + i
  set.seed(n)
  df_train<- df[-i,]
  df_test <- df[i,]
  
  tps <- Tps(df_train[,c(3,4,5)], df_train$estacion10)
  
  
  ###hacer la interpolación##
  
  p3 <- interpolate(ALTITUD, tps, xyOnly = F)
  
  plot(p3, main = paste0('Interpolación iteración ', i, '  ', 'Tenerife', ' ', tecnica))
  
  ###exportar el raster de la interpolación###
  
  # writeRaster(p3, filename = paste0(dst_folder,'TPS_LOOCV','_',i,'.tif'), format = 'Gtiff', overwrite = T)
  
  
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
  df_test_plot[nrow(df_test_plot) + 1,] = c(i,r2_test, rmse_test, mae_test, df_pred_test)
}

df_validacion_loocv <- data.frame(df$estacion10, df_test_plot$Pred)
r2 <- metrica::R2(data = df_validacion_loocv, obs = estacion10, pred = df_test_plot.Pred)

###JLME##

fit<-lm(df_validacion_loocv$df.estacion10 ~ df_validacion_loocv$df_test_plot.Pred, data=df_validacion_loocv)
R2<-format(summary(fit)$adj.r.squared, digits=2)


tiles_plot(data = df_validacion_loocv, obs = df.estacion10, pred = df_test_plot.Pred,
           bins = 10,
           colors = c(low = 'yellow', high = 'green'))


###bucle pero con particiones 70:30###

for (i in seq(1:10)){
  n <- 10 + i
  set.seed(n)
  df_train<- sample_frac(df, 0.7)
  df_test <- subset(df, !(Estación %in% c(df_train$Estación)))
  
  tps <- Tps(df_train[,c(3,4,5)], df_train$estacion10)
  
  
  ###hacer la interpolación##
  
  p3 <- interpolate(ALTITUD, tps, xyOnly = F)
  
  plot(p3, main = paste0('Interpolación iteración ', i, '  ', 'Tenerife', ' ', tecnica))
  
  ###exportar el raster de la interpolación###
  
  # writeRaster(p3, filename = paste0(dst_folder,'TPS_LOOCV','_',i,'.tif'), format = 'Gtiff', overwrite = T)
  
  
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

df_validacion_loocv <- data.frame(df$estacion10, df_test_plot$Pred)
r2 <- metrica::R2(data = df_validacion_loocv, obs = estacion10, pred = df_test_plot.Pred)

###JLME##

fit<-lm(df_validacion_loocv$df.estacion10 ~ df_validacion_loocv$df_test_plot.Pred, data=df_validacion_loocv)
R2<-format(summary(fit)$adj.r.squared, digits=2)


tiles_plot(data = df_validacion_loocv, obs = df.estacion10, pred = df_test_plot.Pred,
           bins = 10,
           colors = c(low = 'yellow', high = 'green'))


####bucle###

for (i in seq(1:10)){
  set.seed(i)
  df_train<- dplyr::sample_frac(df, 0.7)
  df_test <- dplyr::sample_frac(df, 0.3)
  train.control <- trainControl(method = "LOOCV")
  
  step.model1 <- train(estacion10 ~ XUTM + YUTM + ALTITUD , data = df_train,
                       method = "lmStepAIC", 
                       trControl = train.control,
                       trace = FALSE
  )
  
  YUTM <- raster('../RASTERS_ENTRADA/YUTM.tif')
  XUTM <- raster('../RASTERS_ENTRADA/YUTM.tif')
  ALTITUD <- raster('../RASTERS_ENTRADA/DEM_TENERIFE.tif')
  names(ALTITUD) <- 'ALTITUD'
  XUTM <- raster('../RASTERS_ENTRADA/XUTM.tif')
  
  r <- raster::stack(YUTM, ALTITUD, XUTM)
  
  prediction_model1 <- predict(r, step.model1)
  
  plot(prediction_model1)
  
  ###predichos y observados###
  pred <- step.model1$pred
  df_train_val <- cbind(df_train, pred)
  
  df_train_val$Residuos <- df_train_val$estacion10 - df_train_val$pred
  
  coordinates(df_train_val) <- ~ XUTM + YUTM
  
  proj4string(df_train_val) <- CRS("+proj=utm +zone=28 +datum=WGS84 +units=m +no_defs")
  
  crs_stations <- CRS(SRS_string = "EPSG:32628")
  slot(df_train_val, "proj4string") <- crs_stations
  g <- automap::autofitVariogram(formula = Residuos ~ 1, input_data = df_train_val)
  
  plot(g)
  g <- gstat::gstat(formula = Residuos ~ 1, data = df_train_val, model = g$var_model)
  
  z1 <- interpolate(prediction_model1, g) ###el prediction_model1 sólo se usa como grid
  
  plot(z1)
  
  ####MLR + OK
  
  mlr_ok <- z1 + prediction_model1
  
  superficies <- raster::brick(prediction_model1, mlr_ok)
  
  plot(superficies)
  
  writeRaster(mlr_ok, filename = paste0(dst_folder,'mlr_ok','_',i,'.tif'), format = 'Gtiff', overwrite = T)
  
  io_interpolado_train <- raster::extract(mlr_ok, df_train[3:4])
  
  df_validacion_train <- data.frame(df_train$estacion10, io_interpolado_train)
  
  r2_train <- metrica::R2(df_validacion_train, obs = df_validacion_train$df_train.estacion10, pred = df_validacion_train$io_interpolado_train)
  rmse_train <- metrica::RMSE(df_validacion_train, obs = df_validacion_train$df_train.estacion10, pred = df_validacion_train$io_interpolado_train)
  mae_train <- metrica::MAE(df_validacion_train, obs = df_validacion_train$df_train.estacion10, pred = df_validacion_train$io_interpolado_train)
  df_train_plot[nrow(df_train_plot) + 1,] = c(i,r2_train, rmse_train, mae_train)
  
  ###test###
  io_interpolado_test <- raster::extract(mlr_ok, df_test[3:4])
  df_validacion_test <- data.frame(df_test$estacion10, io_interpolado_test)
  r2_test <- metrica::R2(df_validacion_test, obs = df_validacion_test$df_test.estacion10, pred = df_validacion_test$io_interpolado_test)
  rmse_test <- metrica::RMSE(df_validacion_test, obs = df_validacion_test$df_test.estacion10, pred = df_validacion_test$io_interpolado_test)
  mae_test <- metrica::MAE(df_validacion_test, obs = df_validacion_test$df_test.estacion10, pred = df_validacion_test$io_interpolado_test)
  df_test_plot[nrow(df_test_plot) + 1,] = c(i,r2_test, rmse_test, mae_test)
  
  

}

write.csv(df_test_plot, 'C:/BIOCLIMATOLOGIA_TENERIFE/validacion/test.csv')
  
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


####MODEL2 : ITC ###

train.control <- trainControl(method = "LOOCV")

step.model2 <- train(TP ~ XUTM + YUTM + Elev1 , data = df1,
                                    method = "lmStepAIC", 
                                    trControl = train.control,
                                    trace = FALSE)
  
# Model accuracy
step.model2$results
# Final model coefficients
step.model2$finalModel
# Summary of the model
summary(step.model2$finalModel)


prediction_model2 <- predict(r, step.model2)

plot(prediction_model2)

df_validacion <- extract(prediction_model2, df1[2:3])

df_final_val <- cbind(df1, df_validacion)

df_final_val$Residuos <- df_final_val$TP - df_final_val$df_validacion

coordinates(df_final_val) <- ~ XUTM + YUTM

proj4string(df_final_val) <- CRS("+proj=utm +zone=28 +datum=WGS84 +units=m +no_defs")

crs_stations <- CRS(SRS_string = "EPSG:32628")
slot(df_final_val, "proj4string") <- crs_stations
g <- automap::autofitVariogram(formula = Residuos ~ 1, input_data = df_final_val)

plot(g)
g <- gstat::gstat(formula = Residuos ~ 1, data = df_final_val, model = g$var_model)

z1 <- interpolate(prediction_model2, g) ###el prediction_model1 sólo se usa como grid


plot(z1)

####MLR + OK

mlr_ok <- z1 + prediction_model2

superficies <- raster::brick(prediction_model2, mlr_ok)

plot(superficies)

writeRaster(mlr_ok, filename = '../SUPERFICIES_INTERPOLADAS/Tp.tif', format = 'Gtiff')
















####Model 3: Tp

step.model3 <- train(TP ~ XUTM + YUTM + Elev1 , data = df1,
                     method = "lmStepAIC", 
                     trControl = train.control,
                     trace = FALSE)

# Model accuracy
step.model3$results
# Final model coefficients
step.model3$finalModel
# Summary of the model
summary(step.model3$finalModel)