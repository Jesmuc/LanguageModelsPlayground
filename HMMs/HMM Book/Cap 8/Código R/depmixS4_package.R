# /////////////////////////////////////////////////////////////
# //////////////// Usando depmixS4 Package ////////////////////
# /////////////////////////////////////////////////////////////
# Ejemplos del libro HMMs For Time Series, Capitulo 8

library(depmixS4)

options(scipen = 999)   # evitar notacion cientifica
# //////////////// DATOS ////////////////////
# Leer los datos de terremotos desde la URL proporcionada
#quakes <- read.table(
#  "http://www.hmms-for-time-series.de/second/data/earthquakes.txt")$V2

quakes <- c(
  13,14,8,10,16,26,32,27,18,32,36,24,22,23,22,18,25,21,21,14,
  8,11,14,23,18,17,19,20,22,19,13,26,13,14,22,24,21,22,26,21,
  23,24,27,41,31,27,35,26,28,36,39,21,17,22,17,19,15,34,10,15,
  22,18,15,20,15,22,19,16,30,27,29,23,20,16,21,21,25,16,18,15,
  18,14,10,15,8,15,6,11,8,7,18,16,13,12,13,20,15,16,12,18,
  15,16,13,15,16,11,11
)

# //////////////// ESPECIFICAR MODELO ////////////////////
# quakes ~ 1: "~" se usa para definir la formula de emision, 
#             el "1" es para usar solo el intercepto del GLM: logλ_s,t = β_0,s +β_1,s x_1,t + β_2,s x_2,t (ver https://online.stat.psu.edu/stat504/lesson/6/6.1)
#             donde x1,t y x2,t son predictores/covariables/variables independientes
#             en caso de tener dichos predictores, por ejemplo, mediciones de temperatura (x1), humedad (x2),
#             se podrian incluir usando quakes ~ x1 + x2. Tambien es posible usar covariables para la transicion de estados usando transition=
# nstates   : Numero de estados del modelo
# ntimes    : Numero de observaciones en cada serie de tiempo individual, si no se especifica se asume que las observaciones forman una sola secuencia
mod_2 <- depmix(quakes ~ 1, nstates = 2, family = poisson(), ntimes = length(quakes))
mod_3 <- depmix(quakes ~ 1, nstates = 3, family = poisson(), ntimes = length(quakes))
mod_4 <- depmix(quakes ~ 1, nstates = 4, family = poisson(), ntimes = length(quakes))

# //////////////// ESTIMACION DE PARAMETROS ////////////////////
# Calcular MLEs (Maximum Likelihood Estimates)
#set.seed(1)
cat("\nMaximum Likelihood Estimates, m = 2:\n")
fitted_model_2 <- fit(mod_2)                                   # usa EM por default si no hay restricciones en los parametros
summary(fitted_model_2)
theta2 <- getpars(fitted_model_2)                              # guardar todos los parametros en un vector
vals <- unname(unlist(theta1[names(theta2) == "(Intercept)"])) # obtener solo los estimadores de logλ_s
cat("\nLambdas estimadas:\n")
lambda2 <- exp(vals)                                            # convertir intercepts a λ_s
print(lambda2)

cat("\nMaximum Likelihood Estimates, m = 3:\n")
fitted_model_3 <- fit(mod_3)
summary(fitted_model_3)
theta3 <- getpars(fitted_model_3)                              # guardar todos los parametros en un vector
vals <- unname(unlist(theta3[names(theta3) == "(Intercept)"])) # obtener solo los estimadores de logλ_s
cat("\nLambdas estimadas:\n")
lambda3 <- exp(vals)                                            # convertir intercepts a λ_s
print(lambda3)

cat("\nMaximum Likelihood Estimates, m = 4:\n")
fitted_model_4 <- fit(mod_4)
summary(fitted_model_4)
theta4 <- getpars(fitted_model_4)                              # guardar todos los parametros en un vector
vals <- unname(unlist(theta4[names(theta4) == "(Intercept)"])) # obtener solo los estimadores de logλ_s
cat("\nLambdas estimadas:\n")
lambda4 <- exp(vals)                                            # convertir intercepts a λ_s
print(lambda4)

# //////////////// DECODIFICACION ////////////////////
# Obtener probabilidades de estados ocultos usando algoritmo Forward-Backward
FB <- forwardbackward(fitted_model_2)$gamma
FB_maxIdx <- max.col(FB)
#plot(FB[ ,1], type = "o")
plot(FB_maxIdx, col= "brown", type = "o")
#abline(h=0.5, col="red")
#lines(FB[ ,2], col= "blue", type = "o")
# Obtener la secuencia de estados ocultos mas probable usando algoritmo Viterbi

viterbiPath <- posterior(fitted_model_2, type="viterbi")$state # type="smoothing"
points(viterbiPath, pch=20, col= "purple", type = "o")

#comparar FB_maxIdx con viterbiPath

cat("\nComparacion de secuencias de estados ocultos (Forward-Backward vs Viterbi):\n")
cat("\nSon Iguales?\n")
print(all(FB_maxIdx == viterbiPath))

# //////////////// GRAFICA MAP ////////////////////
#graficar serie de tiempo quakes coloreando cada punto segun el estado oculto estimado por Viterbi, grafico de barras

#plot(quakes, col=viterbiPath, type="h", main="Secuencia de estados ocultos estimada por Viterbi", xlab="Tiempo", ylab="Número de terremotos")

#plot(quakes, col=viterbiPath, pch=16, main="Secuencia de estados ocultos estimada por Viterbi", xlab="Tiempo", ylab="Número de terremotos")
#graficar segun FB_maxIdx
#plot(quakes, col=FB_maxIdx, type="h", main="Secuencia de estados ocultos estimada por Forward-Backward", xlab="Tiempo", ylab="Número de terremotos")

# ////////////// SELECCION DE MODELO ////////////////////
plot(2:4,c(BIC(fitted_model_2),BIC(fitted_model_3),BIC(fitted_model_4)),ty="b") # se interpreta como el mejor modelo aquel con el BIC mas bajo
title("Criterio BIC para modelos con diferente número de estados")


#////////////// NOTAS ////////////////////////////////////
#/// OTRAS FORMAS DE OBTENER DATOS DEL MODELO AJUSTADO ///
#slotNames(fitted_model_2)                    # ver slots disponibles
#summary(fitted_model_2, which = "response")  # imprimir solo los response parameters
#resp_list <- fitted_model_2@response          #igual que slot(fitted_model_2, "response")





