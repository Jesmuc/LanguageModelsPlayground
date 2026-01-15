# /////////////////////////////////////////////////////////////
# /////////////////// msm PACKAGE ////////////////////////
# /////////////////////////////////////////////////////////////
#  Ejemplos del libro HMMs For Time Series, Capitulo 8
# ////////////////////////////////////////////////////////////
library(msm)

# ///////////////// DATOS //////////////////////////////
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

# ///////////////// ESTIMACIÓN MLE //////////////////////////////
timeq <- 1:length(quakes)

# Poisson-HMM 2 estados (continuo) - Matriz Q (ver sección 3.3.2)
# notas: f(x; λ) = λe^{-λx}, mean: 1/λ, memoryless property P(X > s + t | X > s) = P(X > t)
#                                                           usar definicion de prob condicional y A⊆B
fitted.mod.msm <- msm(quakes~timeq,
                      qmatrix=rbind(c(-.1,.1), c(.1,-.1)),
                      hmodel=list(hmmPois(rate=10),
                                  hmmPois(rate=20)))
cat("\nModelo Continuo:\n")
print(fitted.mod.msm)

#Γ = exp(Q), 
cat("\nModelo Discreto:\n")
print(pmatrix.msm(fitted.mod.msm))

# ///////////////// DECODIFICACION //////////////////////////////
cat("\nviterbi:\n")
print(viterbi.msm(fitted.mod.msm))