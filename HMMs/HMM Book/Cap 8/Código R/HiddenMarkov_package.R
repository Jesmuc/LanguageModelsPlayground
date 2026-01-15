# /////////////////////////////////////////////////////////////
# /////////////////// HiddenMarkov PACKAGE ////////////////////////
# /////////////////////////////////////////////////////////////
#  Ejemplos del libro HMMs For Time Series, Capitulo 8
# ////////////////////////////////////////////////////////////
library(HiddenMarkov)

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

# ///////////////// INICIALIZACION //////////////////////////////
Pi <- rbind(c(.5,.5),c(.5,.5))          # Matriz de transicion inicial
d <- c(.5,.5)                           # Distribucion inicial de estados
lambda <- c(10,30)                      # Parametros de Poisson iniciales

# Inicializacion 3 estados
Pi3 <- rbind(c(.8,.1,.1),c(.1,.8,.1),c(.1,.1,.8))  # Matriz de transicion inicial
d3 <- c(.33,.33,.34)                               # Distribucion inicial de estados
lambda3 <- c(5,15,30)                              # Parametros de Poisson iniciales

# ///////////////// ESTIMACION - EM ////////////////////////////
mod <- dthmm(quakes,Pi=Pi,delta=d,"pois", list(lambda=lambda))
fitted.mod.HM <- BaumWelch(mod)
summary(fitted.mod.HM)

mod3 <- dthmm(quakes,Pi=Pi3,delta=d3,"pois", list(lambda=lambda3))
fitted.mod.HM.M3 <- BaumWelch(mod3)
summary(fitted.mod.HM.M3)

# ///////////////////////////////////////////////////////////////
# //////////////  ESTIMACION MODELO ESTACIONARIO ////////////////
# ///////////////////////////////////////////////////////////////
# Se hace el fit pero ahora asumiendo que es un HMM estacionario 
# la distribucion inicial se estima como la distribucion estacionaria
#  obtenida de resolver δΓ = Γ
# ////////////////////////////////////////////////////////////// ///////////
fitted.mod.HM.stat <- dthmm(quakes,Pi=Pi,delta=d,"pois", list(lambda=lambda),nonstat=FALSE) # Discrete Time HMM Object
fitted.mod.HM.stat <- BaumWelch(fitted.mod.HM.stat)
print(BaumWelch(fitted.mod.HM.stat))
summary(fitted.mod.HM.stat)

# ///////////////// DECODIFICACION /////////////////////////////////////////
P_c <- fitted.mod.HM$u # Probabilidades de los estados ocultos
viterbiPath <- Viterbi(fitted.mod.HM) # Secuencia mas probable de estados ocultos

# ///////////////// RESIDUALES (pseudo-residuals o quantile residuals) /////
res <- residuals(fitted.mod.HM)    # Para Poisson–HMM en este caso son puntos 
w <- hist(res, main="Poisson HMM: Pseudo Residuales")
n <- length(res)
xfit <- seq(min(res), max(res), length = 40)
yfit <- dnorm(xfit, mean = mean(res), sd = sd(res))
yfit <- yfit * diff(w$mids[1:2]) * n
lines(xfit, yfit, col = "blue", lwd = 2)

# ///////////// Q-Q PLOT ////////////////////////////
qqnorm(res, main="Poisson HMM: Q-Q Plot of Pseudo Residuals")
abline(a=0, b=1, lty=3)
abline(h=seq(-2, 2, 1), lty=3)
abline(v=seq(-2, 2, 1), lty=3)
