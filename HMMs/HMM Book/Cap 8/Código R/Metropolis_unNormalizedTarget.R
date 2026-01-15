## Algoritmo Metropolis-Hastings -- funcion objetivo no normalizada

set.seed(123)

# Densidad objetivo no normalizada (solo la conocemos hasta una constante)
target_unnorm <- function(x) {
  exp( -x^4 / 4 + x^2 / 2 )
}

N <- 200000                 # número de muestras
x <- numeric(N)
x[1] <- 0                   # iniciar en algún punto

prop_sd <- 0.5              # desviación estándar de la propuesta

for (i in 2:N) {
  # proponer desde una normal centrada en el estado actual
  x.c <- rnorm(1, mean = x[i-1], sd = prop_sd)
  
  # evaluar funcion objetivo no normalizada en el punto propuesto y en el actual
  p.c <- target_unnorm(x.c)
  p.p <- target_unnorm(x[i-1])
  
  # probabilidad de aceptación (la constante se cancela)
  alpha <- min(1, p.c / p.p)
  
  if (runif(1) < alpha) {
    x[i] <- x.c
  } else {
    x[i] <- x[i-1]
  }
}

## Descartar primeros elementos (burn-in)
burn <- 5000
x.use <- x[(burn+1):N]

## Graficar histograma como densidad
h <- hist(x.use,
          breaks = 120,
          freq = FALSE,
          col = "gray90",
          border = "gray70",
          main = "Muestreo de Metropolis",
          xlab = "x",
          ylim = c(0, max(yy)))

## Superponer la funcion objetivo (reescalada usando como referencia resultados del algoritmo)
xx <- seq(min(x.use), max(x.use), length.out = 400)
yy <- target_unnorm(xx)

# escalar el objetivo al histograma para visualización
scale_factor <- max(h$density) / max(yy)
lines(xx, yy * scale_factor, col = "blue", lwd = 2)
lines(xx, yy, col = "red", lwd = 2)
