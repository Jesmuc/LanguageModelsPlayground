## ------------------------------------------------------------
## Poisson–HMM con muestreador de Gibbs (basado en Scott / Frühwirth-Schnatter)
## Datos: Tabla 1.1 – Terremotos mayores (Mag ≥ 7), 1900–2006 (107 obs)
## Modelo:
##   C_t | C_{t-1} ~ Markov(Γ) en estados {1,...,m}
##   X_t | C_t=i    ~ Poisson(λ_i), con λ_i = sum_{j=1}^i τ_j, τ_j >= 0
## Priori:
##   Fila r de Γ ~ Dirichlet(ν_r)      (independientes por fila)
##   τ_j        ~ Gamma(a_j, b_j)      (forma a_j, tasa b_j; independientes)
## Pasos del Gibbs (como en las páginas adjuntas):
##   1) Muestrear C_{1:T} | x, θ con mensajes hacia adelante + muestreo hacia atrás
##   2) Descomponer cada x_t en contribuciones (x_{1t},...,x_{it}) vía Multinomial
##   3) Actualizar Γ por filas con Dirichlet(ν_r + conteos de transición)
##   4) Actualizar τ_j con Gamma(a_j + sum_t x_{jt},  b_j + N_j),
##      donde N_j = #{t : C_t ≥ j} (régimen j activo en t)
## ------------------------------------------------------------

set.seed(1)

## -----------------------
## 0) Datos de terremotos (Tabla 1.1)
##    (leer por renglones; 107 conteos anuales)
## -----------------------
x <- c(
  13,14,8,10,16,26,32,27,18,32,36,24,22,23,22,18,25,21,21,14,
  8,11,14,23,18,17,19,20,22,19,13,26,13,14,23,22,24,21,22,26,21,
  23,24,27,41,31,27,35,26,28,36,39,21,17,22,17,19,15,34,10,15,
  22,18,15,20,15,22,19,16,30,27,29,23,20,16,21,25,16,18,15,
  18,14,10,15,8,15,6,11,8,7,18,16,13,12,13,20,15,16,12,18,
  15,16,13,15,16,11,11
)

T <- length(x)  # debería ser 107

## -----------------------
## 1) Tamaño del modelo y priori
## -----------------------
m <- 3  # número de estados ocultos (puedes probar m=2,3,4)

# Hiperparámetros Dirichlet para las filas de Γ (ν_r es vector por fila)
# Usamos priori débilmente informativa y simétrica.
nu <- matrix(1, nrow = m, ncol = m)   # ν_rj = 1 (uniforme en el simplex)

# Priori Gamma(a_j, b_j) para τ_j (forma a_j, tasa b_j).
# Priori ligeras que permiten λ ordenadas pero flexibles.
a <- rep(1.0, m)     # forma
b <- rep(0.1, m)     # tasa  (media a/b = 10; ajusta si quieres)

## -----------------------
## 2) Auxiliares: muestreo y utilidades
## -----------------------

# Muestreador Dirichlet (para una fila)
rdirichlet_row <- function(alpha) {
  g <- rgamma(length(alpha), shape = alpha, rate = 1)
  g / sum(g)
}

# PMF de Poisson vectorizada sobre un conjunto de medias
dpois_vec <- function(y, lambda) {
  # regresa vector de largo length(lambda): P(Y=y | λ=lambda[i])
  dpois(y, lambda = lambda)
}

# De τ a λ (acumulada, impone orden λ_1 ≤ λ_2 ≤ ...)
tau_to_lambda <- function(tau) cumsum(tau)

# Mensajes hacia adelante α_t(i) ∝ P(x_{1:t}, C_t=i | θ), normalizados por t.
# Regresa list(alpha = matriz T x m de α normalizados,
#              loglik = suma de log c_t con c_t las constantes de escalamiento)
forward_messages <- function(x, lambda, Gamma) {
  T <- length(x); m <- length(lambda)
  alpha <- matrix(0, nrow = T, ncol = m)
  cscale <- numeric(T)
  
  # Distribución inicial: estacionaria de Γ o uniforme (usamos uniforme)
  pi0 <- rep(1/m, m)
  
  # t = 1
  emit1 <- dpois_vec(x[1], lambda)
  alpha1_unnorm <- pi0 * emit1
  cscale[1] <- sum(alpha1_unnorm)
  alpha[1, ] <- alpha1_unnorm / cscale[1]
  
  # t >= 2
  for (t in 2:T) {
    emit <- dpois_vec(x[t], lambda)
    pred <- as.numeric(alpha[t-1, ] %*% Gamma)   # longitud m
    at_unnorm <- pred * emit
    cscale[t] <- sum(at_unnorm)
    alpha[t, ] <- at_unnorm / cscale[t]
  }
  list(alpha = alpha, loglik = sum(log(cscale)))
}

# Muestreo hacia atrás de la trayectoria oculta usando α y Γ
backward_sample_states <- function(alpha, Gamma) {
  T <- nrow(alpha); m <- ncol(alpha)
  C <- integer(T)
  
  # Muestrea C_T ~ Categórica(α_T)
  C[T] <- sample.int(m, size = 1, prob = alpha[T, ])
  
  # Para t = T-1,...,1:
  for (t in (T-1):1) {
    # P(C_t = i | C_{t+1}=j, x) ∝ α_t(i) Γ_{ij}
    prob <- alpha[t, ] * Gamma[, C[t+1]]
    prob <- prob / sum(prob)
    C[t] <- sample.int(m, size = 1, prob = prob)
  }
  C
}

# Construye la matriz de conteos de transición a partir de una trayectoria de estados
transition_counts <- function(C, m) {
  M <- matrix(0, m, m)
  for (t in 1:(length(C)-1)) {
    M[C[t], C[t+1]] <- M[C[t], C[t+1]] + 1
  }
  M
}

# Descomponer x_t en contribuciones de regímenes activos cuando C_t = i.
# Muestrea (x_{1t},...,x_{it}) ~ Multinomial(x_t, p_j ∝ τ_{1:i})
# y fija x_{jt}=0 para j>i. Regresa vector de tamaño m.
decompose_xt <- function(xt, i, tau) {
  m <- length(tau)
  contrib <- integer(m)
  if (xt > 0 && i >= 1) {
    probs <- tau[1:i]
    probs <- probs / sum(probs)
    draw <- rmultinom(1, size = xt, prob = probs)[,1]
    contrib[1:i] <- draw
  }
  contrib
}

## -----------------------
## 3) Inicialización
## -----------------------
# Γ inicial con filas Dirichlet aleatorias
Gamma <- t(apply(nu, 1, rdirichlet_row))

# τ inicial desde la priori; calcular λ
tau <- rgamma(m, shape = a, rate = b)
lambda <- tau_to_lambda(tau)

# Trayectoria inicial: aleatoria (podrías usar clustering simple)
C <- sample.int(m, size = T, replace = TRUE)

## -----------------------
## 4) Parámetros del Gibbs
## -----------------------
n_iter   <- 4000
burn_in  <- 1500
thin     <- 2

keep_idx <- seq(burn_in + 1, n_iter, by = thin)
S <- length(keep_idx)

# Almacenamiento
store_lambda <- matrix(NA_real_, nrow = S, ncol = m)
store_tau    <- matrix(NA_real_, nrow = S, ncol = m)
store_Gamma  <- array(NA_real_, dim = c(S, m, m))
store_loglik <- numeric(S)

## -----------------------
## 5) Muestreador de Gibbs
## -----------------------
s <- 0
for (it in 1:n_iter) {
  
  ## ---- (a) Muestrear C_{1:T} | x, θ con FF/BS
  fw <- forward_messages(x, lambda, Gamma)
  C  <- backward_sample_states(fw$alpha, Gamma)
  
  ## ---- (b) Descomponer conteos x_t en contribuciones x_{jt}
  # Construir matriz de contribuciones Xcontrib (m x T) y N_j activos
  Xcontrib <- matrix(0L, nrow = m, ncol = T)
  Nj <- integer(m)  # veces que el régimen j está activo (C_t ≥ j)
  
  for (t in 1:T) {
    i <- C[t]
    # el régimen j está activo sii j ≤ i
    if (i >= 1) {
      contrib <- decompose_xt(x[t], i, tau)
      Xcontrib[, t] <- contrib
      # Contar activos: para j=1..i sumar 1
      if (i >= 1) Nj[1:i] <- Nj[1:i] + 1
    }
  }
  
  ## ---- (c) Actualizar Γ por filas desde Dirichlet(ν_r + conteos de transición)
  Tr <- transition_counts(C, m)  # matriz m x m de conteos
  for (r in 1:m) {
    Gamma[r, ] <- rdirichlet_row(nu[r, ] + Tr[r, ])
  }
  
  ## ---- (d) Actualizar τ_j ~ Gamma(a_j + sum_t x_{jt},  b_j + N_j)
  shape_post <- a + rowSums(Xcontrib)
  rate_post  <- b + Nj
  tau <- rgamma(m, shape = shape_post, rate = rate_post)
  
  ## ---- (e) Recalcular λ (medias ordenadas)
  lambda <- tau_to_lambda(tau)
  
  ## ---- (f) Guardar draws después de burn-in/thinning
  if (it %in% keep_idx) {
    s <- s + 1
    store_lambda[s, ] <- lambda
    store_tau[s, ]    <- tau
    store_Gamma[s, ,] <- Gamma
    store_loglik[s]   <- fw$loglik
  }
  
  if (it %% 500 == 0) {
    cat(sprintf("Iter %4d | loglik %.2f | lambda: %s\n",
                it, fw$loglik, paste(round(lambda,2), collapse=", ")))
  }
}

## -----------------------
## 6) Resúmenes posteriores
## -----------------------
post_lambda_mean <- colMeans(store_lambda)
post_lambda_ci   <- t(apply(store_lambda, 2, quantile, probs = c(0.05,0.95)))

post_Gamma_mean  <- apply(store_Gamma, c(2,3), mean)

cat("\nMedia posterior de lambda (medias por estado):\n")
print(round(post_lambda_mean, 3))
cat("\nIntervalos posteriores 90% para lambda:\n")
print(round(post_lambda_ci, 3))

cat("\nMatriz de transición Γ (media posterior):\n")
print(round(post_Gamma_mean, 3))

## -----------------------
## 7) Opcional: vistazo a estados decodificados
##    (estado MAP puntual usando los α del último promedio)
## -----------------------
fw_last <- forward_messages(x, colMeans(store_lambda), apply(store_Gamma, c(2,3), mean))
map_states <- max.col(fw_last$alpha, ties.method = "first")

# Plots rápidos en base R (puedes cambiar a ggplot si prefieres)
par(mfrow = c(2,2), mar=c(4,4,2,1))

# (i) Datos con bandas por estado MAP
plot(x, type="h", lwd=2, main="Terremotos por año (con estados MAP)",
     xlab="Tiempo (años 1900–2006)", ylab="Conteo")
cols <- c("#56B4E9","#E69F00","#009E73","#CC79A7")[1:m]
rect_x <- 1:T
rect_y0 <- rep(par("usr")[3], T)
rect_y1 <- rep(par("usr")[4], T)
for (t in 1:T) rect(t-0.5, rect_y0[t], t+0.5, rect_y1[t], col=adjustcolor(cols[map_states[t]], 0.08), border=NA)

# (ii) Traza de la verosimilitud logarítmica
plot(store_loglik, type="l", main="Traza: log-verosimilitud", xlab="Iteración guardada", ylab="log p(x|θ)")

# (iii) Trazas de λ
matplot(store_lambda, type="l", lty=1, main="Trazas: λ (medias por estado)", xlab="Iteración guardada", ylab="λ")
legend("topright", legend=paste0("λ",1:m), col=1:m, lty=1, bty="n")

# (iv) Densidades posteriores de λ
plot(density(store_lambda[,1]), main="Densidades posteriores de λ", xlab="λ", ylim=c(0,max(sapply(1:m, function(j) max(density(store_lambda[,j])$y)))))
for (j in 2:m) lines(density(store_lambda[,j]), col=j)
legend("topright", legend=paste0("λ",1:m), col=1:m, lty=1, bty="n")
