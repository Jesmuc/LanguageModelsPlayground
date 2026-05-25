## -----------------------------------------------
## Metropolis–Hastings con historial para gráficos
## -----------------------------------------------

mh_sampler <- function(log_target,
                       propose,
                       x_inicial,
                       n_muestras,
                       burn_in = 0,
                       thin = 1,
                       semilla = NULL) {
  if (!is.null(semilla)) set.seed(semilla)
  
  n_a_guardar <- n_muestras
  n_total <- burn_in + n_a_guardar * thin
  
  # Historial COMPLETO (útil para traza/diagnóstico)
  cadena_cruda <- vector(mode = "list", length = n_total)
  aceptado_flag <- logical(n_total)
  
  # Almacén de muestras "post-procesadas" (después de burn-in y thinning)
  muestras <- vector(mode = "list", length = n_a_guardar)
  
  # Estado inicial
  x_actual <- x_inicial
  log_pi_x <- log_target(x_actual)
  
  k_guardado <- 0L
  aceptas <- 0L
  
  for (t in seq_len(n_total)) {
    # Propuesta
    prop <- propose(x_actual)
    y <- prop$y
    log_q_xy <- prop$log_q_xy
    log_q_yx <- prop$log_q_yx
    
    # Razón de aceptación (en log)
    log_pi_y <- log_target(y)
    log_alpha <- (log_pi_y - log_pi_x) + (log_q_yx - log_q_xy)
    
    # Aceptar / Rechazar
    if (log(runif(1)) < log_alpha) {
      x_actual <- y
      log_pi_x <- log_pi_y
      aceptas <- aceptas + 1L
      aceptado_flag[t] <- TRUE
    }
    
    # Registrar estado actual en la cadena cruda (incluye burn-in)
    cadena_cruda[[t]] <- x_actual
    
    # Guardar según burn-in y thin
    if (t > burn_in && ((t - burn_in) %% thin == 0)) {
      k_guardado <- k_guardado + 1L
      muestras[[k_guardado]] <- x_actual
    }
  }
  
  # Limpiezas: convertir listas a vectores si son escalares
  if (all(vapply(muestras, length, integer(1)) == 1L)) {
    muestras <- unlist(muestras, use.names = FALSE)
  }
  if (all(vapply(cadena_cruda, length, integer(1)) == 1L)) {
    cadena_cruda <- unlist(cadena_cruda, use.names = FALSE)
  }
  
  list(
    muestras = muestras,                  # posterior a burn-in/thin
    cadena_cruda = cadena_cruda,          # toda la cadena (para trazas)
    aceptado_flag = aceptado_flag,        # TRUE/FALSE por iteración
    tasa_aceptacion = mean(aceptado_flag),
    iteraciones_totales = n_total,
    burn_in = burn_in,
    thin = thin,
    inicial = x_inicial
  )
}
## Discreto: caminata en enteros con límites (±1, borde maneja asimetría)
propuesta_entera_rw <- function(paso = 1L, lower = 1L, upper = Inf) {
  force(paso); force(lower); force(upper)
  function(x_actual) {
    cand <- c(x_actual - paso, x_actual + paso)
    ok <- (cand >= lower) & (cand <= upper)
    cand <- cand[ok]
    
    if (length(cand) == 1L) {
      y <- cand
      log_q_xy <- 0
      vec_y <- c(y - paso, y + paso)
      ok_y <- (vec_y >= lower) & (vec_y <= upper)
      log_q_yx <- -log(sum(ok_y))
    } else {
      y <- sample(cand, 1L)
      log_q_xy <- -log(length(cand))
      vec_y <- c(y - paso, y + paso)
      ok_y <- (vec_y >= lower) & (vec_y <= upper)
      log_q_yx <- -log(sum(ok_y))
    }
    list(y = y, log_q_xy = log_q_xy, log_q_yx = log_q_yx)
  }
}

## Continuo: caminata Normal simétrica
propuesta_normal_rw <- function(sigma_prop = 1) {
  force(sigma_prop)
  function(x_actual) {
    y <- rnorm(1L, mean = x_actual, sd = sigma_prop)
    list(y = y, log_q_xy = 0, log_q_yx = 0)  # simétrica
  }
}
## ------------------------------
## Gráficos genéricos de diagnóstico
## ------------------------------

## 3.1 Traza de la cadena (incluye burn-in visualmente)
traza_mh <- function(cadena_cruda, burn_in = 0,
                     main = "Traza de la cadena (MH)",
                     ylab = "Estado", xlab = "Iteración") {
  plot(cadena_cruda, type = "l", col = "gray40", main = main, ylab = ylab, xlab = xlab)
  if (burn_in > 0) abline(v = burn_in, col = "tomato", lty = 2)
}

## 3.2 Media corrida (convergencia de la media muestral)
media_corrida <- function(cadena_cruda, burn_in = 0,
                          main = "Media corrida", ylab = "Media hasta t") {
  mc <- cumsum(cadena_cruda) / seq_along(cadena_cruda)
  plot(mc, type = "l", col = "steelblue", main = main, ylab = ylab, xlab = "Iteración")
  if (burn_in > 0) abline(v = burn_in, col = "tomato", lty = 2)
}

## 3.3 Autocorrelación de las muestras post-proceso
acf_mh <- function(muestras, lag_max = 50, main = "ACF de muestras") {
  acf(muestras, lag.max = lag_max, main = main)
}

## 3.4 Tasa de aceptación acumulada
tasa_aceptacion_acumulada <- function(aceptado_flag,
                                      main = "Tasa de aceptación acumulada") {
  tasa <- cumsum(aceptado_flag) / seq_along(aceptado_flag)
  plot(tasa, type = "l", col = "darkgreen", main = main,
       ylab = "Aceptación acumulada", xlab = "Iteración")
  abline(h = mean(aceptado_flag), col = "gray60", lty = 2)
}

## ------------------------------------------
## Overlays contra la distribución objetivo
## ------------------------------------------

## 3.5 Objetivo DISCRETO: pmf empírica vs teórica (opcionalmente en log-log)
overlay_discreto <- function(muestras, dtarget_pmf, soporte,
                             log_log = FALSE,
                             main = "PMF empírica vs teórica",
                             ylab = "Probabilidad") {
  # pmf empírica (normalizada)
  tab <- prop.table(table(factor(muestras, levels = soporte)))
  pmf_emp <- as.numeric(tab)
  # pmf teórica en el mismo soporte (asegúrate que dtarget_pmf devuelve pmf normalizada)
  pmf_teo <- vapply(soporte, dtarget_pmf, numeric(1))
  
  if (!log_log) {
    plot(soporte, pmf_teo, type = "h", lwd = 2, col = "gray60",
         main = main, xlab = "k", ylab = ylab)
    points(soporte, pmf_emp, pch = 16, col = "steelblue")
    legend("topright", c("Teórica", "Empírica"),
           col = c("gray60", "steelblue"), lwd = c(2, NA), pch = c(NA, 16), bty = "n")
  } else {
    plot(soporte, pmf_teo, log = "xy", type = "l", lwd = 2, col = "gray60",
         main = paste(main, "(log-log)"), xlab = "k (log)", ylab = "PMF (log)")
    points(soporte, pmf_emp, pch = 16, col = "steelblue")
    legend("topright", c("Teórica", "Empírica"),
           col = c("gray60", "steelblue"), lwd = c(2, NA), pch = c(NA, 16), bty = "n")
  }
}

## 3.6 Objetivo CONTINUO: densidad empírica (kernel) vs teórica
overlay_continuo <- function(muestras, dtarget_pdf,
                             main = "Densidad empírica vs teórica",
                             xlab = "x", ylab = "Densidad") {
  dens_emp <- density(muestras)
  plot(dens_emp, main = main, xlab = xlab, ylab = ylab, col = "steelblue", lwd = 2)
  curve(dtarget_pdf(x), add = TRUE, col = "gray40", lwd = 2)
  legend("topright", c("Empírica (kernel)", "Teórica"),
         col = c("steelblue", "gray40"), lwd = 2, bty = "n")
}

# Ejemplo de uso: Ley de Potencias
## Objetivo discreto: log-objetivo (sin constante)
log_target_potencia <- function(k, alpha = 1.5) {
  if (length(k) != 1L || is.na(k) || k < 1 || k != as.integer(k)) return(-Inf)
  -alpha * log(k)
}

## pmf teórica normalizada en {1,...,Kmax} para overlay (truncada)
pmf_power_trunc <- function(alpha, Kmax) {
  z <- sum((1:Kmax)^(-alpha))
  function(k) ifelse(k >= 1 & k <= Kmax, k^(-alpha) / z, 0)
}

set.seed(2025)
fit_pow <- mh_sampler(
  log_target = function(z) log_target_potencia(z, alpha = 1.5),
  propose    = propuesta_entera_rw(paso = 1L, lower = 1L, upper = Inf),
  x_inicial  = 2L,
  n_muestras = 5e4,
  burn_in    = 1000,
  thin       = 1
)

## Gráficos de diagnóstico
par(mfrow = c(2,2))
traza_mh(fit_pow$cadena_cruda, burn_in = fit_pow$burn_in,
         main = "Traza (potencia)")
media_corrida(fit_pow$cadena_cruda, burn_in = fit_pow$burn_in)
acf_mh(fit_pow$muestras, lag_max = 60, main = "ACF (potencia)")
tasa_aceptacion_acumulada(fit_pow$aceptado_flag)
par(mfrow = c(1,1))

## Overlay empírico vs teórico (truncamos soporte para graficar)
# Kmax <- 200
# pmf_teo_fun <- pmf_power_trunc(alpha = 1.5, Kmax = Kmax)
# overlay_discreto(fit_pow$muestras, dtarget_pmf = pmf_teo_fun,
#                  soporte = 1:Kmax, log_log = TRUE,
#                  main = "Power-law α=1.5 (truncada a 200)")
