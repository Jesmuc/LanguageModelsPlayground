##############################################
# 10. Evaluación de la precisión en n oraciones del corpus
##############################################

n_oraciones <- 200   # <-- cambia esto (50, 200, 1000, etc.)

# Aseguramos que no exceda el tamaño del corpus
n_oraciones <- min(n_oraciones, length(spanish_docs))

acc_por_oracion <- numeric(n_oraciones)
tokens_totales <- 0
aciertos_totales <- 0

# (Opcional) guardar algunas comparaciones para inspección
comparaciones_ejemplo <- list()

for (i in seq_len(n_oraciones)) {
  
  sent_text <- spanish_docs[i]
  
  # Etiquetamos con udpipe (referencia)
  sent_annot <- udpipe_annotate(ud_model, x = sent_text, doc_id = paste0("eval_", i))
  sent_df <- as.data.frame(sent_annot)
  
  # Filtrar solo tokens con letras
  sent_df <- subset(sent_df, grepl("[[:alpha:]]", token))
  
  # Si no quedan tokens, saltamos
  if (nrow(sent_df) == 0) {
    acc_por_oracion[i] <- NA
    next
  }
  
  words_sent <- tolower(sent_df$token)
  upos_sent  <- sent_df$upos
  
  # ---- Predicción HMM + Viterbi (manejo rápido de casos cortos) ----
  if (length(words_sent) == 1) {
    # Caso de 1 token: elegimos argmax pi(tag)*P(w|tag)
    w1 <- words_sent[1]
    w_idx <- match(w1, word_levels)
    
    # log P(w|tag)
    log_emit <- if (is.na(w_idx)) rep(log(1e-8), length(tag_levels)) else log(B[, w_idx])
    
    scores <- log(pi) + log_emit
    pred_tags <- tag_levels[which.max(scores)]
    
  } else {
    # Caso normal: usamos tu Viterbi tal cual
    v <- viterbi_hmm(words_sent, tag_levels, word_levels, pi, Gamma, B)
    pred_tags <- v$tags
  }
  
  # Comparación token a token
  match_vec <- (pred_tags == upos_sent)
  
  # Accuracy de esta oración
  acc_i <- mean(match_vec)
  acc_por_oracion[i] <- acc_i
  
  # Acumulado global (por token)
  tokens_totales <- tokens_totales + length(match_vec)
  aciertos_totales <- aciertos_totales + sum(match_vec)
  
  # Guardar algunas oraciones para inspección (por ejemplo las primeras 3)
  if (i <= 3) {
    comparaciones_ejemplo[[i]] <- data.frame(
      oracion_id  = i,
      word        = words_sent,
      HMM_tag     = pred_tags,
      udpipe_tag  = upos_sent,
      match       = match_vec,
      stringsAsFactors = FALSE
    )
  }
}

# Resultados
cat("\n=== Evaluación en", n_oraciones, "oraciones ===\n")

# Accuracy micro (por token)
cat("Accuracy MICRO (por token): ",
    if (tokens_totales > 0) aciertos_totales / tokens_totales else NA,
    "\n")

# Accuracy macro (promedio por oración, ignorando NA)
cat("Accuracy MACRO (promedio por oración): ",
    mean(acc_por_oracion, na.rm = TRUE),
    "\n")

cat("Oraciones válidas (no vacías): ", sum(!is.na(acc_por_oracion)), "\n")

# Mostrar algunas comparaciones
# if (length(comparaciones_ejemplo) > 0) {
#   cat("\n--- Ejemplos (primeras oraciones) ---\n")
#   print(comparaciones_ejemplo[[1]])
#   if (length(comparaciones_ejemplo) >= 2) print(comparaciones_ejemplo[[2]])
#   if (length(comparaciones_ejemplo) >= 3) print(comparaciones_ejemplo[[3]])
# }