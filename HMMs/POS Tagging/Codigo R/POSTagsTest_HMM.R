##############################################
# 1. Carga de librería y modelo de udpipe
##############################################

library(udpipe)
library(ggplot2)
library(matrixcalc)
library(MASS)

#Solo es necesario ejecutar las lineas de esta sección 1 vez

# Descarga de un modelo de udpipe para español.
model <- udpipe_download_model(language = "spanish")

# Cargamos el modelo descargado en memoria para usarlo en las anotaciones
ud_model <- udpipe_load_model(model$file_model)


##############################################
# 2. OBTENCION DE DATOS - Texto de entrenamiento en español
##############################################

#TEST 1 - Solo para probar
spanish_docs <- c(
  "El precio del petróleo subió hoy en los mercados internacionales.",
  "La inflación en México disminuyó ligeramente el mes pasado."
)

#TEST 2 - Se descartó por su estilo muy literario
# Leemos todas las líneas del archivo de texto directamente desde la web
#Unificamos todas las líneas en un solo string (un "documento" grande)

url <- "https://www.gutenberg.org/ebooks/29731.txt.utf-8"
lineas <- readLines(url, encoding = "UTF-8")
spanish_docs <- paste(lineas, collapse = " ")

#TEST 3 - Base de datos seleccionada
# Ajusta según path con la ubicación del archivo
leipzig_file <- "D://spa_wikipedia_2021_100K/spa_wikipedia_2021_100K-sentences.txt"

# Leemos todas las líneas (cada una es "ID \t oración")
lineas <- readLines(leipzig_file, encoding = "UTF-8")

# Separar ID y oración (usamos strsplit por tabulador)
parts <- strsplit(lineas, "\t")

# Extraemos solo el texto de la oración
spanish_docs <- vapply(parts, function(p) p[2], FUN.VALUE = character(1))

length(spanish_docs)      # ~100000 oraciones
head(spanish_docs, 3)     # inspección rápida

#TEST 4 - Al final solo se tomaron x oraciones
set.seed(123)
idx <- sample(seq_along(spanish_docs), size = 30000)  # por ejemplo 30000 oraciones
spanish_docs <- spanish_docs[idx]

##############################################
# 3. Anotación con udpipe (tokens + etiquetas POS)
#    CREACION DE BASE DE DATOS SINTETICA
##############################################

# Anotamos el texto: se añaden las etiquetas
annot <- udpipe_annotate(
  ud_model,
  x      = spanish_docs,
  doc_id = seq_along(spanish_docs)  # identificador de documento: 1, 2, ...
)

# Convertimos el resultado a data.frame
annot_df <- as.data.frame(annot)

# Inspección rápida de algunas columnas relevantes
head(annot_df[, c("doc_id", "sentence_id", "token_id",
                  "token", "lemma", "upos", "xpos")])


##############################################
# 4. POST PROCESAMIENTO DE DATOS SINTETICOS
##############################################

# Construimos un data.frame simplificado con:
# - doc_id: identificador del documento
# - word: la palabra en minúsculas
# - pos: la etiqueta POS universal (upos)
tokens <- data.frame(
  doc_id = annot_df$doc_id,
  word   = tolower(annot_df$token),
  pos    = annot_df$upos,    # o usar 'xpos' para etiquetas más finas
  stringsAsFactors = FALSE
)

# Nos quedamos solo con tokens que contienen letras (descartamos solo signos, etc.)
tokens <- subset(tokens, grepl("[[:alpha:]]", word))


##############################################
# 5. POST PROCESAMIENTO DE DATOS SINTETICOS PT2. 
#        Extracción de etiquetas y vocabulario
##############################################

# Vector de etiquetas POS (una por token)
tags  <- tokens$pos

# Vector de palabras observadas (una por token)
words <- tokens$word

# Conjunto finito de etiquetas (estados ocultos del HMM)
tag_levels  <- sort(unique(tags))

# Conjunto finito de palabras (vocabulario de observaciones del HMM)
word_levels <- sort(unique(words))

# Mapear cada etiqueta/palabra a un índice entero:
# tag_idx: para cada token, el índice de su etiqueta en 'tag_levels'
# word_idx: para cada token, el índice de su palabra en 'word_levels'
tag_idx  <- match(tags,  tag_levels)
word_idx <- match(words, word_levels)


##############################################
# 6. Estimación de π (distribución inicial)
##############################################

# Para la distribución inicial π(tag):
# usamos la primera etiqueta de cada documento como "estado inicial"
first_tag_each_doc <- tapply(tokens$pos, tokens$doc_id, function(v) v[1])

# Conteos de cuántas veces cada etiqueta es inicial
pi_counts <- table(factor(first_tag_each_doc, levels = tag_levels))

# Estimación de π con suavizado add-1:
# π(tag) = (conteo(tag) + 1) / suma_sobre_todos(conteo + 1)
pi <- (pi_counts + 1) / sum(pi_counts + 1)


##############################################
# 7. Estimación de Γ (matriz de transición)
##############################################

# Construimos pares de etiquetas consecutivas dentro de cada documento.
# from_idx: etiqueta en tiempo t
# to_idx:   etiqueta en tiempo t+1
from_idx <- tag_idx[-nrow(tokens)] #todos excepto ultimo
to_idx   <- tag_idx[-1]           #todos excepto primero

# same_doc indica si ambos tokens (t, t+1) pertenecen al mismo documento
same_doc <- tokens$doc_id[-nrow(tokens)] == tokens$doc_id[-1]

# Matriz de conteos de transiciones etiqueta_t -> etiqueta_{t+1}
Gamma_counts <- table(
  factor(from_idx[same_doc], levels = seq_along(tag_levels)),
  factor(to_idx[same_doc],   levels = seq_along(tag_levels))
)

# Suavizado add-1 para evitar ceros:
Gamma <- Gamma_counts + 1

options(scipen = 999)   # evitar notacion cientifica
# Normalizamos por filas para obtener probabilidades:
# Cada fila i representa P(tag_{t+1} | tag_t = i)
Gamma <- Gamma / rowSums(Gamma)


##############################################
# 8. Estimación de B (matriz de emisión)
##############################################

# Matriz de conteos de emisiones etiqueta -> palabra
B_counts <- table(
  factor(tag_idx,  levels = seq_along(tag_levels)),  # filas: etiquetas
  factor(word_idx, levels = seq_along(word_levels))  # columnas: palabras
)

# Suavizado add-1
B <- B_counts + 1

# Normalizamos por filas:
# Cada fila i representa P(palabra | etiqueta = i)
B <- B / rowSums(B)

#imprimir 5 columnas de B
#print(B[, 1:5])


##############################################
# 9. Implementación del algoritmo de Viterbi
##############################################

viterbi_hmm <- function(sentence_tokens, tag_levels, word_levels, pi, Gamma, B) {
  # sentence_tokens: vector de palabras (en minúsculas) de la oración a etiquetar
  # tag_levels:     vector con todas las etiquetas posibles (estados ocultos)
  # word_levels:    vocabulario de palabras visto en entrenamiento
  # pi:             distribución inicial π
  # Gamma:          matriz de transición
  # B:              matriz de emisión
  
  T_len  <- length(sentence_tokens)   # longitud de la secuencia
  N_tags <- length(tag_levels)        # número de estados (etiquetas)
  
  # Función auxiliar para logar probabilidades, manejando ceros como -Inf
  log_safe <- function(x) ifelse(x == 0, -Inf, log(x))
  
  # Logaritmos de las probabilidades de π, Γ y B
  log_pi    <- log_safe(pi)
  log_Gamma <- log_safe(Gamma)
  log_B     <- log_safe(B)
  
  # Mapear cada palabra de la oración a su índice en el vocabulario
  # (tendrá NA si la palabra no fue vista en el entrenamiento)
  w_idx <- match(sentence_tokens, word_levels)
  
  # Función auxiliar para obtener log P(palabra | etiqueta)
  get_emission_logprob <- function(tag_i, w_i) {
    if (is.na(w_i)) {
      # Si la palabra no existe en el vocabulario:
      # le damos una probabilidad muy pequeña pero no cero (palabra desconocida)
      return(log(1e-8))
    }
    # Si la palabra es conocida, usamos la matriz B
    log_B[tag_i, w_i]
  }
  
  # Matriz delta: delta[i, t] = log prob máxima de la mejor ruta que acaba en etiqueta i en tiempo t
  delta <- matrix(-Inf, nrow = N_tags, ncol = T_len)
  
  # Matriz psi: psi[i, t] = índice de la etiqueta anterior que lleva a delta[i, t] máximo
  psi   <- matrix(NA_integer_, nrow = N_tags, ncol = T_len)
  
  # ----------------------------
  # Inicialización (t = 1)
  # ----------------------------
  for (i in seq_len(N_tags)) {
    delta[i, 1] <- log_pi[i] + get_emission_logprob(i, w_idx[1])
    psi[i, 1]   <- 0  # no hay estado previo en t = 1
  }
  
  # ----------------------------
  # Recursión (t = 2, ..., T_len)
  # ----------------------------
  for (t in 2:T_len) {
    for (j in seq_len(N_tags)) {
      # Para llegar a etiqueta j en tiempo t, consideramos todas las etiquetas i en t-1
      # scores[i] = delta[i, t-1] + log P(j | i)
      scores <- delta[, t - 1] + log_Gamma[, j]
      
      # i_star = índice de la etiqueta previa que maximiza scores
      i_star <- which.max(scores)
      
      # Mejor log prob para estar en j en tiempo t
      delta[j, t] <- scores[i_star] + get_emission_logprob(j, w_idx[t])
      
      # Guardamos la "ruta" (de dónde venimos)
      psi[j, t] <- i_star
    }
  }
  
  # ----------------------------
  # Terminación y Backward Step
  # ----------------------------
  
  # 1) Terminación: elegimos la etiqueta final con mayor probabilidad en t = T_len
  best_last <- which.max(delta[, T_len])
  
  # 2) Backward step: reconstruimos la secuencia óptima de etiquetas hacia atrás
  best_path_idx <- integer(T_len)
  best_path_idx[T_len] <- best_last
  
  # Vamos hacia atrás desde t = T_len-1 hasta t = 1
  for (t in (T_len - 1):1) {
    best_path_idx[t] <- psi[best_path_idx[t + 1], t + 1]
  }
  
  # Devolvemos etiquetas (en su forma de texto) y el log-verosimilitud de la ruta óptima
  list(
    tags   = tag_levels[best_path_idx],
    loglik = delta[best_last, T_len]
  )
}


##############################################
# 10. Oración de prueba: udpipe vs HMM
##############################################

# Oración de prueba (en español, consistente con el modelo)
sent_text <- "Una copa de vino al dia no hace daño."
sent_text <- "Hoy voy a escribir diez ensayos sobre la historia de España."
sent_text <- "Hola cómo estás?"
sent_text <- "El modelo predice la etiqueta correcta para cada palabra"

# Etiquetamos la oración con udpipe (esto será nuestra "referencia")
sent_annot <- udpipe_annotate(ud_model, x = sent_text, doc_id = "test1")
sent_df    <- as.data.frame(sent_annot)

# Aplicamos el mismo filtrado: solo tokens con letras
# [[:alpha:]] significa cualquier letra (mayúscula o minúscula) 
sent_df <- subset(sent_df, grepl("[[:alpha:]]", token))

# Palabras y etiquetas POS (upos) de referencia
words_sent <- tolower(sent_df$token)
upos_sent  <- sent_df$upos

# Etiquetamos la oración con el HMM + Viterbi
v <- viterbi_hmm(words_sent, tag_levels, word_levels, pi, Gamma, B)

# Construimos un data.frame para comparar etiqueta a etiqueta
comparison <- data.frame(
  word       = words_sent,       # palabra
  HMM_tag    = v$tags,           # etiqueta asignada por el HMM
  udpipe_tag = upos_sent,        # etiqueta de referencia (udpipe)
  match      = v$tags == upos_sent,  # TRUE si coinciden, FALSE si no
  stringsAsFactors = FALSE
)

# Mostramos la comparación numérica
print(comparison)
cat("Proporción de aciertos del HMM en la oración de prueba:\n")
print(mean(comparison$match))

# VISUALIZAR GAMMA
plot(0,0, ylim=c(0,0.7), xlim=c(0,17))
for (i in 1:17) 
{
  lines(1:17, Gamma[i,], col=i)
}
