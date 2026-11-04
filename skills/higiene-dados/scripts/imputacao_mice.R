# ==============================================================================
# imputacao_mice.R
# Modulo de Diagnostico de Dados Ausentes e Imputacao Multipla MICE
# Skill: higiene-dados (v1.0)
# ==============================================================================

#' Diagnostico Estrutural de Dados Ausentes
#' 
#' @param df Data.frame com os dados a serem avaliados.
#' @param limiar_alerta_item Proporcao maxima toleravel de omissao por item (padrao 0.20).
#' @param limiar_descarte_caso Proporcao maxima toleravel de omissao por participante (padrao 0.20).
#' @return Lista contendo estatisticas agregadas, colunas criticas e linhas candidatas a descarte.
diagnosticar_ausentes <- function(df, limiar_alerta_item = 0.20, limiar_descarte_caso = 0.20) {
  n_total <- nrow(df)
  p_total <- ncol(df)
  
  # Omissao por coluna
  n_miss_col <- colSums(is.na(df))
  pct_miss_col <- round(n_miss_col / n_total, 4)
  
  resumo_colunas <- data.frame(
    variavel = colnames(df),
    n_ausente = n_miss_col,
    prop_ausente = pct_miss_col,
    alerta_omissao_alta = pct_miss_col > limiar_alerta_item,
    stringsAsFactors = FALSE
  )
  rownames(resumo_colunas) <- NULL
  
  # Omissao por linha (participante)
  n_miss_row <- rowSums(is.na(df))
  prop_miss_row <- round(n_miss_row / p_total, 4)
  descarte_sugerido <- prop_miss_row > limiar_descarte_caso
  
  # Padrões de completeza
  n_casos_completos <- sum(stats::complete.cases(df))
  pct_casos_completos <- round(n_casos_completos / n_total, 4)
  
  list(
    n_total_participantes = n_total,
    n_total_variaveis = p_total,
    casos_completos = n_casos_completos,
    pct_casos_completos = pct_casos_completos,
    resumo_colunas = resumo_colunas,
    linhas_descarte_sugerido = which(descarte_sugerido),
    total_linhas_descarte = sum(descarte_sugerido)
  )
}

#' Executa Imputacao Multipla por Equacoes Encadeadas (MICE) para Dados Psicometricos
#' 
#' @param df Data.frame contendo as variaveis a imputar e eventuais variaveis auxiliares.
#' @param m Numero de conjuntos imputados (padrao 20, conforme literatura moderna).
#' @param maxit Numero de iteracoes do algoritmo de Gibbs (padrao 10).
#' @param seed Semente de reproduzibilidade aleatoria.
#' @param metodo_padrao Metodo padrao de imputacao ("pmm" para predicao por pares, seguro para ordinal/continuo).
#' @return Objeto mids do pacote mice ou lista estruturada de bancos imputados.
executar_mice_psicometria <- function(df, m = 20, maxit = 10, seed = 42, metodo_padrao = "pmm") {
  
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("O pacote 'mice' e necessario para executar a imputacao multipla encadeada.")
  }
  
  set.seed(seed)
  
  # Configuracao inicial de metodos
  ini <- mice::mice(df, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  
  # Ajusta metodos por tipo de variavel
  for (col_name in colnames(df)) {
    if (ini$nmis[col_name] > 0) {
      if (is.factor(df[[col_name]])) {
        if (nlevels(df[[col_name]]) == 2) {
          meth[col_name] <- "logreg"
        } else if (is.ordered(df[[col_name]])) {
          meth[col_name] <- "polr"
        } else {
          meth[col_name] <- "polyreg"
        }
      } else if (is.numeric(df[[col_name]])) {
        # pmm preserva a escala original e valores discretos de itens Likert
        meth[col_name] <- metodo_padrao
      }
    } else {
      meth[col_name] <- ""
    }
  }
  
  # Executa o MICE
  imp <- mice::mice(
    data = df,
    m = m,
    maxit = maxit,
    method = meth,
    predictorMatrix = pred,
    seed = seed,
    printFlag = FALSE
  )
  
  imp
}

#' Teste de MCAR de Little (1988) para Dados Ausentes
#' 
#' Avalia a hipotese nula de que os dados sao Missing Completely at Random (MCAR).
#' Se o p-valor for > 0.05, nao se rejeita MCAR. Se p < 0.05, assume-se MAR ou MNAR.
#' 
#' @param df Data.frame com dados numericos continuos ou ordinais.
#' @return Lista contendo estatistica qui-quadrado, graus de liberdade e p-valor.
teste_little_mcar <- function(df, max_iter = 1000, tol = 1e-8) {
  # v1.0.1: mu e Sigma por EM e todos os padroes de ausencia na estatistica, inclusive os de um so caso.
  # A versao 1.0 usava medias e covariancias pairwise e pulava padroes com n_j < 2, o que zerava o teste
  # quando cada omissao era de um participante diferente (chi2 = 0, gl = 1, p = 1).
  df_num <- df[, vapply(df, is.numeric, logical(1)), drop = FALSE]
  df_num <- df_num[rowSums(!is.na(df_num)) > 0, , drop = FALSE]
  if (!anyNA(df_num)) return(list(qui_quadrado = 0, df = 0, p_valor = 1, metodo = "sem ausentes",
                                  conclusao = "Nenhum dado ausente nas variaveis numericas (MCAR trivial)."))
  conclui <- function(p) if (p >= 0.05) "Nao rejeita MCAR (p >= .05)" else "Rejeita MCAR (p < .05, dados provavelmente MAR/MNAR)"
  if (requireNamespace("naniar", quietly = TRUE)) {
    r <- tryCatch(suppressWarnings(naniar::mcar_test(df_num)), error = function(e) NULL)
    if (!is.null(r) && !is.na(r$statistic)) return(list(qui_quadrado = round(r$statistic, 3), df = r$df, p_valor = round(r$p.value, 4),
      padroes = r$missing.patterns, metodo = "naniar::mcar_test (Little, 1988)", conclusao = conclui(r$p.value)))
  }
  X <- as.matrix(df_num); n <- nrow(X); p <- ncol(X)
  mu <- colMeans(X, na.rm = TRUE); S <- diag(apply(X, 2, stats::var, na.rm = TRUE), p); convergiu <- FALSE
  for (it in seq_len(max_iter)) {
    T1 <- numeric(p); T2 <- matrix(0, p, p)
    for (i in seq_len(n)) { x <- X[i, ]; o <- !is.na(x); m <- !o; xf <- x; Cm <- matrix(0, p, p)
      if (any(m)) { B <- S[m, o, drop = FALSE] %*% solve(S[o, o, drop = FALSE])
        xf[m] <- mu[m] + B %*% (x[o] - mu[o]); Cm[m, m] <- S[m, m, drop = FALSE] - B %*% S[o, m, drop = FALSE] }
      T1 <- T1 + xf; T2 <- T2 + tcrossprod(xf) + Cm }
    mu1 <- T1 / n; S1 <- T2 / n - tcrossprod(mu1)
    dif <- max(abs(mu1 - mu), abs(S1 - S)); mu <- mu1; S <- S1
    if (dif < tol) { convergiu <- TRUE; break }
  }
  if (!convergiu) warning("EM do teste de Little nao convergiu em ", max_iter, " iteracoes.", call. = FALSE)
  pad <- apply(!is.na(X), 1, paste, collapse = ""); d2 <- 0; gl <- 0
  for (u in unique(pad)) { r <- which(pad == u); o <- !is.na(X[r[1], ])
    db <- colMeans(X[r, o, drop = FALSE]) - mu[o]
    d2 <- d2 + length(r) * drop(crossprod(db, solve(S[o, o, drop = FALSE], db))); gl <- gl + sum(o) }
  gl <- gl - p
  if (gl < 1) return(list(qui_quadrado = NA_real_, df = gl, p_valor = NA_real_, padroes = length(unique(pad)),
                          metodo = "Little (1988), EM em base R", conclusao = "Teste nao calculavel (gl < 1)."))
  pv <- stats::pchisq(d2, gl, lower.tail = FALSE)
  list(qui_quadrado = round(d2, 3), df = gl, p_valor = round(pv, 4), padroes = length(unique(pad)), convergiu = convergiu,
       metodo = "Little (1988), mu e Sigma por EM em base R", conclusao = conclui(pv))
}

#' Exporta Imputacoes MICE para Lista de Data.Frames Compativel com lavaan.mi / semTools
#'
#' Extrai cada um dos m bancos imputados em um objeto mids do pacote mice e empacota em
#' uma lista de data.frames pronta para semTools::runMI() ou lavaan.mi (Rubin, 1987; Enders, 2022).
#'
#' @param imp Objeto da classe mids gerado por mice::mice() ou executar_mice_psicometria().
#' @param cols Vetor opcional com nomes das colunas de interesse para reter.
#' @return Lista de data.frames completos com atributo 'm' e classe 'lavaan_mi_datalist'.
exportar_lavaan_mi <- function(imp, cols = NULL) {
  if (!inherits(imp, "mids")) {
    stop("imp deve ser um objeto da classe 'mids' gerado pelo pacote mice.")
  }
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("O pacote 'mice' e necessario para extrair os bancos imputados.")
  }
  m <- imp$m
  lista_df <- vector("list", m)
  for (i in seq_len(m)) {
    d <- mice::complete(imp, action = i)
    if (!is.null(cols)) {
      d <- d[, intersect(cols, colnames(d)), drop = FALSE]
    }
    lista_df[[i]] <- d
  }
  attr(lista_df, "m") <- m
  attr(lista_df, "classe_origem") <- "mids"
  class(lista_df) <- c("lavaan_mi_datalist", "list")
  lista_df
}
