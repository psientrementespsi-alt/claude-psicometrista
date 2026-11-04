# ==============================================================================
# checagem_efeitos.R  (relato-estatistico v1.5.0)
# Confere se cada teste relatado traz tamanho de efeito e IC (APA 7; JARS-Quant;
# Lakens, 2013) e reconstitui d e eta2p a partir da estatistica relatada.
# ==============================================================================

.PAD_TESTE <- "(\\bt\\s*\\(|\\bF\\s*\\(|χ²?\\s*\\(|χ2\\s*\\(|chi2\\s*\\(|\\bz\\s*=|\\br\\s*\\()"
.PAD_EFEITO <- "(\\b[dg]\\s*=|η²?p?\\s*=|ηp²\\s*=|ω²?\\s*=|R²\\s*=|R2\\s*=|\\br\\s*[=(]|β\\s*=|\\bOR\\s*=|\\bRR\\s*=|\\bV\\s*=|φ\\s*=|f²\\s*=|\\bB\\s*=)"
.PAD_IC <- "(\\b(IC|CI)\\s*\\d{2}\\s*%|\\[[^]]*[,;][^]]*\\])"

#' Uma linha por frase com teste: tem tamanho de efeito? tem IC?
auditar_efeitos_texto <- function(texto) {
  x <- gsub("\\*", "", paste(texto, collapse = " "))
  frases <- unlist(strsplit(x, "(?<=[.!?])\\s+(?=[A-ZÀ-Ü])", perl = TRUE))
  idx <- which(grepl(.PAD_TESTE, frases, perl = TRUE))
  if (!length(idx)) return(data.frame(frase = integer(0), trecho = character(0), tem_efeito = logical(0), tem_ic = logical(0)))
  data.frame(
    frase = idx,
    trecho = substr(frases[idx], 1, 90),
    tem_efeito = grepl(.PAD_EFEITO, frases[idx], perl = TRUE),
    tem_ic = grepl(.PAD_IC, frases[idx], perl = TRUE),
    stringsAsFactors = FALSE)
}

.intervalo_arred <- function(v, casas) c(v - 0.5 * 10^(-casas), v + 0.5 * 10^(-casas))
.casas <- function(txt) { txt <- as.character(txt); if (!grepl("[.,]", txt)) 0L else nchar(sub("^.*[.,]", "", txt)) }

#' Reconstitui d de Cohen a partir de t e confere o d relatado (Lakens, 2013, eq. 2 e 5)
#'
#' @param t_txt t como relatado (texto, para saber as casas)
#' @param n1,n2 tamanhos dos grupos; n2 = NULL para d_z de amostras pareadas (n1 = n de pares)
#' @param d_txt d relatado (texto)
#' @param hedges TRUE confere g_s de Hedges: d_s * (1 - 3 / (4 * (n1 + n2) - 9)) (Lakens, 2013, eq. 4)
conferir_d_de_t <- function(t_txt, n1, n2 = NULL, d_txt, hedges = FALSE) {
  tt <- .intervalo_arred(as.numeric(sub(",", ".", t_txt)), .casas(t_txt))
  k <- if (is.null(n2)) 1 / sqrt(n1) else sqrt(1 / n1 + 1 / n2)
  d_lim <- sort(abs(tt) * k)
  if (isTRUE(hedges)) { if (is.null(n2)) stop("conferir_d_de_t: hedges = TRUE so para grupos independentes (g_s).", call. = FALSE)
    d_lim <- d_lim * (1 - 3 / (4 * (n1 + n2) - 9)) }
  d_rel <- abs(as.numeric(sub(",", ".", d_txt)))
  dd <- .intervalo_arred(d_rel, .casas(d_txt))
  data.frame(d_relatado = d_rel, d_min = round(d_lim[1], 4), d_max = round(d_lim[2], 4),
             formula = if (is.null(n2)) "d_z = t / sqrt(n)" else if (isTRUE(hedges)) "g_s = t * sqrt(1/n1 + 1/n2) * (1 - 3/(4(n1 + n2) - 9))" else "d_s = t * sqrt(1/n1 + 1/n2)",
             consistente = dd[2] >= d_lim[1] && dd[1] <= d_lim[2])
}

#' Reconstitui eta2 parcial a partir de F(df1, df2) e confere o relatado (Lakens, 2013, eq. 13)
conferir_eta2p_de_F <- function(F_txt, df1, df2, eta_txt) {
  ff <- .intervalo_arred(as.numeric(sub(",", ".", F_txt)), .casas(F_txt))
  e <- ff * df1 / (ff * df1 + df2)
  ee <- .intervalo_arred(as.numeric(sub(",", ".", eta_txt)), .casas(eta_txt))
  data.frame(eta2p_relatado = as.numeric(sub(",", ".", eta_txt)), eta_min = round(e[1], 4), eta_max = round(e[2], 4),
             consistente = ee[2] >= e[1] && ee[1] <= e[2])
}

#' Reconstitui omega2 parcial a partir de F(df1, df2) e N e confere o relatado
#' (Lakens, 2013, eq. 15, com F = MS_efeito / MS_erro: df1 (F - 1) / (df1 F + N - df1))
#' @param N numero total de participantes
conferir_omega2p_de_F <- function(F_txt, df1, N, omega_txt) {
  ff <- .intervalo_arred(as.numeric(sub(",", ".", F_txt)), .casas(F_txt))
  w <- df1 * (ff - 1) / (df1 * ff + N - df1)
  ww <- .intervalo_arred(as.numeric(sub(",", ".", omega_txt)), .casas(omega_txt))
  data.frame(omega2p_relatado = as.numeric(sub(",", ".", omega_txt)), omega_min = round(w[1], 4), omega_max = round(w[2], 4),
             consistente = ww[2] >= w[1] && ww[1] <= w[2])
}
