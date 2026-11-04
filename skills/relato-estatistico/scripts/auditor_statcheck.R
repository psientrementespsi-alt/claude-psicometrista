# ==============================================================================
# auditor_statcheck.R
# Auditor de consistencia de testes NHST (relato-estatistico v1.2.0)
# Regra de decisao portada do pacote statcheck (Nuijten et al., 2016; codigo em
# repositorios-relato-estatistico/statcheck/R: compute_p.R, error_test.R,
# decision_error_test.R), com leitura de texto em portugues (virgula decimal,
# chi-quadrado em Unicode, sinal de menos tipografico, r(gl) e z).
# A v1.0 comparava |p_declarado - p_calculado| com tolerancia fixa de ,015, nao lia
# r(gl) nem virgula decimal; agora o p e recalculado nos dois extremos do intervalo de
# arredondamento da estatistica relatada e comparado nas casas decimais do p relatado.
# ==============================================================================

.sc_compute_p <- function(tipo, est, df1, df2, bicaudal = TRUE) {
  p <- switch(tipo,
    t = stats::pt(-abs(est), df2),
    F = stats::pf(est, df1, df2, lower.tail = FALSE),
    z = stats::pnorm(abs(est), lower.tail = FALSE),
    r = { tt <- est / sqrt((1 - est^2) / df2); stats::pt(-abs(tt), df2) },
    chi2 = stats::pchisq(est, df1, lower.tail = FALSE))
  if (tipo %in% c("t", "z", "r") && isTRUE(bicaudal)) p <- 2 * p
  p
}

# error_test() do statcheck, sem as comparacoes do proprio teste (o auditor so le "=")
.sc_erro <- function(p_rel, tipo, est, df1, df2, p_comp, p_dec, est_dec, bicaudal = TRUE) {
  meia <- 0.5 / 10^est_dec
  baixo <- if (est >= 0) est - meia else est + meia
  alto <- if (est >= 0) est + meia else est - meia
  if (tipo == "r") { baixo <- max(min(baixo, 1), -1); alto <- max(min(alto, 1), -1) }
  p_alto <- .sc_compute_p(tipo, baixo, df1, df2, bicaudal)
  p_baixo <- .sc_compute_p(tipo, alto, df1, df2, bicaudal)
  if (p_rel <= 0) return(TRUE)
  switch(p_comp,
    "=" = p_rel > round(p_alto, p_dec) | p_rel < round(p_baixo, p_dec),
    "<" = p_rel < p_baixo,
    ">" = p_rel > p_alto)
}

# decision_error_test() do statcheck com pEqualAlphaSig = TRUE (padrao do pacote)
.sc_decisao <- function(p_rel, p_calc, p_comp, alfa) {
  switch(p_comp,
    "=" = (p_rel <= alfa & p_calc > alfa) | (p_rel > alfa & p_calc <= alfa),
    "<" = p_rel <= alfa & p_calc > alfa,
    ">" = p_rel >= alfa & p_calc <= alfa)
}

.sc_normalizar <- function(x) {
  x <- gsub("\\\\chi\\^?\\{?2\\}?|χ²|χ\\^?2|χ2|chi\\^2|chi²|qui2", "chi2", x)
  x <- gsub("[−–—]", "-", x)
  x <- gsub("[*$_{}]", "", x)
  x <- gsub("\\\\text|\\\\mathrm|\\\\left|\\\\right", "", x)
  x <- gsub("(?<=[\\d\\s\\[(<=>;-]),(?=\\d)", ".", x, perl = TRUE)   # virgula decimal
  gsub("\\s+", " ", x)
}

#' Audita testes NHST em texto (t, F, chi2, r, z) pela regra do statcheck
#'
#' @param texto String com o texto; ou `arquivo` (.md/.txt).
#' @param nivel_alfa Alfa da decisao (padrao .05).
#' @param bicaudal Testes t, r e z bicaudais (padrao TRUE, como o statcheck).
#' @return data.frame com uma linha por teste: teste, expressao, estatistica, df1, df2, operador, p_declarado,
#'   p_recalculado, p_min, p_max, inconsistente, erro_grosseiro.
auditar_texto_statcheck <- function(texto = NULL, arquivo = NULL, nivel_alfa = 0.05, bicaudal = TRUE) {
  if (is.null(texto) && !is.null(arquivo)) {
    if (!file.exists(arquivo)) stop("Arquivo nao encontrado: ", arquivo)
    texto <- paste(readLines(arquivo, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  }
  if (is.null(texto) || !nzchar(texto)) stop("Forneca texto ou arquivo.")
  x <- .sc_normalizar(texto)
  num <- "-?\\d*\\.?\\d+"
  pp <- "\\s?[,;]\\s?p\\s?([<>=])\\s?(\\d?\\.\\d+)"   # ; separa itens em texto com virgula decimal
  padroes <- list(
    t = paste0("(?<![A-Za-z])t\\s?\\(\\s?(\\d+\\.?\\d*)\\s?\\)\\s?=\\s?(", num, ")", pp),
    r = paste0("(?<![A-Za-z])r\\s?\\(\\s?(\\d+)\\s?\\)\\s?=\\s?(", num, ")", pp),
    F = paste0("(?<![A-Za-z])F\\s?\\(\\s?(\\d+\\.?\\d*)\\s?[,;]\\s?(\\d+\\.?\\d*)\\s?\\)\\s?=\\s?(", num, ")", pp),
    chi2 = paste0("chi2\\s?\\(\\s?(\\d+\\.?\\d*)\\s?(?:,\\s?N\\s?=\\s?\\d+)?\\s?\\)\\s?=\\s?(", num, ")", pp),
    z = paste0("(?<![A-Za-z(])z\\s?=\\s?(", num, ")", pp))
  res <- list()
  for (tipo in names(padroes)) {
    m <- regmatches(x, gregexpr(padroes[[tipo]], x, perl = TRUE))[[1]]
    for (e in m) {
      g <- regmatches(e, regexec(padroes[[tipo]], e, perl = TRUE))[[1]]
      if (tipo == "F") { df1 <- as.numeric(g[2]); df2 <- as.numeric(g[3]); v <- g[4]; comp <- g[5]; pr <- g[6] }
      else if (tipo == "z") { df1 <- NA; df2 <- NA; v <- g[2]; comp <- g[3]; pr <- g[4] }
      else if (tipo == "chi2") { df1 <- as.numeric(g[2]); df2 <- NA; v <- g[3]; comp <- g[4]; pr <- g[5] }
      else { df1 <- NA; df2 <- as.numeric(g[2]); v <- g[3]; comp <- g[4]; pr <- g[5] }
      est <- as.numeric(v); p_rel <- as.numeric(pr)
      est_dec <- max(0, nchar(sub("^[^.]*\\.?", "", v)))
      p_dec <- max(0, nchar(sub("^[^.]*\\.?", "", pr)))
      p_calc <- .sc_compute_p(tipo, est, df1, df2, bicaudal)
      meia <- 0.5 / 10^est_dec
      lim <- sort(c(.sc_compute_p(tipo, est - meia, df1, df2, bicaudal), .sc_compute_p(tipo, est + meia, df1, df2, bicaudal)))
      erro <- .sc_erro(p_rel, tipo, est, df1, df2, comp, p_dec, est_dec, bicaudal)
      dec <- erro && .sc_decisao(p_rel, p_calc, comp, nivel_alfa)
      res[[length(res) + 1]] <- data.frame(teste = tipo, expressao = e, estatistica = est, df1 = df1, df2 = df2,
        operador = comp, p_declarado = p_rel, p_recalculado = round(p_calc, 5), p_min = signif(lim[1], 4),
        p_max = signif(lim[2], 4), inconsistente = isTRUE(erro), erro_grosseiro = isTRUE(dec), stringsAsFactors = FALSE)
    }
  }
  if (!length(res)) { cat("Nenhuma expressao NHST (t, F, chi2, r, z) com p detectada.\n"); return(invisible(data.frame())) }
  df_res <- do.call(rbind, res)
  cat(sprintf("AUDITORIA STATCHECK (regra do pacote statcheck): %d testes | inconsistencias: %d | erros grosseiros: %d\n",
              nrow(df_res), sum(df_res$inconsistente), sum(df_res$erro_grosseiro)))
  if (any(df_res$inconsistente)) print(df_res[df_res$inconsistente, c("teste", "expressao", "p_declarado", "p_min", "p_max", "erro_grosseiro")], row.names = FALSE)
  df_res
}
