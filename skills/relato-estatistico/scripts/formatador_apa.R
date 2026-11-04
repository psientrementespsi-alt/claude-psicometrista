# ==============================================================================
# formatador_apa.R  (relato-estatistico v1.3.0)
# Formatador unico de numeros com perfis por periodico.
# Regra: perfil de periodico so traz regra conferida na norma publicada; o que
# a norma nao diz fica herdado do perfil-base e marcado verificado = FALSE.
# ==============================================================================

.PERFIS_APA <- list(
  apa7_relato = list(
    descricao = "APA 7 na convencao desta skill (regra 2 do SKILL.md)",
    decimal = ".", sep = ", ", sem_zero = c("p", "r", "R2", "alfa", "omega", "cfi", "tli", "carga", "beta", "prop", "kappa", "icc"),
    verificado = TRUE, fonte = "SKILL.md regra 2; fabrica-de-artigo/references/statistical_reporting_standards.md"),
  apa7_estrito = list(
    descricao = "APA 7 secao 6.36 literal: sem zero so quando a estatistica nao pode passar de 1 (papaja::apa_num gt1 = FALSE)",
    decimal = ".", sep = ", ", sem_zero = c("p", "r", "R2", "alfa", "omega", "cfi", "carga", "prop", "kappa", "icc"),
    verificado = TRUE, fonte = "APA (2020) 6.36; papaja apa_num"),
  psicometria_r = list(
    descricao = "Convencao da psicometria-r (relato-apa.md e associacao.md)",
    decimal = ".", sep = ", ", sem_zero = c("p", "r", "R2", "alfa", "omega", "cfi", "tli", "rmsea", "srmr", "carga", "prop", "kappa", "icc"),
    verificado = TRUE, fonte = "psicometria-r/references/relato-apa.md; associacao.md (beta com zero)"),
  apa7_ptbr = list(
    descricao = "Manuscrito em portugues: virgula decimal com zero a esquerda (0,31), itens separados por ponto e virgula",
    decimal = ",", sep = "; ", sem_zero = character(0),
    verificado = FALSE, fonte = "fabrica-de-artigo/references/estatistica-apa.md (virgula em texto portugues); zero com virgula varia entre periodicos"),
  apa7_ptbr_sem_zero = list(
    base = "apa7_relato", descricao = "Virgula decimal sem zero nas grandezas limitadas a 1 (,31), como o decimal = ',' da psicometria-r",
    decimal = ",", sep = "; ",
    verificado = FALSE, fonte = "psicometria-r/references/relato-apa.md (zero_inicial = FALSE com virgula); conferir a norma do periodico"),
  rbtc = list(
    base = "apa7_ptbr", descricao = "Revista Brasileira de Terapias Cognitivas",
    url_norma = "https://pepsic.bvsalud.org/revistas/rbtc/pinstruc.htm",
    verificado = FALSE, fonte = "Instrucoes aos autores consultadas em 29/10/2026 sem regra de zero ou decimal; herda apa7_ptbr"),
  avaliacao_psicologica = list(
    base = "apa7_ptbr", descricao = "Avaliacao Psicologica (IBAP)",
    url_norma = "https://pepsic.bvsalud.org/revistas/avp/pinstruc.htm",
    verificado = FALSE, fonte = "Norma de numeros nao conferida em 29/10/2026; herda apa7_ptbr"),
  prc = list(
    base = "apa7_estrito", descricao = "Psicologia: Reflexao e Critica (SpringerOpen, em ingles)",
    url_norma = "https://prc.springeropen.com/submission-guidelines/preparing-your-manuscript",
    verificado = FALSE, fonte = "Diretrizes consultadas em 29/10/2026 sem regra de zero; herda apa7_estrito")
)

#' Lista os perfis disponiveis
listar_perfis_apa <- function() {
  do.call(rbind, lapply(names(.PERFIS_APA), function(n) {
    p <- perfil_apa(n)
    data.frame(perfil = n, descricao = p$descricao, decimal = p$decimal,
               verificado = p$verificado, url_norma = if (is.null(p$url_norma)) "" else p$url_norma,
               stringsAsFactors = FALSE)
  }))
}

#' Resolve um perfil (com heranca) ou devolve um perfil personalizado ja montado
perfil_apa <- function(perfil = "apa7_relato") {
  if (is.list(perfil)) return(perfil)
  if (!perfil %in% names(.PERFIS_APA))
    stop("Perfil desconhecido: ", perfil, ". Veja listar_perfis_apa() ou use perfil_apa_personalizado().")
  p <- .PERFIS_APA[[perfil]]
  if (!is.null(p$base)) {
    b <- perfil_apa(p$base)
    b[names(p)] <- p
    b$base <- NULL
    b$herdado_de <- p$base
    p <- b
  }
  p$nome <- perfil
  p
}

#' Monta perfil personalizado a partir de um perfil-base
perfil_apa_personalizado <- function(base = "apa7_relato", decimal = NULL, sep = NULL,
                                     sem_zero = NULL, adicionar_sem_zero = NULL,
                                     remover_sem_zero = NULL, fonte = "definido pelo autor") {
  p <- perfil_apa(base)
  if (!is.null(decimal)) p$decimal <- decimal
  if (!is.null(sep)) p$sep <- sep
  if (!is.null(sem_zero)) p$sem_zero <- sem_zero
  if (!is.null(adicionar_sem_zero)) p$sem_zero <- union(p$sem_zero, adicionar_sem_zero)
  if (!is.null(remover_sem_zero)) p$sem_zero <- setdiff(p$sem_zero, remover_sem_zero)
  p$nome <- "personalizado"; p$verificado <- FALSE; p$fonte <- fonte
  p
}

.DIGITOS_PADRAO <- c(p = 3, cfi = 3, tli = 3, rmsea = 3, srmr = 3, omega = 2, alfa = 2,
                     r = 2, R2 = 2, carga = 2, beta = 2, prop = 2, kappa = 2, icc = 2)

#' Formata numeros segundo o tipo de estatistica e o perfil
#'
#' @param x vetor numerico
#' @param tipo "p", "r", "R2", "alfa", "omega", "cfi", "tli", "rmsea", "srmr",
#'   "carga", "beta", "prop", "kappa", "icc" ou "livre" (M, DP, t, F, chi2, d, B...)
#' @param perfil nome do perfil ou lista de perfil_apa_personalizado()
#' @param digitos casas decimais (padrao: 3 para p e indices de ajuste, 2 nos demais)
#' @param menos "-" (ASCII, padrao) ou "−"
fmt_apa <- function(x, tipo = "livre", perfil = "apa7_relato", digitos = NULL, menos = "-") {
  pf <- perfil_apa(perfil)
  if (is.null(digitos)) digitos <- if (tipo %in% names(.DIGITOS_PADRAO)) .DIGITOS_PADRAO[[tipo]] else 2
  vapply(x, function(v) {
    if (is.na(v)) return("—")
    if (tipo == "p") {
      lim <- 10^(-digitos)
      if (v < lim) return(paste0("< ", .fmt_num(lim, digitos, TRUE, pf$decimal, menos)))
      if (v > 1 - lim) return(paste0("> ", .fmt_num(1 - lim, digitos, TRUE, pf$decimal, menos)))
    }
    .fmt_num(v, digitos, tipo %in% pf$sem_zero, pf$decimal, menos)
  }, character(1), USE.NAMES = FALSE)
}

.fmt_num <- function(v, digitos, sem_zero, decimal, menos) {
  txt <- formatC(abs(v), format = "f", digits = digitos)
  if (sem_zero) txt <- sub("^0\\.", ".", txt)
  if (decimal != ".") txt <- sub(".", decimal, txt, fixed = TRUE)
  neg <- v < 0 && as.numeric(formatC(abs(v), format = "f", digits = digitos)) != 0
  paste0(if (neg) menos else "", txt)
}

#' "p = .032" ou "p < .001"
fmt_p_apa <- function(p, perfil = "apa7_relato", digitos = 3, italico = TRUE) {
  s <- fmt_apa(p, "p", perfil, digitos)
  rot <- if (italico) "*p*" else "p"
  ifelse(grepl("^[<>]", s), paste(rot, s), paste(rot, "=", s))
}

#' Intervalo de confianca "[.12, .45]" ou "[0,12; 0,45]"
fmt_ic_apa <- function(li, ls, tipo = "livre", perfil = "apa7_relato", digitos = NULL) {
  pf <- perfil_apa(perfil)
  sep <- if (pf$decimal == ",") "; " else ", "
  paste0("[", fmt_apa(li, tipo, pf, digitos), sep, fmt_apa(ls, tipo, pf, digitos), "]")
}

#' Separador de itens de uma frase de resultado no perfil
sep_apa <- function(perfil = "apa7_relato") perfil_apa(perfil)$sep

#' Aviso de perfil nao verificado (para ir junto da entrega)
aviso_perfil_apa <- function(perfil = "apa7_relato") {
  pf <- perfil_apa(perfil)
  if (isTRUE(pf$verificado)) return(invisible(NULL))
  msg <- sprintf("Perfil '%s' nao verificado: %s%s", pf$nome, pf$fonte,
                 if (!is.null(pf$url_norma)) paste0(". Conferir em ", pf$url_norma) else "")
  message(msg)
  invisible(msg)
}
