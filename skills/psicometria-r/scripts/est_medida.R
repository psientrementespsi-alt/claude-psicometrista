# ==============================================================================
# est_medida.R  (psicometria-r 2.2.0)
# Erro-padrão condicional de medida (CSEM; Cook & Pitoniak, 2026, cap. 5,
# eq. 71-78, pp. 351-354) e equacionamento de escores observados entre formas
# (média, linear e equipercentil; Cook & Pitoniak, 2026, cap. 11, eq. 15-16,
# pp. 789-790; Kolen & Brennan, 2014, para a continuização por posto
# percentil). Em R base; os pacotes equate e kequate não estão instalados.
# ==============================================================================

#' Erro-padrão condicional de medida do escore bruto
#'
#' @param metodo "binomial" (Lord, eq. 78; itens 0/1), "thorndike" (DP da diferença entre metades por escore
#'   total, eq. 73) ou "tri" (1/sqrt(informação do teste), com `fit` do mirt)
#' @param metades lista com os itens de cada metade (thorndike); padrão ímpares contra pares
#' @param n_min escores (ou faixas) com menos casos que isto ficam sem CSEM no método thorndike
csem <- function(dados, itens, metodo = c("binomial", "thorndike", "tri"), metades = NULL, fit = NULL,
                 theta = seq(-3, 3, by = .5), n_min = 5, decimal = ".") {
  metodo <- match.arg(metodo)
  if (metodo == "tri") {
    if (is.null(fit)) stop("metodo = 'tri' exige fit (ajuste do mirt).", call. = FALSE)
    info <- mirt::testinfo(fit, matrix(theta))
    tab <- data.frame(theta = theta, informacao = info, csem = 1 / sqrt(info))
    return(.ee_resultado("CSEM pela TRI", n = NA, tabela = tab, decisoes = "CSEM(theta) = 1 / sqrt(I(theta)), na métrica de theta.",
      apa = sprintf("O erro-padrão condicional variou de %s a %s entre theta = %s e %s.", .pr_num(min(tab$csem), 2, TRUE, decimal),
                    .pr_num(max(tab$csem), 2, TRUE, decimal), min(theta), max(theta))))
  }
  X <- as.matrix(dados[stats::complete.cases(dados[, itens]), itens]); k <- ncol(X); tot <- rowSums(X)
  if (metodo == "binomial") {
    if (!all(X %in% c(0, 1))) stop("O modelo binomial de Lord exige itens 0/1.", call. = FALSE)
    x <- 0:k; cs <- sqrt(x * (k - x) / (k - 1)); fr <- tabulate(tot + 1, nbins = k + 1)
    tab <- data.frame(escore = x, frequencia = fr, csem = cs)
    sem_geral <- sqrt(sum(fr * cs^2) / sum(fr))
    dec <- c("CSEM de Lord (modelo binomial): sqrt(x(n - x)/(n - 1)) (Cook & Pitoniak, 2026, eq. 78), condicionado ao escore observado.",
             "O modelo trata os itens como amostra aleatória de itens de dificuldade igual; com dificuldades muito diferentes, superestima o erro (o capítulo cita o ajuste de Keats).",
             "SEM geral = raiz da média das variâncias condicionais ponderada pela frequência dos escores.")
  } else {
    if (is.null(metades)) metades <- list(itens[seq(1, k, 2)], itens[seq(2, k, 2)])
    dif <- rowSums(X[, metades[[1]], drop = FALSE]) - rowSums(X[, metades[[2]], drop = FALSE])
    g <- split(dif, tot); fr <- lengths(g)
    cs <- vapply(g, function(v) if (length(v) >= n_min) stats::sd(v) else NA_real_, 0)
    tab <- data.frame(escore = as.numeric(names(g)), frequencia = fr, csem = cs, row.names = NULL)
    sem_geral <- stats::sd(dif)
    dec <- c(sprintf("CSEM de Thorndike: DP da diferença entre metades entre quem tem o mesmo escore total (Cook & Pitoniak, 2026, eq. 73); escores com menos de %d casos sem estimativa.", n_min),
             "Supõe metades tau-equivalentes; a divisão padrão é ímpares contra pares.", "SEM geral = DP da diferença entre metades (eq. 72).")
  }
  apa <- sprintf("O erro-padrão condicional de medida (%s) variou de %s a %s pontos; o erro-padrão geral foi %s.",
                 if (metodo == "binomial") "modelo binomial de Lord" else "método de Thorndike",
                 .pr_num(min(tab$csem[tab$frequencia > 0], na.rm = TRUE), 2, TRUE, decimal), .pr_num(max(tab$csem[tab$frequencia > 0], na.rm = TRUE), 2, TRUE, decimal),
                 .pr_num(sem_geral, 2, TRUE, decimal))
  .ee_resultado(paste("CSEM", metodo), n = nrow(X), n_excluidos = nrow(dados) - nrow(X), tabela = tab, decisoes = dec, apa = apa,
                extra = list(sem_geral = sem_geral))
}

.eq_pr <- function(x, F, mn) {   # posto percentil (0-100) de x numa escala inteira mn..mx com F acumulada
  xs <- floor(x + .5); Fm1 <- if (xs - 1 < mn) 0 else F[xs - 1 - mn + 1]; Fx <- F[xs - mn + 1]
  100 * (Fm1 + (x - (xs - .5)) * (Fx - Fm1))
}

#' Equacionamento de escores observados (delineamento de grupos aleatórios ou grupo único)
#'
#' @param x escores na forma X (a converter); y escores na forma Y (escala de referência)
#' @param escala_x,escala_y valores possíveis (inteiros) de cada forma; padrão do mínimo ao máximo observado
#' @param metodo "media", "linear" (eq. 16) ou "equipercentil" (eq. 15, com continuização por posto percentil)
equacionar <- function(x, y, metodo = c("equipercentil", "linear", "media"), escala_x = NULL, escala_y = NULL, decimal = ".") {
  metodo <- match.arg(metodo); x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (is.null(escala_x)) escala_x <- seq(floor(min(x)), ceiling(max(x)))
  if (is.null(escala_y)) escala_y <- seq(floor(min(y)), ceiling(max(y)))
  if (any(x != round(x)) || any(y != round(y))) stop("equacionar() supõe escores inteiros.", call. = FALSE)
  mx <- mean(x); my <- mean(y); sx <- stats::sd(x); sy <- stats::sd(y)
  conv <- switch(metodo,
    media = escala_x - mx + my,
    linear = my + (escala_x - mx) * sy / sx,
    equipercentil = {
      fx <- tabulate(x - min(escala_x) + 1, length(escala_x)) / length(x); Fx <- cumsum(fx)
      fy <- tabulate(y - min(escala_y) + 1, length(escala_y)) / length(y); Fy <- cumsum(fy)
      vapply(escala_x, function(v) {
        P <- .eq_pr(v, Fx, min(escala_x)) / 100
        if (P <= 0) return(min(escala_y) - .5)
        if (P >= 1) return(max(escala_y) + .5)
        u <- which(Fy > P + 1e-12)[1]; if (is.na(u)) return(max(escala_y) + .5)
        yu <- escala_y[u]; Fprev <- if (u == 1) 0 else Fy[u - 1]
        (P - Fprev) / fy[u] + (yu - .5)
      }, 0)
    })
  tab <- data.frame(escore_x = escala_x, equivalente_y = conv)
  eqx <- conv[match(x, escala_x)]
  mom <- data.frame(distribuicao = c("Y observada", "X convertida"), media = c(my, mean(eqx)), dp = c(sy, stats::sd(eqx)))
  dec <- switch(metodo,
    media = "Equacionamento pela média: e(x) = x - média(X) + média(Y).",
    linear = "Equacionamento linear: e(x) = média(Y) + (x - média(X)) DP(Y)/DP(X) (Cook & Pitoniak, 2026, eq. 16).",
    equipercentil = "Equipercentil: e(x) = inversa da acumulada de Y no posto percentil de x, com continuização uniforme no intervalo de cada escore (Cook & Pitoniak, 2026, eq. 15; Kolen & Brennan, 2014). Sem suavização: com amostras pequenas, considerar pré ou pós-suavização.")
  .ee_resultado(paste("Equacionamento", metodo), n = length(x) + length(y), tabela = tab, posthoc = mom,
    decisoes = c(dec, "Delineamento de grupos aleatórios ou grupo único; com âncora (NEAT) a rota é outra e está fora do escopo testado."),
    apa = sprintf("Os escores da forma X foram convertidos para a escala da forma Y pelo método %s (N_X = %d, N_Y = %d); a média convertida foi %s (Y: %s).",
                  metodo, length(x), length(y), .pr_num(mean(eqx), 2, TRUE, decimal), .pr_num(my, 2, TRUE, decimal)))
}
