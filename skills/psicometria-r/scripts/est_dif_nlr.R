# ==============================================================================
# est_dif_nlr.R  (psicometria-r 2.1.0)
# DIF por regressão logística generalizada (modelo 9.9 de Martinková & Hladká,
# 2023, p. 259; Drabinová & Martinková, 2017): assíntotas inferior (c) e
# superior (d) e diferenças de grupo em a, b, c e d, por máxima verossimilhança
# em R base, com teste da razão de verossimilhança. É a rota do difNLR, que não
# está instalado no Mac; com modelo = "2PL" o LRT é idêntico ao da regressão
# logística com glm (conferido nos evals).
# ==============================================================================

.nlr_p <- function(X, G, a, b, c, d, aD, bD, cD, dD) {
  cc <- c + cD * G; dd <- d + dD * G
  cc + (dd - cc) * stats::plogis((a + aD * G) * (X - b - bD * G))
}

.nlr_ajuste <- function(y, X, G, livres, inicios) {
  # livres: nomes dos parâmetros estimados entre a, b, c, d, aD, bD, cD, dD; os outros ficam fixos
  # c e d (e c + cD, d + dD) são estimados na escala logit para ficarem em [0, 1]
  fixos <- c(a = 1, b = 0, c = 0, d = 1, aD = 0, bD = 0, cD = 0, dD = 0)
  monta <- function(th) {
    p <- fixos; p[names(th)] <- th
    if ("c" %in% names(th)) p["c"] <- stats::plogis(th["c"])
    if ("d" %in% names(th)) p["d"] <- stats::plogis(th["d"])
    if ("cD" %in% names(th)) p["cD"] <- stats::plogis(th["cD"]) - p["c"]
    if ("dD" %in% names(th)) p["dD"] <- stats::plogis(th["dD"]) - p["d"]
    p
  }
  nll <- function(th) {
    names(th) <- livres; p <- monta(th)
    pr <- .nlr_p(X, G, p["a"], p["b"], p["c"], p["d"], p["aD"], p["bD"], p["cD"], p["dD"])
    pr <- pmin(pmax(pr, 1e-10), 1 - 1e-10)
    -sum(y * log(pr) + (1 - y) * log(1 - pr))
  }
  melhor <- NULL
  for (ini in inicios) {
    o <- tryCatch(stats::optim(ini[livres], nll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)), error = function(e) NULL)
    if (!is.null(o) && (is.null(melhor) || o$value < melhor$value)) melhor <- o
  }
  if (is.null(melhor)) return(list(logLik = NA_real_, par = rep(NA_real_, 8), convergiu = FALSE))
  th <- melhor$par; names(th) <- livres
  list(logLik = -melhor$value, par = monta(th), convergiu = melhor$convergence == 0)
}

#' DIF por regressão logística generalizada (rota do difNLR)
#'
#' @param dados data.frame com itens 0/1 e a variável de grupo
#' @param itens nomes dos itens dicotômicos
#' @param grupo variável de grupo (2 níveis)
#' @param focal nível focal (codificado 1)
#' @param criterio "escore_z" (padrão do difNLR: escore total padronizado) ou "escore" (total bruto)
#' @param modelo "2PL", "3PLc" (acerto ao acaso), "3PLd" (desatenção) ou "4PL"
#' @param tipo "todos" (todas as diferenças de grupo, padrão do difNLR), "uniforme" (só bD; e cD/dD se houver) ou
#'   "nao_uniforme" (aD, com bD no modelo nulo)
#' @param ajuste_p correção para múltiplos itens ("none", padrão do difNLR, "holm", "BH")
dif_nlr <- function(dados, itens, grupo, focal, criterio = c("escore_z", "escore"), modelo = c("2PL", "3PLc", "3PLd", "4PL"),
                    tipo = c("todos", "uniforme", "nao_uniforme"), ajuste_p = c("none", "holm", "BH"), alfa = .05) {
  criterio <- match.arg(criterio); modelo <- match.arg(modelo); tipo <- match.arg(tipo); ajuste_p <- match.arg(ajuste_p)
  d <- dados[stats::complete.cases(dados[, c(itens, grupo)]), c(itens, grupo)]
  Y <- as.matrix(d[, itens]); if (!all(Y %in% c(0, 1))) stop("dif_nlr() aceita só itens dicotômicos 0/1.", call. = FALSE)
  G <- as.integer(as.character(d[[grupo]]) == as.character(focal)); if (length(unique(G)) != 2) stop("focal não encontrado ou grupo com 1 nível.", call. = FALSE)
  tot <- rowSums(Y); X <- if (criterio == "escore_z") (tot - mean(tot)) / stats::sd(tot) else tot
  base <- c("a", "b"); if (modelo %in% c("3PLc", "4PL")) base <- c(base, "c"); if (modelo %in% c("3PLd", "4PL")) base <- c(base, "d")
  difs <- c(aD = "a", bD = "b", cD = "c", dD = "d"); difs <- names(difs)[difs %in% base]
  testados <- switch(tipo, todos = difs, uniforme = setdiff(difs, "aD"), nao_uniforme = "aD")
  nulo_livres <- c(base, setdiff(difs, testados)); cheio_livres <- c(base, difs)
  res <- lapply(itens, function(it) {
    y <- d[[it]]; g0 <- stats::glm(y ~ X, family = stats::binomial)
    a0 <- unname(stats::coef(g0)[2]); b0 <- unname(-stats::coef(g0)[1] / stats::coef(g0)[2])
    if (!is.finite(a0) || abs(a0) < 1e-3) { a0 <- 1; b0 <- 0 }
    inis <- lapply(list(c(.05, .95), c(.20, .90), c(.01, .99)), function(cd)
      c(a = a0, b = b0, c = stats::qlogis(cd[1]), d = stats::qlogis(cd[2]), aD = 0, bD = 0,
        cD = stats::qlogis(cd[1]), dD = stats::qlogis(cd[2])))
    f0 <- .nlr_ajuste(y, X, G, nulo_livres, inis); f1 <- .nlr_ajuste(y, X, G, cheio_livres, inis)
    LR <- max(0, 2 * (f1$logLik - f0$logLik)); gl <- length(testados)
    data.frame(item = it, LR = LR, gl = gl, p = stats::pchisq(LR, gl, lower.tail = FALSE),
               a = f1$par["a"], b = f1$par["b"], c = f1$par["c"], d = f1$par["d"], aDif = f1$par["aD"], bDif = f1$par["bD"],
               cDif = f1$par["cD"], dDif = f1$par["dD"], convergiu = f0$convergiu && f1$convergiu, row.names = NULL)
  })
  tab <- do.call(rbind, res); tab$p_ajust <- stats::p.adjust(tab$p, ajuste_p); tab$dif <- tab$p_ajust < alfa
  avisos <- if (any(!tab$convergiu)) paste("Sem convergência em:", paste(tab$item[!tab$convergiu], collapse = ", ")) else character()
  dec <- c(sprintf("Modelo %s (Martinková & Hladká, 2023, eq. 9.9); critério de pareamento: %s.", modelo,
                   if (criterio == "escore_z") "escore total padronizado (padrão do difNLR)" else "escore total"),
           sprintf("Teste da razão de verossimilhança das diferenças de grupo (%s), %d gl por item; correção de p: %s.",
                   paste(testados, collapse = ", "), length(testados), ajuste_p),
           "Estimação por máxima verossimilhança com três pontos de partida; c e d restritos a [0, 1]. O difNLR usa por padrão mínimos quadrados não lineares, então os valores podem diferir na terceira casa.",
           "O escore total inclui o item testado e não há purificação; com muitos itens com DIF, refazer sem eles no critério.")
  apa <- sprintf("O DIF por regressão logística generalizada (modelo %s) apontou %d de %d itens (%s).", modelo, sum(tab$dif), nrow(tab),
                 if (any(tab$dif)) paste(sprintf("%s: LR(%d) = %.2f, p %s", tab$item[tab$dif], tab$gl[tab$dif], tab$LR[tab$dif],
                   ifelse(tab$p_ajust[tab$dif] < .001, "< .001", paste("=", sub("^0", "", sprintf("%.3f", tab$p_ajust[tab$dif]))))), collapse = "; ") else "nenhum")
  .ee_resultado("DIF por regressão logística generalizada", n = nrow(d), n_excluidos = nrow(dados) - nrow(d), tabela = tab,
                decisoes = dec, avisos = avisos, apa = apa)
}
