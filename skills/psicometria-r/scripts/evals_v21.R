# evals_v21.R, psicometria-r 2.1.0: fator de Bayes (t e correlação), DIF por regressão logística generalizada,
# Angoff modificado e Bookmark. Conferência com BayesFactor e difNLR quando instalados; sem eles, contra uma
# implementação independente (Monte Carlo, aproximação de Jeffreys, glm) e verdade simulada.
args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
suppressMessages(source(file.path(pasta, "carregar.R")))
resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome)); resultados[[nome]] <<- ok; invisible(ok)
}
set.seed(3); d <- data.frame(y = c(rnorm(30, 0.5), rnorm(30, 0)), g = rep(c("a", "b"), each = 30)); d$x2 <- d$y * .4 + rnorm(60)
r1 <- bf_teste_t(d, "y", "g")
set.seed(9); del <- stats::rcauchy(4e5, 0, sqrt(2) / 2)
mc <- mean(stats::dt(r1$t, r1$gl, ncp = del * sqrt(r1$N_efetivo))) / stats::dt(r1$t, r1$gl)
teste("BF do t independente (JZS) = Monte Carlo independente da verossimilhança marginal (erro relativo < 1%)", abs(r1$tabela$BF10[1] / mc - 1) < .01)
r0 <- bf_teste_t(d, "y", mu = 0); mc0 <- mean(stats::dt(r0$t, r0$gl, ncp = del * sqrt(r0$N_efetivo))) / stats::dt(r0$t, r0$gl)
teste("BF do t de uma amostra = Monte Carlo (erro relativo < 1%)", abs(r0$tabela$BF10[1] / mc0 - 1) < .01)
teste("BF do t: prior mais larga favorece H0 quando o efeito é pequeno (BF10 cai de r = .707 para r = sqrt(2))", r1$tabela$BF10[1] > r1$tabela$BF10[3])
rc <- bf_correlacao(d, "y", "x2"); r <- rc$r; n <- rc$n
L <- function(p) (1 - p^2)^((n - 1) / 2) / (1 - p * r)^(n - 1.5)
jf <- stats::integrate(function(p) 0.5 * L(p), -1, 1)$value / L(0)
teste("BF da correlação (prior uniforme) próximo da aproximação de Jeffreys (erro relativo < 1%)", abs(rc$tabela$BF10[1] / jf - 1) < .01)
if (requireNamespace("BayesFactor", quietly = TRUE)) {
  bt <- BayesFactor::extractBF(BayesFactor::ttestBF(formula = y ~ g, data = transform(d, g = factor(g)), rscale = "medium"))$bf
  bc <- BayesFactor::extractBF(BayesFactor::correlationBF(d$y, d$x2, rscale = 1))$bf
  teste("BayesFactor: ttestBF igual (tol. 1e-4 relativa)", abs(r1$tabela$BF10[1] / bt - 1) < 1e-4)
  teste("BayesFactor: correlationBF (rscale = 1, beta(1, 1) = uniforme) igual (tol. 1e-3 relativa)", abs(rc$tabela$BF10[1] / bc - 1) < 1e-3)
} else cat("PULADO  comparação com o BayesFactor (pacote ausente)\n")
# DIF
set.seed(7); n <- 1000; G <- rep(0:1, each = n / 2); th <- stats::rnorm(n); k <- 12
a <- stats::runif(k, .8, 2); b <- stats::rnorm(k, 0, .7)
Y <- sapply(1:k, function(i) { bb <- b[i] + (if (i == 3) .8 * G else 0); aa <- a[i] * (if (i == 5) 1 - .6 * G else 1)
  cc <- if (i == 7) .30 * G else 0; stats::rbinom(n, 1, cc + (1 - cc) * stats::plogis(aa * (th - bb))) })
colnames(Y) <- paste0("i", 1:k); D <- data.frame(Y, grupo = G)
r2 <- dif_nlr(D, colnames(Y), "grupo", 1, modelo = "2PL")
X <- as.numeric(scale(rowSums(Y)))
lr_glm <- sapply(colnames(Y), function(it) stats::deviance(stats::glm(D[[it]] ~ X, family = stats::binomial)) - stats::deviance(stats::glm(D[[it]] ~ X * G, family = stats::binomial)))
teste("DIF 2PL: LR igual ao da regressão logística com glm (tol. 1e-6)", max(abs(r2$tabela$LR - lr_glm)) < 1e-6)
teste("DIF 2PL: acha os 3 itens plantados (uniforme, não uniforme, acerto ao acaso) sem falso positivo", all(r2$tabela$dif[c(3, 5, 7)]) && !any(r2$tabela$dif[-c(3, 5, 7)]))
r3 <- dif_nlr(D[, c("i1", "i2", "i3", "i4", "i7", "grupo")], c("i1", "i2", "i3", "i4", "i7"), "grupo", 1, modelo = "3PLc")
teste("DIF 3PLc: 3 gl por item, c em [0, 1] e itens plantados detectados", all(r3$tabela$gl == 3) && all(r3$tabela$c >= 0 & r3$tabela$c <= 1) && all(r3$tabela$dif[r3$tabela$item %in% c("i3", "i7")]))
if (requireNamespace("difNLR", quietly = TRUE)) {
  fn <- difNLR::difNLR(Data = as.data.frame(Y), group = G, focal.name = 1, model = "2PL", method = "mle")
  teste("difNLR (method = 'mle', 2PL): mesmas estatísticas LR (tol. 1e-3)", max(abs(fn$Sval - r2$tabela$LR)) < 1e-3)
} else cat("PULADO  comparação com o difNLR (pacote ausente)\n")
# Standard setting
set.seed(2); J <- matrix(pmin(.95, pmax(.05, stats::rnorm(8 * 20, .6, .15))), 8, 20)
ra <- corte_angoff(J)
teste("Angoff: corte = média das somas por juiz e EP = DP/sqrt(J)", abs(ra$tabela$valor[1] - mean(rowSums(J))) < 1e-12 && abs(ra$tabela$valor[4] - stats::sd(rowSums(J)) / sqrt(8)) < 1e-12)
teste("Angoff: julgamento ausente ou fora de [0, 1] é recusado", inherits(tryCatch(corte_angoff(rbind(J[1, ], c(NA, J[2, -1]))), error = function(e) e), "error") &&
  inherits(tryCatch(corte_angoff(rbind(J[1, ], c(1.2, J[2, -1]))), error = function(e) e), "error"))
p <- data.frame(item = paste0("i", 1:10), a = seq(.8, 2, length = 10), b = seq(-2, 2, length = 10), c = c(rep(0, 5), rep(.2, 5)), d = 1)
rb <- corte_bookmark(p, marcadores = c(5, 6, 6, 7, 5), rp = .67); o <- rb$posthoc
teste("Bookmark: P(local) = RP em todos os itens e caderno em ordem crescente", all(abs(o$c + (o$d - o$c) * stats::plogis(o$a * (o$local - o$b)) - .67) < 1e-12) && !is.unsorted(o$local))
teste("Bookmark: corte = mediana dos thetas das páginas marcadas", abs(rb$tabela$valor[1] - stats::median(o$local[c(5, 6, 6, 7, 5)])) < 1e-12)
set.seed(5); thx <- stats::rnorm(600); Yb <- sapply(1:10, function(i) stats::rbinom(600, 1, stats::plogis(p$a[i] * (thx - p$b[i])))); colnames(Yb) <- p$item
fm <- mirt::mirt(as.data.frame(Yb), 1, "2PL", verbose = FALSE)
rm <- corte_bookmark(fm, c(4, 5, 5))
teste("Bookmark com ajuste do mirt: escore esperado no corte entre 0 e o número de itens", rm$tabela$valor[2] > 0 && rm$tabela$valor[2] < 10)
cat(sprintf("\nTotal evals_v21 2.1.0: %d | Passaram: %d | Falharam: %d\n", length(resultados), sum(unlist(resultados)), sum(!unlist(resultados))))
