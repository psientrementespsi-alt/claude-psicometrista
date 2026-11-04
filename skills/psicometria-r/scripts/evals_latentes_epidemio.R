# evals_latentes_epidemio.R, testes de est_latentes.R (classes e perfis latentes) e de
# est_epidemio.R (razão de prevalência, sobrevida), v1.4.
#
# Rodar da pasta da skill:  Rscript scripts/evals_latentes_epidemio.R
# Dados simulados com semente fixa e o banco lung do pacote survival (público).

args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
source(file.path(pasta, "carregar.R"))
resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome))
  resultados[[nome]] <<- ok
  invisible(ok)
}
t_inicio <- Sys.time()
quieto <- function(expr) suppressWarnings(suppressMessages(expr))

cat("== Perfis latentes (mclust, modelo EEI = modelo 1 do tidyLPA)\n")
set.seed(4); n <- 600
cl <- sample(1:3, n, TRUE, prob = c(.5, .3, .2))
mu <- rbind(c(0, 0, 0, 0), c(1.5, 1.5, -1, -1), c(-1.5, 1, 1.5, -1.5))
dl <- as.data.frame(mu[cl, ] + matrix(stats::rnorm(n * 4), n)); names(dl) <- paste0("v", 1:4)
lp <- quieto(perfis_latentes(dl, names(dl), k_max = 4))
teste("LPA: BIC mínimo na solução de 3 perfis simulada", lp$tabela$perfis[which.min(lp$tabela$bic)] == 3)
teste("LPA: AIC = -2LL + 2p, BIC = -2LL + p ln n, SABIC = -2LL + p ln((n + 2)/24), p = (k - 1) + vk + v no EEI", {
  t <- lp$tabela; pp <- (t$perfis - 1) + 4 * t$perfis + 4
  all(t$parametros == pp) && all(abs(t$aic - (-2 * t$ll + 2 * pp)) < 1e-8) && all(abs(t$bic - (-2 * t$ll + pp * log(n))) < 1e-8) &&
    all(abs(t$sabic - (-2 * t$ll + pp * log((n + 2) / 24))) < 1e-8) })
if (requireNamespace("tidyLPA", quietly = TRUE)) {
  teste("LPA igual ao tidyLPA (estimate_profiles, models = 1): LL, BIC, entropia, prob_min, n_min (tol. 1e-6)", {
    g <- as.data.frame(quieto(tidyLPA::get_fit(tidyLPA::estimate_profiles(dl, n_profiles = 1:4, models = 1))))
    t <- lp$tabela
    max(abs(g$LogLik - t$ll)) < 1e-6 && max(abs(g$BIC - t$bic)) < 1e-6 && max(abs(g$Entropy - t$entropia)) < 1e-6 &&
      max(abs(g$prob_min - t$prob_min)) < 1e-6 && max(abs(g$n_min - t$menor_classe)) < 1e-6 })
} else cat("PULADO  comparação com tidyLPA (não instalado; as fórmulas acima cobrem)\n")
teste("LPA: entropia = 1 + Σ z ln z / (n ln K) na solução de 3 perfis (tol. 1e-8)", {
  z <- lp$modelos[[3]]$z; abs(lp$tabela$entropia[3] - (1 + sum(z * log(z + 1e-12)) / (n * log(3)))) < 1e-8 })

cat("== Classes latentes (poLCA)\n")
set.seed(5); n2 <- 800; c2 <- sample(1:2, n2, TRUE, prob = c(.6, .4))
pr <- rbind(c(.9, .8, .85, .2, .15), c(.2, .3, .25, .8, .9))
dc <- as.data.frame(sapply(1:5, function(j) stats::rbinom(n2, 1, pr[c2, j]))); names(dc) <- paste0("i", 1:5)
lc <- quieto(classes_latentes(dc, names(dc), k_max = 3, nrep = 5))
teste("LCA: BIC mínimo em 2 classes; categorias 0/1 recodificadas para 1/2 sem perder o rótulo nas probabilidades",
      lc$tabela$classes[which.min(lc$tabela$bic)] == 2 && all(c("0", "1") %in% names(lc$perfis[[1]])))
teste("LCA: AIC e BIC iguais aos do poLCA para 2 classes (tol. 1e-6)", {
  m <- lc$modelos[[2]]; abs(lc$tabela$aic[2] - m$aic) < 1e-6 && abs(lc$tabela$bic[2] - m$bic) < 1e-6 })
teste("LCA: LMR entre 1 e 2 classes pela fórmula do tidyLPA::calc_lrt (tol. 1e-8)", {
  t <- lc$tabela; lr <- 2 * (t$ll[2] - t$ll[1]); mod <- lr / (1 + ((3 * 2 - 1) - (3 * 1 - 1)) * log(n2))^-1 * 1
  mod <- lr / (1 + (((3 * 2 - 1) - (3 * 1 - 1)) * log(n2))^-1)
  abs(t$lmr[2] - mod) < 1e-8 })
teste("privacidade: impressão sem a classe de cada pessoa (saída curta)", length(capture.output(print(lc))) < 60)

cat("== Razão de prevalência (Poisson robusta) e sobrevida\n")
set.seed(12); n3 <- 800
x <- stats::rbinom(n3, 1, .4); g <- factor(sample(c("A", "B", "C"), n3, TRUE)); idade <- stats::rnorm(n3, 40, 10)
y <- stats::rbinom(n3, 1, pmin(.95, .25 * exp(.4 * x + .2 * (g == "B") + .01 * (idade - 40))))
dp <- data.frame(y = y, x = x, g = g, idade = idade)
rp <- razao_prevalencia(dp, "y", c("x", "g", "idade"), fatores = "x")
teste("RP: EP robusto = GEE independente com um caso por cluster (geepack) ou sanduíche HC0 à mão (tol. 1e-8)", {
  fit <- glm(y ~ factor(x) + g + idade, poisson, dp); X <- model.matrix(fit); mu <- fitted(fit); e <- dp$y - mu
  B <- solve(crossprod(X * sqrt(mu))); V <- B %*% crossprod(X * e) %*% B
  max(abs(sqrt(diag(V)) - rp$tabela$ep_robusto)) < 1e-8 })
teste("RP: estimativa pontual = exp(coeficiente da Poisson) e IC de Wald exponenciado", {
  i <- 2; abs(rp$tabela$rp[i] - exp(rp$tabela$b[i])) < 1e-12 &&
    abs(rp$tabela$ic_inf[i] - exp(rp$tabela$b[i] - qnorm(.975) * rp$tabela$ep_robusto[i])) < 1e-12 })
sv <- quieto(sobrevida(survival::lung, "time", "status", grupo = "sex", evento_codigo = 2, tempos = c(180, 365), fatores = "sex"))
teste("sobrevida (lung, survival): log-rank χ²(1) = 10,33 e HR de sex = 0,588 [0,424; 0,816]", {
  lr <- sv$tabelas_extra$`Log-rank`
  abs(lr$qui2 - 10.3267) < 1e-3 && abs(sv$tabela$hr - .5880) < 1e-3 && abs(sv$tabela$ic_inf - .4237) < 1e-3 && abs(sv$tabela$ic_sup - .8160) < 1e-3 })
teste("sobrevida: medianas de Kaplan-Meier 270 (homens) e 426 (mulheres), e tábua em 180 e 365 dias", {
  m <- sv$tabelas_extra$`Mediana de sobrevida (Kaplan-Meier)`; tb <- sv$tabelas_extra$`Tábua de sobrevida nos tempos pedidos`
  all(m$mediana == c(270, 426)) && nrow(tb) == 4 && all(tb$sobrevida > 0 & tb$sobrevida < 1) })
teste("sobrevida: riscos proporcionais conferidos (cox.zph) e frase com log-rank e HR", {
  !is.null(sv$tabelas_extra$`Riscos proporcionais (Schoenfeld)`) && any(grepl("log-rank", sv$apa)) && any(grepl("HR = ", sv$apa)) })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
print(versoes_pacotes(c("mclust", "poLCA", "tidyLPA", "survival", "geepack")))
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
