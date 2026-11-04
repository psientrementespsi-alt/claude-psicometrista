# evals_naoparametricos.R, testes de est_naoparametricos.R (v1.5): rankFD, nparLD,
# MANOVA.RM e as comparações do nparcomp.
#
# Rodar da pasta da skill:  Rscript scripts/evals_naoparametricos.R
# Dados simulados com semente fixa.

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
quieto <- function(expr) suppressWarnings(suppressMessages(expr))

cat("== np_fatorial (rankFD)\n")
set.seed(11); n <- 90
d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)), h = factor(rep(c("x", "y"), 45)))
d$y <- round(stats::rexp(n) + (d$g == "c") * 1.2 + (d$g == "b" & d$h == "y") * .8, 1)
r2 <- quieto(np_fatorial(d, "y", c("g", "h")))
ref <- quieto(rankFD::rankFD(y ~ g * h, data = d, hypothesis = "H0F"))
teste("ATS e WTS iguais à chamada direta do rankFD (tol. 1e-10)", {
  a <- r2$tabela[r2$tabela$estatistica == "ATS", ]; w <- r2$tabela[r2$tabela$estatistica == "WTS", ]
  max(abs(a$valor - ref$ANOVA.Type.Statistic[, 1])) < 1e-10 && max(abs(a$p - ref$ANOVA.Type.Statistic[, 4])) < 1e-10 &&
    max(abs(w$valor - ref$Wald.Type.Statistic[, 1])) < 1e-10 })
teste("efeito de g detectado (p < .001) e efeitos relativos de c acima de 0,5 com IC fora de 0,5", {
  a <- r2$tabela[r2$tabela$estatistica == "ATS", ]
  a$p[a$efeito == "g"] < .001 && all(r2$efeito$ic_inf[grepl("^c:", r2$efeito$celula)] > .5) })
d2 <- d[d$g %in% c("a", "c"), ]; d2$g <- droplevels(d2$g)
r1 <- quieto(np_fatorial(d2, "y", "g", posthoc = FALSE))
teste("dois grupos (gl1 = 1): ATS = WTS, forma quadrática de posto 1", {
  abs(r1$tabela$valor[1] - r1$tabela$valor[2]) < 1e-8 && r1$tabela$gl1[1] == 1 })
teste("um fator, n igual: efeito relativo = (posto médio - 0,5) / N", {
  rk <- rank(d2$y); pm <- tapply(rk, d2$g, mean)
  max(abs(r1$efeito$valor - (pm - .5) / nrow(d2))) < 1e-4 })  # o rankFD arredonda o efeito a 4 casas
r1p <- quieto(np_fatorial(d, "y", "g"))
teste("comparações (mctp, Fisher): diferença = diferença dos efeitos relativos e IC dentro de [-1, 1]", {
  ph <- r1p$posthoc; ef <- r1p$tabelas_extra[[1]]$efeito_relativo
  abs(ph$diferenca_efeito[ph$comparacao == "3 - 1"] - round(ef[3] - ef[1], 3)) < 1e-3 && all(ph$ic_inf >= -1 & ph$ic_sup <= 1) })
teste("interação significativa abre comparações por subamostra do segundo fator", {
  di <- d; di$y <- round(stats::rexp(n) + ifelse(di$h == "x", (di$g == "c") * 2, -(di$g == "c") * 2), 1)
  ri <- quieto(np_fatorial(di, "y", c("g", "h")))
  all(c("h = x", "h = y") %in% ri$posthoc$subamostra) })

cat("== np_longitudinal (nparLD)\n")
set.seed(12); dl <- expand.grid(id = 1:40, tempo = 1:3); dl$grupo <- factor(ifelse(dl$id <= 20, "C", "E"))
dl$sexo <- factor(ifelse(dl$id %% 2 == 0, "F", "M"))
dl$y <- round(stats::rnorm(nrow(dl)) + (dl$grupo == "E") * (dl$tempo - 1) * .9, 1)
l1 <- quieto(np_longitudinal(dl, "y", "id", "tempo", "grupo"))
dlf <- dl; dlf$tempo <- factor(dlf$tempo); dlf$id <- factor(dlf$id)
refl <- quieto(nparLD::nparLD(y ~ tempo * grupo, data = dlf, subject = "id"))
teste("F1-LD-F1: ATS e WTS iguais à chamada direta do nparLD (tol. 1e-10)", {
  a <- l1$tabela[l1$tabela$estatistica == "ATS", ]
  max(abs(sort(a$valor) - sort(refl$ATS[, 1]))) < 1e-10 && max(abs(sort(l1$tabela$valor[l1$tabela$estatistica == "WTS"]) - sort(refl$WTS[, 1]))) < 1e-10 })
teste("interação grupo × tempo simulada detectada e RTE do grupo E sobe com o tempo", {
  a <- l1$tabela[l1$tabela$estatistica == "ATS", ]
  rte <- l1$efeito$valor[grepl("^E:", l1$efeito$celula)]
  a$p[a$efeito == "grupo:tempo"] < .01 && all(diff(rte) > 0) })
l0 <- quieto(np_longitudinal(dl, "y", "id", "tempo"))
teste("LD-F1 (só tempo): desenho identificado e comparações entre tempos", grepl("LD-F1", l0$analise) && !grepl("F1-LD", l0$analise))
l2 <- quieto(np_longitudinal(dl, "y", "id", "tempo", c("grupo", "sexo")))
teste("F2-LD-F1: dois fatores entre sujeitos aceitos", grepl("F2-LD-F1", l2$analise))
teste("pessoa sem todos os tempos é excluída e contada", {
  dm <- dl[!(dl$id == 3 & dl$tempo == 2), ]
  lm_ <- quieto(np_longitudinal(dm, "y", "id", "tempo", "grupo")); lm_$n == 39 && lm_$n_excluidos == 1 })

cat("== manova_np (MANOVA.RM)\n")
set.seed(13); dm <- data.frame(g = factor(rep(c("a", "b", "c"), each = 25)))
dm$v1 <- stats::rnorm(75) + (dm$g == "c"); dm$v2 <- stats::rnorm(75, sd = 3); dm$v3 <- stats::rnorm(75)
m1 <- quieto(manova_np(dm, c("v1", "v2", "v3"), "g", iter = 2000))
refm <- quieto(MANOVA.RM::MANOVA.wide(cbind(v1, v2, v3) ~ g, data = dm, iter = 2000, seed = 2026))
teste("MANOVA.wide: WTS, MATS e p bootstrap iguais à chamada direta com a mesma semente", {
  abs(m1$tabela$WTS - refm$WTS[1, 1]) < 1e-8 && abs(m1$tabela$MATS - refm$MATS[1, 1]) < 1e-8 &&
    abs(m1$tabela$p_MATS_bootstrap - as.numeric(refm$resampling[1, 2])) < 1e-8 })
teste("MATS e WTS invariantes à escala das VDs (v2 × 100)", {
  dm2 <- dm; dm2$v2 <- dm2$v2 * 100
  m2 <- quieto(manova_np(dm2, c("v1", "v2", "v3"), "g", iter = 2000, univariados = FALSE))
  abs(m2$tabela$MATS - m1$tabela$MATS) < 1e-8 && abs(m2$tabela$WTS - m1$tabela$WTS) < 1e-8 })
teste("univariados com Holm: v1 (com efeito) tem o menor p", {
  u <- m1$tabelas_extra[[1]]; u$vd[which.min(u$p_MATS_bootstrap)] == "v1" && all(u$p_holm >= u$p_MATS_bootstrap) })
dl$y2 <- dl$y + stats::rnorm(nrow(dl))
mr <- quieto(manova_np(dl, c("y", "y2"), "grupo", id = "id", tempo = "tempo", iter = 1000, univariados = FALSE))
teste("multRM (medidas repetidas): três efeitos e interação detectada", nrow(mr$tabela) == 3 && mr$tabela$p_MATS_bootstrap[3] < .01)
teste("privacidade: nenhuma tabela impressa tem coluna de id", {
  out <- utils::capture.output(print(l1), print(m1), print(r2)); !any(grepl("\\bid\\b", out[grepl("^ *[a-z_]+ +[a-z_]", out)])) })

ok <- unlist(resultados)
cat(sprintf("\nevals_naoparametricos: %d/%d passaram\n", sum(ok), length(ok)))
if (!all(ok)) quit(status = 1)
