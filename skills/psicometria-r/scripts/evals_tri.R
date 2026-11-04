# evals_tri.R, testes do módulo de TRI (est_tri.R, v1.4): tri_ajustar(), tri_comparar(),
# tri_dif().
#
# Rodar da pasta da skill:  Rscript scripts/evals_tri.R
# Dados simulados com parâmetros conhecidos e semente fixa; nenhum dado de participante.
# Conferências: recuperação de parâmetros, mirt::M2, fórmula da confiabilidade marginal,
# difR::difLogistic (regressão logística) e eRm (confiabilidade de separação), quando
# instalados.

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

cat("== 2PL dicotômico (parâmetros conhecidos)\n")
set.seed(7); n <- 1000; th <- stats::rnorm(n)
a <- c(1.2, .8, 1.5, 1, 2, .6, 1.3, 1.1); b <- c(-1, -.5, 0, .5, 1, -1.5, .3, .8)
d2 <- as.data.frame(sapply(1:8, function(j) stats::rbinom(n, 1, stats::plogis(a[j] * (th - b[j]))))); names(d2) <- paste0("i", 1:8)
r2 <- quieto(tri_ajustar(d2, names(d2), modelo = "2PL"))
teste("2PL recupera os parâmetros simulados: correlação > .90 para a e > .97 para b",
      stats::cor(r2$parametros$a, a) > .90 && stats::cor(r2$parametros$b, b) > .97)
teste("ajuste global = mirt::M2 do mesmo ajuste (M2 para dicotômicos; tol. 1e-8)", {
  m <- mirt::M2(r2$modelo_mirt); r2$tabela$estatistica == "M2" && abs(r2$tabela$valor - m$M2) < 1e-8 && r2$tabela$gl == m$df })
teste("confiabilidade marginal = var / (var + E[1/I(θ)]) sobre N(0, var), conferida numa grade independente (tol. 1e-4)", {
  f <- r2$modelo_mirt; th_g <- seq(-8, 8, length.out = 6001); I <- mirt::testinfo(f, matrix(th_g)); w <- stats::dnorm(th_g); w <- w / sum(w)
  abs(r2$efeito$valor[1] - 1 / (1 + sum(w / I))) < 1e-4 })
teste("Q3 relativo = Q3 menos a média dos Q3 fora da diagonal (Christensen et al., 2017), com os Q3 do mirt", {
  q <- quieto(mirt::residuals(r2$modelo_mirt, type = "Q3", verbose = FALSE)); fora <- q[upper.tri(q)]
  abs(r2$efeito$valor[3] - max(fora - mean(fora))) < 1e-8 && abs(r2$efeito$valor[4] - mean(fora)) < 1e-8 })

cat("== Rasch: confiabilidade com variância latente diferente de 1\n")
set.seed(9); thr <- stats::rnorm(1500, 0, 1.6); br <- seq(-2, 2, length.out = 12)
dr <- as.data.frame(sapply(br, function(bj) stats::rbinom(1500, 1, stats::plogis(thr - bj)))); names(dr) <- paste0("r", 1:12)
rr <- quieto(tri_ajustar(dr, names(dr), modelo = "rasch"))
teste("Rasch: marginal da skill perto da empírica EAP (|dif| < .03), usando a variância latente estimada", {
  abs(rr$efeito$valor[1] - rr$efeito$valor[2]) < .03 })
if (requireNamespace("eRm", quietly = TRUE)) {
  # A de separação do eRm usa os θ de máxima verossimilhança condicional (extremos
  # extrapolados) e sai um pouco abaixo; a conferência é de ordem de grandeza.
  teste("Rasch: marginal da skill na mesma ordem da confiabilidade de separação do eRm (|dif| < .08)", {
    e <- eRm::RM(dr); pp <- eRm::person.parameter(e); sr <- eRm::SepRel(pp)$sep.rel
    abs(rr$efeito$valor[1] - sr) < .08 })
} else cat("PULADO  comparação com eRm::SepRel (não instalado)\n")
teste("divergência documentada: mirt::marginal_rxx integra em N(0, 1) e se afasta da empírica com variância latente 1,6² (> .10); a da skill não (< .03)", {
  abs(quieto(mirt::marginal_rxx(rr$modelo_mirt)) - rr$efeito$valor[2]) > .10 && abs(rr$efeito$valor[1] - rr$efeito$valor[2]) < .03 })

cat("== DIF e comparação de modelos\n")
set.seed(8); nd <- 1200; g <- rep(c("ref", "foc"), each = nd / 2); thd <- stats::rnorm(nd)
ad <- c(1.2, .8, 1.5, 1, 2, .6, 1.3, 1.1, 1, .9); bd <- c(-1, -.5, 0, .5, 1, -1.5, .3, .8, -.2, .4)
bb <- matrix(bd, nd, 10, byrow = TRUE); bb[g == "foc", 3] <- bb[g == "foc", 3] + .8
dd <- as.data.frame(sapply(1:10, function(j) stats::rbinom(nd, 1, stats::plogis(ad[j] * (thd - bb[, j]))))); names(dd) <- paste0("i", 1:10); dd$g <- g
dl <- quieto(tri_dif(dd, paste0("i", 1:10), "g", metodo = "logistica", purificar = FALSE))
teste("DIF logístico: só o item com deslocamento simulado (i3) marcado, com ΔR² de Nagelkerke na faixa B de Jodoin e Gierl", {
  identical(dl$tabela$item[dl$tabela$dif], "i3") && dl$tabela$classe[3] == "B (moderado)" })
if (requireNamespace("difR", quietly = TRUE)) {
  teste("DIF logístico igual ao difR::difLogistic (type = both, LRT, sem purificar): estatística e ΔR² (tol. 1e-8)", {
    r <- difR::difLogistic(dd[, 1:10], dd$g, focal.name = "foc", type = "both", criterion = "LRT", purify = FALSE)
    max(abs(r$Logistik - dl$tabela$qui2_conjunto)) < 1e-8 && max(abs(r$deltaR2 - dl$tabela$delta_r2_nagelkerke)) < 1e-8 })
} else cat("PULADO  comparação com difR (não instalado)\n")
dm <- quieto(tri_dif(dd, paste0("i", 1:10), "g", modelo = "2PL", metodo = "mirt"))
teste("DIF por mirt (2PL multigrupo, âncoras em duas etapas): i3 marcado, nenhum outro", identical(dm$tabela$item[dm$tabela$dif], "i3"))
rc <- quieto(tri_comparar(dd, paste0("i", 1:10), modelos = c("rasch", "2PL")))
teste("comparação Rasch contra 2PL: razão de verossimilhança = 2ΔlogLik com gl = Δparâmetros", {
  t <- rc$tabela; l <- rc$posthoc
  abs(l$qui2 - 2 * (t$logLik[2] - t$logLik[1])) < 1e-6 && l$gl == t$parametros[2] - t$parametros[1] })

cat("== Graded politômico\n")
set.seed(11); np <- 800; thp <- stats::rnorm(np)
dp <- as.data.frame(sapply(1:6, function(j) { ys <- 1.2 * thp + stats::rnorm(np); as.integer(cut(ys, c(-Inf, -1, 0, 1, Inf))) })); names(dp) <- paste0("p", 1:6)
rg <- quieto(tri_ajustar(dp, names(dp), modelo = "graded"))
teste("graded: ajuste global pelo M2* quando há gl e, sem gl (6 itens de 4 categorias), pelo C2; igual ao mirt::M2 do mesmo tipo", {
  tp <- rg$tabela$estatistica; m <- mirt::M2(rg$modelo_mirt, type = tp)
  tp == "C2" && abs(rg$tabela$valor - m$M2) < 1e-6 && rg$tabela$gl == m$df })
teste("privacidade: a impressão traz só agregados (ajuste de pessoas agregado, sem θ individual)", length(capture.output(print(rg))) < 150)

cat("\n== TRI multidimensional (v1.5)\n")
set.seed(31); n_m <- 800
a_m <- matrix(0, 10, 2); a_m[1:5, 1] <- stats::runif(5, 1.2, 2); a_m[6:10, 2] <- stats::runif(5, 1.2, 2)
d_m <- matrix(rep(c(1, 0, -1), each = 10), 10) + stats::rnorm(10, sd = .3)
th_m <- MASS::mvrnorm(n_m, c(0, 0), matrix(c(1, .4, .4, 1), 2))
dat_m <- as.data.frame(mirt::simdata(a_m, d_m, n_m, itemtype = "graded", Theta = th_m)); names(dat_m) <- paste0("i", 1:10)
me <- suppressWarnings(tri_multidimensional(dat_m, names(dat_m), 2))
teste("multidimensional exploratório: cada item carrega mais na dimensão simulada", {
  cg <- as.matrix(me$tabelas_extra[[1]][, 2:3]); dom <- apply(abs(cg), 1, which.max)
  length(unique(dom[1:5])) == 1 && length(unique(dom[6:10])) == 1 && dom[1] != dom[6] })
teste("multidimensional exploratório: correlação entre dimensões perto de .40 (tol. .10)", abs(abs(me$tabelas_extra[[2]][1, 3]) - .4) < .1)
teste("multidimensional: 2 dimensões vencem o unidimensional no BIC e na razão de verossimilhança", {
  cp <- me$tabelas_extra[["Comparação com o unidimensional"]]; cp$bic[2] < cp$bic[1] && cp$p[2] < .001 })
mc <- suppressWarnings(tri_multidimensional(dat_m, names(dat_m), list(A = paste0("i", 1:5), B = paste0("i", 6:10)), comparar_unidimensional = FALSE))
teste("multidimensional confirmatório igual à chamada direta do mirt (logLik, tol. 1e-4)", {
  set.seed(2026)
  f <- suppressWarnings(mirt::mirt(as.data.frame(lapply(dat_m, function(x) x - min(x))), mirt::mirt.model("A = 1-5\nB = 6-10\nCOV = A*B"),
                                   itemtype = "graded", SE = TRUE, verbose = FALSE))
  abs(mirt::extract.mirt(f, "logLik") - mirt::extract.mirt(mc$modelo, "logLik")) < 1e-4 })
teste("confiabilidade empírica perto do quadrado da correlação entre θ verdadeiro e EAP (tol. .05)", {
  th <- mirt::fscores(mc$modelo, verbose = FALSE)
  r2 <- c(stats::cor(th[, 1], th_m[, 1])^2, stats::cor(th[, 2], th_m[, 2])^2)
  all(abs(mc$efeito$valor - r2) < .05) })
teste("Rasch multidimensional exploratório recusado (só confirmatório)",
      inherits(try(tri_multidimensional(dat_m, names(dat_m), 2, modelo = "pcm"), silent = TRUE), "try-error"))
teste("3 dimensões: MHRM escolhido sozinho, como nas aulas", {
  m3 <- suppressWarnings(tri_multidimensional(dat_m[1:400, 1:9], names(dat_m)[1:9], list(A = paste0("i", 1:3), B = paste0("i", 4:6), C = paste0("i", 7:9)),
                                              comparar_unidimensional = FALSE))
  grepl("MHRM", m3$decisoes[1]) })
teste("privacidade: o objeto da TRI multidimensional não guarda θ por pessoa fora do modelo", !any(vapply(me[setdiff(names(me), "modelo")], function(x) is.matrix(x) && nrow(x) == n_m, TRUE)))

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
print(versoes_pacotes(c("mirt", "difR", "eRm")))
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
