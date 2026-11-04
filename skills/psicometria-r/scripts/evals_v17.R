# evals_v17.R, testes das funções da v1.7: imputacao_multipla, rede_temporal (graphicalVAR
# e mlVAR), escala_mokken, simular_cat, normas, equivalencia_ajuste e invariancia_continua.
#
# Rodar da pasta da skill:  Rscript scripts/evals_v17.R
# Âncoras: os três testes de equivalência de Petersen (2024, 15.10.3.5), as fórmulas de
# escores referidos à norma de Petersen (2024, cap. 1), o banco acl do mokken e
# HolzingerSwineford1939 do lavaan; o resto é verdade conhecida em dados simulados e
# igualdade com a chamada direta do pacote (mice::pool, graphicalVAR, mlVAR, mokken,
# mirtCAT, strucchange::sctest).

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
pular <- function(nome, motivo) cat(sprintf("PULADO  %s (%s)\n", nome, motivo))
quieto <- function(expr) suppressWarnings(suppressMessages(expr))
tem <- function(p) requireNamespace(p, quietly = TRUE)

cat("== imputacao_multipla (mice)\n")
if (tem("mice")) {
  set.seed(10); n <- 500; x <- rnorm(n); z <- .6 * x + rnorm(n, sd = .8); y <- 1 + .5 * x + .3 * z + rnorm(n)
  d <- data.frame(y, x, z); d$y[runif(n) < plogis(-1 + 1.5 * z)] <- NA; d$x[runif(n) < .1] <- NA
  mi <- quieto(imputacao_multipla(d, c("y", "x", "z"), y ~ x, m = 20))
  fit <- with(mi$imputacoes, lm(y ~ x)); pl <- summary(mice::pool(fit)); pp <- mice::pool(fit)$pooled
  teste("estimativa, EP, gl e p iguais ao mice::pool (Rubin com gl de Barnard-Rubin; 1e-8)",
        max(abs(mi$combinado$B - pl$estimate)) < 1e-8 && max(abs(mi$combinado$ep - pl$std.error)) < 1e-8 &&
          max(abs(mi$combinado$gl - pl$df)) < 1e-6 && max(abs(mi$combinado$p - pl$p.value)) < 1e-10)
  teste("FMI igual ao mice::pool e R² combinado igual ao pool.r.squared (1e-8)",
        max(abs(mi$combinado$fmi - pp$fmi)) < 1e-8 && abs(mi$efeito$valor[2] - mice::pool.r.squared(fit)[1, "est"]) < 1e-8)
  teste("ausência MAR ligada a z: a inclinação de y em x (verdadeira .5 + .3 × .6 = .68) fica perto da verdade com a imputação (tol. .08)",
        abs(mi$combinado$B[2] - .68) < .08)
  teste("aviso quando m é menor que a % de casos incompletos (White et al., 2011)", any(grepl("White", mi$avisos)))
  teste("logística: OR combinado = exp(B) combinado", {
    d2 <- d; d2$yb <- as.integer(d2$y > 1); r2 <- quieto(imputacao_multipla(d2[c("yb", "x", "z")], c("yb", "x", "z"), yb ~ x, familia = "binomial", m = 5))
    abs(r2$efeito$valor[1] - exp(r2$combinado$B[2])) < 1e-12 })
} else pular("imputação múltipla", "mice não instalado")

cat("== rede_temporal (Isvoranu et al., 2022, cap. 10)\n")
if (tem("graphicalVAR")) {
  set.seed(1); gm <- graphicalVAR::randomGVARmodel(4, probKappaEdge = .3, probBetaEdge = .3)
  ts <- as.data.frame(graphicalVAR::graphicalVARsim(200, gm$beta, gm$kappa)); names(ts) <- paste0("v", 1:4)
  rt <- quieto(rede_temporal(ts, names(ts), destendenciar = FALSE))
  teste("N = 1: PDC e PCC iguais à chamada direta do graphicalVAR (1e-12)", {
    g <- quieto(graphicalVAR::graphicalVAR(ts, vars = names(ts), verbose = FALSE)); max(abs(g$PDC - rt$temporal)) < 1e-12 && max(abs(g$PCC - rt$contemporanea)) < 1e-12 })
  teste("N = 1: recupera a rede temporal simulada (r > .90 com a PDC verdadeira)", stats::cor(as.vector(gm$PDC), as.vector(rt$temporal)) > .90)
  teste("destendenciar detecta a tendência linear plantada (p < .001) e a retira antes da rede", {
    ts2 <- ts; ts2$v1 <- ts2$v1 + seq(0, 3, length.out = nrow(ts2)); r2 <- quieto(rede_temporal(ts2, names(ts2)))
    any(r2$tendencias$variavel == "v1" & r2$tendencias$p < .001) })
} else pular("rede temporal N = 1", "graphicalVAR não instalado")
if (tem("mlVAR")) {
  set.seed(2); sv <- mlVAR::mlVARsim(nPerson = 40, nNode = 4, nTime = 60, lag = 1)
  rm_ <- quieto(rede_temporal(sv$Data, sv$vars, id = sv$idvar, destendenciar = FALSE))
  mm <- quieto(mlVAR::mlVAR(sv$Data, vars = sv$vars, idvar = sv$idvar, lags = 1, temporal = "orthogonal", contemporaneous = "orthogonal", verbose = FALSE))
  teste("N > 1: efeitos fixos temporais e p iguais ao mlVAR direto (1e-10)", {
    M <- mm$results$Beta; i <- cbind(match(rm_$tabela$para, sv$vars), match(rm_$tabela$de, sv$vars), 1)
    max(abs(rm_$tabela$b_fixo - M$mean[i])) < 1e-10 && max(abs(rm_$tabela$p - M$P[i])) < 1e-10 })
  teste("N > 1: recupera os efeitos fixos simulados (r > .95)", stats::cor(as.vector(sv$model$Beta$mean[, , 1]), as.vector(mm$results$Beta$mean[, , 1])) > .95)
  teste("N > 1: redes contemporânea e entre pessoas presentes", all(c("Rede contemporânea", "Rede entre pessoas") %in% names(rm_$tabelas_extra)))
} else pular("rede temporal multinível", "mlVAR não instalado")

cat("== escala_mokken (Meijer & Tendeiro, 2018, Wiley, cap. 15)\n")
if (tem("mokken")) {
  data(acl, package = "mokken", envir = environment()); A <- as.data.frame(acl[, 1:10])
  mk <- quieto(escala_mokken(A, names(A), iio = TRUE))
  teste("H, Hi e EP iguais ao mokken::coefH (1e-12)", {
    h <- quieto(mokken::coefH(as.matrix(A), se = TRUE, nice.output = FALSE))
    abs(mk$H - as.numeric(h$H)) < 1e-12 && max(abs(mk$tabela$Hi - as.numeric(h$Hi))) < 1e-12 && max(abs(mk$tabela$ep_Hi - as.numeric(h$se.Hi))) < 1e-12 })
  teste("partição da AISP igual ao mokken::aisp", all(mk$tabela$escala_aisp == quieto(mokken::aisp(as.matrix(A)))[, 1]))
  teste("H = .26 no acl (10 primeiros itens): classe 'não escalável' e aviso de Hi < ,30", abs(mk$H - .2635) < 5e-4 && grepl("não escalável", mk$apa) && any(grepl("Hi <", mk$avisos)))
  teste("dados de Guttman perfeitos dão H = 1", {
    set.seed(3); th <- sort(rnorm(200)); b <- c(-1, -.3, .4, 1.1); G <- as.data.frame(sapply(b, function(bb) as.integer(th > bb)))
    abs(quieto(escala_mokken(G, names(G), aisp = FALSE, monotonicidade = FALSE))$H - 1) < 1e-10 })
} else pular("Mokken", "mokken não instalado")

cat("== simular_cat (Martinková & Hladká, 2023, 10.2.6)\n")
if (tem("mirtCAT")) {
  set.seed(5); a <- matrix(stats::rlnorm(30, .2, .3)); dd <- matrix(stats::rnorm(30))
  XC <- as.data.frame(mirt::simdata(a, dd, 150, itemtype = "2PL")); names(XC) <- paste0("i", 1:30)
  sc <- quieto(simular_cat(XC, names(XC)))
  teste("número médio de itens igual ao mirtCAT direto com o mesmo modelo e desenho", {
    s <- quieto(mirtCAT::mirtCAT(mo = sc$modelo, local_pattern = as.matrix(XC), start_item = "MI", method = "MAP", criteria = "MI", design = list(min_SEM = .30)))
    abs(mean(vapply(s, function(x) length(x$items_answered), 0L)) - sc$tabela$itens_media) < 1e-12 })
  teste("CAT usa menos itens que o banco e o θ correlaciona > .90 com o teste completo", sc$tabela$itens_media < 30 && sc$tabela$r_theta_completo > .90)
  teste("max_itens limita o teste", quieto(simular_cat(XC, names(XC), fit = sc$modelo, ep_max = .1, max_itens = 8))$tabela$itens_max == 8)
} else pular("CAT", "mirtCAT não instalado")

cat("== normas (Petersen, 2024, cap. 1)\n")
set.seed(1); dn <- data.frame(x = round(stats::rnorm(200, 20, 4)), g = sample(c("jovem", "idoso"), 200, TRUE))
nm <- normas(dn, "x", percentil = "petersen")
teste("percentil = dplyr::percent_rank do livro: (posto mínimo - 1)/(n - 1)", {
  x <- dn$x; esperado <- vapply(nm$normas$bruto, function(v) 100 * (rank(x, ties.method = "min")[match(v, x)] - 1) / (length(x) - 1), 0)
  isTRUE(all.equal(nm$normas$percentil, esperado)) })
teste("z, T = 50 + 10z, escalar = 10 + 3z e padrão = 100 + 15z", {
  t <- nm$normas; z <- (t$bruto - mean(dn$x)) / sd(dn$x)
  max(abs(t$z - z)) < 1e-12 && max(abs(t$T - (50 + 10 * z))) < 1e-12 && max(abs(t$escalar - (10 + 3 * z))) < 1e-12 && max(abs(t$padrao - (100 + 15 * z))) < 1e-12 })
teste("estanino pelas faixas de Petersen (4, 11, 23, 40, 60, 77, 89, 96%)", {
  t <- nm$normas; e <- findInterval(t$percentil, c(4, 11, 23, 40, 60, 77, 89, 96), left.open = TRUE) + 1
  all(t$estanino == e) })
teste("percentil médio e estratos: cada estrato com a própria média", {
  ng <- normas(dn, "x", "g"); r <- ng$resumo
  abs(r$media[r$grupo == "jovem"] - mean(dn$x[dn$g == "jovem"])) < 1e-12 && all(ng$normas$percentil > 0 & ng$normas$percentil < 100) })

cat("== equivalencia_ajuste (Counsell et al., 2020; Petersen, 2024, 15.10.3.5)\n")
eq <- equivalencia_ajuste(c(70.9, 75.2, 5.146), c(48, 54, 6), n = 259, grupos = 2)
teste("os três exemplos do livro: Fml .2759/.2926/.02002, ε0 .1536/.1728/.0192, λ 39.48/44.41/4.934, p .147/.076/.1375",
      all(abs(eq$tabela$Fml - c(.2759, .2926, .02002)) < 6e-5) && all(abs(eq$tabela$epsilon0 - c(.1536, .1728, .0192)) < 1e-10) &&
        all(abs(eq$tabela$lambda - c(39.48, 44.41, 4.934)) < .006) && all(abs(eq$tabela$p - c(.147, .07603, .1375)) < 6e-4))

cat("== invariancia_continua (Merkle & Zeileis, 2013)\n")
if (tem("strucchange")) {
  HS <- lavaan::HolzingerSwineford1939
  fh <- lavaan::cfa("visual =~ x1 + x2 + x3\n textual =~ x4 + x5 + x6\n speed =~ x7 + x8 + x9", data = HS, meanstructure = TRUE)
  teste("categórico (escola, LMuo): igual ao strucchange::sctest direto (1e-10)", {
    r <- invariancia_continua(fh, HS$school, "cargas")
    s <- strucchange::sctest(fh, order.by = factor(HS$school), parm = 1:6, functional = strucchange::catL2BB(factor(HS$school)))
    abs(r$tabela$estatistica - unname(s$statistic)) < 1e-10 && abs(r$tabela$p - s$p.value) < 1e-10 })
  set.seed(3); n <- 600; idade <- stats::runif(n, 18, 70); f <- stats::rnorm(n)
  mk_d <- function(l4) data.frame(y1 = .7 * f + stats::rnorm(n, 0, .7), y2 = .7 * f + stats::rnorm(n, 0, .7), y3 = .7 * f + stats::rnorm(n, 0, .7), y4 = l4 * f + stats::rnorm(n, 0, .7))
  d1 <- mk_d(.3 + .012 * (idade - 18)); d0 <- mk_d(.7)
  p1 <- invariancia_continua(lavaan::cfa("F =~ y1+y2+y3+y4", data = d1, meanstructure = TRUE), idade, "cargas")$tabela$p
  p0 <- invariancia_continua(lavaan::cfa("F =~ y1+y2+y3+y4", data = d0, meanstructure = TRUE), idade, "cargas")$tabela$p
  teste("carga que cresce com a idade é detectada (p < .001) e a invariante não (p > .05)", p1 < .001 && p0 > .05)
  teste("ordinal usa WDMo e recusa moderador com ausentes", {
    r <- invariancia_continua(fh, ordered(HS$ageyr), "cargas"); r$tabela$funcional == "WDMo" &&
      inherits(try(invariancia_continua(fh, ordered(HS$grade), "cargas"), silent = TRUE), "try-error") })
} else pular("invariância contínua", "strucchange não instalado")

ok <- unlist(resultados)
cat(sprintf("\nevals_v17: %d/%d passaram\n", sum(ok), length(ok)))
if (!all(ok)) quit(status = 1)
