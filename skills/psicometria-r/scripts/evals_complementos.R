# evals_complementos.R, testes das funções da v1.5 fora dos módulos novos: classes_distal
# (BCH), rede_bayesiana, confiabilidade_metades, ausentes, controle_desejabilidade,
# escolha_forcada, tendencia_joinpoint e taxa_padronizada.
#
# Rodar da pasta da skill:  Rscript scripts/evals_complementos.R
# Dados simulados com semente fixa, bancos de exemplo de pacotes (bnlearn) e a incidência
# de tuberculose 2010-2019 da aula de séries temporais (dados públicos agregados).

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

cat("== classes_distal (três passos, BCH)\n")
set.seed(51); n <- 1500; cl <- sample(1:2, n, TRUE, prob = c(.6, .4))
pr <- rbind(c(.25, .3, .2, .25, .3), c(.7, .75, .65, .7, .7))
di <- as.data.frame(sapply(1:5, function(j) 1 + stats::rbinom(n, 1, pr[cl, j]))); names(di) <- paste0("u", 1:5)
di$z <- c(0, 1)[cl] + stats::rnorm(n)
lc <- quieto(classes_latentes(di, paste0("u", 1:5), k_max = 2, nrep = 5))
bd <- quieto(classes_distal(lc, di, "z", k = 2, tipo = "continua"))
ordem <- order(bd$tabela$media_bch)
teste("BCH recupera as médias simuladas (0 e 1; tol. .15) e corrige o viés da atribuição modal", {
  m <- bd$tabela$media_bch[ordem]; mm <- bd$tabela$media_modal[ordem]
  abs(m[1] - 0) < .15 && abs(m[2] - 1) < .15 && abs(diff(m) - 1) < abs(diff(mm) - 1) })
teste("BCH: classificação perfeita (posterior 0/1) = médias da atribuição modal", {
  r2 <- lc; P <- r2$posterior[[2]]; r2$posterior[[2]] <- diag(2)[max.col(P), ]
  b2 <- quieto(classes_distal(r2, di, "z", k = 2, tipo = "continua")); max(abs(b2$tabela$media_bch - b2$tabela$media_modal)) < 1e-10 })
teste("BCH: diferença entre classes detectada (Wald p < .001)", bd$efeito$valor > 10 && grepl("p < .001", bd$apa))
teste("BCH: soma dos pesos por pessoa = 1 (D⁻¹ preserva o total)", { P <- lc$posterior[[2]]; W <- diag(2)[max.col(P), ]
  D <- t(W) %*% P / colSums(W); max(abs(rowSums(W %*% solve(D)) - 1)) < 1e-10 })
di$cat <- ifelse(stats::runif(n) < c(.2, .7)[cl], "sim", "nao")
bc <- quieto(classes_distal(lc, di, "cat", k = 2))
teste("BCH categórica: proporções somam 1 por classe", { t <- bc$tabela; all(abs(tapply(t$proporcao_bch, t$classe, sum) - 1) < 1e-10) })

cat("== rede_bayesiana (bnlearn)\n")
data(gaussian.test, package = "bnlearn")
rb <- quieto(rede_bayesiana(gaussian.test, names(gaussian.test), reamostras = 100))
verdade <- bnlearn::model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
teste("gaussian.test (exemplo do bnlearn usado nas aulas): o esqueleto da rede média é o da rede geradora", {
  par <- function(x) sort(apply(x, 1, function(r) paste(sort(r), collapse = "-")))
  identical(unique(par(as.matrix(rb$tabela[, 1:2]))), unique(par(bnlearn::arcs(verdade)))) })
teste("forças iguais à chamada direta do boot.strength com a mesma semente", {
  set.seed(2026); bs <- bnlearn::boot.strength(gaussian.test, R = 100, algorithm = "hc")
  m <- merge(rb$tabela, bs, by.x = c("de", "para"), by.y = c("from", "to")); max(abs(m$forca - m$strength)) < 1e-12 })
teste("variáveis com até 7 valores viram rede discreta", {
  dd <- as.data.frame(lapply(gaussian.test[1:300, 1:4], function(x) cut(x, 3, labels = FALSE)))
  grepl("discreta", quieto(rede_bayesiana(dd, names(dd), reamostras = 30))$analise) })

cat("== confiabilidade_metades\n")
set.seed(52); n <- 400; f <- stats::rnorm(n)
dm <- as.data.frame(sapply(1:8, function(i) round(.7 * f + stats::rnorm(n, sd = .7 + .05 * i), 1))); names(dm) <- paste0("q", 1:8)
mp <- confiabilidade_metades(dm, names(dm))
teste("ímpares × pares: Spearman-Brown = 2r/(1 + r) e λ4 = 2[1 - (s²A + s²B)/s²X]", {
  a <- rowSums(dm[c(1, 3, 5, 7)]); b <- rowSums(dm[c(2, 4, 6, 8)]); r <- stats::cor(a, b)
  abs(mp$tabela$spearman_brown - 2 * r / (1 + r)) < 1e-12 && abs(mp$tabela$guttman_lambda4 - 2 * (1 - (stats::var(a) + stats::var(b)) / stats::var(a + b))) < 1e-12 })
mt <- quieto(confiabilidade_metades(dm, names(dm), metodo = "todas"))
teste("todas as divisões (8 itens, 35 divisões): λ4 médio = alfa (Cronbach, 1951; tol. 1e-10)", abs(mt$tabela$lambda4[1] - mt$tabela$alfa[1]) < 1e-10 && mt$tabela$n_divisoes[1] == 35)
if (requireNamespace("psych", quietly = TRUE)) teste("todas as divisões: mínimo, máximo e média = psych::splitHalf(covar = TRUE) (tol. 1e-8)", {
  sh <- suppressWarnings(psych::splitHalf(dm, brute = TRUE, covar = TRUE))
  abs(mt$tabela$lambda4[3] - sh$maxrb) < 1e-8 && abs(mt$tabela$lambda4[2] - sh$minrb) < 1e-8 && abs(mt$tabela$lambda4[1] - sh$meanr) < 1e-8 })
teste("todas as divisões: mínimo <= médio <= máximo", mt$tabela$lambda4[2] <= mt$tabela$lambda4[1] && mt$tabela$lambda4[1] <= mt$tabela$lambda4[3])

cat("== ausentes (Little)\n")
set.seed(53); n <- 500; S <- matrix(.5, 4, 4); diag(S) <- 1
x <- as.data.frame(MASS::mvrnorm(n, rep(0, 4), S)); names(x) <- paste0("v", 1:4)
xm <- x; xm$v2[stats::runif(n) < .15] <- NA; xm$v4[stats::runif(n) < .10] <- NA
am <- ausentes(xm, names(xm))
if (requireNamespace("naniar", quietly = TRUE)) {
  teste("Little igual ao naniar::mcar_test (χ², gl, p; tol. 1e-6)", { nt <- naniar::mcar_test(xm)
    abs(am$little$qui2 - nt$statistic) < 1e-6 && am$little$gl == nt$df && abs(am$little$p - nt$p.value) < 1e-6 })
} else pular("comparação com naniar", "não instalado")
xr <- x; xr$v2[x$v1 > .3] <- NA
ar <- ausentes(xr, names(xr))
teste("ausência que depende de outra variável (MAR forte): Little rejeita MCAR (p < .001)", ar$little$p < .001)
teste("porcentagens por variável e padrões com n < 5 fora do print", {
  out <- utils::capture.output(print(am)); abs(am$tabela$pct_ausentes[2] - 100 * mean(is.na(xm$v2))) < 1e-10 && all(am$tabelas_extra[[2]]$n >= 5) })

cat("== controle_desejabilidade (quádruplas)\n")
set.seed(54); n <- 1200; Tr <- stats::rnorm(n); DS <- stats::rnorm(n); Aq <- stats::rnorm(n, sd = .3); crit <- .5 * Tr + stats::rnorm(n, sd = .87)
pt <- c(1, 1, 1, 1, -1, -1, -1, -1); pdes <- c(1, -1, 1, -1, 1, -1, 1, -1)
dq <- as.data.frame(sapply(1:8, function(i) .6 * pt[i] * Tr + .45 * pdes[i] * DS + Aq + stats::rnorm(n, sd = .6))); names(dq) <- paste0("e", 1:8)
dq$crit <- crit
for (i in c(5:8)) dq[[paste0("e", i)]] <- dq[[paste0("e", i)]]
cd <- quieto(controle_desejabilidade(dq, paste0("e", 1:8), setNames(pdes, paste0("e", 1:8)), criterio = "crit", estimador = "ML"))
teste("cargas de DS com o sinal do polo de desejabilidade em todos os itens", all(cd$tabelas_extra[[2]]$sinal_coerente))
teste("com o controle, as cargas no traço ficam mais perto das simuladas (|.6| padronizada pela variância)", {
  cg <- cd$tabelas_extra[[1]]; v <- .6^2 + .45^2 + .3^2 + .6^2; alvo <- .6 / sqrt(v)
  mean(abs(abs(cg$traco_desejabilidade) - alvo)) < mean(abs(abs(cg$traco) - alvo)) })
teste("validade contra o critério: o modelo com controle fica mais perto do β = .5 simulado", {
  v <- cd$tabelas_extra[[3]]; abs(abs(v$beta[v$modelo == "traco_desejabilidade"]) - .5) <= abs(abs(v$beta[v$modelo == "traco"]) - .5) })
teste("o fator de desejabilidade melhora o ajuste (Δχ² p < .001)", tail(cd$tabela$p_dif_qui2, 1) < .001)

cat("== escolha_forcada (TRI thurstoniana)\n")
if (requireNamespace("thurstonianIRT", quietly = TRUE)) {
  set.seed(55); nt <- 5; ni <- 30
  lam <- stats::runif(ni, .65, .95) * sample(c(1, 1, -1), ni, TRUE)
  Phi <- diag(nt); Phi[1, 2] <- Phi[2, 1] <- .3; Phi[3, 4] <- Phi[4, 3] <- -.3
  sim <- quieto(thurstonianIRT::sim_TIRT_data(npersons = 1000, ntraits = nt, nitems_per_block = 3, nblocks_per_trait = 6,
                                              gamma = stats::rnorm(ni, 0, .5), lambda = lam, Phi = Phi))
  ef <- quieto(escolha_forcada(sim, blocos = list(list(itens = "x", tracos = "t", sinais = 1))))
  teste("5 traços, sinais mistos: correlações entre traços perto das simuladas (.3 e -.3; tol. .12), já no sentido declarado", {
    r <- as.matrix(ef$tabelas_extra[[1]][, -1]); abs(r[1, 2] - .3) < .12 && abs(r[3, 4] + .3) < .12 })
  teste("confiabilidade empírica perto de cor(θ verdadeiro, θ estimado)² (tol. .08)", {
    p <- as.data.frame(stats::predict(ef$modelo)); th <- attr(sim, "eta")
    est <- sapply(paste0("trait", 1:nt), function(tt) p$estimate[p$trait == tt][order(p$id[p$trait == tt])])
    r2 <- vapply(seq_len(nt), function(j) stats::cor(est[, j], th[, j])^2, 0)
    all(abs(ef$efeito$valor - r2) < .08) })
  teste("ajuste aceitável no modelo correto (RMSEA < .05)", ef$tabela$rmsea < .05)
} else pular("escolha forçada", "thurstonianIRT não instalado")

cat("== tendencia_joinpoint (incidência de TB 2010-2019, aula de séries temporais)\n")
tb <- data.frame(ano = rep(2010:2019, 3), sexo = rep(c("masc", "fem", "ambos"), each = 10),
                 casos = c(46870, 48299, 47465, 47334, 46916, 47485, 48049, 51193, 52943, 53458, 24044, 24483, 23653, 23850, 23118, 22320, 22468, 22893, 23972, 24426,
                           70915, 72783, 71128, 71192, 70039, 69809, 70519, 74090, 76923, 77891),
                 populacao = c(96706703, 97610297, 98487258, 99336858, 100159507, 100955522, 101726102, 102471274, 103189829, 103881681,
                               98791094, 99786721, 100755204, 101695856, 102609055, 103495127, 104355330, 105189655, 105996973, 106777332,
                               195497797, 197397018, 199242462, 201032714, 202768562, 204450649, 206081432, 207660929, 209186802, 210659013))
# Saída do Joinpoint 4.9.1 no relato da aula (APC e AAPC com uma casa)
jp_ref <- list(masc = list(ano = 2015, v = c(-1.0, 2.7, 0.7)), fem = list(ano = 2016, v = c(-2.4, 2.3, -0.8)), ambos = list(ano = 2015, v = c(-1.6, 2.2, 0.1)))
for (s in names(jp_ref)) {
  r <- tendencia_joinpoint(tb[tb$sexo == s, ], "ano", "casos", "populacao")
  teste(sprintf("%s: ponto em %d e APC/AAPC = Joinpoint 4.9.1 com uma casa (%s)", s, jp_ref[[s]]$ano, paste(jp_ref[[s]]$v, collapse = "; ")), {
    grepl(as.character(jp_ref[[s]]$ano), r$tabela$periodo[1]) && all(round(r$tabela$valor, 1) == jp_ref[[s]]$v) })
}
teste("sem ponto (max_pontos = 0): APC = 100(e^b - 1) da reta em log(taxa)", {
  set.seed(56); a <- 2000:2019; y <- exp(3 + .02 * (a - 2000) + stats::rnorm(20, sd = .01))
  r <- tendencia_joinpoint(data.frame(ano = a, tx = y), "ano", taxa = "tx", max_pontos = 0)
  b <- stats::coef(stats::lm(log(y) ~ a))[2]; nrow(r$tabela) == 2 && abs(r$tabela$valor[1] - 100 * (exp(b) - 1)) < 1e-8 })
teste("série com mudança nítida em 2010: um ponto, em 2010, com APC antes e depois perto das simuladas (+3% e -2%)", {
  set.seed(57); a <- 2000:2019; y <- exp(3 + .03 * pmin(a - 2000, 10) - .02 * pmax(a - 2010, 0) + stats::rnorm(20, sd = .004))
  r <- tendencia_joinpoint(data.frame(ano = a, tx = y), "ano", taxa = "tx")
  any(grepl("2000-2010", r$tabela$periodo)) && abs(r$tabela$valor[1] - 100 * (exp(.03) - 1)) < .3 && abs(r$tabela$valor[2] - 100 * (exp(-.02) - 1)) < .3 })
teste("IC empírico: reproduzível com a mesma semente e contém a estimativa", {
  r1 <- tendencia_joinpoint(tb[tb$sexo == "masc", ], "ano", "casos", "populacao", ic = "empirico", reamostras = 499)
  r2 <- tendencia_joinpoint(tb[tb$sexo == "masc", ], "ano", "casos", "populacao", ic = "empirico", reamostras = 499)
  identical(r1$tabela, r2$tabela) && all(r1$tabela$ic_inf <= r1$tabela$valor & r1$tabela$valor <= r1$tabela$ic_sup) })

cat("== taxa_padronizada (método direto)\n")
dp <- data.frame(idade = c("0-19", "20-39", "40-59", "60+"), casos = c(5, 20, 60, 150), pop = c(20000, 30000, 25000, 10000), padrao = c(30000, 30000, 25000, 15000))
tp <- taxa_padronizada(dp, "casos", "pop", "padrao")
teste("taxa padronizada = Σ (padrão/total) × casos/pop", abs(tp$tabela$taxa_padronizada - 1e5 * sum(dp$padrao / sum(dp$padrao) * dp$casos / dp$pop)) < 1e-8)
if (requireNamespace("epitools", quietly = TRUE)) {
  teste("igual ao epitools::ageadjust.direct (taxa e IC gama de Fay-Feuer; tol. 1e-10)", {
    e <- epitools::ageadjust.direct(dp$casos, dp$pop, stdpop = dp$padrao) * 1e5
    abs(tp$tabela$taxa_padronizada - e[["adj.rate"]]) < 1e-8 && abs(tp$tabela$ic_inf - e[["lci"]]) < 1e-8 && abs(tp$tabela$ic_sup - e[["uci"]]) < 1e-8 })
} else pular("comparação com epitools", "não instalado")

ok <- unlist(resultados)
cat(sprintf("\nevals_complementos: %d/%d passaram\n", sum(ok), length(ok)))
if (!all(ok)) quit(status = 1)
