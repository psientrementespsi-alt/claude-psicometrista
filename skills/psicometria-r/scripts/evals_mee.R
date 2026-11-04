# evals_mee.R, testes de est_mee.R (v1.5): MEE completa, multigrupo, interação latente,
# curva de crescimento e invariância longitudinal.
#
# Rodar da pasta da skill:  Rscript scripts/evals_mee.R
# Dados simulados com semente fixa (lavaan::simulateData e simulação própria).

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

cat("== mee (modelo estrutural com efeito indireto)\n")
pop <- "X =~ 0.7*x1 + 0.7*x2 + 0.7*x3
M =~ 0.7*m1 + 0.7*m2 + 0.7*m3
Y =~ 0.7*y1 + 0.7*y2 + 0.7*y3
M ~ 0.5*X
Y ~ 0.4*M + 0.1*X"
set.seed(21); dm <- lavaan::simulateData(pop, sample.nobs = 2000, standardized = TRUE)
mod <- "X =~ x1 + x2 + x3
M =~ m1 + m2 + m3
Y =~ y1 + y2 + y3
M ~ a*X
Y ~ b*M + c*X
indireto := a*b
total := a*b + c"
r <- quieto(mee(mod, dm, estimador = "ML"))
dir <- quieto(lavaan::sem(mod, data = dm, estimator = "ML", missing = "fiml"))
teste("estimativas iguais à chamada direta do lavaan (tol. 1e-8)", {
  pe <- lavaan::parameterEstimates(dir); pe <- pe[pe$op %in% c("~", ":="), ]
  max(abs(r$tabela$b - pe$est)) < 1e-8 })
teste("caminhos padronizados recuperam a população (a = .5, b = .4, c = .1; tol. .08)", {
  b <- setNames(r$tabela$beta, r$tabela$parametro)
  abs(b["M <- X"] - .5) < .08 && abs(b["Y <- M"] - .4) < .08 && abs(b["Y <- X"] - .1) < .08 })
teste("efeito indireto padronizado = produto a × b padronizados (tol. 1e-8)", {
  b <- setNames(r$tabela$beta, r$tabela$parametro); abs(b["indireto"] - b["M <- X"] * b["Y <- M"]) < 1e-8 })
rb <- quieto(mee(mod, dm, estimador = "ML", bootstrap = 300, semente = 5))
teste("bootstrap: IC percentil contém a estimativa e difere do IC delta", {
  t <- rb$tabela; all(t$ic_inf <= t$beta & t$beta <= t$ic_sup) && max(abs(t$ic_inf - r$tabela$ic_inf)) > 1e-4 })
teste("bootstrap reproduzível com a mesma semente", {
  rb2 <- quieto(mee(mod, dm, estimador = "ML", bootstrap = 300, semente = 5)); max(abs(rb2$tabela$ic_inf - rb$tabela$ic_inf)) < 1e-12 })
dord <- dm; for (v in setdiff(names(dord), character())) dord[[v]] <- cut(dord[[v]], c(-Inf, -1, 0, 1, Inf), labels = FALSE)
ro <- quieto(mee(mod, dord))
teste("itens com 4 categorias: WLSMV com parametrização theta escolhido sozinho", grepl("WLSMV", ro$analise) && grepl("theta", ro$decisoes[1]))
teste("sem regressão no modelo: erro com orientação para ajustar_afc()", inherits(try(mee("X =~ x1 + x2 + x3", dm), silent = TRUE), "try-error"))

cat("== mee_multigrupo (moderação categórica)\n")
set.seed(22)
pg <- function(b) paste0("X =~ 0.7*x1 + 0.7*x2 + 0.7*x3\nY =~ 0.7*y1 + 0.7*y2 + 0.7*y3\nW =~ 0.7*w1 + 0.7*w2 + 0.7*w3\nY ~ ", b, "*X + 0.3*W")
dg <- rbind(cbind(lavaan::simulateData(pg(.6), sample.nobs = 500, standardized = TRUE), g = "A"),
            cbind(lavaan::simulateData(pg(0), sample.nobs = 500, standardized = TRUE), g = "B"))
modg <- "X =~ x1 + x2 + x3\nY =~ y1 + y2 + y3\nW =~ w1 + w2 + w3\nY ~ X + W"
mg <- quieto(mee_multigrupo(modg, dg, "g", estimador = "ML"))
teste("igualar as regressões piora o ajuste quando um caminho difere entre grupos (p < .001)", tail(mg$tabela$p_dif_qui2, 1) < .001)
teste("teste por caminho aponta Y ~ X e não Y ~ W", {
  p <- setNames(mg$posthoc$p_dif_qui2, mg$posthoc$caminho); p["Y ~ X"] < .001 && p["Y ~ W"] > .01 })
teste("medida invariante na simulação: passos métrico e escalar sem piora (p > .01)", all(mg$tabela$p_dif_qui2[2:3] > .01))
mg2 <- quieto(mee_multigrupo(modg, dg, "g", estimador = "ML", igualar_latentes = TRUE))
teste("igualar_latentes = TRUE acrescenta o passo das latentes (rota das aulas)", "latentes" %in% mg2$tabela$modelo)

cat("== interacao_latente (indicadores-produto, centralização residual)\n")
set.seed(23); n <- 1500
X <- stats::rnorm(n); Z <- .3 * X + stats::rnorm(n, sd = sqrt(1 - .09)); Yl <- .3 * X + .2 * Z + .25 * X * Z + stats::rnorm(n, sd = .8)
ind <- function(f, p) as.data.frame(sapply(1:3, function(i) .75 * f + stats::rnorm(n, sd = sqrt(1 - .75^2))))
di <- cbind(setNames(ind(X), paste0("x", 1:3)), setNames(ind(Z), paste0("z", 1:3)), setNames(ind(Yl), paste0("y", 1:3)))
ri <- quieto(interacao_latente(di, paste0("x", 1:3), paste0("z", 1:3), paste0("y", 1:3), estimador = "ML"))
teste("interação simulada (b = .25) detectada com p < .001 e ΔR² > 0", ri$tabela$p[ri$tabela$preditor == "XxZ"] < .001 && ri$efeito$valor[2] > 0)
teste("igual ao código das aulas (indProd residualC + sem ML com std.lv), tol. 1e-8", {
  d2 <- semTools::indProd(di, var1 = paste0("x", 1:3), var2 = paste0("z", 1:3), match = TRUE, meanC = FALSE, residualC = TRUE, doubleMC = FALSE)
  m <- paste("X =~ x1 + x2 + x3\nZ =~ z1 + z2 + z3\nY =~ y1 + y2 + y3\nXxZ =~", paste(setdiff(names(d2), names(di)), collapse = " + "), "\nY ~ X + Z + XxZ")
  f <- quieto(lavaan::sem(m, data = d2, estimator = "ML", meanstructure = TRUE, std.lv = TRUE))
  pe <- lavaan::parameterEstimates(f); abs(pe$est[pe$lhs == "Y" & pe$op == "~" & pe$rhs == "XxZ"] - ri$tabela$b[ri$tabela$preditor == "XxZ"]) < 1e-8 })
teste("inclinações simples crescem com o moderador", all(diff(ri$posthoc$inclinacao) > 0))
di0 <- di; Y0 <- .3 * X + .2 * Z + stats::rnorm(n, sd = .8); di0[paste0("y", 1:3)] <- ind(Y0)
r0 <- quieto(interacao_latente(di0, paste0("x", 1:3), paste0("z", 1:3), paste0("y", 1:3), estimador = "ML"))
teste("sem interação na população: p > .01", r0$tabela$p[r0$tabela$preditor == "XxZ"] > .01)

cat("== curva_crescimento (lavaan::growth)\n")
set.seed(24); n <- 600; tt <- 0:3
b0 <- 2 + stats::rnorm(n, sd = 1); b1 <- .5 + stats::rnorm(n, sd = .3)
yy <- sapply(tt, function(t) b0 + b1 * t + stats::rnorm(n, sd = .6)); colnames(yy) <- paste0("y", tt)
dw <- as.data.frame(yy)
cg <- quieto(curva_crescimento(dw, paste0("y", tt), estimador = "ML", residuos_iguais = TRUE))
teste("recupera médias (2; 0,5) e variâncias (1; 0,09) da população (tol. .1)", {
  e <- cg$tabela$estimativa; abs(e[1] - 2) < .1 && abs(e[2] - .5) < .1 && abs(e[3] - 1) < .15 && abs(e[4] - .09) < .05 })
if (requireNamespace("lme4", quietly = TRUE)) {
  teste("resíduos iguais = modelo misto com intercepto e inclinação aleatórios (lmer ML; tol. 1e-3)", {
    dl <- data.frame(id = rep(seq_len(n), length(tt)), t = rep(tt, each = n), y = as.vector(yy))
    lm_ <- lme4::lmer(y ~ t + (t | id), data = dl, REML = FALSE)
    vc <- as.data.frame(lme4::VarCorr(lm_))
    e <- cg$tabela$estimativa
    max(abs(e[1:2] - lme4::fixef(lm_))) < 1e-3 && abs(e[3] - vc$vcov[1]) < 1e-3 && abs(e[4] - vc$vcov[2]) < 1e-3 && abs(e[5] - vc$vcov[3]) < 1e-3 })
} else cat("PULADO  comparação com lme4\n")
cl <- quieto(curva_crescimento(dw, paste0("y", tt), forma = "livre", estimador = "ML"))
teste("base latente: cargas livres a menos de 2 EP de 2 e 3 (crescimento linear simulado)", {
  cb <- cl$tabelas_extra[[2]]; abs(cb$carga[3] - 2) < 2 * cb$ep[3] && abs(cb$carga[4] - 3) < 2 * cb$ep[4] })
teste("base latente com resíduos adjacentes em 4 ondas: modelo saturado avisado", {
  cs <- quieto(curva_crescimento(dw, paste0("y", tt), forma = "livre", residuos_adjacentes = TRUE, estimador = "ML")); any(grepl("saturado", cs$avisos)) })
teste("igual ao modelo da aula (growth MLR, cargas 0:4) com o mesmo banco", {
  m <- "intercepto =~ 1*y0 + 1*y1 + 1*y2 + 1*y3\ninclinacao =~ 0*y0 + 1*y1 + 2*y2 + 3*y3"
  f <- quieto(lavaan::growth(m, data = dw, estimator = "MLR", missing = "fiml"))
  c2 <- quieto(curva_crescimento(dw, paste0("y", tt)))
  abs(lavaan::parameterEstimates(f)$est[lavaan::parameterEstimates(f)$lhs == "inclinacao" & lavaan::parameterEstimates(f)$op == "~1"] - c2$tabela$estimativa[2]) < 1e-8 })

cat("== invariancia_longitudinal\n")
set.seed(25); n <- 700
f1 <- stats::rnorm(n); f2 <- .6 * f1 + stats::rnorm(n, sd = .8)
it <- function(f, int) sapply(1:4, function(i) .7 * f + int[i] + stats::rnorm(n, sd = .7))
w1 <- it(f1, c(0, 0, 0, 0)); w2 <- it(f2, c(0, 0, 0, 0)); w2b <- it(f2, c(0, .8, 0, 0))
mk <- function(a, b) { d <- as.data.frame(cbind(a, b)); names(d) <- c(paste0("i", 1:4, "_t1"), paste0("i", 1:4, "_t2")); d }
ond <- list(T1 = paste0("i", 1:4, "_t1"), T2 = paste0("i", 1:4, "_t2"))
il <- quieto(invariancia_longitudinal(mk(w1, w2), ond, "Afeto", estimador = "ML"))
teste("dados invariantes: nenhum passo marcado como falha", length(il$avisos) == 0)
il2 <- quieto(invariancia_longitudinal(mk(w1, w2b), ond, "Afeto", estimador = "ML"))
teste("intercepto deslocado na onda 2: o passo escalar falha", grepl("escalar", il2$avisos))
teste("resíduos do mesmo item correlacionados entre ondas no modelo configural", {
  pe <- lavaan::parameterEstimates(il$modelos$configural); any(pe$op == "~~" & pe$lhs == "i1_t1" & pe$rhs == "i1_t2") })

ok <- unlist(resultados)
cat(sprintf("\nevals_mee: %d/%d passaram\n", sum(ok), length(ok)))
if (!all(ok)) quit(status = 1)
