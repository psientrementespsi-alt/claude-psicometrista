# evals_modelos.R: testes do módulo de modelos avançados (GLM, GEE, multinível,
# meta-análise) da skill psicometria-r.
#
# Rodar da pasta da skill:  Rscript scripts/evals_modelos.R
# Uma linha PASSOU/FALHOU por teste e o resumo; código de saída 1 se algum falhar.
# Tipos de teste, nesta ordem: (1) âncoras do curso, com dados que existem em pacotes do
# R (mlmRev::Exam, metadat::dat.bcg) ou reconstruídos a partir dos resumos que o curso
# imprime; o livro do usuário não traz estes modelos, então não há âncora JASP do livro;
# (2) verdade conhecida em dados simulados com semente fixa; (3) equivalência com outra
# implementação (car, sandwich, geepack, lme4/lmerTest, performance, metafor).
# Tolerância declarada em cada teste: metade da última casa impressa na fonte, salvo
# nota em contrário. As âncoras do livro com o banco Viagra (grupo gravado como número)
# leem a biblioteca do JASP de PR_DADOS_JASP e são puladas com aviso sem ela:
#   PR_DADOS_JASP="/caminho/para/Data Library" Rscript scripts/evals_modelos.R

args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
source(file.path(pasta, "carregar.R"))  # chamado direto: carregar.R acha a própria pasta pelo frame do source()
for (p in c("lme4", "lmerTest", "geepack", "metafor", "metadat", "emmeans", "car", "MASS", "mlmRev", "sandwich"))
  if (!requireNamespace(p, quietly = TRUE)) stop("pacote necessário para os evals: ", p)

resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome))
  resultados[[nome]] <<- ok
  invisible(ok)
}
perto <- function(x, alvo, tol) isTRUE(all(abs(as.numeric(x) - alvo) <= tol + 1e-12))
DJ <- Sys.getenv("PR_DADOS_JASP")
jasp <- function(arquivo) {
  if (!nzchar(DJ)) return(NULL)
  f <- list.files(DJ, recursive = TRUE, full.names = TRUE)
  f <- f[basename(f) == arquivo]
  if (!length(f)) return(NULL)
  utils::read.csv(f[1], check.names = FALSE)
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
t_inicio <- Sys.time()

# ============================================================================ GLM
cat("== GLM\n")
# Âncora do curso (T4, 'GLM - Descrição dos Resultados'): M ± EP por grupo AP5 CONT
# 75,25 ± 18,04 (n 8), Salina 186,90 ± 16,14 (n 10), AP5 IPSI 207,50 ± 18,04 (n 8);
# Wald = 32,27; β = 112,65 e 133,25; AIC normal = 286,27. Os EPs do SPSS (escala por ML)
# dão σ²_ML entre 2603,4 e 2605,0; usamos 2604,2. Os β superam as diferenças de média em
# 1,00, o que aponta AP5 CONT = 74,25 (erro de digitação no curso): com 74,25 o Wald
# reproduz 32,27; com 75,25 daria 31,76.
monta_ratos <- function(m_cont, sig2 = 2604.2) {
  n <- c(8, 10, 8); dfs <- n - 1
  sds <- sqrt(26 * sig2 * dfs / sum(dfs) / dfs)
  z <- function(k) as.numeric(scale(seq_len(k)))
  data.frame(grupo = rep(c("AP5 CONT", "Salina", "AP5 IPSI"), n),
             congelamento = c(m_cont + z(8) * sds[1], 186.90 + z(10) * sds[2], 207.50 + z(8) * sds[3]))
}
ratos <- monta_ratos(74.25)
g1 <- modelo_glm(ratos, "congelamento", "grupo", referencia = list(grupo = "AP5 CONT"))
# tolerância: faixa de σ²_ML compatível com os EPs de 2 casas -> W em [32.27, 32.29]
teste("âncora curso: Wald χ²(2) = 32,27 com escala ML do SPSS (tol .02)", perto(g1$tabela$wald_qui2, 32.27, .02) && g1$tabela$gl == 2)
teste("âncora curso: AIC normal-identidade = 286,27 (tol .01)", perto(g1$aic, 286.27, .01))
cf <- g1$coeficientes
teste("âncora curso: β Salina = 112,65 e β AP5 IPSI = 133,25 contra AP5 CONT (tol .005)",
      perto(cf$B[cf$termo == "grupoSalina"], 112.65, .005) && perto(cf$B[cf$termo == "grupoAP5 IPSI"], 133.25, .005))
med <- g1$posthoc[g1$posthoc$tipo == "média marginal", ]
teste("âncora curso: EP das médias marginais 18,04 (n 8) e 16,14 (n 10) (tol .005)",
      perto(med$ep[med$nivel_ou_par == "AP5 CONT"], 18.04, .005) && perto(med$ep[med$nivel_ou_par == "Salina"], 16.14, .005))
g75 <- modelo_glm(monta_ratos(75.25), "congelamento", "grupo")
teste("com o 75,25 impresso no curso o Wald seria 31,76 (confirma o erro de digitação)", abs(g75$tabela$wald_qui2 - 32.27) > .4)
teste("frase APA usa Wald χ²(gl) com um só gl, nunca 'Wald (2, 23)'", grepl("Wald χ²\\(2\\) = 32\\.2[78]", g1$apa[1]) && !grepl("Wald \\(2,", g1$apa[1]))
teste("decimal = ',' troca o separador na frase", grepl("Wald χ²\\(2\\) = 32,2[78]; p < ,001", modelo_glm(ratos, "congelamento", "grupo", decimal = ",")$apa[1]))

set.seed(2026)
fat <- data.frame(A = factor(sample(c("a1", "a2", "a3"), 90, TRUE)), B = factor(sample(c("b1", "b2"), 90, TRUE)))
fat$y <- 10 + 2 * (fat$A == "a2") - 1.5 * (fat$B == "b2") + 1.2 * (fat$A == "a3") * (fat$B == "b2") + rnorm(90, 0, 2)
gp <- modelo_glm(fat, "y", c("A", "B", "A:B"), escala = "pearson")
ref <- local({ op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op))
  car::Anova(glm(y ~ A * B, data = fat), type = 3, test.statistic = "Wald") })
teste("equivalência: Wald Tipo III (escala Pearson) = car::Anova(glm, type = 3, Wald) (tol 1e-6)",
      perto(gp$tabela$wald_qui2, ref[c("A", "B", "A:B"), "Chisq"], 1e-6))
reff <- local({ op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op))
  car::Anova(lm(y ~ A * B, data = fat), type = 3) })
gm <- modelo_glm(fat, "y", c("A", "B", "A:B"))
teste("equivalência: coluna F (normal) = ANOVA Tipo III de car::Anova(lm) também com escala ML (tol 1e-6)",
      perto(gm$tabela$F, reff[c("A", "B", "A:B"), "F value"], 1e-6))
pa <- gm$posthoc[gm$posthoc$termo == "A" & grepl("^par", gm$posthoc$tipo), ]
emA <- suppressMessages(emmeans::emmeans(gm$modelo_soma, ~ A, vcov. = summary(gm$modelo_soma, dispersion = gm$escala_phi)$cov.scaled, df = Inf))
pu <- summary(emmeans::contrast(emA, "pairwise", adjust = "none"))$p.value
teste("Bonferroni: p ajustado = min(1, 3 × p sem correção) (tol 1e-8)", perto(pa$p, pmin(1, 3 * pu), 1e-8))

set.seed(11)
lg <- data.frame(x = rnorm(3000), g = factor(sample(c("c", "t"), 3000, TRUE)))
lg$y <- rbinom(3000, 1, plogis(-0.3 + log(2) * lg$x + 0.5 * (lg$g == "t")))
gl1 <- modelo_glm(lg, "y", c("x", "g"), familia = "binomial")
cx <- gl1$coeficientes[gl1$coeficientes$termo == "x", ]
teste("verdade conhecida (logística): IC do OR de x cobre 2 e coeficientes = glm() (tol 1e-8)",
      cx$exp_ic_inf < 2 && cx$exp_ic_sup > 2 && perto(gl1$coeficientes$B, coef(glm(y ~ x + g, binomial, lg)), 1e-8))
set.seed(12)
ps <- data.frame(g = factor(rep(c("a", "b"), each = 400)))
ps$y <- rpois(800, 3 * ifelse(ps$g == "b", 1.5, 1))
gpo <- modelo_glm(ps, "y", "g", familia = "poisson")
cb <- gpo$coeficientes[2, ]
disp <- gpo$pressupostos[gpo$pressupostos$checagem == "Dispersão de Pearson (χ²/gl)", ]
teste("verdade conhecida (Poisson): IRR cobre 1,5 e não acusa superdispersão (p > .01)",
      cb$exp_ic_inf < 1.5 && cb$exp_ic_sup > 1.5 && disp$p > .01)
set.seed(13)
nbd <- data.frame(g = factor(rep(c("a", "b"), each = 400)))
nbd$y <- MASS::rnegbin(800, mu = 4 * ifelse(nbd$g == "b", 1.5, 1), theta = 1.5)
gsd <- modelo_glm(nbd, "y", "g", familia = "poisson")
dsd <- gsd$pressupostos[gsd$pressupostos$checagem == "Dispersão de Pearson (χ²/gl)", ]
teste("verdade conhecida: contagem superdispersa ajustada como Poisson acusa (p < .001) e recomenda binomial negativa",
      dsd$p < .001 && any(grepl("binomial_negativa", gsd$avisos)))
gnb <- modelo_glm(nbd, "y", "g", familia = "binomial_negativa")
cnb <- gnb$coeficientes[2, ]
teste("verdade conhecida (binomial negativa): IRR cobre 1,5 e θ estimado perto de 1,5 (± .5)",
      cnb$exp_ic_inf < 1.5 && cnb$exp_ic_sup > 1.5 && abs(gnb$modelo$theta - 1.5) < .5)
set.seed(14)
gm2 <- data.frame(g = factor(rep(c("a", "b"), each = 300)), x = rnorm(600))
gm2$y <- rgamma(600, shape = 4, scale = 2 * ifelse(gm2$g == "b", 1.3, 1) / 4)
gga <- modelo_glm(gm2, "y", c("g", "x"), familia = "gama")
cg <- gga$coeficientes[gga$coeficientes$termo == "gb", ]
teste("verdade conhecida (gama, ligação log com preditor contínuo): razão de médias cobre 1,3 e φ_ML perto de 1/4 (± .04)",
      gga$ligacao == "log" && cg$exp_ic_inf < 1.3 && cg$exp_ic_sup > 1.3 && abs(gga$escala_phi - .25) < .04)
cfm <- comparar_familias(gm2, "y", c("g", "x"), familias = c("normal", "gama"), ligacoes = c("identity", "log"))
teste("comparar_familias: gama (verdadeira) com AIC menor que normal no mesmo Y", cfm$tabela$familia[which.min(cfm$tabela$AIC)] == "gama")
cmix <- comparar_familias(nbd, "y", "g", familias = c("normal", "poisson"))
teste("comparar_familias avisa que AIC entre contínua e discreta não é comparável", any(grepl("não é comparável", cmix$avisos)))
teste("gama recusa desfecho com zero e explica que somar constante muda o modelo",
      { r <- tryCatch(modelo_glm(data.frame(y = c(0, 1, 2, 3), g = c("a", "a", "b", "b")), "y", "g", familia = "gama"), error = conditionMessage)
        is.character(r) && grepl("constante muda o modelo", r) })

cat("== Grupo gravado como número (regressão da revisão CoVe)\n")
set.seed(61)
gn <- data.frame(grupo = rep(1:3, each = 30)); gn$y <- rpois(90, c(2, 6, 2)[gn$grupo])
g_num <- modelo_glm(gn, "y", "grupo", familia = "poisson")
g_fat <- modelo_glm(gn, "y", "grupo", familia = "poisson", fatores = "grupo")
g_ref <- modelo_glm(transform(gn, grupo = factor(grupo)), "y", "grupo", familia = "poisson")
teste("GLM: grupo 1/2/3 numérico avisa que entra como contínuo (1 gl) e sugere fatores =",
      g_num$tabela$gl == 1 && any(grepl('fatores = "grupo"', g_num$avisos, fixed = TRUE)))
teste("GLM: com fatores = 'grupo', Wald χ²(2) igual ao de factor(grupo) (tol 1e-8) e sem o aviso",
      g_fat$tabela$gl == 2 && perto(g_fat$tabela$wald_qui2, g_ref$tabela$wald_qui2, 1e-8) && !any(grepl("entra como contínua", g_fat$avisos)))
teste("GLM: decisões só falam em categoria de referência quando há fator", !any(grepl("categoria de referência", g_num$decisoes)) &&
        any(grepl("categoria de referência", g_fat$decisoes)))
v <- jasp("Viagra.csv")
if (is.null(v)) cat("PULADO  âncoras do livro com o banco Viagra (PR_DADOS_JASP vazia ou banco não encontrado)\n") else {
  vn <- modelo_glm(v, "libido", "dose")
  vf <- modelo_glm(v, "libido", "dose", fatores = "dose")
  teste("livro: dose gravada como número entra como contínua e a função avisa", vn$tabela$gl == 1 && any(grepl('fatores = "dose"', vn$avisos, fixed = TRUE)))
  teste("âncora livro: libido por dose com fatores = 'dose', F(2, 27) = 2,416, p = ,108 (tol .0005)",
        vf$tabela$gl == 2 && vf$tabela$gl_res == 27 && perto(vf$tabela$F, 2.416, .0005) && perto(vf$tabela$p_F, .108, .0005))
  va <- modelo_glm(v, "libido", c("dose", "partnerLibido"), fatores = "dose")
  teste("âncora livro: ANCOVA dose F(2, 26) = 4,142 e covariável F(1, 26) = 4,959 (tol .0005)",
        perto(va$tabela$F[va$tabela$termo == "dose"], 4.142, .0005) && perto(va$tabela$F[va$tabela$termo == "partnerLibido"], 4.959, .0005))
  # H-10 (verificação final): references/modelos.md compara as duas codificações na mesma
  # métrica. Conta independente com lm(): Wald com escala ML = n (SQR0 - SQR1) / SQR1 e F da
  # anova() de modelos aninhados. Dose numérica: Wald χ²(1) = 3,89 (p = .049), F(1, 28) =
  # 3,63 (p = .067); dose como fator: Wald χ²(2) = 5,37 (p = .068), F(2, 27) = 2,416 (p = .108).
  teste("H-10: Viagra na mesma métrica, Wald χ² (3,89, p .049 contra 5,37, p .068) e F (3,63, p .067 contra 2,416, p .108) = lm() (tol 1e-8)", {
    nv <- nrow(v); s0 <- sum(resid(lm(libido ~ 1, v))^2); s1 <- sum(resid(lm(libido ~ dose, v))^2); s2 <- sum(resid(lm(libido ~ factor(dose), v))^2)
    a1 <- anova(lm(libido ~ 1, v), lm(libido ~ dose, v)); a2 <- anova(lm(libido ~ 1, v), lm(libido ~ factor(dose), v))
    perto(vn$tabela$wald_qui2, nv * (s0 - s1) / s1, 1e-8) && perto(vf$tabela$wald_qui2, nv * (s0 - s2) / s2, 1e-8) &&
      perto(c(vn$tabela$F, vn$tabela$p_F), c(a1$F[2], a1$`Pr(>F)`[2]), 1e-8) && perto(c(vf$tabela$F, vf$tabela$p_F), c(a2$F[2], a2$`Pr(>F)`[2]), 1e-8) &&
      perto(c(vn$tabela$wald_qui2, vn$tabela$F, vf$tabela$wald_qui2), c(3.89, 3.63, 5.37), .005) &&
      perto(c(vn$tabela$p, vn$tabela$p_F, vf$tabela$p), c(.049, .067, .068), .0005) })
}

# ============================================================================ GEE
cat("== GEE\n")
sim_long <- function(K, T, rho, seed, falta = 0) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(K), each = T), tempo = rep(seq_len(T), K))
  d$grupo <- factor(ifelse(d$id <= K / 2, "A", "B"))
  u <- rnorm(K, 0, sqrt(rho)); d$y <- 10 + 0.5 * (d$grupo == "B") + 0.3 * d$tempo + u[d$id] + rnorm(K * T, 0, sqrt(1 - rho))
  if (falta > 0) d$y[sample(nrow(d), falta)] <- NA
  d
}
dl <- sim_long(60, 4, .5, 21, falta = 12)
gi <- modelo_gee(dl, "y", c("grupo", "tempo"), id = "id", tempo = "tempo", correlacao = "independence")
dlc <- dl[!is.na(dl$y), ]; dlc$tempo <- factor(dlc$tempo)
gref <- glm(y ~ grupo + tempo, data = dlc)
teste("equivalência: GEE independente = glm() nos coeficientes (tol 1e-6)", perto(gi$coeficientes$B, coef(gref), 1e-6))
V0 <- sandwich::vcovCL(gref, cluster = ~id, type = "HC0", cadjust = FALSE)
teste("equivalência: EP robusto (independente) = sandwich::vcovCL HC0 por sujeito (tol rel. 1e-6)", rel(gi$vcov$robusto, V0) < 1e-6)
V3 <- sandwich::vcovCL(gref, cluster = ~id, type = "HC3", cadjust = FALSE)
G <- length(unique(dlc$id))
# o HC3 por cluster do sandwich é a forma jackknife, (G - 1)/G vezes a fórmula de Mancl e DeRouen
teste("equivalência: Mancl e DeRouen (independente, normal) × (G - 1)/G = sandwich::vcovCL HC3 (tol rel. 1e-6)",
      rel(gi$vcov$mancl_derouen * (G - 1) / G, V3) < 1e-6)
reproduz <- vapply(c("exchangeable", "ar1", "unstructured"), function(cs) {
  g <- modelo_gee(dl, "y", c("grupo", "tempo"), id = "id", tempo = "tempo", correlacao = cs)
  ft <- g$modelo; ok <- !is.na(coef(ft))
  g$correlacao == cs && rel(g$vcov$robusto, ft$geese$vbeta[ok, ok]) < 1e-6
}, logical(1))
teste("equivalência: sanduíche recalculado = geepack (intercambiável, AR(1) e não estruturada, com tempos faltando) (tol rel. 1e-6)", all(reproduz))
gu <- modelo_gee(dl, "y", c("grupo", "tempo"), id = "id", tempo = "tempo", correlacao = "unstructured")
teste("não estruturada com tempos faltando no meio roda (sem segfault) e estima 6 correlações para 4 tempos",
      length(gu$modelo$geese$alpha) == 6 && all(abs(gu$modelo$geese$alpha) < 1))
dg <- sim_long(400, 4, .5, 22)
ge <- modelo_gee(dg, "y", c("grupo", "tempo"), id = "id", tempo = "tempo", correlacao = "exchangeable")
# tolerância = 3 EP de Monte Carlo (EP do alfa com 400 sujeitos e 4 tempos ≈ .023)
teste("verdade conhecida: correlação intercambiável estimada perto de .50 (± .07)", abs(ge$modelo$geese$alpha - .5) < .07)
ce <- ge$coeficientes
teste("verdade conhecida: efeito de grupo (0,5) e de tempo 4 vs 1 (0,9) dentro do IC 95%",
      ce$ic_inf[ce$termo == "grupoB"] < .5 && ce$ic_sup[ce$termo == "grupoB"] > .5 &&
        ce$ic_inf[ce$termo == "tempo4"] < .9 && ce$ic_sup[ce$termo == "tempo4"] > .9)
teste("frase APA do GEE com Wald χ²(gl) e sem a forma 'Wald (1, 24)'", grepl("tempo, Wald χ²\\(3\\) = ", ge$apa[1]) && !grepl("Wald \\(", ge$apa[1]))
# com desenho balanceado e preditores só entre sujeitos ou balanceados no tempo, o GEE
# intercambiável coincide com o MQO; a diferença de eficiência aparece com covariável
# que varia no tempo e entre sujeitos
set.seed(27)
dx <- sim_long(150, 4, .6, 27); dx$x <- rnorm(150)[dx$id] + rnorm(600); dx$y <- dx$y + .4 * dx$x
ccg <- comparar_correlacoes_gee(dx, "y", c("grupo", "x"), id = "id", tempo = "tempo")
teste("comparar_correlacoes_gee: com covariável que varia no tempo, QIC e CIC da intercambiável (verdadeira) menores que os da independência",
      ccg$tabela$CIC[ccg$tabela$correlacao == "exchangeable"] < ccg$tabela$CIC[ccg$tabela$correlacao == "independence"] &&
        ccg$tabela$QIC[ccg$tabela$correlacao == "exchangeable"] < ccg$tabela$QIC[ccg$tabela$correlacao == "independence"])
d2 <- sim_long(80, 2, .5, 23)
c2 <- comparar_correlacoes_gee(d2, "y", "grupo", id = "id", tempo = "tempo", correlacoes = c("exchangeable", "ar1", "unstructured"))
teste("2 tempos: intercambiável, AR(1) e não estruturada dão o mesmo QIC (tol 1e-3) e a função avisa",
      diff(range(c2$tabela$QIC)) < 1e-3 && any(grepl("2 tempos", c2$avisos)))
# QIC: equivalência com geepack::QIC quando a escala é fixa em 1 (Poisson)
set.seed(24)
dp <- data.frame(.id = rep(1:50, each = 3), .onda = rep(1:3, 50)); dp$x <- rnorm(150)
dp$y <- rpois(150, exp(.5 + .3 * dp$x + rep(rnorm(50, 0, .4), each = 3)))
fq <- geepack::geeglm(y ~ x, family = poisson, data = dp, id = .id, corstr = "exchangeable", scale.fix = TRUE)
fqi <- geepack::geeglm(y ~ x, family = poisson, data = dp, id = .id, corstr = "independence", scale.fix = TRUE)
teste("equivalência: QIC (Poisson, escala fixa em 1) = geepack::QIC (tol 1e-6)",
      perto(.ee_mod_qic(fq, fqi)[["QIC"]], geepack::QIC(fq)[["QIC"]], 1e-6))
# informativo: a quasi-verossimilhança da gama no geepack 1.3.13
dgm <- dp; dgm$y <- rgamma(150, shape = 3, scale = exp(.5 + .3 * dgm$x) / 3)
fgm <- geepack::geeglm(y ~ x, family = Gamma("log"), data = dgm, id = .id, corstr = "independence")
qg <- geepack::QIC(fgm)[["Quasi Lik"]]; qc <- sum(-fgm$y / fitted(fgm) - log(fitted(fgm)))
cat(sprintf("  info: quasi-verossimilhança da gama, geepack::QIC = %.3f; -y/mu - log(mu) = %.3f\n", qg, qc))
# poucos clusters: erro tipo I do Wald com EP robusto vs Mancl-DeRouen (K = 12, preditor de nível do sujeito)
set.seed(25)
rej <- replicate(200, {
  K <- 12; d <- data.frame(.id = rep(1:K, each = 4), .onda = rep(1:4, K))
  d$g <- rep(c(-1, 1), each = K * 2); d$y <- rep(rnorm(K), each = 4) + rnorm(K * 4)
  f <- geepack::geeglm(y ~ g, data = d, id = .id, corstr = "exchangeable")
  V <- .ee_mod_gee_vcov(f, d, "exchangeable")
  b <- coef(f)[2]
  c(robusto = abs(b) / sqrt(V$robusto[2, 2]) > 1.96, md = abs(b) / sqrt(V$mancl_derouen[2, 2]) > 1.96)
})
taxa <- rowMeans(rej)
cat(sprintf("  info: erro tipo I nominal .05 com 12 clusters: robusto %.3f; Mancl-DeRouen %.3f\n", taxa[1], taxa[2]))
teste("verdade conhecida: com 12 clusters o EP robusto infla o erro tipo I e Mancl-DeRouen reduz", taxa[2] < taxa[1] && taxa[1] > .07)
gk <- modelo_gee(sim_long(12, 4, .5, 26), "y", c("grupo", "tempo"), id = "id", tempo = "tempo", correlacao = "exchangeable")
teste("com menos de 40 sujeitos a função avisa e traz o p com Mancl-DeRouen ao lado",
      any(grepl("Mancl e DeRouen", gk$avisos)) && "p_mancl_derouen" %in% names(gk$tabela))
teste("GEE avisa MCAR quando há linhas excluídas", any(grepl("MCAR", gi$avisos)))
dl$g3 <- (dl$id - 1) %% 3 + 1
gq_n <- modelo_gee(dl, "y", c("g3", "tempo"), id = "id", tempo = "tempo", correlacao = "exchangeable")
gq_f <- modelo_gee(dl, "y", c("g3", "tempo"), id = "id", tempo = "tempo", correlacao = "exchangeable", fatores = "g3")
teste("GEE: grupo numérico avisa (tempo não, que já vira fator); com fatores = 'g3' o termo tem 2 gl",
      any(grepl('fatores = "g3"', gq_n$avisos, fixed = TRUE)) && !any(grepl('fatores = "tempo"', gq_n$avisos, fixed = TRUE)) &&
        gq_n$tabela$gl[gq_n$tabela$termo == "g3"] == 1 && gq_f$tabela$gl[gq_f$tabela$termo == "g3"] == 2)

# ============================================================================ multinível
cat("== Multinível linear (âncora do curso: mlmRev::Exam)\n")
data(Exam, package = "mlmRev")
ex <- multinivel_linear(Exam, "normexam", "school", nivel1 = c("standLRT", "sex"), nivel2 = c("schavg", "schgend"),
                        inclinacao = "standLRT", interacao = "schgend", referencias = list(sex = "F", schgend = "mixed"),
                        estimador_final = "ML")
cp <- ex$posthoc
teste("âncora curso (reporte): ICC do nulo em ML = .166 (tol .0005)", perto(ex$icc_ml, .166, .0005))
teste("âncora curso (prática no JASP): ICC em REML = .168, τ00 = .172, σ² = .848 (tol .0005)",
      perto(ex$icc_reml, .168, .0005) && local({ v <- as.data.frame(lme4::VarCorr(lme4::lmer(normexam ~ 1 + (1 | school), Exam)))
        perto(v$vcov[1], .172, .0005) && perto(v$vcov[2], .848, .0005) }))
fx <- function(m) lme4::fixef(ex$modelos_ml[[m]])
teste("âncora curso: M1 standLRT = .559 (tol .001; o curso trunca .5595) e sexo M = -.171 (tol .0005)",
      perto(fx("M1")[["standLRT"]], .559, .001) && perto(fx("M1")[["sexM"]], -.171, .0005))
teste("âncora curso: Δ-2LL M1 vs M0 = 1680,6 com Δgl = 2 (tol .05)", perto(cp$delta_desvio[2], 1680.6, .05) && cp$delta_gl[2] == 2)
teste("âncora curso: M2 schavg = .348, boys = .194, girls = .151 (tol .0005); Δ-2LL = 14,7 (tol .05), Δgl = 3",
      perto(fx("M2")[c("schavg", "schgendboys", "schgendgirls")], c(.348, .194, .151), .0005) && perto(cp$delta_desvio[3], 14.7, .05) && cp$delta_gl[3] == 3)
v3 <- as.data.frame(lme4::VarCorr(ex$modelos_ml$M3))
teste("âncora curso: M3 Var(u1j) = .015, DP = .122 (tol .0005), r intercepto-inclinação = .48 (tol .005)",
      perto(v3$vcov[2], .015, .0005) && perto(v3$sdcor[2], .122, .0005) && perto(v3$sdcor[3], .48, .005))
fe <- ex$efeito
teste("âncora curso: inclinação média .549 (tol .001) e faixa de 68% de .427 a .672 (tol .0005)",
      perto(fx("M3")[["standLRT"]], .549, .001) &&
        perto(fe$ic_inf[grepl("68%", fe$indice)], .427, .0005) && perto(fe$ic_sup[grepl("68%", fe$indice)], .672, .0005))
# Não reproduzido na última casa: o reporte do curso dá Δ-2LL = 40,199; lme4 2.0-6 dá 40,181
# com três otimizadores e com contr.sum. Tolerância declarada de .03 (ver references/modelos.md).
teste("âncora curso (com ressalva): Δ-2LL M3 vs M2 = 40,2 (curso 40,199; R 40,18; tol .03), Δgl = 2",
      perto(cp$delta_desvio[4], 40.199, .03) && cp$delta_gl[4] == 2)
teste("âncora curso: M4 interações boys×standLRT = -.016 e girls×standLRT = -.021 (tol .0005); Δ-2LL = 0,2 (tol .05), p = .89 (tol .005)",
      perto(fx("M4")[c("standLRT:schgendboys", "standLRT:schgendgirls")], c(-.016, -.021), .0005) &&
        perto(cp$delta_desvio[5], .2, .05) && perto(cp$p_lrt[5], .89, .005))
teste("modelo final automático = M3 (M4 não melhora), como no curso", ex$modelo_final == "M3")
j1 <- multinivel_linear(Exam, "normexam", "school", nivel1 = "standLRT", nivel2 = "schavg")
teste("âncora curso (prática no JASP): standLRT sem sexo = .563 (tol .0005)", perto(lme4::fixef(j1$modelos_ml$M1)[["standLRT"]], .563, .0005))
j2 <- multinivel_linear(Exam, "normexam", "school", nivel1 = c("standLRT", "sex"), nivel2 = "schavg", referencias = list(sex = "F"))
teste("âncora curso (prática no JASP): schavg = .343 com standLRT e sexo (tol .0005)", perto(lme4::fixef(j2$modelos_ml$M2)[["schavg"]], .343, .0005))
teste("desvio = -2 logLik (positivo; o checklist do curso tem erro de sinal) e = deviance() do lme4 (tol 1e-6)",
      all(cp$desvio > 0) && perto(cp$desvio, vapply(ex$modelos_ml, function(m) deviance(m), numeric(1)), 1e-6))
an <- anova(ex$modelos_ml$M2, ex$modelos_ml$M3)
teste("equivalência: p do LRT = anova() do lme4 (tol 1e-8) e mistura de χ² menor que o p ingênuo",
      perto(cp$p_lrt[4], an$`Pr(>Chisq)`[2], 1e-8) &&
        perto(cp$p_mistura[4], .5 * pchisq(cp$delta_desvio[4], 1, lower.tail = FALSE) + .5 * pchisq(cp$delta_desvio[4], 2, lower.tail = FALSE), 1e-12) &&
        cp$p_mistura[4] < cp$p_lrt[4])
if (requireNamespace("performance", quietly = TRUE)) {
  r2p1 <- suppressWarnings(performance::r2_nakagawa(ex$modelos_ml$M1)); r2p3 <- suppressWarnings(performance::r2_nakagawa(ex$modelos_ml$M3))
  teste("equivalência: R² de Nakagawa = performance::r2_nakagawa, intercepto aleatório e inclinação aleatória (tol 1e-4)",
        perto(c(cp$r2_marginal[2], cp$r2_condicional[2]), c(r2p1$R2_marginal, r2p1$R2_conditional), 1e-4) &&
          perto(c(cp$r2_marginal[4], cp$r2_condicional[4]), c(r2p3$R2_marginal, r2p3$R2_conditional), 1e-4))
  teste("equivalência: ICC do nulo = performance::icc (tol 1e-8)",
        perto(ex$icc_ml, performance::icc(ex$modelos_ml$M0)$ICC_adjusted, 1e-8))
} else cat("  aviso: performance ausente, testes de equivalência do R² pulados\n")
exr <- multinivel_linear(Exam, "normexam", "school", nivel1 = c("standLRT", "sex"), nivel2 = c("schavg", "schgend"),
                         inclinacao = "standLRT", referencias = list(sex = "F", schgend = "mixed"))
teste("estimativas finais em REML com gl de Satterthwaite = lmerTest (tol 1e-8)",
      local({ m <- lmerTest::lmer(normexam ~ standLRT + sex + schavg + schgend + (1 + standLRT | school), Exam, REML = TRUE)
        s <- summary(m)$coefficients; perto(exr$tabela$b, s[, 1], 1e-6) && perto(exr$tabela$gl, s[, "df"], 1e-4) }))
teste("decisões registram a inversão M2/M3 da página 'Ajuste' e o erro de sinal do desvio",
      any(grepl("inverte M2 e M3", ex$decisoes)) && any(grepl("erro de sinal", ex$decisoes)))
ex2 <- Exam; ex2$schgend_n <- as.integer(ex2$schgend)  # 1 = mixed, 2 = boys, 3 = girls
mn <- multinivel_linear(ex2, "normexam", "school", nivel1 = c("standLRT", "sex"), nivel2 = c("schavg", "schgend_n"), referencias = list(sex = "F"))
mfz <- multinivel_linear(ex2, "normexam", "school", nivel1 = c("standLRT", "sex"), nivel2 = c("schavg", "schgend_n"),
                         referencias = list(sex = "F"), fatores = "schgend_n")
teste("multinível: nível 2 gravado como número avisa; com fatores reproduz boys .194 e girls .151 do curso (tol .0005)",
      any(grepl('fatores = "schgend_n"', mn$avisos, fixed = TRUE)) &&
        perto(lme4::fixef(mfz$modelos_ml$M2)[c("schgend_n2", "schgend_n3")], c(.194, .151), .0005))
teste("multinível: inclinação aleatória de preditor categórico para na entrada com mensagem clara",
      { r <- tryCatch(multinivel_linear(Exam, "normexam", "school", nivel1 = c("standLRT", "sex"), inclinacao = "sex",
                                        referencias = list(sex = "F")), error = conditionMessage)
        is.character(r) && grepl("categórico", r) })

cat("== Multinível: verdade conhecida e centralização\n")
set.seed(31)
J <- 150; m <- 20
s0 <- data.frame(cl = rep(1:J, each = m)); s0$y <- rnorm(J, 0, sqrt(.25))[s0$cl] + rnorm(J * m, 0, sqrt(.75))
n0 <- multinivel_linear(s0, "y", "cl")
# tolerância = 3 EP de Monte Carlo (EP do ICC com 150 clusters de 20 ≈ .026)
teste("verdade conhecida: ICC = .25 recuperado (± .08) e efeito de desenho = 1 + 19 × ICC", abs(n0$icc_reml - .25) < .08 &&
        perto(n0$efeito$valor[n0$efeito$indice == "Efeito de desenho (M0)"], 1 + 19 * n0$icc_reml, 1e-10))
set.seed(32)
J <- 120; m <- 20
s1 <- data.frame(cl = rep(1:J, each = m)); w <- rnorm(J); u0 <- rnorm(J, 0, .5); u1 <- rnorm(J, 0, .2)
s1$w <- w[s1$cl]; s1$x <- rnorm(J * m) + rnorm(J)[s1$cl]
s1 <- centralizar_multinivel(s1, "x", "cl", tipo = "grupo")
s1$y <- .3 * s1$w + (.5 + u1[s1$cl]) * s1$x_cwc + .2 * s1$x_media_cgm + u0[s1$cl] + rnorm(J * m)
teste("centralizar_multinivel: X_cwc com média zero em cada cluster e média do cluster constante no cluster",
      max(abs(tapply(s1$x_cwc, s1$cl, mean))) < 1e-10 && all(tapply(s1$x_media_cgm, s1$cl, function(z) diff(range(z))) == 0))
teste("centralizar_multinivel recusa centralizar no grupo uma variável de nível 2",
      { r <- tryCatch(centralizar_multinivel(s1, "w", "cl", tipo = "grupo"), error = conditionMessage); is.character(r) && grepl("nível 2", r) })
teste("centralizar_multinivel na média geral: nível 2 centrado na média entre clusters",
      { s2 <- suppressMessages(centralizar_multinivel(s1, "w", "cl", tipo = "geral")); abs(mean(tapply(s2$w_cgm, s2$cl, mean))) < 1e-10 })
cru <- data.frame(cl = rep(1:30, each = 5), x = rnorm(150), w = rep(rnorm(30), each = 5))
saida <- capture.output(centralizar_multinivel(cru, c("x", "w"), "cl"))
cz <- centralizar_multinivel(cru, c("x", "w"), "cl")
teste("centralizar_multinivel: padrão = média geral (curso), com a alternativa CWC criada ao lado para o nível 1",
      all(c("x_cgm", "w_cgm", "x_cwc", "x_media_cgm") %in% names(cz)) && !"w_cwc" %in% names(cz) && abs(mean(cz$x_cgm)) < 1e-10)
pz <- capture.output(print(cz)); ph <- capture.output(print(head(cz)))
teste("centralizar_multinivel: devolve invisível e o print mostra só agregados (nenhuma linha de participante)",
      length(saida) == 0 && length(pz) < 15 && grepl("Centralização multinível", pz[1]) && length(ph) == 1 &&
        !any(grepl(format(signif(cru$x[1], 4)), c(pz, ph), fixed = TRUE)) && is.data.frame(cz))
teste("centralizar_multinivel: o objeto serve de banco para multinivel_linear", {
  cz$y <- cz$x_cwc + cz$w + rnorm(150); r <- multinivel_linear(cz, "y", "cl", nivel1 = "x_cwc"); r$n == 150 })
ks <- multinivel_linear(s1, "y", "cl", nivel1 = "x_cwc", nivel2 = c("w", "x_media_cgm"), inclinacao = "x_cwc")
fk <- ks$tabela
teste("verdade conhecida: γ10 = .5 e γ01 = .3 dentro do IC 95% do modelo final (M3)",
      ks$modelo_final == "M3" && fk$ic_inf[fk$termo == "x_cwc"] < .5 && fk$ic_sup[fk$termo == "x_cwc"] > .5 &&
        fk$ic_inf[fk$termo == "w"] < .3 && fk$ic_sup[fk$termo == "w"] > .3)
teste("verdade conhecida: variância da inclinação (.04) recuperada entre .015 e .08",
      { t11 <- ks$efeito$valor[grepl("τ11", ks$efeito$indice)]; t11 > .015 && t11 < .08 })
teste("multinivel_linear recusa preditor de nível 2 que varia dentro do cluster",
      { r <- tryCatch(multinivel_linear(s1, "y", "cl", nivel2 = "x_cwc"), error = conditionMessage); is.character(r) && grepl("nível 1", r) })

cat("== Multinível logístico\n")
set.seed(41)
J <- 120; m <- 30
lb <- data.frame(cl = rep(1:J, each = m)); lb$y <- rbinom(J * m, 1, plogis(-.3 + rnorm(J, 0, 1)[lb$cl]))
l0 <- multinivel_logistico(lb, "y", "cl")
teste("verdade conhecida: ICC latente = 1/(1 + π²/3) = .233 recuperado (± .05)", abs(l0$icc_latente - 1 / (1 + pi^2 / 3)) < .05)
set.seed(42)
lb2 <- data.frame(cl = rep(1:J, each = m)); lb2$x <- rnorm(J * m); lb2$w <- rnorm(J)[lb2$cl]
lb2$y <- rbinom(J * m, 1, plogis(-.3 + .7 * lb2$x + .4 * lb2$w + rnorm(J, 0, 1)[lb2$cl]))
l2 <- multinivel_logistico(lb2, "y", "cl", nivel1 = "x", nivel2 = "w")
ox <- l2$tabela[l2$tabela$termo == "x", ]
teste("verdade conhecida: OR condicional de x = exp(.7) dentro do IC 95%", ox$OR_ic_inf < exp(.7) && ox$OR_ic_sup > exp(.7))
if (requireNamespace("performance", quietly = TRUE)) {
  r2l <- suppressWarnings(performance::r2_nakagawa(l2$modelo))
  teste("equivalência: R² de Nakagawa (logística, π²/3) = performance::r2_nakagawa (tol 1e-4)",
        perto(c(l2$posthoc$r2_marginal[3], l2$posthoc$r2_condicional[3]), c(r2l$R2_marginal, r2l$R2_conditional), 1e-4))
}
teste("logística multinível declara OR condicionais ao cluster", any(grepl("condicionais ao cluster", l2$decisoes)))
lb2$g3n <- (lb2$cl %% 3) + 1
l2n <- multinivel_logistico(lb2, "y", "cl", nivel2 = "g3n")
teste("logística multinível: nível 2 gravado como número (1/2/3) avisa e sugere fatores =", any(grepl('fatores = "g3n"', l2n$avisos, fixed = TRUE)))
# Regra única dos módulos (associação, processos, modelos): 2 valores não avisam, porque a
# codificação numérica da dicotômica é a que o livro manda.
lb2$g2 <- (lb2$cl %% 2) + 1
l2b <- multinivel_logistico(lb2, "y", "cl", nivel2 = "g2")
teste("logística multinível: preditor de 2 valores (1/2) não avisa (mesma regra da associação)", !any(grepl('fatores = "g2"', l2b$avisos, fixed = TRUE)))
Exam$alto <- as.integer(Exam$normexam > 0)
la <- multinivel_logistico(Exam, "alto", "school", nivel1 = "standLRT", nivel2 = "schavg", inclinacao = "standLRT", nAGQ = 10)
m2l <- lme4::glmer(alto ~ standLRT + schavg + (1 | school), Exam, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa"))
teste("nAGQ > 1: o LRT do M3 usa o M2 reajustado por Laplace (Δ-2LL = 3,52; tol 1e-3) e a decisão registra",
      perto(la$posthoc$delta_desvio[4], -2 * as.numeric(logLik(m2l)) + 2 * as.numeric(logLik(la$modelos$M3)), 1e-3) &&
        perto(la$posthoc$delta_desvio[4], 3.52, .005) && any(grepl("Laplace", la$decisoes)))

# ============================================================================ meta-análise
cat("== Meta-análise (âncoras do curso)\n")
# Curso, 'Modelo de Efeitos Aleatórios' e 'Estimadores de Variância Tau²': DM de mudança
# de 3 estudos, com N, médias e DPs do forest plot da aula 'Anatomia do Forest Plot'.
cur <- data.frame(estudo = c("Cho 2012", "Fritz 2013", "Morone 2014"), m1 = c(4.00, 2.80, 7.50), d1 = c(1.18, 1.84, 1.01), n1 = c(11, 15, 25),
                  m2 = c(2.81, 0.31, 5.66), d2 = c(0.40, 0.81, 2.29), n2 = c(11, 13, 22))
ec <- calcular_efeitos_meta(cur, "MD", estudo = "estudo", media1 = "m1", dp1 = "d1", n1 = "n1", media2 = "m2", dp2 = "d2", n2 = "n2")
teste("âncora curso: DM 1,19 / 2,49 / 1,84 e variâncias .141 / .276 / .279 (tol .0005)",
      perto(ec$yi, c(1.19, 2.49, 1.84), .0005) && perto(ec$vi, c(.141, .276, .279), .0005))
ma <- meta_analise(ec, sensibilidade = c("HS", "ML", "HE", "PM", "EB", "DL", "SJ"))
l1 <- ma$tabela[1, ]
teste("âncora curso: REML DM = 1,7713 [1,0013; 2,5414], τ² = .2407, z = 4,51 (tol meia última casa)",
      perto(l1$estimativa, 1.7713, 5e-5) && perto(c(l1$ic_inf, l1$ic_sup), c(1.0013, 2.5414), 5e-5) && perto(l1$tau2, .2407, 5e-5) && perto(l1$estatistica, 4.51, .005))
hq <- ma$pressupostos
teste("âncora curso: Q = 4,17 (gl 2, p = .1244), I² = 52,0% pela fórmula (Q - gl)/Q (tol meia última casa)",
      perto(hq$estatistica[1], 4.17, .005) && hq$gl[1] == 2 && perto(hq$p[1], .1244, 5e-5) && perto(hq$estatistica[2], 52.0, .05))
lin <- function(rot) ma$tabela[ma$tabela$modelo == paste0("Sensibilidade: ", rot), ]
alvo <- list(HS = c(.0818, 1.7305), ML = c(.0934, 1.7351), HE = c(.1903, 1.7621), PM = c(.2171, 1.7673), EB = c(.2171, 1.7673),
             DL = c(.2420, 1.7715), SJ = c(.2441, 1.7719))
teste("âncora curso: τ² e DM pelos 8 estimadores (HS, ML, HE, PM, EB, REML, DL, SJ) (tol 5e-5)",
      all(vapply(names(alvo), function(m) perto(c(lin(m)$tau2, lin(m)$estimativa), alvo[[m]], 5e-5), logical(1))))
teste("âncora curso: IC de PM [1,0176; 2,5170], DL [1,0003; 2,5428] e SJ [0,9989; 2,5449] (tol 5e-5)",
      perto(c(lin("PM")$ic_inf, lin("PM")$ic_sup, lin("DL")$ic_inf, lin("DL")$ic_sup, lin("SJ")$ic_inf, lin("SJ")$ic_sup),
            c(1.0176, 2.5170, 1.0003, 2.5428, .9989, 2.5449), 5e-5))
ce1 <- lin("efeito comum")
teste("âncora curso: efeito comum DM = 1,68 [1,16; 2,20], z = 6,36 (tol .005)",
      perto(c(ce1$estimativa, ce1$ic_inf, ce1$ic_sup, ce1$estatistica), c(1.68, 1.16, 2.20, 6.36), .005))
es <- ma$estudos
teste("âncora curso: IC por estudo (Cho [0,45; 1,93]) e pesos comum 49,6/25,3/25,1 e aleatório 40,4/29,9/29,7 (tol .05)",
      perto(c(es$ic_inf[1], es$ic_sup[1]), c(.45, 1.93), .005) && perto(es$peso_comum, c(49.6, 25.3, 25.1), .05) &&
        perto(es$peso_aleatorio, c(40.4, 29.9, 29.7), .05))
teste("com 3 estudos: avisa τ² impreciso e não roda testes de assimetria", any(grepl("Só 3 estudos", ma$avisos)) && any(grepl("não rodam", ma$avisos)))

# BCG (metadat::dat.bcg = banco das aulas de subgrupos e meta-regressão; RR, REML)
data("dat.bcg", package = "metadat")
bcg <- dat.bcg; bcg$estudo <- paste(bcg$author, bcg$year); bcg$n1 <- bcg$tpos + bcg$tneg; bcg$n2 <- bcg$cpos + bcg$cneg
eb <- calcular_efeitos_meta(bcg, "RR", estudo = "estudo", eventos1 = "tpos", n1 = "n1", eventos2 = "cpos", n2 = "n2")
mb <- meta_analise(eb)
teste("âncora curso (BCG): RR = .49 [.34; .70], τ² = .3132, I² = 92,1% (tol meia última casa)",
      perto(mb$efeito$valor[1], .49, .005) && perto(c(mb$efeito$ic_inf[1], mb$efeito$ic_sup[1]), c(.34, .70), .005) &&
        perto(mb$tabela$tau2[1], .3132, 5e-5) && perto(mb$i2_q, 92.1, .05))
teste("equivalência (documentação do metafor, dat.bcg): log RR = -.7145, EP = .1798, Q = 152,2330, I² (τ²) = 92,22% (tol meia última casa)",
      perto(mb$tabela$estimativa[1], -.7145, 5e-5) && perto(mb$tabela$ep[1], .1798, 5e-5) && perto(mb$pressupostos$estatistica[1], 152.2330, 5e-5) &&
        perto(mb$pressupostos$estatistica[3], 92.22, .005))
sb <- meta_subgrupos(eb, "alloc")
st <- sb$tabela
alvo_sg <- rbind(random = c(.38, .22, .65, 94.6), alternate = c(.58, .34, 1.01, 82.0), systematic = c(.65, .32, 1.32, 81.9))
teste("âncora curso: subgrupos de alocação (RR e IC, tol .005; I², tol .05)",
      all(vapply(rownames(alvo_sg), function(g) { l <- st[st$modelo == g, ]
        perto(c(l$efeito, l$efeito_ic_inf, l$efeito_ic_sup), alvo_sg[g, 1:3], .005) && perto(l$i2_q, alvo_sg[g, 4], .05) }, logical(1))))
teste("âncora curso: teste entre subgrupos χ² = 1,86, gl = 2, p = .3943 (tol meia última casa)",
      perto(sb$pressupostos$estatistica[1], 1.86, .005) && sb$pressupostos$gl[1] == 2 && perto(sb$pressupostos$p[1], .3943, 5e-5))
mr <- meta_regressao(eb, "ablat")
mt <- mr$tabela[2, ]; mp <- mr$pressupostos
teste("âncora curso: meta-regressão em latitude B = -.0291, EP .0072, z = -4,0444, IC [-.0432; -.0150] (tol meia última casa)",
      perto(mt$B, -.0291, 5e-5) && perto(mt$ep, .0072, 5e-5) && perto(mt$estatistica, -4.0444, 5e-5) && perto(c(mt$ic_inf, mt$ic_sup), c(-.0432, -.0150), 5e-5))
teste("âncora curso: τ² residual .0764, I² 68,39%, R² 75,62%, QE(11) = 30,7331, QM(1) = 16,3571 (tol meia última casa)",
      perto(mp$estatistica[3], .0764, 5e-5) && perto(mp$estatistica[4], 68.39, .005) && perto(mp$estatistica[5], 75.62, .005) &&
        perto(mp$estatistica[2], 30.7331, 5e-5) && mp$gl[2] == 11 && perto(mp$estatistica[1], 16.3571, 5e-5))
teste("âncora curso: RR cai cerca de 2,9% por grau de latitude (1 - exp(-.0291) = .0287)", perto(-mt$mudanca_pct, 2.87, .01))
ma2 <- meta_regressao(eb, "alloc")
teste("âncora curso: alocação QM(2) = 1,7675, p = .4132, R² = 0%, τ² residual = .3615, I² = 88,77%",
      perto(ma2$pressupostos$estatistica[1], 1.7675, 5e-5) && perto(ma2$pressupostos$p[1], .4132, 5e-5) && perto(ma2$pressupostos$estatistica[5], 0, .005) &&
        perto(ma2$pressupostos$estatistica[3], .3615, 5e-5) && perto(ma2$pressupostos$estatistica[4], 88.77, .005))
ma3 <- meta_regressao(eb, c("ablat", "alloc"))
teste("âncora curso: latitude + alocação QM(3) = 11,0605, R² = 53,84%, latitude B = -.0273 (EP .0092, z = -2,9650)",
      perto(ma3$pressupostos$estatistica[1], 11.0605, 5e-5) && perto(ma3$pressupostos$estatistica[5], 53.84, .005) &&
        perto(ma3$tabela$B[2], -.0273, 5e-5) && perto(ma3$tabela$ep[2], .0092, 5e-5) && perto(ma3$tabela$estatistica[2], -2.9650, 5e-5))
teste("13 estudos: roda Peters (RR) como teste principal de assimetria e trim-and-fill só como sensibilidade",
      any(grepl("Peters", mb$pressupostos$checagem)) && any(grepl("trim-and-fill", mb$tabela$modelo)))

cat("== Meta-análise: equivalência e verdade conhecida\n")
set.seed(51)
k <- 40; n1 <- sample(15:120, k, TRUE); n2 <- sample(15:120, k, TRUE); delta <- rnorm(k, .4, .2)
smd <- data.frame(m1 = NA, d1 = NA, m2 = NA, d2 = NA, n1 = n1, n2 = n2)
for (i in 1:k) { a <- rnorm(n1[i], delta[i]); b <- rnorm(n2[i]); smd[i, c("m1", "d1", "m2", "d2")] <- c(mean(a), sd(a), mean(b), sd(b)) }
es1 <- calcular_efeitos_meta(smd, "SMD", media1 = "m1", dp1 = "d1", n1 = "n1", media2 = "m2", dp2 = "d2", n2 = "n2")
mm <- smd$n1 + smd$n2 - 2; J <- exp(lgamma(mm / 2) - log(sqrt(mm / 2)) - lgamma((mm - 1) / 2))
sp <- sqrt(((smd$n1 - 1) * smd$d1^2 + (smd$n2 - 1) * smd$d2^2) / mm); gman <- J * (smd$m1 - smd$m2) / sp
teste("equivalência: g de Hedges com correção exata e variância 1/n1 + 1/n2 + g²/(2(n1+n2)) (tol 1e-10)",
      perto(es1$yi, gman, 1e-10) && perto(es1$vi, 1 / smd$n1 + 1 / smd$n2 + gman^2 / (2 * (smd$n1 + smd$n2)), 1e-10))
ms <- meta_analise(es1)
teste("verdade conhecida: g combinado = .40 (± .08), intervalo de predição cobre .40, τ² < .12",
      abs(ms$efeito$valor[1] - .4) < .08 && ms$efeito$ic_inf[2] < .4 && ms$efeito$ic_sup[2] > .4 && ms$tabela$tau2[1] < .12)
eg_ref <- summary(lm(I(es1$yi / es1$sei) ~ I(1 / es1$sei)))$coefficients
eg <- ms$vies[ms$vies$teste == "Egger, regressão no EP", ]
teste("equivalência: Egger = intercepto da regressão de y/EP em 1/EP (Egger et al., 1997) (tol 1e-8)",
      perto(eg$vies, eg_ref[1, 1], 1e-8) && perto(eg$p, eg_ref[1, 4], 1e-8))
teste("SMD: leitura principal de assimetria pela regressão em 1/√n, com aviso sobre o Egger", any(grepl("1/√n", ms$pressupostos$checagem)) && any(grepl("SMD o Egger", ms$avisos)))
cor_d <- data.frame(r = c(.1, .3, .25, .4, .2), n = c(50, 120, 80, 60, 200))
ezc <- calcular_efeitos_meta(cor_d, "COR", r = "r", n = "n")
mz <- meta_analise(ezc)
teste("equivalência: z de Fisher = atanh(r), vi = 1/(n - 3) e síntese volta para r (tanh) sem zero à esquerda na frase",
      perto(ezc$yi, atanh(cor_d$r), 1e-12) && perto(ezc$vi, 1 / (cor_d$n - 3), 1e-12) &&
        perto(mz$efeito$valor[1], tanh(mz$tabela$estimativa[1]), 1e-12) && grepl("r = \\.[0-9]{2}", mz$apa[1]))
# OR: log OR manual; Harbord e Peters sob ausência e presença de efeito de estudos pequenos
set.seed(52)
k <- 60; ns <- round(exp(runif(k, log(30), log(600))))
sim_or <- function(ns, lor) { p0 <- .3; p1 <- plogis(qlogis(p0) + lor)
  data.frame(e1 = rbinom(length(ns), ns, p1), n1 = ns, e2 = rbinom(length(ns), ns, p0), n2 = ns) }
so <- sim_or(ns, .3)
eo <- calcular_efeitos_meta(so, "OR", eventos1 = "e1", n1 = "n1", eventos2 = "e2", n2 = "n2")
a <- so$e1; b <- so$n1 - so$e1; c <- so$e2; d <- so$n2 - so$e2
teste("equivalência: log OR = log(ad/bc) e vi = 1/a + 1/b + 1/c + 1/d (tol 1e-10)",
      perto(eo$yi, log(a * d / (b * c)), 1e-10) && perto(eo$vi, 1 / a + 1 / b + 1 / c + 1 / d, 1e-10))
mo <- meta_analise(eo)
teste("verdade conhecida: sem efeito de estudos pequenos, Peters e Harbord não rejeitam (p > .01)",
      all(mo$vies$p[mo$vies$teste %in% c("Peters, regressão em 1/N", "Harbord, escore")] > .01))
# seleção: estudos pequenos só "publicados" se significativos e na direção esperada
z <- eo$yi / eo$sei
sel <- eo[eo$n_total > 400 | z > 1.96, ]
attr(sel, "medida") <- "OR"; attr(sel, "contagens") <- attr(eo, "contagens")
msel <- meta_analise(sel)
teste("verdade conhecida: com seleção dos estudos pequenos significativos, Peters e Harbord detectam assimetria (p < .05)",
      nrow(sel) >= 10 && all(msel$vies$p[msel$vies$teste %in% c("Peters, regressão em 1/N", "Harbord, escore")] < .05))
teste("Hartung-Knapp como principal: estatística t com gl = k - 1", { h <- meta_analise(eb, hartung_knapp = TRUE); h$tabela$gl[1] == 12 && grepl("t\\(12\\)", h$apa[1]) })
teste("meta-regressão avisa menos de 10 estudos por covariável", any(grepl("menos de 10 estudos por covariável", ma3$avisos)))
teste("meta-regressão com permutação (semente 2026) devolve p de permutação", {
  pr <- meta_regressao(eb, "ablat", permutacao = TRUE, reamostras = 200); is.numeric(pr$tabela$p_permutacao) && pr$tabela$p_permutacao[2] < .05 })

cat("== Meta-análise: correções da revisão CoVe\n")
d5 <- data.frame(yi = c(.1, .5, .3, .8, .2), vi = c(.02, .03, .04, .02, .05))
r5 <- meta_analise(d5, vies = FALSE); f5 <- metafor::rma(d5$yi, d5$vi)
teste("intervalo de predição com t(k - 2) (Riley et al., 2011): k = 5 dá [-0,51; 1,29] = predict(predtype = 'Riley') (tol 1e-8)",
      perto(c(r5$efeito$ic_inf[2], r5$efeito$ic_sup[2]), unlist(predict(f5, predtype = "Riley")[c("pi.lb", "pi.ub")]), 1e-8) &&
        perto(c(r5$efeito$ic_inf[2], r5$efeito$ic_sup[2]), c(-.51, 1.29), .005))
teste("intervalo de predição do BCG com t(11): RR de .13 a 1,78 (tol .005) e frase APA com o mesmo intervalo",
      perto(c(mb$efeito$ic_inf[2], mb$efeito$ic_sup[2]), c(.134, 1.785), .0005) && grepl("\\[0\\.13, 1\\.78\\]", mb$apa[1]))
eor <- calcular_efeitos_meta(bcg, "OR", estudo = "estudo", eventos1 = "tpos", n1 = "n1", eventos2 = "cpos", n2 = "n2")
mor <- meta_analise(eor); pt_ <- mor$vies[mor$vies$teste == "Peters, regressão em 1/N", ]
Np <- bcg$n1 + bcg$n2; Sp <- bcg$tpos + bcg$cpos
pref <- summary(lm(eor$yi ~ I(1 / Np), weights = Sp * (Np - Sp) / Np))$coefficients
teste("Peters et al. (2006): regressão em 1/N com pesos S·F/N; BCG OR t(11) = -1,362, p = ,2003 (tol 1e-6)",
      perto(c(pt_$estatistica, pt_$p), pref[2, 3:4], 1e-6) && perto(pt_$estatistica, -1.3625, 5e-5) && pt_$gl == 11)
if (requireNamespace("meta", quietly = TRUE)) {
  mbias <- meta::metabias(meta::metabin(tpos, n1, cpos, n2, data = bcg, sm = "OR"), method.bias = "Peters", k.min = 5)
  teste("equivalência: Peters = meta::metabias(method.bias = 'Peters') (tol 1e-6)", perto(c(pt_$estatistica, pt_$p), c(mbias$statistic, mbias$pval), 1e-6))
} else cat("  aviso: pacote meta ausente, equivalência com metabias pulada\n")
mee <- meta_analise(ec, metodo = "EE")
teste("metodo = 'EE': frase diz efeito comum, sem 'efeitos aleatórios', sem intervalo de predição nem NA",
      grepl("efeito comum", mee$apa[1]) && !grepl("aleatórios", mee$apa[1]) && !any(grepl("NA", mee$apa)) && !grepl("predição", mee$apa[1]))
shk <- meta_subgrupos(eb, "alloc", hartung_knapp = TRUE)
teste("subgrupos com Hartung-Knapp avisam que o Q entre fica mais liberal", any(grepl("mais liberal", shk$avisos)))
mrf <- meta_regressao(transform(eb, aloc_n = as.integer(factor(alloc))), "aloc_n")
teste("meta-regressão: moderador numérico com poucos valores avisa e sugere fatores =", any(grepl('fatores = "aloc_n"', mrf$avisos, fixed = TRUE)))

cat("== Correções da verificação final (F3, H-01, G-11)\n")
# F3. Com nivel = .90 o rótulo do IC era "IC 95%" ao lado de um intervalo de 90%. Números
# conferidos fora da skill: metafor::rma(level = 90) para a meta-análise e o subgrupo;
# geepack::geeglm (EP robusto) e glm() com z de .95 para os pares de GEE e GLM.
d90 <- data.frame(yi = c(.2, .4, .3, .5, .1, .9, .7, .6), vi = c(.02, .03, .02, .04, .03, .05, .04, .03), sg = rep(c("A", "B"), each = 4))
m90 <- meta_analise(d90, vies = FALSE, nivel = .90); r90 <- metafor::rma(yi, vi, data = d90, level = 90)
teste("F3 meta_analise(nivel = .90): 'IC 90% [0.29, 0.58]' = metafor::rma(level = 90) [.288; .579], sem 'IC 95%'; em inglês '90% CI'",
      grepl("efeito = 0.43, IC 90% [0.29, 0.58]", m90$apa[1], fixed = TRUE) && !grepl("95%", m90$apa[1]) &&
        perto(c(r90$ci.lb, r90$ci.ub), c(.29, .58), .005) &&
        grepl("90% CI [0.29, 0.58]", meta_analise(d90, vies = FALSE, nivel = .90, idioma = "en")$apa[1], fixed = TRUE))
s90 <- meta_subgrupos(d90, "sg", nivel = .90); rA <- metafor::rma(yi, vi, data = d90[d90$sg == "A", ], level = 90)
teste("F3 meta_subgrupos(nivel = .90): 'A (k = 4): efeito = 0.32, IC 90% [0.19, 0.45]' = metafor::rma(level = 90) do subgrupo",
      grepl("A (k = 4): efeito = 0.32, IC 90% [0.19, 0.45]", s90$apa[1], fixed = TRUE) && !grepl("95%", s90$apa[1]) &&
        perto(c(rA$ci.lb, rA$ci.ub), c(.19, .45), .005))
set.seed(2); n90 <- 60
g90 <- data.frame(id = rep(1:n90, each = 3), tempo = rep(1:3, n90), grupo = rep(sample(c("c", "t"), n90, TRUE), each = 3))
g90$y <- rnorm(3 * n90) + (g90$grupo == "t") * .4
e90 <- modelo_gee(g90, "y", "grupo", id = "id", tempo = "tempo", nivel = .90)
cg90 <- summary(geepack::geeglm(y ~ grupo, data = g90, id = id, waves = tempo, corstr = "unstructured"))$coefficients
ic_g <- -cg90["grupot", "Estimate"] + c(-1, 1) * qnorm(.95) * cg90["grupot", "Std.err"]
teste("F3 modelo_gee(nivel = .90): 'c - t: diferença = -0.28, IC 90% [-0.54, -0.02], p = .076' (geepack: -0,281 ± 1,645 × 0,159)",
      grepl("c - t: diferença = -0.28, IC 90% [-0.54, -0.02], p = .076", e90$apa[2], fixed = TRUE) && perto(ic_g, c(-.54, -.02), .005) &&
        grepl("90% CI [-0.54, -0.02]", modelo_gee(g90, "y", "grupo", id = "id", tempo = "tempo", nivel = .90, idioma = "en")$apa[2], fixed = TRUE))
set.seed(5); p90 <- data.frame(grupo = rep(c("c", "t"), each = 40)); p90$z <- rpois(80, ifelse(p90$grupo == "t", 4, 3))
gp90 <- modelo_glm(p90, "z", "grupo", familia = "poisson", nivel = .90)
cp90 <- summary(glm(z ~ grupo, poisson, p90))$coefficients
ic_p <- exp(-cp90["grupot", 1] + c(-1, 1) * qnorm(.95) * cp90["grupot", 2])
teste("F3 modelo_glm(Poisson, nivel = .90): par c / t com 'IC 90%' e limites = exp(-b ± 1,645 EP) do glm()",
      grepl(paste0("c / t: razão = ", sprintf("%.2f", exp(-cp90["grupot", 1])), ", IC 90% [", sprintf("%.2f", ic_p[1]), ", ", sprintf("%.2f", ic_p[2]), "]"),
            gp90$apa[2], fixed = TRUE) && !grepl("95%", gp90$apa[2]))

# H-01. Médias marginais de GLM e GEE com nível de 1 caso: o print mostrava o escore da
# pessoa. Agora $posthoc tem a coluna n (casos; sujeitos no GEE) e o print apaga média e
# IC de nível com n < 5. Número independente: a média de um nível com um caso é o próprio
# escore (97,531); no GEE intercambiável só com grupo e tempo, a média marginal de um
# grupo de um sujeito é a média das 3 medidas dele (93,456), porque R⁻¹1 é proporcional a 1.
set.seed(7)
p1 <- data.frame(g = c(rep("a", 30), rep("b", 30), "c"), x = round(rnorm(61, 50, 10), 3)); p1$x[61] <- 97.531
p1$z <- rpois(61, 3); p1$z[61] <- 17
h1 <- modelo_glm(p1, "x", "g"); o1 <- capture.output(print(h1))
m1c <- h1$posthoc[h1$posthoc$tipo == "média marginal" & h1$posthoc$nivel_ou_par == "c", ]
teste("H-01 GLM: nível com 1 caso fica no objeto (média 97,531, n = 1; n = 30 nos outros) e sai do print e da frase, com aviso",
      perto(m1c$estimativa, 97.531, 1e-8) && m1c$n == 1 && all(h1$posthoc$n[h1$posthoc$tipo == "média marginal" & h1$posthoc$nivel_ou_par != "c"] == 30) &&
        !any(grepl("97.53", o1, fixed = TRUE)) && !any(grepl(format(signif(m1c$ic_inf, 4)), o1, fixed = TRUE)) &&
        !any(grepl(format(signif(m1c$ic_sup, 4)), o1, fixed = TRUE)) && !any(grepl("97.53", h1$apa, fixed = TRUE)) &&
        any(grepl("g = c, n = 1", h1$avisos, fixed = TRUE)) && any(grepl("n < 5 sem estatísticas no print", o1, fixed = TRUE)) &&
        any(grepl(format(signif(mean(p1$x[p1$g == "a"]), 4)), o1, fixed = TRUE)) && inherits(h1, "pr_resultado"))
h1p <- modelo_glm(p1, "z", "g", familia = "poisson"); o1p <- capture.output(print(h1p))
lin_c <- o1p[grepl("marginal", o1p) & grepl(" c ", o1p, fixed = TRUE)]
teste("H-01 GLM Poisson: média marginal do nível de 1 caso (contagem 17) no objeto e linha do print sem valores",
      perto(h1p$posthoc$estimativa[h1p$posthoc$tipo == "média marginal" & h1p$posthoc$nivel_ou_par == "c"], 17, 1e-6) &&
        length(lin_c) == 1 && grepl("NA", lin_c) && !grepl("17", lin_c, fixed = TRUE))
set.seed(3); K1 <- 31
gs <- data.frame(id = rep(1:K1, each = 3), tempo = rep(1:3, K1), grupo = rep(c(rep(c("a", "b"), 15), "c"), each = 3))
gs$y <- rnorm(3 * K1, 50, 10); gs$y[gs$grupo == "c"] <- c(91.234, 93.456, 95.678)
# o emmeans avisa variância negativa no grupo de um só sujeito (EP sanduíche degenerado com 1 cluster)
hg <- suppressWarnings(modelo_gee(gs, "y", c("grupo", "tempo"), id = "id", tempo = "tempo", correlacao = "exchangeable")); og <- capture.output(print(hg))
mgc <- hg$posthoc[hg$posthoc$tipo == "média marginal" & hg$posthoc$nivel_ou_par == "c", ]
teste("H-01 GEE: grupo de 1 sujeito, média marginal 93,456 (média das 3 medidas) e n = 1 sujeito no objeto, fora do print e da frase",
      perto(mgc$estimativa, mean(c(91.234, 93.456, 95.678)), 1e-6) && mgc$n == 1 &&
        all(hg$posthoc$n[hg$posthoc$termo == "tempo" & hg$posthoc$tipo == "média marginal"] == K1) &&
        !any(grepl("93.46", og, fixed = TRUE)) && !any(grepl("93.46", hg$apa, fixed = TRUE)) &&
        any(grepl("menos de 5 sujeitos (grupo = c, n = 1)", hg$avisos, fixed = TRUE)))
# Cluster único num nível de fator tem alavanca 1: a correção de Mancl e DeRouen não se
# define e a função cai para o EP robusto com aviso (no Mac, o solve() dessa matriz
# quebrava a função; no container, passava com número sem sentido).
teste("GEE com sujeito que determina um parâmetro: Mancl e DeRouen indefinida vira aviso, sem erro e sem colunas de MD",
      any(grepl("Mancl e DeRouen não se define", hg$avisos, fixed = TRUE)) && !("p_mancl_derouen" %in% names(hg$tabela)))

# G-11. "imputed 1 studies" em inglês. Banco com semente 4 (k = 11) em que o trim-and-fill
# imputa 1 estudo: metafor::trimfill(rma()) dá k0 = 1 e efeito ajustado 0,74.
tf_d <- data.frame(yi = c(.624, .68, .63, .679, 1.002, 1.091, .787, 1.101, 1.13, .555, .623),
                   vi = c(.051, .0106, .0306, .0294, .067, .0282, .0607, .0734, .0764, .0151, .0628))
tf_ref <- metafor::trimfill(metafor::rma(tf_d$yi, tf_d$vi))
teste("G-11 trim-and-fill com 1 estudo imputado: 'imputed 1 study' e 'imputou 1 estudo', efeito ajustado 0.74 = metafor::trimfill",
      tf_ref$k0 == 1 && perto(as.numeric(tf_ref$b), .74, .005) &&
        grepl("imputed 1 study, adjusted effect = 0.74", meta_analise(tf_d, idioma = "en")$apa[3], fixed = TRUE) &&
        grepl("imputou 1 estudo, efeito ajustado = 0.74", meta_analise(tf_d)$apa[3], fixed = TRUE))

cat("== v1.4: meta-análise de prevalência e de acurácia diagnóstica\n")
mp_d <- data.frame(est = paste0("S", 1:8), ev = c(12, 30, 8, 45, 22, 5, 60, 17), n = c(100, 210, 90, 300, 150, 80, 400, 120))
teste("prevalência: IC de cada estudo = Clopper-Pearson de binom.test (tol. 1e-12)", {
  r <- meta_proporcao(mp_d, "ev", "n", "est"); e <- r$tabelas_extra$Estudos
  all(abs(e$ic_inf - mapply(function(a, b) binom.test(a, b)$conf.int[1], mp_d$ev, mp_d$n)) < 1e-12) })
if (requireNamespace("meta", quietly = TRUE)) {
  teste("prevalência igual ao meta::metaprop (inverso da variância, REML) em logito, Freeman-Tukey e bruta (tol. 1e-6)", {
    all(vapply(c(PLO = "PLOGIT", PFT = "PFT", PR = "PRAW"), function(sm) {
      tr <- names(which(c(PLO = "PLOGIT", PFT = "PFT", PR = "PRAW") == sm))
      r <- meta_proporcao(mp_d, "ev", "n", "est", transformacao = tr)
      m <- meta::metaprop(ev, n, data = mp_d, sm = sm, method = "Inverse", method.tau = "REML")
      bt <- function(x) meta:::backtransf(x, m$sm, n = 1 / mean(1 / mp_d$n))
      max(abs(unlist(r$tabela[1, c("prevalencia", "ic_inf", "ic_sup")]) - c(bt(m$TE.random), bt(m$lower.random), bt(m$upper.random)))) < 1e-6 }, TRUE)) })
} else cat("PULADO  comparação com meta::metaprop (não instalado)\n")
# AuditC (mada): 14 estudos do AUDIT-C; valores do mada::reitsma 0.4.x: sens .8908781 [.8080775; .9405829], FPR .2203392 [.1672318; .2845490].
ac <- data.frame(TP = c(47, 126, 19, 36, 130, 84, 68, 752, 59, 142, 137, 57, 34, 152), FN = c(9, 51, 10, 3, 19, 2, 0, 0, 5, 50, 24, 3, 1, 51),
                 FP = c(101, 272, 12, 78, 211, 68, 112, 3226, 55, 571, 107, 103, 21, 88), TN = c(738, 1543, 192, 276, 959, 89, 423, 2977, 136, 2788, 358, 437, 56, 264))
ra <- meta_acuracia(ac, "TP", "FP", "FN", "TN")
teste("acurácia (AuditC): sensibilidade .8909 [.8081; .9406] e especificidade .7797 [.7155; .8328], como o mada::reitsma (tol. 1e-5)", {
  abs(ra$tabela$valor[1] - .8908781) < 1e-5 && abs(ra$tabela$ic_inf[1] - .8080775) < 1e-5 && abs(ra$tabela$ic_sup[1] - .9405829) < 1e-5 &&
    abs(ra$tabela$valor[2] - (1 - .2203392)) < 1e-5 && abs(ra$tabela$ic_inf[2] - (1 - .2845490)) < 1e-5 })
teste("acurácia: correção de 0,5 avisada (estudos 7 e 8 têm FN = 0) e RV+ = sens / (1 - espec) combinadas", {
  any(grepl("Célula zero", ra$avisos)) && abs(ra$tabela$valor[3] - ra$tabela$valor[1] / (1 - ra$tabela$valor[2])) < 1e-10 })
if (requireNamespace("mada", quietly = TRUE)) {
  teste("acurácia igual ao mada::reitsma ao vivo (τ² e correlação, com o sinal trocado pela FPR; tol. 1e-4)", {
    mr <- mada::reitsma(ac); h <- ra$tabelas_extra$Heterogeneidade
    abs(h$tau2_logito_sens - mr$Psi[1, 1]) < 1e-4 && abs(h$tau2_logito_espec - mr$Psi[2, 2]) < 1e-4 &&
      abs(h$correlacao + mr$Psi[1, 2] / sqrt(mr$Psi[1, 1] * mr$Psi[2, 2])) < 1e-4 })
} else cat("PULADO  comparação ao vivo com mada (não instalado; os valores publicados acima cobrem)\n")

cat("\n== Meta-análise em rede (v1.5, netmeta)\n")
if (requireNamespace("netmeta", quietly = TRUE)) {
  data(Franchini2012, package = "netmeta")
  fr <- Franchini2012
  longo <- do.call(rbind, lapply(1:3, function(j) data.frame(estudo = fr$Study, trat = fr[[paste0("Treatment", j)]], n = fr[[paste0("n", j)]],
                                                              m = fr[[paste0("y", j)]], s = fr[[paste0("sd", j)]], stringsAsFactors = FALSE)))
  longo <- longo[!is.na(longo$trat), ]
  nr <- suppressWarnings(meta_rede(longo, "estudo", "trat", "n", media = "m", dp = "s", referencia = "Placebo", modelo = "comum"))
  pw_f <- if ("pairwise" %in% getNamespaceExports("netmeta")) getExportedValue("netmeta", "pairwise") else getExportedValue("meta", "pairwise")
  pw <- suppressWarnings(pw_f(list(Treatment1, Treatment2, Treatment3), n = list(n1, n2, n3), mean = list(y1, y2, y3),
                                           sd = list(sd1, sd2, sd3), data = fr, studlab = Study, sm = "MD"))
  nd <- suppressWarnings(netmeta::netmeta(pw, reference.group = "Placebo"))
  teste("Franchini 2012 (banco da aula): efeitos do modelo comum iguais ao netmeta direto no formato largo (tol. 1e-10)", {
    t <- nr$tabela[nr$tabela$modelo == "comum", ]; max(abs(t$efeito - nd$TE.common[t$tratamento, "Placebo"])) < 1e-10 })
  teste("decomposição do Q: total = dentro + entre desenhos", { q <- nr$tabelas_extra[["Decomposição do Q"]]$Q; abs(q[1] - q[2] - q[3]) < 1e-8 })
  teste("P-scores com média 0,5", abs(mean(nr$tabelas_extra[["Ranking (P-score)"]]$p_score) - .5) < 1e-10)
  dois <- longo[longo$trat %in% c("Placebo", "Pramipexole") & longo$estudo %in% longo$estudo[longo$trat == "Pramipexole"], ]
  dois <- dois[dois$estudo %in% names(which(table(dois$estudo) == 2)), ]
  n2 <- suppressWarnings(meta_rede(dois, "estudo", "trat", "n", media = "m", dp = "s", referencia = "Placebo"))
  teste("rede de dois tratamentos = meta-análise par a par do metafor (comum: FE; aleatório: DL; tol. 1e-8)", {
    es <- metafor::escalc("MD", m1i = dois$m[dois$trat == "Pramipexole"], sd1i = dois$s[dois$trat == "Pramipexole"], n1i = dois$n[dois$trat == "Pramipexole"],
                          m2i = dois$m[dois$trat == "Placebo"], sd2i = dois$s[dois$trat == "Placebo"], n2i = dois$n[dois$trat == "Placebo"])
    fe <- metafor::rma(yi, vi, data = es, method = "FE"); re <- metafor::rma(yi, vi, data = es, method = "DL")
    t <- n2$tabela; abs(t$efeito[t$modelo == "comum"] - coef(fe)) < 1e-8 && abs(t$efeito[t$modelo == "aleatório"] - coef(re)) < 1e-8 })
} else cat("PULADO  meta-análise em rede (netmeta não instalado)\n")

cat("\n== Extração de dados para meta-análise (v1.5)\n")
teste("mediana e IIQ: média de Luo = (0,7 + 0,39/n)(q1 + q3)/2 + (0,3 - 0,39/n) mediana (exemplo da aula: n 15; 9,6; 14,7; 30,7)", {
  x <- extrair_media_dp("mediana_iqr", n = 15, mediana = 14.7, q1 = 9.6, q3 = 30.7)
  abs(x$media - ((0.7 + 0.39 / 15) * (9.6 + 30.7) / 2 + (0.3 - 0.39 / 15) * 14.7)) < 1e-10 })
teste("mediana e IIQ: média de Wan = (q1 + mediana + q3)/3", abs(extrair_media_dp("mediana_iqr", n = 15, mediana = 14.7, q1 = 9.6, q3 = 30.7, metodo_media = "Wan")$media - 55 / 3) < 1e-10)
teste("mediana e IIQ, n > 201: DP = IIQ / [2 Φ⁻¹((0,75n - 0,125)/(n + 0,25))]", {
  n <- 300; abs(extrair_media_dp("mediana_iqr", n = n, mediana = 10, q1 = 7, q3 = 14)$dp - 7 / (2 * stats::qnorm((0.75 * n - 0.125) / (n + 0.25)))) < 1e-10 })
teste("mediana e amplitude, n = 15: DP = amplitude / amplitude esperada da normal (integração numérica; tol. 1e-3)", {
  er <- stats::integrate(function(x) 1 - stats::pnorm(x)^15 - (1 - stats::pnorm(x))^15, -Inf, Inf)$value
  abs(extrair_media_dp("mediana_amplitude", n = 15, mediana = 14.7, minimo = 3, maximo = 40)$dp - 37 / er) < 1e-3 })
set.seed(41); a1 <- stats::rnorm(23, 10, 3); a2 <- stats::rnorm(31, 13, 4)
teste("combinar dois grupos = média e DP do banco unido (tol. 1e-10)", {
  x <- extrair_media_dp("combinar", n = 23, media = mean(a1), dp1 = stats::sd(a1), n2 = 31, media2 = mean(a2), dp2 = stats::sd(a2))
  abs(x$media - mean(c(a1, a2))) < 1e-10 && abs(x$dp - stats::sd(c(a1, a2))) < 1e-10 })
pre <- stats::rnorm(40, 20, 5); pos <- .6 * pre + stats::rnorm(40, 10, 3)
teste("DP da mudança com a r observada = DP de pós - pré (tol. 1e-10)", {
  x <- extrair_media_dp("mudanca", dp_pre = stats::sd(pre), dp_pos = stats::sd(pos), r = stats::cor(pre, pos)); abs(x$dp - stats::sd(pos - pre)) < 1e-10 })
teste("r implícita recupera a correlação pré-pós (tol. 1e-10)", {
  abs(extrair_media_dp("correlacao", dp_pre = stats::sd(pre), dp_pos = stats::sd(pos), dp_mudanca = stats::sd(pos - pre))$r - stats::cor(pre, pos)) < 1e-10 })
teste("DP pelo IC: t com n < 60 e z com n >= 60 (Handbook)", {
  x1 <- extrair_media_dp("ic", n = 20, ic_inf = 8, ic_sup = 12); x2 <- extrair_media_dp("ic", n = 100, ic_inf = 8, ic_sup = 12)
  abs(x1$dp - sqrt(20) * 4 / (2 * stats::qt(.975, 19))) < 1e-10 && abs(x2$dp - 10 * 4 / (2 * stats::qnorm(.975))) < 1e-10 })
teste("argumento faltando gera erro claro", inherits(try(extrair_media_dp("ep", n = 10), silent = TRUE), "try-error"))

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok), as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
print(checar_pacotes(c("lme4", "lmerTest", "geepack", "metafor", "emmeans", "car", "MASS", "sandwich", "performance")))
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
