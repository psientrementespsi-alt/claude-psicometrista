# evals_v16.R, testes das funções da v1.6: modelo_ordinal, modelo_multinomial,
# efeito_contextual, concordancia_juizes (AC1/AC2 de Gwet e alfa de Krippendorff),
# confiabilidade(ic = TRUE), spearman_brown, corrigir_correlacao, generalizabilidade,
# analise_distratores, dmacs, tri_dif(metodo = "mh"), regressao_logistica(validar = TRUE),
# risco e NNT no qui_quadrado 2 × 2 e PA-MRFA em afe_numero_fatores.
#
# Rodar da pasta da skill:  Rscript scripts/evals_v16.R
# Âncoras: exemplos dos livros da pasta (Martinková & Hladká, 2023; Petersen, 2024),
# bancos públicos de pacotes (MASS::housing, lavaan::HolzingerSwineford1939, difR::verbal,
# ShinyItemAnalysis::HCIlong/HCItest) e os exemplos da aula de limitações do kappa.
# Implementações de referência quando instaladas: VGAM, ordinal, irrCAC, gtheory e dmacs
# (arquivados no CRAN; instalar do arquivo só para conferir), pinsearch, difR, rms,
# DescTools e EFA.MRFA. Sem elas, o teste correspondente é pulado.

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

cat("== modelo_ordinal (MASS::housing, Venables & Ripley)\n")
h <- MASS::housing; hh <- h[rep(seq_len(nrow(h)), h$Freq), c("Sat", "Infl", "Type", "Cont")]
ro <- quieto(modelo_ordinal(hh, "Sat", c("Infl", "Type", "Cont")))
ref <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = h)
teste("coeficientes e limiares iguais ao polr ponderado (1e-6)", max(abs(ro$coeficientes$B - coef(ref))) < 1e-6 && max(abs(ro$limiares$limiar - ref$zeta)) < 1e-6)
teste("valores publicados: InflHigh 1.2888, TypeTerrace -1.0910, limiar Low|Medium -0.4961 (4 casas)",
      abs(ro$coeficientes$B[2] - 1.2888) < 5e-4 && abs(ro$coeficientes$B[5] + 1.0910) < 5e-4 && abs(ro$limiares$limiar[1] + 0.4961) < 5e-4)
teste("OR = exp(B) e R² de Nagelkerke = Cox-Snell / máximo", abs(ro$efeito$valor[1] - exp(ro$coeficientes$B[1])) < 1e-12 &&
        abs(ro$pseudo_r2$valor[2] - ro$pseudo_r2$valor[1] / (1 - exp(2 * as.numeric(logLik(MASS::polr(Sat ~ 1, data = hh))) / nrow(hh)))) < 1e-8)
if (tem("ordinal")) {
  teste("linhas paralelas por variável iguais ao ordinal::nominal_test (1e-4)", {
    nt <- as.data.frame(ordinal::nominal_test(ordinal::clm(Sat ~ Infl + Type + Cont, weights = Freq, data = h)))[-1, ]
    max(abs(ro$paralelas$qui2_rv - nt$LRT)) < 1e-4 })
  teste("teste global de linhas paralelas não rejeita (p > .05)", ro$pressupostos$p[2] > .05)
} else pular("linhas paralelas", "ordinal não instalado")
teste("probit: sinais iguais aos do logit e razão perto de 1,7", {
  rp <- quieto(modelo_ordinal(hh, "Sat", c("Infl", "Type", "Cont"), ligacao = "probit"))
  all(sign(rp$coeficientes$B) == sign(ro$coeficientes$B)) && abs(median(ro$coeficientes$B / rp$coeficientes$B) - 1.7) < .15 })

cat("== modelo_multinomial\n")
set.seed(1); n <- 600; x <- rnorm(n); g <- factor(sample(c("a", "b"), n, TRUE))
P <- cbind(1, exp(-.2 + .8 * x + .5 * (g == "b")), exp(.3 - .6 * x + .2 * (g == "b"))); P <- P / rowSums(P)
dm <- data.frame(y = apply(P, 1, function(p) sample(c("c1", "c2", "c3"), 1, prob = p)), x = x, g = g)
rm_ <- quieto(modelo_multinomial(dm, "y", c("x", "g")))
if (tem("VGAM")) {
  teste("coeficientes iguais ao VGAM::vglm multinomial (1e-4)", {
    v <- VGAM::vglm(factor(y) ~ x + g, VGAM::multinomial(refLevel = 1), data = dm); cv <- coef(v)
    b <- rm_$coeficientes
    max(abs(c(b$B[b$categoria == "c2"], b$B[b$categoria == "c3"]) - c(cv[c(1, 3, 5)], cv[c(2, 4, 6)]))) < 1e-4 })
} else pular("multinomial vs VGAM", "VGAM não instalado")
teste("RV por termo = 2 × diferença de log-verossimilhança contra o modelo sem o termo", {
  m1 <- nnet::multinom(y ~ x + g, data = dm, trace = FALSE); m0 <- nnet::multinom(y ~ g, data = dm, trace = FALSE)
  abs(rm_$tabela$qui2_rv[1] - 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))) < 1e-4 })
teste("referencia_vd troca a categoria-base e inverte o sinal do contraste c2 vs c1", {
  r2 <- quieto(modelo_multinomial(dm, "y", c("x", "g"), referencia_vd = "c2"))
  b1 <- rm_$coeficientes; b2 <- r2$coeficientes
  abs(b1$B[b1$categoria == "c2" & b1$termo == "x"] + b2$B[b2$categoria == "c1" & b2$termo == "x"]) < 1e-3 })

cat("== efeito_contextual\n")
set.seed(4); J <- 60; nj <- 15; gg <- rep(1:J, each = nj); mu <- rnorm(J)
xx <- mu[gg] + rnorm(J * nj); yy <- .5 * (xx - mu[gg]) + 1.2 * mu[gg] + rnorm(J, 0, .5)[gg] + rnorm(J * nj)
dc <- centralizar_multinivel(data.frame(g = gg, x = xx, y = yy), "x", "g", tipo = "grupo")
ml <- quieto(multinivel_linear(dc, "y", "g", nivel1 = "x_cwc", nivel2 = "x_media_cgm"))
ec <- efeito_contextual(ml, "x")
teste("contextual (b_media - b_cwc) = coeficiente da média do cluster com X centrado na média geral (1e-6)", {
  f2 <- lme4::lmer(y ~ x + x_media_cgm + (1 | g), data = transform(dc, x = x - mean(x)))
  abs(ec$efeito$valor - lme4::fixef(f2)[["x_media_cgm"]]) < 1e-6 })

cat("== concordancia_juizes: AC1/AC2 de Gwet e alfa de Krippendorff\n")
mk <- function(ss, sn, ns, nn) { j1 <- c(rep("sim", ss + sn), rep("não", ns + nn)); j2 <- c(rep("sim", ss), rep("não", sn), rep("sim", ns), rep("não", nn))
  data.frame(j1 = factor(j1, levels = c("sim", "não")), j2 = factor(j2, levels = c("sim", "não"))) }
k_a <- quieto(concordancia_juizes(mk(10, 0, 0, 0), c("j1", "j2")))
k_b <- quieto(concordancia_juizes(mk(9, 1, 0, 0), c("j1", "j2")))
k_c <- quieto(concordancia_juizes(mk(8, 1, 0, 1), c("j1", "j2")))
teste("aula (10 sim-sim): κ incalculável, AC1 = 1 e aviso do paradoxo", is.na(k_a$tabela$kappa) && k_a$tabela$gwet == 1 && any(grepl("Paradoxo", k_a$avisos)))
teste("aula (9 sim-sim, 1 sim-não): κ = 0 com 90% de concordância; AC1 = .8895 e aviso", abs(k_b$tabela$kappa) < 1e-12 && abs(k_b$tabela$gwet - .889503) < 1e-5 && any(grepl("Paradoxo", k_b$avisos)))
teste("aula (8, 1, 0, 1): κ = .6154 = (.9 - .74)/.26", abs(k_c$tabela$kappa - (.9 - .74) / .26) < 1e-4)
set.seed(3); nr <- 40; tr <- sample(1:5, nr, TRUE)
mr <- sapply(1:4, function(j) pmin(5, pmax(1, tr + sample(-1:1, nr, TRUE, prob = c(.2, .6, .2))))); dr <- as.data.frame(mr)
c_nom <- quieto(concordancia_juizes(dr, names(dr), "nominal")); c_ord <- quieto(concordancia_juizes(dr, names(dr), "ordinal")); c_int <- quieto(concordancia_juizes(dr, names(dr), "intervalar"))
if (tem("irrCAC")) {
  teste("AC1 (4 juízes) e AC2 quadrático iguais ao irrCAC::gwet.ac1.raw (1e-5)",
        abs(c_nom$tabela$gwet - irrCAC::gwet.ac1.raw(dr)$est$coeff.val) < 1e-5 && abs(c_ord$tabela$gwet - irrCAC::gwet.ac1.raw(dr, weights = "quadratic")$est$coeff.val) < 1e-5)
  teste("Krippendorff nominal e intervalar iguais ao irrCAC::krippen.alpha.raw (1e-5)",
        abs(c_nom$tabela$krippendorff - irrCAC::krippen.alpha.raw(dr)$est$coeff.val) < 1e-5 &&
          abs(c_int$efeito$valor[3] - irrCAC::krippen.alpha.raw(dr, weights = "quadratic")$est$coeff.val) < 1e-5)
} else pular("Gwet e Krippendorff vs irrCAC", "irrCAC não instalado")
teste("Krippendorff ordinal = intervalar sobre os postos médios das categorias (identidade da métrica ordinal)", {
  nc <- table(factor(as.vector(mr), levels = 1:5)); post <- cumsum(nc) - (nc - 1) / 2
  mr2 <- matrix(post[mr], nrow(mr)); cats <- sort(unique(as.vector(mr2)))
  abs(c_ord$tabela$krippendorff - .ee_krippendorff(mr2, cats, "intervalar")) < 1e-10 })
teste("2 juízes ordinais: AC2 igual ao irrCAC e κ ponderado igual ao psych", {
  c2 <- quieto(concordancia_juizes(dr[, 1:2], c("V1", "V2"), "ordinal"))
  (!tem("irrCAC") || abs(c2$tabela$gwet - irrCAC::gwet.ac1.raw(dr[, 1:2], weights = "quadratic")$est$coeff.val) < 1e-5) &&
    abs(c2$tabela$kappa - psych::cohen.kappa(as.matrix(dr[, 1:2]))$weighted.kappa) < 1e-8 })

cat("== confiabilidade(ic = TRUE)\n")
set.seed(11); nn <- 400; ff <- rnorm(nn); lam <- c(.8, .7, .6, .7, .5, .6)
Xr <- as.data.frame(sapply(lam, function(l) l * ff + rnorm(nn, sd = sqrt(1 - l^2)))); names(Xr) <- paste0("x", 1:6)
fr <- lavaan::cfa("F =~ x1+x2+x3+x4+x5+x6", data = Xr, std.lv = TRUE)
cr <- quieto(confiabilidade(fr, Xr, ic = TRUE, reamostras = 200))
teste("IC do ômega por bootstrap perto do IC delta do parâmetro definido no lavaan (tol. .015)", {
  m2 <- paste("F =~ l1*x1+l2*x2+l3*x3+l4*x4+l5*x5+l6*x6", paste0("x", 1:6, "~~e", 1:6, "*x", 1:6, collapse = "\n"),
              "om := (l1+l2+l3+l4+l5+l6)^2/((l1+l2+l3+l4+l5+l6)^2 + e1+e2+e3+e4+e5+e6)", sep = "\n")
  pe <- lavaan::parameterEstimates(lavaan::cfa(m2, data = Xr, std.lv = TRUE)); o <- pe[pe$lhs == "om", ]
  abs(cr$omega_ic_inf - o$ci.lower) < .015 && abs(cr$omega_ic_sup - o$ci.upper) < .015 })
teste("IC do alfa por bootstrap perto do IC de Feldt (psych::alpha.ci; tol. .015)", {
  fe <- psych::alpha.ci(cr$alfa, nn, n.var = 6); abs(cr$alfa_ic_inf - fe$lower) < .015 && abs(cr$alfa_ic_sup - fe$upper) < .015 })
teste("sem ic = TRUE a saída da v1.5 não muda (colunas fator, omega, alfa, itens)", identical(names(confiabilidade(fr, Xr)), c("fator", "omega", "alfa", "itens")))

cat("== spearman_brown e corrigir_correlacao (Martinková & Hladká, 2023)\n")
teste("4.2.5.1: 20 itens com .75 levados a 30 = .8182; alvo .85 pede 38 itens",
      abs(spearman_brown(.75, 20, k_novo = 30)$tabela$confiabilidade_projetada - .8182) < 5e-5 && spearman_brown(.75, 20, alvo = .85)$tabela$itens_necessarios == 38)
set.seed(987); x1 <- rnorm(1500); e1 <- rnorm(1500); b1 <- .729 / sqrt(1 - .729^2); y1 <- b1 * x1 + e1
xv <- x1 * 10 + 50; yv <- y1 * 2.3 + 20; ad <- data.frame(xv, yv)[xv > 58, ]
rr <- corrigir_correlacao(cor(ad$xv, ad$yv), dp_irrestrito = sd(xv), dp_restrito = sd(ad$xv))
teste("2.4.3: r restrito .3633 corrigido para .6548 (= psych::rangeCorrection)", abs(rr$tabela$r[1] - .3633) < 5e-5 && abs(rr$efeito$valor - .6548) < 5e-5 &&
        abs(rr$efeito$valor - psych::rangeCorrection(cor(ad$xv, ad$yv), sd(xv), sd(ad$xv))) < 1e-12)
teste("atenuação: .30 com rxx = .70 e ryy = .80 dá .30/sqrt(.56) (= psych::correct.cor)",
      abs(corrigir_correlacao(.30, .70, .80)$efeito$valor - psych::correct.cor(matrix(c(1, .3, .3, 1), 2), c(.7, .8))[1, 2]) < 1e-12)
teste("valor corrigido acima de 1 gera aviso", any(grepl("acima de 1", corrigir_correlacao(.6, .5, .6)$avisos)))

cat("== generalizabilidade (Petersen, 2024, cap. 5; Martinková & Hladká, 2023, 4.5)\n")
set.seed(52242)
pio <- data.frame(Person = as.factor(rep(1:6, each = 8)), Item = as.factor(rep(1:4, times = 12)),
                  Score = c(9,9,7,4,9,8,5,5,9,8,4,6,6,5,3,3,8,8,6,2,8,7,3,2,9,8,6,3,9,6,6,2,10,9,8,7,8,8,9,7,6,4,5,1,3,2,3,2),
                  Occasion = as.factor(rep(1:2, each = 4, times = 6)))
pio$Score[30] <- NA
gt <- quieto(generalizabilidade(pio, "Person", "Score", c("Item", "Occasion")))
teste("componentes publicados: pessoa 2.4780, item 3.1238, pessoa×item .7873, resíduo .8374 (4 casas)", {
  v <- stats::setNames(gt$tabela$variancia, gt$tabela$componente)
  abs(v[["Person"]] - 2.47804) < 5e-4 && abs(v[["Item"]] - 3.12380) < 5e-4 && abs(v[["Person:Item"]] - .787305) < 5e-4 && abs(v[[7]] - .837426) < 5e-4 })
teste("Eρ² = .8731 e Φ = .6374, como o gtheory::dstudy no livro", abs(gt$efeito$valor[1] - .8731) < 5e-5 && abs(gt$efeito$valor[2] - .6374) < 5e-5)
if (tem("gtheory")) {
  teste("igual ao gtheory::dstudy (1e-6)", {
    gs <- gtheory::gstudy(data = pio, Score ~ (1 | Person) + (1 | Item) + (1 | Occasion) + (1 | Person:Item) + (1 | Person:Occasion) + (1 | Occasion:Item))
    ds <- gtheory::dstudy(gs, colname.objects = "Person", data = pio, colname.scores = "Score")
    abs(gt$efeito$valor[1] - ds$generalizability) < 1e-6 && abs(gt$efeito$valor[2] - ds$dependability) < 1e-6 })
} else pular("generalizabilidade vs gtheory", "gtheory não instalado (arquivado no CRAN)")
if (tem("ShinyItemAnalysis")) {
  teste("HCI (uma faceta): variâncias .0237/.0270/.1885 e Φ = .6874 com 20 itens", {
    data(HCIlong, package = "ShinyItemAnalysis", envir = environment())
    hc <- quieto(generalizabilidade(HCIlong, "id", "rating", "item"))
    all(abs(hc$tabela$variancia - c(.0237, .0270, .1885)) < 5e-5) && abs(hc$efeito$valor[2] - .6874) < 5e-5 })
} else pular("HCI de Martinková", "ShinyItemAnalysis não instalado")
teste("dados balanceados: componentes iguais aos da ANOVA por quadrados médios (p × i)", {
  set.seed(7); np <- 50; ni <- 6; pe_ <- rnorm(np, 0, 1); ie <- rnorm(ni, 0, .5)
  wide <- as.data.frame(outer(pe_, ie, "+") + matrix(rnorm(np * ni, 0, .8), np)); names(wide) <- paste0("i", 1:ni); wide$id <- 1:np
  gw <- quieto(generalizabilidade(wide, "id", itens = paste0("i", 1:ni)))
  lg <- data.frame(p = factor(rep(1:np, ni)), i = factor(rep(1:ni, each = np)), y = unlist(wide[1:ni]))
  ms <- summary(aov(y ~ p + i, data = lg))[[1]][["Mean Sq"]]
  v_p <- (ms[1] - ms[3]) / ni; v_i <- (ms[2] - ms[3]) / np
  max(abs(gw$tabela$variancia - c(v_p, v_i, ms[3]))) < 1e-6 })
teste("estudo D: Φ cresce com o número de itens e Eρ² >= Φ", {
  gd <- quieto(generalizabilidade(pio, "Person", "Score", c("Item", "Occasion"), estudo_d = list(Item = c(2, 4, 8))))
  t <- gd$estudo_d[gd$estudo_d$Occasion == 2, ]; t <- t[order(t$Item), ]
  all(diff(t$phi) > 0) && all(gd$estudo_d$g >= gd$estudo_d$phi) })

cat("== analise_distratores (Martinková & Hladká, 2023, 5.5)\n")
if (tem("ShinyItemAnalysis")) {
  data(HCItest, HCIkey, package = "ShinyItemAnalysis", envir = environment())
  ad_ <- analise_distratores(HCItest, names(HCItest)[1:20], HCIkey$key)
  teste("proporções por grupo iguais às contagens do ShinyItemAnalysis::DistractorAnalysis (item 1)", {
    sia <- ShinyItemAnalysis::DistractorAnalysis(HCItest[, 1:20], HCIkey, num.groups = 3)[[1]]
    x <- ad_$tabela[ad_$tabela$item == "Item 1", ]
    max(abs(as.matrix(x[, c("grupo1", "grupo2", "grupo3")]) - sweep(unclass(sia), 2, colSums(sia), "/")[x$opcao, ])) < 1e-12 })
  teste("gabarito do item 1 discrimina (superior - inferior > 0) e o distrator C não", {
    x <- ad_$tabela[ad_$tabela$item == "Item 1", ]; x$discriminacao[x$gabarito] > 0 && x$discriminacao[x$opcao == "C"] < 0 })
} else pular("análise de distratores", "ShinyItemAnalysis não instalado")

cat("== dmacs (Nye & Drasgow, 2011)\n")
HS <- lavaan::HolzingerSwineford1939; mhs <- "visual =~ x1 + x2 + x3\n textual =~ x4 + x5 + x6\n speed =~ x7 + x8 + x9"
fc <- lavaan::cfa(mhs, data = HS, group = "school", std.lv = TRUE, meanstructure = TRUE)
dmc <- quieto(dmacs(fc, "Pasteur"))
if (tem("pinsearch")) {
  teste("dMACS igual ao pinsearch::dmacs com os mesmos parâmetros e DP combinado (1e-10)", {
    est <- lavaan::lavInspect(fc, "est")
    lam <- rbind(rowSums(est[[2]]$lambda), rowSums(est[[1]]$lambda)); nu <- rbind(est[[2]]$nu[, 1], est[[1]]$nu[, 1])
    max(abs(pinsearch::dmacs(nu, lam, pooled_item_sd = dmc$tabela$dp_combinado) - dmc$tabela$dmacs)) < 1e-10 })
} else pular("dMACS vs pinsearch", "pinsearch não instalado")
if (tem("dmacs")) {
  teste("igual ao dmacs::lavaan_dmacs do livro de Petersen (dMACS tol. .005; ΔMean por item 1e-5, o pacote integra numericamente)", {
    ref <- quieto(dmacs::lavaan_dmacs(fc, RefGroup = "Pasteur"))
    max(abs(dmc$tabela$dmacs - rowSums(ref$DMACS, na.rm = TRUE))) < .005 && max(abs(dmc$tabela$delta_media - rowSums(ref$ItemDeltaMean, na.rm = TRUE))) < 1e-5 })
} else pular("dMACS vs dmacs", "dmacs não instalado (arquivado no CRAN)")
teste("forma fechada = integral numérica com traço de média .5 e variância 1,7 no grupo focal", {
  eR <- list(l = .8, n = 0); eF <- list(l = .6, n = .3); k <- .5; v <- 1.7; z <- seq(-8, 8, by = .0005)
  num <- sqrt(sum(((eF$n + eF$l * (k + z * sqrt(v))) - (eR$n + eR$l * (k + z * sqrt(v))))^2 * dnorm(z)) * .0005)
  fe <- sqrt(((eF$n - eR$n) + (eF$l - eR$l) * k)^2 + (eF$l - eR$l)^2 * v)
  abs(num - fe) < 1e-6 })
teste("no modelo escalar (cargas e interceptos iguais) o dMACS é zero", {
  inv <- quieto(invariancia(mhs, HS, paste0("x", 1:9), "school")); max(quieto(dmacs(inv$fits$escalar, "Pasteur"))$tabela$dmacs) < 1e-6 })
teste("configural com std.lv: aviso de métrica fixa e proporção de viés não calculada", any(grepl("std.lv", dmc$avisos)) && all(is.na(dmc$escala$proporcao_vies)))

cat("== tri_dif(metodo = \"mh\")\n")
if (tem("difR")) {
  data(verbal, package = "difR", envir = environment()); vb <- verbal[, c(1:24, 26)]; vb$Gender <- factor(vb$Gender, levels = c(0, 1))
  teste("sem purificação: χ² de MH e αMH iguais ao difR::difMH (1e-10)", {
    r <- quieto(tri_dif(vb, names(vb)[1:24], "Gender", metodo = "mh", purificar = FALSE))
    ref <- difR::difMH(vb[, 1:24], group = vb$Gender, focal.name = 1, purify = FALSE)
    max(abs(r$tabela$qui2_mh - ref$MH)) < 1e-10 && max(abs(r$tabela$alfa_mh - ref$alphaMH)) < 1e-10 })
  teste("com purificação (p sem correção): igual ao difR::difMH(purify = TRUE)", {
    r <- quieto(tri_dif(vb, names(vb)[1:24], "Gender", metodo = "mh", purificar = TRUE, corrigir = "none"))
    ref <- difR::difMH(vb[, 1:24], group = vb$Gender, focal.name = 1, purify = TRUE)
    max(abs(r$tabela$qui2_mh - ref$MH)) < 1e-8 })
  teste("ΔMH = -2,35 ln αMH e classes A/B/C pela régua do difR", {
    r <- quieto(tri_dif(vb, names(vb)[1:24], "Gender", metodo = "mh", purificar = FALSE))
    all(abs(r$tabela$delta_mh + 2.35 * log(r$tabela$alfa_mh)) < 1e-12) && all((abs(r$tabela$delta_mh) >= 1.5) == grepl("^C", r$tabela$classe)) })
} else pular("MH vs difR", "difR não instalado")

cat("== regressao_logistica(validar = TRUE)\n")
set.seed(5); nl <- 150; XL <- as.data.frame(matrix(rnorm(nl * 8), nl)); XL$y <- rbinom(nl, 1, plogis(-.3 + .8 * XL$V1 + .4 * XL$V2))
rv <- quieto(regressao_logistica(XL, "y", paste0("V", 1:8), evento = 1, validar = TRUE, reamostras = 300, box_tidwell = FALSE))
teste("no próprio banco a inclinação de calibração vale 1 e o intercepto 0 (máxima verossimilhança)",
      abs(rv$validacao$aparente[3] - 1) < 1e-8 && abs(rv$validacao$aparente[4]) < 1e-8)
teste("com 8 preditores e 150 casos o otimismo é positivo: AUC e inclinação corrigidas menores", rv$validacao$corrigido[1] < rv$validacao$aparente[1] && rv$validacao$corrigido[3] < 1)
if (tem("rms")) {
  teste("perto do rms::validate (AUC tol. .02, inclinação tol. .08; reamostras diferentes)", {
    set.seed(1); f <- rms::lrm(y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8, data = XL, x = TRUE, y = TRUE); v <- rms::validate(f, B = 300)
    abs(rv$validacao$corrigido[1] - (v["Dxy", "index.corrected"] + 1) / 2) < .02 && abs(rv$validacao$corrigido[3] - v["Slope", "index.corrected"]) < .08 })
} else pular("validação vs rms", "rms não instalado")

cat("== qui_quadrado 2 × 2: RR, DR e NNT\n")
d22 <- data.frame(g = rep(c("controle", "tto"), c(100, 100)), y = c(rep(c("sim", "nao"), c(40, 60)), rep(c("sim", "nao"), c(25, 75))))
q22 <- quieto(qui_quadrado(d22, "g", "y", evento = "sim", referencia = "controle"))
rk <- q22$tabelas_extra[[grep("^Risco", names(q22$tabelas_extra))]]
teste("RR = .25/.40, DR = -.15 e NNT = 6,7", abs(rk$valor[1] - .625) < 1e-12 && abs(rk$valor[2] + .15) < 1e-12 && abs(rk$valor[3] - 1 / .15) < 1e-12)
if (tem("DescTools")) {
  teste("IC da DR igual ao Newcombe (DescTools::BinomDiffCI, method = \"score\")", {
    b <- DescTools::BinomDiffCI(25, 100, 40, 100, method = "score"); abs(rk$ic_inf[2] - b[2]) < 1e-6 && abs(rk$ic_sup[2] - b[3]) < 1e-6 })
} else pular("Newcombe vs DescTools", "DescTools não instalado")

cat("== afe_numero_fatores(mrfa = TRUE)\n")
if (tem("EFA.MRFA")) {
  set.seed(1)
  XM <- as.data.frame(cbind(psych::sim.congeneric(N = 300, loads = c(.7, .6, .6, .5, .7, .6), short = FALSE)$observed,
                            psych::sim.congeneric(N = 300, loads = c(.7, .6, .6, .5, .7, .6), short = FALSE)$observed)); names(XM) <- paste0("i", 1:12)
  an <- quieto(afe_numero_fatores(XM, names(XM), mrfa = TRUE))
  teste("PA-MRFA recupera 2 fatores simulados e entra na tabela e na frase", an$n_mrfa == 2 && any(grepl("PA-MRFA", an$apa)))
} else pular("PA-MRFA", "EFA.MRFA não instalado")

ok <- unlist(resultados)
cat(sprintf("\nevals_v16: %d/%d passaram\n", sum(ok), length(ok)))
if (!all(ok)) quit(status = 1)
