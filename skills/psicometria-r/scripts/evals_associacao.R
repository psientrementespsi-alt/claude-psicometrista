# evals_associacao.R, testes do módulo de associação (est_associacao.R)
#
# Rodar da pasta da skill:
#   PR_DADOS_JASP="/caminho/para/Data Library" Rscript scripts/evals_associacao.R
# PR_DADOS_JASP aponta para a biblioteca de dados do JASP (a pasta que contém
# "4. Regression" e "5. Frequencies"; qualquer pasta acima delas também serve, a busca
# é recursiva). Sem a variável, os testes de âncora do livro são pulados com aviso.
# Três tipos de teste, nesta ordem: (1) âncoras JASP 0.95.1 impressas no livro do
# usuário (Vol. 1-3, com a linha mdX:linha); (2) verdade conhecida em dados simulados
# com semente fixa; (3) equivalência com outra implementação. Tolerância declarada em
# cada teste: para âncora, metade da última casa impressa pelo JASP.
# Saída: uma linha PASSOU/FALHOU por teste e o resumo; código de saída 1 se algo falhar.

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
perto <- function(a, b, tol) isTRUE(all(abs(a - b) <= tol + 1e-12))
t_inicio <- Sys.time()

DJ <- Sys.getenv("PR_DADOS_JASP")
jasp <- function(arquivo) {
  if (!nzchar(DJ)) return(NULL)
  f <- list.files(DJ, recursive = TRUE, full.names = TRUE)
  f <- f[basename(f) == arquivo]
  if (!length(f)) return(NULL)
  utils::read.csv(f[1], check.names = FALSE)
}
pular <- function(banco) cat("PULADO  âncoras de", banco, "(PR_DADOS_JASP vazia ou banco não encontrado)\n")

# ============================================================ 1. Âncoras JASP (livro)
ex <- jasp("Exam Anxiety.csv")
if (is.null(ex)) pular("Exam Anxiety") else {
  cat("== Âncoras: Exam Anxiety (Vol. 2, caps. 8 a 10)\n")
  cr <- correlacoes(ex, c("Revise", "Exam", "Anxiety"))
  tb <- cr$tabela
  linha <- function(a, b) tb[tb$var1 == a & tb$var2 == b, ]
  teste("r Revise × Exam = .397, IC [.220; .548]",
        perto(c(linha("Revise", "Exam")$r, linha("Revise", "Exam")$r_ic_inf, linha("Revise", "Exam")$r_ic_sup), c(.397, .220, .548), .0005))
  teste("r Revise × Anxiety = -.709, IC [-.794; -.598]",
        perto(unlist(linha("Revise", "Anxiety")[c("r", "r_ic_inf", "r_ic_sup")]), c(-.709, -.794, -.598), .0005))
  teste("r Exam × Anxiety = -.441, IC [-.585; -.271]",
        perto(unlist(linha("Exam", "Anxiety")[c("r", "r_ic_inf", "r_ic_sup")]), c(-.441, -.585, -.271), .0005))
  teste("r² = 15.7%, 50.3%, 19.4%", perto(100 * tb$r2, c(15.7, 50.3, 19.4), .05))
  teste("Spearman Revise × Anxiety ρ = -.622 ao lado do Pearson -.709",
        perto(linha("Revise", "Anxiety")$rho, -.622, .0005) && "rho" %in% names(tb) && "r" %in% names(tb))
  teste("Shapiro-Wilk W = .804 (Revise) e .822 (Anxiety)",
        perto(cr$pressupostos$estatistica[c(1, 3)], c(.804, .822), .0005))
  teste("frase APA com 3 casas como o parágrafo aprovado do livro: r(101) = -.441, p < .001, IC 95% [-.585, -.271]",
        any(grepl("r(101) = -.441, p < .001, IC 95% [-.585, -.271]", cr$apa, fixed = TRUE)))
  teste("decimal = \",\": r(101) = -,441; p < ,001; IC 95% [-,585; -,271] (idêntico a)",
        any(grepl("r(101) = -,441; p < ,001; IC 95% [-,585; -,271]", correlacoes(ex, c("Exam", "Anxiety"), decimal = ",")$apa, fixed = TRUE)))

  rs <- regressao_linear(ex, "Exam", "Revise")
  ms <- rs$tabelas_extra[[1]]; cs <- rs$tabela
  teste("simples: R = .397, R² = .157, R² aj = .149", perto(c(ms$R, ms$R2, ms$R2_ajustado), c(.397, .157, .149), .0005))
  teste("simples: RMSE 23.93 vs DP do desfecho 25.94; F(1, 101) = 18.87",
        perto(c(ms$RMSE, ms$RMSE_nulo, ms$F), c(23.93, 25.94, 18.87), .005) && ms$gl1 == 1 && ms$gl2 == 101)
  teste("simples: intercepto 45.321; B = .567, EP .130, β .397, t 4.343, IC [.308; .826]",
        perto(cs$B, c(45.321, .567), .0005) && perto(cs$EP[2], .130, .0005) && perto(cs$beta[2], .397, .0005) &&
          perto(cs$t[2], 4.343, .0005) && perto(c(cs$ic_inf[2], cs$ic_sup[2]), c(.308, .826), .0005))

  rm2 <- regressao_linear(ex, "Exam", c("Revise", "Anxiety"))
  mm <- rm2$tabelas_extra[[1]]; cm <- rm2$tabela
  teste("múltipla: F(2, 100) = 13.18; R² = .209; aj = .193; RMSE 23.31; f² = .264",
        perto(mm$F, 13.18, .005) && perto(c(mm$R2, mm$R2_ajustado, mm$f2), c(.209, .193, .264), .0005) && perto(mm$RMSE, 23.31, .005))
  teste("múltipla: Anxiety B = -.485 [-.863; -.107], β = -.321, p = .012",
        perto(unlist(cm[3, c("B", "ic_inf", "ic_sup", "beta", "p")]), c(-.485, -.863, -.107, -.321, .012), .0005))
  teste("múltipla: Revise B = .241 [-.116; .599], β = .169, p = .184",
        perto(unlist(cm[2, c("B", "ic_inf", "ic_sup", "beta", "p")]), c(.241, -.116, .599, .169, .184), .0005))
  teste("múltipla: tolerância .497, VIF 2.01", perto(cm$tolerancia[2], .497, .0005) && perto(cm$vif[2], 2.01, .005))
  pr <- rm2$pressupostos
  teste("múltipla: Durbin-Watson = 2.00", perto(pr$estatistica[grepl("Durbin", pr$checagem)], 2.00, .005))
  teste("múltipla: Shapiro-Wilk dos resíduos W = .965, p = .008",
        perto(pr$estatistica[grepl("Normalidade", pr$checagem)], .965, .0005) && perto(pr$p[grepl("Normalidade", pr$checagem)], .008, .0005))
  teste("frase APA da múltipla: F(2, 100) = 13.18, p < .001, 19.3% (R² ajustado = .193, f² = 0.264)",
        grepl("F(2, 100) = 13.18, p < .001, explicando 19.3% da variância de Exam (R² ajustado = .193, f² = 0.264)", rm2$apa[1], fixed = TRUE))
  teste("β com zero à esquerda e B, IC e β com 3 casas, como e: B = -0.485, IC 95% [-0.863, -0.107], β = -0.321",
        any(grepl("Anxiety: B = -0.485, IC 95% [-0.863, -0.107], β = -0.321, t(100) = -2.55, p = .012", rm2$apa, fixed = TRUE)) &&
          any(grepl("β = -0,321", regressao_linear(ex, "Exam", c("Revise", "Anxiety"), decimal = ",")$apa, fixed = TRUE)))

  hi <- regressao_hierarquica(ex, "Exam", list("Gender", "Revise", "Anxiety"))
  ps <- hi$tabela; cf <- hi$tabelas_extra[[1]]
  teste("hierárquica M1: R² = .000, F(1, 101) = 0.00, p = .963, aj = -.010, RMSE 26.07",
        perto(ps$R2[2], 0, .0005) && perto(ps$F_mudanca[2], 0, .005) && perto(ps$p[2], .963, .0005) &&
          perto(ps$R2_ajustado[2], -.010, .0005) && perto(ps$RMSE[2], 26.07, .005) && perto(ps$RMSE[1], 25.94, .005))
  teste("hierárquica M2: R² = ΔR² = .159, F(1, 100) = 18.89, f² do bloco = .189",
        perto(c(ps$R2[3], ps$delta_R2[3], ps$f2_bloco[3]), c(.159, .159, .189), .0005) && perto(ps$F_mudanca[3], 18.89, .005) && ps$gl2[3] == 100)
  teste("hierárquica M3: ΔR² = .050, F(1, 99) = 6.28, p = .014, f² = .063",
        perto(c(ps$delta_R2[4], ps$p[4], ps$f2_bloco[4]), c(.050, .014, .063), .0005) && perto(ps$F_mudanca[4], 6.28, .005))
  teste("hierárquica final: R² = .209, aj = .185, F(3, 99) = 8.72",
        perto(c(ps$R2[4], ps$R2_ajustado[4]), c(.209, .185), .0005) && grepl("F(3, 99) = 8.72", paste(hi$apa, collapse = " "), fixed = TRUE))
  teste("hierárquica final: Anxiety -.481 [-.863; -.100] β -.319; Revise .246 [-.116; .608] β .172",
        perto(unlist(cf[4, c("B", "ic_inf", "ic_sup", "beta", "p")]), c(-.481, -.863, -.100, -.319, .014), .0005) &&
          perto(unlist(cf[3, c("B", "ic_inf", "ic_sup", "beta", "p")]), c(.246, -.116, .608, .172, .180), .0005))
  teste("hierárquica final: Gender B = 1.040 [-8.183; 10.262], β = .020, p = .823 (recodificado 0/1)",
        perto(unlist(cf[2, c("B", "ic_inf", "ic_sup", "beta", "p")]), c(1.040, -8.183, 10.262, .020, .823), .0005))
  teste("hierárquica: β de Revise no degrau 2 = .400 (só consulta)", perto(hi$coeficientes_degraus[[2]]$beta[3], .400, .0005))
}

bf <- jasp("Big Five Personality Traits.csv")
if (is.null(bf)) pular("Big Five") else {
  cat("== Âncoras: Big Five (Vol. 1, cap. 6)\n")
  c5 <- correlacoes(bf, c("Neuroticism", "Extraversion"))$tabela
  teste("Neuroticism × Extraversion: r = -.350 [-.425; -.271]", perto(unlist(c5[c("r", "r_ic_inf", "r_ic_sup")]), c(-.350, -.425, -.271), .0005))
  teste("Neuroticism × Extraversion: ρ = -.325 [-.402; -.245]",
        perto(unlist(c5[c("rho", "rho_ic_inf", "rho_ic_sup")]), c(-.325, -.402, -.245), .0005))
  ca <- correlacoes(bf, c("Extraversion", "Agreeableness"))$tabela
  teste("Extraversion × Agreeableness: r = .05, p = .223; ρ = .06, p = .166",
        perto(c(ca$r, ca$rho), c(.05, .06), .005) && perto(c(ca$r_p, ca$rho_p), c(.223, .166), .0005))
}

ti <- jasp("Titanic.csv")
if (is.null(ti)) pular("Titanic") else {
  cat("== Âncoras: Titanic (Vol. 3, caps. 11 e 12)\n")
  lg <- regressao_logistica(ti, "Survived", c("Age", "PClass", "Sex"))
  md <- lg$tabelas_extra[["Modelo vs nulo"]]; cl <- lg$tabela; de <- lg$tabelas_extra[["Desempenho"]]
  teste("logística: N = 756, 557 excluídos", lg$n == 756 && lg$n_excluidos == 557)
  teste("logística: deviance nula 1025.57, modelo 695.14, χ²(4) = 330.43",
        perto(c(md$deviance_nulo, md$deviance_modelo, md$qui2), c(1025.57, 695.14, 330.43), .005) && md$gl == 4)
  ef <- lg$efeito
  teste("logística: Nagelkerke .477, Tjur .394",
        perto(ef$valor[ef$indice == "R² de Nagelkerke"], .477, .0005) && perto(ef$valor[ef$indice == "R² de Tjur"], .394, .0005))
  teste("logística: Age B = -.039, z = -5.14, OR = .962 [.947; .976]",
        perto(cl$B[2], -.039, .0005) && perto(cl$z[2], -5.14, .005) && perto(unlist(cl[2, c("OR", "OR_ic_inf", "OR_ic_sup")]), c(.962, .947, .976), .0005))
  teste("logística: 2nd OR = .275 [.165; .457]; 3rd OR = .080 [.047; .138]",
        perto(unlist(cl[3, c("OR", "OR_ic_inf", "OR_ic_sup")]), c(.275, .165, .457), .0005) &&
          perto(unlist(cl[4, c("OR", "OR_ic_inf", "OR_ic_sup")]), c(.080, .047, .138), .0005))
  teste("logística: male z = -13.06, OR = .072 [.048; .107]; inverso 13.9",
        perto(cl$z[5], -13.06, .005) && perto(unlist(cl[5, c("OR", "OR_ic_inf", "OR_ic_sup")]), c(.072, .048, .107), .0005) && perto(1 / cl$OR[5], 13.9, .05))
  cf <- lg$tabelas_extra[["Matriz de confusão"]]
  teste("logística: matriz de confusão 372/71 e 91/222 no corte .5",
        identical(as.numeric(unlist(cf[1, 2:3])), c(372, 71)) && identical(as.numeric(unlist(cf[2, 2:3])), c(91, 222)))
  teste("logística: acurácia .786 vs taxa base .586; sensibilidade .709; especificidade .840; AUC .853",
        perto(c(de$acuracia, de$taxa_base, de$sensibilidade, de$especificidade, de$AUC), c(.786, .586, .709, .840, .853), .0005))
  teste("logística: 3 casas por padrão nas frases, como: OR = 0.275, IC 95% [0.165, 0.457]; OR = 0.072",
        any(grepl("OR = 0.275, IC 95% [0.165, 0.457]", lg$apa, fixed = TRUE)) && any(grepl("OR = 0.072, IC 95% [0.048, 0.107]", lg$apa, fixed = TRUE)))
  teste("logística: frase APA traduz OR < 1 em % menor (72.5% para 2nd), nunca \"vezes menos\"",
        any(grepl("72.5% menor", lg$apa, fixed = TRUE)) && !any(grepl("vezes menos", lg$apa)))
  teste("logística: célula vazia de PClass tratada como ausente (sem categoria fantasma)",
        !any(grepl("PClass:  vs|PClass: vs", cl$termo)) && nrow(cl) == 5)

  q2 <- qui_quadrado(ti, "Sex", "Survived")
  tt <- q2$tabela; ef2 <- q2$efeito
  teste("2×2: χ²(1) = 332.06 lido sem correção; Yates 329.84 ao lado",
        perto(tt$estatistica[1], 332.06, .005) && tt$lido[1] && perto(tt$estatistica[2], 329.84, .005))
  teste("2×2: φ = -.503 (sinal do JASP), magnitude .50", perto(ef2$valor[ef2$indice == "φ"], -.503, .0005))
  teste("2×2: OR homens/mulheres = 0.100 [0.077; 0.130] como o JASP",
        perto(unlist(ef2[grepl("^OR \\(", ef2$indice), c("valor", "ic_inf", "ic_sup")]), c(.100, .077, .130), .0005))
  teste("2×2: OR mulheres/homens = 9.99 [7.67; 13.01] como o livro",
        perto(unlist(ef2[grepl("^OR inverso", ef2$indice), c("valor", "ic_inf", "ic_sup")]), c(9.99, 7.67, 13.01), .005))
  ch <- q2$tabelas_extra[[grep("^Chance", names(q2$tabelas_extra))]]
  teste("2×2: chances 2.000 e 0.200, probabilidades 66.7% e 16.7%; menor esperada 158.34",
        perto(ch$chance, c(2.000, .200), .0005) && perto(100 * ch$probabilidade, c(66.7, 16.7), .05) && perto(q2$pressupostos$estatistica, 158.34, .005))
  teste("2×2: correção do erro (0,200 é a chance dos homens; o OR é 0,100, 90% menor)",
        any(grepl("90.0% menor", q2$apa, fixed = TRUE)))
  teste("2×2: φ com IC bilateral igual ao effectsize::phi (sinal do JASP preservado) e na frase: φ = .503, IC 95% [",
        { ep <- effectsize::phi(table(ti$Sex, ti$Survived), adjust = FALSE, alternative = "two.sided")
          f <- ef2[ef2$indice == "φ", ]
          perto(c(f$ic_inf, f$ic_sup), -c(ep$CI_high, ep$CI_low), 1e-10) && any(grepl("φ = .503, IC 95% [", q2$apa, fixed = TRUE)) })
}

hh <- jasp("Health Habits.csv")
if (is.null(hh)) pular("Health Habits") else {
  cat("== Âncoras: Health Habits (Vol. 3, cap. 12)\n")
  hh[[1]] <- factor(hh[[1]], c("Low", "Moderate", "Vigorous")); hh[[2]] <- factor(hh[[2]], c("Low", "Medium", "High"))
  q3 <- qui_quadrado(hh, "Physical Activity", "Fruit Consumption")
  teste("3×3: χ²(4) = 14.15, p = .007", perto(q3$tabela$estatistica[1], 14.15, .005) && q3$tabela$gl[1] == 4 && perto(q3$tabela$p[1], .007, .0005))
  teste("3×3: V de Cramér = .077; V·√2 = .109 para comparar às faixas",
        perto(q3$efeito$valor[1], .077, .0005) && perto(q3$efeito$valor[2], .109, .0005))
  ce <- q3$tabelas_extra[[1]]
  rl <- function(a, b) ce$residuo_ajustado[ce[[1]] == a & ce[[2]] == b]
  teste("3×3: resíduos ajustados Low/Low +3.45 e Low/High -2.99; demais dentro de |1.96|",
        perto(c(rl("Low", "Low"), rl("Low", "High")), c(3.45, -2.99), .005) && sum(abs(ce$residuo_ajustado) > 1.96) == 2)
  teste("3×3: menor esperada 26.82", perto(q3$pressupostos$estatistica, 26.82, .005))
  teste("3×3: V = .077 com 3 casas e IC bilateral (= effectsize::cramers_v) na frase", {
    ev <- effectsize::cramers_v(table(hh[[1]], hh[[2]]), adjust = FALSE, alternative = "two.sided")
    perto(c(q3$efeito$ic_inf[1], q3$efeito$ic_sup[1]), c(ev$CI_low, ev$CI_high), 1e-10) &&
      any(grepl("V de Cramér = .077, IC 95% [", q3$apa, fixed = TRUE)) && any(grepl("w = V·√2 = .109", q3$apa, fixed = TRUE)) })
  ad <- qui_quadrado_aderencia(hh, "Physical Activity")
  teste("aderência: χ²(2) = 358.07, esperado 394.67 por categoria; observados 108/443/633",
        perto(ad$tabela$estatistica[1], 358.07, .005) && perto(ad$tabelas_extra[[1]]$esperado, rep(394.67, 3), .005) &&
          identical(ad$tabelas_extra[[1]]$observado, c(108L, 443L, 633L)))
}

# ============================================================ 2. Verdade conhecida (simulação)
cat("== Verdade conhecida (dados simulados, semente 2026)\n")
set.seed(2026)
n <- 10000
z <- rnorm(n); x <- .5 * z + sqrt(.75) * rnorm(n)
sim <- data.frame(z = z, x = x, ex = exp(x), u = rnorm(n))
cs <- correlacoes(sim, c("z", "x"))$tabela
# ρ de Spearman populacional na normal bivariada com ρ = .50: (6/π)·asen(.25) = .483
teste("correlação populacional .50: IC de r contém .50; ρ de Spearman perto de .483 (tol .03)",
      cs$r_ic_inf < .5 && cs$r_ic_sup > .5 && abs(cs$rho - 6 / pi * asin(.25)) < .03)
cm <- correlacoes(sim, c("z", "ex"), limiar_divergencia = .05)
teste("transformação monotônica (exp): ρ inalterado, Pearson cai e a divergência é avisada (limiar .05)",
        perto(cm$tabela$rho, cs$rho, 1e-10) && cm$tabela$r < cs$r - .05 && any(grepl("divergem", cm$avisos)))
teste("variável independente: |r| < .06 e magnitude desprezível", { t0 <- correlacoes(sim, c("z", "u"))$tabela; abs(t0$r) < .06 && t0$magnitude == "desprezível" })

set.seed(2026)
zc <- rnorm(800); pc <- data.frame(z = zc, a = zc + rnorm(800), b = zc + rnorm(800))
pp <- correlacao_parcial(pc, "a", "b", "z")
teste("parcial: a e b só compartilham z: bivariada ≈ .50, parcial ≈ 0 (|r| < .08)",
      abs(pp$tabela$coeficiente[1] - .5) < .06 && abs(pp$tabela$coeficiente[2]) < .08)

set.seed(2026)
nr <- 3000; x1 <- rnorm(nr); x2 <- .5 * x1 + rnorm(nr); g3 <- factor(sample(c("a", "b", "c"), nr, TRUE))
yr <- 1 + 2 * x1 - 1 * x2 + c(a = 0, b = 1, c = -1)[as.character(g3)] + rnorm(nr)
dr <- data.frame(y = yr, x1 = x1, x2 = x2, g3 = g3)
rr <- regressao_linear(dr, "y", c("x1", "x2", "g3"))
teste("regressão: B recuperados (2, -1, +1 e -1 para as indicadoras de g3) dentro do IC",
      all(rr$tabela$ic_inf[-1] < c(2, -1, 1, -1) & rr$tabela$ic_sup[-1] > c(2, -1, 1, -1)))
teste("regressão: fator de 3 níveis vira 2 indicadoras contra a referência 'a', sem β", sum(grepl("^g3: ", rr$tabela$termo)) == 2 &&
        all(grepl("vs a$", rr$tabela$termo[grepl("^g3", rr$tabela$termo)])) && all(is.na(rr$tabela$beta[grepl("^g3", rr$tabela$termo)])))
rref <- regressao_linear(dr, "y", c("x1", "g3"), referencias = list(g3 = "c"))
teste("regressão: referencias = list(g3 = \"c\") troca a referência sem mudar o modelo",
      all(grepl("vs c$", rref$tabela$termo[grepl("^g3", rref$tabela$termo)])) &&
        perto(rref$tabelas_extra[[1]]$R2, summary(lm(y ~ x1 + g3, dr))$r.squared, 1e-10))
set.seed(2026)
xh <- runif(1500, 0, 3); yh <- 1 + xh + rnorm(1500, sd = .2 + xh)
dh <- data.frame(y = yh, x = xh)
rh <- regressao_linear(dh, "y", "x", ep_robusto = TRUE)
teste("heterocedasticidade plantada: Breusch-Pagan p < .001 e aviso para EP robusto",
      rh$pressupostos$p[grepl("Breusch", rh$pressupostos$checagem)] < .001 && any(grepl("heterocedasticidade", rh$avisos)))
teste("outlier plantado: 1 caso com |z| > 3 e Cook > 1", {
  set.seed(2026); xo <- rnorm(60); yo <- xo + rnorm(60, sd = .5); xo[60] <- 6; yo[60] <- -8
  ro <- regressao_linear(data.frame(y = yo, x = xo), "y", "x")
  pr <- ro$pressupostos
  pr$estatistica[grepl("\\|z\\| > 3", pr$checagem)] >= 1 && any(grepl("Cook > 1", ro$avisos)) && "60" %in% ro$extra_diagnostico$cook_maior_1 })

set.seed(2026)
nl <- 4000; xl <- rnorm(nl); gl <- factor(sample(c("ctrl", "trat"), nl, TRUE))
eta <- -.5 + log(2) * xl + log(3) * (gl == "trat")
dl <- data.frame(y = rbinom(nl, 1, plogis(eta)), x = xl, g = gl)
lr <- regressao_logistica(dl, "y", c("x", "g"))
teste("logística: OR verdadeiros 2 (x) e 3 (trat vs ctrl) dentro do IC de Wald",
      lr$tabela$OR_ic_inf[2] < 2 && lr$tabela$OR_ic_sup[2] > 2 && lr$tabela$OR_ic_inf[3] < 3 && lr$tabela$OR_ic_sup[3] > 3)
teste("logística: logito linear não é acusado pelo Box-Tidwell (p > .01)",
      lr$pressupostos$p[grepl("Box-Tidwell", lr$pressupostos$checagem)] > .01)
set.seed(2026)
xq <- runif(3000, .1, 10); dq <- data.frame(y = rbinom(3000, 1, plogis(-1 + 2 * log(xq))), x = xq)
teste("logística: logito linear em ln(x), não em x, é acusado pelo Box-Tidwell (p < .001)", {
  pq <- regressao_logistica(dq, "y", "x")$pressupostos; pq$p[grepl("Box-Tidwell", pq$checagem)] < .001 })
teste("logística: separação completa é avisada", {
  ds <- data.frame(y = rep(0:1, each = 20), x = c(rnorm(20, -3, .3), rnorm(20, 3, .3)))
  any(grepl("separação", regressao_logistica(ds, "y", "x")$avisos)) })
teste("logística: evento por rótulo de texto e erro claro para desfecho com 3 categorias", {
  dl2 <- dl; dl2$y <- ifelse(dl$y == 1, "sim", "nao")
  r2 <- regressao_logistica(dl2, "y", "x", evento = "sim")
  dl3 <- dl2; dl3$y[1:5] <- "talvez"
  perto(r2$tabela$B, regressao_logistica(dl, "y", "x")$tabela$B, 1e-8) &&
    inherits(try(regressao_logistica(dl3, "y", "x"), silent = TRUE), "try-error") })

set.seed(2026)
yr0 <- rbinom(500, 1, .4); sr <- yr0 + rnorm(500)
dro <- data.frame(y = yr0, s = sr)
ro <- curva_roc(dro, "y", "s", cortes = c(0, .5))
teste("ROC: AUC próxima da teórica Φ(1/√2) = .760", abs(ro$efeito$valor - pnorm(1 / sqrt(2))) < .04)
teste("ROC: maior_indica_evento = FALSE com escore invertido dá a mesma AUC",
      perto(curva_roc(transform(dro, s = -s), "y", "s", maior_indica_evento = FALSE)$efeito$valor, ro$efeito$valor, 1e-10))
teste("ROC: J de Youden é o máximo da varredura de todos os cortes", {
  cc <- sort(unique(sr)); J <- vapply(cc, function(k) mean(sr[yr0 == 1] >= k) + mean(sr[yr0 == 0] < k) - 1, 0)
  perto(ro$tabela$youden_J[1], max(J), 1e-10) })

set.seed(2026)
pre <- rbinom(200, 1, .3); pos <- ifelse(pre == 1, rbinom(200, 1, .8), rbinom(200, 1, .5))
fu <- ifelse(pos == 1, rbinom(200, 1, .9), rbinom(200, 1, .1))
dmc <- data.frame(pre = pre, pos = pos, fu = fu)
mc <- mcnemar(dmc, "pre", "pos")
teste("McNemar: proporção sobe (.30 → ~.59) e o teste detecta (p < .001)", mc$tabela$p[1] < .001 && mc$efeito$valor[2] > .2)
qc <- q_cochran(dmc, c("pre", "pos", "fu"))
teste("Q de Cochran: mudança plantada detectada; pos vs fu (sem mudança plantada na média) p ajustado > .05",
      qc$tabela$p < .001 && qc$posthoc$p_mcnemar_ajustado[3] > .05)

# ============================================================ 3. Equivalência com outras implementações
cat("== Equivalência com outras implementações\n")
teste("Pearson e Spearman iguais ao cor.test (estimativa, p, IC do Pearson)", {
  ct <- cor.test(sim$z, sim$x); cs2 <- suppressWarnings(cor.test(sim$z, sim$x, method = "spearman"))
  perto(c(cs$r, cs$r_p, cs$r_ic_inf, cs$r_ic_sup), c(ct$estimate, ct$p.value, ct$conf.int), 1e-12) && perto(c(cs$rho, cs$rho_p), c(cs2$estimate, cs2$p.value), 1e-12) })
teste("Kendall: τ e p iguais ao cor.test; IC de Hollander-Wolfe confere com o cálculo direto", {
  dk <- sim[1:150, ]; tk <- correlacoes(dk, c("z", "x"), metodos = "kendall")$tabela
  ck <- suppressWarnings(cor.test(dk$z, dk$x, method = "kendall"))
  nn <- 150; Cv <- sapply(1:nn, function(i) sum(sign(dk$z[i] - dk$z[-i]) * sign(dk$x[i] - dk$x[-i])))
  s2 <- (2 * (nn - 2) * var(Cv) / (nn * (nn - 1)) + 1 - ck$estimate^2) * 2 / (nn * (nn - 1))
  perto(c(tk$tau, tk$tau_p), c(ck$estimate, ck$p.value), 1e-12) && perto(tk$tau_ic_inf, ck$estimate - qnorm(.975) * sqrt(s2), 1e-10) })
teste("ponto-bisserial: t equivalente = t de Student do t.test(var.equal = TRUE)", {
  set.seed(2026); dp <- data.frame(y = rnorm(90), g = rep(c("a", "b"), 45)); dp$y[dp$g == "b"] <- dp$y[dp$g == "b"] + .5
  rp <- correlacao_ponto_bisserial(dp, "y", "g")
  perto(rp$tabela$t_student_equivalente, -unname(t.test(y ~ g, dp, var.equal = TRUE)$statistic), 1e-10) })
teste("parcial: igual ao ppcor e à correlação entre resíduos", {
  rres <- cor(resid(lm(a ~ z, pc)), resid(lm(b ~ z, pc)))
  perto(pp$tabela$coeficiente[2], rres, 1e-10) && perto(pp$tabela$p[2], ppcor::pcor.test(pc$a, pc$b, pc$z)$p.value, 1e-12) })
teste("correlações dependentes: âncora do curso (N = 390, r = -.653 vs -.719, r23 = .518): |z| = 2.057, p = .040", {
  cd <- comparar_correlacoes(-.653, 390, -.719, r12 = .518)
  perto(abs(cd$tabela$z[1]), 2.057, .0005) && perto(cd$tabela$p[1], .040, .0005) })
if (requireNamespace("cocor", quietly = TRUE)) {
  teste("correlações dependentes: Steiger (1980) e Meng et al. (1992) iguais ao cocor", {
    cd <- comparar_correlacoes(.40, 250, .25, r12 = .30)
    co <- cocor::cocor.dep.groups.overlap(.40, .25, .30, 250, test = c("steiger1980", "meng1992"))
    perto(cd$tabela$z, c(co@steiger1980$statistic, co@meng1992$statistic), 1e-8) })
  teste("correlações independentes: z de Fisher e IC de Zou iguais ao cocor", {
    ci <- comparar_correlacoes(.5, 100, .3, 120)
    co <- cocor::cocor.indep.groups(.5, .3, 100, 120, test = c("fisher1925", "zou2007"))
    perto(ci$tabela$z, co@fisher1925$statistic, 1e-8) && perto(c(ci$tabela$ic_inf, ci$tabela$ic_sup), co@zou2007$conf.int, 1e-8) })
} else cat("PULADO  equivalência com cocor (pacote não instalado)\n")
teste("regressão: B, EP, p e R² iguais ao lm(); VIF igual ao car::vif", {
  fit <- lm(y ~ x1 + x2, dr); ra <- regressao_linear(dr, "y", c("x1", "x2"))
  perto(ra$tabela$B, coef(fit), 1e-10) && perto(ra$tabela$p, summary(fit)$coefficients[, 4], 1e-12) &&
    perto(ra$tabela$vif[-1], car::vif(fit), 1e-10) })
teste("regressão: VIF de fator com 3 níveis = GVIF^(1/gl) do car", {
  v <- car::vif(lm(y ~ x1 + x2 + g3, dr)); perto(unique(rr$tabela$vif[grepl("^g3", rr$tabela$termo)]), v["g3", 1]^(1 / v["g3", 2]), 1e-10) })
teste("regressão: Durbin-Watson igual ao car e p ao lmtest::dwtest bilateral; Breusch-Pagan igual ao lmtest", {
  fit <- lm(y ~ x1 + x2, dr); ra <- regressao_linear(dr, "y", c("x1", "x2")); pr <- ra$pressupostos
  perto(pr$estatistica[grepl("Durbin", pr$checagem)], car::durbinWatsonTest(fit, simulate = FALSE)$dw, 1e-10) &&
    perto(pr$p[grepl("Durbin", pr$checagem)], lmtest::dwtest(fit, alternative = "two.sided")$p.value, 1e-10) &&
    perto(pr$estatistica[grepl("Breusch", pr$checagem)], unname(lmtest::bptest(fit)$statistic), 1e-10) })
teste("regressão: EP HC3 igual ao sandwich::vcovHC", perto(rh$tabela$EP_HC3, sqrt(diag(sandwich::vcovHC(lm(y ~ x, dh), type = "HC3"))), 1e-12))
teste("regressão: Mahalanobis igual ao mahalanobis() sobre os preditores", {
  X <- as.matrix(dr[, c("x1", "x2")]); D2 <- mahalanobis(X, colMeans(X), cov(X))
  pr <- regressao_linear(dr, "y", c("x1", "x2"))$pressupostos
  perto(pr$estatistica[grepl("Mahalanobis", pr$checagem)], max(D2), 1e-8) })
teste("hierárquica: F de mudança = anova(m1, m2) de dois modelos (não o anova de três)", {
  ex2 <- dr; h2 <- regressao_hierarquica(ex2, "y", list("x2", "g3", "x1"))
  m1 <- lm(y ~ x2, ex2); m2 <- lm(y ~ x2 + g3, ex2); m3 <- lm(y ~ x2 + g3 + x1, ex2)
  perto(h2$tabela$F_mudanca[3], anova(m1, m2)$F[2], 1e-10) && perto(h2$tabela$gl1[3], 2, 0) &&
    !perto(h2$tabela$F_mudanca[3], anova(m1, m2, m3)$F[2], 1e-6) })
teste("logística: B igual ao glm; IC Wald = confint.default; IC por perfil = confint", {
  fit <- glm(y ~ x + g, binomial, dl); lp <- regressao_logistica(dl, "y", c("x", "g"), ic = "perfil")
  perto(lp$tabela$B, coef(fit), 1e-10) && perto(log(lp$tabela$OR_ic_inf), confint.default(fit)[, 1], 1e-10) &&
    perto(log(lp$tabela$OR_ic_sup_perfil), suppressMessages(confint(fit))[, 2], 1e-8) })
if (requireNamespace("performance", quietly = TRUE)) {
  teste("logística: pseudo-R² de McFadden, Cox-Snell, Nagelkerke e Tjur iguais ao performance", {
    fit <- glm(y ~ x + g, binomial, dl); e <- lr$efeito
    perto(c(e$valor[e$indice == "R² de McFadden"], e$valor[e$indice == "R² de Cox e Snell"], e$valor[e$indice == "R² de Nagelkerke"], e$valor[e$indice == "R² de Tjur"]),
          unname(c(performance::r2_mcfadden(fit)$R2, performance::r2_coxsnell(fit), performance::r2_nagelkerke(fit), performance::r2_tjur(fit))), 1e-8) })
} else cat("PULADO  equivalência com performance (pacote não instalado)\n")
teste("ROC: AUC = U de Mann-Whitney/(n1·n0) e IC = pROC::ci.auc (DeLong)", {
  W <- wilcox.test(sr[yr0 == 1], sr[yr0 == 0])$statistic
  rr0 <- pROC::roc(yr0, sr, levels = c(0, 1), direction = "<", quiet = TRUE)
  perto(ro$efeito$valor, W / (sum(yr0) * sum(1 - yr0)), 1e-10) && perto(c(ro$efeito$ic_inf, ro$efeito$ic_sup), as.numeric(pROC::ci.auc(rr0))[c(1, 3)], 1e-10) })
teste("qui-quadrado: χ², p e resíduos iguais ao chisq.test(correct = FALSE); Yates igual ao padrão do R", {
  tbx <- table(dl$g, dl$y); qq <- qui_quadrado(dl, "g", "y")
  c1 <- chisq.test(tbx, correct = FALSE); c2 <- chisq.test(tbx)
  perto(qq$tabela$estatistica[1:2], c(c1$statistic, c2$statistic), 1e-10) && perto(qq$tabelas_extra[[1]]$residuo_ajustado, as.vector(c1$stdres), 1e-10) })
teste("qui-quadrado: V de Cramér igual ao effectsize::cramers_v(adjust = FALSE)", {
  dv <- data.frame(a = dr$g3, b = cut(dr$x1, 4))
  qv <- qui_quadrado(dv, "a", "b")
  perto(qv$efeito$valor[1], effectsize::cramers_v(table(dv$a, dv$b), adjust = FALSE)$Cramers_v, 1e-10) })
teste("qui-quadrado: esperadas pequenas levam ao Fisher exato (troca visível) e o p = fisher.test", {
  qf <- qui_quadrado(matrix(c(3, 1, 1, 4), 2))
  qf$tabela$lido[qf$tabela$teste == "Fisher exato"] && perto(qf$tabela$p[qf$tabela$teste == "Fisher exato"], fisher.test(matrix(c(3, 1, 1, 4), 2))$p.value, 1e-12) &&
    any(grepl("troca visível", qf$decisoes)) })
teste("qui-quadrado r×c com esperadas pequenas: Fisher exato do R e semente sem efeito colateral na sessão", {
  mt <- matrix(c(2, 0, 3, 1, 4, 0, 0, 2, 1), 3); set.seed(1); a1 <- runif(1)
  set.seed(1); qr <- qui_quadrado(mt); a2 <- runif(1)
  perto(qr$tabela$p[qr$tabela$lido], fisher.test(mt)$p.value, 1e-12) && a1 == a2 })
teste("aderência: igual ao chisq.test com proporções informadas", {
  ad2 <- qui_quadrado_aderencia(dr, "g3", proporcoes = c(a = .5, b = .25, c = .25))
  perto(ad2$tabela$estatistica[1], chisq.test(table(dr$g3), p = c(.5, .25, .25))$statistic, 1e-10) })
teste("McNemar: sem e com correção iguais ao mcnemar.test; exato = binom.test", {
  tm <- table(factor(pre, 0:1), factor(pos, 0:1))
  perto(mc$tabela$estatistica[1:2], c(mcnemar.test(tm, correct = FALSE)$statistic, mcnemar.test(tm)$statistic), 1e-10) &&
    perto(mc$tabela$p[3], binom.test(tm[1, 2], tm[1, 2] + tm[2, 1])$p.value, 1e-12) })
teste("Q de Cochran: igual ao rstatix::cochran_qtest", {
  lg <- data.frame(id = rep(1:200, 3), tempo = rep(c("pre", "pos", "fu"), each = 200), y = c(pre, pos, fu))
  rq <- rstatix::cochran_qtest(lg, y ~ tempo | id)
  perto(qc$tabela$estatistica, rq$statistic, 1e-8) })
teste("Q de Cochran: Σz² do pós-hoc tipo SPSS = k·Q/2 (identidade que os números do curso também cumprem: 53.2 = 3 × 35.471/2)",
      perto(sum(qc$posthoc$z_spss^2), 3 * qc$tabela$estatistica / 2, 1e-8) && perto(4.464^2 + 5.646^2 + 1.182^2, 3 * 35.471 / 2, .01))
teste("Q de Cochran com 2 medidas = McNemar sem correção", perto(q_cochran(dmc, c("pre", "pos"))$tabela$estatistica, mc$tabela$estatistica[1], 1e-10))


# ============================================================ 4. Regressões dos achados (CoVe v1.1)
cat("== Regressões dos achados de revisão (CoVe v1.1)\n")
xa <- factor(c(rep(1, 12), rep(2, 10), rep(3, 8)), levels = 1:4)   # face 4 nunca saiu
teste("B2 aderência: categoria com contagem 0 fica no teste: χ²(3) = 11.07, p = .011 (= chisq.test)", {
  a1 <- qui_quadrado_aderencia(data.frame(face = xa), "face")
  a2 <- qui_quadrado_aderencia(data.frame(face = as.character(xa)), "face", proporcoes = c(`1` = .25, `2` = .25, `3` = .25, `4` = .25))
  a3 <- qui_quadrado_aderencia(data.frame(face = xa), "face", proporcoes = rep(.25, 4))
  perto(a1$tabela$estatistica[1], 11.07, .005) && a1$tabela$gl[1] == 3 && perto(a1$tabela$p[1], chisq.test(table(xa))$p.value, 1e-12) &&
    perto(a2$tabela$estatistica[1], a1$tabela$estatistica[1], 1e-12) && perto(a3$tabela$estatistica[1], a1$tabela$estatistica[1], 1e-12) &&
    any(grepl("sem nenhuma observação", a1$decisoes)) })
set.seed(15)
grp <- sample(c("caso", "controle"), 200, TRUE, prob = c(.3, .7))
dro2 <- data.frame(grp = grp, esc = round(ifelse(grp == "caso", rnorm(200, 20, 5), rnorm(200, 14, 5))))
teste("B3 ROC: desfecho de texto sem evento dá erro que pede evento (não escolhe pela ordem alfabética)",
      grepl("informe evento", tryCatch({ curva_roc(dro2, "grp", "esc"); "" }, error = conditionMessage)))
teste("B3 ROC: fator usa o 2º nível como evento, igual à regressao_logistica", {
  df <- transform(dro2, grp = factor(grp, levels = c("controle", "caso")))
  r1 <- curva_roc(df, "grp", "esc"); r2 <- regressao_logistica(df, "grp", "esc")
  grepl("'caso'", r1$decisoes[1]) && grepl("'caso'", r2$decisoes[1]) && r1$efeito$valor > .7 })
teste("B3 ROC: evento invertido (AUC < .5) gera aviso e nenhuma frase de corte com ±Inf", {
  ri <- curva_roc(dro2, "grp", "esc", evento = "controle")
  ri$efeito$valor < .5 && any(grepl("conferir o evento|Conferir o evento", ri$avisos)) && !grepl("Inf", ri$apa) &&
    !grepl("Youden", ri$apa) && all(is.finite(ri$tabela$corte) | is.na(ri$tabela$corte)) })
set.seed(8)
mq <- data.frame(t1 = rbinom(60, 1, .2), t2 = rbinom(60, 1, .6), t3 = rbinom(60, 1, .7))
teste("B4 Cochran: pós-hoc respeita decimal = \",\" (separador \"; \") e idioma = \"en\" (\"adjusted p\")", {
  qv <- q_cochran(mq, c("t1", "t2", "t3"), decimal = ","); qe <- q_cochran(mq, c("t1", "t2", "t3"), idioma = "en")
  ph_pt <- qv$apa[length(qv$apa)]; ph_en <- qe$apa[length(qe$apa)]
  grepl("; p ajustado = ,", ph_pt, fixed = TRUE) && !grepl("[0-9], p ajustado", ph_pt) &&
    grepl(", adjusted p = .", ph_en, fixed = TRUE) && !grepl("ajustado", ph_en, fixed = TRUE) })
teste("B6 Kendall: rótulo próprio (decisão diferente do Pearson = \"divergiu\") e aviso", {
  set.seed(49); xk <- rnorm(25); yk <- .4 * xk + rnorm(25)
  rk <- correlacoes(data.frame(x = xk, y = yk), c("x", "y"), metodos = c("pearson", "spearman", "kendall"))
  pk <- suppressWarnings(cor.test(xk, yk, method = "kendall"))$p.value; pp <- cor.test(xk, yk)$p.value
  ((pk < .05) != (pp < .05)) && grepl("Kendall divergiu", rk$apa[1]) && any(grepl("^Kendall e Pearson", rk$avisos)) })
teste("B7 referencias aceita vetor nomeado c(g3 = \"c\") igual à lista; nome desconhecido dá erro claro", {
  r1 <- regressao_linear(dr, "y", c("x1", "g3"), referencias = c(g3 = "c")); r2 <- regressao_linear(dr, "y", c("x1", "g3"), referencias = list(g3 = "c"))
  identical(r1$tabela, r2$tabela) &&
    grepl("não é preditor", tryCatch({ regressao_linear(dr, "y", "x1", referencias = c(zz = "a")); "" }, error = conditionMessage)) })
teste("B8 McNemar e Cochran: fator usa o 2º nível como evento (não o alfabético); texto sem evento pede evento", {
  dm <- data.frame(pre = c(rep("sim", 42), rep("nao", 58)), pos = c(rep("sim", 30), rep("nao", 12), rep("sim", 5), rep("nao", 53)))
  df <- dm; df$pre <- factor(df$pre, c("sim", "nao")); df$pos <- factor(df$pos, c("sim", "nao"))
  m1 <- mcnemar(df, "pre", "pos"); q1 <- q_cochran(transform(df, fu = pos), c("pre", "pos", "fu"))
  grepl("'nao'", m1$decisoes[1]) && grepl("evento \\(código 1\\) = 'nao'", q1$decisoes[1]) &&
    grepl("informe evento", tryCatch({ mcnemar(dm, "pre", "pos"); "" }, error = conditionMessage)) &&
    grepl("'sim'", mcnemar(dm, "pre", "pos", evento = "sim")$decisoes[1]) })
teste("B9 separação: frase APA diz \"OR não estimável\" (Firth) e não imprime OR absurdo nem pseudo-R²", {
  set.seed(20); dsp <- data.frame(x = rnorm(60), g = sample(c("a", "b"), 60, TRUE)); dsp$y <- as.integer(dsp$x > 0)
  rs <- regressao_logistica(dsp, "y", c("x", "g"))
  ap <- paste(rs$apa, collapse = " ")
  grepl("OR não estimável", ap) && grepl("Firth", ap) && !grepl("[0-9]{8,}", ap) && !grepl("Nagelkerke = ", ap) })
teste("B11 casas adaptativas: B = -8.9e-06 sai com algarismos significativos, não 0.000", {
  set.seed(3); xr <- rnorm(300, 5000, 1500); dr2 <- data.frame(x = xr, y = 10 - 8.9e-06 * xr * 10 + rnorm(300, sd = .01))
  ap <- regressao_linear(dr2, "y", "x")$apa
  grepl("x: B = -0.0000", ap[2], fixed = TRUE) && !grepl("B = 0.000,", ap[2], fixed = TRUE) &&
    .ee_assoc_casas(-8.9e-06, 3) == 7 && .ee_assoc_casas(.4, 3) == 3 })
teste("célula pequena: qui_quadrado avisa para não publicar cruzamento demográfico com células < 5",
      any(grepl("agrupar categorias", qui_quadrado(matrix(c(3, 1, 1, 4), 2))$avisos)))

# ============================================================ formatação e privacidade
cat("== Formatação e privacidade\n")
teste("idioma = \"en\": rótulo 95% CI e frase em inglês",
      any(grepl("95% CI", correlacoes(sim, c("z", "x"), idioma = "en")$apa, fixed = TRUE)))
teste("nenhuma tabela impressa tem uma linha por participante", {
  objs <- list(rr, lr, ro, mc, qc, cs = correlacoes(sim, c("z", "x", "u")))
  all(vapply(objs, function(o) all(vapply(c(list(o$tabela, o$efeito, o$pressupostos, o$posthoc), o$tabelas_extra),
                                          function(t) is.null(t) || NROW(t) < 50, logical(1))), logical(1))) })

# ============================================================ 5. Correções da verificação final (CoVe final, G e H)
cat("== Correções da verificação final (achados G e H)\n")
# G-03: escolaridade gravada 1/2/3, efeito não monotônico (grupo 2 mais alto). Dados da
# reprodução cove_final/G/g_rep_G03_codigo_numerico.R. Números de referência: lm() e
# glm() com factor(escolaridade), calculados aqui mesmo.
set.seed(7); n_g3 <- 150
d_g3 <- data.frame(escolaridade = sample(1:3, n_g3, TRUE), idade = round(rnorm(n_g3, 35, 8)))
d_g3$escore <- 10 + 2 * (d_g3$escolaridade == 2) - 1 * (d_g3$escolaridade == 3) + .1 * d_g3$idade + rnorm(n_g3, 0, 3)
d_g3$bin <- rbinom(n_g3, 1, .5)
aviso_cod <- function(r, v) any(grepl(paste0(v, " é numérica com só"), r$avisos, fixed = TRUE) &
                                  grepl(paste0("fatores = \"", v, "\""), r$avisos, fixed = TRUE))
r_g3a <- regressao_linear(d_g3, "escore", c("escolaridade", "idade"))
teste("G-03 escolaridade 1/2/3 sem fatores: entra com 1 coeficiente linear e $avisos manda declarar fatores = \"escolaridade\" (livro)",
      aviso_cod(r_g3a, "escolaridade") && sum(grepl("^escolaridade", r_g3a$tabela$termo)) == 1)
r_g3b <- regressao_linear(d_g3, "escore", c("escolaridade", "idade"), fatores = "escolaridade")
lm_g3 <- lm(escore ~ factor(escolaridade) + idade, d_g3)
teste("G-03 com fatores = \"escolaridade\": 2 indicadoras contra o código 1, B iguais ao lm(factor()) (tol 1e-10), 2 vs 1 = 2.786 (G-03) e sem o aviso",
      identical(r_g3b$tabela$termo[2:3], c("escolaridade: 2 vs 1", "escolaridade: 3 vs 1")) &&
        perto(r_g3b$tabela$B, unname(coef(lm_g3)), 1e-10) && perto(r_g3b$tabela$B[2], 2.786, .0005) &&
        !aviso_cod(r_g3b, "escolaridade") && any(grepl("numérica declarada em fatores", r_g3b$decisoes, fixed = TRUE)))
teste("G-03 hierárquica com fatores: F de mudança do bloco escolaridade = anova(lm(idade), lm(factor + idade)), F(2, 146) = 18.10 (G-03)", {
  h_g3 <- regressao_hierarquica(d_g3, "escore", list("idade", "escolaridade"), fatores = "escolaridade")
  a_g3 <- anova(lm(escore ~ idade, d_g3), lm_g3)
  perto(h_g3$tabela$F_mudanca[3], a_g3$F[2], 1e-8) && h_g3$tabela$gl1[3] == 2 && h_g3$tabela$gl2[3] == 146 &&
    perto(h_g3$tabela$F_mudanca[3], 18.10, .005) && !aviso_cod(h_g3, "escolaridade") &&
    aviso_cod(regressao_hierarquica(d_g3, "escore", list("idade", "escolaridade")), "escolaridade") })
teste("G-03 fatores respeita referencias (= relevel do lm); referência de numérica sem fatores e fatores fora dos preditores dão erro claro", {
  r_ref <- regressao_linear(d_g3, "escore", c("escolaridade", "idade"), fatores = "escolaridade", referencias = list(escolaridade = 2))
  lm_ref <- lm(escore ~ relevel(factor(escolaridade), "2") + idade, d_g3)
  identical(r_ref$tabela$termo[2:3], c("escolaridade: 1 vs 2", "escolaridade: 3 vs 2")) && perto(r_ref$tabela$B, unname(coef(lm_ref)), 1e-10) &&
    grepl("fatores = \"escolaridade\"", tryCatch({ regressao_linear(d_g3, "escore", c("escolaridade", "idade"), referencias = list(escolaridade = 2)); "" }, error = conditionMessage), fixed = TRUE) &&
    grepl("não é preditor", tryCatch({ regressao_linear(d_g3, "escore", "idade", fatores = "escolaridade"); "" }, error = conditionMessage)) })
teste("G-03 logística: aviso sem fatores; com fatores, B iguais ao glm(factor()) (tol 1e-8)", {
  set.seed(71); d_g3$ev <- rbinom(n_g3, 1, plogis(-1 + 1.2 * (d_g3$escolaridade == 2)))
  l1 <- regressao_logistica(d_g3, "ev", c("escolaridade", "idade")); l2 <- regressao_logistica(d_g3, "ev", c("escolaridade", "idade"), fatores = "escolaridade")
  aviso_cod(l1, "escolaridade") && !aviso_cod(l2, "escolaridade") &&
    perto(l2$tabela$B, unname(coef(glm(ev ~ factor(escolaridade) + idade, binomial, d_g3))), 1e-8) })
teste("G-03 sem aviso para preditor contínuo de verdade (idade, ~40 valores) nem para 0/1; fronteiras da regra: 3 e 7 inteiros avisam, 2, 8 e valores com decimais não", {
  r_ok <- regressao_linear(d_g3, "escore", c("idade", "bin"))
  !any(grepl("valores inteiros distintos", r_ok$avisos)) && length(unique(d_g3$idade)) > 7 &&
    length(.ee_assoc_aviso_codigo(rep(1:3, 5), "v", "")) == 1 && length(.ee_assoc_aviso_codigo(rep(1:7, 5), "v", "")) == 1 &&
    length(.ee_assoc_aviso_codigo(rep(0:1, 5), "v", "")) == 0 && length(.ee_assoc_aviso_codigo(rep(1:8, 5), "v", "")) == 0 &&
    length(.ee_assoc_aviso_codigo(rep(c(1.5, 2, 3), 5), "v", "")) == 0 })

# G-04: magnitude no idioma da frase. Referência: r pelo cor() e rótulo pelas faixas de
# Cohen (.10/.30/.50) aplicadas aqui, fora da função.
set.seed(2026); x_g4 <- rnorm(300); d_g4 <- data.frame(x = x_g4, y = .45 * x_g4 + rnorm(300))
r_g4 <- cor(d_g4$x, d_g4$y)
rot_g4 <- if (abs(r_g4) < .1) "negligible" else if (abs(r_g4) < .3) "weak" else if (abs(r_g4) < .5) "moderate" else "strong"
teste("G-04 correlacoes(idioma = \"en\"): magnitude em inglês na frase e na tabela (r = cor(), rótulo por Cohen)", {
  c_en <- correlacoes(d_g4, c("x", "y"), idioma = "en")
  grepl(paste0(", ", rot_g4, " in magnitude by Cohen (1988)"), c_en$apa[1], fixed = TRUE) && c_en$tabela$magnitude == rot_g4 &&
    !grepl("moderada|forte|fraca|desprezível", c_en$apa[1]) &&
    correlacoes(d_g4, c("x", "y"))$tabela$magnitude == c(negligible = "desprezível", weak = "fraca", moderate = "moderada", strong = "forte")[[rot_g4]] })
if (!is.null(ex))
  teste("G-04 âncora: Anxiety × Exam em inglês, \"moderate in magnitude\", r(101) = -.441",
        grepl("negatively associated, moderate in magnitude by Cohen (1988), r(101) = -.441", correlacoes(ex, c("Anxiety", "Exam"), idioma = "en")$apa[1], fixed = TRUE))

# G-13: 2 × 2 com célula zero (controle sem nenhum 'sim'). OR corrigido conferido à mão:
# (6,5/1,5)/(0,5/8,5) = 73,667; sem correção, 6·8/(1·0) = infinito.
d_g13 <- data.frame(g = rep(c("controle", "tratado"), c(8, 7)), ev = c(rep("nao", 8), rep("sim", 6), "nao"))
teste("G-13 célula zero: frase traz OR = 73.667 = (6.5/1.5)/(0.5/8.5) e diz que é corrigido (Haldane) e que o observado seria infinito; $decisoes também", {
  q13 <- qui_quadrado(d_g13, "g", "ev"); or_mao <- (6.5 / 1.5) / (.5 / 8.5)
  f13 <- q13$apa[2]
  perto(q13$efeito$valor[grepl("^OR \\(", q13$efeito$indice)], or_mao, 1e-10) && grepl("OR = 73.667", f13, fixed = TRUE) &&
    grepl("correção de Haldane, +0.5 em todas as células", f13, fixed = TRUE) && grepl("seria infinito", f13, fixed = TRUE) &&
    any(grepl("Haldane", q13$decisoes) & grepl("infinito", q13$decisoes)) &&
    grepl("Haldane-corrected OR, +0,5", qui_quadrado(d_g13, "g", "ev", idioma = "en", decimal = ",")$apa[2], fixed = TRUE) &&
    !any(grepl("Haldane", qui_quadrado(dl, "g", "y")$apa)) })

# G-11: rótulo do bloco no idioma; alfa da aderência = 1 − nivel. Contagens 30/18/12 com
# esperado 20: χ²(2) = (100 + 4 + 64)/20 = 8.4 e p = exp(−8.4/2) = .015 (gl = 2).
teste("G-11 hierárquica com idioma = \"en\" e blocos sem nome: \"(block 1)\", nada de \"bloco\"", {
  h_en <- regressao_hierarquica(dr, "y", list("x2", "x1"), idioma = "en")
  grepl("adding x2 (block 1)", h_en$apa[1], fixed = TRUE) && !any(grepl("bloco", h_en$apa)) })
teste("G-11 aderência decide com alfa = 1 − nivel: p = .015 (= exp(−8.4/2)) difere com nivel .95 e não difere com nivel .99", {
  d_ad <- data.frame(cat = rep(c("a", "b", "c"), c(30, 18, 12)))
  a95 <- qui_quadrado_aderencia(d_ad, "cat"); a99 <- qui_quadrado_aderencia(d_ad, "cat", nivel = .99)
  perto(a95$tabela$p[1], exp(-8.4 / 2), 1e-12) && perto(a95$tabela$estatistica[1], 8.4, 1e-12) &&
    grepl("diferiu da esperada", a95$apa) && !grepl("não diferiu", a95$apa) && grepl("não diferiu da esperada", a99$apa) })

# H-01: ponto-bisserial com grupo de 1 caso (escore 97.531, dados de cove_final/H/priv2.R).
set.seed(7)
d_h1 <- data.frame(g = c(rep("b", 30), "c"), x = round(rnorm(31, 50, 10), 3)); d_h1$x[31] <- 97.531
teste("H-01 ponto-bisserial com grupo de n = 1: o print não mostra o escore (97.531), a tabela de grupos vai em tabelas_extra com coluna n e $avisos registra", {
  pb_h1 <- correlacao_ponto_bisserial(d_h1, "x", "g")
  txt <- paste(capture.output(print(pb_h1)), collapse = "\n")
  tg <- pb_h1$tabelas_extra[["Codificação e grupos"]]
  !grepl("97.53", txt, fixed = TRUE) && !grepl("97.53", paste(pb_h1$apa, collapse = " "), fixed = TRUE) &&
    is.data.frame(tg) && "n" %in% names(tg) && identical(as.numeric(tg$n), c(30, 1)) && perto(tg$media[2], 97.531, 1e-12) &&
    any(grepl("n < 5", pb_h1$avisos) & grepl("'c', n = 1", pb_h1$avisos, fixed = TRUE)) &&
    !any(grepl("n < 5", correlacao_ponto_bisserial(rbind(d_h1, data.frame(g = "c", x = c(60, 61, 62, 63))), "x", "g")$avisos)) })

# H-06: índices do modelo sem IC ganham a nota em $decisoes (regra 6).
teste("H-06 R² ajustado, ΔR² e f² sem IC: nota \"não calculado\" em $decisoes da linear e da hierárquica (também em inglês)", {
  any(grepl("R² ajustado e f² sem IC: não calculado", regressao_linear(dr, "y", c("x1", "x2"))$decisoes, fixed = TRUE)) &&
    any(grepl("R² ajustado, ΔR² e f² sem IC: não calculado", regressao_hierarquica(dr, "y", list("x2", "x1"))$decisoes, fixed = TRUE)) &&
    any(grepl("Adjusted R² and f² without CI: not computed", regressao_linear(dr, "y", "x1", idioma = "en")$decisoes, fixed = TRUE)) })

# H-07: nível de fator com 1 caso. Referência da proporção |z| > 2,5: rstandard() do lm
# com na.rm (o caso de alavanca 1 tem resíduo padronizado NaN).
set.seed(7)
d_h7 <- data.frame(g = c(rep("a", 30), rep("b", 30), "c"), x = rnorm(61, 50, 10), cov = rnorm(61))
teste("H-07 nível de fator com n = 1: erro claro com o nome da categoria (não \"missing value where TRUE/FALSE needed\")", {
  msg <- tryCatch({ regressao_linear(d_h7, "x", c("cov", "g")); "" }, error = conditionMessage)
  grepl("'c', n = 1", msg, fixed = TRUE) && grepl("menos de 2 casos", msg) && !grepl("missing value", msg) })
teste("H-07 diagnóstico com caso de alavanca 1 não quebra: proporção |z| > 2,5 = mean(|rstandard| > 2,5, na.rm = TRUE)", {
  fit_h7 <- lm(x ~ cov + g, d_h7); dg_h7 <- .ee_assoc_diag_lm(fit_h7)
  zr <- rstandard(fit_h7)
  any(is.nan(zr)) && perto(dg_h7$pressupostos$estatistica[grepl("2,5", dg_h7$pressupostos$checagem)], mean(abs(zr) > 2.5, na.rm = TRUE), 1e-12) })

# H-08: documentação do zero à esquerda do β igual à base (est_comum.R lista β com zero).
if (!file.exists(file.path(pasta, "..", "references", "associacao.md"))) {
  cat("PULADO  H-08 (references/ fora da pasta dos scripts, como na instalação em ~/.local/share)\n")
} else teste("H-08 associacao.md: célula do zero à esquerda do β diz \"mantido (como a base, est_comum.R)\" e est_comum.R põe β entre os com zero", {
  # leitura em bytes (useBytes): o teste roda também em sessão com locale C
  md <- readLines(file.path(pasta, "..", "references", "associacao.md")); co <- readLines(file.path(pasta, "est_comum.R"))
  tem <- function(txt, x) grepl(txt, x, fixed = TRUE, useBytes = TRUE)
  any(tem("| Zero à esquerda do β |", md) & tem("mantido (como a base, est_comum.R)", md)) &&
    !any(tem("lista β sem zero", md)) && any(tem("β padronizado com zero", co)) })

cat("== Equivalência da correlação e aviso de não significativo (v1.3)\n")
set.seed(1311)
xa <- rnorm(200); ya <- rnorm(200); ya <- stats::resid(stats::lm(ya ~ xa))  # r amostral exatamente 0
ca <- data.frame(x = xa, y = ya, w = xa + rnorm(200, 0, .5))
ns_txt <- "Não significativo não quer dizer igual"
teste("A-01 correlacoes() com r = 0: aviso cita Greenland et al. (2016) e aponta teste_equivalencia(tipo = \"correlacao\")",
      local({ a <- correlacoes(ca, c("x", "y"))$avisos
        any(grepl(ns_txt, a, fixed = TRUE)) && any(grepl("Greenland et al., 2016", a, fixed = TRUE)) &&
          any(grepl("teste_equivalencia(tipo = \"correlacao\")", a, fixed = TRUE)) }))
teste("A-02 correlacoes() com associação forte (x e w): sem o aviso",
      !any(grepl(ns_txt, correlacoes(ca, c("x", "w"))$avisos, fixed = TRUE)))
teste("A-03 com três variáveis, o aviso lista só os pares não significativos",
      local({ a <- correlacoes(ca, c("x", "y", "w"))$avisos
        v <- a[grepl(ns_txt, a, fixed = TRUE)]
        length(v) == 1 && grepl("x × y (r = .000", v, fixed = TRUE) && !grepl("x × w", v, fixed = TRUE) }))
teste("A-03b r numericamente zero (1e-17) sai como .000, sem dezenas de casas na frase e no aviso",
      local({ r <- correlacoes(ca, c("x", "y"))
        !any(grepl("[0-9][.,][0-9]{6,}", c(r$apa, r$avisos))) }))
if (!requireNamespace("TOSTER", quietly = TRUE)) cat("PULADO  equivalência com TOSTER::z_cor_test (pacote ausente)\n") else
  teste("A-04 TOST da correlação = TOSTER::z_cor_test (z de Fisher, Lakens, 2017, eq. 8): p, z do limite decisivo e IC 90%",
        local({ cb <- data.frame(x = rnorm(80)); cb$y <- .15 * cb$x + rnorm(80)
          r <- teste_equivalencia(cb, "correlacao", limites = c(-.3, .3), x = "x", y = "y")
          z <- TOSTER::z_cor_test(cb$x, cb$y, alternative = "equivalence", null = c(-.3, .3))
          perto(r$p_tost, z$p.value, 1e-8) && perto(max(abs(r$tabela$estatistica[1:2][which.max(r$tabela$p[1:2])])), abs(unname(z$statistic)), 1e-8) &&
            perto(c(r$efeito$ic_inf[1], r$efeito$ic_sup[1]), as.numeric(z$conf.int), 1e-8) }))
teste("A-05 verdade conhecida: rho = 0 com n = 2000 e limites ±0,1 dá equivalência; com n = 30 o IC 90% é mais largo que os limites",
      local({ g <- data.frame(x = rnorm(2000)); g$y <- rnorm(2000)
        r1 <- teste_equivalencia(g, "correlacao", limites = .1, x = "x", y = "y")
        r2 <- teste_equivalencia(g[1:30, ], "correlacao", limites = .1, x = "x", y = "y")
        r1$equivalente && !r2$equivalente && any(grepl("mais largo que o intervalo de equivalência", r2$avisos)) }))
teste("A-06 qui-quadrado não significativo: aviso sem TOST (a skill não tem equivalência para V ou φ)",
      local({ tb <- data.frame(a = rep(c("s", "n"), each = 100), b = rep(rep(c("x", "z"), each = 50), 2))
        av <- qui_quadrado(tb, "a", "b")$avisos
        any(grepl(ns_txt, av, fixed = TRUE)) && any(grepl("A skill não tem teste de equivalência para esta estatística", av, fixed = TRUE)) }))
teste("A-07 textos novos sem travessão nem palavras vetadas",
      local({ cb <- data.frame(x = rnorm(80)); cb$y <- .15 * cb$x + rnorm(80)
        rs <- list(correlacoes(ca, c("x", "y", "w")), teste_equivalencia(cb, "correlacao", limites = .3, x = "x", y = "y"),
                   teste_equivalencia(cb, "correlacao", limites = .3, x = "x", y = "y", idioma = "en"))
        txt <- unlist(lapply(rs, function(x) c(x$apa, x$decisoes, x$avisos)))
        !any(grepl(enc2utf8("\u2014|\u2013|genuinamente|de fato|notavelmente|importante notar|intricado|meticuloso|cabe ressaltar|vale destacar|essencialmente|fundamentalmente|abrangente"),
                   enc2utf8(txt), useBytes = TRUE)) }))

cat("== v1.4: ROC com prevalência e custo; viés preditivo de Cleary\n")
set.seed(9); roc_n <- 400; roc_y <- stats::rbinom(roc_n, 1, .25); roc_s <- stats::rnorm(roc_n, 10 + 2 * roc_y, 3)
roc_d <- data.frame(y = roc_y, s = roc_s)
rr <- curva_roc(roc_d, "y", "s", prevalencia = .05, custo_fn_fp = 5, cortes = c(11, 12))
teste("ROC: VPP e VPN por Bayes na prevalência da amostra = contagens da tabela 2 x 2 nos cortes informados (tol. 1e-12)", {
  inf <- rr$tabela[rr$tabela$origem == "informado", ]; cc <- rr$classificacao_cortes
  all(abs(inf$vpp_amostra - cc$vpp) < 1e-12) && all(abs(inf$vpn_amostra - cc$vpn) < 1e-12) })
teste("ROC: VPP na prevalência de 5% = Se·π / [Se·π + (1 − Sp)(1 − π)] e RV+ = Se / (1 − Sp)", {
  t1 <- rr$tabela[1, ]; vpp <- t1$sensibilidade * .05 / (t1$sensibilidade * .05 + (1 - t1$especificidade) * .95)
  abs(t1$vpp_populacao - vpp) < 1e-12 && abs(t1$rv_pos - t1$sensibilidade / (1 - t1$especificidade)) < 1e-12 })
teste("ROC: corte ponderado por custo = máximo de Se + r·Sp, r = (1 − π)/(custo·π), em força bruta sobre os cortes do pROC", {
  ro <- pROC::roc(roc_y, roc_s, levels = c(0, 1), direction = "<", quiet = TRUE); r <- .95 / (5 * .05)
  abs(rr$tabela$corte[grepl("custo", rr$tabela$origem)] - ro$thresholds[which.max(ro$sensitivities + r * ro$specificities)]) < 1e-12 })
teste("ROC sem prevalência: decisão pede prevalencia para a população de uso; prevalência fora de (0, 1) para", {
  r0 <- curva_roc(roc_d, "y", "s")
  any(grepl("informar prevalencia", r0$decisoes)) && grepl("entre 0 e 1", tryCatch({ curva_roc(roc_d, "y", "s", prevalencia = 1.2); "" }, error = conditionMessage)) })
set.seed(10); vp_n <- 300; vp_g <- rep(c("A", "B"), each = vp_n / 2); vp_x <- stats::rnorm(vp_n, 50, 10)
vp_d <- data.frame(g = vp_g, x = vp_x, y0 = 5 + .5 * vp_x + stats::rnorm(vp_n, 0, 5),
                   y1 = 5 + .5 * vp_x + ifelse(vp_g == "B", -4, 0) + stats::rnorm(vp_n, 0, 5),
                   y2 = 5 + ifelse(vp_g == "B", .2, .5) * vp_x + stats::rnorm(vp_n, 0, 5))
teste("Cleary: sem viés simulado -> sem evidência; intercepto deslocado -> interceptos; inclinação diferente -> inclinações", {
  vies_preditivo(vp_d, "y0", "x", "g")$conclusao == "sem_evidencia" && vies_preditivo(vp_d, "y1", "x", "g")$conclusao == "interceptos" &&
    vies_preditivo(vp_d, "y2", "x", "g")$conclusao == "inclinacoes" })
teste("Cleary: F das três etapas = anova() entre os modelos aninhados à mão (tol. 1e-10)", {
  r <- vies_preditivo(vp_d, "y1", "x", "g"); xc <- vp_d$x - mean(vp_d$x)
  m0 <- lm(vp_d$y1 ~ xc); m1 <- lm(vp_d$y1 ~ xc + vp_d$g); m2 <- lm(vp_d$y1 ~ xc * vp_d$g)
  all(abs(r$tabela$F - c(anova(m0, m2)$F[2], anova(m1, m2)$F[2], anova(m0, m1)$F[2])) < 1e-10) })
teste("Cleary: retas por grupo iguais ao lm em cada grupo; resíduo médio da reta comum soma zero ponderado por n", {
  r <- vies_preditivo(vp_d, "y1", "x", "g"); rt <- r$tabelas_extra[[1]]; mB <- lm(y1 ~ x, vp_d[vp_d$g == "B", ])
  abs(rt$inclinacao[rt$grupo == "B"] - coef(mB)[2]) < 1e-10 && abs(sum(rt$residuo_medio_reta_comum * rt$n)) < 1e-8 })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
if (!nzchar(DJ)) cat("Aviso: PR_DADOS_JASP vazia; as âncoras do livro não foram testadas.\n")
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
