# evals_poder.R, testes do módulo de poder e tamanho amostral (est_poder.R)
#
# Rodar da pasta da skill:  Rscript scripts/evals_poder.R
# Saída: uma linha PASSOU/FALHOU por teste e o resumo. Código de saída 1 se algum
# teste falhar. Não usa dado nenhum: os gabaritos são números impressos.
# Ordem de valor: (1) âncoras do livro do usuário (md1, cap. 7); (2) exemplos do curso
# conferidos contra o G*Power (spec_T3, seção 25) e o exemplo de Faul et al. (2007);
# (3) equivalência com pwr e pwrss; a seção (6) traz os testes das correções da revisão
# v1.1.0, ancorados nos exemplos do manual do G*Power 3.1 (versão de 1/6/2023) ou em
# conta independente feita no próprio teste. Tolerância declarada em cada teste: N
# inteiro tem de bater exato; poder impresso com 3 casas, tolerância de meia casa (.0005).

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
campo <- function(res, nome) res$tabela$valor[res$tabela$campo == nome][1]

cat("== (1) Âncoras do livro: Vol. 1, cap. 7\n")
teste("r = .35, poder .80: G*Power exato 61", poder_correlacao(.35)$n == 61)
teste("r = .35, poder .90: G*Power exato 81", poder_correlacao(.35, poder = .90)$n == 81)
teste("r = .10, poder .80: G*Power exato 782", poder_correlacao(.10)$n == 782)
teste("r = .35, poder .80: Fisher 62", poder_correlacao(.35, metodo = "fisher")$n == 62)
teste("r = .35, poder .90: Fisher 82", poder_correlacao(.35, poder = .90, metodo = "fisher")$n == 82)
teste("r = .10, poder .80: Fisher 783", poder_correlacao(.10, metodo = "fisher")$n == 783)
teste("r = .10 com N = 500: poder 61% nos dois métodos", {
  a <- quieto(poder_correlacao(.10, n = 500)); b <- quieto(poder_correlacao(.10, n = 500, metodo = "fisher"))
  abs(as.numeric(campo(a, "Actual power")) - .61) < .005 && abs(as.numeric(campo(b, "Actual power")) - .61) < .005 })
teste("r = .35 com N = 500: poder praticamente 1",
      as.numeric(campo(quieto(poder_correlacao(.35, n = 500)), "Actual power")) > .999)
teste("conta de Fisher do livro: 61,77 antes de arredondar",
      abs(((qnorm(.975) + qnorm(.80)) / atanh(.35))^2 + 3 - 61.77) < .005)
teste("perda pelo livro: inflar 10% a 20% leva 62 a 69 e 75",
      ajustar_perda(62, .10)$n_coleta_inflado == 69 && ajustar_perda(62, .20)$n_coleta_inflado == 75)
teste("frase APA segue o modelo do livro (efeito, α bicaudal, poder, N)", {
  a <- poder_correlacao(.35, metodo = "fisher", justificativa = "magnitude estimada a partir de achado anterior comparável")$apa
  grepl("ρ = .35", a, fixed = TRUE) && grepl("α = .05 (bicaudal)", a, fixed = TRUE) &&
    grepl("poder de .80", a, fixed = TRUE) && grepl("mínimo de 62 participantes", a, fixed = TRUE) })

cat("== (2) Exemplos do curso conferidos contra o G*Power (spec_T3, seção 25)\n")
teste("#1 correlação ρ = .16 -> 304", poder_correlacao(.16)$n == 304)
teste("#2 correlação ρ = .30 -> 84 (pwr::pwr.r.test dá 85: divergência documentada)",
      poder_correlacao(.30)$n == 84 && ceiling(pwr::pwr.r.test(r = .30, power = .80)$n) == 85)
teste("correlação ρ = .50 -> 29 e ρ = .30 unicaudal -> 67", poder_correlacao(.50)$n == 29 && poder_correlacao(.30, caudas = 1)$n == 67)
teste("#3 t independente, d de médias 15 e 18 com DP 6,15 -> 67 + 67 = 134", {
  r <- poder_teste_t(efeito_d_de_medias(15, 18, 6.15))
  r$n == 134 && campo(r, "Sample size group 1") == "67" && campo(r, "Sample size group 2") == "67" })
teste("#3 t independente, d = .488 digitado -> 134", poder_teste_t(.488)$n == 134)
teste("#4 t independente, d = .50 -> 128", poder_teste_t(.50)$n == 128)
teste("#5 t pareado, dz = .30, poder .90 -> 119", poder_teste_t(.30, "pareado", poder = .90)$n == 119)
teste("#6 t pareado, dz = .40, poder .80 -> 52", poder_teste_t(.40, "pareado")$n == 52)
teste("#7 one-way, f de médias 12, 12, 15 e DP 6,15 = .22995 (tol. 5e-5)",
      abs(efeito_f_de_medias(c(12, 12, 15), 6.15) - .22995) < 5e-5)
teste("#7 one-way, f = .23, 3 grupos -> 186 (o 189 da aula não se reproduz)", {
  a <- poder_anova(f = .23, niveis = 3); b <- poder_anova(f = efeito_f_de_medias(c(12, 12, 15), 6.15), niveis = 3)
  a$n == 186 && b$n == 186 && a$n != 189 })
teste("#8 one-way, f = .50, 3 grupos -> 42", poder_anova(f = .50, niveis = 3)$n == 42)
teste("one-way, f = .25, k = 3, 4, 5 -> 159, 180, 200 (G*Power arredonda para grupos iguais; mínimos 158, 179, 196)", {
  ns <- vapply(3:5, function(k) poder_anova(f = .25, niveis = k)$n, 0)
  mins <- vapply(3:5, function(k) as.numeric(campo(poder_anova(f = .25, niveis = k), "N mínimo inteiro")), 0)
  all(ns == c(159, 180, 200)) && all(mins == c(158, 179, 196)) })
teste("#9 MR intra, f = .33, 4 medidas, poder .90, rho = NA -> 132",
      quieto(poder_anova_mr(f = .33, medidas = 4, efeito = "intra", rho = NA, poder = .90))$n == 132)
teste("#9 MR intra com 'Number of groups' = 4 (como na aula) -> também 132",
      quieto(poder_anova_mr(f = .33, medidas = 4, grupos = 4, efeito = "intra", rho = NA, poder = .90))$n == 132)
# Nome e checagem corrigidos na revisão v1.1.0 (C11): o antigo dizia "3 casas" e o grepl de
# "0.333" também aceitava "0.3333". Com f = .333 o N seria 130 (conferido no
# pwrss::power.f.mixed.anova, rho.within = NA) e com .33 seria 132 (curso), então a frase
# precisa das 4 casas para devolver 129.
teste("#9 MR intra com f exato de ηp² = .10 (.3333) -> 129; a frase imprime f com 4 casas (.333 daria 130, .33 daria 132)", {
  r <- quieto(poder_anova_mr(eta2p = .10, medidas = 4, efeito = "intra", rho = NA, poder = .90))
  n333 <- quieto(pwrss::power.f.mixed.anova(eta.squared = efeito_eta2_de_f(.333), factor.levels = c(1, 4), rho.within = NA,
                                            power = .9, effect = "within", verbose = 0))$n.total
  r$n == 129 && grepl("f = 0.3333 ", r$apa, fixed = TRUE) && n333 == 130 &&
    quieto(poder_anova_mr(f = .333, medidas = 4, rho = NA, poder = .90))$n == 130 &&
    quieto(poder_anova_mr(f = .33, medidas = 4, rho = NA, poder = .90))$n == 132 })
teste("MR intra, sensibilidade a ε: .75 -> 162 e .50 -> 217 (curso, seção 17)",
      quieto(poder_anova_mr(f = .33, medidas = 4, rho = NA, epsilon = .75, poder = .90))$n == 162 &&
        quieto(poder_anova_mr(f = .33, medidas = 4, rho = NA, epsilon = .50, poder = .90))$n == 217)
teste("#10 fatorial 3 × 5, ηp² = .06, gl 8, 15 células -> 244", poder_anova(eta2p = .06, niveis = c(3, 5))$n == 244)
teste("#10 fatorial com o f transcrito pela aula (.225) -> 305, logo .225 é erro de transcrição",
      poder_anova(f = .225, niveis = c(3, 5))$n == 305)
teste("fatorial 3 × 5: efeitos principais A -> 155 e B -> 193 (curso, seção 18)",
      poder_anova(eta2p = .06, niveis = c(3, 5), efeito = "A")$n == 155 &&
        poder_anova(eta2p = .06, niveis = c(3, 5), efeito = "B")$n == 193)
teste("#11 MR interação 3 × 2, f = .204, poder .90, rho = NA -> 309 (mínimo 308, múltiplo de 3)", {
  r <- quieto(poder_anova_mr(f = .204, medidas = 2, grupos = 3, efeito = "interacao", rho = NA, poder = .90))
  r$n == 309 && campo(r, "N mínimo inteiro") == "308" })
teste("#12 ANCOVA interação 3 × 5 com 1, 2 e 3 covariáveis -> 244",
      all(vapply(1:3, function(cv) poder_anova(eta2p = .06, niveis = c(3, 5), covariaveis = cv)$n, 0) == 244))
teste("#13 ANCOVA, efeito da covariável (gl = 1) com 1, 2 e 3 covariáveis -> 126",
      all(vapply(1:3, function(cv) poder_anova(eta2p = .06, niveis = c(3, 5), covariaveis = cv, gl_efeito = 1)$n, 0) == 126))
teste("#14 regressão global, R² = .20 (f² = .25), 5 preditores -> 58", poder_regressao("global", 5, r2 = .20)$n == 58)
teste("#15 coeficiente, R² parcial = .0625 (f² = .0667), 5 preditores -> 120",
      poder_regressao("coeficiente", 5, r2_parcial = .0625)$n == 120 && poder_regressao("coeficiente", 5, f2 = .067)$n == 120)
teste("#15 coeficiente com R² parcial = .067 (leitura literal da aula) -> 112",
      poder_regressao("coeficiente", 5, r2_parcial = .067)$n == 112)
teste("#16 ΔR² = .05 sobre R² = .20, 3 testados de 5 -> 168 (poder .90 -> 217)",
      poder_regressao("incremento", 5, r2_reduzido = .20, delta_r2 = .05, testados = 3)$n == 168 &&
        poder_regressao("incremento", 5, r2_reduzido = .20, delta_r2 = .05, testados = 3, poder = .90)$n == 217)
teste("#16 ΔR² com f² = .067 digitado -> 167", poder_regressao("incremento", 5, f2 = .067, testados = 3)$n == 167)
teste("#17 logística, OR = 1.52, Pr H0 = .30, Demidenko com correção -> 227", poder_logistica(1.52, .30)$n == 227)
# Âncora do R² trocada na revisão v1.1.0 (C4): o 284 vinha da "Sensibilidade" do spec_T3
# (linha 821), calculada com o pwrss mantendo a correção (227 / 0,8), e não de tela do
# G*Power. O manual do G*Power (seção 29.4, 395 com R² = .1) mostra que o programa usa o
# caminho sem correção quando R² > 0: 234,19 / 0,8 = 292,74 -> 293 (integração numérica
# de v1 com X ~ N(0, 1), β0 = logit(.30), β1 = ln 1,52, feita à parte; ver o teste da
# seção 6).
quieto_log <- function(...) quieto(poder_logistica(...))
teste("logística: Hsieh -> 214, Demidenko sem correção -> 235 (curso, seção 24); R² other X = .20 -> 293 (regra do manual, C4)",
      poder_logistica(1.52, .30, metodo = "hsieh")$n == 214 && poder_logistica(1.52, .30, metodo = "demidenko")$n == 235 &&
        quieto_log(1.52, .30, r2_outros = .20)$n == 293)
teste("logística: p0 = .05 -> 898 e .50 -> 196 (sensibilidade a p0, curso)",
      poder_logistica(1.52, .05)$n == 898 && poder_logistica(1.52, .50)$n == 196)
teste("dz a partir de d = .30 com r = .30, .50, .70 -> N = 166, 119, 72 (curso, seção 14)",
      all(vapply(c(.3, .5, .7), function(r) poder_teste_t(.30, "pareado", poder = .90, r_pareado = r)$n, 0) == c(166, 119, 72)))
teste("regra 50 + 8k diverge do G*Power com f² = .15: k = 4 -> 85 e k = 10 -> 118",
      poder_regressao("global", 4, f2 = .15)$n == 85 && poder_regressao("global", 10, f2 = .15)$n == 118)

cat("== (2b) Faul et al. (2007, pp. 13-14): exemplo de medidas repetidas do artigo do G*Power\n")
faul <- poder_anova_mr(f = .2572, medidas = 3, grupos = 3, efeito = "intra", rho = .3, epsilon = 1, n = 90)
teste("intra: λ = 25.52 (tol. .005)", abs(as.numeric(campo(faul, "Noncentrality parameter λ")) - 25.52) < .005)
teste("intra: poder .997 (tol. .0005), gl = 2 e 174, F crítico 3.048 (tol. .0005)",
      abs(as.numeric(campo(faul, "Actual power")) - .997) < .0005 && campo(faul, "Numerator df") == "2" &&
        campo(faul, "Denominator df") == "174" && abs(as.numeric(campo(faul, "Critical F")) - 3.048) < .0005)
teste("interação: f = .1532 -> poder .653 (tol. .0005)",
      abs(as.numeric(campo(poder_anova_mr(f = .1532, medidas = 3, grupos = 3, efeito = "interacao", rho = .3, n = 90),
                           "Actual power")) - .653) < .0005)
teste("entre: f = .1719571 -> poder .488 (tol. .0005)",
      abs(as.numeric(campo(poder_anova_mr(f = .1719571, medidas = 3, grupos = 3, efeito = "entre", rho = .3, n = 90),
                           "Actual power")) - .488) < .0005)
teste("Berti et al. no artigo: f = .25, N = 20, 2 grupos, 10 medidas, ρ = .5 -> A .30, B .95, A × B .95 (tol. .005)", {
  p <- function(ef) as.numeric(campo(poder_anova_mr(f = .25, medidas = 10, grupos = 2, efeito = ef, rho = .5, n = 20), "Actual power"))
  abs(p("entre") - .30) < .005 && abs(p("intra") - .95) < .005 && abs(p("interacao") - .95) < .005 })

cat("== (3) Equivalência com pwr e pwrss\n")
teste("t independente e pareado: N igual ao teto do pwr::pwr.t.test numa grade de d e poder", {
  ok <- TRUE
  for (d in c(.2, .35, .5, .8)) for (pw in c(.8, .9)) {
    ok <- ok && poder_teste_t(d, poder = pw)$n == 2 * ceiling(pwr::pwr.t.test(d = d, power = pw)$n) &&
      poder_teste_t(d, "pareado", poder = pw)$n == ceiling(pwr::pwr.t.test(d = d, power = pw, type = "paired")$n)
  }
  ok })
# Âncora trocada na revisão v1.1.0 (C5): antes n2 = teto(n1 × .85) com n1 inteiro (70 + 60);
# a regra do G*Power arredonda cada grupo a partir de n1* contínuo (n1* = 69,38 pela raiz
# de pwr::pwr.t2n.test com n2 = .85·n1 contínuo -> 70 + 59).
teste("t com grupos desiguais: poder igual ao pwr::pwr.t2n.test (tol. 1e-4: a tabela guarda 4 casas)", {
  r <- quieto(poder_teste_t(.5, razao = .85)); n1 <- as.numeric(campo(r, "Sample size group 1")); n2 <- as.numeric(campo(r, "Sample size group 2"))
  pw <- quieto(poder_teste_t(.5, n = c(n1, n2)))
  n1c <- stats::uniroot(function(x) pwr::pwr.t2n.test(n1 = x, n2 = .85 * x, d = .5)$power - .80, c(5, 1e4), tol = 1e-10)$root
  abs(as.numeric(campo(pw, "Actual power")) - pwr::pwr.t2n.test(n1 = n1, n2 = n2, d = .5)$power) < 1e-4 &&
    pwr::pwr.t2n.test(n1 = n1, n2 = n2, d = .5)$power >= .80 && n1 == ceiling(n1c) && n2 == ceiling(.85 * n1c) &&
    n1 == 70 && n2 == 59 })
teste("one-way: N igual a teto(pwr::pwr.anova.test) × k numa grade", {
  ok <- TRUE
  for (k in 2:5) for (f in c(.1, .25, .4)) ok <- ok && poder_anova(f = f, niveis = k)$n == k * ceiling(pwr::pwr.anova.test(k = k, f = f, power = .8)$n)
  ok })
teste("fatorial e ANCOVA: N mínimo igual ao teto do pwrss::power.f.ancova (ceil.n = FALSE)", {
  ok <- TRUE
  for (cv in 0:2) for (ef in c("A", "B", "A:B")) {
    x <- quieto(pwrss::power.f.ancova(eta.squared = .06, factor.levels = c(3, 5), target.effect = ef, k.covariates = cv,
                                      power = .8, ceil.n = FALSE, verbose = 0))
    ok <- ok && poder_anova(eta2p = .06, niveis = c(3, 5), efeito = ef, covariaveis = cv)$n == ceiling(x$n.total)
  }
  ok })
teste("medidas repetidas: poder igual ao pwrss::power.f.mixed.anova (modo padrão e rho = NA; tol. 1e-4)", {
  ok <- TRUE
  for (ef in c("within", "between", "interaction")) for (rho in c(.5, .3, NA)) for (eps in c(1, .7)) {
    if (ef == "between" && eps != 1) next
    e2 <- c(within = "intra", between = "entre", interaction = "interacao")[[ef]]
    x <- quieto(pwrss::power.f.mixed.anova(eta.squared = efeito_eta2_de_f(.25), factor.levels = c(3, 4), rho.within = rho,
                                           epsilon = eps, n.total = 60, effect = ef, verbose = 0))
    y <- as.numeric(campo(quieto(poder_anova_mr(f = .25, medidas = 4, grupos = 3, efeito = e2, rho = rho, epsilon = eps, n = 60)), "Actual power"))
    ok <- ok && abs(x$power - y) < 1e-4
  }
  ok })
teste("medidas repetidas: N igual ao pwrss::power.f.mixed.anova (ceil.n = TRUE)", {
  x <- quieto(pwrss::power.f.mixed.anova(eta.squared = efeito_eta2_de_f(.204), factor.levels = c(3, 2), rho.within = NA,
                                         power = .9, effect = "interaction", verbose = 0))
  y <- quieto(pwrss::power.f.mixed.anova(eta.squared = efeito_eta2_de_f(.25), factor.levels = c(2, 3), rho.within = .5,
                                         power = .8, effect = "within", verbose = 0))
  x$n.total == 309 && y$n.total == quieto(poder_anova_mr(f = .25, medidas = 3, grupos = 2, efeito = "intra", rho = .5))$n })
teste("regressão global: N = teto(v) + u + 1 do pwr::pwr.f2.test", {
  ok <- TRUE
  for (k in c(2, 5, 8)) for (f2 in c(.02, .15, .35)) ok <- ok && poder_regressao("global", k, f2 = f2)$n == ceiling(pwr::pwr.f2.test(u = k, f2 = f2, power = .8)$v) + k + 1
  ok })
teste("regressão ΔR²: N igual ao pwrss::power.f.regression (f² convertido para o argumento dele)", {
  ok <- TRUE
  for (f2 in c(.0667, .15)) for (q in c(1, 3)) {
    x <- pwrss::power.f.regression(r.squared.change = f2 / (1 + f2), k.total = 5, k.tested = q, power = .8, verbose = 0)
    ok <- ok && poder_regressao("incremento", 5, f2 = f2, testados = q)$n == x$n
  }
  ok })
teste("regressão ΔR²: o atalho do pwr.f2.test diverge (f² = .15, 2 de 5: G*Power 68, pwr 71)", {
  v <- pwr::pwr.f2.test(u = 2, f2 = .15, power = .8)$v
  poder_regressao("incremento", 5, f2 = .15, testados = 2)$n == 68 && ceiling(v) + 5 + 1 == 71 })
teste("qui-quadrado: N = teto do pwr::pwr.chisq.test (w = .1, .3, .5; gl 1 a 4)", {
  ok <- TRUE
  for (w in c(.1, .3, .5)) for (g in 1:4) ok <- ok && poder_quiquadrado(w, g)$n == ceiling(pwr::pwr.chisq.test(w = w, df = g, power = .8)$N)
  ok && poder_quiquadrado(.3, 1)$n == 88 })
teste("correlação exata: poder igual ao pwrss::power.exact.onecor com N dado (tol. 1e-4)", {
  x <- pwrss::power.exact.onecor(rho = .25, n = 100, verbose = 0)$power
  abs(as.numeric(campo(quieto(poder_correlacao(.25, n = 100)), "Actual power")) - x) < 1e-4 })

cat("== (4) Conversões de efeito\n")
teste("η² .01, .06, .14 -> f .1005, .2526, .4035 (tol. 5e-5)",
      all(abs(efeito_f_de_eta2(c(.01, .06, .14)) - c(.1005, .2526, .4035)) < 5e-5))
teste("f <-> η² ida e volta", all(abs(efeito_eta2_de_f(efeito_f_de_eta2(c(.02, .1, .3))) - c(.02, .1, .3)) < 1e-12))
teste("f de médias com n desiguais pondera pelas frequências", {
  f1 <- efeito_f_de_medias(c(10, 12, 15), 5, n = c(20, 20, 20)); f2 <- efeito_f_de_medias(c(10, 12, 15), 5, n = c(10, 20, 60))
  abs(f1 - efeito_f_de_medias(c(10, 12, 15), 5)) < 1e-12 && f2 < f1 })
teste("d de médias: DP pelo 'Determine' (média das variâncias) e ponderado com n", {
  abs(efeito_d_de_medias(10, 13, 6, 4) - 3 / sqrt(26)) < 1e-12 &&
    abs(efeito_d_de_medias(10, 13, 6, 4, 30, 10) - 3 / sqrt((29 * 36 + 9 * 16) / 38)) < 1e-12 })
teste("r <-> d ida e volta (grupos iguais) e d = .5 -> r = .2425", {
  abs(efeito_r_de_d(efeito_d_de_r(.3)) - .3) < 1e-12 && abs(efeito_r_de_d(.5) - .2425) < 5e-5 })
teste("f² de R² e de ΔR²: .20 -> .25; ΔR² .05 sobre .20 -> .0667", {
  abs(efeito_f2_de_r2(.2) - .25) < 1e-12 && abs(efeito_f2_de_delta_r2(.2, .05) - .05 / .75) < 1e-12 })
teste("OR de probabilidades e dz de d", abs(efeito_or_de_probabilidades(.6, .5) - 1.5) < 1e-12 &&
        abs(efeito_dz_de_d(.3, .5) - .3) < 1e-12 && abs(efeito_dz_de_d(.3, .3) - .2535) < 5e-5)
teste("d de t e n (rota B do curso): t = 2.5 com 30 e 30 -> .6455", abs(efeito_d_de_t(2.5, 30, 30) - .6455) < 5e-5)

cat("== (5) Comportamento das funções\n")
teste("perda N / (1 - perda): 62 com 20% -> 78; one-way 186 com 20% -> 234 (múltiplo de 3)", {
  r <- quieto(poder_anova(f = .23, niveis = 3, perda = .2))
  ajustar_perda(62, .2)$n_coleta == 78 && campo(r, "N para coleta (N / (1 - perda))") == "234" &&
    any(grepl("não vem das fontes do usuário", r$avisos)) && grepl("meta de coleta é de 234", r$apa) })
teste("post hoc devolve poder e avisa sobre efeito observado", {
  r <- quieto(poder_teste_t(.5, n = 128)); abs(as.numeric(campo(r, "Actual power")) - .8015) < .0005 &&
    any(grepl("post hoc", r$avisos)) && !length(r$apa) })
teste("frase APA com placeholder quando falta a origem do efeito", grepl("[origem do tamanho de efeito: preencher]",
                                                                        poder_anova(f = .25, niveis = 3)$apa, fixed = TRUE))
teste("frase APA em inglês e com vírgula decimal", {
  a <- poder_regressao("global", 5, r2 = .2, decimal = ",", idioma = "en")$apa
  grepl("Sample size was determined", a) && grepl("f² = 0,250", a, fixed = TRUE) && grepl("α = ,05", a, fixed = TRUE) })
teste("configuração do G*Power impressa (família, teste, campos)", {
  r <- poder_anova(eta2p = .06, niveis = c(3, 5), covariaveis = 2)
  campo(r, "Test family") == "F tests" && campo(r, "Statistical test") == "ANCOVA: Fixed effects, main effects and interactions" &&
    campo(r, "Numerator df") == "8" && campo(r, "Number of groups") == "15" && campo(r, "Number of covariates") == "2" })
teste("MR com rho = NA avisa que 'as in Cohen (1988)' não está implementada",
      any(grepl("as in Cohen \\(1988\\)", quieto(poder_anova_mr(f = .3, medidas = 3, rho = NA))$avisos)))
teste("objeto é pr_resultado (relatar() funciona) e o print não falha", {
  r <- poder_correlacao(.3); inherits(r, "pr_resultado") && identical(relatar(r), r$apa) &&
    length(utils::capture.output(print(r))) > 10 })
teste("erros claros: f e eta2p juntos; ε fora do limite; testados > preditores", {
  e1 <- tryCatch({ poder_anova(f = .2, eta2p = .05, niveis = 3); FALSE }, error = function(e) TRUE)
  e2 <- tryCatch({ poder_anova_mr(f = .2, medidas = 4, epsilon = .2); FALSE }, error = function(e) TRUE)
  e3 <- tryCatch({ poder_regressao("incremento", 3, f2 = .1, testados = 4); FALSE }, error = function(e) TRUE)
  e1 && e2 && e3 })

cat("== (6) Correções da revisão v1.1.0 (C4, C5, C6, C14a)\n")
# C4. Logística: exemplos do manual do G*Power 3.1 (versão de 1/6/2023), seção 29.4.
# OR = 1,5 (p1 = .5, p2 = .6), preditor normal, bicaudal, α = .05, poder .95.
teste("C4 logística, manual 29.4 com R² = 0: Demidenko com correção 337, sem correção 355, Hsieh 317", {
  poder_logistica(1.5, .5, poder = .95)$n == 337 && poder_logistica(1.5, .5, poder = .95, metodo = "demidenko")$n == 355 &&
    poder_logistica(1.5, .5, poder = .95, metodo = "hsieh")$n == 317 })
# Mesma seção, preditor binário: p1 = .05, p2 = .10 (OR = 2,111), poder .95; π = .5, .2 e
# .8 com correção: 1437, 2158, 2368; π = .5 sem correção: 1498.
teste("C4 logística binária, manual 29.4: 1437, 2158, 2368 com correção e 1498 sem", {
  or_b <- efeito_or_de_probabilidades(.10, .05)
  ns <- vapply(c(.5, .2, .8), function(pi) poder_logistica(or_b, .05, distribuicao = "binomial", prop_x = pi, poder = .95)$n, 0)
  all(ns == c(1437, 2158, 2368)) &&
    poder_logistica(or_b, .05, distribuicao = "binomial", prop_x = .5, poder = .95, metodo = "demidenko")$n == 1498 })
# Manual 29.4: "the necessary sample size increases from 337 to a value of 395 when using
# the procedure of Demidenko (2007) with variance correction" (R² other X = .1). Conferido
# à parte por integração numérica de v1 (X ~ N(0, 1)): N sem correção = 354,74 / 0,9 =
# 394,15 -> 395; a fórmula documentada (eq. 14-15, correção mantida) dá 336,86 / 0,9 = 374,29 -> 375.
teste("C4 logística, manual 29.4 com R² other X = .1 -> 395 em 'G*Power, saída'; 375 da fórmula documentada no bloco da skill e aviso", {
  r <- quieto_log(1.5, .5, r2_outros = .1, poder = .95)
  sai <- r$tabela[r$tabela$secao == "G*Power: saída", ]
  r$n == 395 && sai$valor[sai$campo == "Total sample size"] == "395" &&
    campo(r, "N pela fórmula documentada (correção mantida)") == "375" &&
    any(grepl("inferida de um só exemplo", r$avisos)) && any(grepl("395", r$decisoes)) })
teste("C4 logística com R² > 0 e preditor binário (sem exemplo publicado): N fora de 'G*Power, saída', aviso e frase sem G*Power", {
  r <- quieto_log(efeito_or_de_probabilidades(.10, .05), .05, distribuicao = "binomial", prop_x = .5, poder = .95, r2_outros = .1)
  sai <- r$tabela$campo[r$tabela$secao == "G*Power: saída"]
  !("Total sample size" %in% sai) && !("Actual power" %in% sai) && campo(r, "N mínimo") == as.character(r$n) &&
    any(grepl("não foi conferido contra o G\\*Power", r$avisos)) && !grepl("G*Power", r$apa, fixed = TRUE) &&
    grepl("Demidenko (2007)", r$apa, fixed = TRUE) })
teste("C4 frase da logística com R² > 0 traz o R² dos demais preditores", {
  grepl("com R² = .10 do preditor com os demais preditores", quieto_log(1.5, .5, r2_outros = .1, poder = .95)$apa, fixed = TRUE) })

# C5. Regra de alocação do G*Power (n1 = teto(n1*), n2 = teto(razão × n1*)).
# Manual seção 17 (diferença entre inclinações, teste t; dados de Dupont e Plummer):
# |Δ slope| = .01592, σ = .5578413, σx1 = 9.02914, σx2 = 11.86779, razão 1,571428, α = .05,
# poder .80 -> 163 + 256, δ = 2.811598, gl 415, poder .800980. δ = Δ·√n2, Δ = |Δ slope| /
# σR, σR = σ·√(1/(m·σx1²) + 1/σx2²), m = n1/n2 (eq. 6 da seção 17.5; m = 1/razão).
teste("C5 alocação, manual seção 17 (t de inclinações, razão 1,571428): 163 + 256, δ = 2.811598, poder .800980", {
  rz <- 1.571428
  pw_incl <- function(n1, n2 = rz * n1) {
    sR <- sqrt(.5578413^2 * (1 / ((1 / rz) * 9.02914^2) + 1 / 11.86779^2)); del <- .01592 / sR * sqrt(n2); gl <- n1 + n2 - 4
    cr <- stats::qt(.975, gl); c(poder = stats::pt(cr, gl, del, lower.tail = FALSE) + stats::pt(-cr, gl, del), delta = del) }
  a <- .ee_poder_alocacao(function(x) pw_incl(x)[["poder"]], rz, .80, lo = 5)
  v <- pw_incl(a$n1, a$n2)
  a$n1 == 163 && a$n2 == 256 && abs(v[["delta"]] - 2.811598) < 5e-6 && abs(v[["poder"]] - .800980) < 5e-6 })
# Manual seção 15 (duas variâncias, F): razão de variâncias 1,5, bicaudal, α = .05, poder
# .80; com alocação 1: 193 + 193 (poder .800105); com alocação 2: 148 + 295.
teste("C5 alocação, manual seção 15 (duas variâncias, alocação 2): 148 + 295; com alocação 1, poder .800105 em 193 + 193", {
  pw_var <- function(n1, n2, rv = 1.5) { lo <- stats::qf(.025, n2 - 1, n1 - 1); hi <- stats::qf(.975, n2 - 1, n1 - 1)
    stats::pf(lo / rv, n2 - 1, n1 - 1) + stats::pf(hi / rv, n2 - 1, n1 - 1, lower.tail = FALSE) }
  a <- .ee_poder_alocacao(function(x) pw_var(x, 2 * x), 2, .80, lo = 5)
  a$n1 == 148 && a$n2 == 295 && abs(pw_var(193, 193) - .800105) < 5e-7 })
# t de duas médias, d = .5, razão 1,571428 (sem exemplo no manual): n1* = 52,13 pela raiz de
# pwr::pwr.t2n.test com n2 contínuo -> 53 + 82 = 135. A regra antiga dava 53 + 84 = 137.
teste("C5 t independente, d = .5, razão 1,571428 -> 53 + 82 = 135 (n1* da raiz do pwr::pwr.t2n.test)", {
  r <- quieto(poder_teste_t(.5, razao = 1.571428))
  n1c <- stats::uniroot(function(x) pwr::pwr.t2n.test(n1 = x, n2 = 1.571428 * x, d = .5)$power - .80, c(5, 1e4), tol = 1e-10)$root
  r$n == 135 && campo(r, "Sample size group 1") == "53" && campo(r, "Sample size group 2") == "82" &&
    ceiling(n1c) == 53 && ceiling(1.571428 * n1c) == 82 && grepl("(53 + 82)", r$apa, fixed = TRUE) })

# C6. α na frase com as casas do valor usado (antes: .00, .02 e .01). N do t e da ANOVA
# conferidos no pwr (pwr.t.test e pwr.anova.test com o mesmo α), para garantir que o α
# impresso é o da conta.
teste("C6 frase com α = .001, .0167 (.05/3) e .005, e N do mesmo α", {
  a <- poder_teste_t(.5, alfa = .001); b <- poder_anova(.25, niveis = 3, alfa = .05 / 3); c3 <- poder_correlacao(.3, alfa = .005)
  grepl("α = .001 (bicaudal)", a$apa, fixed = TRUE) && a$n == 2 * ceiling(pwr::pwr.t.test(d = .5, sig.level = .001, power = .8)$n) &&
    grepl("α = .0167 e", b$apa, fixed = TRUE) && b$n == 3 * ceiling(pwr::pwr.anova.test(k = 3, f = .25, sig.level = .05 / 3, power = .8)$n) &&
    grepl("α = .005 (bicaudal)", c3$apa, fixed = TRUE) })
teste("C6 α comum continua com 2 casas (.05, .01, .10) e vírgula decimal funciona (,001)", {
  grepl("α = .05 (bicaudal)", poder_teste_t(.5)$apa, fixed = TRUE) &&
    grepl("α = .01 (bicaudal)", poder_teste_t(.5, alfa = .01)$apa, fixed = TRUE) &&
    .ee_poder_alfa_txt(.10) == ".10" && grepl("α = ,001", poder_teste_t(.5, alfa = .001, decimal = ",")$apa, fixed = TRUE) })

# C14a. Correlação com metodo = "fisher": 'G*Power, saída' só com o que o G*Power devolve
# (exato: 61, livro) e frase sem G*Power; o 62 de Fisher fica no bloco da skill.
teste("C14a correlação ρ = .35 com metodo = 'fisher': saída do G*Power = 61, skill = 62, frase cita Fisher e não o G*Power", {
  r <- poder_correlacao(.35, metodo = "fisher")
  sai <- r$tabela[r$tabela$secao == "G*Power: saída", ]
  r$n == 62 && sai$valor[sai$campo == "Total sample size"] == "61" && campo(r, "N pela aproximação de Fisher") == "62" &&
    !grepl("G*Power", r$apa, fixed = TRUE) && grepl("aproximação z de Fisher", r$apa, fixed = TRUE) &&
    grepl("G*Power 3.1", poder_correlacao(.35)$apa, fixed = TRUE) })

cat("== (7) Correções da verificação final (G-12, G-02, F1, F2, G-08, G-09, H-15)\n")
# G-12. O efeito impresso na frase precisa devolver o N relatado. Antes, a conferência
# usava round() do R (,625 -> ,62) e a frase .pr_num (,625 -> ,63): "d = 0.63 ... 84",
# mas com d = 0,63 o N é 82. Conferência independente com o pwr: N do valor impresso.
teste("G-12 d = .625 -> frase 'd = 0.625' e 84; f = .345 -> 'f = 0.345' e 87 (pwr: 84 e 87 com o impresso; 82 e 84 com .63 e .35)", {
  a <- poder_teste_t(.625); b <- poder_anova(f = .345, niveis = 3)
  a$n == 84 && grepl("d = 0.625 ", a$apa, fixed = TRUE) && 2 * ceiling(pwr::pwr.t.test(d = .625, power = .8)$n) == 84 &&
    2 * ceiling(pwr::pwr.t.test(d = .63, power = .8)$n) == 82 &&
    b$n == 87 && grepl("f = 0.345 ", b$apa, fixed = TRUE) && 3 * ceiling(pwr::pwr.anova.test(k = 3, f = .345, power = .8)$n) == 87 &&
    3 * ceiling(pwr::pwr.anova.test(k = 3, f = .35, power = .8)$n) == 84 })
# Varredura da revisão (cove_final/G/g_rep_G12_casas_poder.R, 68 divergências antes): efeitos
# de ,105 a 1,495 (passo ,01), poder .80, .90 e .95, t independente, t de uma amostra e
# ANOVA de 3 grupos; o valor lido da frase vai para o pwr e o N tem de ser o mesmo.
teste("G-12 varredura: 1260 combinações, o efeito lido da frase devolve no pwr o N da frase", {
  ok <- TRUE
  le <- function(txt, s) as.numeric(sub(paste0(".*com ", s, " = ([0-9.]+) .*"), "\\1", txt))
  for (pw in c(.8, .9, .95)) for (x in round(seq(.105, 1.495, by = .01), 3)) {
    a <- poder_teste_t(x, poder = pw); b <- poder_teste_t(x, "uma_amostra", poder = pw); f <- poder_anova(f = x, niveis = 3, poder = pw)
    ok <- ok && a$n == 2 * ceiling(pwr::pwr.t.test(d = le(a$apa, "d"), power = pw)$n) &&
      b$n == ceiling(pwr::pwr.t.test(d = le(b$apa, "d"), power = pw, type = "one.sample")$n) &&
      f$n == 3 * ceiling(pwr::pwr.anova.test(k = 3, f = le(f$apa, "f"), power = pw)$n)
  }
  ok })
# Mesma classe na regressão (f² com 3 casas no mínimo): f² de ,0105 a ,3505 (passo ,005),
# 2 e 5 preditores, poder .80 e .90; antes, 19 de 276 divergiam (ex.: f² = .1205, 5
# preditores, poder .90: frase com 0.121 e 143, que no pwr dá outro N). Conferência:
# N = teto(v) + u + 1 do pwr::pwr.f2.test com o f² lido da frase.
teste("G-12 regressão global: 276 combinações, o f² lido da frase devolve no pwr o N da frase", {
  ok <- TRUE
  for (pw in c(.8, .9)) for (k in c(2, 5)) for (x in round(seq(.0105, .3505, by = .005), 4)) {
    a <- poder_regressao("global", k, f2 = x, poder = pw); f <- as.numeric(sub(".*f² = ([0-9.]+) .*", "\\1", a$apa))
    ok <- ok && a$n == ceiling(pwr::pwr.f2.test(u = k, f2 = f, power = pw)$v) + k + 1
  }
  ok })

# G-02. Hsieh com preditor binário parava com erro (o pwrss só aceita "binomial" com
# size = 1 no Hsieh). Manual do G*Power, seção 29.4: p1 = .05, p2 = .10, π = .5, poder .95,
# N = 1437 ("replicates the sample size given in Table I in Hsieh et al."), Actual power
# 0.950068. Conta independente: fórmula de Hsieh et al. (1998) para X binário escrita aqui,
# n = (zα·√(p̄(1 - p̄)/B) + zβ·√(p1(1 - p1) + p2(1 - p2)(1 - B)/B))² / ((p1 - p2)²(1 - B)) = 1436,48.
# Preditor normal, mesma seção: Hsieh 317 com Actual power 0.950486 (antes saía o poder pedido, .95).
teste("G-02 logística Hsieh, preditor binário: 1437 e poder 0.950068 (manual 29.4); normal: 317 e 0.950486", {
  r <- poder_logistica(efeito_or_de_probabilidades(.10, .05), .05, distribuicao = "binomial", prop_x = .5, metodo = "hsieh", poder = .95)
  B <- .5; p1 <- .05; p2 <- .10; pb <- (1 - B) * p1 + B * p2
  n_h <- (qnorm(.975) * sqrt(pb * (1 - pb) / B) + qnorm(.95) * sqrt(p1 * (1 - p1) + p2 * (1 - p2) * (1 - B) / B))^2 / ((p1 - p2)^2 * (1 - B))
  rn <- poder_logistica(1.5, .5, poder = .95, metodo = "hsieh")
  r$n == 1437 && ceiling(n_h) == 1437 && abs(n_h - 1436.48) < .005 && campo(r, "Actual power") == "0.950068" &&
    rn$n == 317 && campo(rn, "Actual power") == "0.950486" })

# F1. Com d grande e razão diferente de 1 o uniroot parava ("f() values at end points not of
# opposite sign"); a versão anterior dava 2 + 6 = 8 e 2 + 4 = 6. Conferência: poder do
# desenho mínimo pelo pwr::pwr.t2n.test (acima de .80) e, nos casos em que o mínimo não
# basta, n1* pela raiz do pwr::pwr.t2n.test com n2 = razão × n1 contínuo.
teste("F1 d = 3, razão 3 -> 2 + 6 = 8; d = 4, razão 2 -> 2 + 4; d = 4, razão .5 -> 4 + 2 (desenho mínimo, poder do pwr >= .80, com aviso)", {
  a <- poder_teste_t(3, razao = 3); b <- poder_teste_t(4, razao = 2); c5 <- poder_teste_t(4, razao = .5)
  a$n == 8 && campo(a, "Sample size group 1") == "2" && campo(a, "Sample size group 2") == "6" &&
    b$n == 6 && campo(b, "Sample size group 2") == "4" && c5$n == 6 && campo(c5, "Sample size group 1") == "4" &&
    pwr::pwr.t2n.test(n1 = 2, n2 = 6, d = 3)$power >= .80 && pwr::pwr.t2n.test(n1 = 2, n2 = 4, d = 4)$power >= .80 &&
    pwr::pwr.t2n.test(n1 = 4, n2 = 2, d = 4)$power >= .80 && any(grepl("menor desenho", a$avisos)) })
teste("F1 efeitos grandes com razão != 1 em que o mínimo não basta: d = 2/razão 4, 2,5/3 (poder .90), 1,8/5 = raiz do pwr::pwr.t2n.test", {
  ok <- TRUE
  for (cfg in list(c(2, 4, .8), c(2.5, 3, .9), c(1.8, 5, .8))) {
    r <- poder_teste_t(cfg[1], razao = cfg[2], poder = cfg[3])
    n1c <- stats::uniroot(function(x) pwr::pwr.t2n.test(n1 = x, n2 = cfg[2] * x, d = cfg[1])$power - cfg[3], c(2, 1e4), tol = 1e-10)$root
    ok <- ok && campo(r, "Sample size group 1") == as.character(ceiling(n1c)) && campo(r, "Sample size group 2") == as.character(ceiling(cfg[2] * n1c)) &&
      !any(grepl("menor desenho", r$avisos))
  }
  ok && poder_teste_t(2, razao = 4)$n == 15 })

# F2. Poder na frase com as casas do valor usado, como o α (antes .975 -> ".98" e .999 ->
# "1.00"). N conferido no pwr::pwr.t.test com o mesmo poder: 248 e 412.
teste("F2 poder .975 -> 'poder de .975' (248) e .999 -> 'power of .999' (412); .80 e .90 continuam com 2 casas", {
  a <- poder_teste_t(.5, poder = .975); b <- poder_teste_t(.5, poder = .999, idioma = "en")
  grepl("poder de .975,", a$apa, fixed = TRUE) && a$n == 2 * ceiling(pwr::pwr.t.test(d = .5, power = .975)$n) && a$n == 248 &&
    grepl("power of .999,", b$apa, fixed = TRUE) && b$n == 2 * ceiling(pwr::pwr.t.test(d = .5, power = .999)$n) && b$n == 412 &&
    grepl("poder de .90,", poder_teste_t(.5, poder = .9)$apa, fixed = TRUE) && grepl("poder de ,975", poder_teste_t(.5, poder = .975, decimal = ",")$apa, fixed = TRUE) })

# G-08. Bloco "G*Power, saída" com 6 casas, como o manual. Seção 11.3 (3 × 3 × 4, efeito A,
# f = .7066856, N = 108): λ 53.935690, F crítico 3.123907, gl 2 e 72, "Power (1-β err prob):
# 0.99999" (antes a skill mostrava "1"). O R dá 0,9999994; a tolerância de 1e-5 cobre a
# diferença de cálculo da F não central na sexta casa.
teste("G-08 manual 11.3: λ = 53.935690, F crítico = 3.123907, gl 72 e poder com 6 casas (0.99999, tol. 1e-5), não '1'", {
  r <- poder_anova(f = .7066856, niveis = c(3, 3, 4), efeito = "A", n = 108); pw <- campo(r, "Actual power")
  campo(r, "Noncentrality parameter λ") == "53.935690" && campo(r, "Critical F") == "3.123907" && campo(r, "Denominator df") == "72" &&
    pw != "1" && grepl("^0\\.[0-9]{6}$", pw) && abs(as.numeric(pw) - .99999) < 1e-5 })

# G-09. t de uma amostra é "Effect size d" no G*Power (manual, seção 20.3); dz é do pareado.
# Exemplos da seção: d = .625 unicaudal, poder .95 -> δ 3.423266, t crítico 1.699127, gl 29,
# N = 30, poder 0.955144; d = .1, α = .01, poder .90 -> δ 3.862642, t crítico 2.579131, N = 1492, poder 0.900169.
teste("G-09 manual 20.3: uma amostra com 'Effect size d' e 'd =' na frase; 30 (δ 3.423266, t 1.699127, 0.955144) e 1492 (0.900169)", {
  a <- poder_teste_t(.625, "uma_amostra", caudas = 1, poder = .95); b <- poder_teste_t(.1, "uma_amostra", alfa = .01, poder = .90)
  pa <- poder_teste_t(.3, "pareado")
  a$n == 30 && "Effect size d" %in% a$tabela$campo && !("Effect size dz" %in% a$tabela$campo) &&
    campo(a, "Noncentrality parameter δ") == "3.423266" && campo(a, "Critical t") == "1.699127" && campo(a, "Df") == "29" &&
    campo(a, "Actual power") == "0.955144" && grepl("com d = 0.625 ", a$apa, fixed = TRUE) && a$efeito$indice[1] == "d" &&
    b$n == 1492 && campo(b, "Noncentrality parameter δ") == "3.862642" && campo(b, "Critical t") == "2.579131" &&
    campo(b, "Actual power") == "0.900169" && "Effect size dz" %in% pa$tabela$campo && grepl("com dz = ", pa$apa, fixed = TRUE) })

# H-15. Faul et al. (2009, p. 1158): análises de compromisso (q = β/α = 1) com N = 1386,
# R² other X = .40 e Demidenko sem correção. Ex. 1: X binário, π = .1, p0 = .30, OR = 5,44444
# (o valor que o G*Power transfere), z crítico 3,454879 e α = β = .00055. Ex. 2: X normal,
# p0 = .10, p1 = .15, OR = (.15/.85)(.90/.10) = 1,588235 (o texto arredonda para 1,588),
# z crítico 2,146971 e α = β = .031796. A skill não tem modo de compromisso: o teste acha o
# α com 1 - poder = α usando o poder post hoc de poder_logistica() (motor pwrss). Tolerância
# de 2e-6 no z: a integração numérica do pwrss difere do G*Power na sexta casa.
teste("H-15 Faul et al. (2009), compromisso com R² other X = .40: z crítico 3.454879 (binário) e 2.146971 (normal), tol. 2e-6", {
  comp <- function(or, p0, dist, pi = .5) {
    g <- function(a) { r <- quieto_log(or, p0, r2_outros = .40, distribuicao = dist, prop_x = pi, metodo = "demidenko", n = 1386, alfa = a)
      (1 - r$poder_real) - a }
    a <- stats::uniroot(g, c(1e-8, .3), tol = 1e-14)$root; c(alfa = a, z = stats::qnorm(1 - a / 2))
  }
  e1 <- comp(5.44444, .30, "binomial", .1); e2 <- comp(efeito_or_de_probabilidades(.15, .10), .10, "normal")
  abs(e1[["z"]] - 3.454879) < 2e-6 && abs(e1[["alfa"]] - .00055) < 5e-6 &&
    abs(e2[["z"]] - 2.146971) < 2e-6 && abs(e2[["alfa"]] - .031796) < 5e-7 })

cat("== v1.4: AFC e MEE (Kelley & Lai, 2018, cap. 5 do Wiley Handbook; exemplos do capítulo)\n")
teste("RMSEA, ajuste exato (ε0 = 0, ε = .045, gl = 5): N = 1.268 (capítulo: 1.269 pelo MBESS, que para a .001 do alvo)", {
  r <- poder_afc_rmsea(0, .045, 5); r$n == 1268 && r$poder_alcancado >= .80 &&
    .ee_rmsea_pw(1267, 0, .045, 5, .05)$poder < .80 })
teste("RMSEA, ajuste não aproximado (H0 ε ≥ .08, ε = .045, gl = 5): N = 1.116 (capítulo: 1.113 pelo MBESS)", {
  r <- poder_afc_rmsea(.08, .045, 5); r$n == 1116 && .ee_rmsea_pw(1115, .08, .045, 5, .05)$poder < .80 })
teste("RMSEA, ajuste aproximado de MacCallum et al. (1996): ε0 = .05 contra .08 com gl = 20 dá N = 434; λ = (N - 1) gl ε²", {
  r <- poder_afc_rmsea(.05, .08, 20)
  r$n == 434 && abs(as.numeric(r$tabela$valor[r$tabela$campo == "λ sob H1"]) - 433 * 20 * .08^2) < 1e-3 })
if (requireNamespace("semTools", quietly = TRUE)) {
  teste("igual a semTools::findRMSEAsamplesize em 6 cenários (exato, aproximado, não aproximado, 2 grupos)", {
    cs <- list(c(0, .045, 5, 1), c(.08, .045, 5, 1), c(.05, .08, 20, 1), c(.05, .01, 20, 1), c(.05, .10, 40, 1), c(.05, .08, 20, 2))
    all(vapply(cs, function(z) {
      st <- suppressMessages(semTools::findRMSEAsamplesize(z[1], z[2], z[3], .80, .05, group = z[4]))
      poder_afc_rmsea(z[1], z[2], z[3], grupos = z[4])$n == round(st * z[4]) }, TRUE)) })
} else cat("PULADO  comparação com semTools (não instalado)\n")
teste("AIPE do RMSEA (ε = .045, gl = 5, IC 90%): larguras .060, .050 e .040 dão N = 802, 1.059 e 1.544, como no capítulo", {
  all(vapply(c(.06, .05, .04), function(w) n_aipe_rmsea(.045, 5, w)$n, 0) == c(802, 1059, 1544)) })
teste("IC do RMSEA com N = 802: largura <= .060 e com N = 801 maior (busca do menor N)", {
  w <- function(N) diff(.ee_ic_rmsea(.045, 5, N, .90)[c("inf", "sup")])
  w(802) <= .06 && w(801) > .06 })
teste("aninhados por RMSEA (MacCallum, Browne & Cai, 2006): restrito .0625 -> 336, .07 -> 268 (capítulo) e .06 -> 365 (capítulo: 364 pelo MBESS)", {
  f <- function(e) poder_afc_aninhados(6, 5, rmsea_restrito = e)$n
  f(.0625) == 336 && f(.07) == 268 && f(.06) == 365 })
teste("aninhados por F (Satorra & Saris, 1985): F_restrito = .02346, Δgl = 1 -> N = 336, como no capítulo",
      poder_afc_aninhados(6, 5, F_restrito = .02346)$n == 336)
# Σ1 da tabela 5.4 do capítulo: correlações (triângulo superior) e DP 1,20; 1,20; 0,70; 1,30; 1,10.
s1_dp <- c(1.2, 1.2, .7, 1.3, 1.1); s1_r <- c(.418, .394, .435, .348, .490, .540, .432, .510, .408, .450)
S1 <- matrix(0, 5, 5); S1[lower.tri(S1)] <- s1_r; S1 <- S1 + t(S1); diag(S1) <- 1
S1 <- diag(s1_dp) %*% S1 %*% diag(s1_dp); dimnames(S1) <- list(paste0("X", 1:5), paste0("X", 1:5))
teste("Σ1 montada confere com as covariâncias da tabela 5.4 (.601, .842, .644; tol. .002)",
      abs(S1[2, 1] - .601) < .002 && abs(S1[4, 2] - .842) < .002 && abs(S1[5, 4] - .644) < .002)
m_completo <- "F1 =~ NA*X1 + l1*X1 + X2 + X3 + X4 + X5\nF1 ~~ 1*F1"
teste("AIPE de λ1 (largura .15, IC 95%) com Σ1: N = 1.007, como o MBESS::ss.aipe.sem.path no capítulo",
      n_aipe_parametro(m_completo, S1, "l1", .15)$n == 1007)
teste("aninhados com Σ e sintaxes: Σ implicada pelo próprio modelo restrito dá ΔF ~ 0 e para com mensagem", {
  Sr <- lavaan::fitted(lavaan::cfa("F1 =~ NA*X1 + a*X1 + a*X2 + X3 + X4 + X5\nF1 ~~ 1*F1", sample.cov = S1, sample.nobs = 1000, sample.cov.rescale = FALSE))$cov
  m <- tryCatch({ poder_afc_aninhados(Sigma = Sr, modelo_restrito = "F1 =~ NA*X1 + a*X1 + a*X2 + X3 + X4 + X5\nF1 ~~ 1*F1",
                                      modelo_completo = m_completo); "" }, error = function(e) conditionMessage(e))
  grepl("discrepância do restrito", m) })
teste("aninhados com Σ1: F de cada modelo = 2 x fmin do lavaan sem reescalar e N pela qui-quadrado não central (conta à mão)", {
  r <- poder_afc_aninhados(Sigma = S1, modelo_restrito = "F1 =~ NA*X1 + a*X1 + a*X2 + X3 + X4 + X5\nF1 ~~ 1*F1", modelo_completo = m_completo)
  Fd <- r$delta_F; pw <- function(N) stats::pchisq(stats::qchisq(.95, 1), 1, ncp = (N - 1) * Fd, lower.tail = FALSE)
  pw(r$n) >= .80 && pw(r$n - 1) < .80 })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
cat("R", paste(R.version$major, R.version$minor, sep = "."), "| pwrss", as.character(utils::packageVersion("pwrss")),
    "| pwr", as.character(utils::packageVersion("pwr")), "\n")
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
