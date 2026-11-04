# evals_processos.R, testes do módulo de moderação, mediação e processo condicional
#
# Rodar da pasta da skill:
#   Rscript scripts/evals_processos.R
# Os testes de âncora leem os bancos da biblioteca do JASP (Presumed Media Influence e
# Climate Change Disasters, pasta "20. Process") a partir de Sys.getenv("PR_DADOS_JASP")
# (a pasta "Data Library"); sem a variável, esses testes são pulados com aviso.
# Saída: uma linha PASSOU/FALHOU por teste e o resumo. Código de saída 1 se algum falhar.
# Tolerâncias: valor impresso no livro com 3 casas = metade da última casa (.0005);
# limites de IC bootstrap = erro de Monte Carlo declarado no teste.

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
perto <- function(obtido, esperado, tol) all(abs(obtido - esperado) <= tol)
t_inicio <- Sys.time()

lin <- function(res, caminho) res$tabela[res$tabela$caminho == caminho, ]
efe <- function(res, padrao) res$efeito[grepl(padrao, res$efeito$indice, fixed = TRUE), ]
ind <- function(res, rotulo) res$efeito$valor[res$efeito$indice == paste0("indireto: ", rotulo)]

# ---------------------------------------------------------------- âncoras JASP (livro, Vol. 6)
raiz <- Sys.getenv("PR_DADOS_JASP")
acha <- function(arquivo) {
  if (!nzchar(raiz)) return(NULL)
  f <- file.path(raiz, "20. Process", arquivo)
  if (file.exists(f)) return(f)
  f <- list.files(raiz, pattern = paste0("^", gsub(".", "\\.", arquivo, fixed = TRUE), "$"), recursive = TRUE, full.names = TRUE)
  if (length(f)) f[1] else NULL
}
f_pmi <- acha("Presumed Media Influence.csv"); f_cc <- acha("Climate Change Disasters.csv")
if (is.null(f_pmi) || is.null(f_cc)) {
  cat("AVISO   PR_DADOS_JASP vazia ou sem a pasta '20. Process': testes de âncora pulados\n")
} else {
  pmi <- read.csv(f_pmi, fileEncoding = "UTF-8-BOM")
  cc <- read.csv(f_cc, fileEncoding = "UTF-8-BOM")

  cat("== Âncora: mediação simples, Presumed Media Influence\n")
  med <- mediacao(pmi, "Condition", "MediaInfluence", "Intention")   # padrão: 10.000 reamostras, semente 2026
  teste("a = 0,477; b = 0,506; c = 0,496; c' = 0,254; ab = 0,241 (tol .0005)",
        perto(c(lin(med, "a")$b, lin(med, "b")$b, lin(med, "c")$b, lin(med, "c'")$b, efe(med, "indireto:")$valor),
              c(.477, .506, .496, .254, .241), .0005))
  teste("p de b < ,001; p de c = ,077; p de c' = ,322 (tol .0005)",
        lin(med, "b")$p < .001 && perto(lin(med, "c")$p, .077, .0005) && perto(lin(med, "c'")$p, .322, .0005))
  teste("sem arredondar: c' = 0,2544, ab = 0,2413, c = 0,4957 (tol .00005) e c = c' + ab exata",
        perto(c(lin(med, "c'")$b, efe(med, "indireto:")$valor, lin(med, "c")$b), c(.2544, .2413, .4957), .00005) &&
          abs(med$identidade[["diferenca"]]) < 1e-12)
  teste("Sobel z = 1,89 (tol .005), p = ,059 (tol .0005), só como comparação",
        perto(med$sobel$z, 1.89, .005) && perto(med$sobel$p, .059, .0005))
  ic <- unlist(efe(med, "indireto:")[c("ic_inf", "ic_sup")])
  cat(sprintf("        IC bootstrap obtido [%.4f; %.4f]; livro [0,006; 0,521] e [0,023; 0,541] em outra semente\n", ic[1], ic[2]))
  teste("IC 95% bootstrap percentílico, 10.000 reamostras, semente 2026: [0,006; 0,521] (tol Monte Carlo .01; zero fora)",
        perto(ic, c(.006, .521), .01) && ic[1] > 0)
  teste("frase APA com IC bootstrap percentílico e sem rótulo completa/parcial",
        any(grepl("bootstrap percentílico", med$apa)) && !any(grepl("complet|parcial", med$apa)) &&
          any(grepl("semente 2026", med$apa)))

  cat("== Âncora: moderação, Climate Change Disasters\n")
  mod <- moderacao(cc, "Justification", "Framing", "Skepticism")
  cf <- setNames(mod$tabela$b, mod$tabela$termo); pv <- setNames(mod$tabela$p, mod$tabela$termo)
  teste("b1 = -0,562; b2 = 0,105; b3 = 0,201 (tol .0005)",
        perto(unname(cf[2:4]), c(-.562, .105, .201), .0005))
  teste("p: b1 = ,011; b2 = ,006 (tol .0005); b3 < ,001",
        perto(unname(pv[2:3]), c(.011, .006), .0005) && pv[[4]] < .001)
  teste("R² = ,246; ΔR² = ,048; f² = ,064 (tol .0005)",
        perto(c(mod$resumo_modelo$r2, mod$teste_interacao$delta_r2, mod$teste_interacao$f2), c(.246, .048, .064), .0005))
  modc <- moderacao(cc, "Justification", "Framing", "Skepticism", centralizar = TRUE)
  teste("centralizado: b1 = 0,117; b3 = 0,201 e R² = ,246 inalterados; EP(b3) idêntico",
        perto(modc$tabela$b[2], .117, .0005) && abs(modc$tabela$b[4] - mod$tabela$b[4]) < 1e-12 &&
          abs(modc$tabela$ep[4] - mod$tabela$ep[4]) < 1e-12 && abs(modc$resumo_modelo$r2 - mod$resumo_modelo$r2) < 1e-12)
  teste("efeitos condicionais -1 DP, média, +1 DP: -0,292, 0,117, 0,526 (tol .0005)",
        perto(mod$posthoc$efeito, c(-.292, .117, .526), .0005))
  teste("p dos efeitos condicionais: ,068; ,297; ,001 (tol .0005)",
        perto(mod$posthoc$p, c(.068, .297, .001), .0005))
  teste("padrão descrito como inversão de sinal", identical(mod$padrao, "inverte"))
  teste("frase APA da moderação traz b3, ΔR² e f² do livro e as fronteiras de Johnson-Neyman",
        grepl("b = 0.201", mod$apa[1], fixed = TRUE) && grepl("ΔR² = .048", mod$apa[1], fixed = TRUE) &&
          grepl("f² = 0.064", mod$apa[1], fixed = TRUE) && length(mod$apa) == 3 && grepl("Johnson-Neyman", mod$apa[3]))

  cat("== Âncora do curso (mesmo banco, valores de Hayes reproduzidos no curso)\n")
  par <- mediacao(pmi, "Condition", c("Importance", "MediaInfluence"), "Intention", reamostras = 2000)
  con <- efe(par, "contraste")
  teste("paralela: via Importance 0,203; via MediaInfluence 0,189; total 0,392; contraste 0,014 (tol .0005)",
        perto(c(ind(par, "Condition → Importance → Intention"), ind(par, "Condition → MediaInfluence → Intention"),
                efe(par, "indireto total")$valor, con$valor),
              c(.203, .189, .392, .014), .0005))
  teste("paralela: IC do contraste inclui zero (curso 13c-5)", con$ic_inf < 0 && con$ic_sup > 0)
  ser <- mediacao(pmi, "Condition", c("Importance", "MediaInfluence"), "Intention", tipo = "serial", reamostras = 2000)
  mm14 <- mediacao_moderada(cc, "Framing", "DonationAttitude", "Justification", "Skepticism", modelo = 14, reamostras = 1000)
  teste("B11: índice de mediação moderada pequeno (-0,00049) sai com sinal e 2 algarismos, não '0.000'",
        grepl("moderada foi -0.00049 (", mm14$apa[2], fixed = TRUE) && !grepl("foi 0.000 ", mm14$apa[2], fixed = TRUE) &&
          any(grepl("algarismos significativos", mm14$decisoes)))
  teste("serial: efeito Importance → MediaInfluence = 0,049 (tol .0005); indireto total = 0,392 como na paralela",
        perto(ind(ser, "Condition → Importance → MediaInfluence → Intention"), .049, .0005) &&
          abs(efe(ser, "indireto total")$valor - efe(par, "indireto total")$valor) < 1e-12)
}

# ---------------------------------------------------------------- verdade conhecida
cat("== Verdade conhecida em dados simulados\n")
set.seed(101)
n <- 2000
s1 <- data.frame(x = rnorm(n)); s1$m <- .4 * s1$x + rnorm(n); s1$y <- .3 * s1$m + .2 * s1$x + rnorm(n)
r1 <- mediacao(s1, "x", "m", "y", reamostras = 2000)
i1 <- efe(r1, "indireto:")
teste("mediação simples: a = .4, b = .3 -> ab perto de .12 (tol .03) e IC cobre .12",
      perto(i1$valor, .12, .03) && i1$ic_inf < .12 && i1$ic_sup > .12)
set.seed(102)
s0 <- data.frame(x = rnorm(300)); s0$m <- rnorm(300); s0$y <- .5 * s0$m + .3 * s0$x + rnorm(300)
r0 <- mediacao(s0, "x", "m", "y", reamostras = 2000)
teste("mediação nula (a = 0): IC do indireto inclui zero", efe(r0, "indireto:")$ic_inf < 0 && efe(r0, "indireto:")$ic_sup > 0)

set.seed(103)
s2 <- data.frame(x = rnorm(n), w = rnorm(n)); s2$y <- 1 + .5 * s2$x + .3 * s2$w + .4 * s2$x * s2$w + rnorm(n)
r2 <- moderacao(s2, "y", "x", "w", valores_w = c(-1, 0, 1))
teste("moderação: b3 = .4 recuperado (tol .06); efeito em W = 1 perto de .9 (tol .08)",
      perto(r2$tabela$b[4], .4, .06) && perto(r2$posthoc$efeito[3], .9, .08))
jn <- r2$johnson_neyman
zero_theta <- -r2$tabela$b[2] / r2$tabela$b[4]
teste("Johnson-Neyman: |t| = t crítico em cada fronteira (tol 1e-8) e o zero de θ(W) = b1 + b3·W fica entre as duas fronteiras",
      !is.null(jn$pontos) && nrow(jn$pontos) == 2 && perto(abs(jn$pontos$t), jn$t_critico, 1e-8) &&
        zero_theta > jn$pontos$w[1] && zero_theta < jn$pontos$w[2])
g <- jn$grade
reg <- jn$regioes
na_reg <- vapply(g$w, function(v) { j <- which(v >= reg$de & v <= reg$ate)[1]; reg$significativo[j] }, logical(1))
fora_front <- !vapply(g$w, function(v) any(abs(v - jn$pontos$w) < 1e-9), logical(1))
teste("Johnson-Neyman: na grade, p < .05 exatamente nas regiões marcadas como significativas",
      all((g$p < .05)[fora_front] == na_reg[fora_front]))
teste("Johnson-Neyman: % de casos das regiões soma 100", perto(sum(reg$pct_casos), 100, 1e-8))

set.seed(104)
s3 <- data.frame(x = rnorm(3000), g = factor(sample(c("A", "B", "C"), 3000, TRUE)))
inc <- c(A = .2, B = .5, C = .8)
s3$y <- inc[as.character(s3$g)] * s3$x + .3 * (s3$g == "B") + rnorm(3000)
r3 <- moderacao(s3, "y", "x", "g")
teste("moderador categórico (3 níveis): inclinações .2/.5/.8 recuperadas (tol .08) e F com gl1 = 2",
      perto(r3$posthoc$efeito, c(.2, .5, .8), .08) && r3$teste_interacao$gl1 == 2)

set.seed(105)
s4 <- data.frame(x = rnorm(3000)); s4$m1 <- .5 * s4$x + rnorm(3000)
s4$m2 <- .2 * s4$x + .4 * s4$m1 + rnorm(3000); s4$y <- .2 * s4$x + .2 * s4$m1 + .3 * s4$m2 + rnorm(3000)
r4 <- mediacao(s4, "x", c("m1", "m2"), "y", tipo = "serial", reamostras = 1000)
teste("serial: a1·d21·b2 = .5 × .4 × .3 = .06 recuperado (tol .02); três caminhos indiretos",
      perto(ind(r4, "x → m1 → m2 → y"), .06, .02) && sum(grepl("^indireto: ", r4$efeito$indice)) == 3)

set.seed(106)
s5 <- data.frame(x = rnorm(3000), w = rnorm(3000))
s5$m <- .3 * s5$x + .1 * s5$w + .2 * s5$x * s5$w + rnorm(3000); s5$y <- .5 * s5$m + .1 * s5$x + rnorm(3000)
r5 <- mediacao_moderada(s5, "x", "m", "y", "w", modelo = 7, valores_w = c(-1, 0, 1), reamostras = 1000)
teste("modelo 7: índice a3·b = .2 × .5 = .1 (tol .03); indireto em W = 1 perto de .25 (tol .04)",
      perto(r5$efeito$valor[1], .1, .03) && perto(r5$efeito$valor[4], .25, .04))
s5$y2 <- .2 * s5$x + .3 * s5$m + .1 * s5$w + .2 * s5$m * s5$w + rnorm(3000)
s5$m2 <- .5 * s5$x + rnorm(3000); s5$y3 <- .2 * s5$x + .3 * s5$m2 + .1 * s5$w + .2 * s5$m2 * s5$w + rnorm(3000)
r6 <- mediacao_moderada(s5, "x", "m2", "y3", "w", modelo = 14, reamostras = 1000)
teste("modelo 14: índice a·b3 = .5 × .2 = .1 (tol .03)", perto(r6$efeito$valor[1], .1, .03))
set.seed(108)
s15 <- data.frame(x = rnorm(4000), w = rnorm(4000)); s15$m <- .5 * s15$x + rnorm(4000)
s15$y <- .2 * s15$x + .3 * s15$m + .1 * s15$w + .25 * s15$x * s15$w + .2 * s15$m * s15$w + rnorm(4000)
r7 <- mediacao_moderada(s15, "x", "m", "y", "w", modelo = 15, valores_w = c(-1, 0, 1), reamostras = 1000)
dir7 <- r7$posthoc[grepl("direto", r7$posthoc$efeito_de), ]
teste("modelo 15: efeito direto condicional c'1 + c'3·W = .2 ± .25 recuperado (tol .07) e índice perto de .1 (tol .03)",
      perto(dir7$efeito, c(-.05, .2, .45), .07) && perto(r7$efeito$valor[1], .1, .03))

# ---------------------------------------------------------------- equivalência e identidades
cat("== Equivalência com outras implementações e identidades\n")
fit <- lm(y ~ x * w, s2)
rh <- moderacao(s2, "y", "x", "w", ep = "HC3")
suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(boot); library(emmeans); library(lavaan) })
teste("EP HC3 igual a sandwich::vcovHC(type = 'HC3') (tol 1e-10)",
      perto(rh$tabela$ep, sqrt(diag(vcovHC(fit, type = "HC3"))), 1e-10))
bp <- bptest(fit)
teste("Breusch-Pagan igual a lmtest::bptest (tol 1e-8)",
      perto(r2$pressupostos$estatistica[1], unname(bp$statistic), 1e-8) && perto(r2$pressupostos$p[1], unname(bp$p.value), 1e-8))
teste("F da interação = t² do produto, clássico e HC3 (tol 1e-8)",
      perto(r2$teste_interacao$F, r2$tabela$t[4]^2, 1e-8) && perto(rh$teste_interacao$F, rh$tabela$t[4]^2, 1e-8))
em <- summary(emtrends(fit, ~ w, var = "x", at = list(w = c(-1, 0, 1))))
teste("efeitos condicionais iguais a emmeans::emtrends (estimativa e EP, tol 1e-6)",
      perto(r2$posthoc$efeito, em$x.trend, 1e-6) && perto(r2$posthoc$ep, em$SE, 1e-6))
pe <- summary(emtrends(lm(y ~ x * g, s3), ~ g, var = "x"))
teste("moderador categórico igual a emmeans::emtrends (tol 1e-6)", perto(r3$posthoc$efeito, pe$x.trend, 1e-6))

s1b <- s1[1:150, ]
rb <- mediacao(s1b, "x", "m", "y", reamostras = 2000, semente = 77)
est_ab <- function(d, i) { d <- d[i, ]; coef(lm(m ~ x, d))[2] * coef(lm(y ~ x + m, d))[3] }
set.seed(77); bb <- boot(s1b, est_ab, R = 2000)
ci_boot <- boot.ci(bb, type = "perc")$percent[4:5]
teste("IC percentílico igual a boot::boot + boot.ci(type = 'perc') com a mesma semente (tol 1e-10)",
      perto(unlist(efe(rb, "indireto:")[c("ic_inf", "ic_sup")]), ci_boot, 1e-10))

lv <- sem("m1 ~ x\n m2 ~ x + m1\n y ~ x + m1 + m2", s4)
pe_lv <- parameterEstimates(lv)
coef_lv <- function(l, r) pe_lv$est[pe_lv$lhs == l & pe_lv$rhs == r & pe_lv$op == "~"]
teste("serial: caminhos iguais aos do lavaan (modelo saturado, MV = MQO; tol 1e-6)",
      perto(c(lin(r4, "a1")$b, lin(r4, "d21")$b, lin(r4, "b2")$b, lin(r4, "c'")$b),
            c(coef_lv("m1", "x"), coef_lv("m2", "m1"), coef_lv("y", "m2"), coef_lv("y", "x")), 1e-6))

set.seed(107)
s6 <- data.frame(x = rbinom(400, 1, .5), cv = rnorm(400), sx = sample(c("F", "M"), 400, TRUE))
s6$m <- .5 * s6$x + .3 * s6$cv + rnorm(400); s6$y <- .4 * s6$m + .2 * s6$x + .2 * s6$cv + rnorm(400)
s6$y[c(3, 9)] <- NA
r8 <- mediacao(s6, "x", "m", "y", covariaveis = c("cv", "sx"), reamostras = 500)
ok6 <- s6[complete.cases(s6), ]
teste("covariáveis (numérica e fator): a, b, c', c iguais aos de lm (tol 1e-10), identidade fecha e n_excluidos = 2",
      perto(c(lin(r8, "a")$b, lin(r8, "b")$b, lin(r8, "c'")$b, lin(r8, "c")$b),
            c(coef(lm(m ~ x + cv + sx, ok6))[["x"]], coef(lm(y ~ x + m + cv + sx, ok6))[["m"]],
              coef(lm(y ~ x + m + cv + sx, ok6))[["x"]], coef(lm(y ~ x + cv + sx, ok6))[["x"]]), 1e-10) &&
        abs(r8$identidade[["diferenca"]]) < 1e-12 && r8$n_excluidos == 2)

ind5 <- r5$efeito$valor[2:4]; idx5 <- r5$efeito$valor[1]
a3b <- lin(r5, "a3")$b * lin(r5, "b")$b
teste("modelo 7: índice = a3·b da tabela e diferença entre indiretos condicionais = índice × ΔW (tol 1e-12)",
      abs(idx5 - a3b) < 1e-12 && perto(diff(ind5), idx5 * diff(c(-1, 0, 1)), 1e-12))
l15 <- coef(lm(y ~ x + m + w + x:w + m:w, s15))
teste("modelo 15: coeficientes da equação de Y iguais aos de lm com os dois produtos (tol 1e-10)",
      perto(r7$tabela$b[r7$tabela$desfecho == "y"], unname(l15), 1e-10))
l14 <- coef(lm(y3 ~ x + m2 + w + m2:w, s5))
teste("modelo 14: índice = a·b3 com a e b3 de lm (tol 1e-10)",
      abs(r6$efeito$valor[1] - coef(lm(m2 ~ x, s5))[["x"]] * l14[["m2:w"]]) < 1e-10)
rq <- moderacao(s2, "y", "x", "w", valores_w = "percentis", johnson_neyman = FALSE)
teste("valores_w = 'percentis': W sondada nos percentis 16, 50 e 84 (quantile tipo 7)",
      perto(rq$posthoc$w, unname(quantile(s2$w, c(.16, .5, .84))), 1e-12) && identical(rq$posthoc$w_nivel, c("P16", "P50", "P84")))
rpad <- mediacao(s6, "x", "m", "y", padronizado = TRUE, reamostras = 500)
teste("X binário: indireto parcialmente padronizado = ab / DP_Y (tol 1e-12)",
      abs(efe(rpad, "parcialmente padronizado")$valor -
            efe(rpad, "indireto:")$valor / sd(s6$y[complete.cases(s6[c("x", "m", "y")])])) < 1e-12)

cat("== Comportamento da função\n")
set.seed(1); u1 <- runif(1)
set.seed(1); invisible(mediacao(s0, "x", "m", "y", reamostras = 200)); u2 <- runif(1)
teste("a semente da análise não altera o gerador aleatório do usuário", identical(u1, u2))
rv <- mediacao(s0, "x", "m", "y", reamostras = 200, decimal = ",")
teste("decimal = ',' gera vírgula decimal e ponto e vírgula no IC", any(grepl("IC 95% \\[-?[0-9]+,[0-9]+; ", rv$apa)))
s7 <- s0; s7$yb <- as.numeric(s7$y > 0)
teste("Y dicotômica é recusada com orientação para logística",
      inherits(tryCatch(mediacao(s7, "x", "m", "yb", reamostras = 200), error = function(e) e), "error"))
teste("objeto pr_resultado com coluna 'principal' no efeito",
      inherits(r1, "pr_resultado") && "principal" %in% names(r1$efeito) && any(r1$efeito$principal))
pmed <- paste(capture.output(print(r1)), collapse = "\n")
teste("print sem linha de participante (nenhum valor bruto de y no print)", !grepl(format(s1$y[1], digits = 7), pmed, fixed = TRUE))
t0 <- Sys.time(); invisible(mediacao(s1[1:300, ], "x", "m", "y")); dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("        mediação simples, n = 300, 10.000 reamostras: %.1f s\n", dt))
teste("bootstrap de 10.000 reamostras em menos de 10 s (n = 300)", dt < 10)

cat("== Regressões da revisão CoVe (B1, B5, B10, B11)\n")
set.seed(12); nv <- 120
dv <- data.frame(x = rnorm(nv), w = sample(c("baixo", "medio", "alto"), nv, TRUE), cv = sample(c("F", "M"), nv, TRUE))
dv$y <- .3 * dv$x + .5 * dv$x * (dv$w == "alto") + rnorm(nv)
dv$w[c(3, 15, 27, 40, 51, 77)] <- ""; dv$cv[c(5, 9)] <- "  "
rv1 <- moderacao(dv, "y", "x", "w", covariaveis = "cv")
okv <- dv[nzchar(trimws(dv$w)) & nzchar(trimws(dv$cv)), ]
lv1 <- coef(lm(y ~ x * w + cv, okv))
teste("B1: célula de texto vazia sai por listwise em moderacao() (n = 112, 8 excluídos, W com 3 níveis, coeficientes = lm)",
      rv1$n == 112 && rv1$n_excluidos == 8 && identical(rv1$posthoc$w_nivel, c("alto", "baixo", "medio")) &&
        perto(rv1$tabela$b, unname(lv1[c("(Intercept)", "x", "wbaixo", "wmedio", "x:wbaixo", "x:wmedio", "cvM")]), 1e-10) &&
        any(grepl("vazia", rv1$decisoes)))
dvf <- dv; dvf$w <- factor(dvf$w); dvf$cv <- factor(dvf$cv)
rv2 <- moderacao(dvf, "y", "x", "w", covariaveis = "cv")
teste("B1: fator com nível vazio tem o mesmo resultado (n = 112, sem nível '')",
      rv2$n == 112 && perto(rv2$tabela$b, rv1$tabela$b, 1e-12) && !any(rv2$posthoc$w_nivel == ""))
dg <- data.frame(g = sample(c("ctrl", "trat"), nv, TRUE), w = rnorm(nv), m = rnorm(nv)); dg$y <- .4 * dg$m + rnorm(nv)
dg$g[c(2, 4)] <- ""
rv3 <- mediacao(dg, "g", "m", "y", reamostras = 200)
rv4 <- mediacao_moderada(transform(dg, gw = ifelse(g == "", "", ifelse(w > 0, "alto", "baixo"))), "w", "m", "y", "gw",
                         modelo = 14, reamostras = 200)
teste("B1: X de texto binário com vazios funciona em mediacao() e mediacao_moderada() (2 excluídos, codificação 0/1)",
      rv3$n_excluidos == 2 && any(grepl("g codificada como 0 = 'ctrl' e 1 = 'trat'", rv3$decisoes, fixed = TRUE)) &&
        rv4$n_excluidos == 2 && any(grepl("gw codificada como 0 = 'alto' e 1 = 'baixo'", rv4$decisoes, fixed = TRUE)))

aviso_inc <- function(r) any(grepl("mediação inconsistente", r$avisos))
justificado <- logical(); nulo_opostos <- logical()
for (sd_ in 1:40) {
  set.seed(sd_); dn <- data.frame(x = rnorm(200)); dn$m <- .4 * dn$x + rnorm(200); dn$y <- .4 * dn$m + rnorm(200)
  rn <- mediacao(dn, "x", "m", "y", reamostras = 200)
  cl <- lin(rn, "c'"); sig_c <- cl$ic_inf > 0 || cl$ic_sup < 0
  if (aviso_inc(rn)) justificado <- c(justificado, sig_c && cl$b < 0)
  if (!sig_c && cl$b < 0) nulo_opostos <- c(nulo_opostos, aviso_inc(rn))
}
teste("B5: com c' verdadeiro = 0, c' negativo e IC com zero nunca gera aviso de inconsistência; todo aviso tem c' com IC sem zero",
      length(nulo_opostos) > 5 && !any(nulo_opostos) && all(justificado))
set.seed(109); di <- data.frame(x = rnorm(400)); di$m <- .5 * di$x + rnorm(400); di$y <- .5 * di$m - .5 * di$x + rnorm(400)
teste("B5: mediação inconsistente verdadeira (ab = .25, c' = -.5) gera o aviso",
      aviso_inc(mediacao(di, "x", "m", "y", reamostras = 500)))

set.seed(110); dq <- data.frame(x = rnorm(2001)); dq$m <- .3 * dq$x + rnorm(2001); dq$y <- .3 * dq$m + rnorm(2001)
rq2 <- mediacao(dq, "x", "m", "y")
teste("B10: bootstrap sequencial (n × reamostras > 2e7) avisa em $avisos, como a referência descreve",
      any(grepl("uma reamostra por vez", rq2$avisos)) && any(grepl("uma reamostra por vez", rq2$decisoes)))

set.seed(4); dr <- data.frame(idade = runif(300, 18, 80), renda = rlnorm(300, log(4000), .6))
dr$y <- 10 + .05 * dr$idade + .0005 * dr$renda - .00001 * dr$idade * dr$renda + rnorm(300)
rr <- moderacao(dr, "y", "idade", "renda", decimal = ",")
teste("B11: b3 = -0,0000089 e EP = 0,0000012 saem com sinal e tamanho (não '0,000'), com nota em $decisoes",
      grepl("b = -0,0000089; EP = 0,0000012", rr$apa[1], fixed = TRUE) && any(grepl("algarismos significativos", rr$decisoes)))
fn <- function(v) .pr_num(v, .ee_proc_casas(v, 3), TRUE, ".")
teste("B11: casas só ampliam quando o valor arredondaria para zero; limite de IC perto de zero mantém o sinal",
      .ee_proc_casas(.241, 3) == 3 && .ee_proc_casas(-.00049, 3) == 5 && .ee_proc_casas(0, 3) == 3 &&
        .ee_proc_ic(-.0003, .52, fn, ".", .95, "pt") == "IC 95% [-0.00030, 0.520]")

cat("== Correções da verificação final (achados G e H)\n")
# G-03: grupo gravado 1/2/3 no W da moderação, no X da mediação, no W da mediação
# moderada e numa covariável. Dados da reprodução cove_final/G/g_rep_G03_codigo_numerico.R;
# referência do F conjunto: anova() de dois lm com factor(), calculada aqui.
set.seed(7); n_g3 <- 150
d_g3 <- data.frame(escolaridade = sample(1:3, n_g3, TRUE), idade = round(rnorm(n_g3, 35, 8)), x01 = rbinom(n_g3, 1, .5))
d_g3$y2 <- 5 + .2 * d_g3$idade + rnorm(n_g3)
d_g3$m <- .5 * d_g3$x01 + rnorm(n_g3); d_g3$y <- .4 * d_g3$m + rnorm(n_g3)
aviso_cod <- function(r, v) any(grepl(paste0(v, " é numérica com só 3 valores inteiros distintos (1, 2, 3)"), r$avisos, fixed = TRUE))
teste("G-03 moderação: W = escolaridade 1/2/3 numérica gera aviso (sugere factor()); W contínua não", {
  mo <- moderacao(d_g3, "y2", "idade", "escolaridade")
  aviso_cod(mo, "escolaridade") && any(grepl("factor(dados$escolaridade)", mo$avisos, fixed = TRUE)) &&
    !any(grepl("valores inteiros distintos", r2$avisos)) })
teste("G-03 moderação: com W em fator o aviso some e o F conjunto (gl1 = 2) = anova(lm(y ~ x + W), lm(y ~ x * W)) (tol 1e-8)", {
  d_f <- transform(d_g3, escolaridade = factor(escolaridade))
  mf <- moderacao(d_f, "y2", "idade", "escolaridade")
  a_f <- anova(lm(y2 ~ idade + escolaridade, d_f), lm(y2 ~ idade * escolaridade, d_f))
  !any(grepl("valores inteiros distintos", mf$avisos)) && mf$teste_interacao$gl1 == 2 && perto(mf$teste_interacao$F, a_f$F[2], 1e-8) })
teste("G-03 mediação: X = escolaridade 1/2/3 gera aviso; X 0/1 continua aceito sem aviso; covariável 1/2/3 também avisa", {
  me <- mediacao(d_g3, "escolaridade", "m", "y", reamostras = 200)
  mb <- mediacao(d_g3, "x01", "m", "y", reamostras = 200)
  mc <- mediacao(d_g3, "x01", "m", "y", covariaveis = "escolaridade", reamostras = 200)
  aviso_cod(me, "escolaridade") && any(grepl("dois grupos por vez", me$avisos)) &&
    !any(grepl("valores inteiros distintos", mb$avisos)) && aviso_cod(mc, "escolaridade") && any(grepl("covariável em fator", mc$avisos)) })
teste("G-03 mediação moderada: W = escolaridade 1/2/3 gera aviso; W contínua (s5$w) não", {
  mm <- mediacao_moderada(d_g3, "idade", "m", "y", "escolaridade", modelo = 7, reamostras = 200)
  aviso_cod(mm, "escolaridade") && !any(grepl("valores inteiros distintos", r5$avisos)) })

# G-11: número de reamostras na frase sem separador de milhar ambíguo. Esperado escrito à
# mão: en com decimal "." = "1,000"; pt com decimal "," = "1.000"; en com "," = "1000".
teste("G-11 reamostras: idioma = \"en\" com decimal = \",\" sai \"1000 resamples\" (não \"1,000\" ao lado de \"0,xxx\"); en/\".\" = \"1,000\"; pt/\",\" = \"1.000\"", {
  a_en_v <- mediacao(s0, "x", "m", "y", reamostras = 1000, idioma = "en", decimal = ",")$apa[1]
  a_en_p <- mediacao(s0, "x", "m", "y", reamostras = 1000, idioma = "en")$apa[1]
  a_pt_v <- mediacao(s0, "x", "m", "y", reamostras = 1000, decimal = ",")$apa[1]
  grepl("percentile bootstrap, 1000 resamples", a_en_v, fixed = TRUE) && !grepl("1,000", a_en_v, fixed = TRUE) &&
    grepl("1,000 resamples", a_en_p, fixed = TRUE) && grepl("1.000 reamostras", a_pt_v, fixed = TRUE) &&
    .ee_proc_int(10000, ",", "en") == "10000" && .ee_proc_int(10000, ".", "pt") == "10000" })

# H-06: ΔR², f² e R² da moderação (e ΔR² da mediação moderada) sem IC ganham a nota.
teste("H-06 moderação e mediação moderada: nota \"sem IC: não calculado\" em $decisoes (também em inglês)", {
  any(grepl("ΔR², f² e R² do modelo sem IC: não calculado", r2$decisoes, fixed = TRUE)) &&
    any(grepl("ΔR² das interações sem IC: não calculado", r5$decisoes, fixed = TRUE)) &&
    any(grepl("without CI: not computed", moderacao(s2, "y", "x", "w", idioma = "en", johnson_neyman = FALSE)$decisoes, fixed = TRUE)) })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
