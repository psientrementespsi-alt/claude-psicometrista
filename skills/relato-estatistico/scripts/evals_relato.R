# evals_relato.R, relato-estatistico 1.5.0 (1.5: g de Hedges, omega2p, rede): auditor statcheck, formatador por perfil, frases do R, efeitos, .docx
# Veredictos de referencia obtidos com o codigo do pacote statcheck (repositorios-relato-estatistico/statcheck/R,
# funcao statcheck(), padroes do pacote) em 29/10/2026; quando o repositorio existe na maquina, a comparacao e refeita ao vivo.
args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else if (file.exists("auditor_statcheck.R")) "." else "scripts"
for (f in c("auditor_statcheck.R", "formatador_apa.R", "frases_apa.R", "checagem_efeitos.R", "docx_apa.R", "gerador_tabelas_apa.R")) source(file.path(pasta, f))
resultados <- list()
teste <- function(nome, expr) { ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome)); resultados[[nome]] <<- ok; invisible(ok) }
aud <- function(s) suppressMessages({ invisible(utils::capture.output(o <- auditar_texto_statcheck(texto = s))); o })
ref <- data.frame(caso = c("t(28) = 2.20, p = .036", "t(28) = 2.20, p = .06", "t(28) = 2.20, p < .05", "F(2, 56) = 1.03, p < .05",
  "F(2, 56) = 3.03, p = .056", "r(112) = .18, p = .062", "r(112) = .53, p < .001", "chi2(40) = 25.92, p = .958", "chi2(3) = 8.05, p = .045",
  "t(110) = 2.15, p = .033", "t(110) = 1.63, p = .106", "F(3, 110) = 27.54, p < .001", "t(95) = 2.66, p = .009", "chi2(186) = 318.84, p < .001",
  "the test z = 1.96, p = .05", "the test z = 2.50, p = .012", "t(40) = -2.31, p = .026", "F(1, 98) = 4.12, p = .045", "r(50) = .25, p = .02",
  "t(15) = 1.10, p < .01", "t(20) = 2.09, p = .05", "F(1, 20) = 4.35, p = .05", "r(30) = .35, p = .049", "chi2(2) = 5.99, p = .05"),
  erro = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  grosseiro = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE)
v <- lapply(ref$caso, aud)
teste("24 casos: um teste extraido por caso", all(vapply(v, nrow, 1L) == 1L))
teste("inconsistencias iguais as do statcheck nos 24 casos", identical(vapply(v, function(o) o$inconsistente[1], TRUE), ref$erro))
teste("erros grosseiros iguais aos do statcheck nos 24 casos", identical(vapply(v, function(o) o$erro_grosseiro[1], TRUE), ref$grosseiro))
pt <- aud("A regressão foi significativa, *F*(3, 110) = 27,54, *p* < ,001. Estresse, *t*(110) = 4,79, *p* < ,001; idade, t(110) = 2,15, p = ,033; χ²(40) = 25,92, p = ,958; r(112) = ,18, p = ,062; t(40) = −2,31, p = ,026.")
teste("portugues: virgula decimal, χ², sinal de menos tipografico e italico de markdown (6 testes, nenhum erro)", nrow(pt) == 6 && !any(pt$inconsistente))
pe <- aud("O efeito foi t(28) = 2,20, p = ,06.")
teste("portugues: erro grosseiro detectado com virgula decimal", nrow(pe) == 1 && pe$erro_grosseiro)
teste("arredondamento: t(110) = 2.15 com p = .034 tambem e consistente (intervalo de arredondamento)", !aud("t(110) = 2.15, p = .034")$inconsistente)
teste("p relatado zero e erro (pZeroError do statcheck)", aud("t(30) = 5.20, p = .000")$inconsistente)

# ---- 1.3.0: formatador por perfil ----
teste("perfis: apa7_relato tira zero de beta e mantem em RMSEA; apa7_estrito mantem em beta; psicometria_r tira em RMSEA",
  fmt_apa(.31, "beta") == ".31" && fmt_apa(.048, "rmsea") == "0.048" && fmt_apa(.31, "beta", "apa7_estrito") == "0.31" && fmt_apa(.048, "rmsea", "psicometria_r") == ".048")
teste("p: limites < .001 e > .999, e zero negativo nao vira -0.00", fmt_apa(.0004, "p") == "< .001" && fmt_apa(.9996, "p") == "> .999" && fmt_apa(-.001, "livre") == "0.00")
teste("ptbr: virgula decimal, IC com ponto e virgula; perfis de periodico herdam e saem nao verificados",
  fmt_ic_apa(.12, .45, "r", "apa7_ptbr") == "[0,12; 0,45]" && fmt_ic_apa(.12, .45, "r", "apa7_ptbr_sem_zero") == "[,12; ,45]" &&
  !perfil_apa("rbtc")$verificado && perfil_apa("rbtc")$decimal == "," && perfil_apa("prc")$herdado_de == "apa7_estrito")
teste("perfil personalizado altera so o pedido", { pp <- perfil_apa_personalizado("apa7_relato", remover_sem_zero = "beta"); fmt_apa(.31, "beta", pp) == "0.31" && fmt_apa(.31, "r", pp) == ".31" })
# ---- 1.3.0: frases geradas do objeto passam no auditor ----
set.seed(11); x <- rnorm(90); g <- rep(1:2, 45); y <- .5 * x + rnorm(90); z <- rnorm(90)
for (pf in c("apa7_relato", "apa7_ptbr", "apa7_ptbr_sem_zero")) {
  fr <- c(frases_apa(t.test(y ~ g, var.equal = TRUE), pf), frases_apa(cor.test(x, y), pf), frases_apa(lm(y ~ x + z), pf), frases_apa(chisq.test(table(g, x > 0)), pf))
  a <- aud(paste("Resultado:", fr, ".", collapse = " "))
  teste(sprintf("frases_apa no perfil %s: 6 testes reconhecidos (t, r, F, 2 coeficientes, chi2), nenhum inconsistente", pf), nrow(a) == 6 && !any(a$inconsistente))
}
teste("frases_apa: p da regressao bate com summary(lm)", { f <- lm(y ~ x + z); grepl(fmt_apa(summary(f)$coefficients["x", 4], "p"), frases_apa(f)["x"], fixed = TRUE) })
if (requireNamespace("lavaan", quietly = TRUE)) {
  fl <- lavaan::cfa("vis =~ x1 + x2 + x3\n txt =~ x4 + x5 + x6", data = lavaan::HolzingerSwineford1939, estimator = "MLR")
  fr <- frases_apa(fl)
  teste("frases_apa lavaan MLR: usa indices robustos e chi2 passa no auditor", attr(fr, "versao_indices") == "robusta" && nrow(aud(paste("Ajuste:", fr))) == 1)
  teste("frases_apa lavaan: modelo saturado (gl = 0) e recusado", inherits(tryCatch(frases_apa(lavaan::cfa("vis =~ x1 + x2 + x3", data = lavaan::HolzingerSwineford1939)), error = function(e) e), "error"))
}
# ---- 1.4.0: frases_apa para glm, lmer, coxph, rma e roc ----
df_reg <- data.frame(y_bin = c(rep(0, 20), rep(1, 20)), x1 = rnorm(40), grp = rep(1:4, each = 10))
fit_glm <- stats::glm(y_bin ~ x1, data = df_reg, family = stats::binomial)
f_glm <- frases_apa(fit_glm)
teste("frases_apa glm: relata OR, EP, z e p para cada preditor",
  "x1" %in% names(f_glm) && grepl("OR =", f_glm["x1"]) && grepl("EP", f_glm["x1"]) && grepl("z", f_glm["x1"]))

if (requireNamespace("lme4", quietly = TRUE)) {
  fit_lmer <- suppressWarnings(lme4::lmer(x1 ~ (1 | grp), data = df_reg))
  f_lmer <- frases_apa(fit_lmer)
  teste("frases_apa lmer: relata ICC, variancia do intercepto e residuo",
    "modelo" %in% names(f_lmer) && grepl("ICC =", f_lmer["modelo"]) && grepl("τ₀₀ =", f_lmer["modelo"]))
} else cat("PULADO  frases_apa lmer (lme4 ausente)\n")

if (requireNamespace("metafor", quietly = TRUE)) {
  fit_meta <- suppressWarnings(metafor::rma(yi = c(0.2, 0.4, 0.3), vi = c(0.04, 0.05, 0.03), method = "FE"))
  f_meta <- frases_apa(fit_meta)
  teste("frases_apa metafor: relata efeito combinado, IC, heterogeneidade Q e I2",
    "resumo" %in% names(f_meta) && grepl("Tamanho de efeito combinado =", f_meta["resumo"]) && grepl("I\\*²", f_meta["resumo"]))
} else cat("PULADO  frases_apa metafor (metafor ausente)\n")

if (requireNamespace("pROC", quietly = TRUE)) {
  roc_ex <- suppressMessages(pROC::roc(df_reg$y_bin, df_reg$x1, quiet = TRUE))
  f_roc <- frases_apa(roc_ex)
  teste("frases_apa roc: relata AUC e ponto de corte de Youden",
    grepl("AUC =", f_roc) && grepl("Ponto de corte ótimo", f_roc))
} else cat("PULADO  frases_apa roc (pROC ausente)\n")
# ---- 1.3.0: checagem de efeitos (Lakens, 2013) ----
ef <- auditar_efeitos_texto("O grupo diferiu, t(78) = 2.40, p = .019, d = 0.54, IC 95% [0.09, 0.98]. A ANOVA foi F(2, 77) = 4.10, p = .020.")
teste("efeitos: frase com d e IC passa; F sem efeito nem IC e marcado", nrow(ef) == 2 && ef$tem_efeito[1] && ef$tem_ic[1] && !ef$tem_efeito[2] && !ef$tem_ic[2])
teste("d reconstituido de t (independentes e pareado) e eta2p de F", conferir_d_de_t("2.40", 40, 40, "0.54")$consistente && !conferir_d_de_t("2.40", 40, 40, "0.80")$consistente &&
  conferir_d_de_t("3.00", 25, NULL, "0.60")$consistente && conferir_eta2p_de_F("4.10", 2, 77, ".10")$consistente && !conferir_eta2p_de_F("4.10", 2, 77, ".20")$consistente)
# ---- 1.3.0: .docx e tabela de correlacoes ----
arq <- tempfile(fileext = ".docx")
escrever_tabelas_docx(arq, list(numero = 1, titulo = "T *it*", headers = c("V", "*M*", "1"), linhas = list(c("A & B", "0.12", "—"), c("C", "1.50", ".30")), nota = "*Nota.* x"))
xml <- paste(readLines(unz(arq, "word/document.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
teste(".docx: partes OOXML presentes, tres linhas (6 bordas em 3 colunas), escape de &, sem borda vertical",
  all(c("[Content_Types].xml", "_rels/.rels", "word/document.xml") %in% utils::unzip(arq, list = TRUE)$Name) &&
  lengths(regmatches(xml, gregexpr('w:val="single"', xml))) == 9 && grepl("A &amp; B", xml, fixed = TRUE) && !grepl('w:left w:val="single"', xml))
dd <- data.frame(a = c(x[1:80], rep(NA, 10)), b = y, c = z)
teste("tabela de correlacoes relata n por par quando ha ausentes", grepl("entre 80 e 90", gerar_tabela_correlacoes_apa(dd, c("a", "b", "c"))))

SC <- Sys.getenv("PR_STATCHECK_R", "~/Desktop/psicomtria rrr/repositorios-relato-estatistico/statcheck/R")
if (dir.exists(SC)) {
  es <- new.env(); for (f in c("regex.R","constants.R","helpers-parse-stats.R","extract_stats.R","extract_1tail.R","process_stats.R","compute_p.R","error_test.R","decision_error_test.R","calc_APA_factor.R","statcheck.R"))
    if (file.exists(file.path(SC, f))) sys.source(file.path(SC, f), envir = es)
  r0 <- as.data.frame(suppressWarnings(suppressMessages(es$statcheck(paste(ref$caso, collapse = ". "), messages = FALSE))))
  teste("ao vivo: o statcheck do repositorio da os mesmos erros que a tabela de referencia", nrow(r0) == 24 && sum(r0$error) == sum(ref$erro) && sum(r0$decision_error) == sum(ref$grosseiro))
} else cat("PULADO  comparacao ao vivo com o statcheck (repositorio ausente)\n")
# ---- v1.5.0: g de Hedges e omega2 parcial (Lakens, 2013); frases de rede (Epskamp, Borsboom & Fried, 2018) ----
teste("Lakens (2013), exemplo dos filmes: t(18) = 2,52 com 10 e 10 da d_s = 1,13 e g_s = 1,08",
  conferir_d_de_t("2.52", 10, 10, "1.13")$consistente && conferir_d_de_t("2.52", 10, 10, "1.08", hedges = TRUE)$consistente &&
  !conferir_d_de_t("2.52", 10, 10, "1.13", hedges = TRUE)$consistente)
teste("Lakens (2013), eq. 20: omega2p = ,21 com MS 5 e 0,789, gl 1 e N = 20 (F = 6,34)",
  conferir_omega2p_de_F(sprintf("%.2f", 5 / .789), 1, 20, ".21")$consistente && !conferir_omega2p_de_F(sprintf("%.2f", 5 / .789), 1, 20, ".26")$consistente)
fr_cs <- frases_rede(structure(list(graph = matrix(0, 2, 2), nPerson = 359, default = "EBICglasso", arguments = list(tuning = .5)), class = c("bootnetResult", "list")),
  cs = c(betweenness = .05, closeness = .05, strength = .44), perfil = "apa7_relato")
teste("frases_rede: CS do exemplo de Epskamp et al. (2018), forca ,44 abaixo de ,50 e intermediacao ,05 abaixo de ,25",
  grepl("for\u00e7a = 0.44 (abaixo de 0.50", fr_cs["estabilidade"], fixed = TRUE) && grepl("intermedia\u00e7\u00e3o = 0.05 (abaixo de 0.25", fr_cs["estabilidade"], fixed = TRUE))
if (requireNamespace("bootnet", quietly = TRUE) && requireNamespace("MASS", quietly = TRUE)) {
  set.seed(1); Sg <- diag(6); Sg[Sg == 0] <- .3; Xr <- MASS::mvrnorm(300, rep(0, 6), Sg); colnames(Xr) <- paste0("v", 1:6)
  net <- suppressMessages(suppressWarnings(bootnet::estimateNetwork(as.data.frame(Xr), default = "EBICglasso", tuning = .5, verbose = FALSE)))
  fr <- frases_apa(net); ne <- sum(net$graph[upper.tri(net$graph)] != 0)
  teste("frases_apa(bootnetResult): nos, N, arestas e gama saem do objeto", grepl("rede de 6 n\u00f3s", fr["rede"], fixed = TRUE) &&
    grepl(sprintf("*N* = 300, teve %d de 15", ne), fr["rede"], fixed = TRUE) && grepl("\u03b3 = 0.50", fr["rede"], fixed = TRUE))
} else cat("PULADO  frases_apa com objeto do bootnet (bootnet ausente)\n")
cat(sprintf("\nTotal relato-estatistico 1.5.0: %d | Passaram: %d | Falharam: %d\n", length(resultados), sum(unlist(resultados)), sum(!unlist(resultados))))
