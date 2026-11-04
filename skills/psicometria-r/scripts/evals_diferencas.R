# evals_diferencas.R, testes do módulo de diferenças (est_diferencas.R)
#
# Rodar da pasta da skill:
#   PR_DADOS_JASP="/caminho/para/Data Library" Rscript scripts/evals_diferencas.R
# PR_DADOS_JASP aponta para a biblioteca de dados do JASP (a pasta que contém
# "2. T-Tests" e "3. ANOVA"; qualquer pasta acima delas também serve, a busca é
# recursiva). Sem a variável, os testes de âncora do livro são pulados com aviso.
# Três tipos de teste, nesta ordem: (1) âncoras JASP 0.95.1 impressas no livro do
# usuário (Vol. 4 e 5, com a linha mdX:linha); (2) verdade conhecida em dados simulados
# com semente fixa; (3) equivalência com outra implementação. Tolerância declarada em
# cada teste: para âncora, metade da última casa impressa pelo JASP.
# Saída: uma linha PASSOU/FALHOU por teste e o resumo; código de saída 1 se algo falhar.

options(rgl.useNULL = TRUE)  # heplots carrega o rgl; sem tela, evita avisos de X11
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
ef <- function(r, ind, termo = NULL) {
  e <- r$efeito
  if (!is.null(termo)) e <- e[e$termo == termo, ]
  e[grepl(ind, e$indice, fixed = TRUE), , drop = FALSE][1, ]
}
pres <- function(r, txt) r$pressupostos[grepl(txt, r$pressupostos$checagem, fixed = TRUE), , drop = FALSE][1, ]

# ============================================================ 1. Âncoras JASP (livro)
dra <- jasp("Directed Reading Activities.csv")
if (is.null(dra)) pular("Directed Reading Activities") else {
  cat("== Âncoras: Directed Reading Activities (Vol. 4, cap. 13)\n")
  r <- teste_t_independente(dra, "drp", "group")
  st <- r$tabela[1, ]; we <- r$tabela[2, ]
  teste("Student t(42) = -2,267, p = ,029",
        perto(st$t, -2.267, .0005) && st$gl == 42 && perto(st$p, .029, .0005))
  teste("Welch t = -2,311, gl = 37,86, p = ,026",
        perto(we$t, -2.311, .0005) && perto(we$gl, 37.86, .005) && perto(we$p, .026, .0005))
  d <- ef(r, "Student")
  teste("d = -0,684, IC 95% [-1,290; -0,071]",
        perto(d$valor, -.684, .0005) && perto(c(d$ic_inf, d$ic_sup), c(-1.290, -.071), .0005) && d$principal)
  lv <- pres(r, "Levene")
  teste("Levene F(1, 42) = 2,362, p = ,132: o do livro é o centrado na média",
        perto(lv$estatistica, 2.362, .0005) && perto(lv$p, .132, .0005))
  teste("Levene centrado na mediana (padrão do painel t do JASP) dá 2,342: não é o valor do livro",
        perto(r$levene_outro_centro[["F"]], 2.342, .0005))
  teste("Shapiro-Wilk dos resíduos W = 0,973, p = ,396",
        perto(pres(r, "resíduos")$estatistica, .973, .0005) && perto(pres(r, "resíduos")$p, .396, .0005))
  teste("Shapiro-Wilk por grupo: Treat W = 0,966, p = ,652; Control W = 0,972, p = ,732",
        perto(pres(r, "em Treat")$estatistica, .966, .0005) && perto(pres(r, "em Treat")$p, .652, .0005) &&
          perto(pres(r, "em Control")$estatistica, .972, .0005) && perto(pres(r, "em Control")$p, .732, .0005))
  ds <- r$descritivas
  teste("Treat M = 51,48, DP = 11,01 (n 21); Control M = 41,52, DP = 17,15 (n 23)",
        perto(ds$media[ds$grupo == "Treat"], 51.48, .005) && perto(ds$dp[ds$grupo == "Treat"], 11.01, .005) &&
          perto(ds$media[ds$grupo == "Control"], 41.52, .005) && perto(ds$dp[ds$grupo == "Control"], 17.15, .005) &&
          all(ds$n == c(23, 21)))
  teste("diferença 9,95, IC 95% [1,09; 18,82] em módulo",
        perto(abs(st$dif_media), 9.95, .005) && perto(sort(abs(c(st$ic_inf, st$ic_sup))), c(1.09, 18.82), .005))
  teste("rota clássica escolhe Student (Levene p = ,132) e rota welch escolhe Welch",
        r$tabela$escolhido[1] && teste_t_independente(dra, "drp", "group", rota = "welch")$tabela$escolhido[2])
  teste("frase APA completa do Directed Reading (texto exato)",
        identical(r$apa, paste0("Um teste t de Student para amostras independentes indicou média maior em Treat (M = 51.48, DP = 11.01) ",
                                "do que em Control (M = 41.52, DP = 17.15), t(42) = 2.27, p = .029; diferença média 9.95, IC 95% [1.09, 18.82], ",
                                "d = 0.68, IC 95% [0.07, 1.29] (efeito moderado).")))
  teste("frase APA com t(42) = 2.27, p = .029 e d = 0.68 (efeito moderado)",
        grepl("t(42) = 2.27", r$apa, fixed = TRUE) && grepl("p = .029", r$apa, fixed = TRUE) &&
          grepl("d = 0.68", r$apa, fixed = TRUE) && grepl("moderado", r$apa))
  # G-05 (CoVe final; reprodução g_rep_G05_unilateral.R): H1 Control > Treat com os dados no
  # sentido oposto. t de Student -2,267; p unilateral =
  # 1 - ,0285/2 = ,986; diferença -9,95. A frase antes dava "t(42) = 2.27, p = .986"
  # sem dizer a hipótese.
  ru <- teste_t_independente(dra, "drp", "group", alternativa = "maior")$apa
  teste("G-05 unilateral oposto no Directed Reading: H1 declarada, t(42) = -2.27, p = .986, diferença Control - Treat = -9.95",
        grepl("unilateral (H1: média de Control maior que a de Treat)", ru, fixed = TRUE) && grepl("t(42) = -2.27, p = .986", ru, fixed = TRUE) &&
          grepl("diferença média Control - Treat = -9.95", ru, fixed = TRUE) && grepl("no sentido oposto ao de H1", ru, fixed = TRUE))
  # G-11: $decisoes com decimal = ",": Levene 2,362 -> 2,36 e a versão na mediana
  # 2,342 -> 2,34 (teste acima), com "; " entre os itens e p sem zero à esquerda.
  dv <- paste(teste_t_independente(dra, "drp", "group", decimal = ",")$decisoes, collapse = " ")
  teste("G-11 decisões do t com decimal = ',': 'Levene F(1; 42) = 2,36; p = ,132' e 'a outra versão dá F = 2,34; p = ,133'",
        grepl("Levene F(1; 42) = 2,36; p = ,132", dv, fixed = TRUE) && grepl("a outra versão dá F = 2,34; p = ,133", dv, fixed = TRUE))
}

wg <- jasp("Weight Gain.csv")
if (is.null(wg)) pular("Weight Gain") else {
  cat("== Âncoras: Weight Gain (Vol. 4, cap. 13)\n")
  r <- teste_t_pareado(wg, "Weight Before", "Weight After")
  tt <- r$tabela
  teste("t(15) = 10,84 (o JASP mostra -10,84: ordem antes - depois), p < ,001",
        perto(tt$t, 10.84, .005) && tt$gl == 15 && tt$p < .001)
  teste("IC 95% da diferença [8,36; 12,46] e diferença média 10,41",
        perto(c(tt$ic_inf, tt$ic_sup), c(8.36, 12.46), .005) && perto(tt$dif_media, 10.41, .005))
  dz <- ef(r, "d_z")
  teste("d_z = 2,710, IC 95% [1,623; 3,778] (o d da tela, sinal trocado; tol .0005)",
        perto(dz$valor, 2.710, .0005) && perto(c(dz$ic_inf, dz$ic_sup), c(1.623, 3.778), .0005))
  dp <- r$efeito[r$efeito$principal, ]
  teste("d principal na métrica original = 0,47 (Caixa 13.2; tol .005)",
        nrow(dp) == 1 && perto(dp$valor, .47, .005) && grepl("métrica original", dp$indice))
  teste("d_rm ('Correct for correlation' do JASP) = 0,45, diferente do 0,47 do livro",
        perto(ef(r, "d_rm")$valor, .445, .0005))
  teste("Shapiro-Wilk da diferença W = 0,938, p = ,325",
        perto(pres(r, "diferença")$estatistica, .938, .0005) && perto(pres(r, "diferença")$p, .325, .0005))
  ds <- r$descritivas
  teste("antes M = 144,64, DP = 22,70; depois M = 155,05, DP = 21,44; diferença DP = 3,84; r = ,99",
        perto(ds$media[1:2], c(144.64, 155.05), .005) && perto(ds$dp, c(22.70, 21.44, 3.84), .005) && perto(r$r, .99, .005))
  teste("frase APA completa do Weight Gain com decimal = ',' (texto exato; M = 155,05 como o JASP arredonda)",
        identical(teste_t_pareado(wg, "Weight Before", "Weight After", decimal = ",")$apa,
                  paste0("Um teste t para amostras dependentes indicou que os escores aumentaram do momento Weight Before (M = 144,64; DP = 22,70) ",
                         "para Weight After (M = 155,05; DP = 21,44), t(15) = 10,84; p < ,001; diferença média 10,41; IC 95% [8,36; 12,46]; ",
                         "d = 0,47; IC 95% [0,39; 0,56] (d na métrica original, DP combinado das duas medidas; efeito pequeno, próximo do moderado).")))
  teste("frase APA com M = 155,05 (arredondamento do JASP), t(15) = 10,84 e d = 0,47 (pequeno, próximo do moderado)",
        { a <- teste_t_pareado(wg, "Weight Before", "Weight After", decimal = ",")$apa
          grepl("155,05", a) && grepl("t(15) = 10,84", a, fixed = TRUE) && grepl("d = 0,47", a, fixed = TRUE) && grepl("próximo do moderado", a) })
}

eye <- jasp("Response to Eye Color.csv")
if (is.null(eye)) pular("Response to Eye Color") else {
  cat("== Âncoras: Response to Eye Color (Vol. 4, cap. 14)\n")
  r <- anova_um_fator(eye, "Score", "Group")
  cl <- r$tabela[1, ]
  teste("F(3, 218) = 2,894, p = ,036, SQ fator = 24,42, SQ resíduo ≈ 613",
        perto(cl$F, 2.894, .0005) && cl$gl1 == 3 && cl$gl2 == 218 && perto(cl$p, .036, .0005) &&
          perto(cl$sq, 24.42, .005) && perto(cl$sq_residuo, 613, .5))
  teste("ω² = 0,025 (principal) e η²p = 0,038",
        perto(ef(r, "ω²")$valor, .025, .0005) && ef(r, "ω²")$principal && perto(ef(r, "η²")$valor, .038, .0005))
  lv <- pres(r, "Levene")
  teste("Levene F(3, 218) = 0,519, p = ,670", perto(lv$estatistica, .519, .0005) && perto(lv$p, .670, .0005))
  ph <- r$posthoc
  dg <- ph[ph$comparacao == "Down - Green", ]; bg <- ph[ph$comparacao == "Blue - Green", ]
  teste("Tukey Down-Green: dif = -0,752, t = -2,321, p aj = ,096, d = -0,449",
        perto(c(dg$dif, dg$t, dg$p_ajustado, dg$d_qme), c(-.752, -2.321, .096, -.449), .0005))
  teste("Tukey Blue-Green: dif = -0,666, t = -2,376, p aj = ,085, d = -0,397",
        perto(c(bg$dif, bg$t, bg$p_ajustado, bg$d_qme), c(-.666, -2.376, .085, -.397), .0005))
  teste("demais pares com p ajustado > ,36", all(ph$p_ajustado[!ph$comparacao %in% c("Down - Green", "Blue - Green")] > .36))
  teste("g de Hedges par a par ,46 e ,39 (o g que o corpo cita; tol .005)",
        perto(abs(c(dg$g_par, bg$g_par)), c(.46, .39), .005))
  ds <- r$descritivas
  teste("médias Blue 3,19, Brown 3,72, Down 3,11, Green 3,86; n 67/37/41/77",
        perto(ds$media, c(3.19, 3.72, 3.11, 3.86), .005) && all(ds$n == c(67, 37, 41, 77)))
  teste("frase APA: F(3, 218) = 2.89, p = .036, ω² = .025",
        grepl("F(3, 218) = 2.89", r$apa[1], fixed = TRUE) && grepl("p = .036", r$apa[1], fixed = TRUE) && grepl("ω² = .025", r$apa[1], fixed = TRUE))
}

stp <- jasp("Stroop.csv")
if (is.null(stp)) pular("Stroop") else {
  cat("== Âncoras: Stroop, bloco Stroop (Vol. 4, cap. 14)\n")
  med <- c("Congruent-Stroop", "Incongruent-Stroop", "Neutral-Stroop")
  r <- anova_medidas_repetidas(stp, med, intra = list(condicao = c("Congruente", "Incongruente", "Neutra")))
  mc <- pres(r, "Mauchly")
  teste("Mauchly W = 0,670, χ²(2) = 6,803, p = ,033",
        perto(mc$estatistica, .670, .0005) && mc$gl == "2" && perto(mc$p, .033, .0005) && grepl("6.803", mc$leitura, fixed = TRUE))
  teste("ε GG = 0,752, HF = 0,806, LB = 0,500",
        grepl("GG = 0.752", mc$leitura, fixed = TRUE) && grepl("HF = 0.806", mc$leitura, fixed = TRUE) && grepl("LB = 0.5", mc$leitura, fixed = TRUE))
  tb <- r$tabela
  teste("sem correção F(2, 36) = 9,18, p < ,001", perto(tb$F, 9.183, .0005) && tb$gl1 == 2 && tb$gl2 == 36 && tb$p < .001)
  teste("GG escolhido na fronteira (ε = ,752): F(1,504; 27,072) = 9,183, p = ,002",
        tb$correcao == "Greenhouse-Geisser" && perto(c(tb$gl1_corr, tb$gl2_corr), c(1.504, 27.072), .0005) && perto(tb$p_corr, .002, .0005))
  teste("η²p = 0,338 (principal) e ω² = 0,020",
        perto(ef(r, "η²p")$valor, .338, .0005) && ef(r, "η²p")$principal && perto(ef(r, "ω²")$valor, .020, .0005))
  teste("HF daria o mesmo p = ,002",
        perto(pf(9.1828, 2 * .805515, 36 * .805515, lower.tail = FALSE), .002, .0005) &&
          anova_medidas_repetidas(stp, med, limiar_epsilon = .60)$tabela$correcao == "Huynh-Feldt")
  ph <- r$posthoc
  teste("Bonferroni N1-N2: dif = -29,158, t(18) = -4,639, p < ,001, d = -0,366",
        perto(c(ph$dif[1], ph$t[1], ph$d_jasp[1]), c(-29.158, -4.639, -.366), .0005) && ph$p_ajustado[1] < .001 && ph$gl[1] == 18)
  teste("Bonferroni N1-N3: dif = -9,316, t = -1,703, p = ,317; N2-N3: 19,842, 2,284, ,104",
        perto(c(ph$dif[2], ph$t[2], ph$p_ajustado[2], ph$dif[3], ph$t[3], ph$p_ajustado[3]), c(-9.316, -1.703, .317, 19.842, 2.284, .104), .0005))
  teste("d do corpo: ,35 na métrica original do par (d_par) e ≈ 1,06 o d bruto das diferenças",
        perto(abs(ph$d_par[1]), .35, .005) && perto(abs(ph$d_z[1]), 1.06, .005))
  teste("médias 440,42 / 469,58 / 449,74 ms", perto(r$descritivas$media, c(440.42, 469.58, 449.74), .005))
  teste("frase APA com F(1,50; 27,07) = 9,18 e p = ,002",
        { a <- anova_medidas_repetidas(stp, med, decimal = ",")$apa[1]; grepl("F(1,50; 27,07) = 9,18", a, fixed = TRUE) && grepl("p = ,002", a, fixed = TRUE) })
  # G-11 (CoVe final): o trecho do Mauchly na frase tinha ", " fixo entre χ² e p com
  # decimal = ","; W = 0,670, χ²(2) = 6,803, p = ,033
  teste("G-11 frase da RM com decimal = ',': 'W = ,67; χ²(2) = 6,80; p = ,033'",
        grepl("após o teste de Mauchly, W = ,67; χ²(2) = 6,80; p = ,033", anova_medidas_repetidas(stp, med, decimal = ",")$apa[1], fixed = TRUE))
}

vi <- jasp("Viagra.csv")
if (is.null(vi)) pular("Viagra") else {
  vi$dose <- factor(vi$dose, levels = 1:3, labels = c("Placebo", "Baixa", "Alta"))
  cat("== Âncoras: Viagra, ANCOVA (Vol. 5, cap. 15)\n")
  r <- ancova(vi, "libido", "dose", "partnerLibido")
  tb <- r$tabela
  teste("covariável F(1, 26) = 4,959, p = ,035; ω²p = 0,117 e η²p = 0,160",
        perto(c(tb$F[1], tb$p[1]), c(4.959, .035), .0005) && tb$gl2[1] == 26 &&
          perto(ef(r, "ω²p", "partnerLibido")$valor, .117, .0005) && perto(ef(r, "η²p", "partnerLibido")$valor, .160, .0005))
  teste("dose F(2, 26) = 4,142, p = ,027; ω²p = 0,173 (principal) e η²p = 0,242",
        perto(c(tb$F[2], tb$p[2]), c(4.142, .027), .0005) && perto(ef(r, "ω²p", "dose")$valor, .173, .0005) &&
          ef(r, "ω²p", "dose")$principal && perto(ef(r, "η²p", "dose")$valor, .242, .0005))
  np <- r$efeito_nao_parcial
  teste("não parciais da dose (não reportar): η² = 0,211, ω² = 0,156",
        perto(np$eta2[np$termo == "dose"], .211, .0005) && perto(np$omega2[np$termo == "dose"], .156, .0005))
  teste("Levene da VD bruta F(2, 27) = 0,71, p = ,500",
        perto(pres(r, "VD bruta")$estatistica, .71, .005) && perto(pres(r, "VD bruta")$p, .500, .0005))
  teste("Levene dos resíduos (painel do JASP) F(2, 27) = 4,618, p = ,019",
        perto(pres(r, "resíduos do modelo")$estatistica, 4.618, .0005) && perto(pres(r, "resíduos do modelo")$p, .019, .0005))
  teste("Shapiro-Wilk dos resíduos W = ,94, p = ,111",
        perto(pres(r, "Normalidade")$estatistica, .94, .005) && perto(pres(r, "Normalidade")$p, .111, .0005))
  md <- r$medias
  teste("médias brutas 3,22/4,88/4,85 e ajustadas 2,926/4,712/5,151",
        perto(md$media_bruta, c(3.22, 4.88, 4.85), .005) && perto(md$media_ajustada, c(2.926, 4.712, 5.151), .0005))
  teste("libido do parceiro média 3,44 no placebo e 2,00 na dose alta",
        perto(md$media_partnerLibido[c(1, 3)], c(3.44, 2.00), .005))
  ph <- r$posthoc
  teste("Placebo-Baixa: dif = -1,786, d = -1,024, p Bonferroni = ,136, IC não ajustado [-3,53; -0,04]",
        perto(c(ph$dif[1], ph$d_qme[1], ph$p_ajustado[1]), c(-1.786, -1.024, .136), .0005) &&
          perto(c(ph$ic_inf_sem_ajuste[1], ph$ic_sup_sem_ajuste[1]), c(-3.53, -.04), .005))
  teste("Placebo-Alta: dif = -2,225, d = -1,276, p = ,031, IC [-3,88; -0,57]; Baixa-Alta: -0,439, d = -0,252, p = 1,000, IC [-2,11; 1,23]",
        perto(c(ph$dif[2], ph$d_qme[2], ph$p_ajustado[2], ph$dif[3], ph$d_qme[3]), c(-2.225, -1.276, .031, -.439, -.252), .0005) &&
          ph$p_ajustado[3] > .9995 && perto(c(ph$ic_inf_sem_ajuste[2:3], ph$ic_sup_sem_ajuste[2:3]), c(-3.88, -2.11, -.57, 1.23), .005))
  inc <- pres(r, "inclinações")
  teste("inclinações: dose × parceiro F(2, 24) = 4,18, p = ,028 e aviso de violação",
        perto(inc$estatistica, 4.18, .005) && inc$gl == "2, 24" && perto(inc$p, .028, .0005) && any(grepl("inclinações violada", r$avisos)))
  ic <- r$inclinacoes
  teste("inclinações por grupo 0,76/0,82/-0,22 e r = ,88/,97/-,17",
        perto(ic$inclinacao, c(.76, .82, -.22), .005) && perto(ic$r, c(.88, .97, -.17), .005))
  teste("frase APA: F(1, 26) = 4.96, p = .035, ω²p = .117 e F(2, 26) = 4.14, p = .027, ω²p = .173",
        grepl("F(1, 26) = 4.96, p = .035, ω²p = .117", r$apa[1], fixed = TRUE) && grepl("F(2, 26) = 4.14, p = .027, ω²p = .173", r$apa[1], fixed = TRUE))
  # G-11 (CoVe final): a frase das inclinações tinha ", " fixo antes do p com decimal = ","
  # ("F(2; 24) = 4,18, p = ,028"); números do livro. Não parciais ,156 e ,211.
  rv <- ancova(vi, "libido", "dose", "partnerLibido", decimal = ",")
  teste("G-11 ANCOVA com decimal = ',': inclinações 'F(2; 24) = 4,18; p = ,028' e decisões 'ω² = ,156; η² = ,211'",
        grepl("(dose × covariável, F(2; 24) = 4,18; p = ,028)", rv$apa[2], fixed = TRUE) && !grepl("4,18, p", rv$apa[2], fixed = TRUE) &&
          grepl("ω² = ,156; η² = ,211", rv$decisoes[2], fixed = TRUE))

  cat("== Âncoras: Viagra, MANOVA (Vol. 5, cap. 17)\n")
  m <- manova_grupos(vi, c("libido", "partnerLibido"), "dose")
  mv <- m$tabela
  lin <- function(s) unlist(mv[mv$estatistica == s, c("valor", "F", "gl1", "gl2", "p")])
  teste("Pillai V = 0,349, F(4, 54) = 2,854, p = ,032 (principal; tol .0005)",
        perto(lin("Pillai"), c(.349, 2.854, 4, 54, .032), .0005) && mv$principal[mv$estatistica == "Pillai"])
  teste("Wilks Λ = 0,661, F(4, 52) = 2,985, p = ,027", perto(lin("Wilks"), c(.661, 2.985, 4, 52, .027), .0005))
  teste("Hotelling-Lawley 0,496, F = 3,101, gl 50,00, p = ,023 (aproximação da tela do JASP; tol .0005)",
        perto(lin("Hotelling-Lawley"), c(.496, 3.101, 4, 50, .023), .0005))
  teste("Roy 0,462, F(2, 27) = 6,235, p = ,006", perto(lin("Roy"), c(.462, 6.235, 2, 27, .006), .0005))
  bx <- pres(m, "Box")
  teste("M de Box χ²(6) = 26,34, p < ,001", perto(bx$estatistica, 26.34, .005) && bx$gl == "6" && bx$p < .001)
  sm <- pres(m, "Shapiro-Wilk multivariado")
  teste("Shapiro-Wilk multivariado W = 0,955, p = ,200 (o teste do JASP; tol .0005)",
        perto(sm$estatistica, .955, .0005) && perto(sm$p, .200, .0005))
  hzr <- pres(m, "Henze-Zirkler")
  teste("Henze-Zirkler dos resíduos HZ = 0,920, p = ,023: rejeita, como conclui o livro",
        perto(hzr$estatistica, .920, .0005) && perto(hzr$p, .023, .0005) && grepl("afastamento", hzr$leitura))
  teste("frase APA do Viagra abre por Box p < .001 com o critério declarado",
        grepl("indicou heterogeneidade das matrizes de covariância, p < .001 (critério p < .001)", m$apa[1], fixed = TRUE))
  teste("ANCOVA do Viagra sem bootstrap: Levene dos resíduos acusa e o aviso sugere bootstrap = TRUE com 1000 reamostragens",
        any(grepl("bootstrap = TRUE com 1000", r$avisos)) && !any(grepl("bootstrap = TRUE com 1000", ancova(vi, "libido", "dose", "partnerLibido", bootstrap = TRUE, n_boot = 200)$avisos)))
  u <- m$univariadas
  teste("univariadas: libido F(2, 27) = 2,416, p = ,108, η²p = ,15; parceiro F = 1,979, p = ,158, η²p = ,13",
        perto(c(u$F, u$p), c(2.416, 1.979, .108, .158), .0005) && perto(u$eta2p, c(.15, .13), .005))
  teste("efeito multivariado ≈ ,17 = V/s", perto(ef(m, "Pillai")$valor, .17, .005))
  teste("um fator: SQ Tipo I (o do JASP) e Tipo III coincidem",
        perto(manova_grupos(vi, c("libido", "partnerLibido"), "dose", tipo_sq = 1)$tabela$F, mv$F, 1e-8))
  teste("frase APA: Box p < .001, Pillai V = .35, F(4, 54) = 2.85, p = .032",
        grepl("p < .001", m$apa[1], fixed = TRUE) && grepl("V = .35, F(4, 54) = 2.85, p = .032", m$apa[1], fixed = TRUE))
}

hr <- jasp("Heart Rate.csv")
if (is.null(hr)) pular("Heart Rate") else {
  cat("== Âncoras: Heart Rate, ANOVA fatorial (Vol. 5, cap. 16)\n")
  r <- anova_fatorial(hr, "Heart Rate", c("Gender", "Group"))
  tb <- r$tabela
  teste("F: Gender 185,980; Group 695,647; interação 7,409 com p = ,007",
        perto(tb$F, c(185.980, 695.647, 7.409), .0005) && perto(tb$p[3], .007, .0005) && all(tb$gl2 == 796))
  teste("ω²p (principal) ,188/,465/,008 e η²p ,189/,466/,009",
        perto(r$efeito$valor[r$efeito$indice == "ω²p"], c(.188, .465, .008), .0005) &&
          perto(r$efeito$valor[r$efeito$indice == "η²p"], c(.189, .466, .009), .0005) && all(r$efeito$principal[r$efeito$indice == "ω²p"]))
  teste("Levene F(3, 796) = 5,56, p < ,001",
        perto(pres(r, "Levene")$estatistica, 5.56, .005) && pres(r, "Levene")$p < .001)
  teste("Shapiro-Wilk dos resíduos W = ,998, p = ,338 e R² = 0,53",
        perto(pres(r, "Normalidade")$estatistica, .998, .0005) && perto(pres(r, "Normalidade")$p, .338, .0005) && perto(r$r2, .53, .005))
  ds <- r$descritivas
  mg <- function(f, l) mean(ds$media[ds[[f]] == l])
  teste("médias: corredores 109,98, sedentários 139,00, mulheres 131,99, homens 116,99",
        perto(c(mg("Group", "Runners"), mg("Group", "Control"), mg("Gender", "Female"), mg("Gender", "Male")), c(109.98, 139.00, 131.99, 116.99), .005))
  teste("células 148,00 / 115,99 / 130,00 / 103,98",
        perto(sort(ds$media, decreasing = TRUE), c(148.00, 130.00, 115.99, 103.98), .005))
  ph <- r$posthoc
  pick <- function(by, cmp) ph[ph$em_by == by & ph$comparacao == cmp, ]
  f1 <- pick("Gender = Female", "Control - Runners"); f2 <- pick("Gender = Male", "Control - Runners")
  f3 <- pick("Group = Control", "Female - Male"); f4 <- pick("Group = Runners", "Female - Male")
  teste("efeito simples mulheres: 32,02, t(796) = 20,57, IC [28,96; 35,07], d = 2,06",
        perto(c(f1$dif, f1$t, f1$ic_inf, f1$ic_sup, f1$d_qme), c(32.02, 20.57, 28.96, 35.07, 2.06), .005) && f1$gl == 796)
  teste("efeito simples homens: 26,03, t = 16,73, IC [22,97; 29,08], d = 1,67",
        perto(c(f2$dif, f2$t, f2$ic_inf, f2$ic_sup, f2$d_qme), c(26.03, 16.73, 22.97, 29.08, 1.67), .005))
  teste("sexo nos sedentários: H - M = -18,00, t = -11,57, IC [-21,05; -14,95], d = -1,16 (sinal invertido: aqui M - H)",
        perto(-c(f3$dif, f3$t, f3$ic_sup, f3$ic_inf, f3$d_qme), c(-18.00, -11.57, -21.05, -14.95, -1.16), .005))
  teste("sexo nos corredores: H - M = -12,01, t = -7,72, IC [-15,06; -8,96], d = -0,77",
        perto(-c(f4$dif, f4$t, f4$ic_sup, f4$ic_inf, f4$d_qme), c(-12.01, -7.72, -15.06, -8.96, -.77), .005))
  teste("frase APA abre pela interação: F(1, 796) = 7.41, p = .007, ω²p = .008",
        grepl("^Interação Gender × Group: F\\(1, 796\\) = 7.41, p = .007, ω²p = .008", r$apa[1]))
}

# ============================================================ 2. Verdade conhecida (simulados)
cat("== Simulados: testes t e não paramétricos\n")
set.seed(2026)
s1 <- data.frame(g = rep(c("A", "B"), each = 300), y = c(rnorm(300, 10.5, 1), rnorm(300, 10, 1)))
r <- teste_t_independente(s1, "y", "g")
teste("t independente: d verdadeiro 0,5 recuperado (|d - 0,5| < ,17, IC contém 0,5)",
      { e <- r$efeito[r$efeito$principal, ]; abs(e$valor - .5) < .17 && e$ic_inf < .5 && e$ic_sup > .5 })
s2 <- data.frame(g = rep(c("A", "B"), c(40, 160)), y = c(rnorm(40, 0, 1), rnorm(160, 0, 4)))
r2 <- teste_t_independente(s2, "y", "g")
teste("variâncias 1 vs 16: Levene acusa, rota clássica vai para Welch e o d principal é o de Welch",
      r2$tabela$escolhido[2] && grepl("Welch", r2$efeito$indice[r2$efeito$principal]))
teste("desbalanceamento NÃO dispara aviso de g (erro do livro não reproduzido); gl < 20 dispara o aviso de amostra pequena",
      !any(grepl("Hedges|desbalance", r2$avisos)) &&
        any(grepl("Amostra pequena", teste_t_independente(s1[c(1:6, 301:306), ], "y", "g")$avisos)))
teste("indice = 'g' torna o g de Hedges principal e menor que o d em módulo",
      { e <- teste_t_independente(s1, "y", "g", indice = "g")$efeito; e$principal[3] && abs(e$valor[3]) < abs(e$valor[1]) })
pre <- rnorm(80, 50, 10); s3 <- data.frame(pre = pre, pos = pre + 5 + rnorm(80, 0, 3))
r3 <- teste_t_pareado(s3, "pre", "pos")
teste("pareado: d na métrica original ≈ 5/10 e d_z bem maior (r alto)",
      { e <- r3$efeito; abs(e$valor[1] - .5) < .1 && e$valor[4] > 1.2 })
s3b <- s3; s3b$pos[1:3] <- NA
teste("pareado: casos incompletos saem por listwise e são contados", teste_t_pareado(s3b, "pre", "pos")$n_excluidos == 3)
teste("uma amostra: diferença e IC centrados em M - mu",
      { u <- teste_t_uma_amostra(s3, "pre", mu = 50)$tabela[1, ]; perto(u$dif_media, mean(s3$pre) - 50, 1e-10) && u$ic_inf < u$dif_media && u$ic_sup > u$dif_media })
teste("teste t com 3 grupos para com mensagem que aponta anova_um_fator()",
      inherits(tryCatch(teste_t_independente(data.frame(y = 1:9, g = rep(1:3, 3)), "y", "g"), error = function(e) e), "error") &&
        grepl("anova_um_fator", tryCatch(teste_t_independente(data.frame(y = 1:9, g = rep(1:3, 3)), "y", "g"), error = function(e) conditionMessage(e))))
mw <- teste_mann_whitney(s1, "y", "g")
teste("Mann-Whitney: A (média maior) tem r_rb positivo no sinal direcional e negativo no sinal do JASP 0.95.1",
      mw$efeito$valor > 0 && teste_mann_whitney(s1, "y", "g", sinal_rrb = "jasp_0951")$efeito$valor < 0)

cat("== Simulados: ANOVA de um fator\n")
set.seed(11)
s4 <- data.frame(g = rep(c("a", "b", "c", "d"), each = 25), y = rnorm(100))
r4 <- anova_um_fator(s4, "y", "g")
teste("médias iguais (semente 11): ômnibus sem significância e pós-hoc não roda", r4$tabela$p[1] > .05 && is.null(r4$posthoc))
teste("exigir_omnibus = FALSE roda o pós-hoc assim mesmo", !is.null(anova_um_fator(s4, "y", "g", exigir_omnibus = FALSE)$posthoc))
set.seed(12)
s5 <- data.frame(g = rep(c("a", "b", "c"), each = 40), y = c(rnorm(40, 0), rnorm(40, 0), rnorm(40, 1.2)))
r5 <- anova_um_fator(s5, "y", "g")
teste("grupo c deslocado 1,2 DP: Tukey acusa a-c e b-c e não a-b",
      { p <- r5$posthoc; p$p_ajustado[p$comparacao == "a - b"] > .05 && all(p$p_ajustado[p$comparacao != "a - b"] < .001) })
set.seed(13)
s6 <- data.frame(g = rep(c("a", "b", "c"), c(20, 40, 80)), y = c(rnorm(20, 0, 1), rnorm(40, 0, 3), rnorm(80, 3, 6)))
r6 <- anova_um_fator(s6, "y", "g")
teste("variâncias muito desiguais: Levene acusa, linha de Welch e Games-Howell", r6$tabela$escolhido[2] && r6$posthoc$metodo[1] == "games_howell")

cat("== Simulados: medidas repetidas e mista\n")
set.seed(21)
n <- 60; suj <- rnorm(n, 0, 2)
s7 <- data.frame(t1 = suj + rnorm(n), t2 = suj + .8 + rnorm(n), t3 = suj + .8 + rnorm(n))
r7 <- anova_medidas_repetidas(s7, c("t1", "t2", "t3"))
teste("simetria composta (esfericidade verdadeira): Mauchly sem evidência e nenhuma correção",
      pres(r7, "Mauchly")$p > .05 && r7$tabela$correcao == "nenhuma")
teste("ε de GG igual à fórmula de Box a partir da covariância dos contrastes ortonormais",
      { C <- contr.poly(3); S <- t(C) %*% cov(as.matrix(s7)) %*% C
        e <- sum(diag(S))^2 / (2 * sum(S^2))
        perto(as.numeric(sub(".*GG = ([0-9.]+),.*", "\\1", pres(r7, "Mauchly")$leitura)), e, .0005) })
set.seed(22)
s8 <- data.frame(t1 = rnorm(40, 0, 1), t2 = rnorm(40, 0, 1), t3 = rnorm(40, .5, 1)); s8$t2 <- s8$t1 + rnorm(40, 0, .2)
r8 <- anova_medidas_repetidas(s8, c("t1", "t2", "t3"))
teste("variâncias das diferenças muito desiguais: Mauchly acusa e a correção é aplicada (gl fracionários)",
      pres(r8, "Mauchly")$p < .05 && r8$tabela$correcao != "nenhuma" && r8$tabela$gl1_corr < 2)
teste("esfericidade = 'epsilon' corrige mesmo sem Mauchly significativo",
      anova_medidas_repetidas(s7, c("t1", "t2", "t3"), esfericidade = "epsilon")$tabela$correcao != "nenhuma")
set.seed(23)
nm <- 50; gm <- rep(c("ctl", "int"), each = nm / 2); sj <- rnorm(nm, 0, 1.5)
s9 <- data.frame(grupo = gm, pre = sj + rnorm(nm), pos = sj + ifelse(gm == "int", 1.2, 0) + rnorm(nm), seg = sj + ifelse(gm == "int", 1, 0) + rnorm(nm))
r9 <- anova_medidas_repetidas(s9, c("pre", "pos", "seg"), intra = list(tempo = c("pre", "pos", "seg")), entre = "grupo")
teste("mista: interação grupo × tempo significativa e efeitos simples com Bonferroni sobre o conjunto",
      { it <- r9$tabela[r9$tabela$termo == "grupo × tempo", ]; it$p_corr < .01 && !is.null(r9$posthoc) &&
          perto(r9$posthoc$p_ajustado, pmin(1, r9$posthoc$p_sem_ajuste * nrow(r9$posthoc)), 1e-12) })
teste("mista: Levene por tempo e M de Box presentes nos pressupostos",
      sum(grepl("Levene", r9$pressupostos$checagem)) == 3 && any(grepl("Box", r9$pressupostos$checagem)))
teste("mista: η²p é o principal", all(r9$efeito$principal == (r9$efeito$indice == "η²p")))

cat("== Simulados: fatorial, ANCOVA, MANOVA, Kruskal-Wallis, Friedman\n")
set.seed(31)
s10 <- data.frame(A = sample(c("a1", "a2"), 150, TRUE, prob = c(.3, .7)), B = sample(c("b1", "b2", "b3"), 150, TRUE))
s10$y <- rnorm(150) + (s10$A == "a2") * .8
r10 <- anova_fatorial(s10, "y", c("A", "B"))
teste("fatorial desbalanceado: efeito de A detectado, B e A × B sem efeito; sem efeitos simples",
      r10$tabela$p[1] < .001 && r10$tabela$p[2] > .05 && r10$tabela$p[3] > .05 && is.null(r10$posthoc))
set.seed(32)
x <- rnorm(90); gg <- rep(c("c", "t1", "t2"), each = 30)
s11 <- data.frame(g = gg, x = x, y = 2 + .8 * x + ifelse(gg == "t2", 1, 0) + rnorm(90, 0, .7))
r11 <- ancova(s11, "y", "g", "x", grupos_aleatorizados = TRUE)
teste("ANCOVA: diferença ajustada t2 - c ≈ 1 (verdade 1; tol .35) e inclinações homogêneas",
      { p <- r11$posthoc; abs(-p$dif[p$comparacao == "c - t2"] - 1) < .35 && pres(r11, "inclinações")$p > .05 })
teste("ANCOVA com grupos aleatorizados declarados não emite o aviso de Miller & Chapman", !any(grepl("Miller", r11$avisos)))
s12 <- s11; s12$y <- s12$y + ifelse(s12$g == "t1", -1.6 * s12$x, 0)
teste("ANCOVA com inclinações heterogêneas: teste acusa e aviso aparece",
      { a <- ancova(s12, "y", "g", "x"); pres(a, "inclinações")$p < .001 && any(grepl("inclinações violada", a$avisos)) })
teste("bootstrap da ANCOVA: mesma semente, mesmo IC; IC contém a estimativa",
      { b1 <- ancova(s11, "y", "g", "x", bootstrap = TRUE, n_boot = 300)$posthoc
        b2 <- ancova(s11, "y", "g", "x", bootstrap = TRUE, n_boot = 300)$posthoc
        identical(b1$ic_boot_inf, b2$ic_boot_inf) && all(b1$ic_boot_inf < b1$dif & b1$ic_boot_sup > b1$dif) })
teste("bootstrap não altera o gerador de números aleatórios do usuário",
      { set.seed(5); a1 <- runif(1); set.seed(5); invisible(ancova(s11, "y", "g", "x", bootstrap = TRUE, n_boot = 50)); a2 <- runif(1); a1 == a2 })
set.seed(33)
s13 <- data.frame(g = rep(c("a", "b", "c"), each = 30))
s13$y1 <- rnorm(90) + (s13$g == "c") * .9; s13$y2 <- rnorm(90)
r13 <- manova_grupos(s13, c("y1", "y2"), "g")
teste("MANOVA: efeito em y1 detectado pelo Pillai e pela univariada de y1, não pela de y2",
      r13$tabela$p[r13$tabela$estatistica == "Pillai"] < .01 && r13$univariadas$p[1] < .01 && r13$univariadas$p[2] > .05)
set.seed(34)
s14 <- data.frame(g = rep(c("a", "b", "c"), each = 30), y = c(rexp(30), rexp(30), rexp(30) + 1.5))
r14 <- teste_kruskal_wallis(s14, "y", "g")
teste("Kruskal-Wallis: grupo c deslocado; Dunn (Holm) acusa a-c e b-c, não a-b",
      { p <- r14$posthoc; p$p_ajustado[p$comparacao == "a - b"] > .05 && all(p$p_ajustado[p$comparacao != "a - b"] < .01) })
teste("Kruskal-Wallis: Holm nunca dá p maior que Bonferroni", all(r14$posthoc$p_holm <= r14$posthoc$p_bonferroni + 1e-15))
set.seed(35)
s15 <- data.frame(t1 = rnorm(30), t2 = rnorm(30) + 1.3, t3 = rnorm(30) + 1.3)
r15 <- teste_friedman(s15, c("t1", "t2", "t3"))
teste("Friedman: t1 diferente de t2 e t3; t2 e t3 iguais (Conover, Holm)",
      { p <- r15$posthoc; p$p_ajustado[p$comparacao == "t2 - t3"] > .05 && all(p$p_ajustado[p$comparacao != "t2 - t3"] < .05) })

# ============================================================ 3. Equivalência com outra implementação
cat("== Equivalências\n")
suppressPackageStartupMessages(requireNamespace("effectsize", quietly = TRUE))
op <- options(contrasts = c("contr.sum", "contr.poly"))
teste("Levene (média e mediana) = car::leveneTest",
      { a <- .ee_dif_levene(s6$y, factor(s6$g), "media"); b <- car::leveneTest(y ~ factor(g), s6, center = mean)
        a2 <- .ee_dif_levene(s6$y, factor(s6$g), "mediana"); b2 <- car::leveneTest(y ~ factor(g), s6, center = median)
        perto(a[["F"]], b[1, "F value"], 1e-10) && perto(a2[["F"]], b2[1, "F value"], 1e-10) })
teste("d de Cohen e IC (Student) = effectsize::cohens_d; g = effectsize::hedges_g",
      { e <- r$efeito; ce <- effectsize::cohens_d(y ~ g, data = s1); he <- effectsize::hedges_g(y ~ g, data = s1)
        perto(c(e$valor[1], e$ic_inf[1], e$ic_sup[1]), c(ce$Cohens_d, ce$CI_low, ce$CI_high), 1e-6) &&
          perto(c(e$valor[3], e$ic_inf[3], e$ic_sup[3]), c(he$Hedges_g, he$CI_low, he$CI_high), 1e-6) })
teste("d de Welch (média das variâncias) = effectsize::cohens_d(pooled_sd = FALSE) no ponto",
      perto(r2$efeito$valor[2], effectsize::cohens_d(y ~ g, data = s2, pooled_sd = FALSE)$Cohens_d, 1e-10))
teste("d na métrica original, d_av, d_rm e d_z do pareado = fórmulas fechadas; IC do principal = effectsize 'av'",
      { e <- r3$efeito; s1x <- sd(s3$pre); s2x <- sd(s3$pos); dd <- s3$pos - s3$pre; rr <- cor(s3$pre, s3$pos)
        es <- effectsize::repeated_measures_d(s3$pos, s3$pre, method = "av", adjust = FALSE)
        perto(e$valor[1:4], c(mean(dd) / sqrt((s1x^2 + s2x^2) / 2), mean(dd) / ((s1x + s2x) / 2), mean(dd) / sd(dd) * sqrt(2 * (1 - rr)), mean(dd) / sd(dd)), 1e-10) &&
          perto(c(e$ic_inf[1], e$ic_sup[1]), c(es$CI_low, es$CI_high), 1e-10) })
teste("r_rb de Mann-Whitney e IC = effectsize::rank_biserial (mesmo sinal no modo direcional)",
      { rb <- effectsize::rank_biserial(y ~ g, data = s1); perto(c(mw$efeito$valor, mw$efeito$ic_inf, mw$efeito$ic_sup), c(rb$r_rank_biserial, rb$CI_low, rb$CI_high), 1e-6) })
teste("r_rb pareada e IC = effectsize::rank_biserial(paired = TRUE)",
      { w <- teste_wilcoxon_pareado(s3, "pre", "pos"); rb <- effectsize::rank_biserial(s3$pos, s3$pre, paired = TRUE)
        perto(c(w$efeito$valor, w$efeito$ic_inf, w$efeito$ic_sup), c(rb$r_rank_biserial, rb$CI_low, rb$CI_high), 1e-6) })
teste("r_rb pareada exclui diferenças nulas (fórmula do JASP: 2V/maxw - 1, maxw sobre as não nulas)",
      { z <- data.frame(a = c(1:10, 5, 5), b = c(1:10 + c(2, -1, 3, 4, -2, 5, 6, 1, 2, 3), 5, 5))
        w <- teste_wilcoxon_pareado(z, "a", "b"); d <- z$b - z$a; d <- d[d != 0]; rk <- rank(abs(d))
        perto(w$efeito$valor, 2 * sum(rk[d > 0]) / (length(d) * (length(d) + 1) / 2) - 1, 1e-12) && w$tabela$n_dif_nao_nulas == 10 })
teste("Mann-Whitney com empates segue a regra do JASP 0.95.1 (normal com correção), não o p exato condicional do R 4.6",
      { set.seed(41); z <- data.frame(g = rep(c("a", "b"), each = 15), y = c(sample(1:6, 15, TRUE), sample(3:8, 15, TRUE)))
        m1 <- teste_mann_whitney(z, "y", "g"); w <- suppressWarnings(wilcox.test(y ~ g, z, exact = FALSE, correct = TRUE))
        perto(m1$tabela$p, w$p.value, 1e-12) && !m1$tabela$exato })
teste("Wilcoxon pareado com zeros: V e p como no R 4.5/JASP (zeros fora antes de ranquear)",
      { z <- data.frame(a = c(1:10, 5, 5), b = c(1:10 + c(2, -1, 3, 4, -2, 5, 6, 1, 2, 3), 5, 5))
        w <- teste_wilcoxon_pareado(z, "a", "b"); ref <- suppressWarnings(wilcox.test((z$b - z$a)[z$b != z$a], exact = FALSE, correct = TRUE))
        w$tabela$V == 49.5 && perto(w$tabela$p, ref$p.value, 1e-12) })
teste("ANOVA de um fator = car::Anova(type = 3) e Welch = oneway.test",
      { a <- car::Anova(aov(y ~ g, s6), type = 3); perto(r6$tabela$F[1], a["g", "F value"], 1e-10) &&
          perto(r6$tabela$F[2], unname(oneway.test(y ~ g, s6)$statistic), 1e-10) })
teste("Brown-Forsythe F* = onewaytests::bf.test",
      { b <- onewaytests::bf.test(y ~ g, data = transform(s6, g = factor(g)), verbose = FALSE)
        perto(c(r6$tabela$F[3], r6$tabela$gl2[3], r6$tabela$p[3]), c(b$statistic, b$parameter[2], b$p.value), 1e-8) })
teste("Games-Howell = rstatix::games_howell_test (diferença, IC e p; o rstatix subtrai 2º - 1º)",
      { gh <- as.data.frame(rstatix::games_howell_test(s6, y ~ g)); p <- r6$posthoc
        perto(-p$dif, gh$estimate, 1e-8) && perto(p$p_ajustado, gh$p.adj, 1e-3) && perto(sort(-c(p$ic_inf, p$ic_sup)), sort(c(gh$conf.low, gh$conf.high)), 1e-8) })
teste("fatorial desbalanceado Tipo III = car::Anova(lm, type = 3); Tipo II = car::Anova(type = 2)",
      { a3 <- car::Anova(lm(y ~ A * B, s10), type = 3); a2 <- car::Anova(lm(y ~ A * B, s10), type = 2)
        perto(r10$tabela$F, a3[c("A", "B", "A:B"), "F value"], 1e-10) &&
          perto(anova_fatorial(s10, "y", c("A", "B"), tipo_sq = 2)$tabela$F, a2[c("A", "B", "A:B"), "F value"], 1e-10) })
teste("ω²p e η²p da fatorial = effectsize sobre a tabela car (bilateral)",
      { a3 <- car::Anova(lm(y ~ A * B, s10), type = 3); o <- suppressMessages(effectsize::omega_squared(a3, partial = TRUE, alternative = "two.sided"))
        perto(r10$efeito$valor[r10$efeito$indice == "ω²p"], o$Omega2_partial, 1e-10) })
teste("medidas repetidas = afex::aov_car (F, gl corrigidos por GG, p)",
      { lg <- reshape(transform(s8, id = seq_len(nrow(s8))), direction = "long", varying = c("t1", "t2", "t3"), v.names = "y", timevar = "t", times = c("t1", "t2", "t3"), idvar = "id")
        a <- suppressMessages(afex::aov_car(y ~ Error(id / t), data = lg, anova_table = list(correction = "GG")))
        at <- a$anova_table; tb <- r8$tabela
        tb$correcao == "Greenhouse-Geisser" && perto(c(tb$F, tb$gl1_corr, tb$gl2_corr, tb$p_corr), c(at$F, at$`num Df`, at$`den Df`, at$`Pr(>F)`), 1e-8) })
teste("ω² de medidas repetidas (fórmula do JASP) = effectsize::omega_squared(afex) no Stroop simulado",
      { lg <- reshape(transform(s7, id = seq_len(nrow(s7))), direction = "long", varying = c("t1", "t2", "t3"), v.names = "y", timevar = "t", times = c("t1", "t2", "t3"), idvar = "id")
        a <- suppressMessages(afex::aov_car(y ~ Error(id / t), data = lg))
        perto(ef(r7, "ω²")$valor, effectsize::omega_squared(a, alternative = "two.sided")[[2]], 1e-8) })
teste("MANOVA com 2 fatores: Tipo III = car::Manova(type = 3); Tipo I = summary.manova",
      { s16 <- s13; s16$h <- rep(c("x", "y", "y"), 30)
        m3 <- manova_grupos(s16, c("y1", "y2"), c("g", "h")); m1 <- manova_grupos(s16, c("y1", "y2"), c("g", "h"), tipo_sq = 1)
        cm <- car::Manova(lm(cbind(y1, y2) ~ g * h, s16), type = 3)
        pc <- summary(cm)$multivariate.tests$g$SSPH
        sm <- summary(manova(cbind(y1, y2) ~ g * h, s16), test = "Pillai")$stats
        V3 <- m3$tabela$valor[m3$tabela$termo == "g" & m3$tabela$estatistica == "Pillai"]
        V1 <- m1$tabela$valor[m1$tabela$termo == "h" & m1$tabela$estatistica == "Pillai"]
        ev <- Re(eigen(solve(cm$SSPE) %*% cm$SSP$g)$values)
        perto(V3, sum(ev / (1 + ev)), 1e-10) && perto(V1, sm["h", "Pillai"], 1e-10) })
teste("M de Box = heplots::boxM", perto(pres(r13, "Box")$estatistica, unname(heplots::boxM(as.matrix(s13[c("y1", "y2")]), factor(s13$g))$statistic), 1e-10))
teste("Kruskal-Wallis: H = kruskal.test; ε²_R = effectsize::rank_epsilon_squared; η²_H = (H - k + 1)/(N - k)",
      { k <- kruskal.test(y ~ g, s14); e <- effectsize::rank_epsilon_squared(y ~ g, data = s14)
        perto(r14$tabela$H, unname(k$statistic), 1e-10) && perto(r14$efeito$valor[1], e$rank_epsilon_squared, 1e-10) &&
          perto(r14$efeito$valor[2], (k$statistic - 2) / (90 - 3), 1e-10) })
teste("Dunn = rstatix::dunn_test (|z| e p Holm; correção de empates)",
      { dn <- as.data.frame(rstatix::dunn_test(s14, y ~ g, p.adjust.method = "holm"))
        perto(abs(r14$posthoc$z), abs(dn$statistic), 1e-8) && perto(r14$posthoc$p_holm, dn$p.adj, 1e-8) })
teste("Friedman: χ² = friedman.test; W = rstatix::friedman_effsize",
      { ft <- friedman.test(as.matrix(s15)); lg <- data.frame(id = rep(1:30, 3), t = rep(c("t1", "t2", "t3"), each = 30), y = unlist(s15))
        fe <- rstatix::friedman_effsize(lg, y ~ t | id)
        perto(r15$tabela$qui2, unname(ft$statistic), 1e-10) && perto(r15$efeito$valor, fe$effsize, 1e-8) })
if (requireNamespace("PMCMRplus", quietly = TRUE)) {
  teste("Conover = PMCMRplus::frdAllPairsConoverTest (t e p sem ajuste)",
        { M <- as.matrix(s15); rownames(M) <- seq_len(nrow(M)); cv <- PMCMRplus::frdAllPairsConoverTest(M, p.adjust.method = "none")
          perto(abs(r15$posthoc$t_conover), abs(c(cv$statistic[1, 1], cv$statistic[2, 1], cv$statistic[2, 2])), 1e-8) &&
            perto(r15$posthoc$p_sem_ajuste, c(cv$p.value[1, 1], cv$p.value[2, 1], cv$p.value[2, 2]), 1e-8) })
} else cat("PULADO  Conover vs PMCMRplus (pacote ausente; não é dependência da skill)\n")
options(op)

cat("== Regressões da revisão CoVe v1.1\n")
for (tp in c(1, 3)) {
  set.seed(1); n <- 40
  dm <- data.frame(g = rep(c("A", "B"), each = 20), g3 = rep(c("A", "B", "C"), length.out = n), y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), y4 = rnorm(n))
  dm$y1[dm$g == "B"] <- dm$y1[dm$g == "B"] + .6
  ref <- function(f, dd) { o <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(o))
    do.call(rbind, lapply(c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"), function(ts) {
      if (tp == 1) { x <- summary(manova(f, dd), test = ts)$stats; unname(x[1, 2:6]) }
      else { cm <- car::Manova(lm(f, dd), type = 3); x <- capture.output(sm <- summary(cm, multivariate = TRUE))
        st <- sm$multivariate.tests[[2]]; eig <- Re(eigen(solve(st$SSPE) %*% st$SSPH)$values)
        unlist(.ee_dif_mv_estat(pmax(eig, 0), st$df, st$df.residual)[match(ts, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")), 2:6]) } })) }
  m2 <- manova_grupos(dm, c("y1", "y2", "y3"), "g", tipo_sq = tp)
  teste(paste0("MANOVA 2 grupos × 3 VDs (tipo ", tp, "): valor, F, gl e p das 4 estatísticas = ", if (tp == 1) "summary.manova" else "car::Manova"),
        perto(as.matrix(m2$tabela[, c("valor", "F", "gl1", "gl2", "p")]), ref(cbind(y1, y2, y3) ~ g, dm), 1e-8) && all(m2$tabela$gl1 == 3))
  m4 <- manova_grupos(dm, c("y1", "y2", "y3", "y4"), "g3", tipo_sq = tp)
  teste(paste0("MANOVA 3 grupos × 4 VDs (tipo ", tp, "): valor, F, gl e p das 4 estatísticas = ", if (tp == 1) "summary.manova" else "car::Manova"),
        perto(as.matrix(m4$tabela[, c("valor", "F", "gl1", "gl2", "p")]), ref(cbind(y1, y2, y3, y4) ~ g3, dm), 1e-8) && m4$tabela$gl1[1] == 8)
}
teste("MANOVA 2 grupos × 3 VDs: Pillai = car::Manova (F(3, 36) = 0.79, p = .508) e a frase não declara efeito",
      { m2 <- manova_grupos(dm, c("y1", "y2", "y3"), "g"); perto(m2$tabela$p[1], .5077, .0001) && grepl("não foi estatisticamente significativo", m2$apa[1]) })
teste("nivel = .90: todos os IC da frase rotulados 90% (t, pareado, um fator, Mann-Whitney, Kruskal-Wallis)",
      { a <- c(teste_t_independente(s1, "y", "g", nivel = .90)$apa, teste_t_pareado(s3, "pre", "pos", nivel = .90)$apa,
               anova_um_fator(s5, "y", "g", nivel = .90)$apa[1], teste_mann_whitney(s1, "y", "g", nivel = .90)$apa,
               teste_kruskal_wallis(s14, "y", "g", nivel = .90)$apa[1])
        !any(grepl("95%", a)) && all(grepl("IC 90%", a)) })
x1w <- c(10, 12, 15, 11, 9, 14, 13, 8, 16, 12.5, 7, 10)
dw <- data.frame(pre = x1w, pos = x1w + c(0, 0, 1.1, 2.3, -0.4, 3.7, 0, 4.2, 5.6, -1.9, 2.8, 6.1))
teste("Wilcoxon pareado com zeros e sem empates: p normal com correção = R 4.5.1/JASP 0.95.1 (,0330), não exato",
      { w <- teste_wilcoxon_pareado(dw, "pre", "pos"); perto(w$tabela$p, .03297, .00001) && !w$tabela$exato && w$tabela$V == 41 })
teste("Wilcoxon pareado exato = TRUE = wilcox.test(exact = TRUE) do R instalado (no R 4.6: V = 62, zeros nos postos)",
      { w <- teste_wilcoxon_pareado(dw, "pre", "pos", exato = TRUE); ref <- suppressWarnings(wilcox.test(dw$pos, dw$pre, paired = TRUE, exact = TRUE))
        perto(c(w$tabela$V, w$tabela$p), c(unname(ref$statistic), ref$p.value), 1e-12) })
teste("Wilcoxon de uma amostra com valores iguais a mu: p = R 4.5.1 (,236)",
      perto(teste_t_uma_amostra(data.frame(v = c(5, 5, 6.2, 3.1, 7.4, 8.9, 5, 2.2, 9.6, 6.6, 4.3, 7.7)), "v", mu = 5)$tabela$p[2], .2361, .0001))
set.seed(8); Af <- sample(c("a1", "a2", "a3"), 90, TRUE, prob = c(.5, .3, .2)); Bf <- sample(c("b2", "b1"), 90, TRUE)
df8 <- data.frame(y = 10 + (Af == "a2") * 2 + (Af == "a3") * 1 + (Bf == "b1") * 1.5 + (Af == "a3" & Bf == "b1") * 4 + rnorm(90, 0, 2), A = Af, B = Bf)
teste("fatorial com interação: posthoc = 'nenhum' não calcula efeitos simples nem os põe na frase",
      { f0 <- anova_fatorial(df8, "y", c("A", "B"), posthoc = "nenhum"); is.null(f0$posthoc) && length(f0$apa) == 1 })
teste("fatorial com interação: posthoc = 'holm' aplica Holm ao conjunto de efeitos simples",
      { fh <- anova_fatorial(df8, "y", c("A", "B"), posthoc = "holm")$posthoc; perto(fh$p_ajustado, p.adjust(fh$p_sem_ajuste, "holm"), 1e-12) })
teste("unilateral: frase com '∞' e 'unilateral', sem 'Inf' (t independente, pareado, uma amostra)",
      { a <- c(teste_t_independente(s1, "y", "g", alternativa = "maior")$apa, teste_t_pareado(s3, "pre", "pos", alternativa = "maior")$apa,
               teste_t_uma_amostra(s3, "pre", mu = 50, alternativa = "menor")$apa)
        !any(grepl("Inf", a)) && all(grepl("unilateral", a)) && all(grepl("\u221e", enc2native(a), useBytes = TRUE) | grepl(rawToChar(as.raw(c(0xe2, 0x88, 0x9e))), a, fixed = TRUE, useBytes = TRUE)) })
teste("unilateral: IC do r_rb unilateral = effectsize::rank_biserial(alternative = 'less'/'greater')",
      { a <- teste_mann_whitney(s1, "y", "g", alternativa = "menor")$efeito; ra <- effectsize::rank_biserial(y ~ g, data = s1, alternative = "less")
        b <- teste_wilcoxon_pareado(s3, "pre", "pos", alternativa = "maior")$efeito; rb <- effectsize::rank_biserial(s3$pos, s3$pre, paired = TRUE, alternative = "greater")
        perto(c(a$ic_inf, a$ic_sup), c(ra$CI_low, ra$CI_high), 1e-6) && perto(c(b$ic_inf, b$ic_sup), c(rb$CI_low, rb$CI_high), 1e-6) })
teste("privacidade: Mann-Whitney e Kruskal-Wallis param com grupo de 1 caso",
      inherits(tryCatch(teste_mann_whitney(data.frame(y = c(9.3, 1:8 + .1), g = c("solo", rep("resto", 8))), "y", "g"), error = function(e) e), "error") &&
        inherits(tryCatch(teste_kruskal_wallis(data.frame(y = c(3.3, 1:12), g = c("u", rep(c("a", "b"), 6))), "y", "g"), error = function(e) e), "error"))
teste("privacidade: descritivas com n < 5 saem sem estatísticas no print (Kruskal-Wallis e mista)",
      { set.seed(3); k3 <- data.frame(y = c(91.7, 88.2, 77.3, rnorm(20)), g = c(rep("raro", 3), rep(c("a", "b"), 10)))
        o1 <- capture.output(print(teste_kruskal_wallis(k3, "y", "g", n_boot = 100)))
        mx <- data.frame(g = c(rep("raro", 3), rep("comum", 20)), t1 = c(91.7, 88.2, 77.3, rnorm(20)), t2 = c(90.1, 85.5, 70.2, rnorm(20)), t3 = c(95.3, 80.4, 72.9, rnorm(20)))
        o2 <- capture.output(print(anova_medidas_repetidas(mx, c("t1", "t2", "t3"), entre = "g", posthoc = "nenhum")))
        any(grepl("n < 5", o1)) && any(grepl("n < 5", o2)) && !any(grepl("85.57|85.6 |88.2|91.7", c(o1, o2))) })
teste("ANCOVA com 2 covariáveis: inclinação por grupo é a parcial (lm com as duas covariáveis) e a leitura diz 'parciais'",
      { set.seed(9); n <- 60; g <- rep(c("a", "b", "c"), each = 20); x <- rnorm(n); x2 <- .6 * x + rnorm(n)
        yy <- 1 + x + ifelse(g == "a", 2, ifelse(g == "b", .5, -1)) * x2 + rnorm(n); da <- data.frame(y = yy, g = g, x = x, x2 = x2)
        a <- ancova(da, "y", "g", c("x", "x2")); ic <- a$inclinacoes
        ref <- unname(coef(lm(y ~ x + x2, da[da$g == "a", ]))["x2"])
        perto(ic$inclinacao[ic$covariavel == "x2" & ic$grupo == "a"], ref, 1e-10) && all(ic$tipo == "parcial (demais covariáveis no modelo)") &&
          grepl("parciais", pres(a, "g × x2")$leitura) })
teste("M de Box com .001 <= p < .05: pressupostos e frase usam o mesmo critério (p < .001) sem contradição",
      { set.seed(16); n <- 60; db <- data.frame(A = sample(c("a1", "a2", "a3"), n, TRUE), B = sample(c("b1", "b2", "b3"), n, TRUE), y1 = rnorm(n), y2 = rnorm(n))
        db$y1 <- db$y1 + (db$A == "a2") * .8
        mb <- suppressWarnings(manova_grupos(db, c("y1", "y2"), c("A", "B"))); bx <- pres(mb, "Box")
        bx$p > .001 && bx$p < .05 && grepl("sem evidência de heterogeneidade pelo critério p < .001", bx$leitura) && grepl("não indicou heterogeneidade", mb$apa[1]) })
teste("pós-hoc de um fator com grupos desbalanceados: a frase usa d_qme (não troca para g por desbalanceamento)",
      { ph <- r6$posthoc; s <- ph[ph$p_ajustado < .05, ][1, ]; grepl(paste0("d = ", formatC(round(s$d_qme, 2), format = "f", digits = 2)), r6$apa[2], fixed = TRUE) })
teste("g par a par fica NA (não -Inf) quando os dois grupos não variam",
      is.na(.ee_dif_hedges_par(c(1, 1, 1), c(2, 2, 2))))
if (requireNamespace("MVN", quietly = TRUE)) {
  teste("Henze-Zirkler = MVN::hz (estatística e p)",
        { X <- as.data.frame(residuals(lm(cbind(y1, y2, y3) ~ g3, dm))); h <- MVN::hz(X); hz <- .ee_dif_hz(X)
          perto(hz, c(h$Statistic, h$p.value), 1e-10) })
} else cat("PULADO  Henze-Zirkler vs MVN (pacote ausente; não é dependência da skill)\n")

cat("== Regressões da revisão CoVe final (achados H-02, G-01, G-05, G-07, G-11)\n")
# H-02: grupo c com 2 casos (91,234 e 93,456). M = 92,345 e DP = 1,571 devolvem os dois
# escores (M ± DP/√2); a tabela de descritivas do print já os apagava, a frase não. No t, M e
# DP do outro grupo também saem da frase: com eles, a diferença e o t de Welch devolvem
# M_c = M_b - dif e DP_c² = 2(EP² - DP_b²/20), EP = dif/t (conferido: 92,345 e 1,571).
set.seed(7)
dp2 <- data.frame(g = c(rep("b", 20), "c", "c"), x = c(round(rnorm(20, 50, 8), 3), 91.234, 93.456))
alvo2 <- "91\\.23|93\\.4[56]|92\\.3[45]|DP = 1\\.57|SD = 1\\.57"
sai <- function(r) c(capture.output(print(r)), r$apa)
teste("H-02 t independente com grupo de 2 casos: nem o print nem $apa trazem M, DP ou escores; frase sem M e DP nos dois grupos; $avisos registra",
      { r2a <- teste_t_independente(dp2, "x", "g"); r2b <- teste_t_independente(dp2, "x", "g", idioma = "en", alternativa = "menor")
        !any(grepl(alvo2, c(sai(r2a), sai(r2b)))) && !any(grepl("M = |DP = |SD = ", c(r2a$apa, r2b$apa))) &&
          grepl("(n < 5, descritivas omitidas nos dois grupos)", r2a$apa, fixed = TRUE) &&
          grepl("(n < 5, descriptives omitted for both groups)", r2b$apa, fixed = TRUE) && any(grepl("n < 5 (c)", r2a$avisos, fixed = TRUE)) })
teste("H-02 Mann-Whitney com grupo de 2 casos: sem Md do grupo na frase nem no print; $avisos registra",
      { m2a <- teste_mann_whitney(dp2, "x", "g")
        !any(grepl(alvo2, sai(m2a))) && grepl("c (n < 5, descritivas omitidas)", m2a$apa, fixed = TRUE) && grepl("b (Md = ", m2a$apa, fixed = TRUE) &&
          any(grepl("sem mediana", m2a$avisos)) })
teste("H-02 ANCOVA com grupo de 2 casos: média ajustada do grupo fora da frase",
      { set.seed(8); da2 <- data.frame(g = c(rep(c("a", "b"), each = 15), "c", "c"), cv = rnorm(32)); da2$x <- 50 + 3 * da2$cv + rnorm(32, 0, 5)
        a2 <- ancova(da2, "x", "g", "cv", posthoc = "nenhum")
        grepl("c (n < 5, descritivas omitidas)", a2$apa[1], fixed = TRUE) && !grepl("c = ", a2$apa[1], fixed = TRUE) && any(grepl("sem média ajustada", a2$avisos)) })
# pareado com 4 pares: M de pre = 10,90 (DP 1,16), de pos = 12,80 (DP 1,40), medianas 10,75 e
# 13,05 saíam na frase; a tabela de descritivas do print já as apagava
p4 <- data.frame(pre = c(10.1, 12.3, 9.8, 11.4), pos = c(12.2, 13.9, 11.1, 14.0))
teste("H-02 pareado, Wilcoxon pareado e uma amostra com n = 4: frase sem M, DP nem Md (10.90, 12.80, 1.16, 1.40, 10.75, 13.05)",
      { tp <- teste_t_pareado(p4, "pre", "pos"); wp <- teste_wilcoxon_pareado(p4, "pre", "pos"); ua <- teste_t_uma_amostra(p4, "pre", mu = 10)
        fr <- c(tp$apa, wp$apa, ua$apa)
        !any(grepl("M = |Md = |DP = ", c(sai(tp), fr))) && !any(grepl("10\\.90|12\\.80|1\\.16|1\\.40|10\\.75|13\\.05", fr)) &&
          all(grepl("n < 5, descritivas omitidas", fr, fixed = TRUE)) &&
          any(grepl("n < 5", tp$avisos)) && sum(!grepl("^Não significativo", wp$avisos)) == 1 && any(grepl("permite recuperar M", ua$avisos)) &&
          !any(grepl("10\\.90|12\\.80|1\\.16|1\\.40|10\\.75|13\\.05", c(tp$avisos, wp$avisos, ua$avisos))) })
teste("H-02 regra não muda frase com grupos de 5 ou mais (s1, 300 por grupo: M e DP presentes, sem aviso)",
      grepl("(M = ", teste_t_independente(s1, "y", "g")$apa, fixed = TRUE) && !length(grep("n < 5", teste_t_independente(s1, "y", "g")$avisos)))

# G-01 (reprodução g_rep_G01_welch_omega.R): rota Welch com F de Welch significativo e ω²
# clássico = 0. Welch F(2; 29,99) = 4,07, p = ,027; clássica F(2; 50) = 0,75, p = ,479.
# ω² pelo F de Welch = (4,07 - 1)·2/(4,07·2 + 29,99 + 1) = ,157, o mesmo que
# effectsize::omega_squared(oneway.test(var.equal = FALSE)).
set.seed(1); dw1 <- data.frame(gr = rep(1:3, c(8, 15, 30)), y = c(rnorm(8, 10, 1), rnorm(15, 11, 3), rnorm(30, 12, 6)))
rw1 <- anova_um_fator(dw1, "y", "gr")
teste("G-01 rota Welch: aviso de contradição (Welch significativo, ω² clássico com IC que inclui zero) e nota em $decisoes",
      rw1$tabela$escolhido[2] && perto(rw1$tabela$p[2], .027, .0005) && rw1$efeito$valor[1] == 0 &&
        any(grepl("Contradição aparente", rw1$avisos)) && any(grepl("F(2, 50) = 0.75, p = .479", rw1$avisos, fixed = TRUE)) &&
        any(grepl("Ômnibus de Welch: o ω² principal vem da ANOVA clássica", rw1$decisoes, fixed = TRUE)))
teste("G-01 ω² aproximado pelo F de Welch = (F - 1)gl1/(F·gl1 + gl2 + 1) = effectsize::omega_squared(oneway.test(var.equal = FALSE)) (,157)",
      { ow <- suppressMessages(effectsize::omega_squared(oneway.test(y ~ factor(gr), dw1, var.equal = FALSE), ci = .95, alternative = "two.sided"))
        e3 <- rw1$efeito[3, ]; Fw <- rw1$tabela$F[2]; g2 <- rw1$tabela$gl2[2]
        grepl("F de Welch", e3$indice) && !e3$principal && perto(e3$valor, (Fw - 1) * 2 / (Fw * 2 + g2 + 1), 1e-10) &&
          perto(c(e3$valor, e3$ic_inf, e3$ic_sup), c(ow[[1]], ow$CI_low, ow$CI_high), 1e-8) && perto(e3$valor, .157, .0005) })
teste("G-01 frase da rota Welch declara a origem do ω² e traz o aproximado pelo F de Welch",
      grepl("(ω² da ANOVA clássica, como no JASP: efeito abaixo do pequeno; pelo F de Welch, ω² aproximado = .157", rw1$apa[1], fixed = TRUE) &&
        grepl("from Welch's F, approximate ω² = .157", anova_um_fator(dw1, "y", "gr", idioma = "en")$apa[1], fixed = TRUE))
teste("G-01 rota clássica sem mudança: duas linhas em $efeito e nenhum aviso de Welch (s5)",
      nrow(r5$efeito) == 2 && !any(grepl("Welch", r5$avisos)) && !grepl("aproximado", r5$apa[1]))

# G-05: unilateral com os dados no sentido oposto ao de H1. Em s1 a média de A é ~0,5 maior
# que a de B; com alternativa = "menor" (H1: A < B), t, diferença e d saem com o sinal de
# A - B (positivos), como $tabela, e a frase diz a direção testada.
teste("G-05 t independente unilateral oposto: H1 declarada, t com o sinal de $tabela, 'no sentido oposto ao de H1'",
      { ru1 <- teste_t_independente(s1, "y", "g", alternativa = "menor"); tl <- ru1$tabela[ru1$tabela$escolhido, ]
        grepl("unilateral (H1: média de A menor que a de B)", ru1$apa, fixed = TRUE) && tl$t > 0 && tl$p > .5 &&
          grepl(.ee_estat("t", tl$t, tl$gl), ru1$apa, fixed = TRUE) && grepl("diferença média A - B = 0.", ru1$apa, fixed = TRUE) &&
          grepl("no sentido oposto ao de H1", ru1$apa, fixed = TRUE) })
teste("G-05 pareado unilateral oposto (H1: diminuição; dados aumentam): t positivo de pos - pre e direção declarada",
      { rp <- teste_t_pareado(s3, "pre", "pos", alternativa = "menor"); rq <- teste_t_pareado(s3, "pre", "pos", alternativa = "maior")
        grepl("unilateral (H1: diminuição de pre para pos)", rp$apa, fixed = TRUE) && rp$tabela$t > 0 &&
          grepl(.ee_estat("t", rp$tabela$t, 79), rp$apa, fixed = TRUE) && grepl("no sentido oposto ao de H1", rp$apa, fixed = TRUE) &&
          grepl("H1: aumento de pre para pos", rq$apa, fixed = TRUE) && !grepl("oposto", rq$apa) })
teste("G-05 Mann-Whitney, Wilcoxon pareado e uma amostra unilaterais opostos: H1 declarada e marca de sentido oposto",
      { mu1 <- teste_mann_whitney(s1, "y", "g", alternativa = "menor"); wu <- teste_wilcoxon_pareado(s3, "pre", "pos", alternativa = "menor")
        uu <- teste_t_uma_amostra(s3, "pos", mu = 50, alternativa = "menor")
        grepl("H1: valores de A tendem a ser menores que os de B", mu1$apa, fixed = TRUE) && mu1$tabela$z > 0 &&
          grepl(paste0("z = ", .pr_num(mu1$tabela$z, 2)), mu1$apa, fixed = TRUE) &&
          grepl("H1: diminuição de pre para pos", wu$apa, fixed = TRUE) && grepl("H1: média de pos menor que 50.00", uu$apa, fixed = TRUE) &&
          all(grepl("no sentido oposto ao de H1", c(mu1$apa, wu$apa, uu$apa), fixed = TRUE)) })
teste("G-05 bilateral sem mudança: nenhuma frase cita H1",
      !any(grepl("H1", c(r$apa, r3$apa, mw$apa, teste_t_uma_amostra(s3, "pre", mu = 50)$apa))))

# G-07 (reprodução g_rep_G07_levene_degenerado.R): y = 1, 2 | 4, 6. Os desvios |y - média|
# são 0,5 e 0,5 em a, 1 e 1 em b: erro do Levene igual a 0, e o F(1, 2) = 8,1e31 que decidia
# a rota era ruído numérico.
d6 <- data.frame(y = c(1, 2, 4, 6), g = c("a", "a", "b", "b"))
teste("G-07 Levene degenerado com 2 casos por grupo: sem F absurdo, rota Welch, motivo em $decisoes e $avisos, sem aviso do R",
      { avR <- character()
        r6d <- withCallingHandlers(teste_t_independente(d6, "y", "g"), warning = function(w) { avR <<- c(avR, conditionMessage(w)); invokeRestart("muffleWarning") })
        lv6 <- pres(r6d, "Levene")
        !length(avR) && is.na(lv6$estatistica) && is.na(lv6$p) && grepl("degenerado", lv6$leitura) && r6d$tabela$escolhido[2] &&
          any(grepl("Levene está degenerado", r6d$decisoes)) && !any(grepl("81129638|e\\+31", c(r6d$decisoes, capture.output(print(r6d))))) &&
          any(grepl("Levene degenerado", r6d$avisos)) })
teste("G-07 Levene com variância zero num grupo (n = 6) ou com grupo de 2 casos na ANOVA: não decide a rota",
      { dz <- data.frame(y = c(rep(5, 6), 1, 3, 5, 7, 9, 4), g = rep(c("a", "b"), each = 6))
        tz <- teste_t_independente(dz, "y", "g")
        dz3 <- data.frame(y = c(1, 2, 3.5, 4, 5, 6.2, 7, 8, 9.1, 10), g = c("a", "a", rep("b", 4), rep("c", 4)))
        az <- anova_um_fator(dz3, "y", "g")
        grepl("variância zero em a", pres(tz, "Levene")$leitura) && tz$tabela$escolhido[2] &&
          grepl("grupo com 2 casos (a)", pres(az, "Levene")$leitura, fixed = TRUE) && az$tabela$escolhido[2] && any(grepl("Levene degenerado", az$avisos)) })
teste("G-07 VD constante (geral ou dentro de todos os grupos) para com mensagem clara, não com erro cru do R",
      { msg <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
        cte <- data.frame(y = rep(5, 20), g = rep(c("a", "b"), each = 10)); cte3 <- data.frame(y = rep(5, 30), g = rep(c("a", "b", "c"), each = 10))
        sep <- data.frame(y = c(rep(5, 10), rep(7, 10)), g = rep(c("a", "b"), each = 10)); sep3 <- data.frame(y = rep(c(5, 7, 9), each = 10), g = rep(c("a", "b", "c"), each = 10))
        m <- c(msg(teste_t_independente(cte, "y", "g")), msg(anova_um_fator(cte3, "y", "g")), msg(teste_mann_whitney(cte, "y", "g")),
               msg(teste_kruskal_wallis(cte3, "y", "g", n_boot = 20)), msg(teste_t_uma_amostra(cte, "y", mu = 3)))
        m2 <- c(msg(teste_t_independente(sep, "y", "g")), msg(anova_um_fator(sep3, "y", "g")), msg(teste_mann_whitney(sep, "y", "g")))
        mp <- msg(teste_t_pareado(data.frame(a = 1:6, b = 2:7), "a", "b"))
        all(grepl("não varia (todos os casos", m, fixed = TRUE)) && all(grepl("não varia dentro de nenhum grupo", m2, fixed = TRUE)) &&
          grepl("não variam", mp) && !any(grepl("essentially constant|residual sum of squares|missing value", c(m, m2, mp))) &&
          nrow(teste_kruskal_wallis(sep3, "y", "g", n_boot = 20)$tabela) == 1 })

# G-11: separador de itens pelo decimal na frase das inclinações da ANCOVA e nos textos de
# $decisoes (antes: round() com ponto e zero à esquerda, "p = 0.133", com decimal = ",").
teste("G-11 ANCOVA simulada com inclinações heterogêneas e decimal = ',': '; p = ,' na frase das inclinações",
      { a12 <- ancova(s12, "y", "g", "x", decimal = ",")
        grepl("; p < ,001); a ANCOVA", a12$apa[2], fixed = TRUE) && !grepl("[0-9], p ", a12$apa[2]) })
teste("G-11 frases com decimal = ',': pós-hoc sem par significativo 'p = ,077; 1 - 2' (cenário do G-01) e Mauchly da RM simulada com '; p '",
      { ph1 <- anova_um_fator(dw1, "y", "gr", decimal = ",")$apa[2]; rm8 <- anova_medidas_repetidas(s8, c("t1", "t2", "t3"), decimal = ",")$apa[1]
        grepl("(menor p ajustado: p = ,077; 1 - 2)", ph1, fixed = TRUE) && grepl("após o teste de Mauchly, W = ,[0-9]{2}; χ²\\(2\\) = [0-9]+,[0-9]{2}; p [<=] ,[0-9]{3}", rm8) &&
          !grepl("[0-9], p ", rm8) })
teste("G-11 leituras de pressupostos com decimal = ',': Mauchly 'ε GG = 0,xxx; HF = ...' e inclinações da ANCOVA 'c: 0,xx; t1: ...'",
      { lm8 <- pres(anova_medidas_repetidas(s8, c("t1", "t2", "t3"), decimal = ","), "Mauchly")$leitura
        li12 <- pres(ancova(s12, "y", "g", "x", decimal = ","), "inclinações")$leitura
        grepl("ε GG = [01],[0-9]{3}; HF = [01],[0-9]{3}; LB = 0,500", lm8) &&
          grepl("\\(c: -?[0-9]+,[0-9]{2}; t1: -?[0-9]+,[0-9]{2}; t2: -?[0-9]+,[0-9]{2}\\)", li12) && !grepl("[0-9]\\.[0-9]", li12) })
teste("G-11 decisões do t, do pareado e da ANOVA com decimal = ',': sem '0.' nem ', p = '",
      { dv1 <- paste(teste_t_independente(s1, "y", "g", decimal = ",")$decisoes, collapse = " ")
        dv2 <- paste(teste_t_pareado(s3, "pre", "pos", decimal = ",")$decisoes, collapse = " ")
        dv3 <- paste(anova_um_fator(s5, "y", "g", decimal = ",")$decisoes, collapse = " ")
        grepl("a outra versão dá F = [0-9]+,[0-9]{2}; p = ,[0-9]{3}", dv1) && grepl("d_rm = [0-9]+,[0-9]{3}", dv2) &&
          !any(grepl("p = 0\\.|[0-9]\\.[0-9]{3}|[0-9], p = ", c(dv1, dv2, dv3))) })
teste("textos novos (avisos, decisões e frases dos achados G e H) sem travessão nem palavras vetadas",
      { txt <- unlist(lapply(list(rw1, teste_t_independente(d6, "y", "g"), teste_t_independente(dp2, "x", "g"), teste_t_pareado(p4, "pre", "pos"),
                                  teste_t_independente(s1, "y", "g", alternativa = "menor"), teste_mann_whitney(s1, "y", "g", alternativa = "menor")),
                             function(x) c(x$apa, x$decisoes, x$avisos, x$pressupostos$leitura)))
        !any(grepl(enc2utf8("\u2014|genuinamente|de fato|notavelmente|importante notar|intricado|meticuloso"), enc2utf8(txt), useBytes = TRUE)) })

cat("== Equivalência (TOST) e aviso de não significativo (v1.3)\n")
# Equivalência com outra implementação: TOSTER::t_TOST (pacote TOSTER, de Lakens, 2017), que
# os evals usam só para conferir. Tolerância numérica 1e-8 (mesma conta, outra implementação).
tem_toster <- requireNamespace("TOSTER", quietly = TRUE)
set.seed(1301)
q1 <- data.frame(y = c(rnorm(40, 10, 2), rnorm(35, 10.3, 2.5)), g = rep(c("A", "B"), c(40, 35)))
q2 <- data.frame(pre = rnorm(30, 20, 4)); q2$pos <- q2$pre + rnorm(30, .4, 2)
compara_tost <- function(r, tt, lim) {
  lo <- tt$TOST["TOST Lower", ]; up <- tt$TOST["TOST Upper", ]; t0 <- tt$TOST["t-test", ]
  perto(r$tabela$estatistica, c(lo$t, up$t, t0$t), 1e-8) && perto(r$tabela$p, c(lo$p.value, up$p.value, t0$p.value), 1e-8) &&
    perto(r$tabela$gl[1], lo$df, 1e-8) && perto(c(r$efeito$ic_inf[1], r$efeito$ic_sup[1]), unlist(tt$effsize["Raw", c("lower.ci", "upper.ci")]), 1e-8) &&
    perto(r$limites, lim, 1e-8)
}
if (!tem_toster) cat("PULADO  equivalência com TOSTER (pacote ausente; está em .ee_pacotes_evals)\n") else {
  teste("E-01 TOST independente (Welch) = TOSTER::t_TOST: t, gl, p dos dois testes unilaterais e do usual, IC 90%",
        local({ r <- teste_equivalencia(q1, "independente", limites = 1, vd = "y", grupo = "g")
          tt <- TOSTER::t_TOST(x = q1$y[q1$g == "A"], y = q1$y[q1$g == "B"], eqb = 1, var.equal = FALSE)
          compara_tost(r, tt, c(-1, 1)) && r$efeito$nivel[1] == .90 }))
  teste("E-02 TOST independente (Student, welch = FALSE) e limites assimétricos = TOSTER",
        local({ r <- teste_equivalencia(q1, "independente", limites = c(-.8, 1.2), vd = "y", grupo = "g", welch = FALSE)
          tt <- TOSTER::t_TOST(x = q1$y[q1$g == "A"], y = q1$y[q1$g == "B"], eqb = c(-.8, 1.2), var.equal = TRUE)
          compara_tost(r, tt, c(-.8, 1.2)) }))
  teste("E-03 TOST pareado (depois - antes) = TOSTER com paired = TRUE",
        local({ r <- teste_equivalencia(q2, "pareado", limites = 1.5, antes = "pre", depois = "pos")
          tt <- TOSTER::t_TOST(x = q2$pos, y = q2$pre, paired = TRUE, eqb = 1.5)
          compara_tost(r, tt, c(-1.5, 1.5)) }))
  # Uma amostra: a skill segue Lakens (2017, eq. 7), com limites relativos a mu; o TOSTER 0.8.6
  # recebe os limites na escala da média (mu ± limite) e dá o IC da média, não o da diferença.
  teste("E-04 TOST de uma amostra contra mu = 20 com limites ±2 = TOSTER com eqb = c(18, 22) (IC deslocado por mu)",
        local({ r <- teste_equivalencia(q2, "uma_amostra", limites = 2, vd = "pre", mu = 20)
          tt <- TOSTER::t_TOST(x = q2$pre, mu = 20, eqb = c(18, 22))
          lo <- tt$TOST["TOST Lower", ]; up <- tt$TOST["TOST Upper", ]; t0 <- tt$TOST["t-test", ]
          perto(r$tabela$estatistica, c(lo$t, up$t, t0$t), 1e-8) && perto(r$tabela$p, c(lo$p.value, up$p.value, t0$p.value), 1e-8) &&
            perto(c(r$efeito$ic_inf[1], r$efeito$ic_sup[1]) + 20, unlist(tt$effsize["Raw", c("lower.ci", "upper.ci")]), 1e-8) }))
  teste("E-05 limites em d (Welch) = TOSTER eqbound_type = 'SMD' (DP de referência = raiz da média das variâncias) e o aviso de viés sai",
        local({ r <- teste_equivalencia(q1, "independente", limites = .4, vd = "y", grupo = "g", escala = "d")
          tt <- suppressWarnings(suppressMessages(TOSTER::t_TOST(x = q1$y[q1$g == "A"], y = q1$y[q1$g == "B"], eqb = .4,
                                                                 var.equal = FALSE, eqbound_type = "SMD")))
          dref <- sqrt((var(q1$y[q1$g == "A"]) + var(q1$y[q1$g == "B"])) / 2)
          compara_tost(r, tt, c(-.4, .4) * dref) && any(grepl("enviesado", r$avisos)) }))
  teste("E-06 pareado com limites em d_z = TOSTER 'SMD' pareado (DP das diferenças)",
        local({ r <- teste_equivalencia(q2, "pareado", limites = .5, antes = "pre", depois = "pos", escala = "d")
          tt <- suppressWarnings(suppressMessages(TOSTER::t_TOST(x = q2$pos, y = q2$pre, paired = TRUE, eqb = .5, eqbound_type = "SMD")))
          compara_tost(r, tt, c(-.5, .5) * sd(q2$pos - q2$pre)) }))
}
cat("== Equivalência: verdade conhecida\n")
set.seed(1302)
# Grupos com os mesmos valores: diferença observada exatamente 0 (p = 1), sem depender da semente.
z0 <- rnorm(500); v0 <- data.frame(y = c(z0, z0), g = rep(c("A", "B"), each = 500))
teste("E-07 diferença real 0, n = 500 por grupo, limites ±0,3: equivalente e não diferente de zero",
      local({ r <- teste_equivalencia(v0, "independente", limites = .3, vd = "y", grupo = "g")
        r$equivalente && !r$diferente_de_zero && r$resultado == "equivalente e não diferente de zero" && grepl("indicou equivalência", r$apa) }))
v1 <- data.frame(y = c(rnorm(500, .5), rnorm(500)), g = rep(c("A", "B"), each = 500))
teste("E-08 diferença real 0,5 com limites ±0,2: diferente de zero e não equivalente",
      local({ r <- teste_equivalencia(v1, "independente", limites = .2, vd = "y", grupo = "g")
        !r$equivalente && r$diferente_de_zero && r$resultado == "diferente de zero e não equivalente" && grepl("não indicou equivalência", r$apa) }))
v2 <- data.frame(y = c(rnorm(5000, .1), rnorm(5000)), g = rep(c("A", "B"), each = 5000))
teste("E-09 diferença real 0,1 com n = 5000 por grupo e limites ±0,3: diferente de zero, mas dentro dos limites",
      local({ r <- teste_equivalencia(v2, "independente", limites = .3, vd = "y", grupo = "g")
        r$equivalente && r$diferente_de_zero && grepl("dentro dos limites", r$resultado) }))
z3 <- rnorm(10); v3 <- data.frame(y = c(z3, z3), g = rep(c("A", "B"), each = 10))
teste("E-10 n = 10 por grupo e limites ±0,2: IC 90% mais largo que os limites gera aviso e resultado indeterminado sem frase de igualdade",
      local({ r <- teste_equivalencia(v3, "independente", limites = .2, vd = "y", grupo = "g")
        any(grepl("mais largo que o intervalo de equivalência", r$avisos)) && grepl("^indeterminado", r$resultado) &&
          any(grepl("não escrever que os grupos são iguais", r$avisos)) }))
teste("E-11 sem limites a função para e pede o menor efeito de interesse",
      inherits(tryCatch(teste_equivalencia(v0, "independente", vd = "y", grupo = "g"), error = function(e) e), "error"))
teste("E-12 frase em inglês e decimal = ','",
      local({ a <- teste_equivalencia(v0, "independente", limites = .3, vd = "y", grupo = "g", idioma = "en")$apa
        b <- teste_equivalencia(v0, "independente", limites = .3, vd = "y", grupo = "g", decimal = ",")$apa
        grepl("^An equivalence test \\(TOST, Welch\\)", a) && grepl("90% CI", a, fixed = TRUE) && grepl("IC 90% [-", b, fixed = TRUE) && grepl("p < ,001", b, fixed = TRUE) }))
cat("== Aviso de não significativo\n")
ns_txt <- "Não significativo não quer dizer igual"
teste("E-13 t independente não significativo: aviso cita Greenland et al. (2016) e aponta teste_equivalencia(tipo = \"independente\")",
      local({ a <- teste_t_independente(v0, "y", "g")$avisos
        any(grepl(ns_txt, a, fixed = TRUE)) && any(grepl("Greenland et al., 2016", a, fixed = TRUE)) &&
          any(grepl("teste_equivalencia(tipo = \"independente\")", a, fixed = TRUE)) }))
teste("E-14 t independente significativo: sem o aviso",
      !any(grepl(ns_txt, teste_t_independente(v1, "y", "g")$avisos, fixed = TRUE)))
teste("E-15 pareado e uma amostra apontam o tipo certo; Mann-Whitney diz que não há TOST para a estatística",
      local({ e3 <- rnorm(40); q3 <- data.frame(pre = rnorm(40, 20, 4)); q3$pos <- q3$pre + e3 - mean(e3)  # mudança média exatamente 0
        u3 <- rnorm(50); u3 <- data.frame(y = u3 - mean(u3))                                          # média exatamente 0
        ap <- teste_t_pareado(q3, "pre", "pos")$avisos
        au <- teste_t_uma_amostra(u3, "y", mu = 0)$avisos
        am <- teste_mann_whitney(v0, "y", "g")$avisos
        any(grepl("tipo = \"pareado\"", ap, fixed = TRUE)) &&
          any(grepl("tipo = \"uma_amostra\"", au, fixed = TRUE)) &&
          any(grepl("A skill não tem teste de equivalência para esta estatística", am, fixed = TRUE)) }))
teste("E-16 ANOVA fatorial: o aviso lista só os termos não significativos",
      local({ set.seed(1303)
        f2 <- expand.grid(a = c("a1", "a2"), b = c("b1", "b2"), rep = 1:40)
        e2 <- rnorm(nrow(f2)); e2 <- e2 - ave(e2, f2$a, f2$b)  # resíduo com média 0 em cada célula
        f2$y <- e2 + 1.2 * (f2$a == "a2")                     # só o efeito de a; b e a × b exatamente 0
        av <- anova_fatorial(f2, "y", c("a", "b"))$avisos
        x <- av[grepl(ns_txt, av, fixed = TRUE)]
        length(x) == 1 && grepl("em b (", x, fixed = TRUE) && !grepl("em a (", x, fixed = TRUE) }))
teste("E-17 textos novos (TOST e aviso) sem travessão nem palavras vetadas",
      local({ rs <- list(teste_equivalencia(v0, "independente", limites = .3, vd = "y", grupo = "g"),
                   teste_equivalencia(v3, "independente", limites = .2, vd = "y", grupo = "g"),
                   teste_equivalencia(v0, "independente", limites = .3, vd = "y", grupo = "g", idioma = "en"),
                   teste_equivalencia(q1, "independente", limites = .4, vd = "y", grupo = "g", escala = "d"),
                   teste_t_independente(v0, "y", "g"), teste_mann_whitney(v0, "y", "g"))
        txt <- unlist(lapply(rs, function(x) c(x$apa, x$decisoes, x$avisos)))
        !any(grepl(enc2utf8("\u2014|\u2013|genuinamente|de fato|notavelmente|importante notar|intricado|meticuloso|cabe ressaltar|vale destacar|essencialmente|fundamentalmente|abrangente"),
                   enc2utf8(txt), useBytes = TRUE)) }))

cat("== Formato e privacidade\n")
teste("decimal = ',' gera p = ,036 e gl separados por ponto e vírgula",
      { a <- anova_um_fator(s5, "y", "g", decimal = ",")$apa[1]; grepl("p < ,001", a) && grepl("F(2; 117)", a, fixed = TRUE) })
teste("idioma = 'en' gera frase em inglês", grepl("^A one-way ANOVA", anova_um_fator(s5, "y", "g", idioma = "en")$apa[1]))
teste("nenhuma tabela impressa tem uma linha por participante",
      { rs <- list(r, r3, r5, r7, r9, r10, r11, r13, r14, r15)
        all(vapply(rs, function(x) all(vapply(x[c("tabela", "efeito", "pressupostos", "posthoc")], function(d) is.null(d) || NROW(d) < 30, TRUE)), TRUE)) })
teste("frases APA sem números crus vazados (nenhum número com 6+ casas decimais)",
      { txt <- unlist(lapply(list(r, r2, r3, mw, r5, r6, r7, r8, r9, r10, r11, r13, r14, r15), function(x) x$apa))
        !any(grepl("[0-9][.,][0-9]{6,}", txt)) })
teste("sem travessão longo nem palavras vetadas nos textos gerados",
      { txt <- unlist(lapply(list(r, r3, r5, r7, r9, r10, r11, r13, r14, r15), function(x) c(x$apa, x$decisoes, x$avisos, x$pressupostos$leitura)))
        !any(grepl(enc2utf8("\u2014|genuinamente|de fato|notavelmente|importante notar|intricado|meticuloso"), enc2utf8(txt), useBytes = TRUE)) })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
