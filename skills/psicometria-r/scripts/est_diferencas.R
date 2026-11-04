# est_diferencas.R: comparação de grupos e condições (módulo "diferenças").
# Carregar por scripts/carregar.R (usa est_comum.R e psicometria.R).
#
# Funções públicas (todas devolvem objeto "pr_resultado"; detalhes em
# references/diferencas.md):
#   teste_t_independente(dados, vd, grupo, rota, centro_levene, indice, alternativa, controle, ...)
#   teste_t_pareado(dados, antes, depois, alternativa, ...)
#   teste_t_uma_amostra(dados, vd, mu, alternativa, ...)
#   teste_mann_whitney(dados, vd, grupo, sinal_rrb, regua_rrb, alternativa, exato, ...)
#   teste_wilcoxon_pareado(dados, antes, depois, sinal_rrb, regua_rrb, alternativa, exato, ...)
#   anova_um_fator(dados, vd, grupo, rota, posthoc, exigir_omnibus, centro_levene, ...)
#   anova_medidas_repetidas(dados, medidas, intra, entre, esfericidade, limiar_epsilon, posthoc, ...)
#   anova_fatorial(dados, vd, fatores, tipo_sq, posthoc, centro_levene, ...)
#   ancova(dados, vd, fator, covariaveis, posthoc, bootstrap, n_boot, semente, ...)
#   manova_grupos(dados, vds, fatores, tipo_sq, alfa_box, ...)
#   teste_kruskal_wallis(dados, vd, grupo, correcao, exigir_omnibus, sinal_rrb, n_boot, semente, ...)
#   teste_friedman(dados, medidas, nomes, correcao, exigir_omnibus, ...)
#   teste_equivalencia(dados, tipo, limites, vd, grupo, antes, depois, x, y, mu, escala,
#                      welch, alfa, decimal, idioma)   TOST; IC de 1 - 2α, sem nivel
# "..." = alfa = .05, nivel = .95, decimal = ".", idioma = "pt" (em todas).
# Teste principal não significativo: nota de .ee_aviso_nao_sig() (est_comum.R) em $avisos.
#
# Regra geral: o padrão segue o livro do usuário (Vol. 4 e 5, JASP 0.95.1); a
# alternativa recomendada pela literatura é calculada e mostrada ao lado quando é
# barata, e a escolha fica escrita em $decisoes. Convenção de sinal das comparações
# par a par: rótulo "A - B" = média (ou postos) de A menos B.

# ================================================================ auxiliares
.ee_dif_alt <- function(alternativa) {
  c(bilateral = "two.sided", maior = "greater", menor = "less")[[alternativa]]
}

.ee_dif_sw <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(x) > 5000 || stats::sd(x) == 0) return(c(W = NA_real_, p = NA_real_))
  s <- stats::shapiro.test(x)
  c(W = unname(s$statistic), p = s$p.value)
}

# Leitura de normalidade por magnitude (assimetria e curtose, Vol. 1 do livro), com o
# Shapiro-Wilk como um indício entre vários; reaproveita ler_normalidade() da base.
.ee_dif_le_normal <- function(x, rotulo) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || stats::sd(x) == 0) return(paste0(rotulo, ": normalidade indeterminada (n < 3 ou sem variabilidade)"))
  sw <- .ee_dif_sw(x)
  d <- data.frame(variavel = rotulo, n = length(x), assimetria = .ee_assimetria(x),
                  curtose = .ee_curtose(x), shapiro_p = sw[["p"]])
  unname(ler_normalidade(d))
}

.ee_dif_pres <- function(checagem, estatistica, gl, p, leitura) {
  data.frame(checagem = checagem, estatistica = estatistica, gl = as.character(gl), p = p,
             leitura = leitura, stringsAsFactors = FALSE)
}

.ee_dif_gltxt <- function(...) paste(.ee_gl(c(...)), collapse = ", ")

# Levene: |y - centro do grupo| numa ANOVA de um fator. centro "media" = Levene
# clássico (o do livro e o do painel de ANOVA do JASP); "mediana" = variante de
# Brown e Forsythe (padrão do painel de teste t do JASP 0.95.1 e do car::leveneTest).
.ee_dif_levene <- function(y, g, centro = c("media", "mediana")) {
  centro <- match.arg(centro)
  fun <- if (centro == "media") mean else stats::median
  cen <- stats::ave(y, g, FUN = fun)
  # ajuste perfeito (Levene degenerado) faz o anova() avisar; quem decide é .ee_dif_levene_degenerado()
  a <- suppressWarnings(stats::anova(stats::lm(abs(y - cen) ~ g)))
  c(F = a[1, "F value"], gl1 = a[1, "Df"], gl2 = a[2, "Df"], p = a[1, "Pr(>F)"])
}

.ee_dif_le_levene <- function(p, alfa, idioma = "pt") {
  if (is.na(p)) return("indeterminado")
  if (p < alfa) "evidência de variâncias desiguais" else
    "sem evidência de variâncias desiguais (o que não prova que sejam iguais)"
}

# VD sem variação. O critério de "constante" é o do t.test do R (desvio-padrão abaixo de
# 10 vezes o épsilon da máquina na escala da média). Com a VD constante, ou constante
# dentro de todos os grupos, o erro é zero e o R para com mensagem crua ("data are
# essentially constant", "residual sum of squares is 0"); aqui a parada diz o motivo.
.ee_dif_constante <- function(x) {
  x <- x[is.finite(x)]
  length(x) < 2 || isTRUE(stats::sd(x) <= 10 * .Machine$double.eps * max(1, abs(mean(x))))
}
.ee_dif_exige_variacao <- function(y, vd, g = NULL, grupo = NULL, sugestao = NULL) {
  if (.ee_dif_constante(y)) stop("a VD ", vd, " não varia (todos os casos têm o mesmo valor): não há diferença a testar", call. = FALSE)
  if (!is.null(g) && all(vapply(split(y, g), .ee_dif_constante, logical(1))))
    stop("a VD ", vd, " não varia dentro de nenhum grupo de ", grupo, ": a variância do erro é zero e o teste não é calculável",
         " (os grupos se separam sem sobreposição; descrever os valores por grupo", if (!is.null(sugestao)) paste0(" ou usar ", sugestao), ")", call. = FALSE)
}

# Levene degenerado. Com 2 casos num grupo, os dois desvios absolutos ao centro são
# iguais e o grupo não contribui para o erro do Levene; com variância zero num grupo,
# todos os desvios dele são 0. Quando nenhum grupo sobra para o erro, o ajuste é perfeito
# e o F vira ruído numérico (F(1, 2) = 8,1e31 com 2 casos por grupo, achado G-07 da
# revisão CoVe final). Devolve list(motivo = "" se não degenerado, ruido = TRUE quando o F
# não tem valor numérico interpretável).
.ee_dif_levene_degenerado <- function(y, g, lev, centro = "media") {
  ns <- table(g)
  fun <- if (centro == "media") mean else stats::median
  z <- abs(y - stats::ave(y, g, FUN = fun))
  ssw <- sum((z - stats::ave(z, g))^2); sst <- sum((z - mean(z))^2)
  ruido <- !is.finite(lev[["F"]]) || !is.finite(lev[["p"]]) || ssw <= 1e-10 * max(sst, .Machine$double.eps)
  zero <- names(ns)[vapply(levels(g), function(l) .ee_dif_constante(y[g == l]), logical(1))]
  mot <- c(if (ruido) "ajuste perfeito: o F é ruído numérico",
           if (any(ns <= 2)) paste0("grupo com 2 casos (", paste(names(ns)[ns <= 2], collapse = ", "), "), cujos desvios ao centro são iguais por construção"),
           if (length(zero)) paste0("variância zero em ", paste(zero, collapse = ", ")))
  list(motivo = paste(mot, collapse = "; "), ruido = ruido)
}

# Privacidade na frase: grupo (ou momento) com n < 5 não tem M, DP nem Md na frase, a
# mesma regra que o print aplica à tabela de descritivas (est_comum.R, .ee_n_min_print).
# Com n = 2, M e DP devolvem os dois escores (achado H-02 da revisão CoVe final).
.ee_dif_peq <- function(n) !is.na(n) & n < .ee_n_min_print
.ee_dif_omit_txt <- function(idioma) if (idioma == "en") "n < 5, descriptives omitted" else "n < 5, descritivas omitidas"
.ee_dif_aviso_peq <- function(rotulos, o_que) paste0("Grupo ou momento com n < 5 (", paste(rotulos, collapse = ", "), "): sem ", o_que,
                                                     " na frase APA, como na tabela de descritivas do print (privacidade); os valores ficam no objeto.")

# Teste unilateral: texto da hipótese alternativa e marca de dados no sentido oposto.
.ee_dif_oposto <- function(alt, estat) (alt == "greater" && isTRUE(estat < 0)) || (alt == "less" && isTRUE(estat > 0))
.ee_dif_txt_oposto <- function(oposto, idioma) if (!oposto) "" else if (idioma == "en") ", opposite to H1" else ", no sentido oposto ao de H1"

# IC de um d pela inversão do t não central (a conf.limits.nct do MBESS, que o JASP
# copia): acha os ncp extremos para o t observado e reescala por 'fator'.
.ee_dif_ic_nct <- function(t, gl, fator, nivel = .95, alt = "two.sided") {
  a <- 1 - nivel
  achar <- function(prob) {
    f <- function(ncp) suppressWarnings(stats::pt(t, gl, ncp)) - prob
    lo <- t - 10 - abs(t); hi <- t + 10 + abs(t)
    while (f(lo) < 0) lo <- lo - 10
    while (f(hi) > 0) hi <- hi + 10
    stats::uniroot(f, c(lo, hi), tol = 1e-10)$root
  }
  lo <- if (alt == "less") -Inf else achar(if (alt == "two.sided") 1 - a / 2 else 1 - a)
  hi <- if (alt == "greater") Inf else achar(if (alt == "two.sided") a / 2 else a)
  c(lo, hi) * fator
}

# Correção de viés de Hedges, forma gama exata (a do JASP e do effectsize).
.ee_dif_J <- function(gl) exp(lgamma(gl / 2) - log(sqrt(gl / 2)) - lgamma((gl - 1) / 2))

.ee_dif_hedges_par <- function(a, b) {
  n1 <- length(a); n2 <- length(b); gl <- n1 + n2 - 2
  sp <- sqrt(((n1 - 1) * stats::var(a) + (n2 - 1) * stats::var(b)) / gl)
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(a) - mean(b)) / sp * .ee_dif_J(gl)
}

# Rótulos verbais. d: pontos de Cohen (1988) .20/.50/.80 usados como limiares
# inferiores (leitura aproximada, não fronteira rígida; Cohen definiu pontos, não
# faixas). Proporções de variância: .01/.06/.14 (derivados de f = .10/.25/.40 de
# Cohen; aplicar a ω² é extrapolação usual). Correlações de postos e W de Kendall:
# convenção de r, .10/.30/.50; regua "grissom_kim" usa os cortes do curso.
.ee_dif_rotulo <- function(valor, tipo = c("d", "var", "r", "grissom_kim"), idioma = "pt") {
  tipo <- match.arg(tipo)
  a <- abs(valor)
  if (is.na(a)) return("")
  cortes <- switch(tipo, d = c(.2, .5, .8), var = c(.01, .06, .14), r = c(.1, .3, .5),
                   grissom_kim = c(.125, .304, .465))
  i <- findInterval(a, cortes) + 1
  pt <- if (tipo == "var") c("abaixo do pequeno", "pequeno", "médio", "grande") else
    c("abaixo do pequeno", "pequeno", "moderado", "grande")
  en <- c("below small", "small", "medium", "large")
  s <- if (idioma == "en") en[i] else pt[i]
  if (tipo == "d" && i %in% 2:3 && a >= cortes[i] - .05) {
    s <- paste0(s, if (idioma == "en") ", near the next benchmark" else
      if (i == 2) ", próximo do moderado" else ", próximo do grande")
  }
  s
}

.ee_dif_efeito_df <- function(indice, valor, lo = NA_real_, hi = NA_real_, principal = FALSE, termo = NULL) {
  d <- data.frame(indice = indice, valor = valor, ic_inf = lo, ic_sup = hi, principal = principal,
                  stringsAsFactors = FALSE)
  if (!is.null(termo)) d <- cbind(termo = termo, d, stringsAsFactors = FALSE)
  d
}

# arredondamento meio-para-longe-do-zero (o do JASP/SPSS), imune à representação binária
.ee_dif_m <- function(x, decimal, casas = 2) .pr_num(x + sign(x) * 1e-9 * max(1, abs(x)), casas, TRUE, decimal)

# estatísticas de postos (U, V): inteiro sem casas, meio com uma casa
.ee_dif_u <- function(x, decimal) .pr_num(x, if (abs(x - round(x)) < 1e-9) 0 else 1, TRUE, decimal)

.ee_dif_desc <- function(y, g) {
  do.call(rbind, lapply(levels(g), function(l) {
    x <- y[g == l]
    data.frame(grupo = l, n = length(x), media = mean(x), dp = stats::sd(x),
               ep = stats::sd(x) / sqrt(length(x)), mediana = stats::median(x), stringsAsFactors = FALSE)
  }))
}

# Regra do p exato do Mann-Whitney no R até 4.5 (o R 4.5.1 do JASP 0.95.1): exato só com
# menos de 50 casos por grupo e sem empates; senão, aproximação normal com correção de
# continuidade e de empates. O R 4.6 passou a dar p exato condicional com empates.
# exato = NULL aplica a regra antiga; FALSE força a normal; TRUE chama o wilcox.test do R
# instalado com exact = TRUE (no R 4.6, o exato condicional).
.ee_dif_exato <- function(exato, n_max, empates) {
  if (is.null(exato)) n_max < 50 && !empates else isTRUE(exato)
}

# Wilcoxon de postos sinalizados exatamente como o R 4.5.1 (o do JASP 0.95.1): zeros saem
# antes de ranquear; p exato (psignrank) só com n < 50 diferenças não nulas, sem empates
# e SEM zeros; havendo zeros ou empates, normal com correção de continuidade e de
# empates. exato = FALSE força a normal; exato = TRUE chama o wilcox.test do R instalado
# com exact = TRUE sobre todas as diferenças (no R 4.6: exato condicional, com os zeros
# nos postos, o que muda V; no R 4.5 com zeros ou empates o próprio R cai na normal).
.ee_dif_wsr <- function(d, alt, exato, nivel) {
  zeros <- sum(d == 0); x <- d[d != 0]; n <- length(x)
  r <- rank(abs(x)); V <- sum(r[x > 0]); empates <- length(r) != length(unique(r))
  if (isTRUE(exato)) {
    w <- suppressWarnings(stats::wilcox.test(d, mu = 0, alternative = alt, exact = TRUE, conf.int = TRUE, conf.level = nivel))
    return(list(V = unname(w$statistic), p = w$p.value, exato = grepl("exact", w$method), hl = unname(w$estimate),
                ic = w$conf.int, n = n, zeros = zeros, metodo = "wilcox.test do R instalado, exact = TRUE"))
  }
  usa_ex <- is.null(exato) && n < 50 && !empates && zeros == 0
  if (usa_ex) {
    p <- switch(alt,
                two.sided = min(1, 2 * (if (V > n * (n + 1) / 4) stats::psignrank(V - 1, n, lower.tail = FALSE) else stats::psignrank(V, n))),
                greater = stats::psignrank(V - 1, n, lower.tail = FALSE), less = stats::psignrank(V, n))
  } else {
    tb <- table(r)
    z <- V - n * (n + 1) / 4
    sg <- sqrt(n * (n + 1) * (2 * n + 1) / 24 - sum(tb^3 - tb) / 48)
    z <- (z - switch(alt, two.sided = sign(z) * .5, greater = .5, less = -.5)) / sg
    p <- switch(alt, two.sided = 2 * min(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE)),
                greater = stats::pnorm(z, lower.tail = FALSE), less = stats::pnorm(z))
  }
  w <- suppressWarnings(stats::wilcox.test(x, mu = 0, alternative = alt, exact = usa_ex, correct = TRUE, conf.int = TRUE, conf.level = nivel))
  list(V = V, p = p, exato = usa_ex, hl = unname(w$estimate), ic = w$conf.int, n = n, zeros = zeros,
       metodo = if (usa_ex) "exato (sem zeros nem empates, n < 50)" else "aproximação normal com correção de continuidade e de empates")
}

# Correlação rank-biserial com IC por z de Fisher (o mesmo EP do JASP e do effectsize).
# Teste unilateral: IC unilateral no mesmo sentido (como o effectsize e como o IC do d
# nos testes t), com o limite aberto em -1 ou 1.
.ee_dif_ic_rrb <- function(r, ep, nivel, alt = "two.sided") {
  if (!is.finite(r) || abs(r) >= 1) return(c(NA_real_, NA_real_))
  if (alt == "two.sided") return(tanh(atanh(r) + c(-1, 1) * stats::qnorm(1 - (1 - nivel) / 2) * ep))
  z <- stats::qnorm(nivel)
  if (alt == "greater") c(tanh(atanh(r) - z * ep), 1) else c(-1, tanh(atanh(r) + z * ep))
}

# Frase de teste unilateral: limite infinito como "∞" e IC rotulado como unilateral.
# (useBytes: funciona em qualquer locale; o símbolo vai como os bytes UTF-8 de "∞")
.ee_dif_unilateral <- function(txt, alt, idioma = "pt") {
  if (alt == "two.sided") return(txt)
  inf <- rawToChar(as.raw(c(0xe2, 0x88, 0x9e)))
  txt <- gsub("(?<![A-Za-z])Inf(?![A-Za-z])", inf, txt, perl = TRUE, useBytes = TRUE)
  if (idioma == "en") gsub("([0-9]+)% CI \\[", "\\1% one-sided CI [", txt, useBytes = TRUE)
  else gsub("IC ([0-9]+)% \\[", "IC \\1% unilateral [", txt, useBytes = TRUE)
}

# ================================================================ testes t
# Teste t de amostras independentes. Student e Welch sempre calculados.
#  rota = "classica": o Levene decide a linha (livro, Vol. 4);
#  rota = "welch": Welch sempre (Delacre, Lakens & Leys, 2017; Caixa 13.3).
#  centro_levene = "media": Levene clássico, o valor que o livro imprime (F = 2,362 no
#  Directed Reading); "mediana": variante de Brown-Forsythe, padrão do painel de teste t
#  do JASP 0.95.1 (rotulado "Brown-Forsythe" lá).
#  indice: "d" (padrão do livro, Caixa 13.1) ou "g" de Hedges (correção de viés de
#  amostra pequena, que depende só dos gl; o livro a associa também a grupos
#  desbalanceados, o que não procede e não é reproduzido).
#  controle: nível usado como denominador do Δ de Glass (padrão: 2º nível, como o JASP).
teste_t_independente <- function(dados, vd, grupo, rota = c("classica", "welch"),
                                 centro_levene = c("media", "mediana"), indice = c("d", "g"),
                                 alternativa = c("bilateral", "maior", "menor"), controle = NULL,
                                 alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  sepa <- .ee_sep(decimal)
  rota <- match.arg(rota); centro_levene <- match.arg(centro_levene); indice <- match.arg(indice)
  alt <- .ee_dif_alt(match.arg(alternativa))
  dd <- .ee_dados(dados, c(vd, grupo), numericas = vd)
  y <- dd$d[[vd]]; g <- .ee_fator(dd$d[[grupo]], grupo)
  if (nlevels(g) != 2) stop("o teste t compara exatamente 2 grupos; ", grupo, " tem ", nlevels(g),
                            " níveis com dados. Para 3 ou mais grupos, use anova_um_fator().")
  nv <- levels(g); y1 <- y[g == nv[1]]; y2 <- y[g == nv[2]]
  n1 <- length(y1); n2 <- length(y2)
  if (min(n1, n2) < 2) stop("cada grupo precisa de pelo menos 2 casos")
  .ee_dif_exige_variacao(y, vd, g, grupo)
  m1 <- mean(y1); m2 <- mean(y2); s1 <- stats::sd(y1); s2 <- stats::sd(y2)
  st <- stats::t.test(y1, y2, var.equal = TRUE, alternative = alt, conf.level = nivel)
  we <- stats::t.test(y1, y2, var.equal = FALSE, alternative = alt, conf.level = nivel)
  lev <- .ee_dif_levene(y, g, centro_levene)
  lev_outro <- .ee_dif_levene(y, g, setdiff(c("media", "mediana"), centro_levene))
  # Levene degenerado (G-07): não decide a rota; vale a linha de Welch, que não supõe
  # variâncias iguais (a recomendação sem Levene, Delacre et al., 2017).
  lev_deg <- .ee_dif_levene_degenerado(y, g, lev, centro_levene)
  degenerado <- nzchar(lev_deg$motivo)
  usa_welch <- rota == "welch" || degenerado || isTRUE(lev[["p"]] < alfa)

  # pressupostos
  sw1 <- .ee_dif_sw(y1); sw2 <- .ee_dif_sw(y2); swr <- .ee_dif_sw(y - stats::ave(y, g))
  pres <- rbind(
    .ee_dif_pres(paste0("Normalidade em ", nv[1], " (Shapiro-Wilk)"), sw1[["W"]], "", sw1[["p"]], .ee_dif_le_normal(y1, nv[1])),
    .ee_dif_pres(paste0("Normalidade em ", nv[2], " (Shapiro-Wilk)"), sw2[["W"]], "", sw2[["p"]], .ee_dif_le_normal(y2, nv[2])),
    .ee_dif_pres("Normalidade dos resíduos (Shapiro-Wilk; o que o painel do JASP 0.95.1 mostra)", swr[["W"]], "", swr[["p"]],
                 "indício complementar; o livro lê a normalidade por grupo e por magnitude"),
    .ee_dif_pres(paste0("Homogeneidade de variâncias (Levene, centro na ", if (centro_levene == "media") "média" else "mediana", ")"),
                 if (lev_deg$ruido) NA_real_ else lev[["F"]], .ee_dif_gltxt(lev[["gl1"]], lev[["gl2"]]), if (lev_deg$ruido) NA_real_ else lev[["p"]],
                 if (degenerado) paste0("degenerado (", lev_deg$motivo, "): não decide a rota") else .ee_dif_le_levene(lev[["p"]], alfa)))

  # tamanhos de efeito
  gl_s <- n1 + n2 - 2
  fator <- sqrt((n1 + n2) / (n1 * n2))
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / gl_s)
  d_s <- (m1 - m2) / sp
  ic_ds <- .ee_dif_ic_nct(d_s / fator, gl_s, fator, nivel, alt)
  d_w <- (m1 - m2) / sqrt((s1^2 + s2^2) / 2)
  ic_dw <- .ee_dif_ic_nct(d_w / fator, unname(we$parameter), fator, nivel, alt)
  J <- .ee_dif_J(gl_s)
  if (is.null(controle)) controle <- nv[2]
  if (!controle %in% nv) stop("controle precisa ser um dos níveis: ", paste(nv, collapse = ", "))
  s_c <- if (controle == nv[1]) s1 else s2
  n_c <- if (controle == nv[1]) n1 else n2
  glass <- (m1 - m2) / s_c
  ic_gl <- .ee_dif_ic_nct(glass / fator, n_c - 1, fator, nivel, alt)
  princ <- if (indice == "g") 3 else if (usa_welch) 2 else 1
  ef <- .ee_dif_efeito_df(
    c("d de Cohen (DP combinado; Student)", "d de Cohen (média das variâncias; Welch)", "g de Hedges",
      paste0("Δ de Glass (DP de ", controle, ")")),
    c(d_s, d_w, d_s * J, glass), c(ic_ds[1], ic_dw[1], ic_ds[1] * J, ic_gl[1]),
    c(ic_ds[2], ic_dw[2], ic_ds[2] * J, ic_gl[2]), seq_len(4) == princ)

  tab <- data.frame(teste = c("t de Student", "t de Welch"),
                    t = c(unname(st$statistic), unname(we$statistic)),
                    gl = c(unname(st$parameter), unname(we$parameter)),
                    p = c(st$p.value, we$p.value), dif_media = m1 - m2,
                    ic_inf = c(st$conf.int[1], we$conf.int[1]), ic_sup = c(st$conf.int[2], we$conf.int[2]),
                    escolhido = c(!usa_welch, usa_welch), stringsAsFactors = FALSE)
  desc <- .ee_dif_desc(y, g)

  # decisões
  ltxt <- if (lev_deg$ruido) paste0("Levene sem valor interpretável (", lev_deg$motivo, ")") else
    paste0("Levene ", .ee_estat("F", lev[["F"]], c(lev[["gl1"]], lev[["gl2"]]), decimal), sepa, .pr_p(lev[["p"]], decimal))
  dec <- c(
    paste0("Sinal: diferença = ", nv[1], " - ", nv[2], " (ordem dos níveis do fator, a mesma do JASP); o sinal indica direção, não mérito.",
           if (alt != "two.sided") paste0(" Teste unilateral (H1: média de ", nv[1], if (alt == "greater") " maior" else " menor", " que a de ", nv[2],
                                          "): a frase dá t, diferença, d e IC com sinal, nessa ordem, e diz quando os dados vão no sentido oposto ao de H1.") else ""),
    if (rota == "classica" && degenerado) paste0("Rota clássica do livro, mas o Levene está degenerado (", lev_deg$motivo,
                                                 "): ele não decide a rota; usada a linha de Welch, que não supõe variâncias iguais (Delacre et al., 2017; Caixa 13.3).")
    else if (rota == "classica") paste0("Rota clássica do livro: ", ltxt, " -> ",
                                   if (usa_welch) "variâncias desiguais, linha de Welch." else "sem evidência de variâncias desiguais, linha de Student.",
                                   " A recomendação de usar Welch sempre (Delacre et al., 2017; Caixa 13.3) está ao lado na tabela.")
    else paste0("Rota Welch por padrão (Delacre et al., 2017), sem consultar o Levene; ", ltxt, " fica só como descrição."),
    paste0("Levene centrado na ", if (centro_levene == "media") "média (o que o livro imprime: F = 2,362 no Directed Reading)" else "mediana (padrão do painel de teste t do JASP 0.95.1, rotulado 'Brown-Forsythe')",
           "; a outra versão ", if (.ee_dif_levene_degenerado(y, g, lev_outro, setdiff(c("media", "mediana"), centro_levene))$ruido) "também não tem valor interpretável (ajuste perfeito)" else
             paste0("dá ", .ee_junta(.ee_estat("F", lev_outro[["F"]], NULL, decimal), .pr_p(lev_outro[["p"]], decimal), decimal = decimal)), "."),
    paste0("d de Cohen como no JASP: DP combinado na linha de Student, média simples das variâncias na de Welch; IC pelo t não central. ",
           "g de Hedges com correção gama exata (IC = IC do d vezes J): corrige o viés do d em amostras pequenas (depende só dos gl, não do desbalanceamento); indice = 'g' o torna principal."))
  avisos <- character()
  if (rota == "classica" && degenerado)
    avisos <- c(avisos, paste0("Levene degenerado (", lev_deg$motivo, "): o teste de homogeneidade não informa nada com esses dados; a rota não foi decidida por ele (linha de Welch)."))
  if (gl_s < 20 && indice == "d") avisos <- c(avisos, paste0("Amostra pequena (gl = ", gl_s, "): o d superestima o efeito populacional em ",
                                                              round(100 * (1 / J - 1)), "%; o g de Hedges corrige esse viés (indice = 'g'; Lakens, 2013)."))
  if (min(n1, n2) < 30 && any(c(sw1[["p"]], sw2[["p"]]) < alfa, na.rm = TRUE))
    avisos <- c(avisos, "Shapiro-Wilk sugere afastamento da normalidade com grupo pequeno: conferir o Q-Q plot; se o afastamento for grave, teste_mann_whitney() ou bootstrap.")

  # frase APA
  linha <- if (usa_welch) 2 else 1
  tt <- tab[linha, ]
  e <- ef[princ, ]
  maior <- if (m1 >= m2) 1 else 2; menor <- 3 - maior
  ms <- c(m1, m2); ss <- c(s1, s2); nn <- c(n1, n2)
  uni <- alt != "two.sided"
  if (uni) {
    # unilateral (G-05): tudo com sinal na ordem dos níveis (nv1 - nv2), como $tabela e $efeito
    ic_dif <- c(tt$ic_inf, tt$ic_sup); ic_e <- c(e$ic_inf, e$ic_sup)
    t_val <- tt$t; dif_val <- m1 - m2; e_val <- e$valor
  } else {
    sinal <- if (maior == 1) 1 else -1
    ic_dif <- sort(sinal * c(tt$ic_inf, tt$ic_sup)); ic_e <- sort(sinal * c(e$ic_inf, e$ic_sup))
    t_val <- abs(tt$t); dif_val <- abs(m1 - m2); e_val <- abs(e$valor)
  }
  oposto <- uni && .ee_dif_oposto(alt, tt$t)
  simb <- if (princ == 3) "g" else if (princ == 4) "Δ" else "d"
  sig <- isTRUE(tt$p < alfa)
  # privacidade (H-02): com algum grupo de n < 5, M e DP saem da frase nos DOIS grupos. Só o
  # grupo pequeno não basta: M e DP do outro grupo, a diferença e o t devolvem por álgebra a M
  # e o DP do pequeno (M_c = M_b - dif; EP = dif/t; DP_c² = n_c(EP² - DP_b²/n_b)), e com
  # n = 2 os dois escores.
  peq <- any(.ee_dif_peq(nn))
  gtxt <- function(i, lab_m, lab_dp) if (peq) nv[i] else
    paste0(nv[i], " (", lab_m, " = ", .ee_dif_m(ms[i], decimal), .ee_sep(decimal), lab_dp, " = ", .ee_dif_m(ss[i], decimal), ")")
  peq_txt <- if (peq) paste0(" (", .ee_dif_omit_txt(idioma), if (idioma == "en") " for both groups" else " nos dois grupos", ")") else ""
  if (peq) avisos <- c(avisos, paste0(.ee_dif_aviso_peq(nv[.ee_dif_peq(nn)], "M e DP"),
                                      if (!all(.ee_dif_peq(nn))) paste0(" Saem também os do outro grupo: com eles, a diferença e o t devolveriam M e DP do grupo pequeno. ",
                                                                        "As estatísticas do teste no print, com as descritivas do outro grupo, ainda permitem essa conta; com grupo tão pequeno, juntar categorias ou não relatar o teste.") else ""))
  stats_txt <- .ee_junta(.ee_estat("t", t_val, tt$gl, decimal), .pr_p(tt$p, decimal), decimal = decimal)
  h1 <- if (!uni) "" else if (idioma == "en") paste0(", one-tailed (H1: mean of ", nv[1], if (alt == "greater") " greater" else " less", " than mean of ", nv[2], "),") else
    paste0(", unilateral (H1: média de ", nv[1], if (alt == "greater") " maior" else " menor", " que a de ", nv[2], "),")
  if (idioma == "en") {
    dif_lab <- if (uni) paste0("mean difference ", nv[1], " - ", nv[2], " = ") else "mean difference "
    apa <- paste0("An independent-samples ", if (usa_welch) "Welch" else "Student", " t test", h1, " ",
                  if (sig) "indicated that " else if (uni) "did not indicate a statistically significant difference in the direction of H1: " else "did not indicate a statistically significant difference: ",
                  gtxt(maior, "M", "SD"), if (sig) " scored higher than " else " and ", gtxt(menor, "M", "SD"), peq_txt, ", ", stats_txt,
                  "; ", dif_lab, .ee_dif_m(dif_val, decimal), sepa, .ee_ic(ic_dif[1], ic_dif[2], 2, TRUE, decimal, nivel, "en"),
                  sepa, .ee_efeito_txt(simb, e_val, ic_e[1], ic_e[2], FALSE, 2, decimal, "en", nivel = nivel),
                  " (", .ee_dif_rotulo(e$valor, "d", "en"), " effect", .ee_dif_txt_oposto(oposto, "en"), ").")
  } else {
    dif_lab <- if (uni) paste0("diferença média ", nv[1], " - ", nv[2], " = ") else "diferença média "
    apa <- paste0("Um teste t de ", if (usa_welch) "Welch" else "Student", " para amostras independentes", h1, " ",
                  if (sig) "indicou média maior em " else if (uni) "não indicou diferença estatisticamente significativa no sentido de H1 entre " else
                    "não indicou diferença estatisticamente significativa entre ",
                  gtxt(maior, "M", "DP"), if (sig) " do que em " else " e ", gtxt(menor, "M", "DP"), peq_txt, ", ", stats_txt,
                  "; ", dif_lab, .ee_dif_m(dif_val, decimal), sepa, .ee_ic(ic_dif[1], ic_dif[2], 2, TRUE, decimal, nivel, "pt"),
                  sepa, .ee_efeito_txt(simb, e_val, ic_e[1], ic_e[2], FALSE, 2, decimal, "pt", nivel = nivel),
                  " (efeito ", .ee_dif_rotulo(e$valor, "d", "pt"), .ee_dif_txt_oposto(oposto, "pt"), ").")
  }
  apa <- .ee_dif_unilateral(apa, alt, idioma)
  avisos <- c(avisos, .ee_aviso_nao_sig(tt$p, alfa, if (!uni) paste0("diferença média ", nv[1], " - ", nv[2], " = ", .ee_dif_m(m1 - m2, decimal), sepa,
                                                                      .ee_ic(tt$ic_inf, tt$ic_sup, 2, TRUE, decimal, nivel, "pt")),
                                        "diferença", "independente", decimal))
  .ee_resultado("Teste t de amostras independentes", n = n1 + n2, n_excluidos = dd$n_excluidos,
                tabela = tab, efeito = ef, pressupostos = pres, decisoes = dec, avisos = avisos, apa = apa,
                extra = list(descritivas = desc, levene_outro_centro = lev_outro))
}

# Teste t pareado. Diferença = depois - antes (positivo = aumento), como na leitura do
# livro; o JASP subtrai na ordem do par (antes - depois) e mostra o sinal
# oposto. d principal: d na métrica original, média das diferenças sobre o DP combinado
# das duas medidas, sqrt((DP1² + DP2²)/2) (Caixa 13.2; 0,47 no Weight Gain).
teste_t_pareado <- function(dados, antes, depois, alternativa = c("bilateral", "maior", "menor"),
                            alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  sepa <- .ee_sep(decimal)
  alt <- .ee_dif_alt(match.arg(alternativa))
  .pr_exige("effectsize")
  dd <- .ee_dados(dados, c(antes, depois), numericas = c(antes, depois))
  x1 <- dd$d[[antes]]; x2 <- dd$d[[depois]]; n <- length(x1)
  if (n < 3) stop("o teste pareado precisa de pelo menos 3 pares completos")
  dif <- x2 - x1
  if (.ee_dif_constante(dif)) stop("as diferenças ", depois, " - ", antes, " não variam (todos os pares mudaram o mesmo valor): o erro-padrão é zero e o t pareado não é calculável; descrever a mudança", call. = FALSE)
  tt <- stats::t.test(x2, x1, paired = TRUE, alternative = alt, conf.level = nivel)
  r <- stats::cor(x1, x2)
  md <- mean(dif); sdd <- stats::sd(dif); s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  dz <- md / sdd
  ic_dz <- .ee_dif_ic_nct(unname(tt$statistic), n - 1, 1 / sqrt(n), nivel, alt)
  rm_d <- function(met) {
    e <- suppressWarnings(suppressMessages(effectsize::repeated_measures_d(x2, x1, method = met, adjust = FALSE,
                                                                           ci = nivel, alternative = alt)))
    c(e[[1]], e$CI_low, e$CI_high)
  }
  # d na métrica original = média das diferenças / sqrt((DP1² + DP2²)/2): é o método "av"
  # do effectsize (média das variâncias), cujo IC leva o pareamento em conta. O d_av de
  # Lakens (2013) usa a média simples dos dois DPs: quase igual, sem IC aqui.
  dcomb <- rm_d("av"); drm <- rm_d("rm")
  dav_l <- md / ((s1 + s2) / 2)
  ef <- .ee_dif_efeito_df(
    c("d na métrica original (média das diferenças / DP combinado das duas medidas)",
      "d_av de Lakens (média das diferenças / média dos dois DPs)",
      "d_rm (corrigido pela correlação; 'Correct for correlation' do JASP)",
      "d_z (média das diferenças / DP das diferenças; padrão do JASP)"),
    c(dcomb[1], dav_l, drm[1], dz), c(dcomb[2], NA, drm[2], ic_dz[1]),
    c(dcomb[3], NA, drm[3], ic_dz[2]), c(TRUE, FALSE, FALSE, FALSE))
  swd <- .ee_dif_sw(dif)
  ct <- stats::cor.test(x1, x2)
  pres <- rbind(
    .ee_dif_pres("Normalidade da diferença (Shapiro-Wilk)", swd[["W"]], "", swd[["p"]], .ee_dif_le_normal(dif, "diferença")),
    .ee_dif_pres("Correlação entre as medidas (r de Pearson)", r, as.character(n - 2), ct$p.value,
                 if (r >= .7) "correlação alta: o d_z fica inflado; relatar o d na métrica original" else "informativo"))
  tab <- data.frame(teste = "t pareado", t = unname(tt$statistic), gl = unname(tt$parameter), p = tt$p.value,
                    dif_media = md, ic_inf = tt$conf.int[1], ic_sup = tt$conf.int[2], stringsAsFactors = FALSE)
  desc <- data.frame(medida = c(antes, depois, "diferença (depois - antes)"), n = n,
                     media = c(mean(x1), mean(x2), md), dp = c(s1, s2, sdd),
                     mediana = c(stats::median(x1), stats::median(x2), stats::median(dif)), stringsAsFactors = FALSE)
  dec <- c(
    paste0("Diferença = ", depois, " - ", antes, " (positivo = aumento), como a leitura do livro; o JASP mostra o t com o sinal oposto (ordem do par)."),
    "d principal: d na métrica original, com o DP combinado das duas medidas (escolha do livro, Caixa 13.2); é o comparável a estudos de grupos independentes.",
    paste0("O 'Correct for correlation' do JASP 0.95.1 devolve d_rm = d_z * sqrt(2(1 - r)) (Dunlap et al., 1996), não o d do livro: aqui d_rm = ",
           .pr_num(drm[1], 3, TRUE, decimal), " e d na métrica original = ", .pr_num(dcomb[1], 3, TRUE, decimal), ". d_z (", .pr_num(dz, 3, TRUE, decimal), ") infla com r alto e serve para cálculo de poder."),
    "Não há homogeneidade de variâncias a testar no pareado; a normalidade é a da diferença.",
    if (alt != "two.sided") paste0("Teste unilateral (H1: ", if (alt == "greater") "aumento" else "diminuição", " de ", antes, " para ", depois,
                                   "): a frase dá t, diferença, d e IC com sinal (", depois, " - ", antes, ") e diz quando os dados vão no sentido oposto ao de H1."))
  avisos <- character()
  if (n < 30 && isTRUE(swd[["p"]] < alfa))
    avisos <- "Shapiro-Wilk da diferença sugere afastamento da normalidade com n pequeno: ver o Q-Q plot; se grave, teste_wilcoxon_pareado()."
  sig <- isTRUE(tt$p.value < alfa)
  sub <- md >= 0
  uni <- alt != "two.sided"
  oposto <- uni && .ee_dif_oposto(alt, md)
  # bilateral: módulos e verbo (aumento/diminuição); unilateral (G-05): sinal de depois - antes, como $tabela
  t_val <- if (uni) unname(tt$statistic) else abs(tt$statistic)
  dif_val <- if (uni) md else abs(md)
  d_val <- if (uni) dcomb[1] else abs(dcomb[1])
  stats_txt <- .ee_junta(.ee_estat("t", t_val, n - 1, decimal), .pr_p(tt$p.value, decimal), decimal = decimal)
  ic <- if (sub || uni) tt$conf.int else -rev(tt$conf.int)
  ice <- if (sub || uni) dcomb[2:3] else -rev(dcomb[2:3])
  # privacidade (H-02): com n < 5 as duas medidas saem sem M e DP
  peq <- .ee_dif_peq(n)
  if (peq) avisos <- c(avisos, .ee_dif_aviso_peq(c(antes, depois), "M e DP"))
  mtxt <- function(x, s, lab_dp) if (peq) "" else paste0(" (M = ", .ee_dif_m(mean(x), decimal), .ee_sep(decimal), lab_dp, " = ", .ee_dif_m(s, decimal), ")")
  if (idioma == "en") {
    h1 <- if (uni) paste0(", one-tailed (H1: ", if (alt == "greater") "increase" else "decrease", " from ", antes, " to ", depois, "),") else ""
    apa <- paste0("A paired-samples t test", h1, " ", if (sig) paste0("indicated that the scores ", if (sub) "increased" else "decreased", " from ") else
                    if (uni) "did not indicate a statistically significant change in the direction of H1 from " else "did not indicate a statistically significant change from ",
                  antes, mtxt(x1, s1, "SD"), " to ", depois, mtxt(x2, s2, "SD"), if (peq) paste0(" (", .ee_dif_omit_txt("en"), ")") else "", ", ", stats_txt,
                  "; ", if (uni) paste0("mean difference ", depois, " - ", antes, " = ") else "mean difference ", .ee_dif_m(dif_val, decimal), sepa,
                  .ee_ic(ic[1], ic[2], 2, TRUE, decimal, nivel, "en"), sepa,
                  .ee_efeito_txt("d", d_val, ice[1], ice[2], FALSE, 2, decimal, "en", nivel = nivel),
                  " (original-metric d, pooled SD of both measures; ", .ee_dif_rotulo(dcomb[1], "d", "en"), " effect", .ee_dif_txt_oposto(oposto, "en"), ").")
  } else {
    h1 <- if (uni) paste0(", unilateral (H1: ", if (alt == "greater") "aumento" else "diminuição", " de ", antes, " para ", depois, "),") else ""
    apa <- paste0("Um teste t para amostras dependentes", h1, " ", if (sig) paste0("indicou que os escores ", if (sub) "aumentaram" else "diminuíram", " do momento ") else
                    if (uni) "não indicou mudança estatisticamente significativa no sentido de H1 do momento " else "não indicou mudança estatisticamente significativa do momento ",
                  antes, mtxt(x1, s1, "DP"), " para ", depois, mtxt(x2, s2, "DP"), if (peq) paste0(" (", .ee_dif_omit_txt("pt"), ")") else "", ", ", stats_txt,
                  "; ", if (uni) paste0("diferença média ", depois, " - ", antes, " = ") else "diferença média ", .ee_dif_m(dif_val, decimal), sepa,
                  .ee_ic(ic[1], ic[2], 2, TRUE, decimal, nivel, "pt"), sepa,
                  .ee_efeito_txt("d", d_val, ice[1], ice[2], FALSE, 2, decimal, "pt", nivel = nivel),
                  " (d na métrica original, DP combinado das duas medidas; efeito ", .ee_dif_rotulo(dcomb[1], "d", "pt"), .ee_dif_txt_oposto(oposto, "pt"), ").")
  }
  apa <- .ee_dif_unilateral(apa, alt, idioma)
  avisos <- c(avisos, .ee_aviso_nao_sig(tt$p.value, alfa, if (!uni) paste0("diferença média ", depois, " - ", antes, " = ", .ee_dif_m(md, decimal), sepa,
                                                                            .ee_ic(tt$conf.int[1], tt$conf.int[2], 2, TRUE, decimal, nivel, "pt")),
                                        "mudança", "pareado", decimal))
  .ee_resultado("Teste t pareado", n = n, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef,
                pressupostos = pres, decisoes = dec, avisos = avisos, apa = apa,
                extra = list(descritivas = desc, r = r))
}

# Teste t de uma amostra (contra um valor de referência mu), com o Wilcoxon de postos
# sinalizados de uma amostra ao lado (análogo não paramétrico usual).
teste_t_uma_amostra <- function(dados, vd, mu = 0, alternativa = c("bilateral", "maior", "menor"),
                                alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  sepa <- .ee_sep(decimal)
  alt <- .ee_dif_alt(match.arg(alternativa))
  dd <- .ee_dados(dados, vd, numericas = vd)
  x <- dd$d[[vd]]; n <- length(x)
  if (n < 3) stop("são necessários pelo menos 3 casos")
  .ee_dif_exige_variacao(x, vd)
  tt <- stats::t.test(x, mu = mu, alternative = alt, conf.level = nivel)
  d <- (mean(x) - mu) / stats::sd(x)
  ic <- .ee_dif_ic_nct(unname(tt$statistic), n - 1, 1 / sqrt(n), nivel, alt)
  J <- .ee_dif_J(n - 1)
  w <- .ee_dif_wsr(x - mu, alt, NULL, nivel); nd <- w$n
  maxw <- nd * (nd + 1) / 2
  rrb <- 2 * w$V / maxw - 1
  ic_r <- .ee_dif_ic_rrb(rrb, sqrt(nd * (nd + 1) * (2 * nd + 1) / 6) / maxw, nivel, alt)
  ef <- .ee_dif_efeito_df(c("d de Cohen ((M - mu)/DP)", "g de Hedges", "rank-biserial pareada (Wilcoxon)"),
                          c(d, d * J, rrb), c(ic[1], ic[1] * J, ic_r[1]), c(ic[2], ic[2] * J, ic_r[2]), c(TRUE, FALSE, FALSE))
  tab <- data.frame(teste = c("t de uma amostra", "Wilcoxon de uma amostra (V)"),
                    estatistica = c(unname(tt$statistic), w$V),
                    gl = c(n - 1, NA), p = c(tt$p.value, w$p), dif_media = c(mean(x) - mu, NA),
                    ic_inf = c(tt$conf.int[1] - mu, NA), ic_sup = c(tt$conf.int[2] - mu, NA), stringsAsFactors = FALSE)
  sw <- .ee_dif_sw(x)
  pres <- .ee_dif_pres("Normalidade (Shapiro-Wilk)", sw[["W"]], "", sw[["p"]], .ee_dif_le_normal(x, vd))
  sig <- isTRUE(tt$p.value < alfa)
  stats_txt <- .ee_junta(.ee_estat("t", tt$statistic, n - 1, decimal), .pr_p(tt$p.value, decimal), decimal = decimal)
  mu_txt <- .pr_num(mu, 2, TRUE, decimal)
  uni <- alt != "two.sided"
  oposto <- uni && .ee_dif_oposto(alt, tt$statistic)
  # privacidade (H-02): com n < 5, M e DP saem da frase (a diferença M - mu, que é a estimativa testada, fica)
  peq <- .ee_dif_peq(n)
  avisos <- if (peq) paste(.ee_dif_aviso_peq(vd, "M e DP"), "A diferença M - mu, que é a estimativa testada, continua na frase e permite recuperar M.") else character()
  if (idioma == "en") {
    h1 <- if (uni) paste0(", one-tailed (H1: mean of ", vd, if (alt == "greater") " greater" else " less", " than ", mu_txt, "),") else ""
    apa <- paste0("A one-sample t test", h1, " compared the mean of ", vd, " (", if (peq) .ee_dif_omit_txt("en") else paste0("M = ", .ee_dif_m(mean(x), decimal), .ee_sep(decimal), "SD = ",
                  .ee_dif_m(stats::sd(x), decimal)), ") with the reference value ", mu_txt, ", ", stats_txt, "; difference ",
                  .ee_dif_m(mean(x) - mu, decimal), sepa, .ee_ic(tt$conf.int[1] - mu, tt$conf.int[2] - mu, 2, TRUE, decimal, nivel, "en"), sepa,
                  .ee_efeito_txt("d", d, ic[1], ic[2], FALSE, 2, decimal, "en", nivel = nivel), " (", .ee_dif_rotulo(d, "d", "en"), " effect", .ee_dif_txt_oposto(oposto, "en"), ").")
  } else {
    h1 <- if (uni) paste0(", unilateral (H1: média de ", vd, if (alt == "greater") " maior" else " menor", " que ", mu_txt, "),") else ""
    apa <- paste0("Um teste t de uma amostra", h1, " comparou a média de ", vd, " (", if (peq) .ee_dif_omit_txt("pt") else paste0("M = ", .ee_dif_m(mean(x), decimal), .ee_sep(decimal), "DP = ",
                  .ee_dif_m(stats::sd(x), decimal)), ") com o valor de referência ", mu_txt, ", ", stats_txt,
                  if (sig) "" else " (sem diferença estatisticamente significativa)", "; diferença ",
                  .ee_dif_m(mean(x) - mu, decimal), sepa, .ee_ic(tt$conf.int[1] - mu, tt$conf.int[2] - mu, 2, TRUE, decimal, nivel, "pt"), sepa,
                  .ee_efeito_txt("d", d, ic[1], ic[2], FALSE, 2, decimal, "pt", nivel = nivel), " (efeito ", .ee_dif_rotulo(d, "d", "pt"), .ee_dif_txt_oposto(oposto, "pt"), ").")
  }
  apa <- .ee_dif_unilateral(apa, alt, idioma)
  avisos <- c(avisos, .ee_aviso_nao_sig(tt$p.value, alfa, if (!uni) paste0("diferença M - mu = ", .ee_dif_m(mean(x) - mu, decimal), sepa,
                                                                            .ee_ic(tt$conf.int[1] - mu, tt$conf.int[2] - mu, 2, TRUE, decimal, nivel, "pt")),
                                        "diferença em relação ao valor de referência", "uma_amostra", decimal))
  .ee_resultado("Teste t de uma amostra", n = n, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef,
                pressupostos = pres,
                decisoes = c(paste0("Diferença = média - ", mu, "; d mantém o sinal (direção) e a frase o declara.",
                                    if (uni) paste0(" Teste unilateral (H1: média ", if (alt == "greater") "maior" else "menor", " que ", mu, "): a frase diz quando os dados vão no sentido oposto ao de H1.") else ""),
                             paste0("O livro trata o t de uma amostra só em nível conceitual; o Wilcoxon de uma amostra sai ao lado como alternativa sem normalidade (",
                                    w$metodo, "; ", w$zeros, " valor(es) igual(is) a mu descartado(s), como no R 4.5.1 do JASP 0.95.1).")),
                avisos = avisos, apa = apa)
}

# ================================================================ equivalência (TOST)
# Dois testes unilaterais (TOST; Schuirmann, 1987, como apresentado em Lakens, 2017):
# H0 inferior Δ <= ΔL e H0 superior Δ >= ΔU. Há equivalência quando os dois são
# rejeitados, o que é o mesmo que o IC de 1 - 2α (90% com α = .05) cair dentro dos
# limites; relata-se o maior dos dois p (Lakens, 2017). O teste usual (H0: Δ = 0) sai ao
# lado, e os dois juntos dão os quatro resultados de Lakens (2017): equivalente e não
# diferente de zero; diferente de zero e equivalente; diferente e não equivalente;
# indeterminado.
#  tipo "independente": Welch por padrão (gl de Satterthwaite, eq. 3 e 4 de Lakens,
#    2017); welch = FALSE dá Student. Diferença = 1º nível - 2º nível, como
#    teste_t_independente();
#  "pareado": diferença = depois - antes, como teste_t_pareado() (eq. 6);
#  "uma_amostra": M - mu (eq. 7);
#  "correlacao": r de Pearson, z de Fisher com EP = 1/√(N - 3) contra a normal (eq. 8,
#    que Lakens, 2017, toma de Goertzen & Cribbie, 2010).
# limites: c(inferior, superior) ou um número só (simétrico, ±). Precisam vir de antes
# da análise, pelo menor efeito de interesse; a função não sugere valor.
# escala = "bruta" (unidade da variável, o recomendado) ou "d" (d de Cohen; d_z no
# pareado), convertido pelo DP da amostra: média das variâncias com Welch, DP combinado
# com Student, DP das diferenças no pareado. A documentação do TOSTER avisa que limites
# padronizados dão resultado enviesado; a função registra o aviso.
teste_equivalencia <- function(dados, tipo = c("independente", "pareado", "uma_amostra", "correlacao"), limites,
                               vd = NULL, grupo = NULL, antes = NULL, depois = NULL, x = NULL, y = NULL, mu = 0,
                               escala = c("bruta", "d"), welch = TRUE, alfa = .05, decimal = ".", idioma = "pt") {
  tipo <- match.arg(tipo); escala <- match.arg(escala)
  sepa <- .ee_sep(decimal)
  if (missing(limites) || !is.numeric(limites) || !length(limites) %in% 1:2 || anyNA(limites))
    stop("informe limites = c(inferior, superior), ou um número só para limites simétricos, fixados antes da análise ",
         "pelo menor efeito de interesse (Lakens, 2017)", call. = FALSE)
  if (length(limites) == 1) {
    if (limites <= 0) stop("com um número só, limites precisa ser positivo (vira -limites e +limites)", call. = FALSE)
    limites <- c(-limites, limites)
  }
  if (limites[1] >= limites[2]) stop("o limite inferior precisa ser menor que o superior", call. = FALSE)
  if (!(alfa > 0 && alfa < .5)) stop("alfa precisa estar entre 0 e .5", call. = FALSE)
  nivel_eq <- 1 - 2 * alfa
  precisa <- function(...) {
    a <- list(...)
    vaz <- names(a)[vapply(a, is.null, logical(1))]
    if (length(vaz)) stop("tipo = \"", tipo, "\" pede ", paste(vaz, collapse = " e "), call. = FALSE)
  }
  dec <- character(); avisos <- character()
  lim_in <- limites
  rot_lim <- if (escala == "d") (if (tipo == "pareado") "d_z" else "d") else "unidade da variável"

  if (tipo == "correlacao") {
    precisa(x = x, y = y)
    if (escala == "d") stop("na correlação os limites são valores de r; use escala = \"bruta\"", call. = FALSE)
    if (any(abs(limites) >= 1)) stop("na correlação os limites precisam estar entre -1 e 1", call. = FALSE)
    dd <- .ee_dados(dados, c(x, y), numericas = c(x, y))
    a <- dd$d[[x]]; b <- dd$d[[y]]; n <- length(a)
    if (n < 5) stop("são necessários pelo menos 5 pares completos")
    .ee_dif_exige_variacao(a, x); .ee_dif_exige_variacao(b, y)
    ct <- stats::cor.test(a, b)
    est <- unname(ct$estimate); ep <- 1 / sqrt(n - 3); zr <- atanh(est)
    est_L <- (zr - atanh(limites[1])) / ep; est_U <- (zr - atanh(limites[2])) / ep
    p_L <- stats::pnorm(est_L, lower.tail = FALSE); p_U <- stats::pnorm(est_U)
    q <- stats::qnorm(1 - alfa)
    ic_eq <- tanh(zr + c(-1, 1) * q * ep)
    gl <- NA_real_; simb <- "z"
    est0 <- unname(ct$statistic); gl0 <- unname(ct$parameter); p0 <- ct$p.value; simb0 <- "t"
    ic0 <- ct$conf.int[1:2]
    largura_ic <- 2 * q * ep; largura_lim <- atanh(limites[2]) - atanh(limites[1])
    rot_est <- "r"; casas <- 3; zero <- FALSE
    obj <- paste0("a correlação entre ", x, " e ", y); obj_en <- paste0("the correlation between ", x, " and ", y)
    dec <- c(dec, paste0("Correlação de Pearson entre ", x, " e ", y, " (deleção listwise, n = ", n, "); TOST pela z de Fisher com EP = 1/√(n - 3), ",
                         "como em Lakens (2017, eq. 8, de Goertzen & Cribbie, 2010). O teste usual ao lado é o t do cor.test, o mesmo de correlacoes()."))
    n_exc <- dd$n_excluidos
  } else {
    if (tipo == "independente") {
      precisa(vd = vd, grupo = grupo)
      dd <- .ee_dados(dados, c(vd, grupo), numericas = vd)
      yv <- dd$d[[vd]]; g <- .ee_fator(dd$d[[grupo]], grupo)
      if (nlevels(g) != 2) stop("o teste de equivalência independente compara exatamente 2 grupos; ", grupo, " tem ", nlevels(g), " níveis com dados")
      nv <- levels(g); y1 <- yv[g == nv[1]]; y2 <- yv[g == nv[2]]
      n1 <- length(y1); n2 <- length(y2); n <- n1 + n2
      if (min(n1, n2) < 2) stop("cada grupo precisa de pelo menos 2 casos")
      .ee_dif_exige_variacao(yv, vd, g, grupo)
      s1 <- stats::sd(y1); s2 <- stats::sd(y2)
      est <- mean(y1) - mean(y2)
      if (welch) {
        ep <- sqrt(s1^2 / n1 + s2^2 / n2)
        gl <- ep^4 / ((s1^2 / n1)^2 / (n1 - 1) + (s2^2 / n2)^2 / (n2 - 1))
        dp_ref <- sqrt((s1^2 + s2^2) / 2)
      } else {
        sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
        ep <- sp * sqrt(1 / n1 + 1 / n2); gl <- n1 + n2 - 2; dp_ref <- sp
      }
      rot_est <- paste0("diferença média ", nv[1], " - ", nv[2])
      obj <- paste0("a diferença entre ", nv[1], " e ", nv[2]); obj_en <- paste0("the difference between ", nv[1], " and ", nv[2])
      dec <- c(dec, paste0("Diferença = ", nv[1], " - ", nv[2], " (ordem dos níveis, como teste_t_independente()); ",
                           if (welch) "t de Welch com gl de Satterthwaite, padrão de Lakens (2017, eq. 3 e 4; Delacre et al., 2017)." else
                             "t de Student (welch = FALSE), que supõe variâncias iguais; Lakens (2017) recomenda Welch por padrão."))
    } else if (tipo == "pareado") {
      precisa(antes = antes, depois = depois)
      dd <- .ee_dados(dados, c(antes, depois), numericas = c(antes, depois))
      dif <- dd$d[[depois]] - dd$d[[antes]]; n <- length(dif)
      if (n < 3) stop("o teste pareado precisa de pelo menos 3 pares completos")
      if (.ee_dif_constante(dif)) stop("as diferenças ", depois, " - ", antes, " não variam: o erro-padrão é zero", call. = FALSE)
      est <- mean(dif); dp_ref <- stats::sd(dif); ep <- dp_ref / sqrt(n); gl <- n - 1
      rot_est <- paste0("diferença média ", depois, " - ", antes)
      obj <- paste0("a mudança de ", antes, " para ", depois); obj_en <- paste0("the change from ", antes, " to ", depois)
      dec <- c(dec, paste0("Diferença = ", depois, " - ", antes, " (positivo = aumento), como teste_t_pareado(); EP = DP das diferenças/√n (Lakens, 2017, eq. 6)."))
    } else {
      precisa(vd = vd)
      dd <- .ee_dados(dados, vd, numericas = vd)
      xv <- dd$d[[vd]]; n <- length(xv)
      if (n < 3) stop("são necessários pelo menos 3 casos")
      .ee_dif_exige_variacao(xv, vd)
      est <- mean(xv) - mu; dp_ref <- stats::sd(xv); ep <- dp_ref / sqrt(n); gl <- n - 1
      rot_est <- "diferença M - mu"
      obj <- paste0("a diferença entre a média de ", vd, " e o valor de referência ", .pr_num(mu, 2, TRUE, decimal))
      obj_en <- paste0("the difference between the mean of ", vd, " and the reference value ", .pr_num(mu, 2, TRUE, decimal))
      dec <- c(dec, paste0("Diferença = média de ", vd, " - ", mu, " (Lakens, 2017, eq. 7)."))
    }
    n_exc <- dd$n_excluidos
    if (escala == "d") {
      limites <- lim_in * dp_ref
      avisos <- c(avisos, paste0("Limites dados em ", rot_lim, " e convertidos pelo DP da amostra (", .pr_num(dp_ref, 3, TRUE, decimal), "): ",
                                 .pr_num(limites[1], 3, TRUE, decimal), " a ", .pr_num(limites[2], 3, TRUE, decimal),
                                 " na unidade da variável. Assim o limite passa a depender dos dados, e a documentação do TOSTER avisa que limites padronizados ",
                                 "dão resultado enviesado; preferir limites na unidade da escala, fixados antes."))
    }
    est_L <- (est - limites[1]) / ep; est_U <- (est - limites[2]) / ep
    p_L <- stats::pt(est_L, gl, lower.tail = FALSE); p_U <- stats::pt(est_U, gl)
    q <- stats::qt(1 - alfa, gl)
    ic_eq <- est + c(-1, 1) * q * ep
    est0 <- est / ep; gl0 <- gl; p0 <- 2 * stats::pt(-abs(est0), gl); simb <- "t"; simb0 <- "t"
    ic0 <- est + c(-1, 1) * stats::qt(1 - alfa / 2, gl) * ep
    largura_ic <- 2 * q * ep; largura_lim <- limites[2] - limites[1]
    casas <- 2; zero <- TRUE
  }

  p_tost <- max(p_L, p_U); qual <- if (p_L >= p_U) "inferior" else "superior"
  est_tost <- if (qual == "inferior") est_L else est_U
  equiv <- p_tost < alfa; difere <- p0 < alfa
  resultado <- if (equiv && !difere) "equivalente e não diferente de zero" else
    if (equiv && difere) "diferente de zero, mas dentro dos limites de equivalência" else
      if (!equiv && difere) "diferente de zero e não equivalente" else
        "indeterminado: nem diferente de zero nem equivalente"
  num <- function(v) .pr_num(v, casas, zero, decimal)
  lim_txt <- paste0(num(limites[1]), " a ", num(limites[2]))
  par_s <- if (tipo == "correlacao") "ρ" else "Δ"
  tab <- data.frame(teste = c(paste0("H0 inferior: ", par_s, " <= ", num(limites[1])), paste0("H0 superior: ", par_s, " >= ", num(limites[2])),
                              paste0("Teste usual: ", par_s, " = 0")),
                    estatistica = c(est_L, est_U, est0), gl = c(gl, gl, gl0), p = c(p_L, p_U, p0),
                    stringsAsFactors = FALSE)
  efe <- data.frame(indice = c(rot_est, rot_est), valor = est, nivel = c(nivel_eq, 1 - alfa),
                    ic_inf = c(ic_eq[1], ic0[1]), ic_sup = c(ic_eq[2], ic0[2]), principal = c(TRUE, FALSE))
  if (tipo != "correlacao" && escala == "d") efe$limites_d <- paste(lim_in, collapse = " a ")
  dec <- c(dec,
           paste0("Limites de equivalência: ", lim_txt, if (tipo == "correlacao") " (valores de r)" else if (escala == "d")
             paste0(" na unidade da variável, convertidos de ", rot_lim, " = ", paste(lim_in, collapse = " a ")) else " na unidade da variável",
             ". A escolha dos limites é a parte que precisa de justificativa (menor efeito de interesse, fixado antes de ver os dados; Lakens, 2017); a função não sugere valor."),
           paste0("TOST com α = ", .pr_num(alfa, 2, FALSE, decimal), " em cada teste unilateral; relatado o maior dos dois p (aqui o do limite ", qual,
                  "), o que equivale a conferir se o IC ", round(100 * nivel_eq), "% fica dentro dos limites (Lakens, 2017)."),
           paste0("Resultado combinado com o teste usual (Lakens, 2017): ", resultado, "."))
  if (largura_ic >= largura_lim)
    avisos <- c(avisos, paste0("O IC ", round(100 * nivel_eq), "% é mais largo que o intervalo de equivalência: com este n o teste não tinha como concluir equivalência, ",
                               "qualquer que fosse a estimativa (Lakens, 2017). Não concluir nada sobre igualdade; planejar n maior."))
  if (!equiv && !difere)
    avisos <- c(avisos, paste0("Nem diferença nem equivalência: os dados não bastam para nenhuma das duas conclusões; não escrever que ",
                               switch(tipo, independente = "os grupos são iguais", pareado = "não houve mudança",
                                      uma_amostra = "a média é igual ao valor de referência", correlacao = "não há associação"), "."))

  # frase
  ic_txt <- .ee_ic(ic_eq[1], ic_eq[2], casas, zero, decimal, nivel_eq, idioma)
  st_tost <- .ee_junta(if (simb == "z") paste0("z = ", .pr_num(est_tost, 2, TRUE, decimal)) else .ee_estat("t", est_tost, gl, decimal),
                       .pr_p(p_tost, decimal), decimal = decimal)
  st0 <- .ee_junta(if (tipo == "correlacao") paste0("r(", gl0, ") = ", num(est)) else .ee_estat("t", est0, gl0, decimal), .pr_p(p0, decimal), decimal = decimal)
  rota <- if (tipo == "correlacao") "z de Fisher" else if (tipo == "independente") (if (welch) "Welch" else "Student") else if (tipo == "pareado") "pareado" else "uma amostra"
  rota_en <- if (tipo == "correlacao") "Fisher's z" else if (tipo == "independente") (if (welch) "Welch" else "Student") else if (tipo == "pareado") "paired" else "one-sample"
  fora <- if (ic_eq[1] <= limites[1] && ic_eq[2] >= limites[2]) "ambos os limites" else if (ic_eq[2] >= limites[2]) "o limite superior" else "o limite inferior"
  if (idioma == "en") {
    fora_en <- c("ambos os limites" = "both bounds", "o limite superior" = "the upper bound", "o limite inferior" = "the lower bound")[[fora]]
    rot_en <- if (tipo == "correlacao") "r" else sub("^diferença média ", "mean difference ", sub("^diferença M - mu", "difference M - mu", rot_est))
    est_en <- paste0(rot_en, " = ", num(est), sepa, ic_txt)
    lim_en <- paste0(num(limites[1]), " to ", num(limites[2]),
                     if (tipo != "correlacao" && escala == "d") paste0(" (", rot_lim, " from ", paste(lim_in, collapse = " to "), ")") else "")
    apa <- paste0("An equivalence test (TOST, ", rota_en, ") with bounds of ", lim_en,
                  if (equiv) paste0(" indicated equivalence: ", obj_en, " fell within the bounds, ") else
                    paste0(" did not indicate equivalence: the CI of ", obj_en, " crosses ", fora_en, ", "),
                  est_en, sepa, st_tost, " (larger of the two one-sided tests). The usual test ",
                  if (difere) "indicated a statistically significant " else "did not indicate a statistically significant ",
                  if (tipo == "correlacao") "association, " else "difference from zero, ", st0, ".")
  } else {
    est_pt <- paste0(if (tipo == "correlacao") "r" else rot_est, " = ", num(est), sepa, ic_txt)
    lim_pt <- paste0(num(limites[1]), " a ", num(limites[2]),
                     if (tipo != "correlacao" && escala == "d") paste0(" (", rot_lim, " de ", paste(lim_in, collapse = " a "), ")") else "")
    apa <- paste0("Um teste de equivalência (TOST, ", rota, ") com limites de ", lim_pt,
                  if (equiv) paste0(" indicou equivalência: ", obj, " ficou dentro dos limites, ") else
                    paste0(" não indicou equivalência: o IC da ", sub("^a ", "", obj), " ultrapassa ", fora, ", "),
                  est_pt, sepa, st_tost, " (o maior dos dois testes unilaterais). O teste usual ",
                  if (difere) "indicou " else "não indicou ",
                  if (tipo == "correlacao") "associação estatisticamente significativa, " else "diferença estatisticamente significativa de zero, ", st0, ".")
  }
  .ee_resultado(paste0("Teste de equivalência (TOST, ", rota, ")"), n = n, n_excluidos = n_exc, tabela = tab, efeito = efe,
                decisoes = dec, avisos = avisos, apa = apa,
                extra = list(limites = limites, equivalente = equiv, diferente_de_zero = difere, resultado = resultado, p_tost = p_tost))
}

# ================================================================ não paramétricos (2)
# Mann-Whitney. U = W do R = U do 1º grupo (o que o JASP mostra). Rank-biserial:
#  sinal_rrb = "direcional": r = 2U/(n1 n2) - 1, positivo quando o 1º grupo tende a
#  valores maiores (effectsize e JASP >= 0.96);
#  "jasp_0951": r = 1 - 2U/(n1 n2), o sinal invertido do JASP 0.95.1 (corrigido depois no
#  próprio JASP, PR #392); use só para conferir uma tela antiga.
teste_mann_whitney <- function(dados, vd, grupo, sinal_rrb = c("direcional", "jasp_0951"),
                               regua_rrb = c("r", "grissom_kim"), alternativa = c("bilateral", "maior", "menor"),
                               exato = NULL, alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  sepa <- .ee_sep(decimal)
  sinal_rrb <- match.arg(sinal_rrb); regua_rrb <- match.arg(regua_rrb)
  alt <- .ee_dif_alt(match.arg(alternativa))
  dd <- .ee_dados(dados, c(vd, grupo), numericas = vd)
  y <- dd$d[[vd]]; g <- .ee_fator(dd$d[[grupo]], grupo)
  if (nlevels(g) != 2) stop("o Mann-Whitney compara 2 grupos; para 3 ou mais, use teste_kruskal_wallis()")
  nv <- levels(g); y1 <- y[g == nv[1]]; y2 <- y[g == nv[2]]; n1 <- length(y1); n2 <- length(y2); N <- n1 + n2
  if (min(n1, n2) < 2) stop("cada grupo precisa de pelo menos 2 casos")
  # constante dentro de todos os grupos: o IC de Hodges-Lehmann do wilcox.test quebra com erro cru
  .ee_dif_exige_variacao(y, vd, g, grupo)
  ex <- .ee_dif_exato(exato, max(n1, n2), anyDuplicated(y) > 0)
  w <- suppressWarnings(stats::wilcox.test(y1, y2, alternative = alt, conf.int = TRUE, conf.level = nivel, exact = ex, correct = TRUE))
  U1 <- unname(w$statistic)
  rrb <- 2 * U1 / (n1 * n2) - 1
  ic <- .ee_dif_ic_rrb(rrb, sqrt((n1 + n2 + 1) / (3 * n1 * n2)), nivel, alt)
  if (sinal_rrb == "jasp_0951") { rrb <- -rrb; ic <- -rev(ic) }
  rk <- rank(y); emp <- table(rk)
  sig_u <- sqrt(n1 * n2 / 12 * ((N + 1) - sum(emp^3 - emp) / (N * (N - 1))))
  z <- (U1 - n1 * n2 / 2) / sig_u
  tab <- data.frame(teste = "Mann-Whitney", U = U1, z = z, p = w$p.value,
                    exato = grepl("exact", w$method),
                    hodges_lehmann = unname(w$estimate), ic_inf = w$conf.int[1], ic_sup = w$conf.int[2],
                    stringsAsFactors = FALSE)
  desc <- data.frame(grupo = nv, n = c(n1, n2), mediana = c(stats::median(y1), stats::median(y2)),
                     media = c(mean(y1), mean(y2)), dp = c(stats::sd(y1), stats::sd(y2)),
                     media_postos = c(mean(rk[g == nv[1]]), mean(rk[g == nv[2]])),
                     soma_postos = c(sum(rk[g == nv[1]]), sum(rk[g == nv[2]])), stringsAsFactors = FALSE)
  ef <- .ee_dif_efeito_df("correlação rank-biserial", rrb, ic[1], ic[2], TRUE)
  dec <- c(paste0("U = W do R = U do grupo ", nv[1], " (1º nível), o mesmo que o JASP mostra; p ", if (ex) "exato" else "pela aproximação normal com correção de continuidade e de empates",
                  " (regra do JASP 0.95.1: exato só com n < 50 por grupo e sem empates; o R 4.6 mudou o padrão do wilcox.test para p exato condicional com empates; exato = TRUE o pede)."),
           if (sinal_rrb == "direcional") paste0("r_rb positivo = ", nv[1], " tende a valores maiores que ", nv[2], ". O JASP 0.95.1 mostra o mesmo valor com sinal trocado (1 - 2U/(n1 n2)); o JASP 0.96 já corrigiu.")
           else "r_rb no sinal do JASP 0.95.1 (1 - 2U/(n1 n2)): positivo quando o 1º grupo tende a valores MENORES.",
           paste0("Rótulo do r_rb pela ", if (regua_rrb == "r") "convenção de r (.10/.30/.50), a do JASP e do effectsize; os cortes de Grissom & Kim citados no curso (.125/.304/.465) não puderam ser conferidos na fonte" else "régua de Grissom & Kim citada no curso (.125/.304/.465; fonte não conferida)", "."),
           "Mann-Whitney testa dominância estocástica, não medianas: relatar n, mediana, média e DP e média dos postos.",
           if (alt != "two.sided") paste0("Teste unilateral (H1: valores de ", nv[1], " tendem a ser ", if (alt == "greater") "maiores" else "menores", " que os de ", nv[2],
                                          "): z e r_rb com sinal; a frase diz quando os dados vão no sentido oposto ao de H1."))
  rot <- .ee_dif_rotulo(rrb, if (regua_rrb == "r") "r" else "grissom_kim", idioma)
  maior <- if (desc$media_postos[1] >= desc$media_postos[2]) 1 else 2
  st_txt <- .ee_junta(paste0("U = ", .ee_dif_u(U1, decimal)), paste0("z = ", .pr_num(z, 2, TRUE, decimal)), .pr_p(w$p.value, decimal), decimal = decimal)
  uni <- alt != "two.sided"
  oposto <- uni && .ee_dif_oposto(alt, z)  # z sempre no sinal direcional (positivo = 1º grupo maior)
  # privacidade (H-02): grupo com n < 5 sai sem mediana
  avisos <- character()
  if (any(.ee_dif_peq(desc$n))) avisos <- .ee_dif_aviso_peq(nv[.ee_dif_peq(desc$n)], "mediana")
  mdn <- function(i, lab) paste0(" (", if (.ee_dif_peq(desc$n[i])) .ee_dif_omit_txt(idioma) else paste0(lab, " = ", .ee_dif_m(desc$mediana[i], decimal)), ")")
  if (idioma == "en") {
    h1 <- if (uni) paste0(", one-tailed (H1: values of ", nv[1], " tend to be ", if (alt == "greater") "greater" else "smaller", " than those of ", nv[2], "),") else ""
    apa <- paste0("A Mann-Whitney test", h1, " compared ", vd, " between ", nv[1], mdn(1, "Mdn"), " and ", nv[2],
                  mdn(2, "Mdn"), ", ", st_txt, sepa, .ee_efeito_txt("r_rb", rrb, ic[1], ic[2], TRUE, 2, decimal, "en", nivel = nivel),
                  " (", rot, .ee_dif_txt_oposto(oposto, "en"), "); higher mean rank in ", nv[maior], ".")
  } else {
    h1 <- if (uni) paste0(", unilateral (H1: valores de ", nv[1], " tendem a ser ", if (alt == "greater") "maiores" else "menores", " que os de ", nv[2], "),") else ""
    apa <- paste0("Um teste de Mann-Whitney", h1, " comparou ", vd, " entre ", nv[1], mdn(1, "Md"), " e ", nv[2],
                  mdn(2, "Md"), ", ", st_txt, sepa, .ee_efeito_txt("r_rb", rrb, ic[1], ic[2], TRUE, 2, decimal, "pt", nivel = nivel),
                  " (", rot, .ee_dif_txt_oposto(oposto, "pt"), "); média de postos maior em ", nv[maior], ".")
  }
  apa <- .ee_dif_unilateral(apa, alt, idioma)
  avisos <- c(avisos, .ee_aviso_nao_sig(w$p.value, alfa, if (alt == "two.sided") paste0("correlação rank-biserial = ", .pr_num(rrb, 2, FALSE, decimal), sepa,
                                                                                         .ee_ic(ic[1], ic[2], 2, FALSE, decimal, nivel, "pt")),
                                        "diferença", NULL, decimal))
  .ee_resultado("Teste U de Mann-Whitney", n = N, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef,
                decisoes = dec, avisos = avisos, apa = apa, extra = list(descritivas = desc))
}

# Wilcoxon de postos sinalizados para duas medidas do mesmo grupo. Diferença =
# depois - antes; zeros descartados (método de Wilcoxon, como o JASP e o R).
# rank-biserial pareada = 2V/maxw - 1 (Kerby, 2014), positiva = aumento. O JASP 0.95.1
# subtrai na ordem do par (antes - depois): mesmo valor, sinal oposto
# (sinal_rrb = "jasp_0951" reproduz).
teste_wilcoxon_pareado <- function(dados, antes, depois, sinal_rrb = c("direcional", "jasp_0951"),
                                   regua_rrb = c("r", "grissom_kim"), alternativa = c("bilateral", "maior", "menor"),
                                   exato = NULL, alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  sepa <- .ee_sep(decimal)
  sinal_rrb <- match.arg(sinal_rrb); regua_rrb <- match.arg(regua_rrb)
  alt <- .ee_dif_alt(match.arg(alternativa))
  dd <- .ee_dados(dados, c(antes, depois), numericas = c(antes, depois))
  x1 <- dd$d[[antes]]; x2 <- dd$d[[depois]]; n <- length(x1)
  dif <- x2 - x1; nd <- sum(dif != 0)
  if (nd < 2) stop("são necessárias pelo menos 2 diferenças não nulas")
  w <- .ee_dif_wsr(dif, alt, exato, nivel)
  V <- w$V; maxw <- nd * (nd + 1) / 2
  V0 <- sum(rank(abs(dif[dif != 0]))[dif[dif != 0] > 0])  # V de Wilcoxon (sem zeros), base de r_rb e z
  rrb <- 2 * V0 / maxw - 1
  ic <- .ee_dif_ic_rrb(rrb, sqrt(nd * (nd + 1) * (2 * nd + 1) / 6) / maxw, nivel, alt)
  if (sinal_rrb == "jasp_0951") { rrb <- -rrb; ic <- -rev(ic) }
  z <- (V0 - nd * (nd + 1) / 4) / sqrt(nd * (nd + 1) * (2 * nd + 1) / 24)
  tab <- data.frame(teste = "Wilcoxon pareado", V = V, z = z, p = w$p, exato = w$exato, n_dif_nao_nulas = nd,
                    hodges_lehmann = w$hl, ic_inf = w$ic[1], ic_sup = w$ic[2], stringsAsFactors = FALSE)
  desc <- data.frame(medida = c(antes, depois, "diferença (depois - antes)"), n = n,
                     mediana = c(stats::median(x1), stats::median(x2), stats::median(dif)),
                     media = c(mean(x1), mean(x2), mean(dif)), dp = c(stats::sd(x1), stats::sd(x2), stats::sd(dif)),
                     stringsAsFactors = FALSE)
  ef <- .ee_dif_efeito_df("correlação rank-biserial pareada", rrb, ic[1], ic[2], TRUE)
  rot <- .ee_dif_rotulo(rrb, if (regua_rrb == "r") "r" else "grissom_kim", idioma)
  dec <- c(paste0("Diferença = ", depois, " - ", antes, "; V = soma dos postos das diferenças positivas; ", n - nd, " diferença(s) nula(s) descartada(s)."),
           if (sinal_rrb == "direcional") "r_rb positivo = aumento. O JASP 0.95.1 (ordem antes, depois) mostra o mesmo valor com sinal oposto."
           else "r_rb no sinal do JASP 0.95.1 (antes - depois): positivo = diminuição.",
           paste0("z sem correção de empates, como o JASP; p: ", w$metodo, ". Regra do R 4.5.1/JASP 0.95.1: zeros fora antes de ranquear; exato só com n < 50, sem empates e sem zeros. ",
                  "O R 4.6 mudou o wilcox.test (exato condicional, zeros nos postos); exato = TRUE pede essa versão, que muda V e p."),
           "A r_rb pareada pondera as diferenças pelos seus postos (Kerby, 2014): não é a proporção de pessoas que mudaram.",
           if (alt != "two.sided") paste0("Teste unilateral (H1: ", if (alt == "greater") "aumento" else "diminuição", " de ", antes, " para ", depois,
                                          "): z e r_rb com sinal; a frase diz quando os dados vão no sentido oposto ao de H1."))
  st_txt <- .ee_junta(paste0("V = ", .ee_dif_u(V, decimal)), paste0("z = ", .pr_num(z, 2, TRUE, decimal)), .pr_p(w$p, decimal), decimal = decimal)
  uni <- alt != "two.sided"
  oposto <- uni && .ee_dif_oposto(alt, z)  # z no sinal de depois - antes (positivo = aumento)
  # privacidade (H-02): com n < 5 as duas medidas saem sem mediana
  peq <- .ee_dif_peq(n)
  avisos <- if (peq) .ee_dif_aviso_peq(c(antes, depois), "mediana") else character()
  mdn <- function(i, lab) if (peq) "" else paste0(" (", lab, " = ", .ee_dif_m(desc$mediana[i], decimal), ")")
  if (idioma == "en") {
    h1 <- if (uni) paste0(", one-tailed (H1: ", if (alt == "greater") "increase" else "decrease", " from ", antes, " to ", depois, "),") else ""
    apa <- paste0("A Wilcoxon signed-rank test", h1, " compared ", antes, mdn(1, "Mdn"), " and ", depois,
                  mdn(2, "Mdn"), if (peq) paste0(" (", .ee_dif_omit_txt("en"), ")") else "", ", ", st_txt, sepa,
                  .ee_efeito_txt("r_rb", rrb, ic[1], ic[2], TRUE, 2, decimal, "en", nivel = nivel), " (", rot, .ee_dif_txt_oposto(oposto, "en"), ").")
  } else {
    h1 <- if (uni) paste0(", unilateral (H1: ", if (alt == "greater") "aumento" else "diminuição", " de ", antes, " para ", depois, "),") else ""
    apa <- paste0("Um teste de Wilcoxon de postos sinalizados", h1, " comparou ", antes, mdn(1, "Md"), " e ", depois,
                  mdn(2, "Md"), if (peq) paste0(" (", .ee_dif_omit_txt("pt"), ")") else "", ", ", st_txt, sepa,
                  .ee_efeito_txt("r_rb", rrb, ic[1], ic[2], TRUE, 2, decimal, "pt", nivel = nivel), " (", rot, .ee_dif_txt_oposto(oposto, "pt"), ").")
  }
  apa <- .ee_dif_unilateral(apa, alt, idioma)
  avisos <- c(avisos, .ee_aviso_nao_sig(w$p, alfa, if (alt == "two.sided") paste0("correlação rank-biserial pareada = ", .pr_num(rrb, 2, FALSE, decimal), sepa,
                                                                                   .ee_ic(ic[1], ic[2], 2, FALSE, decimal, nivel, "pt")),
                                        "mudança", NULL, decimal))
  .ee_resultado("Teste de Wilcoxon de postos sinalizados", n = n, n_excluidos = dd$n_excluidos, tabela = tab,
                efeito = ef, decisoes = dec, avisos = avisos, apa = apa, extra = list(descritivas = desc))
}

# ================================================================ ANOVA: auxiliares
# Brown-Forsythe F* (teste de médias; Brown & Forsythe, 1974), a fórmula do
# onewaytests::bf.test que o JASP usa.
.ee_dif_bf <- function(y, g) {
  n <- tapply(y, g, length); m <- tapply(y, g, mean); v <- tapply(y, g, stats::var); N <- sum(n)
  w <- (1 - n / N) * v
  Fs <- sum(n * (m - mean(y))^2) / sum(w)
  cc <- w / sum(w)
  gl2 <- 1 / sum(cc^2 / (n - 1))
  c(F = Fs, gl1 = length(n) - 1, gl2 = gl2, p = stats::pf(Fs, length(n) - 1, gl2, lower.tail = FALSE))
}

# Games-Howell (a fórmula manual do JASP e do rstatix): t de Welch por par, gl de
# Welch-Satterthwaite, p pela distribuição da amplitude estudentizada.
.ee_dif_games_howell <- function(y, g, nivel = .95) {
  nv <- levels(g); k <- length(nv)
  n <- tapply(y, g, length); m <- tapply(y, g, mean); v <- tapply(y, g, stats::var)
  pr <- utils::combn(k, 2)
  do.call(rbind, lapply(seq_len(ncol(pr)), function(j) {
    a <- pr[1, j]; b <- pr[2, j]
    va <- v[[a]] / n[[a]]; vb <- v[[b]] / n[[b]]
    se <- sqrt(va + vb); dif <- m[[a]] - m[[b]]
    gl <- (va + vb)^2 / (va^2 / (n[[a]] - 1) + vb^2 / (n[[b]] - 1))
    q <- stats::qtukey(nivel, k, gl)
    data.frame(comparacao = paste(nv[a], "-", nv[b]), dif = dif, ep = se, t = dif / se, gl = gl,
               p_ajustado = stats::ptukey(abs(dif / se) * sqrt(2), k, gl, lower.tail = FALSE),
               ic_inf = dif - q * se / sqrt(2), ic_sup = dif + q * se / sqrt(2), stringsAsFactors = FALSE)
  }))
}

# Comparações par a par por emmeans (Tukey, Bonferroni, Holm). Com Holm, o emmeans dá
# ICs de Bonferroni (comportamento documentado do pacote).
.ee_dif_pares_emm <- function(em, ajuste, nivel, traduz = NULL) {
  if (is.null(traduz)) traduz <- function(v, x) x
  ct <- emmeans::contrast(em, method = "pairwise")
  s <- as.data.frame(summary(ct, adjust = ajuste, infer = c(TRUE, TRUE), level = nivel))
  nm <- names(s)
  # rótulos montados dos níveis (o emmeans prefixa o nome da variável quando os níveis
  # parecem números); a ordem dos pares é a de combn dentro de cada grupo 'by'
  lv <- traduz(em@misc$pri.vars[1], em@levels[[em@misc$pri.vars[1]]])
  pr <- utils::combn(length(lv), 2)
  rot <- paste(lv[pr[1, ]], "-", lv[pr[2, ]])
  out <- data.frame(comparacao = rep(rot, length.out = nrow(s)), dif = s$estimate, ep = s$SE, t = s$t.ratio, gl = s$df,
                    p_ajustado = s$p.value, ic_inf = s[[grep("lower|asymp.LCL", nm)[1]]],
                    ic_sup = s[[grep("upper|asymp.UCL", nm)[1]]], stringsAsFactors = FALSE)
  by <- setdiff(nm, c("contrast", "estimate", "SE", "df", "lower.CL", "upper.CL", "t.ratio", "p.value",
                      "asymp.LCL", "asymp.UCL", "z.ratio"))
  if (length(by)) {
    bt <- s[, by, drop = FALSE]
    for (b in by) bt[[b]] <- traduz(b, trimws(as.character(bt[[b]])))
    out <- cbind(em_by = apply(bt, 1, function(r) paste(paste0(by, " = ", r), collapse = ", ")), out, stringsAsFactors = FALSE)
  }
  out
}

.ee_dif_frase_pares <- function(ph, alfa, metodo_txt, simb = "d", col_ef = "d_qme", decimal = ".", idioma = "pt") {
  if (is.null(ph) || !nrow(ph)) return(character())
  sig <- ph[ph$p_ajustado < alfa, , drop = FALSE]
  um <- function(l) paste0(l$comparacao, if (!is.null(l$em_by)) paste0(" [", l$em_by, "]") else "", ": dif = ",
                           .ee_dif_m(l$dif, decimal), .ee_sep(decimal), .pr_p(l$p_ajustado, decimal),
                           if (!is.na(l[[col_ef]])) paste0(.ee_sep(decimal), simb, " = ", .pr_num(l[[col_ef]], 2, TRUE, decimal)) else "")
  if (idioma == "en") {
    if (!nrow(sig)) return(paste0("Pairwise comparisons (", metodo_txt, "): no pair differed after correction (smallest adjusted ",
                                  .pr_p(min(ph$p_ajustado), decimal), .ee_sep(decimal), ph$comparacao[which.min(ph$p_ajustado)], ")."))
    return(paste0("Pairwise comparisons (", metodo_txt, ", adjusted p): ", paste(vapply(seq_len(nrow(sig)), function(i) um(sig[i, ]), ""), collapse = "; "), "."))
  }
  if (!nrow(sig)) return(paste0("Nas comparações par a par (", metodo_txt, "), nenhum par diferiu após a correção (menor p ajustado: ",
                                .pr_p(min(ph$p_ajustado), decimal), .ee_sep(decimal), ph$comparacao[which.min(ph$p_ajustado)], ")."))
  paste0("Comparações par a par (", metodo_txt, ", p ajustado): ", paste(vapply(seq_len(nrow(sig)), function(i) um(sig[i, ]), ""), collapse = "; "), ".")
}

.ee_dif_es_tab <- function(tab, partial, nivel) {
  om <- suppressWarnings(suppressMessages(effectsize::omega_squared(tab, partial = partial, ci = nivel, alternative = "two.sided")))
  et <- suppressWarnings(suppressMessages(effectsize::eta_squared(tab, partial = partial, ci = nivel, alternative = "two.sided")))
  list(om = as.data.frame(om), et = as.data.frame(et))
}

.ee_dif_nomes <- function(termos, mapa) {
  vapply(termos, function(t) paste(vapply(strsplit(t, ":", fixed = TRUE)[[1]], function(p) if (p %in% names(mapa)) mapa[[p]] else p, ""),
                                   collapse = " × "), "", USE.NAMES = FALSE)
}

# ================================================================ ANOVA de um fator
# Rota clássica (livro): Levene sem evidência -> ANOVA padrão + Tukey; Levene
# aponta desigualdade -> ANOVA de Welch + Games-Howell. rota = "welch": Welch sempre
# (Delacre et al., 2019). Brown-Forsythe F* sai na tabela. ω² principal (Caixa 14.1).
# Pós-hoc só com ômnibus significativo, salvo exigir_omnibus = FALSE.
anova_um_fator <- function(dados, vd, grupo, rota = c("classica", "welch"),
                           posthoc = c("auto", "tukey", "bonferroni", "holm", "games_howell", "nenhum"),
                           exigir_omnibus = TRUE, centro_levene = c("media", "mediana"),
                           alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  rota <- match.arg(rota); posthoc <- match.arg(posthoc); centro_levene <- match.arg(centro_levene)
  .pr_exige("car", "emmeans", "effectsize")
  dd <- .ee_dados(dados, c(vd, grupo), numericas = vd)
  y <- dd$d[[vd]]; g <- .ee_fator(dd$d[[grupo]], grupo)
  k <- nlevels(g); N <- length(y); ns <- table(g)
  if (any(ns < 2)) stop("cada grupo precisa de pelo menos 2 casos")
  .ee_dif_exige_variacao(y, vd, g, grupo, "teste_kruskal_wallis(), que roda com separação perfeita")
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op), add = TRUE)
  d0 <- data.frame(y = y, g = g)
  m <- stats::aov(y ~ g, data = d0)
  t3 <- car::Anova(m, type = 3)
  we <- stats::oneway.test(y ~ g, data = d0, var.equal = FALSE)
  bf <- .ee_dif_bf(y, g)
  lev <- .ee_dif_levene(y, g, centro_levene)
  # Levene degenerado (G-07): não decide a rota; vale a de Welch
  lev_deg <- .ee_dif_levene_degenerado(y, g, lev, centro_levene)
  degenerado <- nzchar(lev_deg$motivo)
  usa_welch <- rota == "welch" || degenerado || isTRUE(lev[["p"]] < alfa)
  tab <- data.frame(teste = c("ANOVA clássica (SQ Tipo III)", "ANOVA de Welch", "Brown-Forsythe (F*)"),
                    sq = c(t3["g", "Sum Sq"], NA, NA), sq_residuo = c(t3["Residuals", "Sum Sq"], NA, NA),
                    gl1 = c(t3["g", "Df"], unname(we$parameter[1]), bf[["gl1"]]),
                    gl2 = c(t3["Residuals", "Df"], unname(we$parameter[2]), bf[["gl2"]]),
                    F = c(t3["g", "F value"], unname(we$statistic), bf[["F"]]),
                    p = c(t3["g", "Pr(>F)"], we$p.value, bf[["p"]]),
                    escolhido = c(!usa_welch, usa_welch, FALSE), stringsAsFactors = FALSE)
  es <- .ee_dif_es_tab(t3, FALSE, nivel)
  ef <- .ee_dif_efeito_df(c("ω²", "η² (igual ao η²p com um fator)"), c(es$om[1, 2], es$et[1, 2]),
                          c(es$om$CI_low[1], es$et$CI_low[1]), c(es$om$CI_high[1], es$et$CI_high[1]), c(TRUE, FALSE))
  # Rota Welch (G-01): o ω² principal continua o da tabela clássica (como o JASP), que usa o
  # QM do erro combinado; ao lado, o ω² aproximado pelo F de Welch, (F - 1)gl1/(F·gl1 + gl2 + 1)
  # com os gl de Welch: é o que effectsize::omega_squared(oneway.test(var.equal = FALSE))
  # devolve, marcado pelo próprio effectsize como aproximação.
  om_w <- NULL
  if (usa_welch) {
    om_w <- as.data.frame(suppressWarnings(suppressMessages(effectsize::F_to_omega2(unname(we$statistic), unname(we$parameter[1]), unname(we$parameter[2]),
                                                                                     ci = nivel, alternative = "two.sided"))))
    ef <- rbind(ef, .ee_dif_efeito_df("ω² aproximado pelo F de Welch (effectsize::F_to_omega2 com os gl de Welch)", om_w[1, 1], om_w$CI_low[1], om_w$CI_high[1], FALSE))
  }
  res <- stats::residuals(m); sw <- .ee_dif_sw(res)
  pres <- rbind(
    .ee_dif_pres("Normalidade dos resíduos (Shapiro-Wilk)", sw[["W"]], "", sw[["p"]], .ee_dif_le_normal(res, "resíduos")),
    .ee_dif_pres(paste0("Homogeneidade de variâncias (Levene, centro na ", if (centro_levene == "media") "média" else "mediana", ")"),
                 if (lev_deg$ruido) NA_real_ else lev[["F"]], .ee_dif_gltxt(lev[["gl1"]], lev[["gl2"]]), if (lev_deg$ruido) NA_real_ else lev[["p"]],
                 if (degenerado) paste0("degenerado (", lev_deg$motivo, "): não decide a rota") else .ee_dif_le_levene(lev[["p"]], alfa)))
  omni <- tab[if (usa_welch) 2 else 1, ]
  metodo <- if (posthoc == "auto") (if (usa_welch) "games_howell" else "tukey") else posthoc
  roda <- metodo != "nenhum" && k > 2 && (!exigir_omnibus || isTRUE(omni$p < alfa))
  ph <- NULL
  if (roda) {
    ph <- if (metodo == "games_howell") .ee_dif_games_howell(y, g, nivel) else
      .ee_dif_pares_emm(emmeans::emmeans(m, "g"), metodo, nivel)
    pr <- utils::combn(k, 2)
    ph$d_qme <- ph$dif / stats::sigma(m)
    ph$g_par <- vapply(seq_len(ncol(pr)), function(j) .ee_dif_hedges_par(y[g == levels(g)[pr[1, j]]], y[g == levels(g)[pr[2, j]]]), 1)
    ph$metodo <- metodo
  }
  nome_m <- c(tukey = "Tukey", bonferroni = "Bonferroni", holm = "Holm", games_howell = "Games-Howell", nenhum = "nenhum")[[metodo]]
  ltxt <- if (lev_deg$ruido) paste0("Levene sem valor interpretável (", lev_deg$motivo, ")") else
    paste0("Levene ", .ee_estat("F", lev[["F"]], c(lev[["gl1"]], lev[["gl2"]]), decimal), .ee_sep(decimal), .pr_p(lev[["p"]], decimal))
  dec <- c(
    if (rota == "classica" && degenerado) paste0("Rota clássica do livro, mas o Levene está degenerado (", lev_deg$motivo,
                                                 "): ele não decide a rota; usadas a ANOVA de Welch e o Games-Howell, que não supõem variâncias iguais (Delacre et al., 2019).")
    else if (rota == "classica") paste0("Rota clássica do livro: ", ltxt, if (usa_welch) " -> variâncias desiguais: ANOVA de Welch e Games-Howell." else
      " -> sem evidência de variâncias desiguais: ANOVA padrão e Tukey.", " Welch e Brown-Forsythe ficam na tabela (Delacre et al., 2019, recomendam Welch por padrão).")
    else "Rota Welch por padrão (Delacre et al., 2019): ANOVA de Welch e Games-Howell, sem consultar o Levene.",
    "SQ Tipo III com contrastes de soma zero (padrão do JASP); com um fator os tipos coincidem.",
    "ω² é o tamanho de efeito principal da ANOVA de um fator (Caixa 14.1); η² sai ao lado (superestima o efeito populacional). IC pelo F não central, bilateral.",
    if (usa_welch) paste0("Ômnibus de Welch: o ω² principal vem da ANOVA clássica (como no JASP), que usa o QM do erro combinado e supõe variâncias iguais; ",
                          "por isso pode divergir do F de Welch. Ao lado, o ω² aproximado pelo F de Welch, (F - 1)gl1/(F·gl1 + gl2 + 1) com os gl de Welch ",
                          "(effectsize::F_to_omega2; é o que effectsize::omega_squared() devolve para oneway.test(var.equal = FALSE), marcado como aproximação)."),
    if (k == 2) "Com 2 grupos, F = t² e o resultado é o do teste t independente; não há pós-hoc." else
      if (!roda && metodo != "nenhum") paste0("Ômnibus sem significância (", .pr_p(omni$p, decimal), "): pós-hoc não rodado; exigir_omnibus = FALSE roda mesmo assim.")
    else if (roda) paste0("Pós-hoc: ", nome_m, " (Caixa 14.3)",
                          if (posthoc == "auto" && metodo == "tukey" && ncol(utils::combn(k, 2)) <= 4)
                            paste0("; com ", ncol(utils::combn(k, 2)), " comparações o livro também aceita Bonferroni (posthoc = \"bonferroni\"); o automático fica no Tukey, que cobre todos os pares") else "",
                          ". d_qme = diferença / raiz do QM do erro (a coluna 'Cohen's d' do JASP, a da frase); ",
                          "g_par = g de Hedges só com os dois grupos do par (o índice que o corpo do livro cita; a correção de Hedges é de amostra pequena, não de desbalanceamento)."))
  avisos <- character()
  if (rota == "classica" && degenerado)
    avisos <- c(avisos, paste0("Levene degenerado (", lev_deg$motivo, "): o teste de homogeneidade não informa nada com esses dados; a rota não foi decidida por ele (ANOVA de Welch e Games-Howell)."))
  if (roda && metodo %in% c("tukey", "bonferroni", "holm") && (degenerado || isTRUE(lev[["p"]] < alfa)))
    avisos <- c(avisos, if (degenerado) "O pós-hoc escolhido supõe variâncias iguais e o Levene está degenerado: o livro indica Games-Howell quando a igualdade não se sustenta." else
      "Levene indica variâncias desiguais e o pós-hoc escolhido as supõe iguais: o livro indica Games-Howell nesse caso.")
  sig <- isTRUE(omni$p < alfa)
  if (usa_welch) {
    avisos <- c(avisos, paste0("Ômnibus de Welch com ω² da ANOVA clássica (como no JASP): os dois podem divergir, porque o ω² clássico usa o QM do erro combinado, ",
                               "que supõe variâncias iguais; o ω² aproximado pelo F de Welch (", .pr_num(om_w[1, 1], 3, FALSE, decimal), ") sai ao lado na frase e em $efeito."))
    if (sig && (isTRUE(ef$valor[1] <= 0) || isTRUE(ef$ic_inf[1] <= 0)))
      avisos <- c(avisos, paste0("Contradição aparente: a ANOVA de Welch indica efeito (", .pr_p(omni$p, decimal), "), mas o ω² da ANOVA clássica é ",
                                 .pr_num(ef$valor[1], 3, FALSE, decimal), " com IC que inclui zero (a ANOVA clássica dá ",
                                 .ee_junta(.ee_estat("F", tab$F[1], c(tab$gl1[1], tab$gl2[1]), decimal), .pr_p(tab$p[1], decimal), decimal = decimal),
                                 "). Não ler o ω² clássico como ausência de efeito: relatar os dois índices e a divergência, ou usar o ω² aproximado pelo F de Welch."))
  }
  st_txt <- .ee_junta(.ee_estat("F", omni$F, c(omni$gl1, omni$gl2), decimal), .pr_p(omni$p, decimal),
                      .ee_efeito_txt("ω²", ef$valor[1], ef$ic_inf[1], ef$ic_sup[1], TRUE, 3, decimal, idioma, nivel = nivel), decimal = decimal)
  rot <- .ee_dif_rotulo(ef$valor[1], "var", idioma)
  # rota Welch: declara a origem do ω² e traz o aproximado pelo F de Welch
  w_txt <- if (!usa_welch) "" else paste0("; ", if (idioma == "en") "from Welch's F, approximate " else "pelo F de Welch, ",
                                          .ee_efeito_txt(if (idioma == "en") "ω²" else "ω² aproximado", om_w[1, 1], om_w$CI_low[1], om_w$CI_high[1], TRUE, 3, decimal, idioma, nivel = nivel))
  simb <- "d"; col <- "d_qme"
  if (idioma == "en") {
    apa <- paste0(if (usa_welch) "A Welch one-way ANOVA " else "A one-way ANOVA ", if (sig) "indicated an effect of " else "did not indicate a statistically significant effect of ",
                  grupo, " on ", vd, ", ", st_txt, " (", if (usa_welch) "ω² from the classical ANOVA, as in JASP: " else "", rot, " effect", w_txt, ").")
  } else {
    apa <- paste0(if (usa_welch) "Uma ANOVA de Welch de um fator " else "Uma ANOVA de um fator ", if (sig) "indicou efeito de " else "não indicou efeito estatisticamente significativo de ",
                  grupo, " sobre ", vd, ", ", st_txt, " (", if (usa_welch) "ω² da ANOVA clássica, como no JASP: " else "", "efeito ", rot, w_txt, ").")
  }
  if (roda) apa <- c(apa, .ee_dif_frase_pares(ph, alfa, nome_m, simb, col, decimal, idioma))
  avisos <- c(avisos, .ee_aviso_nao_sig(omni$p, alfa, NULL, "diferença entre os grupos", NULL, decimal))
  .ee_resultado("ANOVA de um fator entre grupos", n = N, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef,
                pressupostos = pres, posthoc = ph, decisoes = dec, avisos = avisos, apa = apa,
                extra = list(descritivas = .ee_dif_desc(y, g), modelo = m))
}

# ================================================================ medidas repetidas e mista
# Delineamentos com fator(es) intrassujeito, via afex::aov_car (SQ Tipo III).
#  medidas: colunas do banco em formato largo, uma por célula intrassujeito.
#  intra: NULL (um fator "condicao" cujos níveis são os nomes das colunas), um nome
#   (um fator com esse nome) ou lista nomeada de níveis para 2+ fatores, ex.:
#   list(bloco = c("A", "B"), condicao = c("C", "I", "N")); as colunas em 'medidas'
#   seguem a ordem das células com o ÚLTIMO fator variando mais rápido (A-C, A-I, A-N,
#   B-C, ...), a mesma ordem das células do JASP.
#  entre: fator(es) entre sujeitos (desenho misto).
#  esfericidade = "classica" (livro): Mauchly decide se corrige; a correção é
#   GG quando o ε de GG arredondado a duas casas é <= limiar_epsilon (o livro adota GG
#   em ε = ,752, "sobre o limiar" e) e HF acima. "epsilon": decide só
#   pelo ε, sem o Mauchly como porta (Blanca et al., 2023).
#  η²p é o principal em todo desenho com fator intra (: o pacote adota η²p fora
#  da família entre grupos); η²G e ω² saem ao lado.
anova_medidas_repetidas <- function(dados, medidas, intra = NULL, entre = NULL,
                                    esfericidade = c("classica", "epsilon"), limiar_epsilon = .75,
                                    posthoc = c("bonferroni", "holm", "nenhum"), exigir_omnibus = TRUE,
                                    alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  esfericidade <- match.arg(esfericidade); posthoc <- match.arg(posthoc)
  .pr_exige("afex", "emmeans", "effectsize", "car")
  if (is.null(intra)) intra <- list(condicao = medidas)
  if (is.character(intra) && length(intra) == 1 && is.null(names(intra))) intra <- stats::setNames(list(medidas), intra)
  if (!is.list(intra) || is.null(names(intra)) || any(!nzchar(names(intra))))
    stop("intra deve ser NULL, um nome ou uma lista nomeada de níveis, ex.: list(condicao = c('A', 'B', 'C'))")
  if (length(medidas) != prod(lengths(intra)))
    stop("o número de colunas em 'medidas' (", length(medidas), ") precisa ser o produto dos níveis de 'intra' (", prod(lengths(intra)), ")")
  dd <- .ee_dados(dados, c(medidas, entre), numericas = medidas)
  w <- dd$d; n <- nrow(w)
  iw <- paste0("w", seq_along(intra)); ib <- if (length(entre)) paste0("b", seq_along(entre)) else character()
  mapa <- stats::setNames(c(names(intra), entre), c(iw, ib))
  grade <- expand.grid(rev(intra), stringsAsFactors = FALSE)[, rev(seq_along(intra)), drop = FALSE]
  names(grade) <- iw
  bet <- lapply(entre, function(b) .ee_fator(w[[b]], b)); names(bet) <- ib
  lg <- do.call(rbind, lapply(seq_along(medidas), function(j) {
    x <- data.frame(.id = factor(seq_len(n)), .y = w[[medidas[j]]])
    for (b in ib) x[[b]] <- bet[[b]]
    for (v in iw) x[[v]] <- factor(rep(grade[j, v], n), levels = intra[[match(v, iw)]])
    x
  }))
  f <- stats::as.formula(paste(".y ~", if (length(ib)) paste(ib, collapse = " * ") else "1",
                               "+ Error(.id/(", paste(iw, collapse = " * "), "))"))
  a <- suppressMessages(afex::aov_car(f, data = lg, type = 3, factorize = FALSE,
                                      anova_table = list(correction = "none", es = c("pes", "ges"))))
  s <- suppressWarnings(summary(a$Anova, multivariate = FALSE))  # "HF eps > 1 treated as 1": truncado como o JASP
  ut <- unclass(s$univariate.tests)
  b_err <- ut["(Intercept)", c("Error SS", "den Df")]
  ut <- ut[rownames(ut) != "(Intercept)", , drop = FALSE]
  termos <- rownames(ut)
  esf <- if (!is.null(s$sphericity.tests) && length(s$sphericity.tests)) unclass(s$sphericity.tests) else NULL
  adj <- if (!is.null(s$pval.adjustments) && length(s$pval.adjustments)) unclass(s$pval.adjustments) else NULL
  gl_err_mult <- a$lm$df.residual
  linhas <- lapply(termos, function(tm) {
    F <- ut[tm, "F value"]; g1 <- ut[tm, "num Df"]; g2 <- ut[tm, "den Df"]
    W <- pW <- egg <- ehf <- chi <- gchi <- NA_real_
    if (!is.null(esf) && tm %in% rownames(esf) && is.finite(esf[tm, 1])) {
      W <- esf[tm, 1]; pW <- esf[tm, 2]
      egg <- min(1, adj[tm, "GG eps"]); ehf <- min(1, adj[tm, "HF eps"])
      pp <- g1
      chi <- -(gl_err_mult - (2 * pp^2 + pp + 2) / (6 * pp)) * log(W); gchi <- pp * (pp + 1) / 2 - 1
    }
    corr <- "nenhuma"
    if (is.finite(egg)) {
      precisa <- if (esfericidade == "classica") isTRUE(pW < alfa) else egg < 1
      if (precisa) corr <- if (round(egg, 2) <= limiar_epsilon) "Greenhouse-Geisser" else "Huynh-Feldt"
    }
    eps <- switch(corr, nenhuma = 1, "Greenhouse-Geisser" = egg, "Huynh-Feldt" = ehf)
    data.frame(termo = tm, sq = ut[tm, "Sum Sq"], gl1 = g1, sq_erro = ut[tm, "Error SS"], gl2 = g2, F = F,
               p = ut[tm, "Pr(>F)"], correcao = corr, epsilon = eps, gl1_corr = g1 * eps, gl2_corr = g2 * eps,
               p_corr = stats::pf(F, g1 * eps, g2 * eps, lower.tail = FALSE),
               mauchly_W = W, mauchly_chi2 = chi, mauchly_gl = gchi, mauchly_p = pW, eps_gg = egg, eps_hf = ehf,
               eps_lb = if (is.finite(egg)) 1 / g1 else NA_real_, stringsAsFactors = FALSE)
  })
  tb <- do.call(rbind, linhas)
  tb$intra <- vapply(strsplit(tb$termo, ":"), function(p) any(p %in% iw), logical(1))
  nomes <- .ee_dif_nomes(tb$termo, mapa)
  # tamanhos de efeito
  pes <- tb$sq / (tb$sq + tb$sq_erro)
  ges <- a$anova_table[tb$termo, "ges"]
  so_um <- length(iw) == 1 && !length(ib)
  if (so_um) {
    MSm <- tb$sq / tb$gl1; MSr <- tb$sq_erro / tb$gl2; MSb <- b_err[[1]] / b_err[[2]]; df <- tb$gl1
    num <- df / (n * (df + 1)) * (MSm - MSr)
    om <- max(0, num / (MSr + (MSb - MSr) / (df + 1) + num))
    om_lab <- "ω² (generalizado, fórmula do JASP para medidas repetidas)"; om_ic <- c(NA, NA)
    om <- rep(om, nrow(tb)); om_lo <- om_hi <- rep(NA_real_, nrow(tb))
  } else {
    oe <- as.data.frame(suppressWarnings(suppressMessages(effectsize::omega_squared(a, partial = TRUE, ci = nivel, alternative = "two.sided"))))
    ix <- match(tb$termo, oe$Parameter)
    om <- oe[[2]][ix]; om_lo <- oe$CI_low[ix]; om_hi <- oe$CI_high[ix]
    om_lab <- "ω²p (effectsize)"
  }
  ef <- do.call(rbind, lapply(seq_len(nrow(tb)), function(i) .ee_dif_efeito_df(c("η²p", "η²G (generalizado)", om_lab),
                                                                             c(pes[i], ges[i], om[i]), c(NA, NA, om_lo[i]), c(NA, NA, om_hi[i]),
                                                                             c(TRUE, FALSE, FALSE), termo = nomes[i])))
  tab <- tb[, c("termo", "sq", "gl1", "sq_erro", "gl2", "F", "p", "correcao", "epsilon", "gl1_corr", "gl2_corr", "p_corr")]
  tab$termo <- nomes
  # pressupostos
  pres <- list()
  for (i in which(is.finite(tb$mauchly_W))) {
    pres[[length(pres) + 1]] <- .ee_dif_pres(paste0("Esfericidade de ", nomes[i], " (Mauchly W)"), tb$mauchly_W[i], .ee_gl(tb$mauchly_gl[i]), tb$mauchly_p[i],
      paste0("χ² = ", .pr_num(tb$mauchly_chi2[i], 3, TRUE, decimal), "; ε GG = ", .pr_num(tb$eps_gg[i], 3, TRUE, decimal), .ee_sep(decimal), "HF = ",
             .pr_num(tb$eps_hf[i], 3, TRUE, decimal), .ee_sep(decimal), "LB = ", .pr_num(tb$eps_lb[i], 3, TRUE, decimal),
             if (tb$mauchly_p[i] < alfa) "; evidência de violação" else "; sem evidência de violação (não prova esfericidade)",
             "; correção usada: ", tb$correcao[i]))
  }
  form_r <- stats::as.formula(paste(".y ~ .id +", paste(c(paste(iw, collapse = "*"), if (length(ib)) paste0("(", paste(iw, collapse = "*"), "):(", paste(ib, collapse = "*"), ")")), collapse = " + ")))
  res <- stats::residuals(stats::lm(form_r, data = lg)); sw <- .ee_dif_sw(res)
  pres[[length(pres) + 1]] <- .ee_dif_pres("Normalidade dos resíduos (Shapiro-Wilk)", sw[["W"]], "", sw[["p"]], .ee_dif_le_normal(res, "resíduos"))
  box <- NULL
  if (length(ib)) {
    cel_b <- interaction(bet, drop = TRUE)
    for (j in seq_along(medidas)) {
      lv <- .ee_dif_levene(w[[medidas[j]]], cel_b, "media")
      pres[[length(pres) + 1]] <- .ee_dif_pres(paste0("Homogeneidade de variâncias em ", medidas[j], " (Levene, centro na média)"), lv[["F"]],
                                               .ee_dif_gltxt(lv[["gl1"]], lv[["gl2"]]), lv[["p"]], .ee_dif_le_levene(lv[["p"]], alfa))
    }
    box <- tryCatch(suppressWarnings(heplots::boxM(as.matrix(w[, medidas]), cel_b)), error = function(e) NULL)
    if (!is.null(box)) pres[[length(pres) + 1]] <- .ee_dif_pres("Homogeneidade das matrizes de covariância (M de Box)", unname(box$statistic),
                                                                .ee_gl(unname(box$parameter)), box$p.value,
                                                                .ee_dif_le_box(box$p.value, .001))
  }
  pres <- do.call(rbind, pres)
  # comparações
  sigma_jasp <- sqrt(mean(stats::sigma(a$lm)^2))
  # o afex torna sintáticos os níveis intra ("1" vira "X1"); volta aos rótulos originais
  traduz <- function(v, x) {
    if (!v %in% iw) return(x)
    L <- intra[[match(v, iw)]]; m <- match(x, make.names(L)); m2 <- match(x, L)
    ifelse(!is.na(m2), x, ifelse(!is.na(m), L[m], x))
  }
  p_usado <- tb$p_corr
  sig_t <- p_usado < alfa
  ordem <- lengths(strsplit(tb$termo, ":"))
  ph <- NULL; dec_ph <- character()
  inter2 <- which(ordem == 2 & (sig_t | !exigir_omnibus))
  if (posthoc != "nenhum") {
    if (length(inter2) && any(sig_t[inter2])) {
      inter2 <- inter2[sig_t[inter2]]
      blocos <- list()
      for (i in inter2) {
        fs <- strsplit(tb$termo[i], ":")[[1]]
        for (o in list(fs, rev(fs))) {
          em <- emmeans::emmeans(a, o[1], by = o[2], model = "multivariate")
          pp <- .ee_dif_pares_emm(em, "none", nivel, traduz)
          pp$efeito_simples <- paste0(mapa[[o[1]]], " em cada nível de ", mapa[[o[2]]])
          blocos[[length(blocos) + 1]] <- pp
        }
      }
      ph <- do.call(rbind, blocos)
      ph$em_by <- .ee_dif_troca_nomes(ph$em_by, mapa)
      ph$p_sem_ajuste <- ph$p_ajustado
      ph$p_ajustado <- stats::p.adjust(ph$p_sem_ajuste, if (posthoc == "holm") "holm" else "bonferroni")
      ph$d_jasp <- ph$dif / sigma_jasp
      dec_ph <- paste0("Interação significativa: efeitos simples (cada fator dentro dos níveis do outro), com ", if (posthoc == "holm") "Holm" else "Bonferroni",
                       " sobre o conjunto de ", nrow(ph), " contrastes (Caixa 16.A). ICs não ajustados. d_jasp = diferença / DP médio das células (o d do JASP).")
    } else {
      principais <- which(ordem == 1 & (sig_t | !exigir_omnibus))
      blocos <- list()
      for (i in principais) {
        fx <- tb$termo[i]
        nlev <- if (fx %in% iw) length(intra[[match(fx, iw)]]) else nlevels(bet[[fx]])
        if (nlev < 3) next
        em <- emmeans::emmeans(a, fx, model = "multivariate")
        pp <- .ee_dif_pares_emm(em, posthoc, nivel, traduz)
        pp$efeito <- mapa[[fx]]
        pp$d_jasp <- pp$dif / sigma_jasp
        if (so_um) {
          pr <- utils::combn(length(medidas), 2)
          pp$d_par <- vapply(seq_len(ncol(pr)), function(j) {
            x1 <- w[[medidas[pr[1, j]]]]; x2 <- w[[medidas[pr[2, j]]]]
            mean(x1 - x2) / sqrt((stats::var(x1) + stats::var(x2)) / 2) }, 1)
          pp$d_z <- vapply(seq_len(ncol(pr)), function(j) {
            dx <- w[[medidas[pr[1, j]]]] - w[[medidas[pr[2, j]]]]; mean(dx) / stats::sd(dx) }, 1)
        }
        blocos[[length(blocos) + 1]] <- pp
      }
      if (length(blocos)) {
        ph <- do.call(rbind, blocos)
        dec_ph <- paste0("Pós-hoc ", if (posthoc == "holm") "Holm" else "Bonferroni", " (Caixa 14.3) sobre o modelo multivariado do afex: erro separado por comparação, gl = n - 1 como o t(18) do livro.",
                         if (so_um) " d_par = d na métrica original do par (DP combinado das duas condições), a versão que o livro interpreta (,35 no Stroop); d_jasp = diferença / DP médio de todas as condições (a coluna do JASP, -0,366); d_z = d das diferenças (inflado)." else " d_jasp = diferença / DP médio das células (o d do JASP).")
      }
    }
  }
  # decisões
  dec <- c(
    "SQ Tipo III (afex::aov_car, contrastes de soma zero), padrão do JASP.",
    if (esfericidade == "classica") paste0("Esfericidade pela rota do livro: Mauchly decide se corrige; GG se o ε de GG arredondado a duas casas for <= ", limiar_epsilon,
                                          " (fronteira inclusiva: o livro adota GG em ε = ,752), HF acima (Caixa 14.2). Blanca et al. (2023) sugerem decidir só pelo ε (esfericidade = 'epsilon') e corte perto de ,60.")
    else paste0("Esfericidade decidida só pelo ε, sem o Mauchly como porta (Blanca et al., 2023): GG se ε de GG <= ", limiar_epsilon, ", HF acima."),
    "ε de HF truncado em 1, como o JASP; ε < 1 multiplica os dois gl.",
    "η²p é o tamanho de efeito principal em desenhos com fator intrassujeito; η²G (Lakens, 2013) e ω² devolvem a variância entre participantes ao denominador e ficam menores.",
    if (length(ib)) "Desenho misto: o curso prefere ω²p como principal; o livro reserva o ω²p para desenhos só entre grupos e usa η²p fora dessa família, e é o que vale aqui." else NULL,
    dec_ph,
    if (posthoc != "nenhum" && is.null(ph)) "Nenhum termo pediu comparações par a par (sem efeito significativo com 3 ou mais níveis, ou fatores de 2 níveis, cujo F já é a comparação)." else NULL,
    "Casos com qualquer medida ausente saem da análise (listwise); modelos mistos ou GEE usam os casos incompletos.")
  avisos <- character()
  if (any(ordem >= 3 & sig_t)) avisos <- c(avisos, "Há interação de 3ª ordem significativa: os efeitos simples automáticos cobrem só interações de 2 fatores; fatiar a interação por um dos fatores.")
  # frase APA
  ordem_apa <- order(-ordem)
  frases <- vapply(ordem_apa, function(i) {
    st <- .ee_junta(.ee_estat("F", tb$F[i], c(tb$gl1_corr[i], tb$gl2_corr[i]), decimal), .pr_p(tb$p_corr[i], decimal),
                    .ee_efeito_txt("η²p", pes[i], NA, NA, TRUE, 3, decimal, idioma, nivel = nivel), decimal = decimal)
    rot <- .ee_dif_rotulo(pes[i], "var", idioma)
    cor_txt <- if (tb$correcao[i] != "nenhuma") {
      # itens do Mauchly separados pelo decimal (G-11: com decimal = ",", "W = ,67; χ²(2) = 6,80; p = ,033")
      mau <- .ee_junta(paste0("W = ", .pr_num(tb$mauchly_W[i], 2, FALSE, decimal)), paste0("χ²(", .ee_gl(tb$mauchly_gl[i]), ") = ", .pr_num(tb$mauchly_chi2[i], 2, TRUE, decimal)),
                       .pr_p(tb$mauchly_p[i], decimal), decimal = decimal)
      if (idioma == "en") paste0("; ", tb$correcao[i], " correction, ε = ", .pr_num(tb$epsilon[i], 2, FALSE, decimal), ", after Mauchly's test, ", mau)
      else paste0("; correção de ", tb$correcao[i], ", ε = ", .pr_num(tb$epsilon[i], 2, FALSE, decimal), ", após o teste de Mauchly, ", mau)
    } else ""
    if (idioma == "en") paste0(if (ordem[i] > 1) "Interaction " else "Effect of ", nomes[i], ": ", st, " (", rot, " effect", cor_txt, ").")
    else paste0(if (ordem[i] > 1) "Interação " else "Efeito de ", nomes[i], ": ", st, " (efeito ", rot, cor_txt, ").")
  }, "")
  tipo <- if (length(ib)) "ANOVA mista" else if (length(iw) > 1) "ANOVA fatorial de medidas repetidas" else "ANOVA de medidas repetidas"
  abre <- if (idioma == "en") paste0(c("ANOVA mista" = "Mixed ANOVA", "ANOVA fatorial de medidas repetidas" = "Factorial repeated-measures ANOVA",
                                       "ANOVA de medidas repetidas" = "Repeated-measures ANOVA")[[tipo]], " (n = ", n, "). ") else paste0(tipo, " (n = ", n, "). ")
  apa <- paste0(abre, paste(frases, collapse = " "))
  if (!is.null(ph)) {
    col <- if (so_um && "d_par" %in% names(ph)) "d_par" else "d_jasp"
    apa <- c(apa, .ee_dif_frase_pares(ph, alfa, if (posthoc == "holm") "Holm" else "Bonferroni", "d", col, decimal, idioma))
  }
  desc <- do.call(rbind, lapply(seq_along(medidas), function(j) {
    if (length(ib)) {
      cb <- interaction(bet, drop = TRUE, sep = " / ")
      do.call(rbind, lapply(levels(cb), function(l) { x <- w[[medidas[j]]][cb == l]
        data.frame(medida = medidas[j], grupo = l, n = length(x), media = mean(x), dp = stats::sd(x), stringsAsFactors = FALSE) }))
    } else data.frame(medida = medidas[j], n = n, media = mean(w[[medidas[j]]]), dp = stats::sd(w[[medidas[j]]]), stringsAsFactors = FALSE)
  }))
  avisos <- c(avisos, .ee_aviso_nao_sig(p_usado, alfa, NULL, "efeito", NULL, decimal,
                                        rotulos = paste0(nomes, " (", vapply(p_usado, .pr_p, "", decimal = decimal), ")")))
  .ee_resultado(tipo, n = n, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef, pressupostos = pres,
                posthoc = ph, decisoes = dec, avisos = avisos, apa = apa,
                extra = list(descritivas = desc, celulas = cbind(medida = medidas, grade), modelo = a, box_m = box))
}

.ee_dif_troca_nomes <- function(txt, mapa) {
  for (k in names(mapa)) txt <- gsub(paste0("\\b", k, " = "), paste0(mapa[[k]], " = "), txt)
  txt
}

# ================================================================ ANOVA fatorial entre grupos
# 2 ou mais fatores entre sujeitos, SQ Tipo III (Caixa 16.B; tipo_sq = 2 como
# alternativa, Langsrud, 2003). ω²p principal e η²p ao lado. Interação lida
# primeiro; se significativa, efeitos simples (cada fator dentro dos níveis do outro)
# com Bonferroni sobre o conjunto e d sobre a raiz do QM do erro (Caixa 16.A).
# Efeitos principais com 3+ níveis fora de interação significativa: Tukey (padrão).
anova_fatorial <- function(dados, vd, fatores, tipo_sq = 3, posthoc = c("tukey", "bonferroni", "holm", "nenhum"),
                           centro_levene = c("media", "mediana"), exigir_omnibus = TRUE,
                           alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  posthoc <- match.arg(posthoc); centro_levene <- match.arg(centro_levene)
  if (!tipo_sq %in% c(2, 3)) stop("tipo_sq deve ser 2 ou 3")
  if (length(fatores) < 2) stop("a ANOVA fatorial pede 2 ou mais fatores; com um fator use anova_um_fator()")
  .pr_exige("car", "emmeans", "effectsize")
  dd <- .ee_dados(dados, c(vd, fatores), numericas = vd)
  ifs <- paste0("f", seq_along(fatores)); mapa <- stats::setNames(fatores, ifs)
  d0 <- data.frame(.y = dd$d[[vd]])
  for (i in seq_along(fatores)) d0[[ifs[i]]] <- .ee_fator(dd$d[[fatores[i]]], fatores[i])
  cel <- interaction(d0[ifs], drop = FALSE)
  if (any(table(cel) == 0)) stop("há células vazias no cruzamento dos fatores; a SQ Tipo III não é estimável")
  .ee_dif_exige_variacao(d0$.y, vd, droplevels(cel), paste(fatores, collapse = " × "))
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op), add = TRUE)
  m <- stats::aov(stats::as.formula(paste(".y ~", paste(ifs, collapse = " * "))), data = d0)
  t3 <- car::Anova(m, type = tipo_sq)
  termos <- setdiff(rownames(t3), c("(Intercept)", "Residuals"))
  es <- .ee_dif_es_tab(t3, TRUE, nivel)
  nomes <- .ee_dif_nomes(termos, mapa)
  tab <- data.frame(termo = nomes, sq = t3[termos, "Sum Sq"], gl1 = t3[termos, "Df"], gl2 = t3["Residuals", "Df"],
                    F = t3[termos, "F value"], p = t3[termos, "Pr(>F)"], stringsAsFactors = FALSE)
  io <- match(termos, es$om$Parameter); ie <- match(termos, es$et$Parameter)
  ef <- do.call(rbind, lapply(seq_along(termos), function(i) .ee_dif_efeito_df(c("ω²p", "η²p"),
                                                                              c(es$om[io[i], 2], es$et[ie[i], 2]), c(es$om$CI_low[io[i]], es$et$CI_low[ie[i]]),
                                                                              c(es$om$CI_high[io[i]], es$et$CI_high[ie[i]]), c(TRUE, FALSE), termo = nomes[i])))
  res <- stats::residuals(m); sw <- .ee_dif_sw(res)
  lev <- if (centro_levene == "media") {
    a <- stats::anova(stats::lm(abs(res) ~ cel)); c(F = a[1, "F value"], gl1 = a[1, "Df"], gl2 = a[2, "Df"], p = a[1, "Pr(>F)"])
  } else .ee_dif_levene(d0$.y, droplevels(cel), "mediana")
  r2 <- summary(stats::lm(stats::formula(m), data = d0))$r.squared
  pres <- rbind(
    .ee_dif_pres("Normalidade dos resíduos (Shapiro-Wilk)", sw[["W"]], "", sw[["p"]], .ee_dif_le_normal(res, "resíduos")),
    .ee_dif_pres(paste0("Homogeneidade de variâncias entre as células (Levene, centro na ", if (centro_levene == "media") "média" else "mediana", ")"),
                 lev[["F"]], .ee_dif_gltxt(lev[["gl1"]], lev[["gl2"]]), lev[["p"]], .ee_dif_le_levene(lev[["p"]], alfa)))
  ordem <- lengths(strsplit(termos, ":"))
  sig_t <- tab$p < alfa
  sigma_m <- stats::sigma(m)
  ph <- NULL; simples_F <- NULL; dec_ph <- character()
  inter2 <- which(ordem == 2 & sig_t)
  # efeitos simples: Bonferroni sobre o conjunto (Caixa 16.A); posthoc = "holm" troca por
  # Holm; "tukey" (padrão, pensado para os efeitos principais) vira Bonferroni aqui
  aj_simples <- if (posthoc == "holm") "holm" else "bonferroni"
  nome_aj <- if (posthoc == "holm") "Holm" else "Bonferroni"
  if (length(inter2) && posthoc == "nenhum") {
    dec_ph <- "Interação significativa, mas posthoc = 'nenhum': efeitos simples não calculados (a Caixa 16.A, os recomenda quando a interação é substantiva)."
  } else if (length(inter2)) {
    blocos <- list(); jt <- list()
    for (i in inter2) {
      fs <- strsplit(termos[i], ":")[[1]]
      for (o in list(fs, rev(fs))) {
        em <- emmeans::emmeans(m, o[1], by = o[2])
        pp <- .ee_dif_pares_emm(em, "none", nivel)
        pp$efeito_simples <- paste0(mapa[[o[1]]], " em cada nível de ", mapa[[o[2]]])
        blocos[[length(blocos) + 1]] <- pp
        j <- as.data.frame(emmeans::joint_tests(m, by = o[2]))
        j <- j[j$`model term` == o[1], , drop = FALSE]
        jt[[length(jt) + 1]] <- data.frame(efeito_simples = paste0(mapa[[o[1]]], " em ", mapa[[o[2]]], " = ", j[[o[2]]]),
                                           gl1 = j$df1, gl2 = j$df2, F = j$F.ratio, p = j$p.value, stringsAsFactors = FALSE)
      }
    }
    ph <- do.call(rbind, blocos)
    ph$em_by <- .ee_dif_troca_nomes(ph$em_by, mapa)
    ph$p_sem_ajuste <- ph$p_ajustado
    ph$p_ajustado <- stats::p.adjust(ph$p_sem_ajuste, aj_simples)
    ph$d_qme <- ph$dif / sigma_m
    simples_F <- do.call(rbind, jt)
    dec_ph <- paste0("Interação significativa: efeitos simples com ", nome_aj, " sobre o conjunto de ", nrow(ph),
                     " contrastes (Caixa 16.A); t com o gl do erro do modelo; IC não ajustado; d = diferença / raiz do QM do erro. p_sem_ajuste fica na tabela (o livro imprime esses).")
  } else if (!length(inter2) && posthoc != "nenhum") {
    blocos <- list()
    for (i in which(ordem == 1 & (sig_t | !exigir_omnibus))) {
      fx <- termos[i]
      if (nlevels(d0[[fx]]) < 3) next
      pp <- .ee_dif_pares_emm(emmeans::emmeans(m, fx), posthoc, nivel)
      pp$efeito <- mapa[[fx]]; pp$d_qme <- pp$dif / sigma_m
      blocos[[length(blocos) + 1]] <- pp
    }
    if (length(blocos)) {
      ph <- do.call(rbind, blocos)
      dec_ph <- paste0("Sem interação significativa: pós-hoc ", c(tukey = "Tukey", bonferroni = "Bonferroni", holm = "Holm")[[posthoc]],
                       " nos efeitos principais com 3 ou mais níveis (médias marginais); d = diferença / raiz do QM do erro.")
    }
  }
  desbal <- length(unique(as.integer(table(cel)))) > 1
  dec <- c(paste0("SQ Tipo ", if (tipo_sq == 3) "III com contrastes de soma zero (o pré-cálculo do livro e o padrão do JASP)" else "II (Langsrud, 2003; adequado sem interação relevante)",
                  if (desbal) "; desenho desbalanceado: o tipo de SQ muda F e p." else "; desenho balanceado: os tipos coincidem."),
           "ω²p é o principal e η²p sai ao lado; as versões não parciais não são as reportadas.",
           "Leia a interação antes dos efeitos principais; com interação significativa, os efeitos principais se interpretam à luz dela.",
           paste0("R² do modelo completo = ", .pr_num(r2, 3, FALSE, decimal), "."), dec_ph)
  avisos <- character()
  if (isTRUE(lev[["p"]] < alfa)) avisos <- c(avisos, if (desbal) "Levene indica variâncias desiguais com células desbalanceadas: a ANOVA é menos robusta nesse caso; considerar erros-padrão robustos ou reamostragem." else
    "Levene indica variâncias desiguais; com células de tamanho igual a ANOVA é robusta a isso: registrar a violação.")
  if (any(ordem >= 3 & sig_t)) avisos <- c(avisos, "Interação de 3ª ordem significativa: fatiar por um dos fatores (os efeitos simples automáticos cobrem interações de 2 fatores).")
  frases <- vapply(order(-ordem), function(i) {
    st <- .ee_junta(.ee_estat("F", tab$F[i], c(tab$gl1[i], tab$gl2[i]), decimal), .pr_p(tab$p[i], decimal),
                    .ee_efeito_txt("ω²p", ef$valor[2 * i - 1], ef$ic_inf[2 * i - 1], ef$ic_sup[2 * i - 1], TRUE, 3, decimal, idioma, nivel = nivel), decimal = decimal)
    rot <- .ee_dif_rotulo(ef$valor[2 * i - 1], "var", idioma)
    if (idioma == "en") paste0(if (ordem[i] > 1) "Interaction " else "Main effect of ", nomes[i], ": ", st, " (", rot, " effect).")
    else paste0(if (ordem[i] > 1) "Interação " else "Efeito principal de ", nomes[i], ": ", st, " (efeito ", rot, ").")
  }, "")
  apa <- paste(frases, collapse = " ")
  if (!is.null(ph)) apa <- c(apa, .ee_dif_frase_pares(ph, alfa, if (length(inter2)) paste0("efeitos simples, ", nome_aj) else c(tukey = "Tukey", bonferroni = "Bonferroni", holm = "Holm")[[posthoc]],
                                                    "d", "d_qme", decimal, idioma))
  medias <- stats::aggregate(d0$.y, d0[ifs], function(x) c(n = length(x), media = mean(x), dp = stats::sd(x)))
  medias <- data.frame(medias[ifs], medias$x); names(medias)[seq_along(ifs)] <- fatores
  avisos <- c(avisos, .ee_aviso_nao_sig(tab$p, alfa, NULL, "efeito", NULL, decimal,
                                        rotulos = paste0(tab$termo, " (", vapply(tab$p, .pr_p, "", decimal = decimal), ")")))
  .ee_resultado(paste0("ANOVA fatorial entre grupos (", paste(fatores, collapse = " × "), ")"), n = nrow(d0),
                n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef, pressupostos = pres, posthoc = ph,
                decisoes = dec, avisos = avisos, apa = apa,
                extra = list(descritivas = medias, efeitos_simples_F = simples_F, r2 = r2, modelo = m))
}

# ================================================================ ANCOVA
# Um fator entre grupos e 1+ covariáveis contínuas, SQ Tipo III. ω²p principal e η²p ao
# lado. Checagens: homogeneidade das inclinações (fator × covariável; Caixa
# 15.A), Levene sobre a VD bruta e sobre os resíduos do modelo (o do JASP;
#), normalidade dos resíduos e diferença dos grupos na covariável (Miller &
# Chapman, 2001). Médias ajustadas (emmeans, covariáveis na média) e comparações
# Bonferroni sobre elas, com IC não ajustado e ajustado e d sobre a raiz do QM
# do erro. bootstrap = TRUE: IC das diferenças ajustadas por reamostragem de
# casos estratificada por grupo, n_boot = 1000 (Caixa 15.C); não conserta
# inclinações heterogêneas.
.ee_dif_com_semente <- function(semente, expr) {
  if (exists(".Random.seed", envir = globalenv())) {
    velho <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", velho, envir = globalenv()), add = TRUE)
  } else on.exit(if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(semente)
  force(expr)
}

ancova <- function(dados, vd, fator, covariaveis, posthoc = c("bonferroni", "holm", "tukey", "nenhum"),
                   exigir_omnibus = TRUE, bootstrap = FALSE, n_boot = 1000, semente = 2026,
                   tipo_ic_boot = c("bca", "perc"), grupos_aleatorizados = NA,
                   alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  posthoc <- match.arg(posthoc); tipo_ic_boot <- match.arg(tipo_ic_boot)
  .pr_exige("car", "emmeans", "effectsize")
  if (length(fator) != 1) stop("ancova() aceita um fator; para mais fatores com covariável, ajustar o modelo à parte")
  dd <- .ee_dados(dados, c(vd, fator, covariaveis), numericas = c(vd, covariaveis))
  d0 <- data.frame(.y = dd$d[[vd]], .g = .ee_fator(dd$d[[fator]], fator))
  ic <- paste0("c", seq_along(covariaveis))
  for (i in seq_along(covariaveis)) d0[[ic[i]]] <- dd$d[[covariaveis[i]]]
  mapa <- stats::setNames(c(fator, covariaveis), c(".g", ic))
  lv <- levels(d0$.g); k <- length(lv); N <- nrow(d0)
  .ee_dif_exige_variacao(d0$.y, vd)
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op), add = TRUE)
  f <- stats::as.formula(paste(".y ~ .g +", paste(ic, collapse = " + ")))
  m <- stats::aov(f, data = d0)
  t3 <- car::Anova(m, type = 3)
  termos <- c(ic, ".g")
  es <- .ee_dif_es_tab(t3, TRUE, nivel); esn <- .ee_dif_es_tab(t3, FALSE, nivel)
  nomes <- unname(mapa[termos])
  tab <- data.frame(termo = nomes, sq = t3[termos, "Sum Sq"], gl1 = t3[termos, "Df"], gl2 = t3["Residuals", "Df"],
                    F = t3[termos, "F value"], p = t3[termos, "Pr(>F)"], stringsAsFactors = FALSE)
  io <- match(termos, es$om$Parameter); ie <- match(termos, es$et$Parameter)
  ef <- do.call(rbind, lapply(seq_along(termos), function(i) .ee_dif_efeito_df(c("ω²p", "η²p"), c(es$om[io[i], 2], es$et[ie[i], 2]),
                                                                              c(es$om$CI_low[io[i]], es$et$CI_low[ie[i]]), c(es$om$CI_high[io[i]], es$et$CI_high[ie[i]]),
                                                                              c(TRUE, FALSE), termo = nomes[i])))
  # inclinações
  fi <- stats::as.formula(paste(".y ~ .g * (", paste(ic, collapse = " + "), ")"))
  mi <- stats::aov(fi, data = d0)
  ti <- car::Anova(mi, type = 3)
  rows_i <- paste0(".g:", ic)
  # inclinação de cada covariável em cada grupo, com as demais covariáveis no modelo (a
  # inclinação parcial que o teste fator × covariável avalia); com uma covariável é a
  # inclinação simples. r = correlação bivariada VD × covariável no grupo.
  fg <- stats::as.formula(paste(".y ~", paste(ic, collapse = " + ")))
  incl <- do.call(rbind, lapply(seq_along(ic), function(j) do.call(rbind, lapply(lv, function(l) {
    s <- d0[d0$.g == l, ]
    data.frame(covariavel = covariaveis[j], grupo = l, n = nrow(s),
               inclinacao = unname(stats::coef(stats::lm(fg, data = s))[ic[j]]),
               tipo = if (length(ic) > 1) "parcial (demais covariáveis no modelo)" else "simples",
               r = stats::cor(s$.y, s[[ic[j]]]), stringsAsFactors = FALSE) }))))
  res <- stats::residuals(m); sw <- .ee_dif_sw(res)
  lev_b <- .ee_dif_levene(d0$.y, d0$.g, "media")
  ar <- stats::anova(stats::lm(abs(res) ~ d0$.g))
  lev_r <- c(F = ar[1, "F value"], gl1 = ar[1, "Df"], gl2 = ar[2, "Df"], p = ar[1, "Pr(>F)"])
  pres <- list()
  for (j in seq_along(ic)) {
    pres[[length(pres) + 1]] <- .ee_dif_pres(paste0("Homogeneidade das inclinações (", fator, " × ", covariaveis[j], ")"), ti[rows_i[j], "F value"],
                                             .ee_dif_gltxt(ti[rows_i[j], "Df"], ti["Residuals", "Df"]), ti[rows_i[j], "Pr(>F)"],
                                             if (ti[rows_i[j], "Pr(>F)"] < alfa) paste0(if (length(ic) > 1) "inclinações parciais diferem entre os grupos (" else "inclinações diferem entre os grupos (", paste(paste0(incl$grupo[incl$covariavel == covariaveis[j]], ": ", vapply(incl$inclinacao[incl$covariavel == covariaveis[j]], .pr_num, "", casas = 2, decimal = decimal)), collapse = "; "), "): a ANCOVA clássica não se sustenta; modelar a interação")
                                             else "sem evidência de inclinações diferentes (não prova que sejam iguais)")
    ag <- stats::anova(stats::lm(d0[[ic[j]]] ~ d0$.g))
    pres[[length(pres) + 1]] <- .ee_dif_pres(paste0("Diferença entre os grupos na covariável ", covariaveis[j], " (ANOVA)"), ag[1, "F value"],
                                             .ee_dif_gltxt(ag[1, "Df"], ag[2, "Df"]), ag[1, "Pr(>F)"],
                                             if (ag[1, "Pr(>F)"] < alfa) "grupos diferem na covariável: se não houve alocação aleatória, o ajuste não corrige a diferença de forma interpretável (Miller & Chapman, 2001)"
                                             else "sem evidência de diferença dos grupos na covariável")
  }
  pres <- c(pres, list(
    .ee_dif_pres("Normalidade dos resíduos (Shapiro-Wilk)", sw[["W"]], "", sw[["p"]], .ee_dif_le_normal(res, "resíduos")),
    .ee_dif_pres("Homogeneidade de variâncias da VD bruta (Levene, centro na média)", lev_b[["F"]], .ee_dif_gltxt(lev_b[["gl1"]], lev_b[["gl2"]]), lev_b[["p"]],
                 paste0(.ee_dif_le_levene(lev_b[["p"]], alfa), "; o valor que o corpo do livro cita")),
    .ee_dif_pres("Homogeneidade de variâncias dos resíduos do modelo (Levene; o do painel do JASP)", lev_r[["F"]], .ee_dif_gltxt(lev_r[["gl1"]], lev_r[["gl2"]]), lev_r[["p"]],
                 .ee_dif_le_levene(lev_r[["p"]], alfa))))
  pres <- do.call(rbind, pres)
  # médias ajustadas e comparações
  em <- emmeans::emmeans(m, ".g")
  ems <- as.data.frame(summary(em, level = nivel))
  medias <- data.frame(grupo = lv, n = as.integer(table(d0$.g)), media_bruta = tapply(d0$.y, d0$.g, mean),
                       media_ajustada = ems$emmean, ep = ems$SE, ic_inf = ems$lower.CL, ic_sup = ems$upper.CL,
                       stringsAsFactors = FALSE, row.names = NULL)
  for (j in seq_along(ic)) medias[[paste0("media_", covariaveis[j])]] <- tapply(d0[[ic[j]]], d0$.g, mean)
  p_fator <- tab$p[tab$termo == fator]
  roda <- posthoc != "nenhum" && k > 2 && (!exigir_omnibus || isTRUE(p_fator < alfa))
  if (posthoc != "nenhum" && k == 2) roda <- TRUE
  ph <- NULL; dec_boot <- character(); avisos <- character()
  if (roda) {
    ph <- .ee_dif_pares_emm(em, posthoc, nivel)
    sa <- .ee_dif_pares_emm(em, "none", nivel)
    ph$ic_inf_sem_ajuste <- sa$ic_inf; ph$ic_sup_sem_ajuste <- sa$ic_sup
    ph$d_qme <- ph$dif / stats::sigma(m)
    if (bootstrap) {
      .pr_exige("boot")
      pr <- utils::combn(k, 2)
      est <- function(dat, idx) {
        b <- dat[idx, , drop = FALSE]
        if (any(table(b$.g) < 2)) return(rep(NA_real_, ncol(pr)))
        mm <- stats::lm(f, data = b)
        nd <- data.frame(.g = factor(lv, levels = lv))
        for (cc in ic) nd[[cc]] <- mean(b[[cc]])
        am <- stats::predict(mm, newdata = nd)
        unname(am[pr[1, ]] - am[pr[2, ]])
      }
      bt <- .ee_dif_com_semente(semente, boot::boot(d0, est, R = n_boot, strata = d0$.g))
      ci_b <- t(vapply(seq_len(ncol(pr)), function(j) {
        r <- tryCatch(suppressWarnings(boot::boot.ci(bt, conf = nivel, type = tipo_ic_boot, index = j)), error = function(e) NULL)
        if (is.null(r) && tipo_ic_boot == "bca") r <- tryCatch(suppressWarnings(boot::boot.ci(bt, conf = nivel, type = "perc", index = j)), error = function(e) NULL)
        if (is.null(r)) return(c(NA_real_, NA_real_))
        v <- if (!is.null(r$bca)) r$bca else r$percent
        v[1, 4:5] }, numeric(2)))
      ph$ic_boot_inf <- ci_b[, 1]; ph$ic_boot_sup <- ci_b[, 2]
      dec_boot <- paste0("Bootstrap (Caixa 15.C): ", n_boot, " reamostragens de casos estratificadas por grupo, semente ", semente,
                         ", IC ", if (tipo_ic_boot == "bca") "BCa (o tipo do JASP; percentil se o BCa falhar)" else "percentil",
                         " das diferenças entre médias ajustadas. Não conserta inclinações heterogêneas (problema de especificação, não de distribuição dos erros).")
    }
  }
  p_incl <- ti[rows_i, "Pr(>F)"]
  if (any(p_incl < alfa)) avisos <- c(avisos, "Homogeneidade das inclinações violada: a ANCOVA clássica é só primeira aproximação; modelar a interação fator × covariável (moderação). O bootstrap não resolve isso.")
  if (isTRUE(lev_r[["p"]] < alfa) && !isTRUE(lev_b[["p"]] < alfa))
    avisos <- c(avisos, "Levene sem evidência na VD bruta, mas com evidência nos resíduos do modelo (o painel do JASP mostra o segundo): dizer qual se cita.")
  viol_dist <- isTRUE(sw[["p"]] < alfa) || isTRUE(max(abs(c(.ee_assimetria(res), .ee_curtose(res)))) > 2) ||
    isTRUE(lev_b[["p"]] < alfa) || isTRUE(lev_r[["p"]] < alfa)
  if (viol_dist && !bootstrap)
    avisos <- c(avisos, "Normalidade dos resíduos ou homogeneidade de variâncias com evidência de violação: a rota de robustez do livro é bootstrap = TRUE com 1000 reamostragens (Caixa 15.C), que dá IC das diferenças ajustadas sem supor normalidade (não conserta inclinações heterogêneas).")
  if (!isTRUE(grupos_aleatorizados)) avisos <- c(avisos, "Sem alocação aleatória declarada (grupos_aleatorizados = TRUE): a ANCOVA reduz um viés, não transforma o desenho em experimento; escrever 'ajustando por' ou 'levando em conta', não 'controlando'.")
  nome_m <- c(bonferroni = "Bonferroni", holm = "Holm", tukey = "Tukey", nenhum = "nenhum")[[posthoc]]
  dec <- c("SQ Tipo III com contrastes de soma zero (padrão do JASP); covariáveis em Covariates, fator em Fixed Factors.",
           paste0("ω²p principal e η²p ao lado. Não parciais do fator, que não se reportam: ω² = ", .pr_num(esn$om[match(".g", esn$om$Parameter), 2], 3, FALSE, decimal),
                  .ee_sep(decimal), "η² = ", .pr_num(esn$et[match(".g", esn$et$Parameter), 2], 3, FALSE, decimal), "."),
           "Leitura: covariável primeiro, depois o fator ajustado, depois médias ajustadas (não as brutas) e comparações.",
           if (roda) paste0("Comparações ", nome_m, " sobre as médias ajustadas (só para o fator, nunca para a covariável); a tabela traz o IC ajustado (o que o JASP mostra) e o não ajustado (o que o corpo do livro cita). d = diferença ajustada / raiz do QM do erro.")
           else if (posthoc != "nenhum") paste0("Fator sem significância após o ajuste (", .pr_p(p_fator, decimal), "): comparações não rodadas."),
           dec_boot,
           "Não extrapolar o ajuste para faixas da covariável em que um grupo quase não tem casos.")
  # frase APA
  i_f <- which(tab$termo == fator)
  efw <- function(i) ef$valor[2 * i - 1]
  st <- function(i) .ee_junta(.ee_estat("F", tab$F[i], c(tab$gl1[i], tab$gl2[i]), decimal), .pr_p(tab$p[i], decimal),
                              .ee_efeito_txt("ω²p", efw(i), NA, NA, TRUE, 3, decimal, idioma, nivel = nivel), decimal = decimal)
  covs_txt <- paste(vapply(seq_along(ic), function(i) paste0(covariaveis[i], " (", st(i), ")"), ""), collapse = if (idioma == "en") " and " else " e ")
  sig <- isTRUE(tab$p[i_f] < alfa)
  # privacidade (H-02): grupo com n < 5 sai sem média ajustada na frase
  peq <- .ee_dif_peq(medias$n)
  if (any(peq)) avisos <- c(avisos, paste0(.ee_dif_aviso_peq(lv[peq], "média ajustada"),
                                           if (!all(peq)) " As diferenças ajustadas das comparações, com a média de outro grupo, ainda permitem deduzi-la." else ""))
  med_txt <- paste(ifelse(peq, paste0(lv, " (", .ee_dif_omit_txt(idioma), ")"),
                          paste0(lv, " = ", vapply(medias$media_ajustada, .ee_dif_m, "", decimal = decimal))), collapse = .ee_sep(decimal))
  # inclinações: separador de itens pelo decimal (G-11; com decimal = ",", "F(2; 24) = 4,18; p = ,028")
  incl_txt <- .ee_junta(.ee_estat("F", ti[rows_i[which.min(p_incl)], "F value"], c(ti[rows_i[which.min(p_incl)], "Df"], ti["Residuals", "Df"]), decimal),
                        .pr_p(min(p_incl), decimal), decimal = decimal)
  if (idioma == "en") {
    apa <- paste0("After statistically adjusting for ", covs_txt, ", ", fator, if (sig) " showed an effect on " else " did not show a statistically significant effect on ", vd, ", ",
                  st(i_f), " (", .ee_dif_rotulo(efw(i_f), "var", "en"), " effect). Adjusted means: ", med_txt, ".")
    if (any(p_incl < alfa)) apa <- c(apa, paste0("The homogeneity-of-regression-slopes assumption did not hold (", fator, " × covariate, ",
                                                incl_txt, "); the classical ANCOVA is a first approximation."))
  } else {
    apa <- paste0("Após ajuste estatístico por ", covs_txt, ", ", fator, if (sig) " apresentou efeito sobre " else " não apresentou efeito estatisticamente significativo sobre ", vd, ", ",
                  st(i_f), " (efeito ", .ee_dif_rotulo(efw(i_f), "var", "pt"), "). Médias ajustadas: ", med_txt, ".")
    if (any(p_incl < alfa)) apa <- c(apa, paste0("O pressuposto de homogeneidade das inclinações não se sustentou (", fator, " × covariável, ",
                                                incl_txt, "); a ANCOVA clássica é lida como primeira aproximação."))
  }
  if (roda) apa <- c(apa, .ee_dif_frase_pares(ph, alfa, nome_m, "d", "d_qme", decimal, idioma))
  avisos <- c(avisos, .ee_aviso_nao_sig(tab$p[i_f], alfa, NULL, "diferença ajustada entre os grupos", NULL, decimal))
  .ee_resultado(paste0("ANCOVA (", fator, "; covariáveis: ", paste(covariaveis, collapse = ", "), ")"), n = N,
                n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef, pressupostos = pres, posthoc = ph,
                decisoes = dec, avisos = avisos, apa = apa,
                extra = list(medias = medias, inclinacoes = incl, modelo = m, modelo_interacao = mi,
                             efeito_nao_parcial = data.frame(termo = unname(mapa[esn$om$Parameter]), omega2 = esn$om[[2]],
                                                             eta2 = esn$et[[2]][match(esn$om$Parameter, esn$et$Parameter)], stringsAsFactors = FALSE)))
}

# ================================================================ MANOVA
# Estatísticas multivariadas a partir dos autovalores de E^-1 H, com as aproximações F
# do stats::summary.manova (as que o JASP mostra).
.ee_dif_mv_estat <- function(eig, q, gl_res) {
  p <- length(eig); s <- min(p, q)
  m <- .5 * (abs(p - q) - 1); nn <- .5 * (gl_res - p - 1)
  V <- sum(eig / (1 + eig)); t1 <- 2 * m + s + 1; t2 <- 2 * nn + s + 1
  pil <- c(V, (t2 / t1 * V) / (s - V), s * t1, s * t2)
  L <- prod(1 / (1 + eig)); a1 <- gl_res - .5 * (p - q + 1); a2 <- (p * q - 2) / 4; a3 <- p^2 + q^2 - 5
  a3 <- if (a3 > 0) sqrt(((p * q)^2 - 4) / a3) else 1
  wil <- c(L, ((L^(-1 / a3) - 1) * (a1 * a3 - 2 * a2)) / p / q, p * q, a1 * a3 - 2 * a2)
  HL <- sum(eig); h2 <- 2 * (s * nn + 1)
  hot <- c(HL, (h2 * HL) / s / s / t1, s * t1, h2)
  R <- max(eig); r1 <- max(p, q); r2 <- gl_res - r1 + q
  roy <- c(R, (r2 * R) / r1, r1, r2)
  out <- rbind(Pillai = pil, Wilks = wil, "Hotelling-Lawley" = hot, Roy = roy)
  data.frame(estatistica = rownames(out), valor = out[, 1], F = out[, 2], gl1 = out[, 3], gl2 = out[, 4],
             p = stats::pf(out[, 2], out[, 3], out[, 4], lower.tail = FALSE), stringsAsFactors = FALSE, row.names = NULL)
}

# Shapiro-Wilk multivariado de Villasenor-Alva & González-Estrada (2009), o teste que o
# JASP 0.95.1 usa (pacote mvShapiroTest), sobre as VDs brutas.
.ee_dif_sw_mv <- function(X) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  if (n < 12 || n > 5000 || n <= p) return(c(W = NA_real_, p = NA_real_))
  x <- scale(X, scale = FALSE)
  ev <- eigen(stats::var(X), symmetric = TRUE)
  raiz <- ev$vectors %*% diag(1 / sqrt(ev$values), ncol = p) %*% t(ev$vectors)
  z <- t(raiz %*% t(x))
  w <- mean(vapply(seq_len(p), function(k) unname(stats::shapiro.test(z[, k])$statistic), 1))
  y <- log(n); w1 <- log(1 - w)
  mm <- -1.5861 - .31082 * y - 0.083751 * y^2 + .0038915 * y^3
  ss <- exp(-.4803 - .082676 * y + .0030302 * y^2)
  s2 <- log((p - 1 + exp(ss^2)) / p)
  c(W = w, p = stats::pnorm(w1, mm + ss^2 / 2 - s2 / 2, sqrt(s2), lower.tail = FALSE))
}

# Henze-Zirkler (a fórmula do pacote MVN, função hz): estatística e p pela aproximação
# lognormal. Aplicado aos resíduos da MANOVA, que são o objeto do pressuposto.
.ee_dif_hz <- function(X) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  if (n <= p + 1) return(c(HZ = NA_real_, p = NA_real_))
  dif <- scale(X, scale = FALSE)
  S <- ((n - 1) / n) * stats::cov(X)
  Si <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(Si)) return(c(HZ = NA_real_, p = NA_real_))
  Dj <- diag(dif %*% Si %*% t(dif))
  Y <- X %*% Si %*% t(X)
  Djk <- -2 * t(Y) + matrix(diag(t(Y))) %*% matrix(rep(1, n), 1, n) + matrix(rep(1, n), n, 1) %*% matrix(diag(t(Y)), 1, n)
  b <- 1 / sqrt(2) * ((2 * p + 1) / 4)^(1 / (p + 4)) * n^(1 / (p + 4))
  hz <- n * (1 / n^2 * sum(exp(-(b^2) / 2 * Djk)) - 2 * ((1 + b^2)^(-p / 2)) * (1 / n) * sum(exp(-((b^2) / (2 * (1 + b^2))) * Dj)) +
               ((1 + 2 * b^2)^(-p / 2)))
  wb <- (1 + b^2) * (1 + 3 * b^2); a <- 1 + 2 * b^2
  mu <- 1 - a^(-p / 2) * (1 + p * b^2 / a + (p * (p + 2) * b^4) / (2 * a^2))
  si2 <- 2 * (1 + 4 * b^2)^(-p / 2) + 2 * a^(-p) * (1 + (2 * p * b^4) / a^2 + (3 * p * (p + 2) * b^8) / (4 * a^4)) -
    4 * wb^(-p / 2) * (1 + (3 * p * b^4) / (2 * wb) + (p * (p + 2) * b^8) / (2 * wb^2))
  pmu <- log(sqrt(mu^4 / (si2 + mu^2))); psi <- sqrt(log((si2 + mu^2) / mu^2))
  c(HZ = hz, p = stats::plnorm(hz, pmu, psi, lower.tail = FALSE))
}

# Leitura do M de Box com um critério único (o usual é p < .001, pela sensibilidade do
# teste com n grande; o livro não declara critério e só mostra o caso p < ,001).
.ee_dif_le_box <- function(p, alfa_box, decimal = ".") {
  crit <- paste0("critério p < ", .pr_num(alfa_box, 3, FALSE, decimal))
  if (!is.finite(p)) return("não calculável: algum grupo tem n menor ou igual ao número de variáveis")
  if (p < alfa_box) paste0("heterogeneidade das matrizes (", crit, ")") else
    paste0("sem evidência de heterogeneidade pelo ", crit, " (o M de Box é muito sensível com n grande)")
}

# MANOVA. Pillai é a estatística principal, fixada antes de ver o resultado (Olson, 1976;
#); Wilks, Hotelling-Lawley e Roy em tabela. Tamanho de efeito
# multivariado: η²p = V/s (o "aproximadamente ,17" do livro). ANOVAs
# univariadas por VD com η²p, interpretadas só se o multivariado for significativo
#. tipo_sq = 1 reproduz o JASP (summary.manova é sequencial); 3 (padrão) usa
# car::Manova; com um fator os dois coincidem.
manova_grupos <- function(dados, vds, fatores, tipo_sq = 3, alfa_box = .001, alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  if (length(vds) < 2) stop("a MANOVA pede 2 ou mais VDs")
  if (!tipo_sq %in% c(1, 3)) stop("tipo_sq deve ser 1 (sequencial, o do JASP) ou 3")
  .pr_exige("car", "heplots")
  dd <- .ee_dados(dados, c(vds, fatores), numericas = vds)
  ifs <- paste0("f", seq_along(fatores)); mapa <- stats::setNames(fatores, ifs)
  d0 <- dd$d[, vds, drop = FALSE]; names(d0) <- paste0("y", seq_along(vds))
  for (i in seq_along(fatores)) d0[[ifs[i]]] <- .ee_fator(dd$d[[fatores[i]]], fatores[i])
  Y <- as.matrix(d0[paste0("y", seq_along(vds))]); N <- nrow(d0)
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op), add = TRUE)
  rhs <- paste(ifs, collapse = " * ")
  fit <- stats::lm(stats::as.formula(paste("Y ~", rhs)), data = d0)
  if (tipo_sq == 1) {
    sm <- summary(stats::manova(fit), test = "Pillai")
    SS <- sm$SS; E <- SS$Residuals; termos <- setdiff(names(SS), "Residuals")
    gls <- stats::setNames(sm$stats[termos, "Df"], termos); gl_res <- fit$df.residual
    H <- SS[termos]
  } else {
    mv <- car::Manova(fit, type = 3)
    termos <- setdiff(names(mv$SSP), "(Intercept)")
    H <- mv$SSP[termos]; E <- mv$SSPE; gls <- unlist(mv$df[termos]); gl_res <- mv$error.df
  }
  termos <- gsub(" ", "", termos)
  names(H) <- termos; names(gls) <- termos
  nomes <- .ee_dif_nomes(termos, mapa)
  mvt <- do.call(rbind, lapply(seq_along(termos), function(i) {
    eig <- Re(eigen(solve(E) %*% H[[i]], only.values = TRUE)$values)
    eig <- pmax(eig, 0)  # todos os p autovalores (os nulos entram): p = nº de VDs nas aproximações F
    cbind(termo = nomes[i], .ee_dif_mv_estat(eig, gls[[i]], gl_res), stringsAsFactors = FALSE)
  }))
  mvt$principal <- mvt$estatistica == "Pillai"
  s_t <- vapply(termos, function(t) min(length(vds), gls[[t]]), 1)
  pil <- mvt[mvt$principal, ]
  wil <- mvt[mvt$estatistica == "Wilks", ]
  ef <- rbind(.ee_dif_efeito_df(rep("η²p multivariado (Pillai, V/s)", length(termos)), pil$valor / s_t, NA, NA, TRUE, termo = nomes),
              .ee_dif_efeito_df(rep("η²p multivariado (Wilks, 1 - Λ^(1/s))", length(termos)), 1 - wil$valor^(1 / s_t), NA, NA, FALSE, termo = nomes))
  # univariadas
  uni <- do.call(rbind, lapply(seq_along(vds), function(j) {
    mj <- stats::aov(stats::as.formula(paste0("y", j, " ~ ", rhs)), data = d0)
    tj <- if (tipo_sq == 1) { a <- stats::anova(mj); a } else car::Anova(mj, type = 3)
    tr <- setdiff(gsub(" ", "", rownames(tj)), c("(Intercept)", "Residuals"))
    rn <- gsub(" ", "", rownames(tj))
    sqr <- tj[rn == "Residuals", "Sum Sq"]; glr <- tj[rn == "Residuals", "Df"]
    data.frame(vd = vds[j], termo = .ee_dif_nomes(tr, mapa), sq = tj[match(tr, rn), "Sum Sq"], gl1 = tj[match(tr, rn), "Df"], gl2 = glr,
               F = tj[match(tr, rn), "F value"], p = tj[match(tr, rn), "Pr(>F)"],
               eta2p = tj[match(tr, rn), "Sum Sq"] / (tj[match(tr, rn), "Sum Sq"] + sqr), stringsAsFactors = FALSE)
  }))
  ef <- rbind(ef, .ee_dif_efeito_df(rep("η²p univariado", nrow(uni)), uni$eta2p, NA, NA, FALSE, termo = paste0(uni$vd, ": ", uni$termo)))
  cel <- interaction(d0[ifs], drop = TRUE)
  box <- tryCatch(suppressWarnings(heplots::boxM(Y, cel)), error = function(e) NULL)
  swm <- .ee_dif_sw_mv(Y)
  hz <- .ee_dif_hz(stats::residuals(fit))
  if (!is.null(box) && !is.finite(box$p.value)) box <- NULL  # M de Box não calculável (grupo com n <= nº de VDs)
  box_viol <- !is.null(box) && isTRUE(box$p.value < alfa_box)
  pres <- list()
  if (!is.null(box)) pres[[1]] <- .ee_dif_pres("Homogeneidade das matrizes de covariância (M de Box, χ² aproximado)", unname(box$statistic), .ee_gl(unname(box$parameter)), box$p.value,
                                               paste0(.ee_dif_le_box(box$p.value, alfa_box),
                                                      if (box_viol) ": manter o traço de Pillai, o mais robusto, e registrar a limitação" else ""))
  pres[[length(pres) + 1]] <- .ee_dif_pres("Normalidade multivariada dos resíduos (Henze-Zirkler)", hz[["HZ"]], "", hz[["p"]],
                                           paste0(if (isTRUE(hz[["p"]] < alfa)) "evidência de afastamento da normalidade multivariada" else "sem evidência de afastamento (não prova normalidade)",
                                                  "; teste sobre os resíduos, o objeto do pressuposto; com n pequeno nenhum teste é palavra final"))
  pres[[length(pres) + 1]] <- .ee_dif_pres("Normalidade multivariada (Shapiro-Wilk multivariado sobre as VDs brutas: o que o painel do JASP mostra)", swm[["W"]], "", swm[["p"]],
                                           "referência para conferir a tela; sobre as VDs brutas, mistura a diferença entre os grupos com a forma da distribuição")
  for (j in seq_along(vds)) {
    lv <- .ee_dif_levene(d0[[paste0("y", j)]], cel, "media")
    pres[[length(pres) + 1]] <- .ee_dif_pres(paste0("Homogeneidade de variâncias em ", vds[j], " (Levene, centro na média)"), lv[["F"]], .ee_dif_gltxt(lv[["gl1"]], lv[["gl2"]]), lv[["p"]], .ee_dif_le_levene(lv[["p"]], alfa))
  }
  pres <- do.call(rbind, pres)
  sig_mv <- pil$p < alfa
  dec <- c("Pillai é a estatística principal, fixada antes de ver o resultado (Olson, 1976); nunca escolher a estatística pelo menor p (Roy é a mais propensa a erro tipo I).",
           paste0("SQ ", if (tipo_sq == 1) "sequencial (Tipo I), a do summary.manova e a do JASP" else "Tipo III (car::Manova)", if (length(fatores) > 1)
             "; com mais de um fator e desenho desbalanceado, Tipo I e Tipo III diferem (o JASP usa Tipo I)." else "; com um fator, Tipo I e Tipo III coincidem."),
           "F de Hotelling-Lawley pela aproximação do summary.manova (a do JASP: F = 3,101, gl 50 no Viagra); o corpo do livro cita a de McKeon (3,19; 50,19). O traço é o número estável entre programas.",
           "Tamanho de efeito multivariado: η²p = V/s (s = mín(nº de VDs, gl do efeito)); com Pillai = ,349 e s = 2 dá o ,17 do livro.",
           paste0("M de Box lido com o critério p < ", alfa_box, " (alfa_box), o usual pela sensibilidade do teste; o livro não declara critério (só mostra p < ,001)."),
           "Normalidade multivariada: Henze-Zirkler sobre os resíduos (o objeto do pressuposto); o Shapiro-Wilk multivariado do JASP, sobre as VDs brutas, fica para conferir a tela.",
           if (any(sig_mv)) "Multivariado significativo: as ANOVAs univariadas podem ser lidas; multivariado significativo com univariadas nulas indica efeito na combinação linear das VDs (passo seguinte: análise discriminante)."
           else "Multivariado sem significância: as ANOVAs univariadas ficam na saída só por transparência e não se interpretam.",
           "A 'proteção' das univariadas pela MANOVA só vale se a hipótese nula multivariada for toda verdadeira; quem quer as VDs isoladas pode rodar ANOVAs com correção (curso, seção MANOVA).")
  avisos <- character()
  if (box_viol && length(unique(table(cel))) > 1)
    avisos <- "M de Box violado com grupos de tamanhos diferentes: a robustez do Pillai cai nessa combinação; registrar a limitação."
  frases <- vapply(seq_len(nrow(pil)), function(i) {
    st <- .ee_junta(paste0("V = ", .pr_num(pil$valor[i], 2, pil$valor[i] >= 1, decimal)),
                    .ee_estat("F", pil$F[i], c(pil$gl1[i], pil$gl2[i]), decimal), .pr_p(pil$p[i], decimal),
                    .ee_efeito_txt("η²p", pil$valor[i] / s_t[i], NA, NA, TRUE, 3, decimal, idioma, nivel = nivel), decimal = decimal)
    if (idioma == "en") paste0("The multivariate effect of ", nomes[i], if (sig_mv[i]) " was statistically significant, " else " was not statistically significant, ", st, ".")
    else paste0("O efeito multivariado de ", nomes[i], if (sig_mv[i]) " foi estatisticamente significativo, " else " não foi estatisticamente significativo, ", st, ".")
  }, "")
  crit_box <- .pr_num(alfa_box, 3, FALSE, decimal)
  abre <- if (is.null(box)) {
    if (idioma == "en") "Pillai's trace was used. " else "Adotou-se o traço de Pillai. "
  } else if (box_viol) {
    if (idioma == "en") paste0("Because Box's M indicated heterogeneous covariance matrices, ", .pr_p(box$p.value, decimal), " (criterion p < ", crit_box, "), Pillai's trace was used. ")
    else paste0("Dado que o teste M de Box indicou heterogeneidade das matrizes de covariância, ", .pr_p(box$p.value, decimal), " (critério p < ", crit_box, "), adotou-se o traço de Pillai. ")
  } else {
    if (idioma == "en") paste0("Box's M did not indicate heterogeneous covariance matrices by the p < ", crit_box, " criterion (", .pr_p(box$p.value, decimal), "); Pillai's trace, chosen a priori, was used. ")
    else paste0("O teste M de Box não indicou heterogeneidade das matrizes de covariância pelo critério p < ", crit_box, " (", .pr_p(box$p.value, decimal), "); adotou-se o traço de Pillai, fixado a priori. ")
  }
  apa <- paste0(abre, paste(frases, collapse = " "))
  if (any(sig_mv)) {
    u <- uni[uni$termo %in% nomes[sig_mv], ]
    ut <- paste(vapply(seq_len(nrow(u)), function(i) paste0(u$vd[i], if (length(fatores) > 1) paste0(" (", u$termo[i], ")") else "", ": ",
                                                           .ee_junta(.ee_estat("F", u$F[i], c(u$gl1[i], u$gl2[i]), decimal), .pr_p(u$p[i], decimal),
                                                                     .ee_efeito_txt("η²p", u$eta2p[i], NA, NA, TRUE, 3, decimal, idioma, nivel = nivel), decimal = decimal)), ""), collapse = "; ")
    apa <- c(apa, paste0(if (idioma == "en") "Univariate ANOVAs: " else "ANOVAs univariadas: ", ut, "."))
    if (all(u$p >= alfa)) apa <- c(apa, if (idioma == "en") "No univariate ANOVA was significant: the effect lies in the linear combination of the outcomes." else
      "Nenhuma ANOVA univariada foi significativa: o efeito está na combinação linear dos desfechos, não em um deles isolado.")
  }
  avisos <- c(avisos, .ee_aviso_nao_sig(pil$p, alfa, NULL, "efeito multivariado", NULL, decimal,
                                        rotulos = paste0(pil$termo, " (Pillai, ", vapply(pil$p, .pr_p, "", decimal = decimal), ")")))
  .ee_resultado(paste0("MANOVA (", paste(vds, collapse = ", "), " por ", paste(fatores, collapse = " × "), ")"), n = N,
                n_excluidos = dd$n_excluidos, tabela = mvt, efeito = ef, pressupostos = pres, decisoes = dec, avisos = avisos,
                apa = apa, extra = list(univariadas = uni, box_m = box, modelo = fit))
}

# ================================================================ Kruskal-Wallis e Friedman
# Kruskal-Wallis com ε²_R = H/(N - 1) (principal; Tomczak & Tomczak) e η²_H = (H - k + 1)/(N - k)
# (o "rank η²" do JASP, usado no curso); IC do ε² por bootstrap percentil estratificado
# (n_boot = 1000, semente). Pós-hoc de Dunn (o do JASP) só com ômnibus significativo,
# com p de Bonferroni e de Holm na tabela; correcao decide qual vale (Holm por padrão:
# controla o mesmo erro de família e nunca tem menos poder; a regra do curso, Bonferroni
# até 3 grupos, fica como opção). r_rb por par vem do Mann-Whitney do par.
teste_kruskal_wallis <- function(dados, vd, grupo, correcao = c("holm", "bonferroni"), exigir_omnibus = TRUE,
                                 sinal_rrb = c("direcional", "jasp_0951"),
                                 n_boot = 1000, semente = 2026, alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  correcao <- match.arg(correcao); sinal_rrb <- match.arg(sinal_rrb)
  dd <- .ee_dados(dados, c(vd, grupo), numericas = vd)
  y <- dd$d[[vd]]; g <- .ee_fator(dd$d[[grupo]], grupo)
  k <- nlevels(g); N <- length(y); nv <- levels(g)
  if (any(table(g) < 2)) stop("cada grupo precisa de pelo menos 2 casos")
  .ee_dif_exige_variacao(y, vd)  # VD constante: H indefinido (0/0); constante só dentro dos grupos é separação perfeita e o teste roda
  kt <- stats::kruskal.test(y, g)
  H <- unname(kt$statistic)
  eps <- H / (N - 1); eta <- (H - k + 1) / (N - k)
  boot_eps <- .ee_dif_com_semente(semente, replicate(n_boot, {
    yb <- unlist(lapply(nv, function(l) { x <- y[g == l]; x[sample.int(length(x), replace = TRUE)] }))
    gb <- factor(rep(nv, as.integer(table(g))), levels = nv)
    unname(stats::kruskal.test(yb, gb)$statistic) / (N - 1)
  }))
  ic_e <- stats::quantile(boot_eps, c((1 - nivel) / 2, 1 - (1 - nivel) / 2), na.rm = TRUE, names = FALSE)
  ef <- .ee_dif_efeito_df(c("ε²_R (de postos)", "η²_H (rank η² do JASP)"), c(eps, eta), c(ic_e[1], NA), c(ic_e[2], NA), c(TRUE, FALSE))
  rk <- rank(y); emp <- table(rk)
  desc <- data.frame(grupo = nv, n = as.integer(table(g)), mediana = tapply(y, g, stats::median), media = tapply(y, g, mean),
                     dp = tapply(y, g, stats::sd), media_postos = tapply(rk, g, mean), stringsAsFactors = FALSE, row.names = NULL)
  tab <- data.frame(teste = "Kruskal-Wallis", H = H, gl = unname(kt$parameter), p = kt$p.value, stringsAsFactors = FALSE)
  roda <- k > 2 && (!exigir_omnibus || kt$p.value < alfa)
  ph <- NULL
  if (roda) {
    s2 <- N * (N + 1) / 12 - sum(emp^3 - emp) / (12 * (N - 1))
    pr <- utils::combn(k, 2)
    ph <- do.call(rbind, lapply(seq_len(ncol(pr)), function(j) {
      a <- pr[1, j]; b <- pr[2, j]
      ya <- y[g == nv[a]]; yb <- y[g == nv[b]]
      z <- (desc$media_postos[a] - desc$media_postos[b]) / sqrt(s2 * (1 / desc$n[a] + 1 / desc$n[b]))
      U <- unname(suppressWarnings(stats::wilcox.test(ya, yb))$statistic)
      r <- 2 * U / (length(ya) * length(yb)) - 1
      data.frame(comparacao = paste(nv[a], "-", nv[b]), z = z, p_sem_ajuste = 2 * stats::pnorm(-abs(z)),
                 rrb = if (sinal_rrb == "jasp_0951") -r else r, stringsAsFactors = FALSE)
    }))
    ph$p_bonferroni <- stats::p.adjust(ph$p_sem_ajuste, "bonferroni")
    ph$p_holm <- stats::p.adjust(ph$p_sem_ajuste, "holm")
    ph$p_ajustado <- if (correcao == "holm") ph$p_holm else ph$p_bonferroni
  }
  dec <- c("ε²_R é o tamanho de efeito principal (Tomczak & Tomczak); η²_H, o 'rank η²' que o JASP mostra e o curso usa, sai ao lado. O curso prefere não rotular esses índices, e a frase segue isso.",
           paste0("IC do ε² por bootstrap percentil estratificado por grupo, ", n_boot, " reamostragens, semente ", semente, " (o IC do effectsize também é bootstrap e varia entre execuções)."),
           if (roda) paste0("Pós-hoc de Dunn (o do JASP: z com variância agrupada dos postos e correção de empates) com p de Bonferroni e de Holm; vale ", if (correcao == "holm") "Holm" else "Bonferroni",
                            ". Holm controla o mesmo erro de família e é sempre ao menos tão poderoso quanto Bonferroni (Holm, 1979); a regra do curso (Bonferroni até 3 grupos) fica em correcao = 'bonferroni'.")
           else if (k > 2) paste0("Ômnibus sem significância (", .pr_p(kt$p.value, decimal), "): sem pós-hoc."),
           if (roda) paste0("r_rb por par vem do Mann-Whitney do par (método diferente do z de Dunn, como no JASP); ", if (sinal_rrb == "direcional") "positivo = 1º grupo do par com valores maiores (o JASP 0.95.1 mostra o sinal oposto)." else "sinal do JASP 0.95.1."),
           "Kruskal-Wallis testa dominância estocástica, não medianas: relatar n, mediana, média, DP e média dos postos.")
  sig <- kt$p.value < alfa
  st <- .ee_junta(.ee_estat("H", H, unname(kt$parameter), decimal), .pr_p(kt$p.value, decimal),
                  .ee_efeito_txt("ε²_R", eps, ic_e[1], ic_e[2], TRUE, 3, decimal, idioma, nivel = nivel), decimal = decimal)
  apa <- if (idioma == "en") paste0("A Kruskal-Wallis test ", if (sig) "indicated" else "did not indicate", " differences in ", vd, " across ", grupo, ", ", st, ".")
  else paste0("Um teste de Kruskal-Wallis ", if (sig) "indicou" else "não indicou", " diferença em ", vd, " entre os níveis de ", grupo, ", ", st, ".")
  if (roda) {
    ph2 <- ph; ph2$dif <- ph2$z
    sgn <- ph2[ph2$p_ajustado < alfa, , drop = FALSE]
    nome <- if (correcao == "holm") "Holm" else "Bonferroni"
    apa <- c(apa, if (!nrow(sgn)) {
      if (idioma == "en") paste0("Dunn comparisons (", nome, "): no pair differed after correction.") else paste0("Nas comparações de Dunn (", nome, "), nenhum par diferiu após a correção.")
    } else paste0(if (idioma == "en") "Dunn comparisons (" else "Comparações de Dunn (", nome, "): ",
                  paste(vapply(seq_len(nrow(sgn)), function(i) paste0(sgn$comparacao[i], ": z = ", .pr_num(sgn$z[i], 2, TRUE, decimal), .ee_sep(decimal),
                                                                       .pr_p(sgn$p_ajustado[i], decimal), .ee_sep(decimal), "r_rb = ", .pr_num(sgn$rrb[i], 2, FALSE, decimal)), ""), collapse = "; "), "."))
  }
  .ee_resultado("Teste de Kruskal-Wallis", n = N, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef, posthoc = ph,
                decisoes = dec, avisos = .ee_aviso_nao_sig(kt$p.value, alfa, NULL, "diferença entre os grupos", NULL, decimal),
                apa = apa, extra = list(descritivas = desc))
}

# Friedman para 3+ medidas do mesmo grupo. W de Kendall = χ²/(n(k - 1)) principal (o
# tamanho de efeito usual do Friedman: concordância dos participantes na ordenação das
# condições). Pós-hoc de Conover (o do JASP), só com ômnibus significativo, p de
# Bonferroni e de Holm; r_rb por par do Wilcoxon pareado ("A - B": positivo = A maior).
teste_friedman <- function(dados, medidas, nomes = NULL, correcao = c("holm", "bonferroni"), exigir_omnibus = TRUE,
                           alfa = .05, nivel = .95, decimal = ".", idioma = "pt") {
  correcao <- match.arg(correcao)
  if (length(medidas) < 3) stop("o Friedman pede 3 ou mais medidas; com 2, use teste_wilcoxon_pareado()")
  if (is.null(nomes)) nomes <- medidas
  dd <- .ee_dados(dados, medidas, numericas = medidas)
  M <- as.matrix(dd$d[, medidas]); colnames(M) <- nomes
  n <- nrow(M); k <- ncol(M)
  ft <- stats::friedman.test(M)
  chi <- unname(ft$statistic); W <- chi / (n * (k - 1))
  R <- t(apply(M, 1, rank)); Rj <- colSums(R); A1 <- sum(R^2)
  den <- sqrt(2 * (n * A1 - sum(Rj^2)) / ((n - 1) * (k - 1)))
  desc <- data.frame(medida = nomes, n = n, mediana = apply(M, 2, stats::median), media = colMeans(M), dp = apply(M, 2, stats::sd),
                     media_postos = Rj / n, stringsAsFactors = FALSE, row.names = NULL)
  tab <- data.frame(teste = "Friedman", qui2 = chi, gl = unname(ft$parameter), p = ft$p.value, stringsAsFactors = FALSE)
  ef <- .ee_dif_efeito_df("W de Kendall", W, NA, NA, TRUE)
  roda <- !exigir_omnibus || ft$p.value < alfa
  ph <- NULL
  if (roda) {
    pr <- utils::combn(k, 2); glc <- (n - 1) * (k - 1)
    ph <- do.call(rbind, lapply(seq_len(ncol(pr)), function(j) {
      a <- pr[1, j]; b <- pr[2, j]
      tt <- (Rj[[a]] - Rj[[b]]) / den
      dv <- M[, a] - M[, b]; dv <- dv[dv != 0]; nd <- length(dv)
      rr <- if (nd) { rkd <- rank(abs(dv)); 2 * sum(rkd[dv > 0]) / (nd * (nd + 1) / 2) - 1 } else NA_real_
      data.frame(comparacao = paste(nomes[a], "-", nomes[b]), t_conover = tt, gl = glc,
                 p_sem_ajuste = 2 * stats::pt(-abs(tt), glc), rrb = rr, stringsAsFactors = FALSE)
    }))
    ph$p_bonferroni <- stats::p.adjust(ph$p_sem_ajuste, "bonferroni")
    ph$p_holm <- stats::p.adjust(ph$p_sem_ajuste, "holm")
    ph$p_ajustado <- if (correcao == "holm") ph$p_holm else ph$p_bonferroni
  }
  dec <- c("W de Kendall = χ²/(n(k - 1)) é o tamanho de efeito do Friedman: mede a concordância dos participantes na ordenação das condições (W alto = mudança consistente). O curso recomenda não interpretá-lo em intervenções; essa leitura foi revista na extração do curso e o W é mantido, com os r_rb por par ao lado.",
           if (roda) paste0("Pós-hoc de Conover (o do JASP): t = diferença das somas de postos / sqrt(2(n·A1 - ΣR²)/((n - 1)(k - 1))), gl = (n - 1)(k - 1); vale ", if (correcao == "holm") "Holm" else "Bonferroni",
                            " (Holm nunca tem menos poder; o curso aceita Bonferroni com 3 tempos). r_rb do Wilcoxon pareado de cada par, positivo = 1ª medida do par maior (a convenção do JASP).")
           else paste0("Ômnibus sem significância (", .pr_p(ft$p.value, decimal), "): sem pós-hoc."),
           "Casos com qualquer medida ausente saem (listwise).")
  sig <- ft$p.value < alfa
  st <- .ee_junta(.ee_estat("χ²", chi, unname(ft$parameter), decimal), .pr_p(ft$p.value, decimal),
                  .ee_efeito_txt("W", W, NA, NA, TRUE, 2, decimal, idioma, nivel = nivel), decimal = decimal)
  rotW <- .ee_dif_rotulo(W, "r", idioma)
  apa <- if (idioma == "en") paste0("A Friedman test ", if (sig) "indicated" else "did not indicate", " differences among ", paste(nomes, collapse = ", "), " (n = ", n, "), ", st, " (Kendall's W; ", rotW, ").")
  else paste0("Um teste de Friedman ", if (sig) "indicou" else "não indicou", " diferença entre ", paste(nomes, collapse = ", "), " (n = ", n, "), ", st, " (W de Kendall; ", rotW, ").")
  if (roda) {
    sgn <- ph[ph$p_ajustado < alfa, , drop = FALSE]
    nome <- if (correcao == "holm") "Holm" else "Bonferroni"
    apa <- c(apa, if (!nrow(sgn)) { if (idioma == "en") paste0("Conover comparisons (", nome, "): no pair differed after correction.") else paste0("Nas comparações de Conover (", nome, "), nenhum par diferiu após a correção.") }
             else paste0(if (idioma == "en") "Conover comparisons (" else "Comparações de Conover (", nome, "): ",
                         paste(vapply(seq_len(nrow(sgn)), function(i) paste0(sgn$comparacao[i], ": t(", sgn$gl[i], ") = ", .pr_num(sgn$t_conover[i], 2, TRUE, decimal), .ee_sep(decimal),
                                                                            .pr_p(sgn$p_ajustado[i], decimal), .ee_sep(decimal), "r_rb = ", .pr_num(sgn$rrb[i], 2, FALSE, decimal)), ""), collapse = "; "), "."))
  }
  .ee_resultado("Teste de Friedman", n = n, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef, posthoc = ph,
                decisoes = dec, avisos = .ee_aviso_nao_sig(ft$p.value, alfa, NULL, "diferença entre as medidas", NULL, decimal),
                apa = apa, extra = list(descritivas = desc))
}
