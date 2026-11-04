# est_poder.R: poder e tamanho amostral a priori, reproduzindo o G*Power 3.1, e (v1.4) o
# planejamento para AFC e MEE do cap. 5 do Wiley Handbook (RMSEA, aninhados e AIPE).
# Carregar por scripts/carregar.R (usa .pr_num, .pr_exige, .ee_resultado e afins).
#
# Uma função por família de teste, todas com a mesma lógica:
#   - sem `n`: resolve o N mínimo para o `poder` pedido (análise a priori);
#   - com `n`: devolve o poder alcançado (post hoc; com aviso, ver livro);
#   - imprime as configurações para refazer a conta no G*Power (Test family,
#     Statistical test, campos de entrada e de saída), o N mínimo inteiro e, quando o
#     desenho pede, o N arredondado para grupos iguais, que é o que o G*Power mostra
#     na ANOVA de um fator e nas medidas repetidas;
#   - `perda` (proporção prevista de perda) acrescenta N / (1 - perda), ajuste que não
#     vem das fontes do usuário (ver avisos e references/poder-amostra.md).
# Motor: distribuições não centrais do R (stats::pf, stats::pt, stats::pchisq) com os
# parâmetros do manual do G*Power; correlação exata e logística vêm do pwrss 1.3.1
# (power.exact.onecor e power.z.logistic). A equivalência com pwr e pwrss está nos evals.

# ---------------------------------------------------------------- motor
.ee_poder_f <- function(lambda, gl1, gl2, alfa) {
  if (gl2 <= 0) return(list(poder = 0, critico = NA_real_))
  crit <- stats::qf(1 - alfa, gl1, gl2)
  list(poder = stats::pf(crit, gl1, gl2, ncp = lambda, lower.tail = FALSE), critico = crit)
}

.ee_poder_t <- function(delta, gl, alfa, caudas) {
  if (gl <= 0) return(list(poder = 0, critico = NA_real_))
  if (caudas == 2) {
    crit <- stats::qt(1 - alfa / 2, gl)
    p <- stats::pt(crit, gl, ncp = delta, lower.tail = FALSE) + stats::pt(-crit, gl, ncp = delta)
  } else {
    crit <- stats::qt(1 - alfa, gl)
    p <- stats::pt(crit, gl, ncp = delta, lower.tail = FALSE)
  }
  list(poder = p, critico = crit)
}

# Menor N = inicio + passo * j (j >= 0) com poder >= alvo. O poder cresce com N em
# todos os desenhos daqui, então busca exponencial seguida de bisseção nos inteiros.
.ee_poder_busca <- function(fp, alvo, inicio, passo = 1, n_max = 1e7) {
  if (fp(inicio) >= alvo) return(inicio)
  lo <- 0; hi <- 1
  while (fp(inicio + passo * hi) < alvo) {
    lo <- hi; hi <- hi * 2
    if (inicio + passo * hi > n_max) stop("nenhum N até ", format(n_max, big.mark = "."),
                                          " alcança o poder pedido: conferir o tamanho de efeito")
  }
  while (hi - lo > 1) {
    mid <- (lo + hi) %/% 2
    if (fp(inicio + passo * mid) >= alvo) hi <- mid else lo <- mid
  }
  inicio + passo * hi
}

.ee_poder_checa <- function(poder, alfa, perda, caudas = 2, n = NULL) {
  if (is.null(n) && (is.null(poder) || poder <= 0 || poder >= 1)) stop("poder precisa estar entre 0 e 1 (ex.: 0.80)")
  if (alfa <= 0 || alfa >= 1) stop("alfa precisa estar entre 0 e 1 (ex.: 0.05)")
  if (perda < 0 || perda >= 1) stop("perda é uma proporção entre 0 e 1 (ex.: 0.20 para 20%)")
  if (!caudas %in% c(1, 2)) stop("caudas precisa ser 1 ou 2")
  invisible(TRUE)
}

.ee_poder_fmt <- function(x, casas = 4) {
  if (is.character(x)) return(x)
  vapply(x, function(v) if (is.na(v)) "NA" else format(round(v, casas), nsmall = 0, scientific = FALSE, trim = TRUE), "")
}

# Bloco "G*Power, saída" com as casas do programa: λ, δ, valores críticos e poder com 6
# casas fixas (manual, seção 11.3: "Power (1-β err prob): 0.99999", que com 4 casas
# virava "1"); gl e N, inteiros, sem casas.
.ee_poder_fmt_saida <- function(x, campo) {
  if (is.character(x) || is.na(x)) return(.ee_poder_fmt(x))
  if (grepl("df|sample size", campo, ignore.case = TRUE)) return(.ee_poder_fmt(x))
  formatC(x, format = "f", digits = 6)
}

# Menor número de casas (de 2 a 4) com que o efeito arredondado devolve o mesmo N.
# Evita o problema do curso em que f = 0,3333 vira 0,33 e o N muda de 129 para 132.
# Confere com o mesmo arredondamento da frase (.pr_num, meio para cima): round() do R
# levaria ,625 a ,62, a frase imprimiria ,63 e o N do valor impresso seria outro
# (verificação final, G-12).
.ee_poder_casas <- function(valor, calc_n, n_ref, de = 2) {
  for (k in de:4) {
    n_k <- tryCatch(calc_n(as.numeric(.pr_num(valor, k, TRUE, "."))), error = function(e) NA)
    if (isTRUE(n_k == n_ref)) return(k)
  }
  4
}

# Ajuste por perda. Padrão: N / (1 - perda), que garante o N planejado depois da perda.
# O livro fala em inflar o N em 10% a 20%, o que, feito como
# N * (1 + perda), fica abaixo do planejado depois da perda; sai ao lado para comparação.
ajustar_perda <- function(n, perda, passo = 1) {
  if (perda < 0 || perda >= 1) stop("perda é uma proporção entre 0 e 1 (ex.: 0.20)")
  n_div <- ceiling(ceiling(n / (1 - perda) - 1e-9) / passo) * passo
  n_livro <- ceiling(ceiling(n * (1 + perda) - 1e-9) / passo) * passo
  data.frame(n_planejado = n, perda = perda, n_coleta = n_div,
             n_coleta_inflado = n_livro, restam_inflado = floor(n_livro * (1 - perda)))
}

# ---------------------------------------------------------------- conversões de efeito
# Fórmulas do G*Power (botão "Determine") e de conversão corrente; ver
# references/poder-amostra.md. Devolvem número simples para usar direto nas funções.
efeito_f_de_eta2 <- function(eta2) {
  if (any(eta2 < 0 | eta2 >= 1)) stop("eta2 precisa estar em [0, 1)")
  sqrt(eta2 / (1 - eta2))
}
efeito_eta2_de_f <- function(f) f^2 / (1 + f^2)

# d pelo "Determine" do G*Power: DP = raiz da média das variâncias. Com n1 e n2, DP
# combinado ponderado pelos graus de liberdade (a fonte do curso não diz qual usou).
efeito_d_de_medias <- function(m1, m2, dp1, dp2 = dp1, n1 = NULL, n2 = NULL) {
  dp <- if (is.null(n1) || is.null(n2)) sqrt((dp1^2 + dp2^2) / 2) else
    sqrt(((n1 - 1) * dp1^2 + (n2 - 1) * dp2^2) / (n1 + n2 - 2))
  abs(m2 - m1) / dp
}

# f pelo "Determine" da ANOVA de um fator: DP das médias (ponderado por n quando dado)
# sobre o DP dentro dos grupos.
efeito_f_de_medias <- function(medias, dp, n = NULL) {
  if (is.null(n)) n <- rep(1, length(medias))
  if (length(n) != length(medias)) stop("n precisa ter um valor por média")
  mw <- sum(n * medias) / sum(n)
  sqrt(sum(n * (medias - mw)^2) / sum(n)) / dp
}

efeito_d_de_r <- function(r) 2 * r / sqrt(1 - r^2)
efeito_r_de_d <- function(d, n1 = NULL, n2 = NULL) {
  a <- if (is.null(n1) || is.null(n2)) 4 else (n1 + n2)^2 / (n1 * n2)
  d / sqrt(d^2 + a)
}
efeito_d_de_t <- function(t, n1, n2) t * sqrt(1 / n1 + 1 / n2)
efeito_f2_de_r2 <- function(r2) {
  if (any(r2 < 0 | r2 >= 1)) stop("R² precisa estar em [0, 1)")
  r2 / (1 - r2)
}
# f² do incremento: ΔR² / (1 - R² do modelo completo), com R² completo = R²0 + ΔR².
efeito_f2_de_delta_r2 <- function(r2_reduzido, delta_r2) {
  r2c <- r2_reduzido + delta_r2
  if (any(r2c >= 1)) stop("R² reduzido + ΔR² precisa ficar abaixo de 1")
  delta_r2 / (1 - r2c)
}
efeito_or_de_probabilidades <- function(p1, p0) (p1 / (1 - p1)) / (p0 / (1 - p0))
# dz (t pareado do G*Power) a partir de d entre momentos e da correlação entre eles.
efeito_dz_de_d <- function(d, r) {
  if (any(r <= -1 | r >= 1)) stop("r precisa estar entre -1 e 1")
  d / sqrt(2 * (1 - r))
}

# ---------------------------------------------------------------- montagem do resultado
.ee_poder_tipo <- function(n) if (is.null(n)) "A priori: Compute required sample size" else
  "Post hoc: Compute achieved power"

# α (e poder) com casas suficientes para ser o mesmo valor usado na conta: .05, .001,
# .005; poder .80, .975, .999; valor sem representação curta (.05/3) sai com 3
# algarismos significativos (.0167).
.ee_poder_alfa_txt <- function(alfa, decimal = ".") {
  k <- which(vapply(2:6, function(k) abs(round(alfa, k) - alfa) < 1e-10, TRUE))
  casas <- if (length(k)) (2:6)[k[1]] else max(2, ceiling(-log10(alfa)) + 2)
  .pr_num(alfa, casas, FALSE, decimal)
}

.ee_poder_caudas_txt <- function(caudas, idioma = "pt") {
  if (idioma == "en") (if (caudas == 2) "two-tailed" else "one-tailed") else (if (caudas == 2) "bicaudal" else "unicaudal")
}

.ee_poder_apa <- function(analise_pt, analise_en, efeito_txt, alfa, caudas, poder, n_final, detalhe,
                          perda, n_perda, justificativa, decimal, idioma, abertura = NULL) {
  sep <- .ee_sep(decimal)
  just <- if (is.null(justificativa)) {
    if (idioma == "en") "[source of the effect size: fill in]" else "[origem do tamanho de efeito: preencher]"
  } else justificativa
  a_txt <- paste0("α = ", .ee_poder_alfa_txt(alfa, decimal))
  if (!is.null(caudas)) a_txt <- paste0(a_txt, " (", .ee_poder_caudas_txt(caudas, idioma), ")")
  # probabilidade: sem zero à esquerda, com as
  # casas do valor usado na conta, como o α (.975 não vira .98, nem .999 vira 1.00)
  p_txt <- .ee_poder_alfa_txt(poder, decimal)
  if (idioma == "en") {
    s <- paste0(if (is.null(abertura)) "Sample size was determined a priori with G*Power 3.1 settings (Faul et al., 2007) " else abertura[["en"]], analise_en,
                ", assuming ", efeito_txt, " (", just, "), ", a_txt, " and power of ", p_txt,
                ", which required a minimum of ", n_final, " participants", detalhe, ".")
    if (perda > 0) s <- paste0(s, " Anticipating ", round(100 * perda), "% attrition, the recruitment target is ",
                               n_perda, " participants.")
  } else {
    s <- paste0(if (is.null(abertura)) "O tamanho amostral foi determinado a priori com os parâmetros do G*Power 3.1 (Faul et al., 2007) " else abertura[["pt"]],
                analise_pt, ", com ", efeito_txt, " (", just, "), ", a_txt, " e poder de ", p_txt,
                ", resultando em um mínimo de ", n_final, " participantes", detalhe, ".")
    if (perda > 0) s <- paste0(s, " Prevendo perda de ", round(100 * perda), "%, a meta de coleta é de ",
                               n_perda, " participantes.")
  }
  s
}

# Achata listas e vetores nomeados em vetor de texto (campo -> valor). Números com até
# 4 casas; NULL some (permite `if (...) list(...)` dentro das listas de campos).
.ee_poder_achata <- function(x, saida = FALSE) {
  out <- character(); nms <- character()
  fmt <- function(v, nm) if (saida) .ee_poder_fmt_saida(v, nm) else .ee_poder_fmt(v)
  rec <- function(obj, nm = NULL) {
    if (is.null(obj)) return(invisible())
    if (is.list(obj)) { for (i in seq_along(obj)) rec(obj[[i]], names(obj)[i]); return(invisible()) }
    if (length(obj) == 1 && !is.null(nm) && nzchar(nm)) {
      if (!is.na(obj)) { out <<- c(out, fmt(obj, nm)); nms <<- c(nms, nm) }
      return(invisible())
    }
    for (j in seq_along(obj)) if (!is.na(obj[[j]])) { out <<- c(out, fmt(obj[[j]], names(obj)[j])); nms <<- c(nms, names(obj)[j]) }
  }
  rec(x)
  stats::setNames(out, nms)
}

# entradas, saidas e skill: listas nomeadas (nome = campo do G*Power; skill = linhas que
# a skill acrescenta, como N mínimo, N com grupos iguais e perda).
.ee_poder_resultado <- function(analise, familia, teste, n, entradas, saidas, skill, n_res,
                                efeito, decisoes, avisos, apa, extra = list()) {
  bloco <- function(secao, v) {
    v <- .ee_poder_achata(v, saida = secao == "G*Power: saída")
    if (!length(v)) return(NULL)
    data.frame(secao = secao, campo = names(v), valor = unname(v), stringsAsFactors = FALSE)
  }
  tab <- rbind(bloco("G*Power: entrada", list(`Test family` = familia, `Statistical test` = teste,
                                              `Type of power analysis` = .ee_poder_tipo(n), entradas)),
               bloco("G*Power: saída", saidas), bloco("skill", skill))
  rownames(tab) <- NULL
  if (!is.null(n)) {
    avisos <- c(avisos, paste0("Poder post hoc: se o tamanho de efeito informado veio da própria amostra, o poder ",
                               "é função direta do p e não acrescenta informação (livro). Usar ",
                               "efeito definido antes da coleta (menor efeito de interesse, meta-análise) ou ",
                               "análise de sensibilidade."))
  }
  res <- .ee_resultado(analise, n = n_res, tabela = tab, efeito = efeito, decisoes = decisoes,
                       avisos = avisos, apa = apa, extra = extra)
  class(res) <- c("pr_poder", class(res))
  res
}

# Impressão em lista (campo e valor), que cabe na largura do console; o resto segue o
# formato de print.pr_resultado. O objeto continua sendo um pr_resultado.
print.pr_poder <- function(x, ...) {
  cat("== ", x$analise, " (N = ", x$n, ")\n", sep = "")
  tab <- x$tabela
  for (s in unique(tab$secao)) {
    t <- tab[tab$secao == s, ]
    cat("\n", c(`G*Power: entrada` = "G*Power, entrada:", `G*Power: saída` = "G*Power, saída:",
                skill = "Skill:")[[s]], "\n", sep = "")
    cat(paste0("  ", format(t$campo), "  ", t$valor), sep = "\n")
  }
  .ee_print_df(x$efeito, "Tamanho de efeito:")
  if (length(x$decisoes)) cat("\nDecisões:\n", paste0("  - ", x$decisoes, collapse = "\n"), "\n", sep = "")
  if (length(x$avisos)) cat("\nAvisos:\n", paste0("  ! ", x$avisos, collapse = "\n"), "\n", sep = "")
  if (length(x$apa)) cat("\nAPA: ", paste(x$apa, collapse = "\n     "), "\n", sep = "")
  invisible(x)
}

.ee_poder_efeito <- function(indices, valores, principal) {
  data.frame(indice = indices, valor = valores, ic_inf = NA_real_, ic_sup = NA_real_,
             principal = indices == principal, stringsAsFactors = FALSE)
}

.ee_poder_perda <- function(n_final, perda, passo, decisoes, avisos) {
  if (perda <= 0) return(list(n_perda = NA, skill = NULL, avisos = avisos))
  p <- ajustar_perda(n_final, perda, passo)
  avisos <- c(avisos, paste0(
    "Ajuste por perda: N / (1 - perda) = ", p$n_coleta, " é acréscimo da skill, não vem das fontes do usuário ",
    "(o curso não ensina ajuste por perda). O livro sugere inflar o N em 10% a 20% como decisão do pesquisador ",
    "; feito como N × (1 + perda) daria ", p$n_coleta_inflado, ", e depois da perda restariam ",
    p$restam_inflado, if (p$restam_inflado < n_final) paste0(", abaixo dos ", n_final, " planejados") else "", "."))
  list(n_perda = p$n_coleta,
       skill = list(`N para coleta (N / (1 - perda))` = p$n_coleta), avisos = avisos)
}

# ---------------------------------------------------------------- correlação
# G*Power: Exact > Correlation: Bivariate normal model (distribuição exata de r).
# metodo = "fisher": aproximação z de Fisher, a conta do texto do livro (;
#): N = ((z[1-α/2] + z[poder]) / atanh(r))² + 3. As duas saem sempre.
.ee_poder_r_fisher <- function(r, n, alfa, caudas) {
  z <- atanh(abs(r)) * sqrt(n - 3)
  zc <- stats::qnorm(1 - alfa / caudas)
  if (caudas == 2) stats::pnorm(z - zc) + stats::pnorm(-z - zc) else stats::pnorm(z - zc)
}
.ee_poder_n_fisher <- function(r, poder, alfa, caudas) {
  ceiling(((stats::qnorm(1 - alfa / caudas) + stats::qnorm(poder)) / atanh(abs(r)))^2 + 3 - 1e-9)
}
.ee_poder_r_exato <- function(r, n = NULL, poder = NULL, alfa, caudas, n_max = 1e4) {
  .pr_exige("pwrss")
  alt <- if (caudas == 2) "two.sided" else "one.sided"
  x <- pwrss::power.exact.onecor(rho = abs(r), n = n, power = poder, alpha = alfa, alternative = alt,
                                 n.max = n_max, verbose = 0)
  list(n = x$n, poder = x$power)
}

poder_correlacao <- function(r, n = NULL, poder = 0.80, alfa = 0.05, caudas = 2,
                             metodo = c("exato", "fisher"), perda = 0, justificativa = NULL,
                             decimal = ".", idioma = "pt") {
  metodo <- match.arg(metodo)
  .ee_poder_checa(poder, alfa, perda, caudas, n)
  if (abs(r) <= 0 || abs(r) >= 1) stop("r precisa estar entre 0 e 1 em módulo")
  decisoes <- character(); avisos <- character()
  n_fisher <- .ee_poder_n_fisher(r, if (is.null(n)) poder else .8, alfa, caudas)
  n_max <- max(1e4, 2 * n_fisher)
  if (is.null(n)) {
    ex <- .ee_poder_r_exato(r, poder = poder, alfa = alfa, caudas = caudas, n_max = n_max)
    nf <- .ee_poder_n_fisher(r, poder, alfa, caudas)
    n_gp <- ex$n; p_gp <- ex$poder; p_fi <- .ee_poder_r_fisher(r, nf, alfa, caudas)
    n_res <- if (metodo == "exato") ex$n else nf
    decisoes <- c(decisoes, paste0("N pelo método exato (o do G*Power) = ", ex$n, "; pela aproximação z de Fisher = ", nf,
                                   ". Método do resultado ($n e frase): ", metodo, "."))
  } else {
    n_gp <- n; nf <- NA
    p_gp <- .ee_poder_r_exato(r, n = n, alfa = alfa, caudas = caudas, n_max = n_max)$poder
    p_fi <- .ee_poder_r_fisher(r, n, alfa, caudas)
    n_res <- n
    decisoes <- c(decisoes, paste0("Poder com N = ", n, ": exato (G*Power) = ", round(p_gp, 3), "; Fisher = ", round(p_fi, 3), "."))
  }
  decisoes <- c(decisoes, paste0(
    "O livro faz a conta pela aproximação de Fisher no texto (62 para r = 0,35) e mostra o G*Power exato na tela ",
    "(61); o padrão aqui é o exato, que é o que o G*Power devolve. pwr::pwr.r.test usa ",
    "outra correção e dá 85 onde o G*Power dá 84 (r = 0,30)."))
  if (metodo == "fisher") decisoes <- c(decisoes, paste0(
    "Com metodo = 'fisher' o bloco 'G*Power, saída' continua mostrando o exato (", n_gp, "), que é o que a tela devolve; ",
    "o N de Fisher fica no bloco da skill e a frase não cita o G*Power."))
  crit_t <- stats::qt(1 - alfa / caudas, n_gp - 2)
  r_crit <- crit_t / sqrt(crit_t^2 + n_gp - 2)
  pd <- .ee_poder_perda(n_res, perda, 1, decisoes, avisos); avisos <- pd$avisos
  calc <- function(rr) if (metodo == "exato") .ee_poder_r_exato(rr, poder = poder, alfa = alfa, caudas = caudas, n_max = n_max)$n else
    .ee_poder_n_fisher(rr, poder, alfa, caudas)
  k <- if (is.null(n)) .ee_poder_casas(abs(r), calc, n_res) else 2
  ef_txt <- paste0("ρ = ", .pr_num(abs(r), k, FALSE, decimal))
  abertura <- if (metodo == "fisher") c(pt = "O tamanho amostral foi determinado a priori pela aproximação z de Fisher ",
                                        en = "Sample size was determined a priori with Fisher's z approximation ") else NULL
  apa <- if (is.null(n)) .ee_poder_apa("para detectar uma correlação de Pearson", "to detect a Pearson correlation",
                                       ef_txt, alfa, caudas, poder, n_res, "", perda, pd$n_perda, justificativa,
                                       decimal, idioma, abertura = abertura) else character()
  .ee_poder_resultado(
    "Poder e tamanho amostral: correlação de Pearson", "Exact", "Correlation: Bivariate normal model", n,
    entradas = list(`Tail(s)` = if (caudas == 2) "Two" else "One", `Correlation ρ H1` = abs(r), `Correlation ρ H0` = 0,
                 `α err prob` = .ee_poder_fmt(alfa, 6), if (is.null(n)) list(`Power (1-β err prob)` = poder) else list(`Total sample size` = n)),
    saidas = list(`Lower critical r` = -r_crit, `Upper critical r` = r_crit,
               if (is.null(n)) list(`Total sample size` = n_gp), `Actual power` = p_gp),
    skill = list(`Método do resultado` = metodo, `N mínimo do resultado` = if (is.null(n)) n_res else NA,
                 `N pela aproximação de Fisher` = nf, `Poder pela aproximação de Fisher` = p_fi, pd$skill),
    n_res = n_res, efeito = .ee_poder_efeito("rho", abs(r), "rho"), decisoes = decisoes, avisos = avisos, apa = apa,
    extra = list(n_exato = if (is.null(n)) n_gp else NA, n_fisher = nf))
}

# ---------------------------------------------------------------- teste t
# G*Power: t tests > Means: Difference between two independent means (two groups),
# Difference between two dependent means (matched pairs) ou Difference from constant
# (one sample case). δ = d·√(n1n2/(n1+n2)) ou dz·√N; poder bilateral soma as duas caudas
# (igual a pwr::pwr.t.test e pwr.t2n.test).
# Alocação N2/N1 diferente de 1, regra do G*Power: resolve n1* contínuo com n2 = razão ×
# n1 e devolve n1 = teto(n1*) e n2 = teto(razão × n1*), cada grupo arredondado por conta
# própria. Conferida no manual (2023) em dois exemplos: seção 15, duas variâncias com
# razão 2, 148 + 295 (n1* = 147,44); seção 17, diferença entre inclinações (teste t) com
# razão 1,571428, 163 + 256 (n1* = 162,50). Buscar o menor n1 inteiro com n2 = teto(razão
# × n1) daria 148 + 296 e 163 + 257. `pw_n1(n1)` devolve o poder com n2 = razão × n1.
# Se o menor desenho (n1 = lo) já alcança o poder pedido (efeito muito grande), n1* = lo:
# o uniroot pararia com "f() values at end points not of opposite sign" (verificação
# final, F1; a busca inteira anterior devolvia esse mesmo desenho mínimo).
.ee_poder_alocacao <- function(pw_n1, razao, poder, lo = 2) {
  minimo <- isTRUE(pw_n1(lo) >= poder)
  n1c <- if (minimo) lo else stats::uniroot(function(x) pw_n1(x) - poder, c(lo + 1e-6, 1e7), tol = 1e-10)$root
  list(n1c = n1c, n1 = ceiling(n1c - 1e-9), n2 = ceiling(razao * n1c - 1e-9), minimo = minimo)
}
.ee_poder_t_razao <- function(d, razao, poder, alfa, caudas) {
  pw <- function(n1) { n2 <- razao * n1; .ee_poder_t(d * sqrt(n1 * n2 / (n1 + n2)), n1 + n2 - 2, alfa, caudas)$poder }
  .ee_poder_alocacao(pw, razao, poder, lo = max(2, 2 / razao))$n1c
}

poder_teste_t <- function(d, desenho = c("independente", "pareado", "uma_amostra"), n = NULL,
                          poder = 0.80, alfa = 0.05, caudas = 2, razao = 1, r_pareado = NULL,
                          perda = 0, justificativa = NULL, decimal = ".", idioma = "pt") {
  desenho <- match.arg(desenho)
  .ee_poder_checa(poder, alfa, perda, caudas, n)
  if (d == 0) stop("d precisa ser diferente de zero")
  decisoes <- character(); avisos <- character()
  d_in <- abs(d)
  # dz só no pareado; no de uma amostra o G*Power chama de "Effect size d" (manual, seção 20)
  ef_nome <- if (desenho == "pareado") "dz" else "d"
  if (desenho == "pareado") {
    if (!is.null(r_pareado)) {
      d_in <- efeito_dz_de_d(abs(d), r_pareado)
      decisoes <- c(decisoes, paste0("d entre momentos = ", round(abs(d), 4), " convertido em dz = d / √(2(1 - r)) = ",
                                     round(d_in, 4), " com r = ", r_pareado, "."))
    } else {
      avisos <- c(avisos, paste0("O G*Power usa dz (diferença média / DP das diferenças). Um d, g ou Δ de literatura ",
                                 "calculado com o DP dos escores só é igual a dz quando a correlação entre as medidas é ",
                                 "0,50; se for o caso, informar r_pareado para converter (curso, seção 14)."))
    }
  }
  if (desenho == "independente") {
    pw <- function(n1, n2) .ee_poder_t(d_in * sqrt(n1 * n2 / (n1 + n2)), n1 + n2 - 2, alfa, caudas)
    if (is.null(n)) {
      if (razao == 1) {
        n1 <- .ee_poder_busca(function(k) pw(k, k)$poder, poder, 2); n2 <- n1
      } else {
        n1c <- .ee_poder_t_razao(d_in, razao, poder, alfa, caudas)
        n1 <- ceiling(n1c - 1e-9); n2 <- ceiling(razao * n1c - 1e-9)
        if (n1c <= max(2, 2 / razao) + 1e-9) avisos <- c(avisos, paste0(
          "O menor desenho com essa razão (", n1, " + ", n2, ", pelo menos 2 casos em cada grupo) já passa do poder pedido: ",
          "o efeito é tão grande que o N vem do mínimo para estimar a variância, não do poder. Conferir se o d é plausível."))
        decisoes <- c(decisoes, paste0("Razão N2/N1 = ", razao, ": n1* contínuo = ", round(n1c, 2), "; n1 = teto(n1*) = ", n1,
                                       " e n2 = teto(razão × n1*) = ", n2, ", cada grupo arredondado para cima, como o G*Power ",
                                       "(regra conferida no manual: diferença entre inclinações, teste t com razão 1,571428, 163 + 256; ",
                                       "duas variâncias com razão 2, 148 + 295). Nenhum exemplo do manual é do t de duas médias com ",
                                       "razão diferente de 1: conferir na tela antes de relatar."))
      }
    } else {
      if (length(n) == 2) { n1 <- n[1]; n2 <- n[2] } else { n1 <- round(n / (1 + razao)); n2 <- n - n1 }
    }
    r <- pw(n1, n2); N <- n1 + n2
    delta <- d_in * sqrt(n1 * n2 / (n1 + n2)); gl <- N - 2
    teste <- "Means: Difference between two independent means (two groups)"
    ent <- list(`Effect size d` = d_in, `Allocation ratio N2/N1` = razao)
    calc <- function(dd) {
      if (razao == 1) return(2 * .ee_poder_busca(function(k) .ee_poder_t(dd * sqrt(k / 2), 2 * k - 2, alfa, caudas)$poder, poder, 2))
      a <- .ee_poder_t_razao(dd, razao, poder, alfa, caudas)
      ceiling(a - 1e-9) + ceiling(razao * a - 1e-9)
    }
    sai_n <- list(`Sample size group 1` = n1, `Sample size group 2` = n2, `Total sample size` = N)
    detalhe <- if (n1 == n2) paste0(" (", n1, if (idioma == "en") " per group)" else " por grupo)") else
      paste0(" (", n1, " + ", n2, ")")
    an_pt <- "para comparar duas médias independentes (teste t)"; an_en <- "to compare two independent means (t test)"
  } else {
    pw <- function(k) .ee_poder_t(d_in * sqrt(k), k - 1, alfa, caudas)
    N <- if (is.null(n)) .ee_poder_busca(function(k) pw(k)$poder, poder, 2) else n
    r <- pw(N); delta <- d_in * sqrt(N); gl <- N - 1
    teste <- if (desenho == "pareado") "Means: Difference between two dependent means (matched pairs)" else
      "Means: Difference from constant (one sample case)"
    ent <- if (desenho == "pareado") list(`Effect size dz` = d_in) else list(`Effect size d` = d_in)
    calc <- function(dd) .ee_poder_busca(function(k) .ee_poder_t(dd * sqrt(k), k - 1, alfa, caudas)$poder, poder, 2)
    sai_n <- list(`Total sample size` = N)
    detalhe <- if (desenho == "pareado") (if (idioma == "en") " (pairs)" else " (pares de medidas)") else ""
    an_pt <- if (desenho == "pareado") "para comparar duas medidas pareadas (teste t pareado)" else "para o teste t de uma amostra"
    an_en <- if (desenho == "pareado") "to compare two paired measurements (paired t test)" else "for a one-sample t test"
  }
  pd <- .ee_poder_perda(N, perda, if (desenho == "independente" && razao == 1) 2 else 1, decisoes, avisos); avisos <- pd$avisos
  k <- if (is.null(n)) .ee_poder_casas(d_in, calc, N) else 2
  apa <- if (is.null(n)) .ee_poder_apa(an_pt, an_en, paste0(ef_nome, " = ", .pr_num(d_in, k, TRUE, decimal)),
                                       alfa, caudas, poder, N, detalhe, perda, pd$n_perda, justificativa, decimal, idioma) else character()
  ef <- .ee_poder_efeito(ef_nome, d_in, ef_nome)
  if (desenho == "pareado" && !is.null(r_pareado)) ef <- rbind(ef, .ee_poder_efeito("d_entre_momentos", abs(d), ef_nome))
  .ee_poder_resultado(
    paste0("Poder e tamanho amostral: teste t (", desenho, ")"), "t tests", teste, n,
    entradas = list(`Tail(s)` = if (caudas == 2) "Two" else "One", ent, `α err prob` = .ee_poder_fmt(alfa, 6),
                 if (is.null(n)) list(`Power (1-β err prob)` = poder) else sai_n),
    saidas = list(`Noncentrality parameter δ` = delta, `Critical t` = r$critico, Df = gl,
               if (is.null(n)) sai_n, `Actual power` = r$poder),
    skill = list(`N mínimo` = if (is.null(n)) N else NA, pd$skill),
    n_res = N, efeito = ef, decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- ANOVA entre sujeitos
# Um fator: ANOVA: Fixed effects, omnibus, one-way (λ = f²N; gl = k - 1, N - k); o
# G*Power arredonda N para grupos iguais (128, 159, 180, 200 para f = 0,25 e k = 2 a 5).
# Fatorial: ANOVA: Fixed effects, special, main effects and interactions (λ = f²N; gl1 =
# gl do efeito; gl2 = N - células); o G*Power não força múltiplo de células (244 no
# curso). ANCOVA: gl2 = N - células - covariáveis. Termos: letras A, B, C... na ordem de
# `niveis`; efeito = NULL usa a interação de ordem mais alta (o efeito mais exigente na
# regra do curso, seção 18; livro).
.ee_poder_termos <- function(niveis) {
  letras <- LETTERS[seq_along(niveis)]
  termos <- unlist(lapply(seq_along(niveis), function(m) apply(utils::combn(letras, m), 2, paste, collapse = ":")))
  gl <- vapply(termos, function(t) prod(niveis[match(strsplit(t, ":")[[1]], letras)] - 1), numeric(1))
  stats::setNames(gl, termos)
}

poder_anova <- function(f = NULL, eta2p = NULL, niveis, efeito = NULL, gl_efeito = NULL, covariaveis = 0,
                        n = NULL, poder = 0.80, alfa = 0.05, perda = 0, justificativa = NULL,
                        decimal = ".", idioma = "pt") {
  .ee_poder_checa(poder, alfa, perda, 2, n)
  if (is.null(f) == is.null(eta2p)) stop("informe f ou eta2p (um dos dois)")
  if (is.null(f)) f <- efeito_f_de_eta2(eta2p)
  if (any(niveis < 2)) stop("cada fator precisa de pelo menos 2 níveis")
  decisoes <- character(); avisos <- character()
  celulas <- prod(niveis)
  termos <- .ee_poder_termos(niveis)
  if (!is.null(gl_efeito)) {
    gl1 <- gl_efeito; rotulo <- paste0("efeito com gl = ", gl_efeito)
  } else {
    if (is.null(efeito)) efeito <- names(termos)[length(termos)]
    efeito <- toupper(efeito)
    if (!efeito %in% names(termos)) stop("efeito precisa ser um de: ", paste(names(termos), collapse = ", "))
    gl1 <- termos[[efeito]]; rotulo <- efeito
  }
  um_fator <- length(niveis) == 1 && covariaveis == 0 && is.null(gl_efeito)
  fp <- function(N, ff = f) .ee_poder_f(ff^2 * N, gl1, N - celulas - covariaveis, alfa)
  inicio <- celulas + covariaveis + 1
  if (is.null(n)) {
    n_min <- .ee_poder_busca(function(N) fp(N)$poder, poder, inicio)
    n_iguais <- ceiling(n_min / celulas) * celulas
    N <- if (um_fator) n_iguais else n_min
  } else {
    n_min <- NA; n_iguais <- NA; N <- n
  }
  r <- fp(N)
  if (um_fator) {
    teste <- "ANOVA: Fixed effects, omnibus, one-way"
    ent <- list(`Effect size f` = f, `α err prob` = .ee_poder_fmt(alfa, 6), if (is.null(n)) list(`Power (1-β err prob)` = poder),
             `Number of groups` = celulas)
    if (is.null(n)) decisoes <- c(decisoes, paste0("N mínimo inteiro = ", n_min, "; o G*Power arredonda para grupos iguais: ",
                                                   N, " (", N / celulas, " por grupo)."))
  } else {
    teste <- if (covariaveis > 0) "ANCOVA: Fixed effects, main effects and interactions" else
      "ANOVA: Fixed effects, special, main effects and interactions"
    ent <- list(`Effect size f` = f, `α err prob` = .ee_poder_fmt(alfa, 6), if (is.null(n)) list(`Power (1-β err prob)` = poder),
             `Numerator df` = gl1, `Number of groups` = celulas, if (covariaveis > 0) list(`Number of covariates` = covariaveis))
    if (is.null(n)) decisoes <- c(decisoes, paste0("O G*Power devolve o N mínimo (", n_min, ") sem forçar múltiplo do número de ",
                                                   "células; para células iguais seriam ", n_iguais, " (", n_iguais / celulas, " por célula)."))
    if (length(niveis) > 1 && is.null(gl_efeito) && is.null(n)) {
      nt <- vapply(names(termos), function(t) .ee_poder_busca(function(N) .ee_poder_f(f^2 * N, termos[[t]], N - celulas - covariaveis, alfa)$poder,
                                                              poder, inicio), numeric(1))
      decisoes <- c(decisoes, paste0("Com o mesmo f, N por termo: ", paste(paste0(names(nt), " = ", nt), collapse = "; "),
                                     ". O mais exigente é ", names(nt)[which.max(nt)], " (", max(nt), "); a regra do curso ",
                                     "(seção 18) e do livro é planejar pelo efeito que responde à hipótese principal, ",
                                     "em geral a interação."))
    }
    if (covariaveis > 0) avisos <- c(avisos, paste0(
      "ANCOVA: o N só fica igual ao da ANOVA se o ηp² assumido já for o ajustado pela covariável. Se o efeito veio de ",
      "estudo sem a covariável, o f ajustado é f / √(1 - r²), com r a correlação covariável-desfecho (curso, seção 20)."))
  }
  if (!is.null(eta2p)) decisoes <- c(decisoes, paste0("ηp² = ", eta2p, " convertido em f = √(ηp² / (1 - ηp²)) = ", round(f, 4),
                                                      ". O ηp² precisa ser do mesmo desenho (mesmos fatores e covariáveis) para a conversão valer."))
  pd <- .ee_poder_perda(N, perda, if (um_fator) celulas else 1, decisoes, avisos); avisos <- pd$avisos
  calc <- function(ff) { m <- .ee_poder_busca(function(NN) fp(NN, ff)$poder, poder, inicio); if (um_fator) ceiling(m / celulas) * celulas else m }
  k <- if (is.null(n)) .ee_poder_casas(f, calc, N) else 2
  ef_txt <- if (!is.null(eta2p)) paste0("ηp² = ", .pr_num(eta2p, 3, FALSE, decimal), .ee_sep(decimal), "f = ", .pr_num(f, k, TRUE, decimal)) else
    paste0("f = ", .pr_num(f, k, TRUE, decimal))
  desenho <- if (length(niveis) == 1) paste0("de um fator com ", niveis, " grupos") else paste(niveis, collapse = " × ")
  desenho_en <- if (length(niveis) == 1) paste0("one-way (", niveis, " groups)") else paste(niveis, collapse = " × ")
  liga <- if (length(niveis) == 1) " e " else " com "
  cov_pt <- if (covariaveis == 1) paste0(liga, "1 covariável") else paste0(liga, covariaveis, " covariáveis")
  cov_en <- if (covariaveis == 1) " with 1 covariate" else paste0(" with ", covariaveis, " covariates")
  termo_pt <- if (!is.null(gl_efeito)) paste0("o efeito com ", gl_efeito, " gl no numerador") else
    if (length(niveis) == 1) "o efeito do fator" else
    if (grepl(":", rotulo)) paste0("a interação ", gsub(":", " × ", rotulo)) else paste0("o efeito principal de ", rotulo)
  termo_en <- if (!is.null(gl_efeito)) paste0("the effect with ", gl_efeito, " numerator df") else
    if (length(niveis) == 1) "the factor effect" else
    if (grepl(":", rotulo)) paste0("the ", gsub(":", " × ", rotulo), " interaction") else paste0("the main effect of ", rotulo)
  an_pt <- if (um_fator) paste0("para uma ANOVA de um fator com ", celulas, " grupos") else
    paste0("para ", termo_pt, " de uma ", if (covariaveis > 0) "ANCOVA" else "ANOVA fatorial", " ", desenho,
           if (covariaveis > 0) cov_pt else "")
  an_en <- if (um_fator) paste0("for a one-way ANOVA with ", celulas, " groups") else
    paste0("for ", termo_en, " of a ", desenho_en, if (covariaveis > 0) " ANCOVA" else " factorial ANOVA",
           if (covariaveis > 0) cov_en else "")
  detalhe <- if (um_fator) paste0(" (", N / celulas, if (idioma == "en") " per group)" else " por grupo)") else ""
  apa <- if (is.null(n)) .ee_poder_apa(an_pt, an_en, ef_txt, alfa, NULL, poder, N, detalhe, perda, pd$n_perda,
                                       justificativa, decimal, idioma) else character()
  ef <- .ee_poder_efeito(c("f", "eta2p"), c(f, efeito_eta2_de_f(f)), "f")
  .ee_poder_resultado(
    paste0("Poder e tamanho amostral: ", if (um_fator) "ANOVA de um fator" else if (covariaveis > 0) "ANCOVA" else "ANOVA fatorial"),
    "F tests", teste, n, entradas = c(ent, if (!is.null(n)) list(`Total sample size` = n)),
    saidas = list(`Noncentrality parameter λ` = f^2 * N, `Critical F` = r$critico, `Numerator df` = gl1,
               `Denominator df` = N - celulas - covariaveis, if (is.null(n)) list(`Total sample size` = N), `Actual power` = r$poder),
    skill = list(`N mínimo inteiro` = n_min, `N com células iguais` = n_iguais,
              `Poder com o N mínimo` = if (is.null(n)) fp(n_min)$poder else NA, pd$skill),
    n_res = N, efeito = ef, decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- ANOVA de medidas repetidas
# Modo padrão do G*Power ("as in GPower 3.0"), Faul et al. (2007, pp. 13-14) e
# pwrss::power.f.mixed.anova: λ = f²·N·m·ε/(1 - ρ) para intra e interação, λ = f²·N·m/(1 +
# (m - 1)ρ) para entre; gl intra = (m - 1)ε e (N - k)(m - 1)ε; interação = (k - 1)(m - 1)ε e
# (N - k)(m - 1)ε; entre = k - 1 e N - k. rho = NA: o f informado já incorpora a correlação
# (rota rho.within = NA do pwrss), λ = f²·N·ε; é a rota que reproduz os exemplos do curso
# (132 e 309). A opção "as in Cohen (1988)" do G*Power não está implementada (fórmula não
# confirmada em fonte primária; ver aviso e references/poder-amostra.md). N em múltiplos
# do número de grupos, como o G*Power (309 no curso).
.ee_poder_mr_pars <- function(f, N, k, m, rho, epsilon, efeito) {
  ajuste <- if (is.na(rho)) 1 else if (efeito == "entre") m / (1 + (m - 1) * rho) else m / (1 - rho)
  eps <- if (efeito == "entre") 1 else epsilon
  gl1 <- switch(efeito, intra = (m - 1) * epsilon, entre = k - 1, interacao = (k - 1) * (m - 1) * epsilon)
  gl2 <- if (efeito == "entre") N - k else (N - k) * (m - 1) * epsilon
  list(lambda = f^2 * ajuste * N * eps, gl1 = gl1, gl2 = gl2)
}

poder_anova_mr <- function(f = NULL, eta2p = NULL, medidas, grupos = 1,
                           efeito = c("intra", "entre", "interacao"), rho = 0.5, epsilon = 1,
                           n = NULL, poder = 0.80, alfa = 0.05, perda = 0, justificativa = NULL,
                           decimal = ".", idioma = "pt") {
  efeito <- match.arg(efeito)
  .ee_poder_checa(poder, alfa, perda, 2, n)
  if (is.null(f) == is.null(eta2p)) stop("informe f ou eta2p (um dos dois)")
  if (is.null(f)) f <- efeito_f_de_eta2(eta2p)
  if (medidas < 2 && efeito != "entre") stop("efeito intra ou de interação pede pelo menos 2 medidas")
  if (efeito %in% c("entre", "interacao") && grupos < 2) stop("efeito entre ou de interação pede pelo menos 2 grupos")
  if (medidas > 1 && (epsilon < 1 / (medidas - 1) || epsilon > 1)) stop("epsilon precisa ficar entre 1/(medidas - 1) e 1")
  if (!is.na(rho) && (rho <= -1 || rho >= 1)) stop("rho precisa estar entre -1 e 1 (ou NA)")
  decisoes <- character(); avisos <- character()
  fp <- function(N, ff = f) { p <- .ee_poder_mr_pars(ff, N, grupos, medidas, rho, epsilon, efeito)
    c(.ee_poder_f(p$lambda, p$gl1, p$gl2, alfa), p) }
  inicio <- grupos + 1
  if (is.null(n)) {
    n_min <- .ee_poder_busca(function(N) fp(N)$poder, poder, inicio)
    N <- ceiling(n_min / grupos) * grupos
    decisoes <- c(decisoes, paste0("N mínimo inteiro = ", n_min, "; arredondado para múltiplo do número de grupos (",
                                   grupos, "), como o G*Power: ", N, "."))
  } else { n_min <- NA; N <- n }
  r <- fp(N)
  teste <- switch(efeito, intra = "ANOVA: Repeated measures, within factors",
                  entre = "ANOVA: Repeated measures, between factors",
                  interacao = "ANOVA: Repeated measures, within-between interaction")
  if (is.na(rho)) {
    opcao <- "sem equivalente confirmado (rho = NA; ver aviso)"
    avisos <- c(avisos, paste0(
      "rho = NA: λ = f²·N·ε, sem ajuste pela correlação entre medidas (a rota rho.within = NA do pwrss). Os dois ",
      "exemplos do curso feitos no G*Power com Options = 'as in Cohen (1988)' (132 e 309) saem iguais por esta rota, ",
      "mas a fórmula dessa opção não foi encontrada no manual do G*Power (2023) nem em Faul et al. (2007), e a ",
      "skill não a implementa como tal. Antes de escrever 'as in Cohen (1988)' no método, conferir o N na tela do G*Power."))
  } else {
    opcao <- "as in GPower 3.0 (padrão)"
    decisoes <- c(decisoes, paste0("Modo padrão do G*Power: o f é o do efeito sem a correlação, que entra pelo campo ",
                                   "'Corr among rep measures' (ρ = ", rho, "). Um ηp² publicado de ANOVA de medidas repetidas ",
                                   "(SPSS, JASP) já embute a correlação: com ele, usar rho = NA ou converter antes."))
  }
  if (medidas == 2 && efeito != "entre" && epsilon == 1)
    decisoes <- c(decisoes, "Com 2 medidas não há esfericidade a corrigir: ε = 1 e nada a relatar sobre ela.")
  if (medidas > 2 && efeito != "entre" && epsilon == 1)
    avisos <- c(avisos, paste0("ε = 1 assume esfericidade (escolha liberal, N menor). Com ε menor o N cresce ",
                               "(curso, seção 17: f = 0,33 e poder 0,90 dão 132 com ε = 1, 162 com 0,75 e 217 com 0,50)."))
  if (!is.null(eta2p)) decisoes <- c(decisoes, paste0("ηp² = ", eta2p, " convertido em f = ", round(f, 4), "."))
  pd <- .ee_poder_perda(N, perda, grupos, decisoes, avisos); avisos <- pd$avisos
  calc <- function(ff) ceiling(.ee_poder_busca(function(NN) fp(NN, ff)$poder, poder, inicio) / grupos) * grupos
  k <- if (is.null(n)) .ee_poder_casas(f, calc, N) else 2
  ef_txt <- if (!is.null(eta2p)) paste0("ηp² = ", .pr_num(eta2p, 3, FALSE, decimal), .ee_sep(decimal), "f = ", .pr_num(f, k, TRUE, decimal)) else
    paste0("f = ", .pr_num(f, k, TRUE, decimal))
  nome_ef <- c(intra = "intra-sujeitos", entre = "entre grupos", interacao = "de interação grupo × tempo")[[efeito]]
  nome_en <- c(intra = "within-subjects", entre = "between-subjects", interacao = "group × time interaction")[[efeito]]
  an_pt <- paste0("para o efeito ", nome_ef, " de uma ANOVA de medidas repetidas (", grupos, if (grupos == 1) " grupo, " else " grupos, ", medidas, " medidas",
                  if (!is.na(rho)) paste0(", ρ = ", .pr_num(rho, 2, FALSE, decimal)) else "",
                  if (efeito != "entre" && medidas > 2) paste0(", ε = ", .pr_num(epsilon, 2, FALSE, decimal)) else "", ")")
  an_en <- paste0("for the ", nome_en, " effect of a repeated-measures ANOVA (", grupos, if (grupos == 1) " group, " else " groups, ", medidas, " measurements",
                  if (!is.na(rho)) paste0(", ρ = ", .pr_num(rho, 2, FALSE, decimal)) else "",
                  if (efeito != "entre" && medidas > 2) paste0(", ε = ", .pr_num(epsilon, 2, FALSE, decimal)) else "", ")")
  detalhe <- if (grupos > 1) paste0(" (", N / grupos, if (idioma == "en") " per group)" else " por grupo)") else ""
  apa <- if (is.null(n)) .ee_poder_apa(an_pt, an_en, ef_txt, alfa, NULL, poder, N, detalhe, perda, pd$n_perda,
                                       justificativa, decimal, idioma) else character()
  ent <- list(`Effect size f` = f, `α err prob` = .ee_poder_fmt(alfa, 6), if (is.null(n)) list(`Power (1-β err prob)` = poder) else list(`Total sample size` = n),
           `Number of groups` = grupos, `Number of measurements` = medidas,
           `Corr among rep measures` = if (is.na(rho)) "não usado" else rho,
           if (efeito != "entre") list(`Nonsphericity correction ε` = epsilon), `Options` = opcao)
  .ee_poder_resultado(
    paste0("Poder e tamanho amostral: ANOVA de medidas repetidas (", efeito, ")"), "F tests", teste, n, entradas = ent,
    saidas = list(`Noncentrality parameter λ` = r$lambda, `Critical F` = r$critico, `Numerator df` = r$gl1,
               `Denominator df` = r$gl2, if (is.null(n)) list(`Total sample size` = N), `Actual power` = r$poder),
    skill = list(`N mínimo inteiro` = n_min, `Poder com o N mínimo` = if (is.null(n)) fp(n_min)$poder else NA, pd$skill),
    n_res = N, efeito = .ee_poder_efeito(c("f", "eta2p"), c(f, efeito_eta2_de_f(f)), "f"),
    decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- regressão linear
# global: Linear multiple regression: Fixed model, R² deviation from zero (λ = f²N;
#   gl = p, N - p - 1). incremento: R² increase (λ = f²N; gl = q, N - p - 1; pwr.f2.test
#   não serve aqui porque usa λ = f²(u + v + 1)). coeficiente: t tests > Linear multiple
#   regression: Fixed model, single regression coefficient (δ = √(f²N); gl = N - p - 1).
poder_regressao <- function(tipo = c("global", "incremento", "coeficiente"), preditores,
                            f2 = NULL, r2 = NULL, r2_reduzido = NULL, delta_r2 = NULL, r2_parcial = NULL,
                            testados = NULL, n = NULL, poder = 0.80, alfa = 0.05, caudas = 2,
                            perda = 0, justificativa = NULL, decimal = ".", idioma = "pt") {
  tipo <- match.arg(tipo)
  .ee_poder_checa(poder, alfa, perda, caudas, n)
  decisoes <- character(); avisos <- character()
  p <- preditores
  if (is.null(f2)) {
    f2 <- switch(tipo,
                 global = if (!is.null(r2)) efeito_f2_de_r2(r2) else stop("informe f2 ou r2"),
                 incremento = if (!is.null(delta_r2) && !is.null(r2_reduzido)) efeito_f2_de_delta_r2(r2_reduzido, delta_r2) else
                   stop("informe f2 ou r2_reduzido e delta_r2"),
                 coeficiente = if (!is.null(r2_parcial)) efeito_f2_de_r2(r2_parcial) else stop("informe f2 ou r2_parcial"))
    origem <- switch(tipo, global = paste0("R² = ", r2), incremento = paste0("R²0 = ", r2_reduzido, " e ΔR² = ", delta_r2),
                     coeficiente = paste0("R² parcial = ", r2_parcial))
    decisoes <- c(decisoes, paste0("f² = ", round(f2, 4), " calculado de ", origem, "."))
  }
  if (tipo == "incremento") {
    if (is.null(testados)) stop("informe testados (preditores do bloco testado)")
    if (testados > p) stop("testados não pode passar do total de preditores")
    decisoes <- c(decisoes, "No G*Power, 'Number of tested predictors' é o bloco novo e 'Total number of predictors' é o modelo completo; trocar um pelo outro muda o N (curso, seção 23).")
  }
  if (tipo == "coeficiente") avisos <- c(avisos, paste0(
    "R² parcial = β² (como na aula) só vale com preditores não correlacionados. Com preditores correlacionados, ",
    "usar o r parcial de estudo anterior (r²p = t² / (t² + gl_erro)) ou f² = sr² / (1 - R² total) (curso, seção 22)."))
  if (tipo == "coeficiente") {
    fp <- function(N, ff = f2) { r <- .ee_poder_t(sqrt(ff * N), N - p - 1, alfa, caudas); r$gl <- N - p - 1; r$nc <- sqrt(ff * N); r }
  } else {
    u <- if (tipo == "global") p else testados
    fp <- function(N, ff = f2) { r <- .ee_poder_f(ff * N, u, N - p - 1, alfa); r$gl <- N - p - 1; r$nc <- ff * N; r }
  }
  inicio <- p + 2
  N <- if (is.null(n)) .ee_poder_busca(function(NN) fp(NN)$poder, poder, inicio) else n
  r <- fp(N)
  if (tipo == "global" && is.null(n)) decisoes <- c(decisoes, paste0(
    "Regra de bolso 50 + 8k (citada pelo curso, seção 4) daria ", 50 + 8 * p, "; serve só de checagem, não substitui o cálculo."))
  calc <- function(ff) .ee_poder_busca(function(NN) fp(NN, ff)$poder, poder, inicio)
  k <- if (is.null(n)) .ee_poder_casas(f2, calc, N, de = 3) else 3  # f² com 3 casas no mínimo (est_comum.R)
  ef_txt <- switch(tipo,
                   global = paste0(if (!is.null(r2)) paste0("R² = ", .pr_num(r2, 3, FALSE, decimal), .ee_sep(decimal)) else "",
                                   "f² = ", .pr_num(f2, k, TRUE, decimal)),
                   incremento = paste0(if (!is.null(delta_r2)) paste0("ΔR² = ", .pr_num(delta_r2, 3, FALSE, decimal),
                                                                      if (idioma == "en") " over R² = " else " sobre R² = ",
                                                                      .pr_num(r2_reduzido, 3, FALSE, decimal), .ee_sep(decimal)) else "",
                                       "f² = ", .pr_num(f2, k, TRUE, decimal)),
                   coeficiente = paste0(if (!is.null(r2_parcial)) paste0("R² parcial = ", .pr_num(r2_parcial, 3, FALSE, decimal), .ee_sep(decimal)) else "",
                                        "f² = ", .pr_num(f2, k, TRUE, decimal)))
  an_pt <- switch(tipo, global = paste0("para o R² de uma regressão linear múltipla com ", p, " preditores"),
                  incremento = paste0("para o incremento no R² de ", testados, " preditores num modelo com ", p, " preditores"),
                  coeficiente = paste0("para um coeficiente de regressão num modelo com ", p, " preditores"))
  an_en <- switch(tipo, global = paste0("for the R² of a multiple linear regression with ", p, " predictors"),
                  incremento = paste0("for the R² increase of ", testados, " predictors in a model with ", p, " predictors"),
                  coeficiente = paste0("for a single regression coefficient in a model with ", p, " predictors"))
  pd <- .ee_poder_perda(N, perda, 1, decisoes, avisos); avisos <- pd$avisos
  apa <- if (is.null(n)) .ee_poder_apa(an_pt, an_en, ef_txt, alfa, if (tipo == "coeficiente") caudas else NULL, poder, N, "",
                                       perda, pd$n_perda, justificativa, decimal, idioma) else character()
  teste <- switch(tipo, global = "Linear multiple regression: Fixed model, R² deviation from zero",
                  incremento = "Linear multiple regression: Fixed model, R² increase",
                  coeficiente = "Linear multiple regression: Fixed model, single regression coefficient")
  ent <- c(if (tipo == "coeficiente") list(`Tail(s)` = if (caudas == 2) "Two" else "One"), `Effect size f²` = f2, `α err prob` = .ee_poder_fmt(alfa, 6),
           if (is.null(n)) list(`Power (1-β err prob)` = poder) else list(`Total sample size` = n),
           if (tipo == "incremento") list(`Number of tested predictors` = testados, `Total number of predictors` = p) else
             list(`Number of predictors` = p))
  sai <- if (tipo == "coeficiente") list(`Noncentrality parameter δ` = r$nc, `Critical t` = r$critico, Df = r$gl) else
    list(`Noncentrality parameter λ` = r$nc, `Critical F` = r$critico, `Numerator df` = if (tipo == "global") p else testados,
      `Denominator df` = r$gl)
  ef <- .ee_poder_efeito("f2", f2, "f2")
  if (!is.null(r2)) ef <- rbind(ef, .ee_poder_efeito("R2", r2, "f2"))
  if (!is.null(delta_r2)) ef <- rbind(ef, .ee_poder_efeito(c("R2_reduzido", "delta_R2"), c(r2_reduzido, delta_r2), "f2"))
  if (!is.null(r2_parcial)) ef <- rbind(ef, .ee_poder_efeito("R2_parcial", r2_parcial, "f2"))
  .ee_poder_resultado(
    paste0("Poder e tamanho amostral: regressão linear (", tipo, ")"), if (tipo == "coeficiente") "t tests" else "F tests",
    teste, n, entradas = ent, saidas = c(sai, if (is.null(n)) list(`Total sample size` = N), `Actual power` = r$poder),
    skill = list(`N mínimo` = if (is.null(n)) N else NA, pd$skill),
    n_res = N, efeito = ef, decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- regressão logística
# Poder pelo procedimento de Hsieh et al. (1998) na forma fechada (a fórmula do N, com o
# fator 1/(1 - R²), resolvida para z de 1 - β), como o G*Power mostra em "Actual power":
# manual, seção 29.4, 0.950486 com N = 317 (preditor normal) e 0.950068 com N = 1437
# (preditor binário, π = .5). O pwrss 1.3.1 devolve o poder pedido no cálculo a priori e,
# com n dado, acha o poder por uniroot com tolerância de cerca de 1e-4.
.ee_poder_hsieh <- function(or, p0, n, r2, distribuicao, prop_x, alfa, caudas) {
  p2 <- or * p0 / (1 - p0) / (1 + or * p0 / (1 - p0))
  za <- stats::qnorm(1 - alfa / caudas)
  ne <- n * (1 - r2)
  zb <- if (distribuicao == "normal") sqrt(ne * p0 * (1 - p0)) * abs(log(or)) - za else {
    B <- prop_x; pbar <- (1 - B) * p0 + B * p2
    (sqrt(ne * (p0 - p2)^2 * (1 - B)) - za * sqrt(pbar * (1 - pbar) / B)) / sqrt(p0 * (1 - p0) + p2 * (1 - p2) * (1 - B) / B)
  }
  stats::pnorm(zb)
}

# z tests > Logistic regression. Padrão do G*Power: Demidenko (2007) com correção de
# variância, que o pwrss implementa como method = "demidenko(vc)" (227 no curso).
# p0 = Pr(Y = 1 | X = 1) H0 do G*Power: probabilidade do evento com o preditor na média
# (efeito nulo). r2_outros = "R² other X". Preditor normal (contínuo) ou binomial.
poder_logistica <- function(or, p0, r2_outros = 0, distribuicao = c("normal", "binomial"), prop_x = 0.5,
                            metodo = c("demidenko_vc", "demidenko", "hsieh"), n = NULL, poder = 0.80,
                            alfa = 0.05, caudas = 2, perda = 0, justificativa = NULL, decimal = ".", idioma = "pt") {
  .pr_exige("pwrss")
  distribuicao <- match.arg(distribuicao); metodo <- match.arg(metodo)
  .ee_poder_checa(poder, alfa, perda, caudas, n)
  if (or <= 0 || or == 1) stop("or precisa ser positivo e diferente de 1")
  if (p0 <= 0 || p0 >= 1) stop("p0 precisa estar entre 0 e 1")
  decisoes <- character(); avisos <- character()
  met <- c(demidenko_vc = "demidenko(vc)", demidenko = "demidenko", hsieh = "hsieh")[[metodo]]
  # R² other X > 0 com correção de variância. A fórmula documentada no manual (seção 29.6.2,
  # eq. 14-15: m1 = √(N(1 - R²)/v1)·β1 e s1 = √(v0/v1)) dá N/(1 - R²) com a correção
  # mantida: 375 no exemplo da seção 29.4 (OR = 1,5, p0 = .5, poder .95). O programa
  # imprime 395 nesse exemplo, que é o N de Demidenko SEM correção dividido por (1 - R²)
  # (354,74 / 0,9 = 394,15). Regra da skill: com R² > 0, o caminho sem correção (s1 = 1).
  # Única âncora com correção e R² > 0: esse exemplo, com preditor normal. Com preditor
  # binomial não há exemplo publicado: o N sai no bloco da skill, não no da saída do
  # G*Power, e a frase descreve o procedimento sem atribuí-lo ao G*Power.
  vc_r2 <- metodo == "demidenko_vc" && r2_outros > 0
  sem_ancora <- vc_r2 && distribuicao == "binomial"
  met_doc <- met
  if (vc_r2) met <- "demidenko"
  # Preditor binário como binomial com size = 1: mesma densidade da "bernoulli" (e a = 0,85 na
  # correção de variância), mas é a única forma que o Hsieh do pwrss 1.3.1 aceita; com
  # "bernoulli" o metodo = "hsieh" parava com erro (verificação final, G-02).
  dist <- if (distribuicao == "normal") "normal" else list(dist = "binomial", size = 1, prob = prop_x)
  alt <- if (caudas == 2) "two.sided" else "one.sided"
  roda <- function(oo, nn = NULL, ceil = TRUE, mm = met) pwrss::power.z.logistic(
    base.prob = p0, odds.ratio = oo, n = nn, power = if (is.null(nn)) poder else NULL, r.squared.predictor = r2_outros,
    alpha = alfa, alternative = alt, method = mm, distribution = dist, ceil.n = ceil, verbose = 0)
  x <- roda(or, n)
  x_doc <- if (vc_r2) roda(or, n, mm = met_doc) else NULL
  if (vc_r2) {
    doc_txt <- if (is.null(n)) paste0(x_doc$n, " participantes") else paste0("poder ", round(x_doc$power, 4))
    decisoes <- c(decisoes, paste0(
      "Regra para R² other X > 0 com correção de variância: N de Demidenko (2007) sem a correção, dividido por (1 - R²). ",
      "É o que reproduz o único exemplo do manual do G*Power com R² > 0 (seção 29.4: OR = 1,5, p0 = .5, poder .95, R² = .1, ",
      "preditor normal: 395); a fórmula documentada no mesmo manual (eq. 14-15, correção mantida) dá 375 nesse exemplo ",
      "(neste cálculo: ", doc_txt, ")."))
    avisos <- c(avisos, if (sem_ancora) paste0(
      "R² other X > 0 com preditor binário: nenhum exemplo publicado do G*Power foi conferido para este caso. O N ",
      "(regra do exemplo com preditor normal: sem correção de variância, dividido por 1 - R²) não foi conferido contra ",
      "o G*Power e sai no bloco da skill, não em 'G*Power, saída'; pela fórmula documentada seriam ", doc_txt,
      ". Conferir na tela do G*Power antes de citar o programa.") else paste0(
      "R² other X > 0 com correção de variância: a regra (sem correção, N dividido por 1 - R²) foi inferida de um só ",
      "exemplo do manual do G*Power (seção 29.4, 395) e contradiz a fórmula documentada, que daria ", doc_txt,
      ". Conferir o N na tela do G*Power antes de relatar."))
  }
  N <- if (is.null(n)) x$n else n
  if (metodo == "hsieh") {
    x$power <- .ee_poder_hsieh(or, p0, N, r2_outros, distribuicao, prop_x, alfa, caudas)
    decisoes <- c(decisoes, "Poder alcançado do Hsieh pela fórmula de Hsieh et al. (1998) resolvida para o poder, como o G*Power (o pwrss devolve o poder pedido).")
  }
  if (is.null(n)) decisoes <- c(decisoes, paste0("N contínuo = ", round(roda(or, ceil = FALSE)$n, 2), ", arredondado para cima: ", N, "."))
  decisoes <- c(decisoes, paste0(
    "Método: ", c(demidenko_vc = "Demidenko (2007) com correção de variância, o padrão do G*Power 3.1",
                  demidenko = "Demidenko (2007) sem correção de variância", hsieh = "Hsieh et al. (1998)")[[metodo]],
    ". O N é sensível a p0: com OR = 1,52, p0 = 0,05 dá 898 e 0,30 dá 227 (curso, seção 24); p0 deve refletir a ",
    "prevalência esperada do desfecho na amostra."))
  avisos <- c(avisos, paste0(
    "Os cortes de OR de Chen et al. (2010; 1,52, 2,74 e 4,72 com incidência de 5%) são para exposição binária; ",
    "com preditor contínuo o OR é por 1 DP do preditor e a tabela não se aplica direto (curso, seção 24)."))
  if (distribuicao == "binomial" && metodo == "demidenko_vc") decisoes <- c(decisoes, "Preditor binomial: correção de variância com a = 0,85, conferida contra os exemplos do manual do G*Power (π = .5, .2 e .8: 1437, 2158 e 2368).")
  pd <- .ee_poder_perda(N, perda, 1, decisoes, avisos); avisos <- pd$avisos
  ef_txt <- paste0("OR = ", .pr_num(or, 2, TRUE, decimal), " para ", if (distribuicao == "normal") "um preditor contínuo (por 1 DP)" else "um preditor binário",
                   " e probabilidade do evento sob H0 de ", .pr_num(p0, 2, FALSE, decimal))
  if (idioma == "en") ef_txt <- paste0("OR = ", .pr_num(or, 2, TRUE, decimal), " for a ", if (distribuicao == "normal") "continuous predictor (per 1 SD)" else "binary predictor",
                                       " and event probability under H0 of ", .pr_num(p0, 2, FALSE, decimal))
  if (r2_outros > 0) ef_txt <- paste0(ef_txt, if (idioma == "en") ", with R² = " else ", com R² = ", .pr_num(r2_outros, 2, FALSE, decimal),
                                      if (idioma == "en") " for the predictor on the other predictors" else " do preditor com os demais preditores")
  abertura <- if (sem_ancora) c(
    pt = "O tamanho amostral foi determinado a priori pela aproximação de Demidenko (2007) para o teste de Wald, sem correção de variância e com o N dividido por 1 - R² dos demais preditores (Hsieh et al., 1998), ",
    en = "Sample size was determined a priori with Demidenko's (2007) approximation for the Wald test, without variance correction and with N divided by 1 - R² for the other predictors (Hsieh et al., 1998), ") else NULL
  apa <- if (is.null(n)) .ee_poder_apa("para uma regressão logística binária", "for a binary logistic regression", ef_txt, alfa, caudas,
                                       poder, N, "", perda, pd$n_perda, justificativa, decimal, idioma, abertura = abertura) else character()
  ent <- list(`Tail(s)` = if (caudas == 2) "Two" else "One", `Odds ratio` = or, `Pr(Y=1|X=1) H0` = p0, `α err prob` = .ee_poder_fmt(alfa, 6),
           if (is.null(n)) list(`Power (1-β err prob)` = poder) else list(`Total sample size` = n), `R² other X` = r2_outros,
           `X distribution` = if (distribuicao == "normal") "Normal" else "Binomial",
           if (distribuicao == "normal") list(`X parm μ` = 0, `X parm σ` = 1) else list(`X parm π` = prop_x),
           `Options` = c(demidenko_vc = "Demidenko (2007) with var corr", demidenko = "Demidenko (2007) without var corr",
                         hsieh = "Hsieh et al. (1998)")[[metodo]])
  .ee_poder_resultado(
    "Poder e tamanho amostral: regressão logística", "z tests", "Logistic regression", n, entradas = ent,
    saidas = list(`Critical z` = stats::qnorm(1 - alfa / caudas),
                  if (!sem_ancora) list(if (is.null(n)) list(`Total sample size` = N), `Actual power` = x$power)),
    skill = list(`N mínimo` = if (is.null(n)) N else NA,
                 if (sem_ancora) list(`Poder (regra não conferida no G*Power)` = x$power),
                 if (vc_r2) list(`N pela fórmula documentada (correção mantida)` = if (is.null(n)) x_doc$n else NA,
                                 `Poder pela fórmula documentada` = if (is.null(n)) NA else x_doc$power), pd$skill),
    n_res = N, efeito = .ee_poder_efeito(c("OR", "p0"), c(or, p0), "OR"), decisoes = decisoes, avisos = avisos, apa = apa,
    extra = list(poder_real = x$power))
}

# ---------------------------------------------------------------- qui-quadrado
# χ² tests > Goodness-of-fit tests: Contingency tables (λ = w²N; gl informado). Igual a
# pwr::pwr.chisq.test. Para tabela r × c, gl = (r - 1)(c - 1); w = V·√(min(r, c) - 1).
poder_quiquadrado <- function(w, gl, n = NULL, poder = 0.80, alfa = 0.05, perda = 0, justificativa = NULL,
                              decimal = ".", idioma = "pt") {
  .ee_poder_checa(poder, alfa, perda, 2, n)
  if (w <= 0) stop("w precisa ser positivo")
  fp <- function(N, ww = w) { crit <- stats::qchisq(1 - alfa, gl)
    list(poder = stats::pchisq(crit, gl, ncp = ww^2 * N, lower.tail = FALSE), critico = crit) }
  N <- if (is.null(n)) .ee_poder_busca(function(NN) fp(NN)$poder, poder, 2) else n
  r <- fp(N)
  decisoes <- character(); avisos <- character()
  pd <- .ee_poder_perda(N, perda, 1, decisoes, avisos); avisos <- pd$avisos
  calc <- function(ww) .ee_poder_busca(function(NN) fp(NN, ww)$poder, poder, 2)
  k <- if (is.null(n)) .ee_poder_casas(w, calc, N) else 2
  apa <- if (is.null(n)) .ee_poder_apa(paste0("para um teste qui-quadrado com ", gl, " gl"), paste0("for a chi-square test with ", gl, " df"),
                                       paste0("w = ", .pr_num(w, k, TRUE, decimal)), alfa, NULL, poder, N, "", perda, pd$n_perda,
                                       justificativa, decimal, idioma) else character()
  .ee_poder_resultado(
    "Poder e tamanho amostral: qui-quadrado", "χ² tests", "Goodness-of-fit tests: Contingency tables", n,
    entradas = list(`Effect size w` = w, `α err prob` = .ee_poder_fmt(alfa, 6), if (is.null(n)) list(`Power (1-β err prob)` = poder) else list(`Total sample size` = n), Df = gl),
    saidas = list(`Noncentrality parameter λ` = w^2 * N, `Critical χ²` = r$critico, Df = gl, if (is.null(n)) list(`Total sample size` = N),
               `Actual power` = r$poder),
    skill = list(`N mínimo` = if (is.null(n)) N else NA, pd$skill),
    n_res = N, efeito = .ee_poder_efeito("w", w, "w"), decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- AFC e MEE (v1.4)
# Kelley e Lai (2018, cap. 5 do Wiley Handbook) dão quatro abordagens; as quatro estão
# aqui, com os exemplos do capítulo nos evals:
#   1. poder para o RMSEA do modelo (MacCallum, Browne & Sugawara, 1996):
#      poder_afc_rmsea(); λ = (N - 1) gl ε² / grupos, como semTools::findRMSEApower;
#   2. AIPE do RMSEA (Kelley & Lai, 2011): n_aipe_rmsea(), N para a largura esperada
#      do IC do RMSEA não passar de `largura`;
#   3. poder para um parâmetro por comparação de modelos aninhados (Satorra & Saris,
#      1985, com Σ e os modelos; MacCallum, Browne & Cai, 2006, com os RMSEA ou os F):
#      poder_afc_aninhados();
#   4. AIPE de um parâmetro (largura do IC de Wald): n_aipe_parametro().
# Todas devolvem o menor N inteiro que atinge o alvo. O MBESS (usado no capítulo) para a
# busca quando o poder fica a .001 do alvo; por isso o capítulo dá 1.269, 1.113 e 364
# onde o N mínimo exato é 1.268, 1.116 e 365 (o semTools dá os exatos).
.ee_rmsea_pw <- function(N, e0, e1, gl, alfa, grupos = 1) {
  l0 <- (N - 1) * gl * e0^2 / grupos; l1 <- (N - 1) * gl * e1^2 / grupos
  if (e0 < e1) {
    cr <- stats::qchisq(1 - alfa, gl, ncp = l0)
    list(poder = stats::pchisq(cr, gl, ncp = l1, lower.tail = FALSE), critico = cr, ncp0 = l0, ncp1 = l1)
  } else {
    cr <- stats::qchisq(alfa, gl, ncp = l0)
    list(poder = stats::pchisq(cr, gl, ncp = l1), critico = cr, ncp0 = l0, ncp1 = l1)
  }
}

.ee_tabela_campos <- function(...) {
  x <- list(...); x <- x[!vapply(x, is.null, TRUE)]
  data.frame(campo = names(x), valor = vapply(x, function(v) if (is.numeric(v)) format(signif(v, 6), scientific = FALSE, trim = TRUE) else as.character(v), ""),
             stringsAsFactors = FALSE, row.names = NULL)
}

poder_afc_rmsea <- function(rmsea0 = 0.05, rmsea1 = 0.08, gl, n = NULL, poder = 0.80, alfa = 0.05, grupos = 1,
                            perda = 0, justificativa = NULL, decimal = ".", idioma = "pt") {
  .ee_poder_checa(poder, alfa, perda, 2, n)
  if (!is.numeric(gl) || gl < 1) stop("gl (graus de liberdade do modelo) precisa ser >= 1")
  if (rmsea0 < 0 || rmsea1 < 0 || rmsea0 == rmsea1) stop("rmsea0 e rmsea1 precisam ser >= 0 e diferentes")
  teste <- if (rmsea0 == 0) "ajuste exato" else if (rmsea0 < rmsea1) "ajuste aproximado (close fit)" else "ajuste não aproximado (not-close fit)"
  h0 <- if (rmsea0 == 0) "ε = 0" else if (rmsea0 < rmsea1) paste0("ε ≤ ", rmsea0) else paste0("ε ≥ ", rmsea0)
  fp <- function(NN) .ee_rmsea_pw(NN, rmsea0, rmsea1, gl, alfa, grupos)$poder
  N <- if (is.null(n)) .ee_poder_busca(fp, poder, 2) else n
  r <- .ee_rmsea_pw(N, rmsea0, rmsea1, gl, alfa, grupos)
  avisos <- character(); decisoes <- c(
    paste0("Teste de ", teste, " de MacCallum, Browne e Sugawara (1996): H0 ", h0, ", RMSEA verdadeiro ", rmsea1,
           ", gl = ", gl, ", α = ", alfa, if (grupos > 1) paste0(", ", grupos, " grupos (λ dividido pelo número de grupos, como no semTools; N é o total, ",
                                                              "e o semTools devolve N / grupos)") else "", "."),
    "N mínimo exato (menor inteiro com poder >= alvo), igual a semTools::findRMSEAsamplesize; o MBESS::ss.power.sem para a busca a .001 do alvo.",
    "O gl é o do modelo que vai ser testado: contar antes (p(p + 1)/2 menos os parâmetros livres) ou ler de um ajuste com dados simulados.")
  if (rmsea0 == 0) avisos <- c(avisos, "Rejeitar ε = 0 só diz que o ajuste não é perfeito, o que já se sabe antes de coletar; Kelley e Lai (2018) preferem H0 de ajuste mínimo aceitável (por exemplo ε ≥ .08).")
  pd <- .ee_poder_perda(N, perda, 1, decisoes, avisos); avisos <- pd$avisos
  tab <- .ee_tabela_campos(`Teste` = teste, `H0` = h0, `RMSEA verdadeiro` = rmsea1, gl = gl, `α` = alfa,
                           `Poder pedido` = if (is.null(n)) poder, `Grupos` = if (grupos > 1) grupos,
                           `λ sob H0` = r$ncp0, `λ sob H1` = r$ncp1, `Qui-quadrado crítico` = r$critico,
                           `N` = N, `N com grupos iguais` = if (grupos > 1) grupos * ceiling(N / grupos), `Poder alcançado` = r$poder,
                           `N para coleta (perda)` = if (perda > 0) pd$n_perda)
  just <- if (is.null(justificativa)) (if (idioma == "en") "[source of the RMSEA values: fill in]" else "[origem dos valores de RMSEA: preencher]") else justificativa
  apa <- if (!is.null(n)) character() else if (idioma == "en") {
    paste0("Sample size was planned for the test of ", c(`ajuste exato` = "exact fit", `ajuste aproximado (close fit)` = "close fit",
           `ajuste não aproximado (not-close fit)` = "not-close fit")[[teste]], " (MacCallum, Browne, & Sugawara, 1996), with H0: ", sub("ε", "RMSEA", h0),
           ", a true RMSEA of ", rmsea1, " (", just, "), df = ", gl, ", α = ", .ee_poder_alfa_txt(alfa, decimal), " and power of ",
           .ee_poder_alfa_txt(poder, decimal), ", which required a minimum of ", N, " participants.",
           if (perda > 0) paste0(" Anticipating ", round(100 * perda), "% attrition, the recruitment target is ", pd$n_perda, " participants.") else "")
  } else {
    paste0("O tamanho amostral foi planejado para o teste de ", teste, " (MacCallum, Browne & Sugawara, 1996), com H0: ", sub("ε", "RMSEA", h0),
           ", RMSEA verdadeiro de ", .pr_num(rmsea1, 3, FALSE, decimal), " (", just, "), gl = ", gl, ", α = ", .ee_poder_alfa_txt(alfa, decimal),
           " e poder de ", .ee_poder_alfa_txt(poder, decimal), ", resultando em um mínimo de ", N, " participantes.",
           if (perda > 0) paste0(" Prevendo perda de ", round(100 * perda), "%, a meta de coleta é de ", pd$n_perda, " participantes.") else "")
  }
  .ee_resultado("Poder e tamanho amostral: RMSEA da AFC", n = N, tabela = tab,
                efeito = data.frame(indice = c("RMSEA sob H0", "RMSEA verdadeiro"), valor = c(rmsea0, rmsea1), ic_inf = NA_real_,
                                    ic_sup = NA_real_, principal = c(FALSE, TRUE)),
                decisoes = decisoes, avisos = avisos, apa = apa, extra = list(poder_alcancado = r$poder, n_minimo = N))
}

# IC do RMSEA tomando o valor populacional como estimativa (T = ε² gl (N - 1) + gl),
# como MBESS::ci.rmsea; limites pela inversão da qui-quadrado não central.
.ee_ic_rmsea <- function(rmsea, gl, N, nivel) {
  a <- (1 - nivel) / 2; T <- rmsea^2 * gl * (N - 1) + gl
  acha <- function(alvo) {
    if (stats::pchisq(T, gl, ncp = 0) < alvo) return(0)
    hi <- max(10, T); while (stats::pchisq(T, gl, ncp = hi) > alvo) hi <- hi * 2
    stats::uniroot(function(l) stats::pchisq(T, gl, ncp = l) - alvo, c(0, hi), tol = 1e-10)$root
  }
  lo <- acha(1 - a); up <- acha(a)
  c(inf = sqrt(lo / (gl * (N - 1))), sup = sqrt(up / (gl * (N - 1))))
}

n_aipe_rmsea <- function(rmsea, gl, largura, nivel = 0.90, decimal = ".", idioma = "pt") {
  if (rmsea < 0 || largura <= 0 || gl < 1) stop("rmsea >= 0, largura > 0 e gl >= 1")
  w <- function(N) { ic <- .ee_ic_rmsea(rmsea, gl, N, nivel); ic[["sup"]] - ic[["inf"]] }
  N <- .ee_poder_busca(function(NN) -w(NN), -largura, gl + 2)
  ic <- .ee_ic_rmsea(rmsea, gl, N, nivel)
  tab <- .ee_tabela_campos(`RMSEA populacional` = rmsea, gl = gl, `Nível do IC` = nivel, `Largura desejada` = largura,
                           N = N, `IC esperado: inferior` = ic[["inf"]], `IC esperado: superior` = ic[["sup"]],
                           `Largura esperada` = ic[["sup"]] - ic[["inf"]])
  decisoes <- c("AIPE do RMSEA (Kelley & Lai, 2011): menor N com largura esperada do IC <= largura, com o RMSEA populacional no lugar da estimativa (como MBESS::ss.aipe.rmsea).",
                "Largura esperada: em cerca de metade dos estudos com esse N o IC sai mais largo (Kelley & Lai, 2018); para garantia maior, planejar com largura menor.")
  apa <- if (idioma == "en") paste0("Sample size was planned for accuracy in parameter estimation of the RMSEA (Kelley & Lai, 2011): assuming a population RMSEA of ", rmsea,
                                    " and df = ", gl, ", N = ", N, " yields an expected ", round(100 * nivel), "% CI no wider than ", largura, ".") else
    paste0("O tamanho amostral foi planejado pela acurácia da estimativa do RMSEA (Kelley & Lai, 2011): com RMSEA populacional de ", .pr_num(rmsea, 3, FALSE, decimal),
           " e gl = ", gl, ", N = ", N, " dá IC ", round(100 * nivel), "% com largura esperada de até ", .pr_num(largura, 3, FALSE, decimal), ".")
  .ee_resultado("Tamanho amostral: AIPE do RMSEA", n = N, tabela = tab, decisoes = decisoes, apa = apa,
                extra = list(n_minimo = N, ic_esperado = ic))
}

# Modelos aninhados. Entrada por RMSEA (MacCallum, Browne & Cai, 2006), por F de máxima
# verossimilhança, ou por Σ populacional mais as sintaxes lavaan dos dois modelos (Satorra
# & Saris, 1985: F = mínimo da discrepância ML de cada modelo ajustado a Σ). λ = (N - 1)
# (F_restrito - F_completo), gl = diferença de gl, valor crítico da qui-quadrado central.
poder_afc_aninhados <- function(gl_restrito = NULL, gl_completo = NULL, rmsea_restrito = NULL, rmsea_completo = 0,
                                F_restrito = NULL, F_completo = 0, Sigma = NULL, modelo_restrito = NULL, modelo_completo = NULL,
                                n = NULL, poder = 0.80, alfa = 0.05, perda = 0, justificativa = NULL, decimal = ".", idioma = "pt") {
  .ee_poder_checa(poder, alfa, perda, 2, n)
  via <- if (!is.null(Sigma)) "Sigma" else if (!is.null(F_restrito)) "F" else if (!is.null(rmsea_restrito)) "RMSEA" else
    stop("informar Sigma com as duas sintaxes, ou F_restrito, ou rmsea_restrito (com gl_restrito e gl_completo)")
  if (via == "Sigma") {
    .pr_exige("lavaan")
    if (is.null(modelo_restrito) || is.null(modelo_completo)) stop("com Sigma, informar modelo_restrito e modelo_completo (sintaxe lavaan)")
    aj <- function(m) lavaan::cfa(m, sample.cov = Sigma, sample.nobs = 1000, sample.cov.rescale = FALSE)
    fr <- aj(modelo_restrito); fc <- aj(modelo_completo)
    F_restrito <- 2 * unname(lavaan::fitMeasures(fr, "fmin")); F_completo <- 2 * unname(lavaan::fitMeasures(fc, "fmin"))
    gl_restrito <- unname(lavaan::fitMeasures(fr, "df")); gl_completo <- unname(lavaan::fitMeasures(fc, "df"))
  }
  if (is.null(gl_restrito) || is.null(gl_completo)) stop("informar gl_restrito e gl_completo")
  gd <- gl_restrito - gl_completo
  if (gd < 1) stop("o modelo restrito precisa ter mais gl que o completo")
  Fd <- if (via == "RMSEA") gl_restrito * rmsea_restrito^2 - gl_completo * rmsea_completo^2 else F_restrito - F_completo
  if (Fd < 1e-6) stop("a discrepância do restrito precisa ser maior que a do completo (diferença de F = ", signif(Fd, 4),
                      "): com ΔF nulo a restrição vale em Σ e o poder fica igual a α")
  cr <- stats::qchisq(1 - alfa, gd)
  fp <- function(NN) stats::pchisq(cr, gd, ncp = (NN - 1) * Fd, lower.tail = FALSE)
  N <- if (is.null(n)) .ee_poder_busca(fp, poder, 2) else n
  pw <- fp(N)
  avisos <- character()
  decisoes <- c(paste0(c(Sigma = "Satorra e Saris (1985): F de cada modelo ajustado a Σ por máxima verossimilhança (lavaan, sem reescalar Σ)",
                         F = "Diferença de F informada (Satorra & Saris, 1985)",
                         RMSEA = "MacCallum, Browne e Cai (2006): F = gl × RMSEA² de cada modelo")[[via]],
                       "; λ = (N - 1) × ΔF, gl = ", gd, ", qui-quadrado crítico central."),
                "Teste da razão de verossimilhança entre os modelos: H0 é a restrição valer na população.")
  pd <- .ee_poder_perda(N, perda, 1, decisoes, avisos); avisos <- pd$avisos
  tab <- .ee_tabela_campos(`Entrada` = via, `F restrito` = if (via != "RMSEA") F_restrito, `F completo` = if (via != "RMSEA") F_completo,
                           `RMSEA restrito` = if (via == "RMSEA") rmsea_restrito, `RMSEA completo` = if (via == "RMSEA") rmsea_completo,
                           `gl restrito` = gl_restrito, `gl completo` = gl_completo, `ΔF` = Fd, `Δgl` = gd, `α` = alfa,
                           `Poder pedido` = if (is.null(n)) poder, `λ` = (N - 1) * Fd, `Qui-quadrado crítico` = cr, N = N,
                           `Poder alcançado` = pw, `N para coleta (perda)` = if (perda > 0) pd$n_perda)
  just <- if (is.null(justificativa)) (if (idioma == "en") "[source of the population values: fill in]" else "[origem dos valores populacionais: preencher]") else justificativa
  apa <- if (!is.null(n)) character() else if (idioma == "en")
    paste0("Sample size for the likelihood-ratio test between nested models (Δdf = ", gd, ") was planned with the ",
           if (via == "RMSEA") "MacCallum, Browne, and Cai (2006)" else "Satorra and Saris (1985)", " method (", just, "), α = ",
           .ee_poder_alfa_txt(alfa, decimal), " and power of ", .ee_poder_alfa_txt(poder, decimal), ", which required a minimum of ", N, " participants.") else
    paste0("O tamanho amostral para o teste da razão de verossimilhança entre modelos aninhados (Δgl = ", gd, ") foi planejado pelo método de ",
           if (via == "RMSEA") "MacCallum, Browne e Cai (2006)" else "Satorra e Saris (1985)", " (", just, "), com α = ", .ee_poder_alfa_txt(alfa, decimal),
           " e poder de ", .ee_poder_alfa_txt(poder, decimal), ", resultando em um mínimo de ", N, " participantes.")
  .ee_resultado("Poder e tamanho amostral: modelos aninhados (AFC/MEE)", n = N, tabela = tab, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(n_minimo = N, poder_alcancado = pw, delta_F = Fd))
}

# AIPE de um parâmetro: ajusta o modelo a Σ (N de referência 1000, sem reescalar) e usa
# EP(N) = EP(1000) × raiz(1000 / N), a escala do erro padrão assintótico da ML; menor N
# com largura do IC de Wald (2 z EP) <= largura. Reproduz o N = 1.007 do capítulo
# (λ1, largura .15, IC 95%). A garantia (assurance) do MBESS não está implementada.
n_aipe_parametro <- function(modelo, Sigma, parametro, largura, nivel = 0.95, decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  fit <- lavaan::cfa(modelo, sample.cov = Sigma, sample.nobs = 1000, sample.cov.rescale = FALSE)
  pe <- lavaan::parameterEstimates(fit)
  rot <- paste0(pe$lhs, pe$op, pe$rhs)
  i <- which(pe$label == parametro | rot == gsub(" ", "", parametro))
  if (length(i) != 1) stop("parâmetro '", parametro, "' não encontrado: usar o rótulo da sintaxe ou a forma 'F1=~X1'")
  ep0 <- pe$se[i]; z <- stats::qnorm(1 - (1 - nivel) / 2)
  N <- ceiling(1000 * (2 * z * ep0 / largura)^2)
  while (2 * z * ep0 * sqrt(1000 / (N - 1)) <= largura && N > 2) N <- N - 1
  tab <- .ee_tabela_campos(`Parâmetro` = parametro, `Valor populacional` = pe$est[i], `EP com N = 1000` = ep0, `Nível do IC` = nivel,
                           `Largura desejada` = largura, N = N, `Largura esperada` = 2 * z * ep0 * sqrt(1000 / N))
  apa <- if (idioma == "en") paste0("Sample size was planned for accuracy in parameter estimation of ", parametro, ": N = ", N, " yields an expected ",
                                    round(100 * nivel), "% confidence interval no wider than ", largura, " (Kelley & Lai, 2018).") else
    paste0("O tamanho amostral foi planejado pela acurácia da estimativa de ", parametro, ": N = ", N, " dá IC ", round(100 * nivel),
           "% com largura esperada de até ", .pr_num(largura, 2, FALSE, decimal), " (Kelley & Lai, 2018).")
  .ee_resultado("Tamanho amostral: AIPE de um parâmetro", n = N, tabela = tab,
                decisoes = c("Σ populacional informado pelo pesquisador (correlações e desvios-padrão plausíveis, ou o implicado por cargas propostas).",
                             "Largura esperada: cerca de metade dos estudos com esse N terá IC mais largo; a versão com garantia (assurance) está no MBESS::ss.aipe.sem.path."),
                apa = apa, extra = list(n_minimo = N))
}
