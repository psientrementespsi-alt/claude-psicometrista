# est_associacao.R: módulo de associação da skill psicometria-r (v1.3.0).
# Correlação (Pearson, Spearman, Kendall, ponto-bisserial, parcial, comparação de
# correlações), regressão linear simples, múltipla e hierárquica, regressão logística
# binária, curva ROC, qui-quadrado (aderência e independência), McNemar e Q de Cochran.
#
# Carregar por scripts/carregar.R (depende de est_comum.R e psicometria.R).
# Padrões seguem o livro do usuário (JASP 0.95.1); alternativas da literatura entram
# como argumento e aparecem ao lado quando custam pouco. Documentação para a IA:
# references/associacao.md.
#
# Funções públicas (assinaturas resumidas; todas com decimal = "." e idioma = "pt").
# casas = 3 nas frases (B, β, OR, r, V, φ, AUC), como os parágrafos-modelo do livro.
#   correlacoes(dados, vars, metodos = c("pearson", "spearman"), ausentes = "pairwise",
#               nivel = .95, regua = "cohen", limiar_divergencia = .10, casas = 3)
#   correlacao_ponto_bisserial(dados, continua, dicotomica, codigo_1 = NULL, nivel, casas = 3)
#   correlacao_parcial(dados, x, y, controles, metodo = "pearson", nivel, casas = 3)
#   comparar_correlacoes(r1, n1, r2, n2 = NULL, r12 = NULL, nivel, casas = 3)
#   regressao_linear(dados, vd, preditores, referencias = NULL, fatores = NULL, nivel = .95,
#                    tipo_residuo = "studentizado", ep_robusto = FALSE, casas = 3)
#   regressao_hierarquica(dados, vd, blocos, referencias = NULL, fatores = NULL, nivel = .95, ...)
#   regressao_logistica(dados, vd, preditores, evento = NULL, referencias = NULL, fatores = NULL,
#                       corte = .5, ic = "wald", box_tidwell = TRUE, nivel = .95, casas = 3)
#   curva_roc(dados, desfecho, preditor, evento = NULL, maior_indica_evento = TRUE, prevalencia, custo_fn_fp,
#             cortes = NULL, nivel = .95, casas = 3)
#   qui_quadrado_aderencia(dados, variavel, proporcoes = NULL, semente = 2026, simulacoes = 10000, nivel, casas = 3)
#   qui_quadrado(dados, linhas, colunas, evento = NULL, referencia = NULL, exato = "auto",
#                correcao_continuidade = FALSE, alfa_residuos = .05, semente = 2026,
#                simulacoes = 10000, nivel = .95, casas = 3)
#   mcnemar(dados, antes, depois, evento = NULL, teste = "assintotico", nivel = .95, casas = 3)
#   vies_preditivo(dados, criterio, preditor, grupo, referencia = NULL, alfa = .05)   (Cleary, v1.4)
#   q_cochran(dados, vars, evento = NULL, ajuste = "bonferroni", nivel = .95)
# `referencias`: list(fator = "nível") ou c(fator = "nível"). Evento sem `evento`: 0/1 -> 1,
# fator -> 2º nível, texto -> erro pedindo `evento` (ver .ee_assoc_binaria).
# `fatores`: preditores numéricos que são grupo (1/2/3 vindo do SPSS ou de planilha) e
# devem virar indicadoras; sem isso, preditor numérico com 3 a 7 valores inteiros entra
# como contínuo, com aviso (ver .ee_assoc_aviso_codigo).

# ================================================================ auxiliares internas
.ee_assoc_tr <- function(pt, en, idioma) if (identical(idioma, "en")) en else pt

# Células de texto vazias ("") viram NA antes da deleção listwise (o JASP trata célula
# vazia como ausente; read.csv sem na.strings = "" deixaria uma categoria fantasma).
.ee_assoc_vazios <- function(dados, vars) {
  if (!is.data.frame(dados)) stop("dados precisa ser um data.frame")
  falta <- setdiff(vars, names(dados))
  if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "))
  n_vazios <- 0
  for (v in vars) {
    x <- dados[[v]]
    if (is.character(x)) {
      vaz <- !is.na(x) & !nzchar(trimws(x))
      n_vazios <- n_vazios + sum(vaz); x[vaz] <- NA; dados[[v]] <- x
    } else if (is.factor(x)) {
      vaz <- !is.na(x) & !nzchar(trimws(as.character(x)))
      if (any(vaz) || any(!nzchar(trimws(levels(x))))) {
        # tira só o nível vazio; níveis sem casos continuam (a aderência precisa deles)
        n_vazios <- n_vazios + sum(vaz); x[vaz] <- NA
        dados[[v]] <- factor(x, levels = levels(x)[nzchar(trimws(levels(x)))])
      }
    }
  }
  attr(dados, "n_vazios") <- n_vazios
  dados
}

.ee_assoc_pct <- function(x, casas = 1, decimal = ".") paste0(.pr_num(100 * x, casas, TRUE, decimal), "%")

# Efeito com IC no nível pedido (.ee_efeito_txt fixa o rótulo em 95%).
.ee_assoc_ef <- function(simbolo, valor, lo = NA, hi = NA, limitado = TRUE, casas = 2,
                         decimal = ".", idioma = "pt", nivel = .95) {
  casas <- .ee_assoc_casas(valor, casas)
  s <- paste0(simbolo, " = ", .pr_num(valor, casas, !limitado, decimal))
  if (!is.na(lo) && !is.na(hi))
    s <- paste0(s, .ee_sep(decimal), .ee_ic(lo, hi, casas, !limitado, decimal, nivel, idioma))
  s
}

.ee_assoc_novo <- function(..., tabelas_extra = list()) {
  r <- .ee_resultado(...)
  r$tabelas_extra <- tabelas_extra
  class(r) <- c("pr_associacao", class(r))
  r
}

# As tabelas extras (modelo global, classificação, frequências) são impressas pelo
# print.pr_resultado da base comum, entre o resultado e as decisões.

# Réguas de magnitude para |r| (livro: réguas de partida, não leis).
.ee_assoc_regua_r <- function(r, regua = "cohen", idioma = "pt") {
  a <- abs(r)
  if (is.na(a)) return(NA_character_)
  rot <- switch(regua,
    cohen = if (a < .10) c("desprezível", "negligible") else if (a < .30) c("fraca", "weak") else
      if (a < .50) c("moderada", "moderate") else c("forte", "strong"),
    schober = if (a < .10) c("desprezível", "negligible") else if (a < .40) c("fraca", "weak") else
      if (a < .70) c("moderada", "moderate") else if (a < .90) c("forte", "strong") else c("muito forte", "very strong"),
    gignac = if (a < .10) c("muito pequena", "very small") else if (a < .20) c("relativamente pequena", "relatively small") else
      if (a < .30) c("típica", "typical") else c("relativamente grande", "relatively large"))
  .ee_assoc_tr(rot[1], rot[2], idioma)
}

.ee_assoc_nome_regua <- function(regua, idioma = "pt") {
  switch(regua,
    cohen = .ee_assoc_tr("régua de Cohen (1988): .10/.30/.50", "Cohen's (1988) benchmarks: .10/.30/.50", idioma),
    schober = .ee_assoc_tr("régua de Schober et al. (2018): .10/.40/.70/.90", "Schober et al.'s (2018) benchmarks: .10/.40/.70/.90", idioma),
    gignac = .ee_assoc_tr("régua de Gignac e Szodorai (2016): .10/.20/.30", "Gignac and Szodorai's (2016) benchmarks: .10/.20/.30", idioma))
}

.ee_assoc_autor_regua <- function(regua, idioma = "pt") {
  switch(regua, cohen = "Cohen (1988)", schober = "Schober et al. (2018)",
         gignac = .ee_assoc_tr("Gignac e Szodorai (2016)", "Gignac and Szodorai (2016)", idioma))
}

# IC de uma correlação pela transformação z de Fisher; k = número de variáveis
# controladas (correlação parcial). Mesmo método do cor.test e do JASP para Spearman.
.ee_assoc_ic_fisher <- function(r, n, nivel = .95, k = 0) {
  if (is.na(r) || n - 3 - k <= 0 || abs(r) >= 1) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - nivel) / 2)
  tanh(atanh(r) + c(-1, 1) * z / sqrt(n - 3 - k))
}

# IC assintótico do tau-b de Kendall (Hollander, Wolfe & Chicken, seções 8.3-8.4),
# como o JASP 0.95.1 (concordanceVector_cpp + correlation.R).
.ee_assoc_ic_kendall <- function(x, y, tau, nivel = .95) {
  n <- length(x)
  if (n < 4 || is.na(tau)) return(c(NA_real_, NA_real_))
  C <- rowSums(sign(outer(x, x, "-")) * sign(outer(y, y, "-")))
  s2 <- (2 * (n - 2) * stats::var(C) / (n * (n - 1)) + 1 - tau^2) * 2 / (n * (n - 1))
  z <- stats::qnorm(1 - (1 - nivel) / 2)
  tau + c(-1, 1) * z * sqrt(s2)
}

# Variável binária: devolve o vetor 0/1, os rótulos e a frase da codificação.
# Regra do módulo para o evento (código 1) quando `evento` não é informado:
#   0/1 numérico -> 1; lógico -> TRUE; outro numérico -> o maior valor; fator -> o
#   segundo nível (a ordem dos níveis é do usuário, como no glm); texto -> exige `evento`
#   (a ordem alfabética não diz qual categoria é o evento: "controle" viria depois de
#   "caso"). Com exige_texto = FALSE (ponto-bisserial, qui-quadrado), texto segue a ordem
#   alfabética e a frase diz isso.
.ee_assoc_binaria <- function(x, evento = NULL, nome = "desfecho", exige_texto = TRUE) {
  if (is.logical(x)) x <- as.integer(x)
  vals <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(x[!is.na(x)]))
  vals <- as.character(vals)
  if (length(vals) != 2) stop("a variável ", nome, " precisa ter exatamente 2 categorias com dados (tem ", length(vals), ")")
  regra <- NULL
  if (is.null(evento)) {
    if (all(vals %in% c("0", "1"))) { evento <- "1"; regra <- "código 1" } else if (is.numeric(x)) {
      evento <- vals[2]; regra <- "maior valor numérico"
    } else if (is.factor(x)) { evento <- vals[2]; regra <- "segundo nível do fator" } else if (exige_texto) {
      stop("a variável ", nome, " é texto ('", vals[1], "', '", vals[2], "'): informe evento = '...' com a categoria do evento; ",
           "a ordem alfabética não diz qual é o evento (ou transforme em fator com o evento no segundo nível)", call. = FALSE)
    } else { evento <- vals[2]; regra <- "ordem alfabética do texto" }
  } else regra <- "informado pelo usuário"
  evento <- as.character(evento)
  if (!evento %in% vals) stop("evento '", evento, "' não é uma categoria de ", nome, " (categorias: ", paste(vals, collapse = ", "), ")")
  outro <- setdiff(vals, evento)
  y <- ifelse(is.na(x), NA_integer_, as.integer(as.character(x) == evento))
  list(y = y, evento = evento, outro = outro,
       frase = paste0(nome, ": evento (código 1) = '", evento, "' (", regra, "); categoria de comparação (código 0) = '", outro, "'"))
}

# Casas adaptativas: se o valor arredonda para zero nas casas pedidas, usa casas
# suficientes para 2 algarismos significativos (sem isso, -0,0000089 sairia "0,000").
.ee_assoc_casas <- function(x, casas) {
  # |x| < 1e-10 é zero numérico (por exemplo, r de variáveis residualizadas): sem casas extras
  if (length(x) != 1 || !is.finite(x) || abs(x) < 1e-10 || round(abs(x), casas) > 0) return(casas)
  max(casas, ceiling(-log10(abs(x))) + 1)
}
.ee_assoc_n <- function(x, casas, zero = TRUE, decimal = ".") .pr_num(x, .ee_assoc_casas(x, casas), zero, decimal)

.ee_assoc_bq <- function(v) ifelse(make.names(v) == v, v, paste0("`", v, "`"))

.ee_assoc_formula <- function(vd, preds) {
  stats::as.formula(paste(.ee_assoc_bq(vd), "~", if (length(preds)) paste(.ee_assoc_bq(preds), collapse = " + ") else "1"))
}

# ================================================================ correlações
correlacoes <- function(dados, vars, metodos = c("pearson", "spearman"),
                        ausentes = c("pairwise", "listwise"), nivel = .95,
                        regua = c("cohen", "schober", "gignac"), limiar_divergencia = .10, casas = 3,
                        decimal = ".", idioma = "pt") {
  metodos <- match.arg(metodos, c("pearson", "spearman", "kendall"), several.ok = TRUE)
  ausentes <- match.arg(ausentes); regua <- match.arg(regua)
  if (length(vars) < 2) stop("informe pelo menos 2 variáveis")
  dados <- .ee_assoc_vazios(dados, vars)
  chk <- .ee_dados(dados, vars, numericas = vars)
  base <- if (ausentes == "listwise") chk$d else as.data.frame(dados)[, vars, drop = FALSE]
  pares <- utils::combn(vars, 2, simplify = FALSE)
  linhas <- list(); efeito <- list(); apa <- character(); avisos <- character()
  p_par <- numeric(); rot_par <- character()
  usados <- rep(FALSE, nrow(base))
  principal <- if ("pearson" %in% metodos) "pearson" else metodos[1]
  for (pr in pares) {
    ok <- stats::complete.cases(base[, pr])
    usados <- usados | ok
    x <- base[[pr[1]]][ok]; y <- base[[pr[2]]][ok]; n <- length(x)
    l <- data.frame(var1 = pr[1], var2 = pr[2], n = n, gl = n - 2)
    if (n < 4 || stats::sd(x) == 0 || stats::sd(y) == 0) {
      avisos <- c(avisos, paste0(pr[1], " × ", pr[2], ": n < 4 ou variável constante; correlação não calculada."))
      linhas[[length(linhas) + 1]] <- l; next
    }
    res <- list()
    if ("pearson" %in% metodos) {
      ct <- stats::cor.test(x, y, conf.level = nivel)
      res$pearson <- c(est = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2], p = ct$p.value)
      l <- cbind(l, r = res$pearson[["est"]], r_ic_inf = res$pearson[["lo"]], r_ic_sup = res$pearson[["hi"]],
                 r_p = res$pearson[["p"]], r2 = res$pearson[["est"]]^2)
    }
    if ("spearman" %in% metodos) {
      ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = NULL))
      ic <- .ee_assoc_ic_fisher(unname(ct$estimate), n, nivel)
      res$spearman <- c(est = unname(ct$estimate), lo = ic[1], hi = ic[2], p = ct$p.value)
      l <- cbind(l, rho = res$spearman[["est"]], rho_ic_inf = ic[1], rho_ic_sup = ic[2], rho_p = ct$p.value)
    }
    if ("kendall" %in% metodos) {
      ct <- suppressWarnings(stats::cor.test(x, y, method = "kendall", exact = NULL))
      ic <- .ee_assoc_ic_kendall(x, y, unname(ct$estimate), nivel)
      res$kendall <- c(est = unname(ct$estimate), lo = ic[1], hi = ic[2], p = ct$p.value)
      l <- cbind(l, tau = res$kendall[["est"]], tau_ic_inf = ic[1], tau_ic_sup = ic[2], tau_p = ct$p.value)
    }
    divergente <- FALSE
    if (all(c("pearson", "spearman") %in% metodos)) {
      d_rr <- res$pearson[["est"]] - res$spearman[["est"]]
      l$dif_r_rho <- d_rr
      decisao_dif <- (res$pearson[["p"]] < 1 - nivel) != (res$spearman[["p"]] < 1 - nivel)
      divergente <- abs(d_rr) >= limiar_divergencia || decisao_dif
      if (divergente)
        avisos <- c(avisos, paste0("Pearson e Spearman divergem em ", pr[1], " × ", pr[2], " (r = ",
                                   round(res$pearson[["est"]], 3), "; ρ = ", round(res$spearman[["est"]], 3),
                                   if (decisao_dif) "; a significância muda entre os dois" else "",
                                   "): voltar ao gráfico de dispersão (extremos influentes ou relação monotônica não linear) e relatar os dois; não escolher o resultado mais conveniente (livro)."))
    }
    pp <- res[[principal]]
    l$magnitude <- .ee_assoc_regua_r(pp[["est"]], regua, idioma)
    linhas[[length(linhas) + 1]] <- l
    simb <- c(pearson = "r", spearman = "ρ", kendall = "τ")
    for (m in names(res)) {
      efeito[[length(efeito) + 1]] <- data.frame(indice = paste0(simb[[m]], " ", pr[1], " × ", pr[2]),
                                                 valor = res[[m]][["est"]], ic_inf = res[[m]][["lo"]],
                                                 ic_sup = res[[m]][["hi"]], principal = m == principal)
    }
    # frase APA do par
    txt_coef <- function(m) {
      g <- if (m == "kendall") NULL else n - 2
      cs <- .ee_assoc_casas(res[[m]][["est"]], casas)
      paste0(simb[[m]], if (!is.null(g)) paste0("(", g, ")") else "",
             " = ", .pr_num(res[[m]][["est"]], cs, FALSE, decimal), .ee_sep(decimal),
             .pr_p(res[[m]][["p"]], decimal), .ee_sep(decimal),
             .ee_ic(res[[m]][["lo"]], res[[m]][["hi"]], cs, FALSE, decimal, nivel, idioma))
    }
    sig <- pp[["p"]] < 1 - nivel
    dir <- if (pp[["est"]] > 0) .ee_assoc_tr("positivamente", "positively", idioma) else .ee_assoc_tr("negativamente", "negatively", idioma)
    frase <- if (sig) {
      .ee_assoc_tr(paste0(pr[1], " e ", pr[2], " associaram-se ", dir, ", com magnitude ", l$magnitude, " pela régua de ",
                          .ee_assoc_autor_regua(regua, idioma), ", ", txt_coef(principal)),
                   paste0(pr[1], " and ", pr[2], " were ", dir, " associated, ", l$magnitude, " in magnitude by ",
                          .ee_assoc_autor_regua(regua, idioma), ", ", txt_coef(principal)), idioma)
    } else {
      .ee_assoc_tr(paste0(pr[1], " e ", pr[2], " não se associaram de forma estatisticamente significativa, ", txt_coef(principal)),
                   paste0(pr[1], " and ", pr[2], " were not significantly associated, ", txt_coef(principal)), idioma)
    }
    if (principal == "pearson")
      frase <- paste0(frase, .ee_sep(decimal), "r² = ", .ee_assoc_pct(pp[["est"]]^2, 1, decimal),
                      .ee_assoc_tr(" de variância compartilhada", " shared variance", idioma))
    outros <- setdiff(names(res), principal)
    if (length(outros)) {
      # cada método comparado com o principal: Spearman por magnitude (|Δ| >= limiar) e
      # significância; Kendall só pela significância (τ é menor que r por construção)
      diverge_m <- function(m) {
        dec <- (res[[m]][["p"]] < 1 - nivel) != (pp[["p"]] < 1 - nivel)
        if (m == "kendall" || principal == "kendall") dec else dec || abs(res[[m]][["est"]] - pp[["est"]]) >= limiar_divergencia
      }
      nomes <- c(pearson = "Pearson", spearman = "Spearman", kendall = "Kendall")
      frase <- paste0(frase, "; ", paste(vapply(outros, function(m) {
        conv <- if (diverge_m(m)) .ee_assoc_tr("divergiu", "diverged", idioma) else .ee_assoc_tr("convergiu", "converged", idioma)
        paste0(nomes[[m]], " ", conv, ", ", txt_coef(m))
      }, ""), collapse = "; "))
      if ("kendall" %in% outros && diverge_m("kendall"))
        avisos <- c(avisos, paste0("Kendall e ", nomes[[principal]], " chegam a decisões diferentes em ", pr[1], " × ", pr[2],
                                   " (p = ", round(res$kendall[["p"]], 3), " vs ", round(pp[["p"]], 3), "): relatar os dois e não escolher o mais conveniente."))
    }
    apa <- c(apa, paste0(frase, "."))
    cs_p <- .ee_assoc_casas(pp[["est"]], casas)
    p_par <- c(p_par, pp[["p"]])
    rot_par <- c(rot_par, paste0(pr[1], " × ", pr[2], " (", simb[[principal]], " = ", .pr_num(pp[["est"]], cs_p, FALSE, decimal), .ee_sep(decimal),
                                 .pr_p(pp[["p"]], decimal), .ee_sep(decimal), .ee_ic(pp[["lo"]], pp[["hi"]], cs_p, FALSE, decimal, nivel, "pt"), ")"))
  }
  tab <- do.call(rbind, lapply(linhas, function(l) { l[setdiff(Reduce(union, lapply(linhas, names)), names(l))] <- NA; l }))
  tab <- tab[, Reduce(union, lapply(linhas, names))]
  # pressupostos: normalidade de cada variável como grau (Vol. 1, cap. 5)
  desc <- descritivas(base[usados, , drop = FALSE], vars)
  leit <- ler_normalidade(desc)
  press <- data.frame(checagem = paste0("Normalidade de ", vars), estatistica = NA_real_, gl = NA_real_,
                      p = desc$shapiro_p, leitura = leit)
  for (i in seq_along(vars)) {
    x <- base[[vars[i]]][usados]; x <- x[!is.na(x)]
    if (length(x) >= 3 && length(x) <= 5000 && stats::sd(x) > 0) press$estatistica[i] <- unname(stats::shapiro.test(x)$statistic)
  }
  press$checagem <- paste0(press$checagem, " (W de Shapiro-Wilk)")
  n_par <- tab$n
  decisoes <- c(
    paste0("Coeficientes calculados: ", paste(metodos, collapse = ", "),
           if (all(c("pearson", "spearman") %in% metodos)) "; Pearson e Spearman lado a lado, como o livro manda quando há dúvida. Nenhuma troca silenciosa: os dois ficam no resultado." else "."),
    paste0("Ausentes: ", if (ausentes == "pairwise") "par a par (pairwise, padrão do JASP 0.95.1); n por par de " else "listwise (mesmos casos em todos os pares); n de ",
           min(n_par), " a ", max(n_par), "."),
    paste0("IC ", round(nivel * 100), "%: z de Fisher (Pearson pelo cor.test; Spearman com EP = 1/√(n − 3), como o JASP)",
           if ("kendall" %in% metodos) "; Kendall pelo IC assintótico de Hollander, Wolfe e Chicken (como o JASP)" else "", "."),
    paste0("Magnitude pela ", .ee_assoc_nome_regua(regua), ". Faixas são réguas de partida, não leis (livro); interpretar no contexto e na escala do r²."),
    paste0(length(pares), " par(es) testado(s), sem correção para comparações múltiplas (como no JASP). Comparar duas correlações (\"X se associa mais a Y do que a Z\") exige teste próprio: use comparar_correlacoes() (livro)."))
  if (length(pares) >= 10)
    avisos <- c(avisos, paste0(length(pares), " pares testados: com muitos testes, algumas correlações significativas surgem por acaso; declarar quais pares eram hipótese."))
  if (!"kendall" %in% metodos && any(n_par < 30))
    avisos <- c(avisos, "Algum par tem n < 30: o livro reserva o tau de Kendall para amostras pequenas ou com muitos empates; considere metodos = c(\"pearson\", \"spearman\", \"kendall\").")
  avisos <- c(avisos, .ee_aviso_nao_sig(p_par, 1 - nivel, NULL, "associação", if (principal == "pearson") "correlacao" else NULL, decimal,
                                        rotulos = rot_par))
  apa <- c(apa, .ee_assoc_tr("Por se tratar de desenho correlacional, as associações não autorizam inferência causal.",
                             "Because the design is correlational, these associations do not support causal inference.", idioma))
  .ee_assoc_novo(analise = "Correlações", n = if (ausentes == "listwise") nrow(chk$d) else sum(usados),
                 n_excluidos = if (ausentes == "listwise") chk$n_excluidos else nrow(base) - sum(usados),
                 tabela = tab, efeito = do.call(rbind, efeito), pressupostos = press,
                 decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- ponto-bisserial
# Pearson entre a contínua e a dicotômica codificada 0/1. Mede associação; a mesma
# estrutura de dados do t independente (o t de Student equivalente sai ao lado).
correlacao_ponto_bisserial <- function(dados, continua, dicotomica, codigo_1 = NULL, nivel = .95, casas = 3,
                                       decimal = ".", idioma = "pt") {
  dados <- .ee_assoc_vazios(dados, c(continua, dicotomica))
  dd <- .ee_dados(dados, c(continua, dicotomica), numericas = continua)
  d <- dd$d
  b <- .ee_assoc_binaria(d[[dicotomica]], codigo_1, dicotomica, exige_texto = FALSE)
  y <- d[[continua]]; g <- b$y; n <- length(y)
  ct <- stats::cor.test(y, g, conf.level = nivel)
  r <- unname(ct$estimate)
  t_eq <- r * sqrt((n - 2) / (1 - r^2))
  # Médias e DP por grupo vão em tabelas_extra com a coluna n: o print da base apaga as
  # estatísticas de grupo com n < 5 (com n = 1 a média seria o escore da pessoa). O
  # código fica como texto para continuar visível na linha suprimida.
  grupos <- data.frame(categoria = c(b$outro, b$evento), codigo = c("0", "1"),
                       n = c(sum(g == 0), sum(g == 1)),
                       media = c(mean(y[g == 0]), mean(y[g == 1])), dp = c(stats::sd(y[g == 0]), stats::sd(y[g == 1])))
  pequenos <- grupos$categoria[grupos$n < .ee_n_min_print]
  avisos_pb <- if (length(pequenos))
    paste0("Grupo(s) com n < ", .ee_n_min_print, " (", paste0("'", pequenos, "', n = ", grupos$n[grupos$n < .ee_n_min_print], collapse = "; "),
           "): média e DP desse(s) grupo(s) ficam fora do print, e a frase APA não traz média nem DP por grupo, por privacidade (com n tão pequeno a média se aproxima do escore individual); os valores ficam no objeto. ",
           "Com grupo tão pequeno, r_pb e o t equivalente são frágeis e, junto com média e DP do outro grupo, permitem recalcular a média do grupo pequeno: não publicar a análise assim (agrupar categorias ou excluir o grupo).") else character()
  tab <- data.frame(n = n, gl = n - 2, r_pb = r, ic_inf = ct$conf.int[1], ic_sup = ct$conf.int[2],
                    p = ct$p.value, r2 = r^2, t_student_equivalente = t_eq)
  dir <- if (r > 0) .ee_assoc_tr("maiores", "higher", idioma) else .ee_assoc_tr("menores", "lower", idioma)
  apa <- paste0(.ee_assoc_tr(paste0("A correlação ponto-bisserial entre ", continua, " e ", dicotomica, " (1 = ", b$evento, ") foi "),
                             paste0("The point-biserial correlation between ", continua, " and ", dicotomica, " (1 = ", b$evento, ") was "), idioma),
                "r_pb(", n - 2, ") = ",
                .ee_assoc_n(r, casas, FALSE, decimal), .ee_sep(decimal), .pr_p(ct$p.value, decimal), .ee_sep(decimal),
                .ee_ic(ct$conf.int[1], ct$conf.int[2], .ee_assoc_casas(r, casas), FALSE, decimal, nivel, idioma), .ee_sep(decimal),
                "r² = ", .ee_assoc_pct(r^2, 1, decimal),
                if (ct$p.value < 1 - nivel) .ee_assoc_tr(paste0("; o grupo '", b$evento, "' teve escores ", dir, " que o grupo '", b$outro, "'."),
                                                         paste0("; group '", b$evento, "' scored ", dir, " than group '", b$outro, "'."), idioma) else ".")
  .ee_assoc_novo(analise = "Correlação ponto-bisserial", n = n, n_excluidos = dd$n_excluidos,
                 tabela = tab,
                 efeito = data.frame(indice = "r ponto-bisserial", valor = r, ic_inf = ct$conf.int[1], ic_sup = ct$conf.int[2], principal = TRUE),
                 decisoes = c(b$frase, "O sinal de r_pb depende de quem recebeu o código 1: conferir a codificação antes de ler a direção.",
                              "t equivalente = r·√((n − 2)/(1 − r²)), igual ao t de Student (variâncias iguais) da comparação das médias; r_pb pode servir de tamanho de efeito desse t, em escala de correlação.",
                              "Dicotomia verdadeira (ter ou não ter). Se a dicotomia for artificial (corte numa contínua), o coeficiente adequado é o bisserial, não este."),
                 avisos = c(avisos_pb, .ee_aviso_nao_sig(ct$p.value, 1 - nivel, paste0("r_pb = ", .ee_assoc_n(r, casas, FALSE, decimal), .ee_sep(decimal),
                                                                         .ee_ic(ct$conf.int[1], ct$conf.int[2], .ee_assoc_casas(r, casas), FALSE, decimal, nivel, "pt")),
                                                         "associação", "independente", decimal)),
                 apa = apa, tabelas_extra = list(`Codificação e grupos` = grupos))
}

# ---------------------------------------------------------------- correlação parcial
# Estimativa e p pelo ppcor::pcor.test (como o JASP: t com gl = n − 2 − k). O JASP
# 0.95.1 não dá IC para a parcial; aqui o IC de Pearson sai pela z de Fisher com
# EP = 1/√(n − 3 − k). Para Spearman parcial o IC fica em branco.
correlacao_parcial <- function(dados, x, y, controles, metodo = c("pearson", "spearman"), nivel = .95, casas = 3,
                               decimal = ".", idioma = "pt") {
  metodo <- match.arg(metodo)
  vars <- c(x, y, controles)
  dados <- .ee_assoc_vazios(dados, vars)
  dd <- .ee_dados(dados, vars, numericas = vars)
  d <- dd$d; n <- nrow(d); k <- length(controles)
  .pr_exige("ppcor")
  pc <- ppcor::pcor.test(d[[x]], d[[y]], d[, controles, drop = FALSE], method = metodo)
  biv <- if (metodo == "pearson") stats::cor.test(d[[x]], d[[y]], conf.level = nivel) else
    suppressWarnings(stats::cor.test(d[[x]], d[[y]], method = "spearman", exact = NULL))
  rb <- unname(biv$estimate); rp <- pc$estimate
  ic_b <- if (metodo == "pearson") biv$conf.int[1:2] else .ee_assoc_ic_fisher(rb, n, nivel)
  ic_p <- if (metodo == "pearson") .ee_assoc_ic_fisher(rp, n, nivel, k) else c(NA, NA)
  tab <- data.frame(correlacao = c("bivariada (ordem zero)", paste0("parcial (controlando ", paste(controles, collapse = ", "), ")")),
                    coeficiente = c(rb, rp), ic_inf = c(ic_b[1], ic_p[1]), ic_sup = c(ic_b[2], ic_p[2]),
                    gl = c(n - 2, n - 2 - k), p = c(biv$p.value, pc$p.value), r2 = c(rb^2, rp^2))
  queda <- 100 * (rb^2 - rp^2)
  simb <- if (metodo == "pearson") "r" else "ρ"
  apa <- paste0(.ee_assoc_tr(paste0("A correlação entre ", x, " e ", y, " foi "), paste0("The correlation between ", x, " and ", y, " was "), idioma),
                simb, "(", n - 2, ") = ", .ee_assoc_n(rb, casas, FALSE, decimal), .ee_sep(decimal), .pr_p(biv$p.value, decimal),
                .ee_assoc_tr(paste0("; controlando ", paste(controles, collapse = ", "), ", a correlação parcial foi "),
                             paste0("; controlling for ", paste(controles, collapse = ", "), ", the partial correlation was "), idioma),
                simb, "(", n - 2 - k, ") = ", .ee_assoc_n(rp, casas, FALSE, decimal), .ee_sep(decimal), .pr_p(pc$p.value, decimal),
                if (!is.na(ic_p[1])) paste0(.ee_sep(decimal), .ee_ic(ic_p[1], ic_p[2], .ee_assoc_casas(rp, casas), FALSE, decimal, nivel, idioma)) else "",
                .ee_assoc_tr(paste0(" (variância compartilhada de ", .ee_assoc_pct(rb^2, 1, decimal), " para ", .ee_assoc_pct(rp^2, 1, decimal), ")."),
                             paste0(" (shared variance from ", .ee_assoc_pct(rb^2, 1, decimal), " to ", .ee_assoc_pct(rp^2, 1, decimal), ")."), idioma))
  .ee_assoc_novo(analise = paste0("Correlação parcial (", metodo, ")"), n = n, n_excluidos = dd$n_excluidos,
                 tabela = tab,
                 efeito = data.frame(indice = c(paste0(simb, " bivariada"), paste0(simb, " parcial")), valor = c(rb, rp),
                                     ic_inf = c(ic_b[1], ic_p[1]), ic_sup = c(ic_b[2], ic_p[2]), principal = c(FALSE, TRUE)),
                 decisoes = c(paste0("Parcial = correlação entre os resíduos de ", x, " e de ", y, " depois de regredidos em ", paste(controles, collapse = ", "),
                                     "; gl = n − 2 − k (k = ", k, "). Estimativa e p pelo ppcor (como o JASP)."),
                              if (metodo == "pearson") "IC da parcial pela z de Fisher com EP = 1/√(n − 3 − k); o JASP 0.95.1 não exibe IC para a parcial." else
                                "IC não calculado para Spearman parcial (sem método consolidado no JASP).",
                              paste0("Variância compartilhada caiu ", round(queda, 1), " pontos percentuais com o controle."),
                              "A parcial não distingue confundidor de mediador: o papel da variável controlada é decisão teórica, não resultado da conta."),
                 avisos = .ee_aviso_nao_sig(pc$p.value, 1 - nivel,
                                            if (!is.na(ic_p[1])) paste0(simb, " parcial = ", .ee_assoc_n(rp, casas, FALSE, decimal), .ee_sep(decimal),
                                                                        .ee_ic(ic_p[1], ic_p[2], .ee_assoc_casas(rp, casas), FALSE, decimal, nivel, "pt")),
                                            "associação parcial", NULL, decimal),
                 apa = apa)
}

# ---------------------------------------------------------------- comparar correlações
# Independentes (dois grupos): z de Fisher da diferença e IC de Zou (2007) para r1 − r2.
# Dependentes sobrepostas (mesma amostra, uma variável em comum: r_jk vs r_jh, com
# r_kh entre as não comuns): z de Steiger (1980) com r médio e z de Meng, Rosenthal e
# Rubin (1992), os mesmos do pacote cocor. Para o caso dependente informe n1 = n
# (tamanho da amostra) e r12 = correlação entre as duas variáveis não comuns.
comparar_correlacoes <- function(r1, n1, r2, n2 = NULL, r12 = NULL, nivel = .95, casas = 3, decimal = ".", idioma = "pt") {
  if (any(abs(c(r1, r2, r12)) >= 1)) stop("correlações precisam estar entre -1 e 1")
  zc <- stats::qnorm(1 - (1 - nivel) / 2)
  if (is.null(r12)) {
    if (is.null(n2)) stop("para correlações independentes informe n2; para dependentes informe r12 (correlação entre as variáveis não comuns)")
    z <- (atanh(r1) - atanh(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    p <- 2 * stats::pnorm(-abs(z))
    ic1 <- .ee_assoc_ic_fisher(r1, n1, nivel); ic2 <- .ee_assoc_ic_fisher(r2, n2, nivel)
    lo <- r1 - r2 - sqrt((r1 - ic1[1])^2 + (ic2[2] - r2)^2)
    hi <- r1 - r2 + sqrt((ic1[2] - r1)^2 + (r2 - ic2[1])^2)
    tab <- data.frame(cenario = "independentes", r1 = r1, n1 = n1, r2 = r2, n2 = n2, diferenca = r1 - r2,
                      ic_inf = lo, ic_sup = hi, teste = "z de Fisher", z = z, p = p)
    efe <- data.frame(indice = "r1 − r2", valor = r1 - r2, ic_inf = lo, ic_sup = hi, principal = TRUE)
    dec <- c("Cenário: correlações de amostras independentes (dois grupos). Teste: transformação r → z de Fisher e z da diferença (livro).",
             "IC da diferença r1 − r2 pelo método de Zou (2007), que combina os IC de Fisher de cada r.",
             "Se as duas correlações vêm da mesma amostra, este teste está errado: informe r12 para o teste de correlações dependentes.")
    n_tot <- n1 + n2
  } else {
    n <- n1
    zjk <- atanh(r1); zjh <- atanh(r2)
    rm <- (r1 + r2) / 2
    psi <- r12 * (1 - 2 * rm^2) - .5 * rm^2 * (1 - 2 * rm^2 - r12^2)
    cc <- psi / (1 - rm^2)^2
    z_s <- (zjk - zjh) * sqrt(n - 3) / sqrt(2 - 2 * cc)
    r2m <- (r1^2 + r2^2) / 2; f <- min(1, (1 - r12) / (2 * (1 - r2m))); h <- (1 - f * r2m) / (1 - r2m)
    z_m <- (zjk - zjh) * sqrt((n - 3) / (2 * (1 - r12) * h))
    tab <- data.frame(cenario = "dependentes sobrepostas", r1 = r1, r2 = r2, r12 = r12, n = n, diferenca = r1 - r2,
                      teste = c("z de Steiger (1980)", "z de Meng, Rosenthal e Rubin (1992)"),
                      z = c(z_s, z_m), p = 2 * stats::pnorm(-abs(c(z_s, z_m))))
    z <- z_s; p <- tab$p[1]
    efe <- data.frame(indice = "r1 − r2", valor = r1 - r2, ic_inf = NA, ic_sup = NA, principal = TRUE)
    dec <- c("Cenário: correlações dependentes sobrepostas (mesma amostra, uma variável em comum). O teste leva em conta a correlação entre as próprias correlações (r12), como o livro pede.",
             "Leitura pelo z de Steiger (1980); Meng et al. (1992) ao lado (os dois coincidem com o pacote cocor). Decidir por |z| e p; o sinal só diz qual r é maior.",
             "Conferir que n é o N efetivo dos pares usados nas três correlações (listwise); com outro N o resultado muda.")
    n_tot <- n
  }
  rt1 <- .ee_assoc_n(r1, casas, FALSE, decimal); rt2 <- .ee_assoc_n(r2, casas, FALSE, decimal)
  apa <- paste0(.ee_assoc_tr(paste0("A diferença entre as correlações (r = ", rt1, " e r = ", rt2, ") "),
                             paste0("The difference between the correlations (r = ", rt1, " and r = ", rt2, ") "), idioma),
                if (p < 1 - nivel) .ee_assoc_tr("foi estatisticamente significativa, ", "was statistically significant, ", idioma) else
                  .ee_assoc_tr("não foi estatisticamente significativa, ", "was not statistically significant, ", idioma),
                "z = ", .pr_num(z, 2, TRUE, decimal), .ee_sep(decimal), .pr_p(p, decimal), ".")
  av <- .ee_aviso_nao_sig(p, 1 - nivel,
                          if (!is.na(efe$ic_inf[1])) paste0("r1 - r2 = ", .ee_assoc_n(r1 - r2, casas, FALSE, decimal), .ee_sep(decimal),
                                                            .ee_ic(efe$ic_inf[1], efe$ic_sup[1], .ee_assoc_casas(r1 - r2, casas), FALSE, decimal, nivel, "pt")),
                          "diferença entre as correlações", NULL, decimal)
  .ee_assoc_novo(analise = "Comparação de duas correlações", n = n_tot, tabela = tab, efeito = efe, decisoes = dec, avisos = av, apa = apa)
}

# ================================================================ regressão linear
# Grupo gravado como número (1/2/3, formato comum em banco vindo do SPSS ou de planilha)
# entra como contínuo, com um só coeficiente linear. Regra do aviso: preditor numérico
# com 3 a 7 valores distintos, todos inteiros. Com 2 valores (0/1) o código numérico é o
# que o livro manda para a dicotômica; com decimais ou mais de 7 valores, a
# variável é tratada como escala. `dica` diz o que fazer em cada função.
.ee_assoc_cod_k <- c(3, 7)
.ee_assoc_aviso_codigo <- function(x, v, dica) {
  if (!is.numeric(x) || is.factor(x)) return(character())
  u <- sort(unique(x[!is.na(x)])); k <- length(u)
  if (k < .ee_assoc_cod_k[1] || k > .ee_assoc_cod_k[2] || any(abs(u - round(u)) > 1e-8)) return(character())
  paste0(v, " é numérica com só ", k, " valores inteiros distintos (", paste(u, collapse = ", "),
         ") e entrou como contínua (um coeficiente linear, 1 gl). Se for grupo ou categoria sem ordem e distância, ", dica,
         ": códigos numéricos diretos para 3 ou mais categorias impõem ordem e distância que elas não têm (livro). ",
         "Se for escala (item Likert, contagem), o aviso pode ser ignorado.")
}

# Preditores categóricos (fator, texto, lógico e os numéricos listados em `fatores`) viram
# indicadoras k − 1 contra a referência (primeiro nível ou a indicada em `referencias`,
# lista nomeada); nunca entram como códigos numéricos (livro). Numérico
# fora de `fatores` com 3 a 7 valores inteiros gera aviso (.ee_assoc_aviso_codigo).
.ee_assoc_preparar <- function(d, preds, referencias = NULL, fatores = NULL) {
  if (length(fatores)) {
    fatores <- as.character(fatores)
    fora <- setdiff(fatores, preds)
    if (length(fora)) stop("fatores cita variável que não é preditor: ", paste(fora, collapse = ", "))
  }
  if (!is.null(referencias)) {
    referencias <- as.list(referencias)  # aceita list(g = "C") ou c(g = "C")
    if (is.null(names(referencias)) || any(!nzchar(names(referencias))))
      stop("referencias precisa ser nomeada: list(fator = \"nível\") ou c(fator = \"nível\")")
    fora <- setdiff(names(referencias), preds)
    if (length(fora)) stop("referencias cita variável que não é preditor: ", paste(fora, collapse = ", "))
  }
  info <- list(); avisos <- character()
  for (v in preds) {
    x <- d[[v]]
    if (is.numeric(x) && !is.factor(x) && !v %in% fatores) {
      if (stats::sd(x) == 0) stop("o preditor ", v, " não varia nos casos analisados (livro)")
      if (!is.null(referencias[[v]]))
        stop("referencias pede categoria de referência para ", v, ", que é numérica e entra como contínua; se for grupo, declare fatores = \"", v, "\"")
      avisos <- c(avisos, .ee_assoc_aviso_codigo(x, v, paste0("declare fatores = \"", v, "\" (indicadoras k − 1, menor código como referência)")))
      info[[v]] <- list(tipo = "numerica"); next
    }
    f <- droplevels(as.factor(if (is.logical(x)) as.character(x) else x))
    if (nlevels(f) < 2) stop("o preditor ", v, " tem só uma categoria com dados")
    cont <- table(f)
    if (any(cont < 2))
      stop("o preditor categórico ", v, " tem categoria com menos de 2 casos (", paste0("'", names(cont)[cont < 2], "', n = ", cont[cont < 2], collapse = "; "),
           "): com 1 caso a indicadora reproduz o valor da pessoa e o resíduo dela fica indefinido. Agrupe categorias ou retire o caso antes de rodar.", call. = FALSE)
    ref <- if (!is.null(referencias[[v]])) as.character(referencias[[v]]) else levels(f)[1]
    if (!ref %in% levels(f)) stop("referência '", ref, "' não é categoria de ", v, " (", paste(levels(f), collapse = ", "), ")")
    f <- stats::relevel(f, ref)
    stats::contrasts(f) <- stats::contr.treatment(levels(f))
    d[[v]] <- f
    info[[v]] <- list(tipo = "fator", ref = ref, niveis = levels(f), declarado = v %in% fatores)
  }
  list(d = d, info = info, avisos = avisos)
}

# Frase de $decisoes para os preditores categóricos (com a origem: fator/texto ou `fatores`).
.ee_assoc_txt_fatores <- function(info) {
  fat <- names(Filter(function(i) i$tipo == "fator", info))
  if (!length(fat)) return(character())
  paste(vapply(fat, function(v) paste0(v, " (referência = ", info[[v]]$ref,
                                       if (isTRUE(info[[v]]$declarado)) "; numérica declarada em fatores" else "", ")"), ""), collapse = "; ")
}

.ee_assoc_nomes_coef <- function(X, preds, info) {
  asg <- attr(X, "assign")
  vapply(seq_len(ncol(X)), function(j) {
    if (asg[j] == 0) return("(Intercepto)")
    v <- preds[asg[j]]
    if (info[[v]]$tipo == "numerica") return(v)
    k <- sum(asg[seq_len(j)] == asg[j])
    paste0(v, ": ", info[[v]]$niveis[k + 1], " vs ", info[[v]]$ref)
  }, character(1))
}

# VIF por termo. Termo com mais de 1 gl (fator com 3+ níveis): GVIF^(1/gl), na mesma
# escala do VIF (Fox & Monette). O JASP exibe GVIF^(1/(2gl)) para todos os termos quando
# há fator com 3+ níveis; aqui não (ver references/associacao.md).
.ee_assoc_vif <- function(fit) {
  nt <- length(attr(stats::terms(fit), "term.labels"))
  if (nt < 2) return(data.frame(vif = rep(1, nt), gl = rep(1, nt)))
  .pr_exige("car")
  v <- car::vif(fit)
  if (is.matrix(v)) data.frame(vif = v[, 1]^(1 / v[, 2]), gl = v[, 2]) else data.frame(vif = unname(v), gl = 1)
}

.ee_assoc_coef_lm <- function(fit, d, vd, preds, info, nivel, ep_robusto = FALSE) {
  if (any(is.na(stats::coef(fit)))) stop("coeficiente não estimável: preditores perfeitamente colineares ou categoria sem casos")
  X <- stats::model.matrix(fit); asg <- attr(X, "assign")
  sm <- summary(fit)$coefficients
  ci <- stats::confint(fit, level = nivel)
  beta <- stats::coef(fit) * apply(X, 2, stats::sd) / stats::sd(d[[vd]])
  beta[asg == 0] <- NA
  for (i in seq_along(preds)) if (info[[preds[i]]]$tipo == "fator" && length(info[[preds[i]]]$niveis) > 2) beta[asg == i] <- NA
  vf <- .ee_assoc_vif(fit)
  vif_col <- c(NA, vf$vif[asg[asg > 0]])
  tab <- data.frame(termo = .ee_assoc_nomes_coef(X, preds, info), B = sm[, 1], EP = sm[, 2],
                    ic_inf = ci[, 1], ic_sup = ci[, 2], beta = unname(beta), t = sm[, 3], p = sm[, 4],
                    tolerancia = 1 / vif_col, vif = vif_col, row.names = NULL)
  if (length(preds) < 2) tab$tolerancia <- tab$vif <- NULL
  if (ep_robusto) {
    .pr_exige("sandwich")
    ep <- sqrt(diag(sandwich::vcovHC(fit, type = "HC3")))
    tq <- stats::qt(1 - (1 - nivel) / 2, stats::df.residual(fit))
    tab$EP_HC3 <- ep; tab$ic_inf_HC3 <- tab$B - tq * ep; tab$ic_sup_HC3 <- tab$B + tq * ep
    tab$p_HC3 <- 2 * stats::pt(-abs(tab$B / ep), stats::df.residual(fit))
  }
  list(tabela = tab, vif = vf)
}

.ee_assoc_leitura_norm <- function(x, n) {
  a <- .ee_assimetria(x); k <- .ee_curtose(x); m <- max(abs(c(a, k)), na.rm = TRUE)
  faixa <- if (!is.finite(m)) "indeterminado" else if (m < 1) "afastamento leve" else if (m <= 2) "afastamento moderado" else "afastamento que preocupa"
  paste0(faixa, " (assimetria ", round(a, 2), ", curtose ", round(k, 2), ")")
}

# Diagnósticos do modelo linear (livro; curso, T3 seção 4).
.ee_assoc_diag_lm <- function(fit, tipo_residuo = "studentizado") {
  n <- stats::nobs(fit); e <- stats::residuals(fit)
  X <- stats::model.matrix(fit)[, -1, drop = FALSE]
  k <- ncol(X)
  vf <- .ee_assoc_vif(fit)
  linhas <- list(); avisos <- character()
  add <- function(checagem, estatistica = NA, gl = NA, p = NA, leitura) linhas[[length(linhas) + 1]] <<-
    data.frame(checagem = checagem, estatistica = estatistica, gl = gl, p = p, leitura = leitura)
  termos <- attr(stats::terms(fit), "term.labels")
  if (length(termos) >= 2) {
    i <- which.max(vf$vif); v <- vf$vif[i]
    leit <- if (v > 10) "acima de 10 (tolerância < .10): colinearidade forte; coeficientes instáveis, EP inflados, significâncias que vão e vêm" else
      if (v > 5) "entre 5 e 10 (tolerância entre .10 e .20): alerta de colinearidade; ler os coeficientes com cautela" else
        "abaixo das réguas de alerta (VIF > 5 ou > 10; tolerância < .20 ou < .10)"
    add("Multicolinearidade (maior VIF; tolerância = 1/VIF)", v, NA, NA,
        paste0("VIF = ", sprintf("%.2f", v), " (tolerância ", sprintf("%.3f", 1 / v), ") em ", gsub("`", "", termos[i]), ": ", leit,
               ". Réguas de partida, não vereditos (livro; O'Brien, 2007)."))
    if (v > 5) avisos <- c(avisos, paste0("Multicolinearidade: VIF máximo = ", round(v, 2), " (", gsub("`", "", termos[i]), ")."))
  }
  # Durbin-Watson: estatística como o car e p bilateral do lmtest::dwtest (JASP 0.95.1)
  dw <- sum(diff(e)^2) / sum(e^2); ac <- sum(e[-1] * e[-n]) / sum(e^2)
  .pr_exige("lmtest")
  dwp <- tryCatch(lmtest::dwtest(fit, alternative = "two.sided")$p.value, error = function(err) NA_real_)
  add("Independência dos resíduos (Durbin-Watson)", dw, NA, dwp,
      paste0("DW = ", sprintf("%.2f", dw), " (autocorrelação lag 1 = ", sprintf("%.3f", ac), "); perto de 2 = sem autocorrelação (livro). ",
             "Só informa se a ordem das linhas tem sentido (tempo, ordem de coleta); em corte transversal sem ordem natural tem pouco a dizer."))
  # normalidade dos resíduos brutos (o livro reporta o SW dos resíduos: W = .965 no Exam Anxiety)
  if (n >= 3 && n <= 5000) {
    sw <- stats::shapiro.test(e)
    add("Normalidade dos resíduos (W de Shapiro-Wilk)", unname(sw$statistic), NA, sw$p.value,
        paste0("Resíduos, não o desfecho: ", .ee_assoc_leitura_norm(e, n),
               ". Avaliação de grau: olhar o Q-Q dos resíduos", if (n >= 100) "; com este n a inferência tolera afastamentos modestos e o teste acusa desvios pequenos" else "", "."))
  }
  bp <- lmtest::bptest(fit)
  add("Homocedasticidade (Breusch-Pagan studentizado)", unname(bp$statistic), unname(bp$parameter), bp$p.value,
      if (bp$p.value < .05) "Indício de heterocedasticidade: conferir o funil no gráfico resíduos × previstos. Coeficientes seguem não viesados; EP e IC ficam distorcidos. Alternativa: ep_robusto = TRUE (EP HC3)." else
        "Sem indício formal de heterocedasticidade; confirmar no gráfico resíduos × previstos (padrão de funil).")
  if (any(bp$p.value < .05)) avisos <- c(avisos, "Breusch-Pagan sugere heterocedasticidade: ver o gráfico resíduos × previstos e considerar ep_robusto = TRUE.")
  # resíduos padronizados
  z <- if (tipo_residuo == "zresid") e / stats::sigma(fit) else stats::rstandard(fit)
  rot <- if (tipo_residuo == "zresid") "ZRESID = resíduo/σ, como o SPSS" else "studentizado interno, rstandard(), como o JASP 0.95.1"
  # Caso com alavanca 1 tem resíduo padronizado e Cook indefinidos (NaN): fica fora das
  # contagens (na.rm) em vez de derrubar a função.
  n3 <- sum(abs(z) > 3, na.rm = TRUE); p25 <- mean(abs(z) > 2.5, na.rm = TRUE); p2 <- mean(abs(z) > 2, na.rm = TRUE)
  add("Resíduos padronizados |z| > 3 (casos)", n3, NA, NA,
      paste0(n3, " caso(s) (", rot, "); regra do curso: nenhum caso acima de 3."))
  add("Resíduos padronizados |z| > 2,5 (proporção)", p25, NA, NA,
      paste0(round(100 * p25, 1), "% dos casos; esperado até 1%", if (isTRUE(p25 > .01)) ": acima do esperado, rever a especificação do modelo" else "", "."))
  add("Resíduos padronizados |z| > 2 (proporção)", p2, NA, NA,
      paste0(round(100 * p2, 1), "% dos casos; esperado em torno de 5% sob normalidade", if (isTRUE(p2 > .05)) ": acima, tratar como alerta" else "", "."))
  ck <- stats::cooks.distance(fit)
  add("Distância de Cook (máximo)", max(ck, na.rm = TRUE), NA, NA,
      paste0(sum(ck > 1, na.rm = TRUE), " caso(s) com D > 1 (influente; Cook, 1977). Alternativa mais sensível, 4/n = ", round(4 / n, 3), ": ",
             sum(ck > 4 / n, na.rm = TRUE), " caso(s), só como alerta. Caso influente não sai automaticamente: corrigir se for erro; se real, decidir por critério declarado."))
  # v2.2: DFFITS e DFBETAS (Belsley, Kuh & Welsch), usados nos scripts de regressão do curso (dffits, dfbeta)
  p_par <- length(stats::coef(fit)); dff <- stats::dffits(fit); dfb <- stats::dfbetas(fit)
  c_dff <- 2 * sqrt(p_par / n); c_dfb <- 2 / sqrt(n)
  add("DFFITS (máximo |DFFITS|)", max(abs(dff), na.rm = TRUE), NA, NA,
      paste0(sum(abs(dff) > c_dff, na.rm = TRUE), " caso(s) acima de 2*sqrt(p/n) = ", round(c_dff, 3), "; régua de alerta, não de exclusão."))
  add("DFBETAS (máximo |DFBETAS|)", max(abs(dfb), na.rm = TRUE), NA, NA,
      paste0(sum(apply(abs(dfb) > c_dfb, 1, any), na.rm = TRUE), " caso(s) com algum |DFBETAS| > 2/sqrt(n) = ", round(c_dfb, 3), "; mostra qual coeficiente o caso desloca."))
  if (k >= 1) {
    D2 <- stats::mahalanobis(X, colMeans(X), stats::cov(X))
    crit <- stats::qchisq(.999, k)
    add("Distância de Mahalanobis (máximo)", max(D2), k, stats::pchisq(max(D2), k, lower.tail = FALSE),
        paste0(sum(D2 > crit), " caso(s) acima do crítico χ²(", k, ") a p < .001 = ", round(crit, 2),
               ". A tabela do curso (11, 15 e 25 para 2, 3 e 5 preditores) aproxima esse critério com folga (críticos exatos 13.82, 16.27 e 20.52)."))
  } else D2 <- numeric()
  if (n3 > 0) avisos <- c(avisos, paste0(n3, " caso(s) com resíduo padronizado |z| > 3 (índices em $extra_diagnostico)."))
  if (sum(ck > 1, na.rm = TRUE) > 0) avisos <- c(avisos, paste0(sum(ck > 1, na.rm = TRUE), " caso(s) com distância de Cook > 1."))
  list(pressupostos = do.call(rbind, linhas), avisos = avisos,
       indices = list(residuo_maior_3 = names(z)[which(abs(z) > 3)], cook_maior_1 = names(ck)[which(ck > 1)],
                      mahalanobis_acima = names(D2)[D2 > if (k >= 1) stats::qchisq(.999, k) else Inf]),
       bp = bp$p.value)
}

.ee_assoc_resumo_lm <- function(fit, y) {
  s <- summary(fit); n <- stats::nobs(fit)
  r2 <- s$r.squared; k <- length(stats::coef(fit)) - 1
  fst <- if (k > 0) s$fstatistic else c(value = NA, numdf = 0, dendf = n - 1)
  data.frame(n = n, R = sqrt(r2), R2 = r2, R2_ajustado = s$adj.r.squared, RMSE = stats::sigma(fit),
             RMSE_nulo = stats::sd(y), F = unname(fst[1]), gl1 = unname(fst[2]), gl2 = unname(fst[3]),
             p = if (k > 0) stats::pf(fst[1], fst[2], fst[3], lower.tail = FALSE) else NA, f2 = r2 / (1 - r2), row.names = NULL)
}

.ee_assoc_frase_coef <- function(tab, gl, nivel, casas, decimal, idioma, robusto = FALSE) {
  tab <- tab[tab$termo != "(Intercepto)", , drop = FALSE]
  vapply(seq_len(nrow(tab)), function(i) {
    l <- tab[i, ]
    lo <- if (robusto) l$ic_inf_HC3 else l$ic_inf; hi <- if (robusto) l$ic_sup_HC3 else l$ic_sup
    p <- if (robusto) l$p_HC3 else l$p
    cb <- .ee_assoc_casas(l$B, casas)
    # β com zero à esquerda: pode passar de 1 (livro)
    paste0(l$termo, ": B = ", .pr_num(l$B, cb, TRUE, decimal), .ee_sep(decimal),
           .ee_ic(lo, hi, cb, TRUE, decimal, nivel, idioma),
           if (!is.na(l$beta)) paste0(.ee_sep(decimal), "β = ", .ee_assoc_n(l$beta, casas, TRUE, decimal)) else "",
           .ee_sep(decimal), .ee_estat("t", if (robusto) l$B / l$EP_HC3 else l$t, gl, decimal), .ee_sep(decimal), .pr_p(p, decimal))
  }, character(1))
}

.ee_assoc_causal <- function(idioma) .ee_assoc_tr("Os coeficientes descrevem associações condicionais; sem desenho experimental, não autorizam leitura causal.",
                                                   "Coefficients describe conditional associations; without an experimental design they do not support causal claims.", idioma)

# Nota de $decisoes para os índices do modelo que saem sem IC (regra 6 da skill: o relato
# diz que o IC não foi calculado).
.ee_assoc_nota_sem_ic <- function(pt, en, idioma)
  .ee_assoc_tr(paste0(pt, " sem IC: não calculado (dizer isso no relato, ao lado do valor)."),
               paste0(en, " without CI: not computed (state this in the report, next to the value)."), idioma)

regressao_linear <- function(dados, vd, preditores, referencias = NULL, fatores = NULL, nivel = .95,
                             tipo_residuo = c("studentizado", "zresid"), ep_robusto = FALSE, casas = 3,
                             decimal = ".", idioma = "pt") {
  tipo_residuo <- match.arg(tipo_residuo)
  vars <- c(vd, preditores)
  dados <- .ee_assoc_vazios(dados, vars)
  dd <- .ee_dados(dados, vars, numericas = vd)
  y <- dd$d[[vd]]
  if (length(unique(y)) == 2) stop("o desfecho ", vd, " é dicotômico: use regressao_logistica() (livro)")
  prep <- .ee_assoc_preparar(dd$d, preditores, referencias, fatores)
  d <- prep$d
  fit <- stats::lm(.ee_assoc_formula(vd, preditores), data = d)
  cf <- .ee_assoc_coef_lm(fit, d, vd, preditores, prep$info, nivel, ep_robusto)
  mod <- .ee_assoc_resumo_lm(fit, y)
  dg <- .ee_assoc_diag_lm(fit, tipo_residuo)
  efe <- data.frame(indice = c("R²", "R² ajustado", "f²"), valor = c(mod$R2, mod$R2_ajustado, mod$f2),
                    ic_inf = NA, ic_sup = NA, principal = c(FALSE, TRUE, TRUE))
  txt_fat <- .ee_assoc_txt_fatores(prep$info)
  decisoes <- c(paste0("Método Enter (todos os preditores juntos); modelo é decisão de pesquisa, não do algoritmo (livro)."),
                paste0("B com IC ", round(nivel * 100), "% (t); β = B·DP(x)/DP(y), como o JASP; comparar preditores pelos β, nunca pelos B brutos."),
                "Tamanho de efeito do modelo: R² ajustado (o que se relata) e f² = R²/(1 − R²) (Cohen, 1988; .02/.15/.35 como ponto de partida).",
                .ee_assoc_nota_sem_ic("R² ajustado e f²", "Adjusted R² and f²", idioma),
                paste0("Resíduos padronizados: ", if (tipo_residuo == "zresid") "ZRESID (resíduo/σ, SPSS e curso)" else "studentizados internos (rstandard, como o JASP 0.95.1); tipo_residuo = \"zresid\" reproduz o SPSS", "."))
  if (length(txt_fat))
    decisoes <- c(decisoes, paste0("Preditor(es) categórico(s) em indicadoras k − 1: ", txt_fat,
                                   ". β exibido só para fator de 2 níveis (equivale a recodificar 0/1 como escala, livro); para 3+ níveis, ler B."))
  if (ep_robusto) decisoes <- c(decisoes, "EP robusto HC3 (sandwich) calculado ao lado; a frase APA usa o EP HC3 porque ep_robusto = TRUE.")
  sig <- !is.na(mod$p) && mod$p < 1 - nivel
  frase_mod <- paste0(
    if (sig) .ee_assoc_tr("O modelo foi estatisticamente significativo, ", "The model was statistically significant, ", idioma) else
      .ee_assoc_tr("O modelo não foi estatisticamente significativo, ", "The model was not statistically significant, ", idioma),
    .ee_estat("F", mod$F, c(mod$gl1, mod$gl2), decimal), .ee_sep(decimal), .pr_p(mod$p, decimal),
    .ee_assoc_tr(paste0(", explicando ", .ee_assoc_pct(max(mod$R2_ajustado, 0), 1, decimal), " da variância de ", vd, " ("),
                 paste0(", explaining ", .ee_assoc_pct(max(mod$R2_ajustado, 0), 1, decimal), " of the variance in ", vd, " ("), idioma),
    .ee_assoc_tr("R² ajustado", "adjusted R²", idioma), " = ", .pr_num(mod$R2_ajustado, 3, FALSE, decimal), .ee_sep(decimal),
    "f² = ", .pr_num(mod$f2, 3, TRUE, decimal), ").")
  apa <- c(frase_mod, paste0(.ee_assoc_frase_coef(cf$tabela, mod$gl2, nivel, casas, decimal, idioma, ep_robusto), "."), .ee_assoc_causal(idioma))
  avisos <- c(prep$avisos, dg$avisos)
  if (!sig) avisos <- c(avisos, "Modelo global não significativo: em geral não se aprofunda a leitura dos coeficientes (curso, roteiro de 6 passos).")
  .ee_assoc_novo(analise = paste0("Regressão linear ", if (length(preditores) == 1) "simples" else "múltipla", " (", vd, ")"),
                 n = nrow(d), n_excluidos = dd$n_excluidos, tabela = cf$tabela, efeito = efe,
                 pressupostos = dg$pressupostos, decisoes = decisoes, avisos = avisos, apa = apa,
                 extra = list(modelo = fit, extra_diagnostico = dg$indices),
                 tabelas_extra = list(`Modelo (ANOVA e resumo)` = mod))
}

# ================================================================ regressão hierárquica
# Blocos definidos ANTES de ver os resultados (livro). Mesma amostra em todos
# os degraus (listwise em todas as variáveis). F de mudança por anova(m_anterior, m_atual),
# dois modelos de cada vez, como o JASP; anova(m1, m2, m3) com três modelos usaria o
# resíduo do maior modelo e daria outro F. f² do bloco = ΔR²/(1 − R² do degrau).
regressao_hierarquica <- function(dados, vd, blocos, referencias = NULL, fatores = NULL, nivel = .95,
                                  tipo_residuo = c("studentizado", "zresid"), ep_robusto = FALSE, casas = 3,
                                  decimal = ".", idioma = "pt") {
  tipo_residuo <- match.arg(tipo_residuo)
  if (!is.list(blocos) || length(blocos) < 2) stop("blocos precisa ser uma lista com pelo menos 2 vetores de preditores, na ordem da hipótese")
  todos <- unlist(blocos)
  if (anyDuplicated(todos)) stop("um preditor aparece em mais de um bloco")
  vars <- c(vd, todos)
  dados <- .ee_assoc_vazios(dados, vars)
  dd <- .ee_dados(dados, vars, numericas = vd)
  prep <- .ee_assoc_preparar(dd$d, todos, referencias, fatores)
  d <- prep$d; y <- d[[vd]]
  nomes_bl <- if (!is.null(names(blocos)) && all(nzchar(names(blocos)))) names(blocos) else
    paste0(.ee_assoc_tr("bloco ", "block ", idioma), seq_along(blocos))
  modelos <- list(stats::lm(.ee_assoc_formula(vd, character()), data = d))
  for (b in seq_along(blocos)) modelos[[b + 1]] <- stats::lm(.ee_assoc_formula(vd, unlist(blocos[seq_len(b)])), data = d)
  passos <- list()
  for (b in seq_along(modelos)) {
    m <- modelos[[b]]; s <- .ee_assoc_resumo_lm(m, y)
    l <- data.frame(modelo = paste0("M", b - 1), entrou = if (b == 1) "(nulo)" else paste(blocos[[b - 1]], collapse = ", "),
                    R = s$R, R2 = s$R2, R2_ajustado = s$R2_ajustado, RMSE = s$RMSE,
                    delta_R2 = NA, F_mudanca = NA, gl1 = NA, gl2 = NA, p = NA, f2_bloco = NA)
    if (b > 1) {
      a <- stats::anova(modelos[[b - 1]], m)
      l$delta_R2 <- s$R2 - summary(modelos[[b - 1]])$r.squared
      l$F_mudanca <- a$F[2]; l$gl1 <- a$Df[2]; l$gl2 <- a$Res.Df[2]; l$p <- a$`Pr(>F)`[2]
      l$f2_bloco <- l$delta_R2 / (1 - s$R2)
    }
    passos[[b]] <- l
  }
  passos <- do.call(rbind, passos)
  final <- modelos[[length(modelos)]]
  cf <- .ee_assoc_coef_lm(final, d, vd, todos, prep$info, nivel, ep_robusto)
  mod <- .ee_assoc_resumo_lm(final, y)
  dg <- .ee_assoc_diag_lm(final, tipo_residuo)
  pb <- passos[-1, ]
  efe <- data.frame(indice = c(paste0("ΔR² M", seq_len(nrow(pb)), " (", nomes_bl, ")"), paste0("f² do bloco M", seq_len(nrow(pb))), "R² ajustado do modelo final"),
                    valor = c(pb$delta_R2, pb$f2_bloco, mod$R2_ajustado), ic_inf = NA, ic_sup = NA,
                    principal = c(rep(TRUE, nrow(pb)), rep(TRUE, nrow(pb)), TRUE))
  # coeficientes dos degraus intermediários (só para consulta; a leitura é do modelo final)
  coef_degraus <- lapply(seq_along(blocos)[-length(blocos)], function(b)
    .ee_assoc_coef_lm(modelos[[b + 1]], d, vd, unlist(blocos[seq_len(b)]), prep$info, nivel)$tabela)
  apa <- character()
  for (b in seq_len(nrow(pb))) {
    l <- pb[b, ]
    sig <- l$p < 1 - nivel
    apa <- c(apa, paste0(
      .ee_assoc_tr(paste0("No degrau ", b, ", a inclusão de ", l$entrou, " (", nomes_bl[b], ") "),
                   paste0("At step ", b, ", adding ", l$entrou, " (", nomes_bl[b], ") "), idioma),
      if (sig) .ee_assoc_tr("acrescentou explicação significativa", "significantly increased the explained variance", idioma) else
        .ee_assoc_tr("não acrescentou explicação significativa", "did not significantly increase the explained variance", idioma),
      if (b > 1) .ee_assoc_tr(" acima dos blocos anteriores", " above the previous blocks", idioma) else "", ", ",
      "ΔR² = ", .pr_num(l$delta_R2, 3, FALSE, decimal), .ee_sep(decimal),
      .ee_estat("F", l$F_mudanca, c(l$gl1, l$gl2), decimal), .ee_sep(decimal), .pr_p(l$p, decimal), .ee_sep(decimal),
      "f² = ", .pr_num(l$f2_bloco, 3, TRUE, decimal), "."))
  }
  apa <- c(apa, paste0(.ee_assoc_tr("O modelo final ", "The final model ", idioma),
                       .ee_assoc_tr(paste0("explicou ", .ee_assoc_pct(max(mod$R2_ajustado, 0), 1, decimal), " da variância de ", vd, " ("),
                                    paste0("explained ", .ee_assoc_pct(max(mod$R2_ajustado, 0), 1, decimal), " of the variance in ", vd, " ("), idioma),
                       .ee_assoc_tr("R² ajustado", "adjusted R²", idioma), " = ", .pr_num(mod$R2_ajustado, 3, FALSE, decimal), "), ",
                       .ee_estat("F", mod$F, c(mod$gl1, mod$gl2), decimal), .ee_sep(decimal), .pr_p(mod$p, decimal), "."),
           paste0(c(.ee_assoc_tr("Coeficientes do modelo final. ", "Final-model coefficients. ", idioma), rep("", nrow(cf$tabela) - 2)),
                  .ee_assoc_frase_coef(cf$tabela, mod$gl2, nivel, casas, decimal, idioma, ep_robusto), "."),
           .ee_assoc_causal(idioma))
  decisoes <- c(paste0("Ordem dos blocos: ", paste(paste0(nomes_bl, " = ", vapply(blocos, paste, "", collapse = ", ")), collapse = " → "),
                       ". A ordem é parte da hipótese e precisa ter sido fixada antes de rodar; reordenar depois de ver os resultados invalida os p (livro). Exploração de outras ordens só declarada como exploração."),
                "Mesma amostra em todos os degraus (listwise em todas as variáveis dos blocos). M0 = modelo nulo (só intercepto), numeração do JASP; no texto, Modelo 1 = M1.",
                "F de mudança por anova(modelo anterior, modelo atual), dois de cada vez (= JASP). f² do bloco = ΔR²/(1 − R² do degrau que inclui o bloco), conta que reproduz os valores do livro.",
                "Coeficientes lidos só no modelo final; os dos degraus intermediários ficam em $coeficientes_degraus para consulta, não para relato.",
                .ee_assoc_nota_sem_ic("R² ajustado, ΔR² e f²", "Adjusted R², ΔR² and f²", idioma),
                "Pressupostos checados no modelo final.")
  txt_fat <- .ee_assoc_txt_fatores(prep$info)
  if (length(txt_fat)) decisoes <- c(decisoes, paste0("Preditor(es) categórico(s) em indicadoras k − 1: ", txt_fat, "."))
  .ee_assoc_novo(analise = paste0("Regressão linear hierárquica (", vd, ")"), n = nrow(d), n_excluidos = dd$n_excluidos,
                 tabela = passos, efeito = efe, pressupostos = dg$pressupostos, decisoes = decisoes, avisos = c(prep$avisos, dg$avisos),
                 apa = apa, extra = list(modelos = modelos, coeficientes_degraus = coef_degraus, extra_diagnostico = dg$indices),
                 tabelas_extra = list(`Coeficientes do modelo final` = cf$tabela))
}

# ================================================================ curva ROC
.ee_assoc_faixa_auc <- function(auc, idioma = "pt") {
  if (is.na(auc)) return(NA_character_)
  r <- if (auc < .70) c("abaixo do aceitável", "below acceptable") else if (auc < .80) c("aceitável", "acceptable") else
    if (auc < .90) c("excelente", "excellent") else c("excepcional", "outstanding")
  .ee_assoc_tr(r[1], r[2], idioma)
}

# AUC com IC de DeLong e ponto de corte de Youden (pROC). Direção fixada: com
# maior_indica_evento = TRUE, escores maiores apontam o evento (pROC direction "<").
.ee_assoc_roc <- function(y, score, nivel = .95, maior_indica_evento = TRUE) {
  .pr_exige("pROC")
  r <- pROC::roc(response = y, predictor = score, levels = c(0, 1),
                 direction = if (maior_indica_evento) "<" else ">", quiet = TRUE)
  ci <- as.numeric(suppressWarnings(pROC::ci.auc(r, conf.level = nivel, method = "delong")))
  yd <- pROC::coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
  list(roc = r, auc = ci[2], ic = ci[c(1, 3)], youden = as.data.frame(yd))
}

.ee_assoc_classifica <- function(y, score, corte, maior_indica_evento = TRUE) {
  prev <- if (maior_indica_evento) as.integer(score >= corte) else as.integer(score <= corte)
  vp <- sum(prev == 1 & y == 1); vn <- sum(prev == 0 & y == 0); fp <- sum(prev == 1 & y == 0); fn <- sum(prev == 0 & y == 1)
  data.frame(corte = corte, sensibilidade = vp / (vp + fn), especificidade = vn / (vn + fp),
             youden_J = vp / (vp + fn) + vn / (vn + fp) - 1, vpp = if (vp + fp) vp / (vp + fp) else NA,
             vpn = if (vn + fn) vn / (vn + fn) else NA, VP = vp, FN = fn, VN = vn, FP = fp)
}

# v1.4: prevalencia (proporção de casos na população onde o teste vai ser usado) recalcula
# VPP e VPN por Bayes, VPP = Se·π / [Se·π + (1 − Sp)(1 − π)], porque os da amostra só
# valem para a prevalência da amostra; custo_fn_fp (quanto custa um falso negativo em
# relação a um falso positivo) com a prevalência escolhe o corte que maximiza
# Se + r·Sp, r = (1 − π) / (custo·π) (pROC, best.weights). Razões de verossimilhança
# RV+ = Se / (1 − Sp) e RV− = (1 − Se) / Sp saem para cada corte.
curva_roc <- function(dados, desfecho, preditor, evento = NULL, maior_indica_evento = TRUE, cortes = NULL,
                      prevalencia = NULL, custo_fn_fp = NULL, nivel = .95, casas = 3, decimal = ".", idioma = "pt") {
  if (!is.null(prevalencia) && (prevalencia <= 0 || prevalencia >= 1)) stop("prevalencia deve estar entre 0 e 1")
  if (!is.null(custo_fn_fp) && custo_fn_fp <= 0) stop("custo_fn_fp deve ser positivo")
  dados <- .ee_assoc_vazios(dados, c(desfecho, preditor))
  dd <- .ee_dados(dados, c(desfecho, preditor), numericas = preditor)
  b <- .ee_assoc_binaria(dd$d[[desfecho]], evento, desfecho)
  y <- b$y; s <- dd$d[[preditor]]
  if (sum(y == 1) < 2 || sum(y == 0) < 2) stop("cada categoria do desfecho precisa de pelo menos 2 casos")
  ro <- .ee_assoc_roc(y, s, nivel, maior_indica_evento)
  yd <- ro$youden
  # cortes ±Inf (classificar todos num grupo, J = 0) não são pontos de corte utilizáveis
  yd <- yd[is.finite(yd$threshold), , drop = FALSE]
  if (!nrow(yd)) yd <- data.frame(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_)
  tab <- data.frame(origem = paste0("Youden", if (nrow(yd) > 1) paste0(" (", seq_len(nrow(yd)), ")") else ""),
                    corte = yd$threshold, sensibilidade = yd$sensitivity, especificidade = yd$specificity,
                    youden_J = yd$sensitivity + yd$specificity - 1)
  prev_amostra <- mean(y == 1)
  pi_alvo <- if (is.null(prevalencia)) prev_amostra else prevalencia
  if (!is.null(custo_fn_fp)) {
    wc <- pROC::coords(ro$roc, "best", best.method = "youden", best.weights = c(custo_fn_fp, pi_alvo),
                       ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
    wc <- as.data.frame(wc)[is.finite(as.data.frame(wc)$threshold), , drop = FALSE]
    if (nrow(wc)) tab <- rbind(tab, data.frame(origem = paste0("custo FN/FP = ", custo_fn_fp, ", prevalência ", signif(pi_alvo, 3)),
                                               corte = wc$threshold, sensibilidade = wc$sensitivity, especificidade = wc$specificity,
                                               youden_J = wc$sensitivity + wc$specificity - 1))
  }
  if (length(cortes)) {
    ct <- do.call(rbind, lapply(cortes, function(cc) .ee_assoc_classifica(y, s, cc, maior_indica_evento)))
    tab <- rbind(tab, data.frame(origem = "informado", ct[, c("corte", "sensibilidade", "especificidade", "youden_J")]))
    extra_cortes <- ct
  } else extra_cortes <- NULL
  bayes <- function(se, sp, p) c(vpp = se * p / (se * p + (1 - sp) * (1 - p)), vpn = sp * (1 - p) / (sp * (1 - p) + (1 - se) * p))
  vv <- t(mapply(function(se, sp) bayes(se, sp, prev_amostra), tab$sensibilidade, tab$especificidade))
  tab$vpp_amostra <- vv[, "vpp"]; tab$vpn_amostra <- vv[, "vpn"]
  if (!is.null(prevalencia)) {
    vp <- t(mapply(function(se, sp) bayes(se, sp, prevalencia), tab$sensibilidade, tab$especificidade))
    tab$vpp_populacao <- vp[, "vpp"]; tab$vpn_populacao <- vp[, "vpn"]
  }
  tab$rv_pos <- ifelse(tab$especificidade < 1, tab$sensibilidade / (1 - tab$especificidade), Inf)
  tab$rv_neg <- ifelse(tab$especificidade > 0, (1 - tab$sensibilidade) / tab$especificidade, NA)
  faixa <- .ee_assoc_faixa_auc(ro$auc, idioma)
  avisos <- character()
  invertida <- ro$auc <= .5
  if (invertida)
    avisos <- c(avisos, paste0("AUC = ", round(ro$auc, 3), " (≤ .50): ", preditor, " não discrimina ou discrimina na direção oposta à declarada. ",
                               "Conferir o evento (", b$frase, ") e maior_indica_evento; a frase do ponto de corte não foi gerada."))
  if (nrow(yd) > 1) avisos <- c(avisos, "Mais de um corte empata no índice de Youden: todos listados; a escolha é da pergunta clínica.")
  sinal <- if (maior_indica_evento) "≥" else "≤"
  frase_auc <- paste0(.ee_assoc_tr(paste0("A área sob a curva ROC de ", preditor, " para identificar ", desfecho, " = '", b$evento, "' foi "),
                                   paste0("The area under the ROC curve of ", preditor, " for identifying ", desfecho, " = '", b$evento, "' was "), idioma),
                      .ee_assoc_ef("AUC", ro$auc, ro$ic[1], ro$ic[2], TRUE, casas, decimal, idioma, nivel),
                      .ee_assoc_tr(paste0(" (DeLong), discriminação ", faixa, " pela régua de Hosmer et al. (2013)."),
                                   paste0(" (DeLong), ", faixa, " discrimination by Hosmer et al.'s (2013) benchmarks."), idioma))
  frase_corte <- if (!invertida && is.finite(yd$threshold[1]))
    paste0(.ee_assoc_tr("O corte que maximiza o índice de Youden (positivo se ", "The cutoff maximizing Youden's index (positive if ", idioma),
           preditor, " ", sinal, " ", .ee_assoc_n(yd$threshold[1], 2, TRUE, decimal), ") ",
           .ee_assoc_tr("teve sensibilidade = ", "had sensitivity = ", idioma), .pr_num(yd$sensitivity[1], casas, FALSE, decimal),
           .ee_assoc_tr(" e especificidade = ", " and specificity = ", idioma), .pr_num(yd$specificity[1], casas, FALSE, decimal), ".") else NULL
  apa <- paste(c(frase_auc, frase_corte), collapse = " ")
  .ee_assoc_novo(analise = paste0("Curva ROC (", preditor, " → ", desfecho, ")"), n = length(y), n_excluidos = dd$n_excluidos,
                 tabela = tab,
                 efeito = data.frame(indice = "AUC", valor = ro$auc, ic_inf = ro$ic[1], ic_sup = ro$ic[2], principal = TRUE),
                 decisoes = c(b$frase,
                              "Regra do evento (a mesma da regressao_logistica): desfecho 0/1 = 1; lógico = TRUE; fator = segundo nível; texto exige evento = '...'.",
                              paste0("Direção fixada: ", if (maior_indica_evento) "valores maiores" else "valores menores", " de ", preditor,
                                     " indicam o evento (o pROC escolheria sozinho com direction = \"auto\" e poderia inverter a AUC)."),
                              paste0("IC ", round(nivel * 100), "% da AUC pelo método de DeLong (padrão do pROC). AUC = probabilidade de um caso com evento sorteado ter escore mais extremo que um sem evento; .50 = acaso."),
                              "Régua de Hosmer et al. (2013): .70 aceitável, .80 excelente, .90 excepcional (o curso usa rótulos diferentes para as mesmas faixas). Réguas de partida, não leis.",
                              "Corte de Youden = maximiza sensibilidade + especificidade − 1 (o pROC usa pontos médios entre valores observados; cortes ±Inf são descartados). Qual erro custa mais e onde pôr o corte é decisão da pergunta, não da estatística (livro).",
                              paste0("VPP e VPN pela prevalência da amostra (", .pr_num(100 * prev_amostra, 1, TRUE, decimal), "%)",
                                     if (!is.null(prevalencia)) paste0(" e pela prevalência informada (", .pr_num(100 * prevalencia, 1, TRUE, decimal),
                                                                      "%), por Bayes; os da amostra só valem onde a prevalência for a mesma") else
                                       "; para a população de uso, informar prevalencia: com prevalência baixa o VPP cai mesmo com sensibilidade e especificidade altas",
                                     ". RV+ e RV− não dependem da prevalência."),
                              if (!is.null(custo_fn_fp)) paste0("Corte ponderado por custo: maximiza Se + r·Sp com r = (1 − π)/(custo·π) = ",
                                                                 .pr_num((1 - pi_alvo) / (custo_fn_fp * pi_alvo), 3, TRUE, decimal), " (pROC, best.weights).")),
                 avisos = avisos, apa = apa,
                 extra = list(roc = ro$roc, classificacao_cortes = extra_cortes),
                 tabelas_extra = list(`Casos por categoria` = data.frame(categoria = c(b$outro, b$evento), codigo = 0:1, n = c(sum(y == 0), sum(y == 1)))))
}

# ================================================================ regressão logística binária
.ee_assoc_mudanca_chance <- function(or, termo_ref, idioma = "pt", decimal = ".", casas = 3) {
  if (is.na(or)) return("")
  if (or >= 2) return(.ee_assoc_tr(paste0("chance ", .ee_assoc_n(or, casas, TRUE, decimal), " vezes a ", termo_ref),
                                   paste0("odds ", .ee_assoc_n(or, casas, TRUE, decimal), " times those ", termo_ref), idioma))
  pct <- .pr_num(abs(or - 1) * 100, 1, TRUE, decimal)
  if (or >= 1) .ee_assoc_tr(paste0("chance ", pct, "% maior ", termo_ref), paste0("odds ", pct, "% higher ", termo_ref), idioma) else
    .ee_assoc_tr(paste0("chance ", pct, "% menor ", termo_ref), paste0("odds ", pct, "% lower ", termo_ref), idioma)
}

# validar = TRUE (v1.6): validação interna por bootstrap com correção de otimismo (Harrell,
# Lee & Mark, 1996; Steyerberg et al., 2001): reajusta o modelo em cada reamostra, mede
# AUC, Brier e a inclinação e o intercepto de calibração na reamostra e no banco original,
# e desconta a média da diferença. No próprio banco a inclinação de calibração vale 1 e o
# intercepto 0 por construção (máxima verossimilhança); o valor corrigido mostra o
# encolhimento esperado em dados novos (Petersen, 2024, cap. 9, sobreajuste).
.ee_assoc_desempenho_lp <- function(y, lp) {
  p <- stats::plogis(lp)
  auc <- { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
  cal <- tryCatch(stats::coef(suppressWarnings(stats::glm(y ~ lp, family = stats::binomial()))), error = function(e) c(NA, NA))
  citl <- tryCatch(stats::coef(suppressWarnings(stats::glm(y ~ 1 + offset(lp), family = stats::binomial())))[[1]], error = function(e) NA_real_)
  c(AUC = auc, Brier = mean((p - y)^2), inclinacao = unname(cal[2]), intercepto = citl)
}
.ee_assoc_validar <- function(f, d, reamostras, semente) {
  ap <- { m <- suppressWarnings(stats::glm(f, family = stats::binomial(), data = d)); .ee_assoc_desempenho_lp(d$.y_evento, stats::predict(m)) }
  oti <- .ee_assoc_com_semente(semente, vapply(seq_len(reamostras), function(b) {
    db <- d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]
    mb <- tryCatch(suppressWarnings(stats::glm(f, family = stats::binomial(), data = db)), error = function(e) NULL)
    if (is.null(mb) || any(is.na(stats::coef(mb)))) return(rep(NA_real_, 4))
    .ee_assoc_desempenho_lp(db$.y_evento, stats::predict(mb)) -
      .ee_assoc_desempenho_lp(d$.y_evento, stats::predict(mb, newdata = d))
  }, numeric(4)))
  o <- rowMeans(oti, na.rm = TRUE)
  data.frame(indice = names(ap), aparente = unname(ap), otimismo = unname(o), corrigido = unname(ap - o),
             reamostras_validas = sum(stats::complete.cases(t(oti))))
}

regressao_logistica <- function(dados, vd, preditores, evento = NULL, referencias = NULL, fatores = NULL, corte = .5,
                                ic = c("wald", "perfil"), box_tidwell = TRUE, validar = FALSE, reamostras = 500, semente = 2026,
                                nivel = .95, casas = 3, decimal = ".", idioma = "pt") {
  ic <- match.arg(ic)
  vars <- c(vd, preditores)
  dados <- .ee_assoc_vazios(dados, vars)
  n_total <- nrow(dados)
  dd <- .ee_dados(dados, vars)
  b <- .ee_assoc_binaria(dd$d[[vd]], evento, vd)
  prep <- .ee_assoc_preparar(dd$d, preditores, referencias, fatores)
  d <- prep$d; d$.y_evento <- b$y
  f <- stats::as.formula(paste(".y_evento ~", paste(.ee_assoc_bq(preditores), collapse = " + ")))
  avisos_glm <- character()
  ajusta <- function(fm) withCallingHandlers(stats::glm(fm, family = stats::binomial(), data = d),
                                             warning = function(w) { avisos_glm <<- c(avisos_glm, conditionMessage(w)); invokeRestart("muffleWarning") })
  fit <- ajusta(f)
  if (any(is.na(stats::coef(fit)))) stop("coeficiente não estimável: preditores perfeitamente colineares ou categoria sem casos")
  n <- nrow(d); y <- b$y; p_hat <- stats::fitted(fit)
  X <- stats::model.matrix(fit)
  sm <- summary(fit)$coefficients
  ci_w <- stats::confint.default(fit, level = nivel)
  tab <- data.frame(termo = .ee_assoc_nomes_coef(X, preditores, prep$info), B = sm[, 1], EP = sm[, 2], z = sm[, 3],
                    wald = sm[, 3]^2, gl = 1, p = sm[, 4], OR = exp(sm[, 1]), OR_ic_inf = exp(ci_w[, 1]),
                    OR_ic_sup = exp(ci_w[, 2]), mudanca_chance_pct = 100 * (exp(sm[, 1]) - 1), row.names = NULL)
  tab$mudanca_chance_pct[1] <- NA
  if (ic == "perfil") {
    ci_p <- suppressMessages(stats::confint(fit, level = nivel))
    tab$OR_ic_inf_perfil <- exp(ci_p[, 1]); tab$OR_ic_sup_perfil <- exp(ci_p[, 2])
  }
  vf <- .ee_assoc_vif(fit); asg <- attr(X, "assign")
  if (length(preditores) >= 2) { tab$vif <- c(NA, vf$vif[asg[asg > 0]]) }
  # modelo vs nulo
  dev0 <- fit$null.deviance; dev <- fit$deviance; gl <- fit$df.null - fit$df.residual
  qui <- dev0 - dev; pq <- stats::pchisq(qui, gl, lower.tail = FALSE)
  l0 <- -dev0 / 2; lm_ <- as.numeric(stats::logLik(fit))
  mcf <- 1 - dev / dev0; cs <- 1 - exp(2 * (l0 - lm_) / n); nk <- cs / (1 - exp(2 * l0 / n))
  tjur <- abs(mean(p_hat[y == 1]) - mean(p_hat[y == 0]))
  mod <- data.frame(n = n, deviance_nulo = dev0, deviance_modelo = dev, qui2 = qui, gl = gl, p = pq,
                    AIC = stats::AIC(fit), BIC = stats::BIC(fit))
  # classificação no corte (JASP: classe 1 se probabilidade > corte)
  prev <- as.integer(p_hat > corte)
  vp <- sum(prev == 1 & y == 1); vn <- sum(prev == 0 & y == 0); fp <- sum(prev == 1 & y == 0); fn <- sum(prev == 0 & y == 1)
  conf <- data.frame(observado = c(b$outro, b$evento), previsto_0 = c(vn, fn), previsto_1 = c(fp, vp),
                     pct_correto = 100 * c(vn / (vn + fp), vp / (vp + fn)))
  names(conf)[2:3] <- paste0("previsto ", c(b$outro, b$evento))
  ro <- .ee_assoc_roc(y, p_hat, nivel, TRUE)
  base <- max(mean(y), 1 - mean(y)); brier <- mean((p_hat - y)^2); brier0 <- mean(y) * (1 - mean(y))
  prec <- if (vp + fp) vp / (vp + fp) else NA; sens <- vp / (vp + fn); esp <- vn / (vn + fp)
  desemp <- data.frame(corte = corte, acuracia = (vp + vn) / n, taxa_base = base, sensibilidade = sens, especificidade = esp,
                       precisao = prec, vpn = if (vn + fn) vn / (vn + fn) else NA, F1 = 2 * prec * sens / (prec + sens),
                       Brier = brier, Brier_escalado = 1 - brier / brier0, AUC = ro$auc, AUC_ic_inf = ro$ic[1], AUC_ic_sup = ro$ic[2])
  # pressupostos
  linhas <- list(); avisos <- character()
  add <- function(checagem, estatistica = NA, gl = NA, p = NA, leitura) linhas[[length(linhas) + 1]] <<-
    data.frame(checagem = checagem, estatistica = estatistica, gl = gl, p = p, leitura = leitura)
  if (length(preditores) >= 2) {
    i <- which.max(vf$vif); v <- vf$vif[i]
    add("Multicolinearidade (maior VIF)", v, NA, NA,
        paste0("VIF = ", sprintf("%.2f", v), " em ", preditores[i], if (v > 10) ": acima de 10, colinearidade forte" else if (v > 5) ": entre 5 e 10, alerta" else ": abaixo das réguas de alerta (5 e 10)",
               "; mesma checagem da linear (livro)."))
  }
  ev <- min(sum(y == 1), sum(y == 0)); npar <- ncol(X) - 1; epv <- ev / npar
  add("Eventos por parâmetro (EPV)", epv, NA, NA,
      paste0(ev, " casos na categoria menos frequente para ", npar, " parâmetro(s)", if (epv < 10) ": abaixo de 10, estimativas instáveis e OR inflados; reduzir preditores ou ampliar a amostra" else ": acima da regra usual de 10", "."))
  if (epv < 10) avisos <- c(avisos, paste0("EPV = ", round(epv, 1), " (< 10): risco de coeficientes enviesados."))
  sep <- any(grepl("fitted probabilities numerically 0 or 1|did not converge", avisos_glm)) || any(abs(sm[-1, 1]) > 10)
  add("Separação (quase) completa", NA, NA, NA,
      if (sep) "Há sinais de separação (probabilidades ajustadas 0 ou 1, não convergência ou |B| > 10): OR e IC não são confiáveis; rever categorias raras ou usar regressão penalizada (Firth)." else
        "Sem sinais de separação (o ajuste convergiu sem probabilidades 0 ou 1).")
  if (sep) avisos <- c(avisos, "Possível separação (quase) completa: OR e IC não confiáveis; as frases APA trazem \"OR não estimável\" no lugar dos números (alternativas: Firth, métodos exatos).")
  if (box_tidwell) {
    cont <- preditores[vapply(preditores, function(v) prep$info[[v]]$tipo == "numerica" && length(unique(d[[v]])) > 2, logical(1))]
    if (length(cont)) {
      d_bt <- d; termos_bt <- character(); desloc <- character()
      for (v in cont) {
        x <- d[[v]]
        if (min(x) <= 0) { x <- x - min(x) + 1; desloc <- c(desloc, v) }
        nm <- paste0(".bt_", make.names(v)); d_bt[[nm]] <- x * log(x); termos_bt <- c(termos_bt, nm)
      }
      fbt <- stats::as.formula(paste(".y_evento ~", paste(c(.ee_assoc_bq(preditores), termos_bt), collapse = " + ")))
      fit_bt <- suppressWarnings(stats::glm(fbt, family = stats::binomial(), data = d_bt))
      sbt <- summary(fit_bt)$coefficients
      for (j in seq_along(cont)) {
        pj <- if (termos_bt[j] %in% rownames(sbt)) sbt[termos_bt[j], 4] else NA
        add(paste0("Linearidade do logito: ", cont[j], " (Box-Tidwell)"), if (termos_bt[j] %in% rownames(sbt)) sbt[termos_bt[j], 3] else NA, NA, pj,
            paste0("Termo x·ln(x) ", if (!is.na(pj) && pj < .05) "significativo: a relação com o logito pode não ser linear; considerar transformação ou termo polinomial" else "não significativo: sem indício contra a linearidade do logito",
                   if (cont[j] %in% desloc) " (x deslocado para valores positivos: x − mín + 1)" else "",
                   ". Com n grande o teste acusa desvios pequenos; ler junto com o gráfico."))
      }
    }
  }
  rp <- stats::rstandard(fit, type = "pearson"); ck <- stats::cooks.distance(fit)
  add("Casos atípicos e influentes", sum(abs(rp) > 3), NA, NA,
      paste0(sum(abs(rp) > 3), " caso(s) com resíduo de Pearson padronizado |z| > 3; ", sum(ck > 1), " com distância de Cook > 1 (critérios do JASP)."))
  # efeitos
  tt <- tab[-1, ]
  efe <- rbind(data.frame(indice = paste0("OR ", tt$termo), valor = tt$OR, ic_inf = tt$OR_ic_inf, ic_sup = tt$OR_ic_sup, principal = TRUE),
               data.frame(indice = c("R² de Nagelkerke", "R² de Tjur", "R² de McFadden", "R² de Cox e Snell", "AUC"),
                          valor = c(nk, tjur, mcf, cs, ro$auc), ic_inf = c(NA, NA, NA, NA, ro$ic[1]), ic_sup = c(NA, NA, NA, NA, ro$ic[2]),
                          principal = c(TRUE, TRUE, FALSE, FALSE, FALSE)))
  # APA
  ref_txt <- function(termo) {
    v <- preditores[vapply(preditores, function(p) startsWith(termo, paste0(p, ": ")) || termo == p, logical(1))][1]
    if (prep$info[[v]]$tipo == "numerica") .ee_assoc_tr(paste0("por unidade de ", v), paste0("per unit of ", v), idioma) else
      .ee_assoc_tr(paste0("que em ", prep$info[[v]]$ref), paste0("than in ", prep$info[[v]]$ref), idioma)
  }
  mud <- function(l) {
    v <- preditores[vapply(preditores, function(p) startsWith(l$termo, paste0(p, ": ")) || l$termo == p, logical(1))][1]
    if (l$OR < 2) return(.ee_assoc_mudanca_chance(l$OR, ref_txt(l$termo), idioma, decimal, casas))
    o <- .ee_assoc_n(l$OR, casas, TRUE, decimal)
    if (prep$info[[v]]$tipo == "numerica") .ee_assoc_tr(paste0("chance multiplicada por ", o, " a cada unidade de ", v), paste0("odds multiplied by ", o, " per unit of ", v), idioma) else
      .ee_assoc_tr(paste0("chance ", o, " vezes a de ", prep$info[[v]]$ref), paste0("odds ", o, " times those of ", prep$info[[v]]$ref), idioma)
  }
  # separação: o OR dos termos afetados não é estimável por máxima verossimilhança
  sem_or <- if (sep) { f <- abs(tt$B) > 10 | tt$EP > 10 | !is.finite(tt$OR_ic_sup) | tt$OR_ic_inf == 0; if (any(f)) f else rep(TRUE, nrow(tt)) } else rep(FALSE, nrow(tt))
  frases_or <- vapply(seq_len(nrow(tt)), function(i) {
    l <- tt[i, ]
    if (sem_or[i]) return(paste0(l$termo, .ee_assoc_tr(": OR não estimável por máxima verossimilhança (separação (quase) completa); usar regressão logística penalizada de Firth ou métodos exatos.",
                                                       ": OR not estimable by maximum likelihood ((quasi-)complete separation); use Firth's penalized logistic regression or exact methods.", idioma)))
    paste0(l$termo, ": ", .ee_assoc_ef("OR", l$OR, l$OR_ic_inf, l$OR_ic_sup, FALSE, casas, decimal, idioma, nivel),
           .ee_sep(decimal), .pr_p(l$p, decimal), " (",
           mud(l), ").")
  }, character(1))
  apa <- c(
    if (n_total > n) .ee_assoc_tr(paste0("Dos ", n_total, " casos do banco, ", n_total - n, " foram excluídos por ausência de dados nas variáveis do modelo, resultando em N = ", n, "."),
                                  paste0("Of the ", n_total, " cases, ", n_total - n, " were excluded for missing data on model variables, leaving N = ", n, "."), idioma),
    paste0(if (pq < 1 - nivel) .ee_assoc_tr("O modelo foi estatisticamente significativo, ", "The model was statistically significant, ", idioma) else
      .ee_assoc_tr("O modelo não foi estatisticamente significativo, ", "The model was not statistically significant, ", idioma),
      "χ²(", gl, ") = ", .pr_num(qui, 2, TRUE, decimal), .ee_sep(decimal), .pr_p(pq, decimal),
      if (sep) .ee_assoc_tr(" (pseudo-R² não interpretáveis sob separação).", " (pseudo-R² not interpretable under separation).", idioma) else
        paste0(.ee_assoc_tr(", com R² de Nagelkerke = ", ", with Nagelkerke R² = ", idioma), .pr_num(nk, 3, FALSE, decimal),
               .ee_assoc_tr(" e R² de Tjur = ", " and Tjur R² = ", idioma), .pr_num(tjur, 3, FALSE, decimal), ".")),
    frases_or,
    paste0(.ee_assoc_tr("O modelo classificou corretamente ", "The model correctly classified ", idioma), .ee_assoc_pct(desemp$acuracia, 1, decimal),
           .ee_assoc_tr(" dos casos (corte = ", " of cases (cutoff = ", idioma), .pr_num(corte, 2, FALSE, decimal),
           .ee_assoc_tr("), ante ", "), versus ", idioma), .ee_assoc_pct(base, 1, decimal),
           .ee_assoc_tr(" da classificação pela categoria mais frequente; sensibilidade = ", " by the most frequent category; sensitivity = ", idioma),
           .pr_num(sens, casas, FALSE, decimal), .ee_assoc_tr(", especificidade = ", ", specificity = ", idioma), .pr_num(esp, casas, FALSE, decimal),
           .ee_sep(decimal), .ee_assoc_ef("AUC", ro$auc, ro$ic[1], ro$ic[2], TRUE, casas, decimal, idioma, nivel), "."))
  txt_fat <- .ee_assoc_txt_fatores(prep$info)
  decisoes <- c(b$frase,
                if (length(txt_fat)) paste0("Fatores em indicadoras k − 1: ", txt_fat,
                                            ". Antes de ler qualquer OR: evento é o quê, e comparado com quem? (livro)."),
                "Teste do modelo: diferença de deviance (−2LL) entre nulo e modelo, χ² com gl = parâmetros acrescentados (dummies contam um a um).",
                paste0("Coeficientes: z de Wald do glm e Wald = z² (gl = 1), como o JASP. IC ", round(nivel * 100), "% do OR por Wald (confint.default, = JASP)",
                       if (ic == "perfil") "; IC por perfil de verossimilhança ao lado (confint)." else "; ic = \"perfil\" acrescenta o IC por perfil de verossimilhança."),
                "OR < 1 lê-se como \"a chance se multiplica por OR\" ou \"(1 − OR) × 100% menor\"; nunca \"OR vezes menos chance\" (livro). OR é razão de chances, não de probabilidades.",
                "Pseudo-R² em relação ao modelo nulo, fórmulas do JASP: McFadden, Cox e Snell, Nagelkerke e Tjur. Nenhum é proporção de variância explicada; o livro e o curso destacam Nagelkerke e Tjur.",
                paste0("Classificação com corte ", corte, " na probabilidade prevista; acurácia comparada à taxa base (categoria mais frequente = ", .ee_assoc_pct(base, 1), "). AUC com IC de DeLong (pROC)."),
                if (box_tidwell) "Linearidade do logito por Box-Tidwell (termos x·ln(x) no modelo); o JASP não oferece esse teste (livro)." else NULL,
                "Hosmer-Lemeshow não é calculado (não existe no JASP 0.95.1 nem no livro).",
                "Na discussão, a ressalva de desenho observacional: OR ajustado é associação condicional, não mecanismo (livro).")
  val <- NULL
  if (isTRUE(validar)) {
    val <- .ee_assoc_validar(f, d, reamostras, semente)
    decisoes <- c(decisoes, paste0("Validação interna por bootstrap (", reamostras, " reamostras, semente ", semente, "): desempenho corrigido pelo otimismo (Harrell et al., 1996). ",
                                   "Inclinação de calibração corrigida < 1 indica sobreajuste; o fator de encolhimento dos coeficientes é aproximadamente essa inclinação."))
    sl <- val$corrigido[val$indice == "inclinacao"]
    if (is.finite(sl) && sl < .90) avisos <- c(avisos, paste0("Inclinação de calibração corrigida = ", .pr_num(sl, 2, TRUE, decimal),
      ": o modelo promete mais do que vai entregar em dados novos (eventos por parâmetro: ", round(min(sum(y), n - sum(y)) / (ncol(X) - 1), 1),
      "). Reduzir preditores, encolher coeficientes (regressão penalizada) ou validar externamente."))
    apa <- c(apa, paste0(.ee_assoc_tr("Na validação interna por bootstrap (", "In bootstrap internal validation (", idioma), reamostras,
                         .ee_assoc_tr(" reamostras), a AUC corrigida pelo otimismo foi ", " resamples), the optimism-corrected AUC was ", idioma),
                         .pr_num(val$corrigido[1], casas, FALSE, decimal), .ee_assoc_tr(" e a inclinação de calibração, ", " and the calibration slope ", idioma),
                         .pr_num(sl, 2, TRUE, decimal), "."))
  }
  .ee_assoc_novo(analise = paste0("Regressão logística binária (", vd, " = '", b$evento, "')"), n = n, n_excluidos = dd$n_excluidos,
                 tabela = tab, efeito = efe, pressupostos = do.call(rbind, linhas), decisoes = decisoes,
                 avisos = c(prep$avisos, avisos, if (length(avisos_glm)) paste0("Avisos do glm: ", paste(unique(avisos_glm), collapse = "; "))),
                 apa = apa, extra = list(modelo = fit, roc = ro$roc, validacao = val),
                 tabelas_extra = c(list(`Modelo vs nulo` = mod, `Matriz de confusão` = conf, `Desempenho` = desemp),
                                   if (!is.null(val)) list(`Validação interna (bootstrap)` = val)))
}

# ================================================================ qui-quadrado
# Semente local: não altera o gerador aleatório da sessão do usuário.
.ee_assoc_com_semente <- function(semente, expr) {
  velho <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(velho)) { if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()) } else
    assign(".Random.seed", velho, envir = globalenv()))
  set.seed(semente)
  expr
}

.ee_assoc_faixa_w <- function(w, idioma = "pt") {
  r <- if (w < .10) c("muito pequeno (abaixo de .10)", "very small (below .10)") else if (w < .30) c("pequeno", "small") else
    if (w < .50) c("médio", "medium") else c("grande", "large")
  .ee_assoc_tr(r[1], r[2], idioma)
}

# Regra das frequências esperadas (livro; Cochran, 1954): no máximo 20% das
# células com esperado < 5 e nenhuma < 1.
.ee_assoc_esperados <- function(E) {
  list(prop_menor5 = mean(E < 5), minimo = min(E), algum_menor1 = any(E < 1),
       violada = mean(E < 5) > .20 || any(E < 1))
}

.ee_assoc_tabela_cont <- function(dados, linhas, colunas) {
  if (is.table(dados) || is.matrix(dados)) {
    mm <- as.matrix(dados)
    if (is.null(rownames(mm))) rownames(mm) <- paste0("L", seq_len(nrow(mm)))
    if (is.null(colnames(mm))) colnames(mm) <- paste0("C", seq_len(ncol(mm)))
    tab <- as.table(mm)
    nm <- names(dimnames(tab))
    if (is.null(nm) || any(!nzchar(nm))) names(dimnames(tab)) <- c(if (is.null(linhas)) "linhas" else linhas, if (is.null(colunas)) "colunas" else colunas)
    return(list(tab = tab, n_excluidos = 0))
  }
  dados <- .ee_assoc_vazios(dados, c(linhas, colunas))
  dd <- .ee_dados(dados, c(linhas, colunas))
  fl <- droplevels(as.factor(dd$d[[linhas]])); fc <- droplevels(as.factor(dd$d[[colunas]]))
  tab <- table(fl, fc, dnn = c(linhas, colunas))
  list(tab = tab, n_excluidos = dd$n_excluidos)
}

qui_quadrado <- function(dados, linhas = NULL, colunas = NULL, evento = NULL, referencia = NULL,
                         exato = c("auto", "sim", "nao"), correcao_continuidade = FALSE, alfa_residuos = .05,
                         semente = 2026, simulacoes = 10000, nivel = .95, casas = 3, decimal = ".", idioma = "pt") {
  exato <- match.arg(exato)
  tc <- .ee_assoc_tabela_cont(dados, linhas, colunas)
  tab <- tc$tab; nl <- nrow(tab); nc <- ncol(tab); N <- sum(tab)
  vl <- names(dimnames(tab))[1]; vc <- names(dimnames(tab))[2]
  if (nl < 2 || nc < 2) stop("a tabela precisa de pelo menos 2 linhas e 2 colunas com dados")
  ct <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  E <- ct$expected; esp <- .ee_assoc_esperados(E)
  dois <- nl == 2 && nc == 2
  testes <- data.frame(teste = "χ² de Pearson", estatistica = unname(ct$statistic), gl = unname(ct$parameter), p = ct$p.value)
  if (dois) {
    cy <- suppressWarnings(stats::chisq.test(tab, correct = TRUE))
    testes <- rbind(testes, data.frame(teste = "χ² com correção de continuidade (Yates)", estatistica = unname(cy$statistic), gl = 1, p = cy$p.value))
  }
  usa_exato <- exato == "sim" || (exato == "auto" && esp$violada)
  fisher <- NULL
  if (dois || usa_exato) {
    fisher <- tryCatch(stats::fisher.test(tab, conf.level = nivel, workspace = 2e7), error = function(err) NULL)
    if (!is.null(fisher)) {
      testes <- rbind(testes, data.frame(teste = "Fisher exato", estatistica = NA, gl = NA, p = fisher$p.value))
    } else {
      fisher <- .ee_assoc_com_semente(semente, stats::fisher.test(tab, simulate.p.value = TRUE, B = simulacoes))
      testes <- rbind(testes, data.frame(teste = paste0("Fisher por Monte Carlo (B = ", simulacoes, ", semente ", semente, ")"),
                                         estatistica = NA, gl = NA, p = fisher$p.value))
    }
  }
  lido <- if (usa_exato && !is.null(fisher)) nrow(testes) else if (dois && correcao_continuidade) 2 else 1
  testes$lido <- seq_len(nrow(testes)) == lido
  # efeitos
  k <- min(nl, nc)
  V <- sqrt(unname(ct$statistic) / (N * (k - 1))); w <- sqrt(unname(ct$statistic) / N)
  # IC bilateral do V (e do φ no 2×2) pelo qui-quadrado não central (effectsize, sem correção de viés)
  .pr_exige("effectsize")
  v_ic <- tryCatch({ e <- effectsize::cramers_v(tab, adjust = FALSE, ci = nivel, alternative = "two.sided"); c(e$CI_low, e$CI_high) },
                   error = function(err) c(NA_real_, NA_real_))
  efe <- data.frame(indice = "V de Cramér", valor = V, ic_inf = v_ic[1], ic_sup = v_ic[2], principal = !dois)
  extra_tab <- list()
  cel <- as.data.frame(tab, stringsAsFactors = FALSE); names(cel) <- c(vl, vc, "observado")
  cel$esperado <- as.vector(E)
  cel$pct_na_linha <- as.vector(100 * prop.table(tab, 1))
  cel$residuo_ajustado <- as.vector(ct$stdres)
  m_cel <- nl * nc
  z1 <- stats::qnorm(1 - alfa_residuos / 2); zb <- stats::qnorm(1 - alfa_residuos / (2 * m_cel))
  cel$acima_196 <- abs(cel$residuo_ajustado) > z1
  cel$acima_bonferroni <- abs(cel$residuo_ajustado) > zb
  names(cel)[names(cel) == "acima_196"] <- paste0("acima_", formatC(z1, format = "f", digits = 2))
  extra_tab[["Células (observado, esperado, % na linha, resíduo ajustado)"]] <- cel
  apa <- character(); decisoes <- character(); avisos <- character()
  sig <- testes$p[lido] < 1 - nivel
  est_txt <- if (grepl("Fisher", testes$teste[lido])) paste0(.ee_assoc_tr("teste exato de Fisher, ", "Fisher's exact test, ", idioma), .pr_p(testes$p[lido], decimal)) else
    paste0("χ²(", testes$gl[lido], ", N = ", N, ") = ", .pr_num(testes$estatistica[lido], 2, TRUE, decimal), .ee_sep(decimal), .pr_p(testes$p[lido], decimal))
  if (dois) {
    a <- tab[1, 1]; bb <- tab[1, 2]; cc <- tab[2, 1]; dd_ <- tab[2, 2]
    phi <- sign(a * dd_ - bb * cc) * w
    phi_ic <- tryCatch({ e <- effectsize::phi(tab, adjust = FALSE, ci = nivel, alternative = "two.sided"); c(e$CI_low, e$CI_high) },
                       error = function(err) c(NA_real_, NA_real_))
    phi_ic_s <- if (phi < 0) -rev(phi_ic) else phi_ic
    efe <- rbind(data.frame(indice = "φ", valor = phi, ic_inf = phi_ic_s[1], ic_sup = phi_ic_s[2], principal = TRUE), efe)
    cols <- colnames(tab); rows <- rownames(tab)
    ev <- if (!is.null(evento)) as.character(evento) else if (all(cols %in% c("0", "1"))) "1" else cols[2]
    regra_ev <- if (!is.null(evento)) "informado" else if (all(cols %in% c("0", "1"))) "código 1" else "segunda coluna da tabela: 2º nível do fator ou, em texto, ordem alfabética"
    if (!ev %in% cols) stop("evento '", ev, "' não é categoria de ", vc)
    rf <- if (!is.null(referencia)) as.character(referencia) else rows[1]
    if (!rf %in% rows) stop("referência '", rf, "' não é categoria de ", vl)
    cp <- setdiff(rows, rf); nev <- setdiff(cols, ev)
    x1 <- tab[cp, ev]; y1 <- tab[cp, nev]; x0 <- tab[rf, ev]; y0 <- tab[rf, nev]
    corr <- if (any(c(x1, y1, x0, y0) == 0)) .5 else 0
    or <- ((x1 + corr) / (y1 + corr)) / ((x0 + corr) / (y0 + corr))
    # Com célula zero o OR observado é infinito, zero ou indeterminado; o que sai é o
    # corrigido (Haldane, +0,5), e a frase e $decisoes dizem isso.
    or_bruto <- (x1 * y0) / (y1 * x0)
    bruto_txt <- if (is.nan(or_bruto)) c("indeterminado (0/0)", "undefined (0/0)") else if (is.infinite(or_bruto)) c("infinito", "infinite") else c("zero", "zero")
    meio <- .pr_num(.5, 1, TRUE, decimal)
    haldane_frase <- if (corr > 0) .ee_assoc_tr(paste0(" (OR com correção de Haldane, +", meio, " em todas as células, por haver célula zero; sem a correção o OR seria ", bruto_txt[1], ")"),
                                               paste0(" (Haldane-corrected OR, +", meio, " added to every cell because of a zero cell; uncorrected, the OR would be ", bruto_txt[2], ")"), idioma) else ""
    se <- sqrt(1 / (x1 + corr) + 1 / (y1 + corr) + 1 / (x0 + corr) + 1 / (y0 + corr))
    zc <- stats::qnorm(1 - (1 - nivel) / 2)
    or_ic <- exp(log(or) + c(-1, 1) * zc * se)
    efe <- rbind(efe, data.frame(indice = paste0("OR (", ev, ": ", cp, " vs ", rf, ")"), valor = or, ic_inf = or_ic[1], ic_sup = or_ic[2], principal = TRUE),
                 data.frame(indice = paste0("OR inverso (", ev, ": ", rf, " vs ", cp, ")"), valor = 1 / or, ic_inf = 1 / or_ic[2], ic_sup = 1 / or_ic[1], principal = FALSE))
    chances <- data.frame(grupo = c(rf, cp), n = c(x0 + y0, x1 + y1), com_evento = c(x0, x1),
                          probabilidade = c(x0 / (x0 + y0), x1 / (x1 + y1)), chance = c(x0 / y0, x1 / y1))
    names(chances)[1] <- vl
    extra_tab[[paste0("Chance e probabilidade de ", vc, " = ", ev)]] <- chances
    # v1.6: razão de riscos, diferença de riscos (IC de Newcombe, 1998, híbrido de Wilson) e
    # NNT = 1/|DR| com o IC invertido como em Altman (1998): se o IC da DR cruza o zero, o do
    # NNT passa pelo infinito (NNT de benefício de um lado, de dano do outro).
    n1 <- x1 + y1; n0 <- x0 + y0; p1 <- x1 / n1; p0 <- x0 / n0
    wil <- function(x, n) { z <- zc; ctr <- (x + z^2 / 2) / (n + z^2); h <- z * sqrt(x * (n - x) / n + z^2 / 4) / (n + z^2); c(ctr - h, ctr + h) }
    w1 <- wil(x1, n1); w0 <- wil(x0, n0); rd <- p1 - p0
    rd_ic <- c(rd - sqrt((p1 - w1[1])^2 + (w0[2] - p0)^2), rd + sqrt((w1[2] - p1)^2 + (p0 - w0[1])^2))
    rr <- if (x0 > 0 && x1 > 0) p1 / p0 else NA_real_
    rr_ic <- if (is.finite(rr)) exp(log(rr) + c(-1, 1) * zc * sqrt(1 / x1 - 1 / n1 + 1 / x0 - 1 / n0)) else c(NA_real_, NA_real_)
    nnt <- if (rd != 0) 1 / abs(rd) else Inf
    nnt_txt <- if (rd_ic[1] > 0 || rd_ic[2] < 0) paste0(.pr_num(min(1 / abs(rd_ic)), 1, TRUE, decimal), " a ", .pr_num(max(1 / abs(rd_ic)), 1, TRUE, decimal)) else
      paste0(.pr_num(1 / abs(rd_ic[2]), 1, TRUE, decimal), " (", if (rd_ic[2] > 0) "aumento" else "redução", ") a ∞ a ", .pr_num(1 / abs(rd_ic[1]), 1, TRUE, decimal), " (", if (rd_ic[1] < 0) "redução" else "aumento", ")")
    risco <- data.frame(medida = c("razão de riscos (RR)", "diferença de riscos (DR)", "NNT = 1/|DR|"), valor = c(rr, rd, nnt),
                        ic_inf = c(rr_ic[1], rd_ic[1], NA), ic_sup = c(rr_ic[2], rd_ic[2], NA), ic_texto = c(NA, NA, nnt_txt))
    extra_tab[[paste0("Risco de ", vc, " = ", ev, ": ", cp, " vs ", rf)]] <- risco
    efe <- rbind(efe, data.frame(indice = c(paste0("RR (", ev, ": ", cp, " vs ", rf, ")"), paste0("DR (", ev, ": ", cp, " - ", rf, ")")),
                                 valor = c(rr, rd), ic_inf = c(rr_ic[1], rd_ic[1]), ic_sup = c(rr_ic[2], rd_ic[2]), principal = FALSE))
    decisoes <- c(decisoes, paste0("Risco: RR com IC log-normal; DR com IC de Newcombe (1998, híbrido de Wilson); NNT = 1/|DR| com IC invertido (Altman, 1998). ",
                                   "Em delineamento transversal o RR é razão de prevalências e o NNT não tem leitura de tratamento; NNT pede ensaio ou coorte."))
    if (!is.null(fisher) && !is.null(fisher$estimate))
      extra_tab[["OR condicional de Fisher (máxima verossimilhança condicional, IC exato)"]] <-
        data.frame(OR = unname(fisher$estimate), ic_inf = fisher$conf.int[1], ic_sup = fisher$conf.int[2],
                   nota = "orientação: linha 1 × coluna 1 da tabela; mesmo valor que o JASP mostra ao lado do OR")
    apa <- c(apa, paste0(
      if (sig) .ee_assoc_tr(paste0("Houve associação estatisticamente significativa entre ", vl, " e ", vc, ", "),
                            paste0("There was a statistically significant association between ", vl, " and ", vc, ", "), idioma) else
        .ee_assoc_tr(paste0("Não houve associação estatisticamente significativa entre ", vl, " e ", vc, ", "),
                     paste0("There was no statistically significant association between ", vl, " and ", vc, ", "), idioma),
      est_txt, .ee_sep(decimal), .ee_assoc_ef("φ", abs(phi), phi_ic[1], phi_ic[2], TRUE, casas, decimal, idioma, nivel),
      .ee_assoc_tr(paste0(" (efeito ", .ee_assoc_faixa_w(w, idioma), " pela régua de Cohen, 1988)."), paste0(" (", .ee_assoc_faixa_w(w, idioma), " effect by Cohen's, 1988, benchmarks)."), idioma)),
      paste0(.ee_assoc_tr(paste0("A chance de ", vc, " = ", ev, " em ", cp, " foi "), paste0("The odds of ", vc, " = ", ev, " in ", cp, " were "), idioma),
             .ee_assoc_n(or, casas, TRUE, decimal), .ee_assoc_tr(paste0(" vezes a de ", rf), paste0(" times those in ", rf), idioma),
             if (or < 2) paste0(" (", .pr_num(abs(or - 1) * 100, 1, TRUE, decimal), "% ",
                                if (or < 1) .ee_assoc_tr("menor", "lower", idioma) else .ee_assoc_tr("maior", "higher", idioma), ")") else "", ", ",
             .ee_assoc_ef("OR", or, or_ic[1], or_ic[2], FALSE, casas, decimal, idioma, nivel), haldane_frase,
             .ee_assoc_tr("; em probabilidades, ", "; as probabilities, ", idioma),
             .ee_assoc_pct(chances$probabilidade[2], 1, decimal), .ee_assoc_tr(" em ", " in ", idioma), cp, .ee_assoc_tr(" e ", " and ", idioma),
             .ee_assoc_pct(chances$probabilidade[1], 1, decimal), .ee_assoc_tr(" em ", " in ", idioma), rf, "."))
    decisoes <- c(decisoes,
                  paste0("OR = chance de ", vc, " = '", ev, "' (evento: ", regra_ev, ") em '", cp, "' dividida pela chance em '", rf, "' (referência: linha ", if (is.null(referencia)) "1" else "informada", "); IC de Woolf (log OR), que reproduz o livro (9.99 [7.67; 13.01]) e o JASP (0.100 [0.077; 0.130] na orientação inversa)",
                         if (corr > 0) paste0("; célula zero: OR, IC e OR inverso com correção de Haldane (+0,5 em todas as células), declarada também na frase APA; sem a correção o OR seria ", bruto_txt[1],
                                              ", e o OR condicional de Fisher (IC exato) fica ao lado como alternativa sem correção") else "", "."),
                  "Chance (parte/parte) não é probabilidade (parte/todo): as duas saem lado a lado. A chance de um grupo sozinha não é OR (o erro do livro em chama de OR a chance dos homens, 0,200; o OR homens/mulheres é 0,100, 90% menor).",
                  "φ com sinal pela ordem das categorias (como o JASP: negativo se ad − bc < 0); a magnitude |φ| é o que se compara às faixas .10/.30/.50. IC bilateral do φ pelo effectsize (o JASP não dá IC).")
  } else {
    apa <- c(apa, paste0(
      if (sig) .ee_assoc_tr(paste0("Houve associação estatisticamente significativa entre ", vl, " e ", vc, ", "),
                            paste0("There was a statistically significant association between ", vl, " and ", vc, ", "), idioma) else
        .ee_assoc_tr(paste0("Não houve associação estatisticamente significativa entre ", vl, " e ", vc, ", "),
                     paste0("There was no statistically significant association between ", vl, " and ", vc, ", "), idioma),
      est_txt, .ee_sep(decimal), .ee_assoc_ef(.ee_assoc_tr("V de Cramér", "Cramér's V", idioma), V, v_ic[1], v_ic[2], TRUE, casas, decimal, idioma, nivel),
      " (", if (k > 2) paste0("w = V·√", k - 1, " = ", .ee_assoc_n(w, casas, FALSE, decimal), ", ") else "",
      .ee_assoc_tr("efeito ", "", idioma), .ee_assoc_faixa_w(w, idioma),
      .ee_assoc_tr(" pela régua de Cohen, 1988).", " effect by Cohen's, 1988, benchmarks).", idioma)))
    efe <- rbind(efe, data.frame(indice = paste0("w = V·√(min(l, c) − 1), para comparar com .10/.30/.50"), valor = w,
                                 ic_inf = v_ic[1] * sqrt(k - 1), ic_sup = v_ic[2] * sqrt(k - 1), principal = FALSE))
    dest <- cel[cel[[paste0("acima_", formatC(z1, format = "f", digits = 2))]], , drop = FALSE]
    if (nrow(dest)) {
      itens <- paste0(dest[[vl]], "/", dest[[vc]], " (", ifelse(dest$residuo_ajustado > 0, "+", ""), .ee_assoc_tr("", "", idioma),
                      vapply(dest$residuo_ajustado, function(v) .pr_num(v, 2, TRUE, decimal), ""), ")")
      sobrevive <- all(dest$acima_bonferroni)
      apa <- c(apa, paste0(.ee_assoc_tr("Os resíduos ajustados localizaram a associação nas células ", "Adjusted residuals located the association in cells ", idioma),
                           paste(itens, collapse = .ee_sep(decimal)), " (|z| > ", .pr_num(z1, 2, TRUE, decimal), ")",
                           if (sobrevive) .ee_assoc_tr(paste0(", todas acima também do limiar de Bonferroni para ", m_cel, " células (|z| > ", .pr_num(zb, 2, TRUE, decimal), ")."),
                                                       paste0(", all also beyond the Bonferroni threshold for ", m_cel, " cells (|z| > ", .pr_num(zb, 2, TRUE, decimal), ")."), idioma) else
                             .ee_assoc_tr(paste0("; com o limiar de Bonferroni para ", m_cel, " células (|z| > ", .pr_num(zb, 2, TRUE, decimal), ") ficam ", sum(dest$acima_bonferroni), "."),
                                          paste0("; with the Bonferroni threshold for ", m_cel, " cells (|z| > ", .pr_num(zb, 2, TRUE, decimal), "), ", sum(dest$acima_bonferroni), " remain."), idioma)))
    }
    decisoes <- c(decisoes,
                  paste0("V de Cramér sem correção de viés (como o JASP), IC bilateral pelo qui-quadrado não central (effectsize). Antes de comparar às faixas de Cohen (.10/.30/.50, definidas para o 2×2), multiplicar V por √(min(l, c) − 1) (livro): o valor comparável é w = √(χ²/N) = ", round(w, 3), "."),
                  paste0("Pós-teste: resíduos padronizados ajustados (Haberman, 1973), limiar |z| > ", round(z1, 2), " como o livro e, ao lado, o limiar de Bonferroni para ", m_cel,
                         " células (|z| > ", round(zb, 2), "). Com muitas células o limiar sem correção não protege o conjunto (livro; Beasley & Schumacker, 1995)."))
  }
  decisoes <- c(paste0("Frequências esperadas: ", round(100 * esp$prop_menor5, 1), "% das células abaixo de 5; menor esperada = ", round(esp$minimo, 2),
                       if (esp$violada) ". Regra violada (no máximo 20% < 5 e nenhuma < 1; livro; Cochran, 1954)." else ". Regra atendida (no máximo 20% < 5 e nenhuma < 1; livro)."),
                paste0("p lido: ", testes$teste[lido],
                       if (usa_exato && exato == "auto") " (troca visível: a regra das esperadas não se sustenta)." else
                         if (dois && !correcao_continuidade) ". χ² sem correção de continuidade por padrão, como o JASP; Yates sai ao lado (conservador; Haviland, 1990; Campbell, 2007): declarar qual foi lido." else "."),
                decisoes,
                "Associação entre categóricas: falar em associação, não em correlação. Observações independentes: cada pessoa uma vez; medidas repetidas pedem mcnemar() ou q_cochran().")
  if (esp$violada && exato == "nao") avisos <- c(avisos, "Regra das frequências esperadas violada e exato = \"nao\": o p do χ² pode não ser confiável.")
  if (!dois && is.null(fisher) && esp$violada) avisos <- c(avisos, "Teste exato não calculado.")
  if (any(tab < 5))
    avisos <- c(avisos, paste0("Há célula(s) com contagem observada abaixo de 5 (menor = ", min(tab), "): a tabela de células mostra contagens. ",
                               "Cruzamento de variáveis demográficas com células < 5 não deve ser publicado sem agrupar categorias (risco de identificar participantes)."))
  ep_q <- efe[efe$principal, , drop = FALSE]
  avisos <- c(avisos, .ee_aviso_nao_sig(testes$p[lido], 1 - nivel,
                                        if (nrow(ep_q) && !is.na(ep_q$ic_inf[1])) paste0(ep_q$indice[1], " = ", .pr_num(ep_q$valor[1], 2, FALSE, decimal), .ee_sep(decimal),
                                                                                         .ee_ic(ep_q$ic_inf[1], ep_q$ic_sup[1], 2, FALSE, decimal, nivel, "pt")),
                                        "associação", NULL, decimal))
  tab_testes <- testes
  .ee_assoc_novo(analise = paste0("Qui-quadrado de independência (", vl, " × ", vc, ", ", nl, "×", nc, ")"), n = N, n_excluidos = tc$n_excluidos,
                 tabela = tab_testes, efeito = efe,
                 pressupostos = data.frame(checagem = "Frequências esperadas (≤ 20% abaixo de 5; nenhuma abaixo de 1)",
                                           estatistica = esp$minimo, gl = NA, p = NA,
                                           leitura = paste0(round(100 * esp$prop_menor5, 1), "% das células < 5; menor esperada = ", round(esp$minimo, 2),
                                                            if (esp$violada) ": violada; ler o teste exato" else ": atendida")),
                 decisoes = decisoes, avisos = avisos, apa = apa,
                 extra = list(tabela_contingencia = tab, esperados = E, residuos_ajustados = ct$stdres),
                 tabelas_extra = extra_tab)
}

qui_quadrado_aderencia <- function(dados, variavel, proporcoes = NULL, semente = 2026, simulacoes = 10000,
                                   nivel = .95, casas = 3, decimal = ".", idioma = "pt") {
  dados <- .ee_assoc_vazios(dados, variavel)
  dd <- .ee_dados(dados, variavel)
  x0 <- dd$d[[variavel]]
  obs <- unique(as.character(x0))
  # Na aderência a hipótese define as categorias: categoria esperada sem nenhuma
  # observação conta (contagem 0) e não pode sair por droplevels.
  if (!is.null(proporcoes) && !is.null(names(proporcoes))) {
    cats <- names(proporcoes)
    fora <- setdiff(obs, cats)
    if (length(fora)) stop("há categorias observadas sem proporção esperada: ", paste(fora, collapse = ", "))
    p <- proporcoes
  } else {
    cats <- if (is.factor(x0)) levels(x0) else sort(obs)
    if (is.null(proporcoes)) p <- rep(1 / length(cats), length(cats)) else {
      if (length(proporcoes) != length(cats))
        stop("proporcoes precisa ter uma entrada por categoria (", length(cats), ": ", paste(cats, collapse = ", "),
             "); para incluir categoria sem observações, use proporcoes nomeadas ou um fator com esse nível")
      p <- proporcoes
    }
  }
  if (length(cats) < 2) stop("a aderência precisa de pelo menos 2 categorias")
  p <- as.numeric(p) / sum(p)
  x <- table(factor(as.character(x0), levels = cats))
  N <- sum(x); vazias <- cats[as.vector(x) == 0]
  ct <- suppressWarnings(stats::chisq.test(x, p = p))
  esp <- .ee_assoc_esperados(ct$expected)
  testes <- data.frame(teste = "χ² de aderência", estatistica = unname(ct$statistic), gl = unname(ct$parameter), p = ct$p.value)
  if (esp$violada) {
    mc <- .ee_assoc_com_semente(semente, suppressWarnings(stats::chisq.test(x, p = p, simulate.p.value = TRUE, B = simulacoes)))
    testes <- rbind(testes, data.frame(teste = paste0("p por Monte Carlo (B = ", simulacoes, ", semente ", semente, ")"),
                                       estatistica = unname(mc$statistic), gl = NA, p = mc$p.value))
  }
  lido <- nrow(testes); testes$lido <- seq_len(nrow(testes)) == lido
  w <- sqrt(unname(ct$statistic) / N)
  .pr_exige("effectsize")
  w_ic <- tryCatch({ e <- effectsize::cohens_w(as.vector(x), p = p, ci = nivel, alternative = "two.sided"); c(e$CI_low, e$CI_high) },
                   error = function(err) c(NA_real_, NA_real_))
  cel <- data.frame(categoria = cats, observado = as.vector(x), esperado = as.vector(ct$expected),
                    prop_observada = as.vector(x) / N, prop_esperada = p, residuo_padronizado = as.vector(ct$stdres))
  sig <- testes$p[lido] < 1 - nivel   # alfa = 1 − nivel, como no resto do módulo
  apa <- paste0(if (sig) .ee_assoc_tr(paste0("A distribuição de ", variavel, " diferiu da esperada, "), paste0("The distribution of ", variavel, " differed from the expected one, "), idioma) else
    .ee_assoc_tr(paste0("A distribuição de ", variavel, " não diferiu da esperada, "), paste0("The distribution of ", variavel, " did not differ from the expected one, "), idioma),
    "χ²(", testes$gl[1], ", N = ", N, ") = ", .pr_num(testes$estatistica[1], 2, TRUE, decimal), .ee_sep(decimal), .pr_p(testes$p[lido], decimal),
    .ee_sep(decimal), .ee_assoc_ef("w", w, w_ic[1], w_ic[2], TRUE, casas, decimal, idioma, nivel), ".")
  .ee_assoc_novo(analise = paste0("Qui-quadrado de aderência (", variavel, ")"), n = N, n_excluidos = dd$n_excluidos,
                 tabela = testes, efeito = data.frame(indice = "w de Cohen = √(χ²/N)", valor = w, ic_inf = w_ic[1], ic_sup = w_ic[2], principal = TRUE),
                 pressupostos = data.frame(checagem = "Frequências esperadas (≤ 20% abaixo de 5; nenhuma abaixo de 1)", estatistica = esp$minimo, gl = NA, p = NA,
                                           leitura = paste0("menor esperada = ", round(esp$minimo, 2), if (esp$violada) ": violada; p lido por Monte Carlo" else ": atendida")),
                 decisoes = c(paste0("Distribuição esperada: ", if (is.null(proporcoes)) "proporções iguais" else "proporções informadas", " (", paste(paste0(cats, " = ", round(p, 3)), collapse = "; "), ")."),
                              if (length(vazias)) paste0("Categoria(s) sem nenhuma observação mantida(s) no teste com contagem 0: ", paste(vazias, collapse = ", "),
                                                         ". A hipótese define as categorias; tirá-las mudaria o teste.") else NULL,
                              "Tamanho de efeito w = √(χ²/N), índice de Cohen (1988) para qui-quadrado (.10/.30/.50), com IC bilateral (effectsize); o livro não traz efeito para a aderência.",
                              if (esp$violada) "Regra das esperadas violada: p por Monte Carlo com semente fixa (troca visível)." else NULL),
                 apa = apa, tabelas_extra = list(`Categorias` = cel))
}

# ================================================================ McNemar e Q de Cochran
# Medidas binárias pareadas: junta as colunas preservando o tipo, para o evento seguir a
# regra do módulo (fator -> 2º nível, na ordem de níveis da primeira medida; 0/1 -> 1;
# texto -> exige evento), e não a ordem alfabética.
.ee_assoc_binarias_pareadas <- function(d, vars, evento) {
  cols <- lapply(vars, function(v) { x <- d[[v]]; if (is.logical(x)) as.integer(x) else x }); names(cols) <- vars
  todos <- if (all(vapply(cols, is.factor, logical(1)))) {
    niv <- unique(unlist(lapply(cols, function(x) levels(droplevels(x)))))
    factor(unlist(lapply(cols, as.character)), levels = niv)
  } else if (all(vapply(cols, is.numeric, logical(1)))) unlist(cols) else unlist(lapply(cols, as.character))
  b <- .ee_assoc_binaria(todos, evento, paste(vars, collapse = "/"))
  m <- vapply(cols, function(x) as.integer(as.character(x) == b$evento), integer(nrow(d)))
  if (!is.matrix(m)) m <- matrix(m, nrow = nrow(d), dimnames = list(NULL, vars))
  list(m = m, b = b)
}

mcnemar <- function(dados, antes, depois, evento = NULL, teste = c("assintotico", "corrigido", "exato", "midp"),
                    nivel = .95, casas = 3, decimal = ".", idioma = "pt") {
  teste <- match.arg(teste)
  dados <- .ee_assoc_vazios(dados, c(antes, depois))
  dd <- .ee_dados(dados, c(antes, depois))
  bp <- .ee_assoc_binarias_pareadas(dd$d, c(antes, depois), evento)
  m <- bp$m; n <- nrow(m); ev <- bp$b$evento; ou <- bp$b$outro
  b <- sum(m[, 1] == 0 & m[, 2] == 1); c_ <- sum(m[, 1] == 1 & m[, 2] == 0)
  a <- sum(m[, 1] == 0 & m[, 2] == 0); d_ <- sum(m[, 1] == 1 & m[, 2] == 1)
  nd <- b + c_
  if (nd == 0) stop("nenhum par discordante (todos mantiveram a categoria): o teste de McNemar não se aplica")
  q <- (b - c_)^2 / nd; qc <- (max(abs(b - c_) - 1, 0))^2 / nd
  pe <- stats::binom.test(b, nd, .5)$p.value
  kmin <- min(b, c_); pmid <- min(1, 2 * (stats::pbinom(kmin, nd, .5) - .5 * stats::dbinom(kmin, nd, .5)))
  testes <- data.frame(teste = c("χ² de McNemar (sem correção)", "χ² de McNemar com correção de continuidade", "exato (binomial)", "mid-p exato"),
                       estatistica = c(q, qc, NA, NA), gl = c(1, 1, NA, NA), p = c(stats::pchisq(q, 1, lower.tail = FALSE), stats::pchisq(qc, 1, lower.tail = FALSE), pe, pmid))
  lido <- match(teste, c("assintotico", "corrigido", "exato", "midp")); testes$lido <- seq_len(4) == lido
  zc <- stats::qnorm(1 - (1 - nivel) / 2)
  or <- if (c_ > 0) b / c_ else Inf
  or_ic <- if (b > 0 && c_ > 0) exp(log(b / c_) + c(-1, 1) * zc * sqrt(1 / b + 1 / c_)) else c(NA, NA)
  p1 <- mean(m[, 1]); p2 <- mean(m[, 2]); dif <- p2 - p1
  se <- sqrt(nd - (b - c_)^2 / n) / n
  efe <- data.frame(indice = c(paste0("OR pareado (", ou, "→", ev, " / ", ev, "→", ou, ")"), paste0("diferença de proporções de '", ev, "' (", depois, " − ", antes, ")")),
                    valor = c(or, dif), ic_inf = c(or_ic[1], dif - zc * se), ic_sup = c(or_ic[2], dif + zc * se), principal = c(TRUE, TRUE))
  freq <- data.frame(antes = c(ou, ou, ev, ev), depois = c(ou, ev, ou, ev), n = c(a, b, c_, d_),
                     padrao = c("manteve", "mudou para o evento", "saiu do evento", "manteve"))
  names(freq)[1:2] <- c(antes, depois)
  sig <- testes$p[lido] < 1 - nivel
  est_txt <- if (lido <= 2) paste0("χ²(1, N = ", n, ") = ", .pr_num(testes$estatistica[lido], 2, TRUE, decimal), .ee_sep(decimal), .pr_p(testes$p[lido], decimal)) else
    paste0(testes$teste[lido], ", ", .pr_p(testes$p[lido], decimal))
  apa <- paste0(.ee_assoc_tr("Teste de McNemar: ", "McNemar test: ", idioma), est_txt, ". ",
                .ee_assoc_tr(paste0("A proporção de '", ev, "' passou de ", .ee_assoc_pct(p1, 1, decimal), " (", antes, ") para ", .ee_assoc_pct(p2, 1, decimal), " (", depois, "); dos ", nd, " pares que mudaram, ",
                                    b, " foram de '", ou, "' para '", ev, "' e ", c_, " de '", ev, "' para '", ou, "'"),
                             paste0("The proportion of '", ev, "' went from ", .ee_assoc_pct(p1, 1, decimal), " (", antes, ") to ", .ee_assoc_pct(p2, 1, decimal), " (", depois, "); of the ", nd, " pairs that changed, ",
                                    b, " went from '", ou, "' to '", ev, "' and ", c_, " from '", ev, "' to '", ou, "'"), idioma),
                if (is.finite(or)) paste0(" (", .ee_assoc_ef(.ee_assoc_tr("OR pareado", "paired OR", idioma), or, or_ic[1], or_ic[2], FALSE, casas, decimal, idioma, nivel), ")") else "", ".")
  avisos <- if (nd < 25) paste0("Só ", nd, " pares discordantes: a aproximação do χ² fica frágil; preferir teste = \"exato\" ou \"midp\" (ambos já na tabela).") else character()
  .ee_assoc_novo(analise = paste0("Teste de McNemar (", antes, " → ", depois, ")"), n = n, n_excluidos = dd$n_excluidos,
                 tabela = testes, efeito = efe, decisoes = c(bp$b$frase,
                   "Só os pares discordantes entram no teste; χ² = (b − c)²/(b + c), gl = 1.",
                   paste0("p lido: ", testes$teste[lido], ". Os quatro (sem correção, com correção, exato binomial e mid-p) saem sempre lado a lado; o SPSS aplica correção por padrão, por isso os valores podem diferir."),
                   "O teste diz se houve mudança; a direção se lê na tabela de pares (quem mudou para o evento vs quem saiu dele).",
                   "Efeitos: OR pareado = b/c (IC pelo log) e diferença entre as proporções dos dois momentos (IC de Wald para pares). O curso não traz efeito para o McNemar."),
                 avisos = avisos, apa = apa, tabelas_extra = list(`Pares` = freq))
}

q_cochran <- function(dados, vars, evento = NULL, ajuste = c("bonferroni", "holm"), nivel = .95,
                      decimal = ".", idioma = "pt") {
  ajuste <- match.arg(ajuste)
  if (length(vars) < 2) stop("informe pelo menos 2 medidas")
  dados <- .ee_assoc_vazios(dados, vars)
  dd <- .ee_dados(dados, vars)
  bp <- .ee_assoc_binarias_pareadas(dd$d, vars, evento)
  m <- bp$m; n <- nrow(m); k <- ncol(m); ev <- bp$b$evento
  Cj <- colSums(m); Ri <- rowSums(m); Tt <- sum(m)
  den <- k * Tt - sum(Ri^2)
  if (den == 0) stop("sem variação dentro dos participantes (todos com o mesmo valor em todas as medidas): o Q de Cochran não se aplica")
  Q <- (k - 1) * (k * sum(Cj^2) - Tt^2) / den
  pQ <- stats::pchisq(Q, k - 1, lower.tail = FALSE)
  props <- data.frame(medida = vars, n = n, com_evento = Cj, proporcao = Cj / n, row.names = NULL)
  pares <- utils::combn(seq_len(k), 2, simplify = FALSE)
  se_z <- sqrt(2 * (k * Tt - sum(Ri^2)) / (n^2 * k * (k - 1)))
  ph <- do.call(rbind, lapply(pares, function(pr) {
    i <- pr[1]; j <- pr[2]
    b <- sum(m[, i] == 0 & m[, j] == 1); c_ <- sum(m[, i] == 1 & m[, j] == 0)
    q <- if (b + c_ > 0) (b - c_)^2 / (b + c_) else 0
    z <- (Cj[j] - Cj[i]) / n / se_z
    data.frame(comparacao = paste0(vars[i], " vs ", vars[j]), dif_proporcao = (Cj[j] - Cj[i]) / n, b = b, c = c_,
               qui2_mcnemar = q, p_mcnemar = stats::pchisq(q, 1, lower.tail = FALSE), z_spss = unname(z),
               p_z = 2 * stats::pnorm(-abs(unname(z))))
  }))
  ph$p_mcnemar_ajustado <- stats::p.adjust(ph$p_mcnemar, ajuste)
  ph$p_z_ajustado <- stats::p.adjust(ph$p_z, ajuste)
  sig <- pQ < 1 - nivel
  apa <- c(paste0(.ee_assoc_tr("Q de Cochran: ", "Cochran's Q: ", idioma), "χ²(", k - 1, ", N = ", n, ") = ", .pr_num(Q, 2, TRUE, decimal), .ee_sep(decimal), .pr_p(pQ, decimal), ". ",
                  .ee_assoc_tr(paste0("Proporções de '", ev, "': "), paste0("Proportions of '", ev, "': "), idioma),
                  paste(paste0(vars, " = ", vapply(Cj / n, function(v) .ee_assoc_pct(v, 1, decimal), "")), collapse = .ee_sep(decimal)), "."))
  if (sig && k > 2) {
    rot_p <- .ee_assoc_tr("p ajustado", "adjusted p", idioma)
    itens <- vapply(seq_len(nrow(ph)), function(i)
      paste0(ph$comparacao[i], " (", .ee_estat("χ²", ph$qui2_mcnemar[i], 1, decimal), .ee_sep(decimal),
             sub("^p", rot_p, .pr_p(ph$p_mcnemar_ajustado[i], decimal)), ")"), "")
    lista <- if (length(itens) > 1) paste0(paste(itens[-length(itens)], collapse = ", "), .ee_assoc_tr(" e ", " and ", idioma), itens[length(itens)]) else itens
    apa <- c(apa, paste0(.ee_assoc_tr(paste0("Comparações par a par (McNemar, p ajustado por ", ajuste, "): "), paste0("Pairwise comparisons (McNemar, ", ajuste, "-adjusted p): "), idioma),
                         lista, "."))
  }
  .ee_assoc_novo(analise = paste0("Q de Cochran (", k, " medidas)"), n = n, n_excluidos = dd$n_excluidos,
                 tabela = data.frame(teste = "Q de Cochran", estatistica = Q, gl = k - 1, p = pQ),
                 efeito = data.frame(indice = paste0("proporção de '", ev, "' em ", vars), valor = Cj / n, ic_inf = NA, ic_sup = NA, principal = TRUE),
                 posthoc = ph,
                 decisoes = c(bp$b$frase, "Mesmos participantes em todas as medidas (listwise); quem falta em alguma medida sai do teste.",
                              paste0("Pós-hoc: McNemar par a par (sem correção de continuidade) com p ajustado por ", ajuste,
                                     ". Ao lado, o z das comparações par a par do SPSS (diferença de proporções / EP comum), que o curso usa; os dois podem diferir."),
                              "Ler pré-pós (efeito imediato), pós-seguimento (manutenção ou recaída) e pré-seguimento (efeito duradouro) pela direção das proporções.",
                              "O livro e o curso não trazem tamanho de efeito para o Q; as proporções por medida e as diferenças par a par fazem esse papel."),
                 apa = apa, tabelas_extra = list(`Proporções por medida` = props))
}

# ================================================================ viés preditivo (Cleary)
# Cleary (1968): o teste tem viés preditivo para um grupo quando a reta que prevê o
# critério a partir do escore difere entre os grupos, em inclinação, intercepto ou ambos
# (Petersen, 2024, cap. 15). Distinto do viés de medida (invariância, DIF): aqui a
# pergunta é se o mesmo escore prevê o mesmo critério nos grupos. Sequência step-down de
# Lautenschlager e Mendoza (1986):
#   1. ômnibus: critério ~ escore contra critério ~ escore * grupo;
#   2. inclinações: critério ~ escore + grupo contra critério ~ escore * grupo;
#   3. interceptos, só quando as inclinações não diferem: critério ~ escore contra
#      critério ~ escore + grupo.
# O escore entra centrado na média geral: a diferença de intercepto é a diferença de
# critério prevista no escore médio. A tabela extra traz a reta de cada grupo (na escala
# original) e o resíduo médio de cada grupo em relação à reta comum: positivo, a reta
# comum subestima o grupo; negativo, superestima. Poder para inclinações costuma ser
# baixo (Aguinis, Culpepper & Pierce, 2010): não rejeitar não prova ausência de viés.
vies_preditivo <- function(dados, criterio, preditor, grupo, referencia = NULL, alfa = .05, nivel = .95, casas = 3,
                           decimal = ".", idioma = "pt") {
  dados <- .ee_assoc_vazios(dados, c(criterio, preditor, grupo))
  dd <- .ee_dados(dados, c(criterio, preditor, grupo), numericas = c(criterio, preditor))
  d <- dd$d
  g <- factor(d[[grupo]])
  if (!is.null(referencia)) { if (!referencia %in% levels(g)) stop("referencia não é nível de ", grupo); g <- stats::relevel(g, referencia) }
  if (nlevels(g) < 2) stop("o grupo precisa de pelo menos 2 níveis com dados")
  ng <- table(g)
  if (any(ng < 10)) stop("cada grupo precisa de pelo menos 10 casos com critério e escore (menor: ", min(ng), ")")
  y <- d[[criterio]]; x <- d[[preditor]]; xc <- x - mean(x)
  m0 <- stats::lm(y ~ xc); m1 <- stats::lm(y ~ xc + g); m2 <- stats::lm(y ~ xc * g)
  r2 <- function(m) summary(m)$r.squared
  passo <- function(a, b, nome) {
    av <- stats::anova(a, b)
    data.frame(teste = nome, F = av$F[2], gl1 = av$Df[2], gl2 = av$Res.Df[2], p = av$`Pr(>F)`[2],
               delta_r2 = r2(b) - r2(a), f2 = (r2(b) - r2(a)) / (1 - r2(b)), stringsAsFactors = FALSE)
  }
  tab <- rbind(passo(m0, m2, "1. ômnibus (inclinações e interceptos)"), passo(m1, m2, "2. inclinações"), passo(m0, m1, "3. interceptos"))
  omni <- tab$p[1] < alfa; incl <- tab$p[2] < alfa; inter <- tab$p[3] < alfa
  retas <- do.call(rbind, lapply(levels(g), function(k) {
    i <- g == k; mk <- stats::lm(y[i] ~ x[i]); ci <- stats::confint(mk, level = nivel)
    data.frame(grupo = k, n = sum(i), intercepto = unname(stats::coef(mk)[1]), inclinacao = unname(stats::coef(mk)[2]),
               incl_ic_inf = ci[2, 1], incl_ic_sup = ci[2, 2], r = stats::cor(x[i], y[i]),
               residuo_medio_reta_comum = mean(stats::resid(m0)[i]), dp_residual = stats::sigma(mk), stringsAsFactors = FALSE)
  }))
  conclusao <- if (!omni) "sem_evidencia" else if (incl) "inclinacoes" else if (inter) "interceptos" else "omnibus_sem_componente"
  avisos <- character()
  razao_var <- max(retas$dp_residual^2) / min(retas$dp_residual^2)
  if (razao_var > 1.5) avisos <- c(avisos, paste0("Variância residual difere entre grupos (razão maior/menor = ", .pr_num(razao_var, 2, TRUE, decimal),
                                                  "): a regressão moderada pressupõe variância de erro homogênea; com grupos de tamanhos diferentes, o teste das inclinações fica distorcido."))
  if (!incl) avisos <- c(avisos, "Inclinações não diferiram, mas o teste costuma ter pouco poder (Aguinis, Culpepper & Pierce, 2010); relatar a diferença estimada e o IC das inclinações.")
  if (conclusao == "interceptos") {
    sob <- retas$grupo[retas$residuo_medio_reta_comum > 0]; sobre <- retas$grupo[retas$residuo_medio_reta_comum < 0]
    avisos <- c(avisos, paste0("Interceptos diferem: a reta comum subestima o critério de ", paste(sob, collapse = ", "),
                               " e superestima o de ", paste(sobre, collapse = ", "), " (resíduo médio na tabela de retas)."))
  }
  frase_res <- c(sem_evidencia = "não houve evidência de viés preditivo", inclinacoes = "as inclinações diferiram entre os grupos (validade diferencial)",
                 interceptos = "as inclinações não diferiram, mas os interceptos sim (o mesmo escore prevê critério diferente)",
                 omnibus_sem_componente = "o teste ômnibus rejeitou, sem que inclinações ou interceptos isoladamente o fizessem")[[conclusao]]
  frase_res_en <- c(sem_evidencia = "there was no evidence of predictive bias", inclinacoes = "slopes differed across groups (differential validity)",
                    interceptos = "slopes did not differ, but intercepts did (the same score predicts different criterion levels)",
                    omnibus_sem_componente = "the omnibus test was significant, but neither slopes nor intercepts alone were")[[conclusao]]
  fr <- function(i) paste0("F(", tab$gl1[i], ", ", tab$gl2[i], ") = ", .pr_num(tab$F[i], 2, TRUE, decimal), ", ", .pr_p(tab$p[i], decimal),
                           ", ΔR² = ", .pr_num(tab$delta_r2[i], 3, FALSE, decimal))
  apa <- if (identical(idioma, "en"))
    paste0("Predictive bias of ", preditor, " for ", criterio, " across ", grupo, " was tested with step-down moderated regression (Cleary, 1968; Lautenschlager & Mendoza, 1986): ",
           "omnibus ", fr(1), "; slopes ", fr(2), if (!incl) paste0("; intercepts ", fr(3)) else "", ". Thus, ", frase_res_en, ".") else
    paste0("O viés preditivo de ", preditor, " para ", criterio, " entre os níveis de ", grupo, " foi testado por regressão moderada em etapas (Cleary, 1968; Lautenschlager & Mendoza, 1986): ",
           "ômnibus ", fr(1), "; inclinações ", fr(2), if (!incl) paste0("; interceptos ", fr(3)) else "", ". Assim, ", frase_res, ".")
  .ee_assoc_novo(analise = paste0("Viés preditivo de Cleary (", preditor, " → ", criterio, ", por ", grupo, ")"), n = nrow(d), n_excluidos = dd$n_excluidos,
                 tabela = tab,
                 efeito = data.frame(indice = c("ΔR² ômnibus", "ΔR² inclinações", "ΔR² interceptos"), valor = tab$delta_r2,
                                     ic_inf = NA_real_, ic_sup = NA_real_, principal = c(TRUE, FALSE, FALSE)),
                 decisoes = c(paste0("Escore centrado na média geral (", .pr_num(mean(x), 2, TRUE, decimal), "); referência do grupo: ", levels(g)[1], "."),
                              "Sequência step-down: ômnibus; se rejeitar, inclinações; se as inclinações não diferem, interceptos. Com inclinações diferentes, a diferença de intercepto depende do ponto do escore e não é interpretada isolada.",
                              "Viés preditivo não é viés de medida: um teste pode ser invariante (sem DIF) e ainda prever o critério de forma diferente, e o contrário; ver references/invariancia.md."),
                 avisos = avisos, apa = apa, extra = list(conclusao = conclusao, modelos = list(comum = m0, interceptos = m1, completo = m2)),
                 tabelas_extra = list(`Reta por grupo (escala original)` = retas))
}
