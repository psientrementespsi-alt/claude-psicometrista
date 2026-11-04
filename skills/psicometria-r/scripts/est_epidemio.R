# est_epidemio.R: análises de estudos epidemiológicos (v1.4; joinpoint e padronização na v1.5).
# Carregar por carregar.R.
#
#   razao_prevalencia()  Poisson com erro padrão robusto (sanduíche) para desfecho
#                        binário: razão de prevalência (Zou, 2004; Barros & Hirakata,
#                        2003), a do curso de estudos epidemiológicos
#   sobrevida()          Kaplan-Meier (mediana e tábua nos tempos pedidos), log-rank e
#                        regressão de Cox com checagem de riscos proporcionais
#
# Desfecho binário em estudo transversal ou de coorte com desfecho comum: a razão de
# chances da logística superestima a razão de prevalência quando a prevalência passa de
# uns 10%; a Poisson com ligação log e erro robusto estima a RP direto. O EP robusto é o
# sanduíche HC0, o mesmo da GEE com correlação independente e um caso por cluster.

razao_prevalencia <- function(dados, desfecho, preditores, evento = NULL, referencias = NULL, fatores = NULL, nivel = .95,
                              decimal = ".", idioma = "pt") {
  dados <- .ee_assoc_vazios(dados, c(desfecho, preditores))
  dd <- .ee_dados(dados, c(desfecho, preditores))
  b <- .ee_assoc_binaria(dd$d[[desfecho]], evento, desfecho)
  prep <- .ee_assoc_preparar(dd$d, preditores, referencias, fatores)
  d <- prep$d; d$.y_rp <- b$y
  f <- stats::as.formula(paste(".y_rp ~", paste(.ee_assoc_bq(preditores), collapse = " + ")))
  fit <- suppressWarnings(stats::glm(f, family = stats::poisson(link = "log"), data = d))
  X <- stats::model.matrix(fit); mu <- stats::fitted(fit); e <- d$.y_rp - mu
  B <- solve(crossprod(X * sqrt(mu)))
  M <- crossprod(X * e)
  V <- B %*% M %*% B
  est <- stats::coef(fit); ep <- sqrt(diag(V)); z <- stats::qnorm(1 - (1 - nivel) / 2)
  tab <- data.frame(termo = names(est), b = unname(est), ep_robusto = unname(ep), z = unname(est / ep),
                    p = unname(2 * stats::pnorm(-abs(est / ep))), rp = unname(exp(est)),
                    ic_inf = unname(exp(est - z * ep)), ic_sup = unname(exp(est + z * ep)), stringsAsFactors = FALSE)
  avisos <- prep$avisos
  if (any(mu > 1)) avisos <- c(avisos, paste0(sum(mu > 1), " prevalência(s) prevista(s) acima de 1 (limitação conhecida da Poisson para desfecho binário): ",
                                              "a RP segue válida, mas a previsão individual não; o binomial com ligação log é a alternativa."))
  ef <- tab[-1, ]
  efeito <- data.frame(indice = paste0("RP ", ef$termo), valor = ef$rp, ic_inf = ef$ic_inf, ic_sup = ef$ic_sup,
                       principal = seq_len(nrow(ef)) == 1, stringsAsFactors = FALSE)
  frase <- vapply(seq_len(nrow(ef)), function(i) paste0(ef$termo[i], ": RP = ", .pr_num(ef$rp[i], 2, TRUE, decimal), ", ",
                                                       .ee_ic(ef$ic_inf[i], ef$ic_sup[i], 2, TRUE, decimal, nivel, idioma), ", ", .pr_p(ef$p[i], decimal)), "")
  apa <- if (identical(idioma, "en"))
    paste0("Prevalence ratios for ", desfecho, " = '", b$evento, "' were estimated by Poisson regression with robust (sandwich) standard errors (Zou, 2004): ",
           paste(frase, collapse = "; "), ".") else
    paste0("As razões de prevalência de ", desfecho, " = '", b$evento, "' foram estimadas por regressão de Poisson com erro padrão robusto (sanduíche; Zou, 2004): ",
           paste(frase, collapse = "; "), ".")
  .ee_resultado(paste0("Razão de prevalência (Poisson robusta): ", desfecho), n = nrow(d), n_excluidos = dd$n_excluidos, tabela = tab, efeito = efeito,
                decisoes = c(b$frase, paste0("Poisson com ligação log e EP sanduíche HC0 (igual à GEE independente com um caso por cluster); IC de Wald ",
                                             round(100 * nivel), "% exponenciado."),
                             paste0("Categóricos: ", paste(.ee_assoc_txt_fatores(prep$info), collapse = "; ")),
                             "RP, não RC: com desfecho comum a razão de chances da logística se afasta da RP (Barros & Hirakata, 2003)."),
                avisos = avisos, apa = apa, extra = list(modelo = fit, vcov_robusta = V))
}

# Sobrevida. tempo: tempo até o evento ou até a censura; evento: 1 (ou TRUE, ou o nível
# indicado em `evento_codigo`) = evento, o resto = censura. grupo: comparação de curvas
# (log-rank) e termo do Cox; covariaveis entram só no Cox. IC da sobrevida pelo padrão
# do survival (transformação log); a mediana com IC sai do survfit.
sobrevida <- function(dados, tempo, evento, grupo = NULL, covariaveis = NULL, evento_codigo = NULL, tempos = NULL,
                      referencias = NULL, fatores = NULL, nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("survival")
  vars <- c(tempo, evento, grupo, covariaveis)
  dados <- .ee_assoc_vazios(dados, vars)
  dd <- .ee_dados(dados, vars, numericas = tempo)
  d <- dd$d
  ev <- d[[evento]]
  st <- if (!is.null(evento_codigo)) as.integer(ev == evento_codigo) else if (is.logical(ev)) as.integer(ev) else {
    if (!all(ev %in% c(0, 1))) stop("evento precisa ser 0/1, lógico, ou informe evento_codigo"); as.integer(ev) }
  if (any(d[[tempo]] < 0)) stop("tempo negativo no banco")
  if (sum(st) < 1) stop("nenhum evento observado")
  d$.t <- d[[tempo]]; d$.s <- st
  preds <- c(grupo, covariaveis)
  if (length(preds)) { prep <- .ee_assoc_preparar(d, preds, referencias, fatores); d <- prep$d } else prep <- list(info = list(), avisos = character())
  sem_contraste <- function(expr) withCallingHandlers(expr, warning = function(w)
    if (grepl("contrasts dropped", conditionMessage(w))) invokeRestart("muffleWarning"))
  fk <- sem_contraste(if (is.null(grupo)) survival::survfit(survival::Surv(.t, .s) ~ 1, data = d, conf.int = nivel) else
    survival::survfit(stats::as.formula(paste("survival::Surv(.t, .s) ~", .ee_assoc_bq(grupo))), data = d, conf.int = nivel))
  tb <- summary(fk)$table
  if (is.null(dim(tb))) tb <- t(as.matrix(tb))
  lim <- paste0(c("0.", ""), round(100 * nivel), c("LCL", "UCL"))
  med <- data.frame(estrato = if (is.null(grupo)) "todos" else sub("^[^=]*=", "", rownames(tb)), n = tb[, "records"], eventos = tb[, "events"],
                    mediana = tb[, "median"], mediana_ic_inf = tb[, grep("LCL$", colnames(tb))[1]], mediana_ic_sup = tb[, grep("UCL$", colnames(tb))[1]],
                    row.names = NULL, stringsAsFactors = FALSE)
  extras <- list(`Mediana de sobrevida (Kaplan-Meier)` = med)
  if (length(tempos)) {
    sm <- summary(fk, times = tempos, extend = TRUE)
    extras$`Tábua de sobrevida nos tempos pedidos` <- data.frame(estrato = if (is.null(grupo)) "todos" else sub("^[^=]*=", "", as.character(sm$strata)),
                                                                 tempo = sm$time, em_risco = sm$n.risk, eventos = sm$n.event,
                                                                 sobrevida = sm$surv, ic_inf = sm$lower, ic_sup = sm$upper, stringsAsFactors = FALSE)
  }
  tab <- NULL; efeito <- NULL; decisoes <- character(); avisos <- prep$avisos; apa <- character()
  if (!is.null(grupo)) {
    lr <- sem_contraste(survival::survdiff(stats::as.formula(paste("survival::Surv(.t, .s) ~", .ee_assoc_bq(grupo))), data = d))
    glr <- length(lr$n) - 1; plr <- stats::pchisq(lr$chisq, glr, lower.tail = FALSE)
    extras$`Log-rank` <- data.frame(qui2 = lr$chisq, gl = glr, p = plr)
    decisoes <- c(decisoes, "Log-rank (Mantel-Cox) compara as curvas inteiras; pesa igual todos os tempos.")
  }
  if (length(preds)) {
    fc <- survival::coxph(stats::as.formula(paste("survival::Surv(.t, .s) ~", paste(.ee_assoc_bq(preds), collapse = " + "))), data = d)
    sc <- summary(fc, conf.int = nivel)
    cf <- sc$coefficients; ci <- sc$conf.int
    tab <- data.frame(termo = rownames(cf), b = cf[, "coef"], ep = cf[, "se(coef)"], z = cf[, "z"], p = cf[, "Pr(>|z|)"],
                      hr = ci[, "exp(coef)"], ic_inf = ci[, 3], ic_sup = ci[, 4], row.names = NULL, stringsAsFactors = FALSE)
    efeito <- data.frame(indice = paste0("HR ", tab$termo), valor = tab$hr, ic_inf = tab$ic_inf, ic_sup = tab$ic_sup,
                         principal = seq_len(nrow(tab)) == 1, stringsAsFactors = FALSE)
    zph <- survival::cox.zph(fc)
    zt <- as.data.frame(zph$table); zt$termo <- rownames(zph$table); rownames(zt) <- NULL
    extras$`Riscos proporcionais (Schoenfeld)` <- zt[, c("termo", "chisq", "df", "p")]
    if (any(zt$p < .05)) avisos <- c(avisos, paste0("Resíduos de Schoenfeld indicam violação de riscos proporcionais em ", paste(zt$termo[zt$p < .05], collapse = ", "),
                                                    ": HR médio no tempo; considerar estratificar ou interação com o tempo."))
    epv <- sum(st) / length(stats::coef(fc))
    if (epv < 10) avisos <- c(avisos, paste0("Eventos por parâmetro = ", .pr_num(epv, 1, TRUE, decimal), " (< 10): estimativas do Cox instáveis; reduzir covariáveis."))
    decisoes <- c(decisoes, paste0("Cox com empates pelo método de Efron (padrão do survival); HR com IC de Wald ", round(100 * nivel), "%; riscos proporcionais pelos resíduos de Schoenfeld (cox.zph)."))
    frase <- vapply(seq_len(nrow(tab)), function(i) paste0(tab$termo[i], ": HR = ", .pr_num(tab$hr[i], 2, TRUE, decimal), ", ",
                                                          .ee_ic(tab$ic_inf[i], tab$ic_sup[i], 2, TRUE, decimal, nivel, idioma), ", ", .pr_p(tab$p[i], decimal)), "")
    apa <- if (identical(idioma, "en")) paste0("In the Cox proportional hazards model, ", paste(frase, collapse = "; "), ".") else
      paste0("No modelo de riscos proporcionais de Cox, ", paste(frase, collapse = "; "), ".")
  }
  if (!is.null(grupo)) {
    lrt <- extras$`Log-rank`
    apa <- c(paste0(if (identical(idioma, "en")) (if (lrt$p < .05) "Kaplan-Meier curves differed across " else "Kaplan-Meier curves did not differ across ") else
                      "As curvas de Kaplan-Meier de ", grupo,
                    if (identical(idioma, "en")) " (log-rank test), " else if (lrt$p < .05) " diferiram pelo teste log-rank, " else " não diferiram pelo teste log-rank, ",
                    "χ²(", lrt$gl, ") = ", .pr_num(lrt$qui2, 2, TRUE, decimal), ", ", .pr_p(lrt$p, decimal), "."), apa)
  }
  decisoes <- c(paste0("Evento = ", if (!is.null(evento_codigo)) paste0(evento, " == ", evento_codigo) else paste0(evento, " = 1"), "; o resto é censura. ",
                       sum(st), " eventos em ", nrow(d), " casos."),
                "Kaplan-Meier: mediana com IC pela transformação log (padrão do survival); mediana NA quando a curva não cai abaixo de 50%.", decisoes)
  .ee_resultado("Análise de sobrevida", n = nrow(d), n_excluidos = dd$n_excluidos, tabela = tab, efeito = efeito, decisoes = decisoes, avisos = avisos,
                apa = apa, extra = list(km = fk, cox = if (length(preds)) fc else NULL, tabelas_extra = extras))
}

# ---------------------------------------------------------------- v1.5: tendência temporal
# tendencia_joinpoint(): regressão por pontos de inflexão (joinpoint; Kim et al., 2000) em
# R, a análise que o curso de séries temporais faz no programa Joinpoint do NCI. Modelo
# log-linear contínuo: log(taxa) = b0 + b1 t + Σ δj (t - τj)+, ajustado por mínimos
# quadrados com τj procurados nos anos observados (grade do Joinpoint), com pelo menos
# `min_extremo` observações entre um ponto e a ponta da série e `min_entre` entre dois
# pontos (ponto excluído; padrão 2 e 2, como no Joinpoint). Número de pontos pelo BIC do
# Joinpoint: ln(SQR/n) + 2(k + 1) ln(n)/n. Máximo de pontos pela regra do Joinpoint (até
# 11 observações, 1; 12 a 16, 2; 17 a 21, 3; 22 a 26, 4; mais, 5).
# APC = 100 (e^b - 1) por segmento; AAPC = média das inclinações ponderada pelo comprimento
# dos segmentos no período. IC paramétrico (t com n - 2k - 2 gl) ou por quantis empíricos
# de reamostragem dos resíduos, reajustando a localização dos pontos (ic = "empirico").
# Validação (evals_latentes_epidemio.R): com a incidência de tuberculose 2010-2019 da aula,
# os pontos (2015 e 2016), as APC e as AAPC batem com a saída do Joinpoint 4.9.1 no relato
# da aula; os IC do Joinpoint saíram mais largos que os dois métodos daqui, e a
# configuração da sessão (método do IC, variância) não é legível no arquivo .jpo. Para IC
# idêntico ao do programa, rodar o Joinpoint; aqui, declarar o método.
.jp_ajusta <- function(t, y, tau, w) {
  X <- if (length(tau)) cbind(1, t, vapply(tau, function(s) pmax(t - s, 0), numeric(length(t)))) else cbind(1, t)
  m <- stats::lm.wfit(X, y, w)
  list(sse = sum(w * m$residuals^2), coef = m$coefficients, X = X, ajustado = m$fitted.values, residuos = m$residuals, tau = tau)
}
.jp_candidatos <- function(t, k, me, mb) {
  n <- length(t)
  if (k == 0) return(list(numeric(0)))
  idx <- utils::combn(seq_len(n), k, simplify = FALSE)
  ok <- vapply(idx, function(i) i[1] - 1 >= me && n - i[k] >= me && (k == 1 || all(diff(i) - 1 >= mb)), TRUE)
  lapply(idx[ok], function(i) t[i])
}
.jp_melhor <- function(t, y, k, w, me, mb) {
  cs <- .jp_candidatos(t, k, me, mb)
  if (!length(cs)) return(NULL)
  fs <- lapply(cs, function(tau) .jp_ajusta(t, y, tau, w))
  fs[[which.min(vapply(fs, `[[`, 0, "sse"))]]
}
.jp_inclinacoes <- function(f) { b <- f$coef; cumsum(c(b[2], if (length(f$tau)) b[3:(2 + length(f$tau))])) }

tendencia_joinpoint <- function(dados, tempo, casos = NULL, populacao = NULL, taxa = NULL, por = 1e5, max_pontos = NULL,
                                min_extremo = 2, min_entre = 2, pesos = c("constante", "poisson"), ic = c("parametrico", "empirico"),
                                reamostras = 4999, semente = 2026, nivel = .95, decimal = ".", idioma = "pt") {
  pesos <- match.arg(pesos); ic <- match.arg(ic)
  if (is.null(taxa) && (is.null(casos) || is.null(populacao))) stop("informe taxa, ou casos e populacao")
  if (pesos == "poisson" && is.null(casos)) stop("pesos = 'poisson' pede casos")
  d <- dados[order(dados[[tempo]]), ]
  t <- as.numeric(d[[tempo]]); n <- length(t)
  if (anyDuplicated(t)) stop("tempo repetido: uma linha por período (para estratos, rodar por estrato)")
  tx <- if (!is.null(taxa)) d[[taxa]] else d[[casos]] / d[[populacao]] * por
  if (any(tx <= 0)) stop("taxa zero ou negativa: o modelo log-linear não se aplica")
  y <- log(tx); w <- if (pesos == "poisson") d[[casos]] else rep(1, n)
  if (is.null(max_pontos)) max_pontos <- if (n <= 11) 1 else if (n <= 16) 2 else if (n <= 21) 3 else if (n <= 26) 4 else 5
  fits <- list(); bic <- c()
  for (k in 0:max_pontos) { f <- .jp_melhor(t, y, k, w, min_extremo, min_entre); if (is.null(f)) break
    fits[[k + 1]] <- f; bic[k + 1] <- log(f$sse / n) + 2 * (k + 1) * log(n) / n }
  k <- which.min(bic) - 1; f <- fits[[k + 1]]
  gl <- n - 2 * k - 2
  sl <- .jp_inclinacoes(f)
  limites <- c(t[1], f$tau, t[n])
  comp <- diff(limites)
  W <- diag(w); s2 <- f$sse / gl
  V <- s2 * solve(t(f$X) %*% W %*% f$X)
  L <- t(vapply(seq_along(sl), function(i) c(0, 1, rep(1, i - 1), rep(0, k - i + 1)), numeric(2 + k)))
  if (k == 0) L <- matrix(c(0, 1), 1)
  se <- sqrt(diag(L %*% V %*% t(L)))
  la <- colSums(L * comp) / sum(comp); a <- sum(la * f$coef); sea <- sqrt(drop(t(la) %*% V %*% la))
  q <- stats::qt(1 - (1 - nivel) / 2, gl)
  lo <- sl - q * se; hi <- sl + q * se; alo <- a - q * sea; ahi <- a + q * sea
  p <- 2 * stats::pt(-abs(sl / se), gl); pa <- 2 * stats::pt(-abs(a / sea), gl)
  if (ic == "empirico") {
    set.seed(semente)
    fator <- sqrt(n / gl); res <- f$residuos * fator
    bt <- t(vapply(seq_len(reamostras), function(i) {
      ys <- f$ajustado + sample(res, n, TRUE); fb <- .jp_melhor(t, ys, k, w, min_extremo, min_entre)
      sb <- .jp_inclinacoes(fb); cb <- diff(c(t[1], fb$tau, t[n])); c(sb, sum(sb * cb) / sum(cb)) }, numeric(k + 2)))
    qq <- apply(bt, 2, stats::quantile, probs = c((1 - nivel) / 2, 1 - (1 - nivel) / 2))
    lo <- qq[1, seq_len(k + 1)]; hi <- qq[2, seq_len(k + 1)]; alo <- qq[1, k + 2]; ahi <- qq[2, k + 2]
    pe <- function(v) min(1, 2 * min((sum(v <= 0) + 1) / (reamostras + 1), (sum(v >= 0) + 1) / (reamostras + 1)))
    p <- vapply(seq_len(k + 1), function(j) pe(bt[, j]), 0); pa <- pe(bt[, k + 2])
  }
  apc <- function(b) 100 * (exp(b) - 1)
  alfa <- 1 - nivel
  tend <- function(lo_, hi_) ifelse(lo_ > 0, "crescente", ifelse(hi_ < 0, "decrescente", "estacionária"))
  tab <- data.frame(periodo = c(paste(limites[-length(limites)], limites[-1], sep = "-"), paste(t[1], t[n], sep = "-")),
                    medida = c(rep("APC", k + 1), "AAPC"), valor = apc(c(sl, a)), ic_inf = apc(c(lo, alo)), ic_sup = apc(c(hi, ahi)),
                    p = c(p, pa), stringsAsFactors = FALSE)
  tab$tendencia <- tend(tab$ic_inf, tab$ic_sup)
  sel <- data.frame(pontos = seq_along(bic) - 1, bic = bic, sqr = vapply(fits, `[[`, 0, "sse"),
                    localizacao = vapply(fits, function(x) if (length(x$tau)) paste(x$tau, collapse = ", ") else "-", ""), stringsAsFactors = FALSE)
  efeito <- data.frame(indice = paste0(tab$medida, " ", tab$periodo), valor = tab$valor, ic_inf = tab$ic_inf, ic_sup = tab$ic_sup,
                       principal = tab$medida == "AAPC", stringsAsFactors = FALSE)
  en <- identical(idioma, "en")
  fr <- paste(vapply(seq_len(nrow(tab)), function(i) paste0(tab$medida[i], " ", tab$periodo[i], " = ", .pr_num(tab$valor[i], 1, TRUE, decimal), "%, ",
                                                            .ee_ic(tab$ic_inf[i], tab$ic_sup[i], 1, TRUE, decimal, nivel, idioma)), ""), collapse = .ee_sep(decimal))
  apa <- if (en) paste0("Joinpoint regression (log-linear model, ", k, " joinpoint", if (k != 1) "s", " selected by BIC; ", if (ic == "parametrico") "parametric" else "empirical quantile",
                        " CIs) for ", if (is.null(taxa)) "the rate" else taxa, ": ", fr, ".") else
    paste0("A regressão por pontos de inflexão (modelo log-linear, ", k, " ponto", if (k != 1) "s", " escolhido", if (k != 1) "s", " pelo BIC; IC ",
           if (ic == "parametrico") "paramétrico" else "por quantis empíricos", ") para ", if (is.null(taxa)) "a taxa" else taxa, ": ", fr, ".")
  .ee_resultado(paste0("Tendência temporal por pontos de inflexão (", k, " ponto", if (k != 1) "s", ")"), n = n, tabela = tab, efeito = efeito,
                decisoes = c(paste0("Modelo log-linear contínuo, grade nos anos observados, mínimo de ", min_extremo, " observações da ponta e ", min_entre,
                                    " entre pontos; máximo de ", max_pontos, " ponto(s) pela regra do Joinpoint para ", n, " observações."),
                             paste0("Número de pontos pelo BIC do Joinpoint; variância ", if (pesos == "constante") "constante (homocedástica)" else "proporcional a 1/casos (Poisson)", "."),
                             if (ic == "parametrico") paste0("IC e p paramétricos (t com ", gl, " gl), localização dos pontos tratada como fixa.") else
                               paste0("IC por quantis empíricos de ", reamostras, " reamostras dos resíduos, reajustando a localização dos pontos (semente ", semente, ")."),
                             "Tendência: crescente ou decrescente quando o IC não inclui zero; senão, estacionária, como no relato das aulas.",
                             "O IC pode não bater com o do programa Joinpoint (o método do IC dele depende da sessão); declarar o método usado."),
                avisos = c(if (n < 7) "Menos de 7 observações: pontos de inflexão ficam instáveis; relatar a APC do período inteiro." else character(),
                           .ee_aviso_nao_sig(tab$p, alfa, o_que = "tendência", decimal = decimal, rotulos = paste(tab$medida, tab$periodo))),
                apa = apa, extra = list(tabelas_extra = list(`Seleção do número de pontos` = sel), modelo = f))
}

# taxa_padronizada(): padronização direta por idade (ou outra variável de estrato), com IC
# gama de Fay e Feuer (1997), o padrão dos órgãos de saúde e do epitools::ageadjust.direct.
# dados: uma linha por estrato, com casos, população e a população padrão do estrato;
# `grupo` (opcional) devolve uma taxa por grupo (por exemplo, por ano e sexo).
taxa_padronizada <- function(dados, casos, populacao, padrao, grupo = NULL, por = 1e5, nivel = .95, decimal = ".", idioma = "pt") {
  vars <- c(casos, populacao, padrao, grupo)
  falta <- setdiff(vars, names(dados)); if (length(falta)) stop("colunas ausentes: ", paste(falta, collapse = ", "))
  um <- function(d) {
    x <- d[[casos]]; nn <- d[[populacao]]; w <- d[[padrao]] / sum(d[[padrao]])
    tb <- sum(x) / sum(nn); tp <- sum(w * x / nn); v <- sum(w^2 * x / nn^2)
    wm <- max(w / nn); a <- (1 - nivel) / 2
    li <- stats::qgamma(a, tp^2 / v, scale = v / tp)
    ls <- stats::qgamma(1 - a, (tp + wm)^2 / (v + wm^2), scale = (v + wm^2) / (tp + wm))
    data.frame(casos = sum(x), taxa_bruta = tb * por, taxa_padronizada = tp * por, ic_inf = li * por, ic_sup = ls * por)
  }
  if (is.null(grupo)) tab <- um(dados) else {
    g <- interaction(dados[grupo], drop = TRUE, sep = " | ")
    tab <- do.call(rbind, lapply(levels(g), function(l) cbind(grupo = l, um(dados[g == l, ]))))
  }
  rownames(tab) <- NULL
  .ee_resultado(paste0("Taxa padronizada (método direto, por ", formatC(por, format = "d", big.mark = "."), ")"), n = nrow(dados), tabela = tab,
                decisoes = c("Padronização direta: média das taxas específicas ponderada pela distribuição da população padrão.",
                             "IC gama de Fay e Feuer (1997), conservador com poucos casos por estrato; declarar a população padrão usada (por exemplo, Censo 2010 ou OMS)."),
                apa = if (identical(idioma, "en")) "Directly standardized rates with Fay-Feuer gamma confidence intervals." else
                  "Taxas padronizadas pelo método direto, com IC gama de Fay e Feuer.",
                extra = list())
}
