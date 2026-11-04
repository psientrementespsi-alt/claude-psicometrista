# est_processos.R: moderação, mediação e processo condicional da skill psicometria-r.
# Modelos do tipo PROCESS (numeração de Hayes, 2022) estimados por MQO, com bootstrap
# percentílico próprio (reamostragem de casos, .lm.fit). Não usa a macro PROCESS.
#
# Funções públicas (todas devolvem um objeto "pr_resultado", ver est_comum.R):
#   moderacao(dados, y, x, w, ...)                       modelo 1: X × W, W contínua ou categórica
#   mediacao(dados, x, m, y, ...)                        modelo 4 (simples e paralela) e 6 (serial)
#   mediacao_moderada(dados, x, m, y, w, modelo, ...)    modelos 7 (W modera a), 14 (W modera b)
#                                                        e 15 (W modera b e c')
#
# Padrões do livro do usuário (Vol. 6): coeficientes brutos; IC bootstrap
# percentílico com 10.000 reamostras e semente 2026; efeitos condicionais na
# média e a ±1 DP do moderador; centralização só para leitura de b1 e b2
#; sem rótulos de mediação "completa" ou "parcial"; c não
# precisa ser significativo. Detalhes e divergências: references/processos.md.

# ---------------------------------------------------------------- auxiliares gerais
.ee_proc_tr <- function(idioma, pt, en) if (identical(idioma, "en")) en else pt

# Inteiro da frase (número de reamostras). Separador de milhar só quando não se confunde
# com o decimal da frase: vírgula em inglês com decimal = ".", ponto em português com
# decimal = ","; nos casos mistos, nenhum ("1000 resamples" ao lado de "ab = 0,241",
# nunca "1,000", que com decimal = "," se leria como um).
.ee_proc_int <- function(n, decimal = ".", idioma = "pt") {
  en <- identical(idioma, "en")
  marca <- if (en && decimal == ".") "," else if (!en && decimal == ",") "." else ""
  formatC(n, format = "d", big.mark = marca, decimal.mark = if (marca == ".") "," else ".")
}

# Grupo gravado como número (1/2/3, comum em banco vindo do SPSS ou de planilha) entra
# como contínuo, com um só coeficiente linear. Mesma regra de est_associacao.R
# (.ee_assoc_aviso_codigo): 3 a 7 valores distintos, todos inteiros, geram aviso; 0/1 é
# a codificação usual de X ou W binário e não avisa. `dica` diz o que fazer no papel.
.ee_proc_aviso_codigo <- function(v, nome, dica) {
  if (!is.numeric(v)) return(character())
  u <- sort(unique(v[!is.na(v)])); k <- length(u)
  if (k < 3 || k > 7 || any(abs(u - round(u)) > 1e-8)) return(character())
  paste0(nome, " é numérica com só ", k, " valores inteiros distintos (", paste(u, collapse = ", "),
         ") e entrou como contínua (um coeficiente linear). ", dica,
         " Códigos numéricos diretos para 3 ou mais categorias impõem ordem e distância que elas não têm (livro). ",
         "Se for escala (item Likert, contagem), o aviso pode ser ignorado.")
}
.ee_proc_dica_x <- "Se for grupo, esta função não cobre X com 3 ou mais categorias: rode com dois grupos por vez (X de dois níveis vira 0/1)."
.ee_proc_dica_cov <- "Se for grupo, transforme a covariável em fator antes (factor()): covariável fator entra como indicadores."
.ee_proc_aviso_covs <- function(d, covariaveis)
  unlist(lapply(covariaveis, function(v) .ee_proc_aviso_codigo(d[[v]], v, .ee_proc_dica_cov)))

# Nota de $decisoes para índices que saem sem IC (regra 6 da skill).
.ee_proc_nota_sem_ic <- function(pt, en, idioma)
  .ee_proc_tr(idioma, paste0(pt, " sem IC: não calculado (dizer isso no relato, ao lado do valor)."),
              paste0(en, " without CI: not computed (state this in the report, next to the value)."))

# Seleção, célula de texto vazia como ausente e deleção listwise. Célula "" ou só com
# espaços (o read.csv gera isso em coluna de texto) vira NA aqui mesmo, antes do
# .ee_dados, para que nunca vire categoria nem nível de referência.
.ee_proc_dados <- function(dados, vars, numericas = character()) {
  if (!is.data.frame(dados)) stop("dados precisa ser um data.frame")
  falta <- setdiff(vars, names(dados))
  if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "))
  d <- as.data.frame(dados)[, vars, drop = FALSE]
  n_vazios <- 0
  for (v in vars) {
    x <- d[[v]]
    if (!(is.character(x) || is.factor(x))) next
    xc <- as.character(x)
    vaz <- !is.na(xc) & !nzchar(trimws(xc))
    n_vazios <- n_vazios + sum(vaz)
    xc[vaz] <- NA
    d[[v]] <- if (is.factor(x)) droplevels(factor(xc, levels = levels(x)[nzchar(trimws(levels(x)))])) else xc
  }
  pre <- .ee_dados(d, vars, numericas)
  pre$n_vazios <- n_vazios
  pre$decisao <- if (n_vazios > 0) paste0(n_vazios, " célula(s) de texto vazia(s) tratada(s) como ausente(s) e excluída(s) por listwise.") else character()
  pre
}

# Casas da frase: o padrão 'casas', ampliado só quando um valor não nulo arredondaria
# para zero (passa a 2 algarismos significativos, até 10 casas), para não perder sinal e
# tamanho de coeficientes e índices pequenos.
.ee_proc_casas <- function(valores, casas) {
  v <- abs(valores[is.finite(valores) & valores != 0])
  v <- v[v < .5 * 10^-casas]
  if (!length(v)) return(casas)
  min(10, max(casas, floor(-log10(min(v))) + 2))
}

# IC com casas por limite (ver .ee_proc_casas): um limite perto de zero não vira "0.000".
.ee_proc_ic <- function(lo, hi, num, decimal, nivel, idioma) {
  rot <- if (identical(idioma, "en")) paste0(round(nivel * 100), "% CI") else paste0("IC ", round(nivel * 100), "%")
  paste0(rot, " [", num(lo), .ee_sep(decimal), num(hi), "]")
}

.ee_proc_aviso_esquema <- function(bt) if (identical(attr(bt, "esquema"), "sequencial"))
  "n × reamostras acima de 2e7: índices sorteados uma reamostra por vez para poupar memória; o IC continua reprodutível pela semente, mas não coincide com o do boot::boot sob a mesma semente." else character()

.ee_proc_nota_casas <- function(casas) paste0("Valor que arredondaria para zero com casas = ", casas,
  " saiu na frase com 2 algarismos significativos, para manter sinal e tamanho. Coeficiente muito pequeno costuma indicar variável em escala grande (ex.: renda em reais): reescalar (ex.: em milhares) deixa o relato mais legível.")

.ee_proc_checa_papeis <- function(...) {
  v <- c(...)
  if (anyDuplicated(v)) stop("a mesma variável aparece em mais de um papel: ",
                             paste(unique(v[duplicated(v)]), collapse = ", "))
  invisible(TRUE)
}

.ee_proc_checa_y <- function(v, nome) {
  if (!is.numeric(v)) stop("a variável ", nome, " precisa ser numérica")
  k <- length(unique(v))
  if (k < 2) stop("a variável ", nome, " não varia")
  if (k == 2) stop("a variável ", nome, " é dicotômica: modelos por MQO não se aplicam; ",
                   "use regressão logística (este módulo cobre desfechos contínuos)")
  invisible(TRUE)
}

# Codifica X (ou W de dois níveis): numérica fica como está; lógica, fator ou texto com
# dois níveis vira 0/1 (primeiro nível = 0). Mais de dois níveis: erro com orientação.
# Numérica com 3 a 7 valores inteiros (grupo gravado como 1/2/3) sai com aviso.
.ee_proc_x <- function(v, nome, dica = .ee_proc_dica_x) {
  if (is.logical(v)) v <- as.integer(v)
  if (is.numeric(v)) {
    k <- length(unique(v))
    if (k < 2) stop("a variável ", nome, " não varia")
    return(list(valor = as.numeric(v), binaria = k == 2, decisao = character(),
                aviso = .ee_proc_aviso_codigo(v, nome, dica)))
  }
  f <- droplevels(as.factor(v))
  if (nlevels(f) != 2) stop("a variável ", nome, " é categórica com ", nlevels(f),
                            " níveis; aqui ela precisa ser numérica ou ter dois níveis. ",
                            "Em moderação, um fator de vários níveis pode entrar como W.")
  list(valor = as.numeric(f) - 1, binaria = TRUE,
       decisao = paste0(nome, " codificada como 0 = '", levels(f)[1], "' e 1 = '", levels(f)[2], "'."),
       aviso = character())
}

# Covariáveis: numéricas entram como estão; fatores viram indicadores (primeiro nível =
# referência). Nas predições, as covariáveis ficam na média (como no PROCESS).
.ee_proc_cov <- function(d, covariaveis) {
  if (!length(covariaveis)) return(NULL)
  out <- list()
  for (v in covariaveis) {
    z <- d[[v]]
    if (is.logical(z)) z <- as.integer(z)
    if (is.numeric(z)) {
      if (length(unique(z)) < 2) stop("a covariável ", v, " não varia")
      out[[v]] <- as.numeric(z); next
    }
    f <- droplevels(as.factor(z))
    if (nlevels(f) < 2) stop("a covariável ", v, " não varia")
    for (lv in levels(f)[-1]) out[[paste0(v, "[", lv, "]")]] <- as.numeric(f == lv)
  }
  m <- do.call(cbind, out)
  colnames(m) <- names(out)
  m
}

# ---------------------------------------------------------------- MQO e testes
# Ajuste por MQO com EP clássico ou HC3 (Long & Ervin; o mesmo de sandwich::vcovHC
# type = "HC3"). X já traz a coluna do intercepto.
.ee_proc_mqo <- function(X, y, ep = "classico") {
  n <- nrow(X); k <- ncol(X)
  qx <- qr(X)
  if (qx$rank < k) stop("colinearidade perfeita entre os preditores: ",
                        paste(colnames(X)[qx$pivot[(qx$rank + 1):k]], collapse = ", "))
  gl <- n - k
  if (gl < 1) stop("casos insuficientes para o número de preditores")
  b <- qr.coef(qx, y)
  e <- as.vector(y - X %*% b)
  XtXi <- solve(crossprod(X))
  if (ep == "HC3") {
    h <- rowSums((X %*% XtXi) * X)
    if (any(h > 1 - 1e-10)) stop("há caso com alavancagem 1; o EP HC3 não é calculável")
    V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  } else {
    V <- XtXi * sum(e^2) / gl
  }
  names(b) <- colnames(X); dimnames(V) <- list(colnames(X), colnames(X))
  r2 <- 1 - sum(e^2) / sum((y - mean(y))^2)
  aj <- list(b = b, V = V, gl = gl, n = n, k = k, r2 = r2,
             r2_aj = 1 - (1 - r2) * (n - 1) / gl, ep = ep, e = e, X = X)
  aj$F <- if (k > 1) .ee_proc_wald(aj, 2:k) else c(F = NA, gl1 = 0, gl2 = gl, p = NA)
  aj
}

# Teste de Wald F de um bloco de coeficientes. Com EP clássico é idêntico ao F da
# mudança de R² (modelo com e sem o bloco); com HC3 é o F robusto que o PROCESS usa.
.ee_proc_wald <- function(aj, idx) {
  bq <- aj$b[idx]; Vq <- aj$V[idx, idx, drop = FALSE]
  Fv <- as.numeric(crossprod(bq, solve(Vq, bq))) / length(idx)
  c(F = Fv, gl1 = length(idx), gl2 = aj$gl, p = stats::pf(Fv, length(idx), aj$gl, lower.tail = FALSE))
}

.ee_proc_coefs <- function(aj, nivel, desfecho, caminhos = NULL) {
  ep <- sqrt(diag(aj$V)); t <- aj$b / ep
  tc <- stats::qt(1 - (1 - nivel) / 2, aj$gl)
  termo <- names(aj$b); termo[termo == "(Intercept)"] <- "constante"
  out <- data.frame(desfecho = desfecho, termo = termo, b = unname(aj$b), ep = unname(ep),
                    t = unname(t), gl = aj$gl, p = 2 * stats::pt(-abs(unname(t)), aj$gl),
                    ic_inf = unname(aj$b - tc * ep), ic_sup = unname(aj$b + tc * ep),
                    stringsAsFactors = FALSE)
  if (!is.null(caminhos)) {
    cam <- rep("", nrow(out)); cam[out$termo == "constante"] <- "constante"
    for (nm in names(caminhos)) cam[names(aj$b) == nm] <- caminhos[[nm]]
    out <- cbind(out[1:2], caminho = cam, out[-(1:2)])
  }
  out
}

# Combinação linear L'b com EP e IC pela distribuição t (efeitos condicionais).
.ee_proc_comb <- function(aj, L, nivel) {
  est <- sum(L * aj$b); ep <- sqrt(as.numeric(t(L) %*% aj$V %*% L))
  t <- est / ep; tc <- stats::qt(1 - (1 - nivel) / 2, aj$gl)
  c(efeito = est, ep = ep, t = t, gl = aj$gl, p = 2 * stats::pt(-abs(t), aj$gl),
    ic_inf = est - tc * ep, ic_sup = est + tc * ep)
}

# Breusch-Pagan estudentizado (Koenker): n·R² da regressão de e² nos preditores; igual a
# lmtest::bptest(modelo) com os padrões.
.ee_proc_bp <- function(aj) {
  e2 <- aj$e^2; X <- aj$X
  if (ncol(X) < 2) return(c(estatistica = NA, gl = 0, p = NA))
  f <- stats::lm.fit(X, e2)
  r2 <- 1 - sum(f$residuals^2) / sum((e2 - mean(e2))^2)
  bp <- aj$n * r2
  c(estatistica = bp, gl = ncol(X) - 1, p = stats::pchisq(bp, ncol(X) - 1, lower.tail = FALSE))
}

.ee_proc_pressupostos <- function(ajustes, ep, idioma = "pt") {
  linhas <- list()
  for (nm in names(ajustes)) {
    aj <- ajustes[[nm]]
    bp <- .ee_proc_bp(aj)
    leit <- if (is.na(bp[["p"]])) "não calculável" else if (bp[["p"]] < .05) {
      if (ep == "HC3") "variância dos resíduos não constante; o EP HC3 já está em uso" else
        "variância dos resíduos não constante: considerar ep = \"HC3\" para os testes t (o IC bootstrap não depende disso)"
    } else "sem indício forte de heterocedasticidade"
    linhas[[length(linhas) + 1]] <- data.frame(checagem = paste0("Breusch-Pagan (", nm, ")"),
                                               estatistica = bp[["estatistica"]], gl = bp[["gl"]],
                                               p = bp[["p"]], leitura = leit, stringsAsFactors = FALSE)
    n <- aj$n
    if (n >= 3 && n <= 5000 && stats::sd(aj$e) > 0) {
      sw <- stats::shapiro.test(aj$e)
      leit <- if (sw$p.value >= .05) "resíduos compatíveis com normalidade" else if (n >= 300)
        "desvio da normalidade nos resíduos; com n grande o teste acusa desvios irrelevantes, olhar o Q-Q plot" else
          "desvio da normalidade nos resíduos: conferir Q-Q plot e casos atípicos; o IC bootstrap não supõe normalidade"
      linhas[[length(linhas) + 1]] <- data.frame(checagem = paste0("Shapiro-Wilk dos resíduos (", nm, ")"),
                                                 estatistica = unname(sw$statistic), gl = NA, p = sw$p.value,
                                                 leitura = leit, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, linhas)
}

# ---------------------------------------------------------------- bootstrap
# Reamostragem de casos com reposição. Com n × B até 2e7 os índices são sorteados de
# uma vez, na mesma ordem do boot::boot (sample.int(n, n * B) arrumado em B × n), de
# modo que o resultado coincide com boot::boot sob a mesma semente; acima disso, uma
# reamostra por vez (economiza memória). A semente do usuário é restaurada no fim.
.ee_proc_bootstrap <- function(n, B, semente, f, nomes) {
  if (!is.numeric(B) || B < 100) stop("reamostras precisa ser um número >= 100")
  B <- as.integer(B)
  if (!is.null(semente)) {
    ge <- globalenv()
    tinha <- exists(".Random.seed", envir = ge, inherits = FALSE)
    if (tinha) antiga <- get(".Random.seed", envir = ge)
    on.exit({
      if (tinha) assign(".Random.seed", antiga, envir = ge) else
        if (exists(".Random.seed", envir = ge, inherits = FALSE)) rm(".Random.seed", envir = ge)
    }, add = TRUE)
    set.seed(semente)
  }
  em_bloco <- as.double(n) * B <= 2e7
  if (em_bloco) {
    ind <- sample.int(n, n * B, replace = TRUE)
    dim(ind) <- c(B, n)
    ind <- t(ind)
  }
  res <- matrix(NA_real_, B, length(nomes), dimnames = list(NULL, nomes))
  for (r in seq_len(B)) {
    i <- if (em_bloco) ind[, r] else sample.int(n, n, replace = TRUE)
    v <- f(i)
    if (!is.null(v)) res[r, ] <- v
  }
  attr(res, "esquema") <- if (em_bloco) "bloco" else "sequencial"
  res
}

# Coeficientes de uma reamostra; NULL se a reamostra ficou sem posto completo.
.ee_proc_coef_boot <- function(X, y, i) {
  f <- .lm.fit(X[i, , drop = FALSE], y[i])
  if (f$rank < ncol(X)) return(NULL)
  f$coefficients
}

# IC percentílico com a mesma interpolação de boot.ci(type = "perc") (norm.inter).
.ee_proc_percentil <- function(t, nivel) {
  t <- sort(t[is.finite(t)]); R <- length(t)
  if (R < 2) return(c(NA_real_, NA_real_))
  alfa <- c((1 - nivel) / 2, 1 - (1 - nivel) / 2)
  rk <- (R + 1) * alfa; k <- trunc(rk)
  vapply(seq_along(alfa), function(j) {
    if (k[j] == rk[j] && k[j] >= 1 && k[j] <= R) return(t[k[j]])
    if (k[j] <= 0) return(t[1])
    if (k[j] >= R) return(t[R])
    z <- stats::qnorm(alfa[j]); z1 <- stats::qnorm(k[j] / (R + 1)); z2 <- stats::qnorm((k[j] + 1) / (R + 1))
    t[k[j]] + (z - z1) / (z2 - z1) * (t[k[j] + 1] - t[k[j]])
  }, numeric(1))
}

.ee_proc_resumo_boot <- function(est, bt, nivel) {
  ic <- t(apply(bt, 2, .ee_proc_percentil, nivel = nivel))
  data.frame(indice = names(est), valor = unname(est), ep_boot = apply(bt, 2, stats::sd, na.rm = TRUE),
             ic_inf = ic[, 1], ic_sup = ic[, 2], row.names = NULL, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------- sondagem da interação
# Valores do moderador para sondar: "media_dp" (livro), "percentis" 16/50/84
# (PROCESS; quantile tipo 7 do R) ou vetor numérico de valores teóricos.
.ee_proc_valores_w <- function(wv, valores_w, nome, idioma = "pt") {
  avisos <- character()
  if (length(unique(wv)) == 2) {
    u <- sort(unique(wv))
    return(list(valores = u, rotulos = format(u), metodo = "binaria", avisos = avisos))
  }
  if (is.numeric(valores_w)) {
    v <- as.numeric(valores_w); rot <- format(v)
    metodo <- "valores informados"
  } else {
    metodo <- match.arg(valores_w[1], c("media_dp", "percentis"))
    if (metodo == "media_dp") {
      mu <- mean(wv); s <- stats::sd(wv)
      v <- c(mu - s, mu, mu + s)
      rot <- .ee_proc_tr(idioma, c("-1 DP", "média", "+1 DP"), c("-1 SD", "mean", "+1 SD"))
    } else {
      v <- unname(stats::quantile(wv, c(.16, .5, .84), type = 7)); rot <- c("P16", "P50", "P84")
    }
  }
  fora <- v < min(wv) | v > max(wv)
  if (any(fora)) avisos <- c(avisos, paste0("Valor de ", nome, " fora da faixa observada (",
                                            paste(round(v[fora], 3), collapse = ", "),
                                            "): o efeito condicional ali é extrapolação",
                                            if (metodo == "media_dp") "; com W assimétrica, os percentis 16/50/84 evitam isso (valores_w = \"percentis\")" else "", "."))
  list(valores = v, rotulos = rot, metodo = metodo, avisos = avisos)
}

# Johnson-Neyman: fronteiras analíticas em que o efeito condicional θ(w) = b1 + b3·w tem
# |t| igual ao t crítico, isto é, raízes de (b3² - t²·V33)·w² + 2(b1·b3 - t²·V13)·w +
# (b1² - t²·V11) = 0, mantidas as que caem na faixa observada de W. Tudo na escala do
# modelo; 'centro' devolve à escala original.
.ee_proc_jn <- function(b1, b3, v11, v13, v33, gl, nivel, w_modelo, centro = 0, prop_min = .05) {
  tc <- stats::qt(1 - (1 - nivel) / 2, gl)
  A <- b3^2 - tc^2 * v33; Bq <- 2 * (b1 * b3 - tc^2 * v13); C <- b1^2 - tc^2 * v11
  raizes <- if (abs(A) < 1e-14) {
    if (abs(Bq) < 1e-14) numeric() else -C / Bq
  } else {
    disc <- Bq^2 - 4 * A * C
    if (disc < 0) numeric() else sort((-Bq + c(-1, 1) * sqrt(disc)) / (2 * A))
  }
  lo <- min(w_modelo); hi <- max(w_modelo)
  raizes <- unique(raizes[raizes > lo & raizes < hi])
  cond <- function(w) {
    est <- b1 + b3 * w; ep <- sqrt(v11 + 2 * w * v13 + w^2 * v33); t <- est / ep
    data.frame(w = w + centro, efeito = est, ep = ep, t = t, gl = gl, p = 2 * stats::pt(-abs(t), gl),
               ic_inf = est - tc * ep, ic_sup = est + tc * ep)
  }
  pontos <- if (length(raizes)) cbind(cond(raizes),
                                       pct_abaixo = vapply(raizes, function(r) 100 * mean(w_modelo < r), numeric(1)),
                                       pct_acima = vapply(raizes, function(r) 100 * mean(w_modelo > r), numeric(1))) else NULL
  cortes <- c(lo, raizes, hi)
  regioes <- do.call(rbind, lapply(seq_len(length(cortes) - 1), function(j) {
    meio <- (cortes[j] + cortes[j + 1]) / 2; cm <- cond(meio)
    dentro <- if (j == 1) w_modelo <= cortes[j + 1] & w_modelo >= cortes[j] else w_modelo > cortes[j] & w_modelo <= cortes[j + 1]
    data.frame(de = cortes[j] + centro, ate = cortes[j + 1] + centro, significativo = cm$p < (1 - nivel),
               sinal = if (cm$efeito >= 0) "positivo" else "negativo", n_casos = sum(dentro),
               pct_casos = 100 * mean(dentro))
  }))
  grade <- cond(sort(unique(c(seq(lo, hi, length.out = 20), raizes))))
  list(pontos = pontos, regioes = regioes, grade = grade, t_critico = tc,
       pequenas = regioes$significativo & regioes$pct_casos < 100 * prop_min,
       faixa = c(lo, hi) + centro)
}

.ee_proc_jn_txt <- function(jn, nome_x, nome_w, casas, decimal, idioma) {
  sig <- jn$regioes[jn$regioes$significativo, , drop = FALSE]
  if (!nrow(sig)) return(.ee_proc_tr(idioma,
    paste0("Pela técnica de Johnson-Neyman, o efeito de ", nome_x, " não foi significativo em nenhum ponto da faixa observada de ", nome_w, "."),
    paste0("The Johnson-Neyman technique showed no value of ", nome_w, " within the observed range at which the effect of ", nome_x, " was significant.")))
  if (nrow(sig) == 1 && nrow(jn$regioes) == 1) return(.ee_proc_tr(idioma,
    paste0("Pela técnica de Johnson-Neyman, o efeito de ", nome_x, " foi significativo em toda a faixa observada de ", nome_w, "."),
    paste0("The Johnson-Neyman technique showed that the effect of ", nome_x, " was significant across the observed range of ", nome_w, ".")))
  pedacos <- vapply(seq_len(nrow(sig)), function(j) {
    r <- sig[j, ]
    ini <- isTRUE(all.equal(r$de, jn$faixa[1])); fim <- isTRUE(all.equal(r$ate, jn$faixa[2]))
    faixa <- if (ini) paste0(nome_w, " ≤ ", .pr_num(r$ate, casas, TRUE, decimal)) else
      if (fim) paste0(nome_w, " ≥ ", .pr_num(r$de, casas, TRUE, decimal)) else
        .ee_proc_tr(idioma, paste0(nome_w, " entre ", .pr_num(r$de, casas, TRUE, decimal), " e ", .pr_num(r$ate, casas, TRUE, decimal)),
                    paste0(nome_w, " between ", .pr_num(r$de, casas, TRUE, decimal), " and ", .pr_num(r$ate, casas, TRUE, decimal)))
    sinal <- .ee_proc_tr(idioma, r$sinal, if (r$sinal == "positivo") "positive" else "negative")
    paste0(faixa, " (", sinal, "; ", .pr_num(r$pct_casos, 1, TRUE, decimal), "% ",
           .ee_proc_tr(idioma, "dos casos", "of cases"), ")")
  }, character(1))
  .ee_proc_tr(idioma,
    paste0("Pela técnica de Johnson-Neyman, o efeito de ", nome_x, " foi significativo quando ", paste(pedacos, collapse = .ee_proc_tr(idioma, " ou ", " or ")), "."),
    paste0("The Johnson-Neyman technique showed that the effect of ", nome_x, " was significant when ", paste(pedacos, collapse = " or "), "."))
}

.ee_proc_padrao <- function(ef) {
  e1 <- ef[1]; e2 <- ef[length(ef)]
  if (sign(e1) != sign(e2) && e1 != 0 && e2 != 0) return("inverte")
  if (abs(e2) > abs(e1)) "amplifica" else "atenua"
}

# ---------------------------------------------------------------- moderação (modelo 1)
moderacao <- function(dados, y, x, w, covariaveis = NULL, centralizar = FALSE,
                      valores_w = c("media_dp", "percentis"), johnson_neyman = TRUE,
                      ep = c("classico", "HC3"), nivel = .95, jn_prop_min = .05,
                      casas = 3, decimal = ".", idioma = "pt") {
  ep <- match.arg(ep)
  .ee_proc_checa_papeis(y, x, w, covariaveis)
  pre <- .ee_proc_dados(dados, c(y, x, w, covariaveis), numericas = y)
  d <- pre$d; n <- nrow(d)
  .ee_proc_checa_y(d[[y]], y)
  yv <- d[[y]]
  xx <- .ee_proc_x(d[[x]], x, "Se for grupo, X com 3 ou mais categorias não é coberto aqui: troque os papéis (o grupo como W, em fator) ou rode com dois grupos por vez.")
  decisoes <- c(paste0("Modelo: ", y, " = b0 + b1·", x, " + b2·", w, " + b3·", x, "×", w,
                       if (length(covariaveis)) " + covariáveis" else "",
                       ", por MQO; X, W e o produto entram juntos. Equivale ao modelo 1 do PROCESS."),
                pre$decisao, xx$decisao)
  avisos <- c(xx$aviso, .ee_proc_aviso_covs(d, covariaveis))

  wbruta <- d[[w]]
  w_cat <- !(is.numeric(wbruta) || is.logical(wbruta))
  if (w_cat) {
    wf <- droplevels(as.factor(wbruta))
    if (nlevels(wf) < 2) stop("o moderador ", w, " precisa de pelo menos 2 níveis com dados")
    niveis <- levels(wf)
    Wd <- vapply(niveis[-1], function(lv) as.numeric(wf == lv), numeric(n))
    if (!is.matrix(Wd)) Wd <- matrix(Wd, ncol = 1)
    colnames(Wd) <- paste0(w, "[", niveis[-1], "]")
    decisoes <- c(decisoes, paste0(w, " categórica com ", length(niveis), " níveis; codificação por indicadores com '",
                                   niveis[1], "' como referência. O teste da moderação é o F conjunto dos ",
                                   ncol(Wd), " produtos."))
    w_cont <- FALSE; wnum <- NULL
  } else {
    wnum <- as.numeric(wbruta)
    if (length(unique(wnum)) < 2) stop("o moderador ", w, " não varia")
    w_cont <- length(unique(wnum)) > 2
    avisos <- c(avisos, .ee_proc_aviso_codigo(wnum, w, paste0(
      "Se for grupo, transforme em fator antes (dados$", w, " <- factor(dados$", w,
      ")): a moderacao() aceita W categórica, com indicadores e F conjunto dos produtos; sondar W em média ± 1 DP de códigos não tem leitura.")))
  }
  C <- .ee_proc_cov(d, covariaveis)

  monta <- function(cx, cw) {
    xc <- xx$valor - cx
    Wm <- if (w_cat) Wd else matrix(wnum - cw, ncol = 1, dimnames = list(NULL, w))
    P <- xc * Wm; colnames(P) <- paste0(x, " × ", colnames(Wm))
    X <- cbind(1, xc, Wm, P, C)
    colnames(X) <- c("(Intercept)", x, colnames(Wm), colnames(P), colnames(C))
    list(X = X, iw = 2 + seq_len(ncol(Wm)), ip = 2 + ncol(Wm) + seq_len(ncol(Wm)))
  }
  pode_x <- !xx$binaria; pode_w <- !w_cat && w_cont
  cx <- if (centralizar && pode_x) mean(xx$valor) else 0
  cw <- if (centralizar && pode_w) mean(wnum) else 0
  des <- monta(cx, cw)
  aj <- .ee_proc_mqo(des$X, yv, ep)
  aj0 <- .ee_proc_mqo(des$X[, -des$ip, drop = FALSE], yv, ep)
  coefs <- .ee_proc_coefs(aj, nivel, y)

  # teste da interação: ΔR² e F (com HC3, F de Wald robusto, como no PROCESS)
  dr2 <- aj$r2 - aj0$r2
  wt <- .ee_proc_wald(aj, des$ip)
  f2 <- dr2 / (1 - aj$r2)
  teste <- data.frame(termo = paste0(x, " × ", w), delta_r2 = dr2, F = wt[["F"]], gl1 = wt[["gl1"]],
                      gl2 = wt[["gl2"]], p = wt[["p"]], f2 = f2)
  sig_int <- wt[["p"]] < 1 - nivel
  decisoes <- c(decisoes, paste0("Teste da moderação: ΔR² = ", round(dr2, 4), ", F(", wt[["gl1"]], ", ", wt[["gl2"]], ") = ",
                                 round(wt[["F"]], 3), ", p = ", signif(wt[["p"]], 3), "; f² = ΔR²/(1 - R²) = ", round(f2, 4),
                                 ".", if (ep == "HC3") " Com HC3, o F é o de Wald robusto; o ΔR² é descritivo." else
                                   if (length(des$ip) == 1) " Com EP clássico, esse F é idêntico ao t² do produto." else " Com EP clássico, esse F é o da mudança de R²."),
                .ee_proc_nota_sem_ic("ΔR², f² e R² do modelo", "ΔR², f² and model R²", idioma))

  # alternativa de centralização (só leitura de b1 e b2; b3, EP(b3) e R² não mudam)
  alt <- NULL
  if (pode_x || pode_w) {
    cx2 <- if (!centralizar && pode_x) mean(xx$valor) else 0
    cw2 <- if (!centralizar && pode_w) mean(wnum) else 0
    des2 <- monta(cx2, cw2)
    aj2 <- .ee_proc_mqo(des2$X, yv, ep)
    alt <- .ee_proc_coefs(aj2, nivel, y)
    igual <- isTRUE(all.equal(unname(aj$b[des$ip]), unname(aj2$b[des2$ip]), tolerance = 1e-8)) &&
      isTRUE(all.equal(unname(sqrt(diag(aj$V))[des$ip]), unname(sqrt(diag(aj2$V))[des2$ip]), tolerance = 1e-8)) &&
      isTRUE(all.equal(aj$r2, aj2$r2, tolerance = 1e-10))
    cen_txt <- paste(c(if (pode_x) x, if (pode_w) w), collapse = " e ")
    b1a <- aj2$b[[x]]
    decisoes <- c(decisoes, paste0(
      if (centralizar) paste0("Centralizadas pela média: ", cen_txt, " (variáveis binárias e categóricas não são centralizadas). ") else
        "Variáveis na escala original (centralizar = FALSE): b1 é o efeito de X quando W = 0. ",
      "No modelo ", if (centralizar) "sem centralizar" else paste0("com ", cen_txt, " centralizada(s)"), ", b1 = ", round(b1a, 4),
      "; b3, EP(b3) e R² ", if (igual) "são idênticos" else "DIFEREM (conferir)", " nos dois modelos: centralizar muda só a leitura de b1 e b2. Tabela em $coeficientes_alternativos."))
  } else {
    decisoes <- c(decisoes, "X e W binárias ou categóricas: nada a centralizar.")
  }

  # efeitos condicionais de X
  ix <- 2
  if (w_cat) {
    L0 <- numeric(ncol(des$X)); L0[ix] <- 1
    ef <- lapply(seq_along(niveis), function(j) {
      L <- L0; if (j > 1) L[des$ip[j - 1]] <- 1
      .ee_proc_comb(aj, L, nivel)
    })
    cond <- data.frame(w_nivel = niveis, w = NA_real_, do.call(rbind, ef), stringsAsFactors = FALSE)
    vw <- list(metodo = "niveis", avisos = character())
  } else {
    vw <- .ee_proc_valores_w(wnum, valores_w, w, idioma)
    avisos <- c(avisos, vw$avisos)
    ef <- lapply(vw$valores, function(v) { L <- numeric(ncol(des$X)); L[ix] <- 1; L[des$ip] <- v - cw; .ee_proc_comb(aj, L, nivel) })
    cond <- data.frame(w_nivel = vw$rotulos, w = vw$valores, do.call(rbind, ef), stringsAsFactors = FALSE)
  }
  decisoes <- c(decisoes, switch(vw$metodo,
    media_dp = "Efeitos condicionais na média e a ±1 DP de W (padrão do livro); percentis 16/50/84 do PROCESS com valores_w = \"percentis\".",
    percentis = "Efeitos condicionais nos percentis 16, 50 e 84 de W (padrão do PROCESS; quantile tipo 7 do R).",
    binaria = "W com dois valores: efeito condicional de X em cada um.",
    niveis = "Efeito condicional de X em cada nível de W (b1 no nível de referência, b1 + b3j nos demais).",
    "Efeitos condicionais nos valores de W informados (valores teóricos ou clínicos)."))
  padrao <- if (nrow(cond) >= 2 && !w_cat) .ee_proc_padrao(cond$efeito) else NA
  if (!is.na(padrao)) decisoes <- c(decisoes, paste0("Padrão descritivo nos valores sondados: W ",
    switch(padrao, inverte = "inverte o sinal do efeito de X", amplifica = "amplifica o efeito de X (módulo maior com W alto)",
           atenua = "atenua o efeito de X (módulo menor com W alto)"), ". É leitura dos pontos, não teste."))
  if (!sig_int) decisoes <- c(decisoes, "Interação não significativa: os efeitos condicionais ficam no objeto só como descrição, não como achado.")

  # Johnson-Neyman
  jn <- NULL
  if (johnson_neyman && w_cont) {
    jn <- .ee_proc_jn(aj$b[[ix]], aj$b[[des$ip]], aj$V[ix, ix], aj$V[ix, des$ip], aj$V[des$ip, des$ip],
                      aj$gl, nivel, wnum - cw, cw, jn_prop_min)
    decisoes <- c(decisoes, paste0("Johnson-Neyman por solução analítica (θ(W) = b1 + b3·W com |t| = t crítico = ",
                                   round(jn$t_critico, 4), "): ",
                                   if (is.null(jn$pontos)) "nenhuma fronteira dentro da faixa observada de W." else
                                     paste0("fronteira(s) em W = ", paste(round(jn$pontos$w, 4), collapse = " e "), " (",
                                            paste0(round(jn$pontos$pct_abaixo, 1), "% dos casos abaixo", collapse = "; "), ")."),
                                   " Grade de 20 pontos em $johnson_neyman$grade."))
    if (any(jn$pequenas)) avisos <- c(avisos, paste0("Região de significância de Johnson-Neyman com menos de ",
                                                     round(100 * jn_prop_min), "% dos casos: não interpretar sem olhar quantos casos há ali."))
  } else if (johnson_neyman && !w_cont) {
    decisoes <- c(decisoes, "Johnson-Neyman não se aplica: W binária ou categórica.")
  }

  # grade para o gráfico de retas simples (covariáveis na média)
  xs <- if (xx$binaria) sort(unique(xx$valor)) else .ee_proc_valores_w(xx$valor, if (is.numeric(valores_w)) "media_dp" else valores_w, x, idioma)$valores
  wlev <- if (w_cat) niveis else cond$w
  grade <- expand.grid(x = xs, w = wlev, stringsAsFactors = FALSE)
  cm <- if (is.null(C)) numeric() else colMeans(C)
  grade$y_previsto <- vapply(seq_len(nrow(grade)), function(r) {
    xc <- grade$x[r] - cx
    wm <- if (w_cat) as.numeric(niveis[-1] == grade$w[r]) else grade$w[r] - cw
    sum(c(1, xc, wm, xc * wm, cm) * aj$b)
  }, numeric(1))
  names(grade)[1:2] <- c(x, w)

  if (ep == "HC3") decisoes <- c(decisoes, "EP robusto HC3 em coeficientes, efeitos condicionais e Johnson-Neyman (opção do curso e do PROCESS; o livro usa o EP clássico).")
  pres <- .ee_proc_pressupostos(list(modelo = aj), ep, idioma)
  if (any(grepl("^Breusch", pres$checagem) & pres$p < .05) && ep == "classico")
    avisos <- c(avisos, "Breusch-Pagan sugere heterocedasticidade: considerar ep = \"HC3\".")

  efeito <- data.frame(indice = c("ΔR² da interação", "f² da interação", "R² do modelo"),
                       valor = c(dr2, f2, aj$r2), ic_inf = NA_real_, ic_sup = NA_real_,
                       principal = c(TRUE, FALSE, FALSE), stringsAsFactors = FALSE)

  # frase APA
  sep <- .ee_sep(decimal); eps <- .ee_proc_tr(idioma, "EP", "SE")
  ampliou <- FALSE
  num <- function(v, grupo = v) { cs <- .ee_proc_casas(grupo, casas); if (cs > casas) ampliou <<- TRUE; .pr_num(v, cs, TRUE, decimal) }
  ip1 <- des$ip
  int_txt <- if (length(ip1) == 1) {
    g1 <- c(aj$b[[ip1]], sqrt(aj$V[ip1, ip1]))
    .ee_junta(paste0("b = ", num(g1[1], g1)), paste0(eps, " = ", num(g1[2], g1)),
              .ee_estat("t", aj$b[[ip1]] / sqrt(aj$V[ip1, ip1]), aj$gl, decimal), .pr_p(wt[["p"]], decimal),
              .ee_efeito_txt("ΔR²", dr2, casas = 3, decimal = decimal, idioma = idioma, nivel = nivel),
              .ee_efeito_txt("f²", f2, limitado = FALSE, casas = 3, decimal = decimal, idioma = idioma, nivel = nivel), decimal = decimal)
  } else {
    .ee_junta(.ee_estat("F", wt[["F"]], c(wt[["gl1"]], wt[["gl2"]]), decimal), .pr_p(wt[["p"]], decimal),
              .ee_efeito_txt("ΔR²", dr2, casas = 3, decimal = decimal, idioma = idioma, nivel = nivel),
              .ee_efeito_txt("f²", f2, limitado = FALSE, casas = 3, decimal = decimal, idioma = idioma, nivel = nivel), decimal = decimal)
  }
  frase1 <- .ee_proc_tr(idioma,
    paste0("A interação entre ", x, " e ", w, " na predição de ", y, if (sig_int) " foi significativa, " else " não foi significativa, ",
           int_txt, if (ep == "HC3") " (EP robusto HC3)" else "", "."),
    paste0("The ", x, " × ", w, " interaction predicting ", y, if (sig_int) " was significant, " else " was not significant, ",
           int_txt, if (ep == "HC3") " (HC3 robust SE)" else "", "."))
  apa <- frase1
  if (sig_int) {
    itens <- vapply(seq_len(nrow(cond)), function(j) {
      r <- cond[j, ]
      onde <- if (w_cat) .ee_proc_tr(idioma, paste0("em ", w, " = ", r$w_nivel), paste0("at ", w, " = ", r$w_nivel)) else
        paste0(.ee_proc_tr(idioma, "em ", "at "), w, " = ", .pr_num(r$w, casas, TRUE, decimal),
               if (vw$metodo %in% c("media_dp", "percentis")) paste0(" (", r$w_nivel, ")") else "")
      paste0("b = ", num(r$efeito), " (", .pr_p(r$p, decimal), ") ", onde)
    }, character(1))
    apa <- c(apa, .ee_proc_tr(idioma,
      paste0("O efeito condicional de ", x, " sobre ", y, " foi ", paste(itens, collapse = "; "), "."),
      paste0("The conditional effect of ", x, " on ", y, " was ", paste(itens, collapse = "; "), ".")))
    if (!is.null(jn)) apa <- c(apa, .ee_proc_jn_txt(jn, x, w, casas, decimal, idioma))
  }
  if (ampliou) decisoes <- c(decisoes, .ee_proc_nota_casas(casas))

  tab_mod <- data.frame(r2 = aj$r2, r2_aj = aj$r2_aj, F = aj$F[["F"]], gl1 = aj$F[["gl1"]], gl2 = aj$F[["gl2"]], p = aj$F[["p"]])
  .ee_resultado(paste0("Moderação (modelo 1): ", x, " × ", w, " → ", y), n = n, n_excluidos = pre$n_excluidos,
                tabela = coefs, efeito = efeito, pressupostos = pres, posthoc = cond,
                decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(teste_interacao = teste, resumo_modelo = tab_mod,
                             johnson_neyman = jn, grade_grafico = grade,
                             coeficientes_alternativos = alt, centros = c(x = cx, w = cw),
                             padrao = padrao, modelo = list(b = aj$b, vcov = aj$V, gl = aj$gl, ep = ep)))
}

# ---------------------------------------------------------------- mediação (modelos 4 e 6)
.ee_proc_frase_delineamento <- function(delineamento, idioma = "pt") {
  switch(delineamento,
    transversal = .ee_proc_tr(idioma,
      "Delineamento transversal: os caminhos são associações condicionais, e a ordem X → M → Y vem da teoria, não dos dados. Estimativas transversais de mediação podem divergir muito das longitudinais (Maxwell & Cole, 2007). Redigir em termos associativos, por exemplo 'a associação entre X e Y foi explicada, em termos estatísticos, por M'.",
      "Cross-sectional design: paths are conditional associations and the X → M → Y order comes from theory, not from the data; cross-sectional mediation estimates can be badly biased relative to longitudinal ones (Maxwell & Cole, 2007)."),
    longitudinal = .ee_proc_tr(idioma,
      "Delineamento longitudinal: a precedência temporal fica mais forte, mas confundidores não medidos continuam possíveis.",
      "Longitudinal design: temporal precedence is stronger, but unmeasured confounders remain possible."),
    experimental = .ee_proc_tr(idioma,
      "X manipulado: os caminhos a e c têm apoio causal; M não foi manipulado, então o caminho b continua correlacional e sujeito a confundimento.",
      "Manipulated X supports the a and c paths causally; M was not manipulated, so the b path remains correlational."))
}

mediacao <- function(dados, x, m, y, covariaveis = NULL, tipo = c("paralela", "serial"),
                     reamostras = 10000, semente = 2026, nivel = .95, padronizado = FALSE,
                     contrastes = TRUE, ep = c("classico", "HC3"), proporcao_mediada = FALSE,
                     delineamento = c("transversal", "longitudinal", "experimental"),
                     casas = 3, decimal = ".", idioma = "pt") {
  tipo <- match.arg(tipo); ep <- match.arg(ep); delineamento <- match.arg(delineamento)
  k <- length(m)
  if (k < 1) stop("informe pelo menos um mediador em m")
  if (k > 6) stop("no máximo 6 mediadores (o mesmo limite do PROCESS)")
  .ee_proc_checa_papeis(x, m, y, covariaveis)
  pre <- .ee_proc_dados(dados, c(x, m, y, covariaveis), numericas = c(m, y))
  d <- pre$d; n <- nrow(d)
  .ee_proc_checa_y(d[[y]], y)
  for (mm in m) if (length(unique(d[[mm]])) < 3) stop("o mediador ", mm, " precisa ser contínuo (mais de dois valores)")
  xx <- .ee_proc_x(d[[x]], x)
  xv <- xx$valor; yv <- d[[y]]; M <- as.matrix(d[m])
  C <- .ee_proc_cov(d, covariaveis)
  serial <- tipo == "serial" && k > 1
  modelo_nome <- if (k == 1) "simples (modelo 4)" else if (serial) "serial (modelo 6)" else "paralela (modelo 4)"

  # equações
  Xm <- lapply(seq_len(k), function(j) {
    ant <- if (serial && j > 1) M[, seq_len(j - 1), drop = FALSE] else NULL
    X <- cbind(1, xv, ant, C); colnames(X) <- c("(Intercept)", x, colnames(ant), colnames(C)); X
  })
  Xy <- cbind(1, xv, M, C); colnames(Xy) <- c("(Intercept)", x, m, colnames(C))
  Xt <- cbind(1, xv, C); colnames(Xt) <- c("(Intercept)", x, colnames(C))
  ajm <- lapply(seq_len(k), function(j) .ee_proc_mqo(Xm[[j]], M[, j], ep))
  ajy <- .ee_proc_mqo(Xy, yv, ep)
  ajt <- .ee_proc_mqo(Xt, yv, ep)

  # rótulos dos caminhos
  sufixo <- function(j) if (k == 1) "" else as.character(j)
  tabs <- lapply(seq_len(k), function(j) {
    cam <- list(); cam[[x]] <- paste0("a", sufixo(j))
    if (serial && j > 1) for (l in seq_len(j - 1)) cam[[m[l]]] <- paste0("d", j, l)
    if (!is.null(C)) for (cc in colnames(C)) cam[[cc]] <- "covariável"
    .ee_proc_coefs(ajm[[j]], nivel, m[j], cam)
  })
  camy <- list(); camy[[x]] <- "c'"
  for (j in seq_len(k)) camy[[m[j]]] <- paste0("b", sufixo(j))
  if (!is.null(C)) for (cc in colnames(C)) camy[[cc]] <- "covariável"
  camt <- list(); camt[[x]] <- "c"
  if (!is.null(C)) for (cc in colnames(C)) camt[[cc]] <- "covariável"
  tabela <- rbind(do.call(rbind, tabs), .ee_proc_coefs(ajy, nivel, y, camy),
                  transform(.ee_proc_coefs(ajt, nivel, y, camt), desfecho = paste0(y, " (efeito total)")))

  # efeitos indiretos: caminhos específicos
  caminhos <- if (serial) {
    unlist(lapply(seq_len(k), function(s) utils::combn(k, s, simplify = FALSE)), recursive = FALSE)
  } else as.list(seq_len(k))
  rot_cam <- vapply(caminhos, function(s) paste(c(x, m[s], y), collapse = " → "), character(1))
  ind_de <- function(bm, by) {
    vapply(caminhos, function(s) {
      v <- bm[[s[1]]][[2]]
      if (length(s) > 1) for (q in 2:length(s)) v <- v * bm[[s[q]]][[2 + s[q - 1]]]
      v * by[[2 + s[length(s)]]]
    }, numeric(1))
  }
  n_esp <- length(caminhos)
  pares <- if (contrastes && n_esp > 1) utils::combn(n_esp, 2, simplify = FALSE) else list()
  rot_con <- vapply(pares, function(p) paste0("contraste: (", rot_cam[p[1]], ") - (", rot_cam[p[2]], ")"), character(1))
  pad_tipo <- if (xx$binaria) "parcialmente padronizado (ab/DP_Y)" else "completamente padronizado (ab·DP_X/DP_Y)"
  nomes_ef <- c(paste0("indireto: ", rot_cam), if (n_esp > 1) "indireto total", rot_con,
                if (padronizado) paste0(pad_tipo, ": ", c(rot_cam, if (n_esp > 1) "indireto total")))
  efeitos_de <- function(bm, by, sx, sy) {
    esp <- ind_de(bm, by)
    v <- c(esp, if (n_esp > 1) sum(esp), vapply(pares, function(p) esp[p[1]] - esp[p[2]], numeric(1)))
    if (padronizado) {
      base <- c(esp, if (n_esp > 1) sum(esp))
      v <- c(v, if (xx$binaria) base / sy else base * sx / sy)
    }
    v
  }
  est <- efeitos_de(lapply(ajm, `[[`, "b"), ajy$b, stats::sd(xv), stats::sd(yv))
  names(est) <- nomes_ef

  bt <- .ee_proc_bootstrap(n, reamostras, semente, function(i) {
    bm <- lapply(seq_len(k), function(j) .ee_proc_coef_boot(Xm[[j]], M[, j], i))
    by <- .ee_proc_coef_boot(Xy, yv, i)
    if (is.null(by) || any(vapply(bm, is.null, logical(1)))) return(NULL)
    efeitos_de(bm, by, if (padronizado) stats::sd(xv[i]) else 1, if (padronizado) stats::sd(yv[i]) else 1)
  }, nomes_ef)
  falhas <- sum(!stats::complete.cases(bt))
  rb <- .ee_proc_resumo_boot(est, bt, nivel)

  # efeitos direto e total (teoria t)
  cdir <- tabela[tabela$caminho == "c'", ]; ctot <- tabela[tabela$caminho == "c", ]
  efeito <- rbind(
    data.frame(rb[c("indice", "valor", "ic_inf", "ic_sup")], ep_boot = rb$ep_boot,
               principal = seq_len(nrow(rb)) <= n_esp, stringsAsFactors = FALSE),
    data.frame(indice = c("direto (c')", "total (c)"), valor = c(cdir$b, ctot$b),
               ic_inf = c(cdir$ic_inf, ctot$ic_inf), ic_sup = c(cdir$ic_sup, ctot$ic_sup),
               ep_boot = NA_real_, principal = FALSE, stringsAsFactors = FALSE))

  decisoes <- c(paste0("Mediação ", modelo_nome, " por MQO: ", k + 2, " equações (", paste(c(m, y), collapse = ", "),
                       " e o modelo do efeito total). Coeficientes brutos, na métrica original."),
                pre$decisao, xx$decisao,
                paste0("IC ", round(100 * nivel), "% bootstrap percentílico dos efeitos indiretos, ", reamostras,
                       " reamostras de casos, semente ", if (is.null(semente)) "não fixada" else semente,
                       " (livro: 10.000 e semente 2026; o PROCESS e o curso usam 5.000). Inferência pelo IC: se ele não inclui o zero, há evidência de efeito indireto."),
                if (attr(bt, "esquema") == "sequencial") "n × reamostras grande: índices sorteados uma reamostra por vez (o resultado não coincide com boot::boot)." else
                  "Índices sorteados como no boot::boot (mesma semente dá o mesmo IC).",
                "Efeitos direto e total com IC pela distribuição t (como no PROCESS).",
                "Sem rótulos 'mediação completa' ou 'parcial': relatar magnitude e precisão dos efeitos direto e indireto. O efeito total não precisa ser significativo para haver efeito indireto.")
  if (serial) decisoes <- c(decisoes, paste0("Ordem serial: ", paste(m, collapse = " → "), ". A ordem vem da teoria e muda o modelo; d_jl é o caminho de M_l para M_j controlando X e os mediadores anteriores."))
  if (k > 1 && !serial) decisoes <- c(decisoes, "Mediadores em paralelo: equações separadas por MQO (covariância residual entre os mediadores livre). Foco nos efeitos indiretos específicos; o total pode esconder sinais opostos.")
  if (length(pares)) decisoes <- c(decisoes, "Contrastes = diferença bruta entre efeitos indiretos específicos, com IC bootstrap. Dois efeitos significativos não implicam diferença significativa entre eles.")
  if (padronizado) decisoes <- c(decisoes, paste0("Efeito indireto ", pad_tipo, ", com DP recalculados em cada reamostra; X ", if (xx$binaria) "binário: só Y é padronizado." else "contínuo."))
  if (falhas) decisoes <- c(decisoes, paste0(falhas, " reamostras sem posto completo foram descartadas."))

  # identidade c = c' + soma dos indiretos
  tot_ind <- sum(est[seq_len(n_esp)])
  dif <- ctot$b - (cdir$b + tot_ind)
  decisoes <- c(decisoes, paste0("Identidade c = c' + efeito indireto total: ", round(ctot$b, 4), " = ", round(cdir$b, 4), " + ",
                                 round(tot_ind, 4), " (diferença ", signif(dif, 2), ")."))
  avisos <- c(xx$aviso, .ee_proc_aviso_covs(d, covariaveis))
  if (abs(dif) > 1e-8 * max(1, abs(ctot$b))) avisos <- c(avisos, "A identidade c = c' + indireto não fechou: conferir casos e covariáveis entre as equações.")
  avisos <- c(avisos, .ee_proc_aviso_esquema(bt))
  # Inconsistência só com evidência: indireto e direto com IC sem o zero e sinais opostos.
  # Sinal oposto só na estimativa pontual (c' compatível com zero) é o caso comum.
  esp <- est[seq_len(n_esp)]
  sig_esp <- !(rb$ic_inf[seq_len(n_esp)] <= 0 & rb$ic_sup[seq_len(n_esp)] >= 0)
  sig_dir <- cdir$ic_inf > 0 || cdir$ic_sup < 0
  if (sig_dir && any(sig_esp & sign(esp) != sign(cdir$b)))
    avisos <- c(avisos, "Efeito indireto e efeito direto com sinais opostos, ambos com IC sem o zero (mediação inconsistente): |c'| pode superar |c|, c pode ser não significativo e a razão ab/c não é interpretável.")
  if (n_esp > 1 && any(sig_esp & esp > 0) && any(sig_esp & esp < 0))
    avisos <- c(avisos, "Efeitos indiretos específicos com sinais opostos e IC sem o zero: o indireto total pode ficar perto de zero com mecanismos ativos, e o contraste bruto mistura direção e magnitude.")
  avisos <- c(avisos, .ee_proc_frase_delineamento(delineamento, idioma))

  sobel <- NULL
  if (k == 1) {
    a <- ajm[[1]]$b[[2]]; b <- ajy$b[[3]]
    sa <- sqrt(ajm[[1]]$V[2, 2]); sb <- sqrt(ajy$V[3, 3])
    z <- a * b / sqrt(b^2 * sa^2 + a^2 * sb^2)
    sobel <- data.frame(z = z, p = 2 * stats::pnorm(-abs(z)), ep = sqrt(b^2 * sa^2 + a^2 * sb^2))
    decisoes <- c(decisoes, paste0("Teste de Sobel (z = ", round(z, 3), ", p = ", signif(sobel$p, 3),
                                   ") só como comparação histórica em $sobel; não é a inferência."))
    if (proporcao_mediada) {
      pm <- est[[1]] / ctot$b
      if (sign(est[[1]]) != sign(ctot$b)) {
        pm <- NA_real_
        avisos <- c(avisos, "Proporção mediada não calculada: efeito indireto e efeito total com sinais opostos.")
      }
      efeito <- rbind(efeito, data.frame(indice = "proporção mediada (ab/c)", valor = pm, ic_inf = NA, ic_sup = NA,
                                         ep_boot = NA, principal = FALSE))
      avisos <- c(avisos, "Proporção mediada = ab/c (nunca c'/c): instável com c pequeno e pede n grande (MacKinnon et al., 1995); relatar como descrição, sem IC.")
    }
  } else if (proporcao_mediada) {
    avisos <- c(avisos, "Proporção mediada só é calculada na mediação simples.")
  }

  modelos <- do.call(rbind, lapply(c(ajm, list(ajy, ajt)), function(a) data.frame(r2 = a$r2, r2_aj = a$r2_aj,
    F = a$F[["F"]], gl1 = a$F[["gl1"]], gl2 = a$F[["gl2"]], p = a$F[["p"]])))
  modelos <- cbind(desfecho = c(m, y, paste0(y, " (efeito total)")), modelos)
  pres <- .ee_proc_pressupostos(stats::setNames(c(ajm, list(ajy)), c(m, y)), ep, idioma)

  # frase APA
  sep <- .ee_sep(decimal); eps <- .ee_proc_tr(idioma, "EP", "SE")
  boot_txt <- .ee_proc_tr(idioma,
    paste0("bootstrap percentílico, ", .ee_proc_int(reamostras, decimal, idioma), " reamostras", if (!is.null(semente)) paste0(", semente ", semente) else ""),
    paste0("percentile bootstrap, ", .ee_proc_int(reamostras, decimal, idioma), " resamples", if (!is.null(semente)) paste0(", seed ", semente) else ""))
  ampliou <- FALSE
  num <- function(v, grupo = v) { cs <- .ee_proc_casas(grupo, casas); if (cs > casas) ampliou <<- TRUE; .pr_num(v, cs, TRUE, decimal) }
  ic_f <- function(lo, hi) .ee_proc_ic(lo, hi, num, decimal, nivel, idioma)
  ef_txt <- function(j, simb = "ab", extra = "") paste0(simb, " = ", num(rb$valor[j]), " (", ic_f(rb$ic_inf[j], rb$ic_sup[j]), extra, ")")
  cam_txt <- function(lin, simb) .ee_junta(paste0(simb, " = ", num(lin$b, c(lin$b, lin$ep))),
                                           paste0(eps, " = ", num(lin$ep, c(lin$b, lin$ep))),
                                           .ee_estat("t", lin$t, lin$gl, decimal), .pr_p(lin$p, decimal), decimal = decimal)
  if (k == 1) {
    la <- tabela[tabela$caminho == "a", ]; lb <- tabela[tabela$caminho == "b", ]
    apa <- .ee_proc_tr(idioma,
      paste0("O efeito indireto de ", x, " sobre ", y, " por meio de ", m, " foi ", ef_txt(1, extra = paste0("; ", boot_txt)), ", com a = ",
             num(la$b), " (", .pr_p(la$p, decimal), ") e b = ", num(lb$b), " (",
             .pr_p(lb$p, decimal), ")."),
      paste0("The indirect effect of ", x, " on ", y, " through ", m, " was ", ef_txt(1, extra = paste0("; ", boot_txt)), ", with a = ",
             num(la$b), " (", .pr_p(la$p, decimal), ") and b = ", num(lb$b), " (",
             .pr_p(lb$p, decimal), ")."))
  } else {
    itens <- vapply(seq_len(n_esp), function(j) paste0(rot_cam[j], ": ", ef_txt(j)), character(1))
    apa <- .ee_proc_tr(idioma,
      paste0("Os efeitos indiretos específicos de ", x, " sobre ", y, " foram: ", paste(itens, collapse = "; "),
             ". O efeito indireto total foi ", ef_txt(n_esp + 1, "ab total"), ". Intervalos por ", boot_txt, "."),
      paste0("The specific indirect effects of ", x, " on ", y, " were: ", paste(itens, collapse = "; "),
             ". The total indirect effect was ", ef_txt(n_esp + 1, "total ab"), ". Intervals by ", boot_txt, "."))
    if (length(pares)) {
      con <- vapply(seq_along(pares), function(q) { j <- n_esp + 1 + q
        paste0(sub("^contraste: ", "", rot_con[q]), ": ", .ee_proc_tr(idioma, "diferença", "difference"), " = ",
               num(rb$valor[j]), " (", ic_f(rb$ic_inf[j], rb$ic_sup[j]), ")") },
        character(1))
      apa <- c(apa, .ee_proc_tr(idioma, paste0("Contrastes entre efeitos indiretos: ", paste(con, collapse = "; "), "."),
                                paste0("Contrasts between indirect effects: ", paste(con, collapse = "; "), ".")))
    }
  }
  apa <- c(apa, .ee_proc_tr(idioma,
    paste0("O efeito direto foi ", cam_txt(cdir, "c'"), ". O efeito total foi ", cam_txt(ctot, "c"), "."),
    paste0("The direct effect was ", cam_txt(cdir, "c'"), ". The total effect was ", cam_txt(ctot, "c"), ".")))
  if (ampliou) decisoes <- c(decisoes, .ee_proc_nota_casas(casas))

  .ee_resultado(paste0("Mediação ", modelo_nome, ": ", x, " → ", paste(m, collapse = if (serial) " → " else ", "), " → ", y),
                n = n, n_excluidos = pre$n_excluidos, tabela = tabela, efeito = efeito, pressupostos = pres,
                decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(modelos = modelos, efeitos_bootstrap = rb, sobel = sobel, reamostras = bt,
                             identidade = c(total = ctot$b, direto = cdir$b, indireto_total = tot_ind, diferenca = dif),
                             config = list(tipo = if (k == 1) "simples" else tipo, reamostras = reamostras, semente = semente,
                                           nivel = nivel, ep = ep, delineamento = delineamento)))
}

# ---------------------------------------------------------------- mediação moderada (7, 14, 15)
mediacao_moderada <- function(dados, x, m, y, w, modelo = 7, covariaveis = NULL,
                              valores_w = c("media_dp", "percentis"), centralizar = FALSE,
                              reamostras = 10000, semente = 2026, nivel = .95,
                              ep = c("classico", "HC3"), johnson_neyman = TRUE, jn_prop_min = .05,
                              delineamento = c("transversal", "longitudinal", "experimental"),
                              casas = 3, decimal = ".", idioma = "pt") {
  ep <- match.arg(ep); delineamento <- match.arg(delineamento)
  modelo <- as.integer(modelo)
  if (!modelo %in% c(7L, 14L, 15L)) stop("modelo precisa ser 7 (W modera a), 14 (W modera b) ou 15 (W modera b e c')")
  if (length(m) != 1) stop("mediacao_moderada aceita um mediador")
  .ee_proc_checa_papeis(x, m, y, w, covariaveis)
  pre <- .ee_proc_dados(dados, c(x, m, y, w, covariaveis), numericas = c(m, y))
  d <- pre$d; n <- nrow(d)
  .ee_proc_checa_y(d[[y]], y)
  if (length(unique(d[[m]])) < 3) stop("o mediador ", m, " precisa ser contínuo (mais de dois valores)")
  xx <- .ee_proc_x(d[[x]], x)
  ww <- .ee_proc_x(d[[w]], w, "Se for grupo, esta função não cobre W com 3 ou mais categorias: rode com dois grupos por vez (W de dois níveis vira 0/1).")
  C <- .ee_proc_cov(d, covariaveis)
  mv <- d[[m]]; yv <- d[[y]]
  cx <- if (centralizar && !xx$binaria && modelo != 14L) mean(xx$valor) else 0
  cm <- if (centralizar && modelo != 7L) mean(mv) else 0
  cw <- if (centralizar && !ww$binaria) mean(ww$valor) else 0
  xc <- xx$valor - cx; mc <- mv - cm; wc <- ww$valor - cw
  pxw <- paste0(x, " × ", w); pmw <- paste0(m, " × ", w)

  if (modelo == 7L) {
    Xm <- cbind(1, xc, wc, xc * wc, C); colnames(Xm) <- c("(Intercept)", x, w, pxw, colnames(C))
    Xy <- cbind(1, xc, mc, C); colnames(Xy) <- c("(Intercept)", x, m, colnames(C))
  } else {
    Xm <- cbind(1, xc, C); colnames(Xm) <- c("(Intercept)", x, colnames(C))
    Xy <- if (modelo == 14L) cbind(1, xc, mc, wc, mc * wc, C) else cbind(1, xc, mc, wc, xc * wc, mc * wc, C)
    colnames(Xy) <- c("(Intercept)", x, m, w, if (modelo == 15L) pxw, pmw, colnames(C))
  }
  ajm <- .ee_proc_mqo(Xm, mv, ep); ajy <- .ee_proc_mqo(Xy, yv, ep)
  cov_cam <- function(cam) { if (!is.null(C)) for (cc in colnames(C)) cam[[cc]] <- "covariável"; cam }
  camm <- if (modelo == 7L) list(a1 = "a1", a2 = "a2", a3 = "a3") else list(a = "a")
  names(camm) <- if (modelo == 7L) c(x, w, pxw) else x
  camy <- if (modelo == 7L) list("c'", "b") else if (modelo == 14L) list("c'", "b1", "b2", "b3") else list("c'1", "b1", "b2", "c'3", "b3")
  names(camy) <- colnames(Xy)[2:(1 + length(camy))]
  tabela <- rbind(.ee_proc_coefs(ajm, nivel, m, cov_cam(camm)), .ee_proc_coefs(ajy, nivel, y, cov_cam(camy)))

  vw <- .ee_proc_valores_w(ww$valor, valores_w, w, idioma)
  wv_mod <- vw$valores - cw
  # índice e indiretos condicionais a partir dos coeficientes
  efeitos_de <- function(bm, by) {
    if (modelo == 7L) {
      ind <- (bm[[2]] + bm[[4]] * wv_mod) * by[[3]]; idx <- bm[[4]] * by[[3]]
    } else {
      bmw <- by[[if (modelo == 14L) 5 else 6]]
      ind <- bm[[2]] * (by[[3]] + bmw * wv_mod); idx <- bm[[2]] * bmw
    }
    c(idx, ind)
  }
  rot_w <- paste0(w, " = ", format(round(vw$valores, 4)), " (", vw$rotulos, ")")
  nomes_ef <- c("índice de mediação moderada", paste0("indireto condicional em ", rot_w))
  est <- stats::setNames(efeitos_de(ajm$b, ajy$b), nomes_ef)
  bt <- .ee_proc_bootstrap(n, reamostras, semente, function(i) {
    bm <- .ee_proc_coef_boot(Xm, mv, i); by <- .ee_proc_coef_boot(Xy, yv, i)
    if (is.null(bm) || is.null(by)) return(NULL)
    efeitos_de(bm, by)
  }, nomes_ef)
  falhas <- sum(!stats::complete.cases(bt))
  rb <- .ee_proc_resumo_boot(est, bt, nivel)
  efeito <- data.frame(rb[c("indice", "valor", "ic_inf", "ic_sup")], ep_boot = rb$ep_boot, principal = TRUE,
                       stringsAsFactors = FALSE)

  # sondagem do caminho moderado (teoria t) e Johnson-Neyman
  if (modelo == 7L) { aj_s <- ajm; i1 <- 2; i3 <- 4; cam_mod <- paste0("a (", x, " → ", m, ")") } else {
    aj_s <- ajy; i1 <- 3; i3 <- if (modelo == 14L) 5 else 6; cam_mod <- paste0("b (", m, " → ", y, ")") }
  sond <- function(aj, a1, a3, rotulo) {
    out <- lapply(wv_mod, function(v) { L <- numeric(length(aj$b)); L[a1] <- 1; L[a3] <- v; .ee_proc_comb(aj, L, nivel) })
    data.frame(efeito_de = rotulo, w_nivel = vw$rotulos, w = vw$valores, do.call(rbind, out), stringsAsFactors = FALSE)
  }
  posthoc <- sond(aj_s, i1, i3, paste0("caminho ", cam_mod))
  if (modelo == 15L) posthoc <- rbind(posthoc, sond(ajy, 2, 5, paste0("efeito direto c' (", x, " → ", y, ")")))
  testes <- list(.ee_proc_teste_produto(aj_s, i3, colnames(aj_s$X)[i3], if (modelo == 7L) Xm else Xy, if (modelo == 7L) mv else yv, ep))
  if (modelo == 15L) testes <- c(testes, list(.ee_proc_teste_produto(ajy, 5, pxw, Xy, yv, ep)))
  testes <- do.call(rbind, testes)

  avisos <- c(xx$aviso, ww$aviso, .ee_proc_aviso_covs(d, covariaveis), vw$avisos, .ee_proc_aviso_esquema(bt))
  jn <- NULL
  if (johnson_neyman && !ww$binaria) {
    jn <- .ee_proc_jn(aj_s$b[[i1]], aj_s$b[[i3]], aj_s$V[i1, i1], aj_s$V[i1, i3], aj_s$V[i3, i3], aj_s$gl, nivel, wc, cw, jn_prop_min)
    if (any(jn$pequenas)) avisos <- c(avisos, paste0("Região de Johnson-Neyman do caminho ", cam_mod, " com menos de ",
                                                     round(100 * jn_prop_min), "% dos casos."))
  }

  sig_idx <- !(rb$ic_inf[1] <= 0 && rb$ic_sup[1] >= 0)
  desc_mod <- switch(as.character(modelo), `7` = paste0("W modera o caminho a (", x, " → ", m, ")"),
                     `14` = paste0("W modera o caminho b (", m, " → ", y, ")"),
                     `15` = paste0("W modera o caminho b (", m, " → ", y, ") e o efeito direto (", x, " → ", y, ")"))
  formula_ind <- if (modelo == 7L) "(a1 + a3·W)·b; índice = a3·b" else "a·(b1 + b3·W); índice = a·b3"
  decisoes <- c(paste0("Modelo ", modelo, " do PROCESS (Hayes, 2022): ", desc_mod, ". Duas equações por MQO; efeito indireto condicional = ",
                       formula_ind, "."),
                pre$decisao, xx$decisao, ww$decisao,
                paste0("Índice de mediação moderada e efeitos indiretos condicionais com IC ", round(100 * nivel),
                       "% bootstrap percentílico, ", reamostras, " reamostras, semente ", if (is.null(semente)) "não fixada" else semente,
                       "; os valores de W ficam fixos nas reamostras (como no PROCESS). O livro pede os indiretos condicionais com incerteza e um índice da variação ao longo de W."),
                paste0("O índice é a inclinação do efeito indireto em W: a diferença entre dois indiretos condicionais é o índice vezes a distância entre os valores de W. ",
                       if (sig_idx) "IC do índice sem o zero: o efeito indireto varia com W." else "IC do índice com o zero: não há evidência de que o efeito indireto varie com W, mesmo que algum indireto condicional seja significativo."),
                "Interação significativa num caminho é condição necessária, mas não suficiente, para mediação moderada.",
                switch(vw$metodo, media_dp = "W sondada na média e a ±1 DP (livro); percentis 16/50/84 com valores_w = \"percentis\".",
                       percentis = "W sondada nos percentis 16, 50 e 84 (PROCESS; quantile tipo 7 do R).",
                       binaria = "W binária: efeitos nos dois valores.", "W sondada nos valores informados."),
                if (centralizar) paste0("Centralizadas pela média as variáveis contínuas que formam produtos (", paste(c(if (!xx$binaria && modelo != 14L) x, if (modelo != 7L) m, if (!ww$binaria) w), collapse = ", "),
                                        "); só muda a leitura dos termos de primeira ordem, não índice nem indiretos condicionais.") else
                  "Variáveis na escala original: termos de primeira ordem são efeitos com o moderador em zero.")
  if (modelo == 15L) decisoes <- c(decisoes, "Efeitos diretos condicionais (c'1 + c'3·W) pela distribuição t, em $posthoc.")
  decisoes <- c(decisoes, .ee_proc_nota_sem_ic("ΔR² das interações", "ΔR² of the interactions", idioma))
  if (!is.null(jn)) decisoes <- c(decisoes, paste0("Johnson-Neyman do caminho ", cam_mod, ": ",
    if (is.null(jn$pontos)) "nenhuma fronteira na faixa observada de W." else paste0("fronteira(s) em W = ", paste(round(jn$pontos$w, 4), collapse = " e "), ".")))
  if (falhas) decisoes <- c(decisoes, paste0(falhas, " reamostras sem posto completo foram descartadas."))
  if (ep == "HC3") decisoes <- c(decisoes, "EP robusto HC3 nos coeficientes e nos testes das interações (o bootstrap não usa EP).")
  avisos <- c(avisos, .ee_proc_frase_delineamento(delineamento, idioma),
              "Processo condicional encadeia mais suposições causais que a mediação simples e pede amostra maior.")

  pres <- .ee_proc_pressupostos(stats::setNames(list(ajm, ajy), c(m, y)), ep, idioma)
  modelos <- do.call(rbind, lapply(list(ajm, ajy), function(a) data.frame(r2 = a$r2, r2_aj = a$r2_aj,
    F = a$F[["F"]], gl1 = a$F[["gl1"]], gl2 = a$F[["gl2"]], p = a$F[["p"]])))
  modelos <- cbind(desfecho = c(m, y), modelos)

  # frase APA
  sep <- .ee_sep(decimal)
  boot_txt <- .ee_proc_tr(idioma,
    paste0("bootstrap percentílico, ", .ee_proc_int(reamostras, decimal, idioma), " reamostras", if (!is.null(semente)) paste0(", semente ", semente) else ""),
    paste0("percentile bootstrap, ", .ee_proc_int(reamostras, decimal, idioma), " resamples", if (!is.null(semente)) paste0(", seed ", semente) else ""))
  ampliou <- FALSE
  num <- function(v, grupo = v) { cs <- .ee_proc_casas(grupo, casas); if (cs > casas) ampliou <<- TRUE; .pr_num(v, cs, TRUE, decimal) }
  ic_txt <- function(j) .ee_proc_ic(rb$ic_inf[j], rb$ic_sup[j], num, decimal, nivel, idioma)
  frase_int <- function(t1, alvo) {
    int_txt <- .ee_junta(paste0("b = ", num(t1$b)), .ee_estat("t", t1$t, t1$gl, decimal), .pr_p(t1$p, decimal),
                         .ee_efeito_txt("ΔR²", t1$delta_r2, casas = 3, decimal = decimal, idioma = idioma, nivel = nivel), decimal = decimal)
    .ee_proc_tr(idioma,
      paste0("A interação ", t1$termo, " na predição de ", alvo, if (t1$p < 1 - nivel) " foi significativa, " else " não foi significativa, ", int_txt, "."),
      paste0("The ", t1$termo, " interaction predicting ", alvo, if (t1$p < 1 - nivel) " was significant, " else " was not significant, ", int_txt, "."))
  }
  frases_int <- vapply(seq_len(nrow(testes)), function(j) frase_int(testes[j, ], if (modelo == 7L) m else y), character(1))
  itens <- vapply(seq_along(vw$valores), function(j) paste0(num(rb$valor[j + 1]), " (", ic_txt(j + 1), ")",
    .ee_proc_tr(idioma, " em ", " at "), w, " = ", .pr_num(vw$valores[j], casas, TRUE, decimal),
    if (vw$metodo %in% c("media_dp", "percentis")) paste0(" (", vw$rotulos[j], ")") else ""), character(1))
  apa <- .ee_proc_tr(idioma,
    c(paste0("O índice de mediação moderada foi ", num(rb$valor[1]), " (", ic_txt(1), "; ", boot_txt, ")",
             if (sig_idx) ", indicando que o efeito indireto de " else ", sem evidência de que o efeito indireto de ",
             x, " sobre ", y, " por meio de ", m, if (sig_idx) " varia com " else " varie com ", w, "."),
      paste0("Os efeitos indiretos condicionais foram ", paste(itens, collapse = "; "), ".")),
    c(paste0("The index of moderated mediation was ", num(rb$valor[1]), " (", ic_txt(1), "; ", boot_txt, ")",
             if (sig_idx) ", indicating that the indirect effect of " else ", providing no evidence that the indirect effect of ",
             x, " on ", y, " through ", m, if (sig_idx) " depends on " else " depends on ", w, "."),
      paste0("The conditional indirect effects were ", paste(itens, collapse = "; "), ".")))
  apa <- c(frases_int, apa)
  if (ampliou) decisoes <- c(decisoes, .ee_proc_nota_casas(casas))

  .ee_resultado(paste0("Mediação moderada (modelo ", modelo, "): ", x, " → ", m, " → ", y, ", moderador ", w),
                n = n, n_excluidos = pre$n_excluidos, tabela = tabela, efeito = efeito, pressupostos = pres,
                posthoc = posthoc, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(modelos = modelos, testes_interacao = testes, johnson_neyman = jn,
                             efeitos_bootstrap = rb, reamostras = bt, centros = c(x = cx, m = cm, w = cw),
                             config = list(modelo = modelo, reamostras = reamostras, semente = semente,
                                           nivel = nivel, ep = ep, delineamento = delineamento)))
}

# Teste de um produto: b, t, p (EP do ajuste), ΔR² ao retirar o produto e F (Wald; com
# EP clássico igual ao F da mudança de R²).
.ee_proc_teste_produto <- function(aj, i, termo, X, yv, ep) {
  aj0 <- .ee_proc_mqo(X[, -i, drop = FALSE], yv, ep)
  wt <- .ee_proc_wald(aj, i)
  se <- sqrt(aj$V[i, i])
  data.frame(termo = termo, b = aj$b[[i]], ep = se, t = aj$b[[i]] / se, gl = aj$gl,
             p = wt[["p"]], delta_r2 = aj$r2 - aj0$r2, F = wt[["F"]], gl1 = 1, gl2 = aj$gl,
             stringsAsFactors = FALSE)
}
