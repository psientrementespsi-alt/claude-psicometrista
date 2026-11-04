# est_modelos.R: modelos avançados da skill psicometria-r. GLM (normal, gama, binomial,
# Poisson, binomial negativa), GEE, multinível linear e logístico e meta-análise.
# Carregar por scripts/carregar.R (usa est_comum.R e as auxiliares .pr_* de psicometria.R).
#
# O livro do usuário não cobre estes modelos: o GEE aparece como "horizonte, fora do
# escopo" e multinível e meta-análise não aparecem. As regras vêm das
# notas do curso (GLM e GEE; multinível; meta-análise), conferidas contra a literatura
# com DOI verificado; os erros apontados nas extrações do curso não são reproduzidos.
# Detalhes, âncoras reproduzidas e divergências: references/modelos.md.
#
# Funções públicas (todas as de modelo aceitam fatores = c("var") para grupo gravado como
# número; sem isso, preditor numérico com até 5 valores entra como contínuo, com aviso)
#   escolher_familia(dados, vd)
#   modelo_glm(dados, vd, preditores, familia, ligacao, referencia, fatores, escala,
#              exposicao, comparar, ajuste, nivel, decimal, idioma)
#   comparar_familias(dados, vd, preditores, familias, ligacoes, fatores, decimal, idioma)
#   modelo_gee(dados, vd, preditores, id, tempo, familia, ligacao, correlacao, ep,
#              referencia, fatores, comparar, ajuste, nivel, decimal, idioma)
#   comparar_correlacoes_gee(dados, vd, preditores, id, tempo, familia, ligacao,
#              correlacoes, fatores, decimal, idioma)
#   centralizar_multinivel(dados, vars, cluster, tipo = c("geral", "grupo"))
#   multinivel_linear(dados, vd, cluster, nivel1, nivel2, inclinacao, interacao,
#              referencias, fatores, modelo_final, estimador_final, nivel, decimal, idioma)
#   multinivel_logistico(dados, vd, cluster, nivel1, nivel2, inclinacao, interacao,
#              referencias, fatores, evento, modelo_final, nAGQ, nivel, decimal, idioma)
#   calcular_efeitos_meta(dados, medida, estudo, media1, dp1, n1, media2, dp2, n2,
#              eventos1, eventos2, r, n)
#   meta_analise(dados, yi, vi, estudo, medida, metodo, hartung_knapp, sensibilidade,
#              vies, k_min_vies, nivel, decimal, idioma)
#   meta_subgrupos(dados, subgrupo, yi, vi, estudo, medida, metodo, hartung_knapp,
#              nivel, decimal, idioma)
#   meta_regressao(dados, moderadores, yi, vi, medida, metodo, hartung_knapp, fatores,
#              permutacao, semente, reamostras, nivel, decimal, idioma)
#   modelo_ordinal(dados, vd, preditores, ligacao, referencia, fatores, nivel, decimal, idioma)   (v1.6)
#   efeito_contextual(resultado, var, nivel, decimal, idioma)   (v1.6)
#   modelo_multinomial(dados, vd, preditores, referencia_vd, referencia, fatores, nivel,
#              decimal, idioma)   (v1.6)

# ================================================================ auxiliares gerais
.ee_mod_t <- function(idioma, pt, en) if (identical(idioma, "en")) en else pt

.ee_mod_vars <- function(preditores) {
  if (!length(preditores)) return(character())
  unique(all.vars(stats::as.formula(paste("~", paste(preditores, collapse = " + ")))))
}

.ee_mod_nome_ok <- function(v) {
  ruins <- v[make.names(v) != v]
  if (length(ruins)) stop("nomes de variável precisam ser sintáticos no R (sem espaço ou acento no início): ",
                          paste(ruins, collapse = ", "), ". Renomeie a coluna antes.")
}

# Converte texto, lógico e as variáveis em `fatores` (grupo gravado como 1/2/3) em fator,
# aplica as categorias de referência pedidas e tira níveis vazios.
# referencia: lista ou vetor nomeado, ex.: list(grupo = "controle").
.ee_mod_fatores <- function(d, vars, referencia = NULL, fatores = NULL) {
  if (length(fatores)) {
    fora <- setdiff(fatores, vars)
    if (length(fora)) stop("fatores cita variável que não está no modelo: ", paste(fora, collapse = ", "))
    for (v in fatores) d[[v]] <- factor(d[[v]])
  }
  for (v in vars) {
    if (is.character(d[[v]]) || is.logical(d[[v]])) d[[v]] <- factor(d[[v]])
    if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  }
  if (length(referencia)) {
    for (v in names(referencia)) {
      if (!v %in% names(d) || !is.factor(d[[v]])) stop("referência pedida para ", v, ", que não é fator no modelo")
      if (!referencia[[v]] %in% levels(d[[v]])) stop("o nível '", referencia[[v]], "' não existe em ", v)
      d[[v]] <- stats::relevel(d[[v]], ref = as.character(referencia[[v]]))
    }
  }
  d
}

# Preditor numérico com poucos valores distintos (grupo exportado do SPSS/JASP como 1/2/3)
# entra como contínuo, com um só coeficiente linear: avisa.
# Mesma regra dos módulos de associação e processos: preditor numérico com 3 a 7
# valores distintos, todos inteiros, avisa (código de grupo vindo do SPSS ou do JASP);
# 0/1 é a codificação que o livro manda para a dicotômica e não avisa.
.ee_mod_aviso_numericos <- function(d, vars, isentos = character()) {
  out <- character()
  for (v in setdiff(vars, isentos)) {
    x <- d[[v]]
    if (!is.numeric(x)) next
    u <- sort(unique(x[!is.na(x)])); k <- length(u)
    if (k >= 3 && k <= 7 && all(abs(u - round(u)) < 1e-8))
      out <- c(out, paste0(v, " é numérica com só ", k, " valores inteiros distintos (", paste(u, collapse = ", "),
                           ") e entra como contínua (um coeficiente linear, 1 gl). Se for grupo, use fatores = \"", v,
                           "\"; se for escala, o aviso pode ser ignorado."))
  }
  out
}

.ee_mod_formula <- function(vd, preditores, extra = NULL) {
  rhs <- if (length(preditores)) paste(preditores, collapse = " + ") else "1"
  if (length(extra)) rhs <- paste(rhs, extra, sep = " + ")
  stats::as.formula(paste(vd, "~", rhs))
}

.ee_mod_rot_familia <- function(f, idioma) {
  switch(f, normal = "normal",
         gama = .ee_mod_t(idioma, "gama", "gamma"),
         binomial = "binomial", poisson = "Poisson",
         binomial_negativa = .ee_mod_t(idioma, "binomial negativa", "negative binomial"))
}
.ee_mod_rot_ligacao <- function(l, idioma) {
  switch(l, identity = .ee_mod_t(idioma, "identidade", "identity"),
         log = "log", logit = "logit", inverse = .ee_mod_t(idioma, "inversa", "inverse"), l)
}
.ee_mod_rot_exp <- function(familia, ligacao, idioma) {
  if (ligacao == "logit") return("OR")
  if (ligacao != "log") return(NA_character_)
  switch(familia, poisson = "IRR", binomial_negativa = "IRR", binomial = "RR",
         .ee_mod_t(idioma, "razão de médias", "ratio of means"))
}

# Teste de Wald Tipo III por termo: W = b' (L V L')^-1 b sobre as colunas do termo.
# Com contrastes soma-zero (contr.sum), é o "Tests of Model Effects" Tipo III do SPSS.
.ee_mod_wald3 <- function(b, V, assign, rotulos) {
  ok <- !is.na(b)
  b <- b[ok]; assign <- assign[ok]; V <- V[ok, ok, drop = FALSE]
  linhas <- lapply(seq_along(rotulos), function(j) {
    idx <- which(assign == j)
    if (!length(idx)) return(NULL)
    W <- as.numeric(crossprod(b[idx], solve(V[idx, idx, drop = FALSE], b[idx])))
    data.frame(termo = rotulos[j], wald_qui2 = W, gl = length(idx),
               p = stats::pchisq(W, length(idx), lower.tail = FALSE))
  })
  do.call(rbind, linhas)
}

# Tabela de coeficientes de Wald (z), com exp(B) e IC para ligação log ou logit.
.ee_mod_coefs <- function(b, V, nivel = .95, rot_exp = NA) {
  ok <- !is.na(b); b <- b[ok]; V <- V[ok, ok, drop = FALSE]
  ep <- sqrt(diag(V)); z <- stats::qnorm(1 - (1 - nivel) / 2)
  tab <- data.frame(termo = names(b), B = unname(b), ep = unname(ep),
                    ic_inf = unname(b - z * ep), ic_sup = unname(b + z * ep),
                    wald_qui2 = unname((b / ep)^2), gl = 1, p = unname(2 * stats::pnorm(-abs(b / ep))))
  if (!is.na(rot_exp)) {
    tab$exp_B <- exp(tab$B); tab$exp_ic_inf <- exp(tab$ic_inf); tab$exp_ic_sup <- exp(tab$ic_sup)
  }
  tab
}

# Normaliza o summary do emmeans (os nomes das colunas mudam com a escala). Médias:
# nível do primeiro fator, "por" = níveis dos demais. Pares: contraste e "por".
.ee_mod_emm_df <- function(s, fatores) {
  s <- as.data.frame(s)
  est <- intersect(c("emmean", "response", "prob", "rate", "estimate", "ratio", "odds.ratio"), names(s))[1]
  lo <- intersect(c("asymp.LCL", "lower.CL"), names(s))[1]
  hi <- intersect(c("asymp.UCL", "upper.CL"), names(s))[1]
  zc <- intersect(c("z.ratio", "t.ratio"), names(s))[1]
  par <- "contrast" %in% names(s)
  rot <- if (par) as.character(s$contrast) else as.character(s[[fatores[1]]])
  outros <- intersect(fatores[-1], names(s))
  por <- if (length(outros)) apply(s[, outros, drop = FALSE], 1, function(x) paste(trimws(x), collapse = " ")) else rep("", nrow(s))
  data.frame(nivel_ou_par = rot, por = por,
             estimativa = s[[est]], ep = s$SE,
             ic_inf = if (!is.na(lo)) s[[lo]] else NA_real_, ic_sup = if (!is.na(hi)) s[[hi]] else NA_real_,
             z = if (!is.na(zc)) s[[zc]] else NA_real_, p = if ("p.value" %in% names(s)) s$p.value else NA_real_,
             stringsAsFactors = FALSE)
}

# Casos por linha de médias marginais: linhas do banco com o nível (e os níveis do "por")
# daquela linha; no GEE, `id` conta sujeitos distintos, não observações.
.ee_mod_emm_n <- function(s, d, partes, id = NULL) {
  s <- as.data.frame(s)
  cols <- intersect(partes, intersect(names(s), names(d)))
  if (!length(cols)) return(rep(NA_real_, nrow(s)))
  vapply(seq_len(nrow(s)), function(i) {
    ok <- rep(TRUE, nrow(d))
    for (v in cols) ok <- ok & as.character(d[[v]]) == as.character(s[[v]][i])
    ok <- ok & !is.na(ok)
    if (is.null(id)) sum(ok) else length(unique(d[[id]][ok]))
  }, numeric(1))
}

# Médias marginais estimadas (escala da resposta) e pares com correção, por termo.
# Termo "A" = médias de A (média sobre os demais fatores); termo "A:B" = A dentro de B.
# A coluna n (casos do nível ou da célula; sujeitos no GEE; NA nos pares) serve à
# supressão de grupo pequeno no print (print.pr_modelo): com n = 1 a média marginal do
# nível é o escore da pessoa.
.ee_mod_emm <- function(fit, termos, V, d, ajuste = "bonferroni", nivel = .95, id = NULL) {
  out <- list()
  for (tt in termos) {
    partes <- strsplit(tt, ":", fixed = TRUE)[[1]]
    specs <- if (length(partes) == 1) stats::as.formula(paste("~", partes)) else
      stats::as.formula(paste("~", partes[1], "|", paste(partes[-1], collapse = " * ")))
    em <- suppressMessages(emmeans::emmeans(fit, specs, vcov. = V, data = d, df = Inf, type = "response"))
    sm <- summary(em, infer = c(TRUE, FALSE), level = nivel)
    med <- .ee_mod_emm_df(sm, partes)
    n_med <- tryCatch(.ee_mod_emm_n(sm, d, partes, id), error = function(e) rep(NA_real_, nrow(med)))
    med <- cbind(med[, c("nivel_ou_par", "por")], n = n_med, med[, setdiff(names(med), c("nivel_ou_par", "por"))])
    pr <- emmeans::contrast(em, "pairwise", adjust = ajuste)
    par <- .ee_mod_emm_df(summary(pr, infer = c(TRUE, TRUE), level = nivel, adjust = ajuste), partes)
    par <- cbind(par[, c("nivel_ou_par", "por")], n = NA_real_, par[, setdiff(names(par), c("nivel_ou_par", "por"))])
    out[[tt]] <- rbind(cbind(termo = tt, tipo = "média marginal", med),
                       cbind(termo = tt, tipo = paste0("par (", ajuste, ")"), par))
  }
  if (!length(out)) return(NULL)
  r <- do.call(rbind, out); rownames(r) <- NULL; r
}

# Níveis (ou células) com menos de .ee_n_min_print casos nas médias marginais: texto do
# aviso. As comparações e o coeficiente desse nível continuam no print e na frase.
.ee_mod_aviso_pequenos <- function(ph, unidade = "casos") {
  if (is.null(ph) || !"n" %in% names(ph)) return(character())
  peq <- ph[ph$tipo == "média marginal" & !is.na(ph$n) & ph$n < .ee_n_min_print, , drop = FALSE]
  if (!nrow(peq)) return(character())
  rot <- paste0(peq$termo, " = ", peq$nivel_ou_par, ifelse(nzchar(peq$por), paste0(" [", peq$por, "]"), ""), ", n = ", peq$n)
  paste0("Nível com menos de ", .ee_n_min_print, " ", unidade, " (", paste(unique(rot), collapse = "; "), "): a média marginal ",
         "e o IC desse nível saem do print, por privacidade (os valores ficam em $posthoc). As comparações e o coeficiente ",
         "desse nível continuam no print e na frase, e, somados à média do outro nível, devolvem o escore de quem está ",
         "sozinho no grupo. Antes de relatar, agrupar a categoria ou omitir essas comparações.")
}

# Print de modelo_glm() e modelo_gee(): a tabela de comparações passa pela supressão de
# grupo pequeno da base (via tabelas_extra, que tem a coluna n), e o resto segue
# print.pr_resultado. O objeto não muda: $posthoc guarda todos os valores.
print.pr_modelo <- function(x, ...) {
  y <- x
  class(y) <- setdiff(class(x), "pr_modelo")
  if (is.data.frame(x$posthoc) && "n" %in% names(x$posthoc)) {
    rot <- if (is.character(x$titulo_posthoc)) sub(":\\s*$", "", x$titulo_posthoc) else "Comparações"
    y$tabelas_extra <- c(stats::setNames(list(x$posthoc), rot), if (is.list(x$tabelas_extra)) x$tabelas_extra)
    y$posthoc <- NULL
  }
  print(y, ...)
  invisible(x)
}

.ee_mod_classe <- function(r) { class(r) <- c("pr_modelo", class(r)); r }

# Termos para pós-hoc: efeitos principais dos fatores e interações só entre fatores.
.ee_mod_termos_fator <- function(rotulos, d) {
  rotulos[vapply(rotulos, function(tt) all(vapply(strsplit(tt, ":", fixed = TRUE)[[1]],
                                                    function(v) v %in% names(d) && is.factor(d[[v]]), logical(1))),
                 logical(1))]
}

.ee_mod_frase_wald <- function(tab, decimal, idioma) {
  vapply(seq_len(nrow(tab)), function(i) {
    termo <- gsub(":", " × ", tab$termo[i], fixed = TRUE)
    paste0(termo, ", ", .ee_estat("Wald χ²", tab$wald_qui2[i], tab$gl[i], decimal), .ee_sep(decimal),
           .pr_p(tab$p[i], decimal))
  }, character(1))
}

# Só os pares entram na frase; as médias marginais por nível ficam fora dela (com nível
# de n < 5 seriam o escore de quem está no grupo). O rótulo do IC segue `nivel`.
.ee_mod_frase_pares <- function(ph, decimal, idioma, rot_dif, nivel = .95) {
  if (is.null(ph)) return(character())
  pares <- ph[grepl("^par", ph$tipo), , drop = FALSE]
  if (!nrow(pares)) return(character())
  por_termo <- split(pares, factor(pares$termo, levels = unique(pares$termo)))
  vapply(names(por_termo), function(tt) {
    x <- por_termo[[tt]]
    itens <- vapply(seq_len(nrow(x)), function(i) {
      cond <- if (nzchar(x$por[i])) paste0(" [", x$por[i], "]") else ""
      paste0(x$nivel_ou_par[i], cond, ": ", rot_dif, " = ", .pr_num(x$estimativa[i], 2, TRUE, decimal), .ee_sep(decimal),
             .ee_ic(x$ic_inf[i], x$ic_sup[i], 2, TRUE, decimal, nivel = nivel, idioma = idioma), .ee_sep(decimal), .pr_p(x$p[i], decimal))
    }, character(1))
    paste0(.ee_mod_t(idioma, "Comparações de ", "Comparisons for "), gsub(":", " × ", tt, fixed = TRUE),
           " (", sub("^par \\((.*)\\)$", "\\1", x$tipo[1]), "): ", paste(itens, collapse = "; "), ".")
  }, character(1))
}

# ================================================================ escolha da família
# Lê o tipo do desfecho e sugere distribuição e ligação (curso, T4 §3 a §6), sem ajustar
# modelo. A decisão final é teórica; AIC entre candidatas só na mesma escala de Y.
escolher_familia <- function(dados, vd) {
  dd <- .ee_dados(dados, vd)
  y <- dd$d[[vd]]
  if (is.factor(y) || is.character(y) || is.logical(y)) {
    k <- length(unique(y))
    if (k == 2) {
      sug <- data.frame(familia = "binomial", ligacao = "logit", leitura = "desfecho binário: OR = exp(B)")
    } else stop("desfecho categórico com ", k, " categorias: modelo_ordinal() se houver ordem, modelo_multinomial() se não houver")
    tab <- data.frame(n = length(y), categorias = k)
  } else {
    n <- length(y); inteiro <- all(abs(y - round(y)) < 1e-8)
    tab <- data.frame(n = n, minimo = min(y), maximo = max(y), prop_zeros = mean(y == 0),
                      prop_negativos = mean(y < 0), assimetria = .ee_assimetria(y),
                      razao_var_media = if (mean(y) > 0) stats::var(y) / mean(y) else NA_real_)
    if (all(y %in% c(0, 1))) {
      sug <- data.frame(familia = "binomial", ligacao = "logit", leitura = "0/1: regressão logística, OR = exp(B)")
    } else if (inteiro && min(y) >= 0) {
      sug <- data.frame(familia = c("poisson", "binomial_negativa"), ligacao = "log",
                        leitura = c("contagem: IRR = exp(B); checar superdispersão após o ajuste",
                                    if (isTRUE(tab$razao_var_media > 1.5)) "variância bem acima da média: binomial negativa provável" else
                                      "alternativa se houver superdispersão"))
    } else if (min(y) > 0) {
      sug <- data.frame(familia = c("normal", "gama", "gama"), ligacao = c("identity", "identity", "log"),
                        leitura = c("contínua: diferença média (B)",
                                    "positiva e assimétrica: diferença média na escala original (opção do curso)",
                                    "positiva e assimétrica: razão de médias = exp(B); mais estável com preditor contínuo"))
    } else {
      sug <- data.frame(familia = "normal", ligacao = "identity",
                        leitura = if (min(y) == 0 && max(y) > 0) "há zeros: gama impossível; normal, ou modelo em duas partes/Tweedie (fora deste módulo)" else
                          "valores negativos: só a normal entre as famílias deste módulo")
    }
  }
  avisos <- character()
  if (is.numeric(y) && min(y) == 0 && !all(y %in% c(0, 1)) && !all(abs(y - round(y)) < 1e-8)) {
    avisos <- c(avisos, "Há zeros num desfecho contínuo: somar uma constante para usar a gama muda o modelo (com ligação log os coeficientes mudam); não é neutro como sugere o curso. Preferir modelo em duas partes ou Tweedie.")
  }
  .ee_resultado("Escolha da distribuição e da ligação (GLM/GEE)", n = nrow(dd$d), n_excluidos = dd$n_excluidos,
                tabela = tab, efeito = NULL,
                decisoes = c("Sugestão pelo tipo de Y (curso, T4): contínua simétrica, normal-identidade; positiva assimétrica, gama (identidade para diferença média, log para razão de médias); 0/1, binomial-logit; contagem, Poisson-log, binomial negativa se superdispersa.",
                             "A escolha final é teórica; comparar_familias() confere o AIC entre candidatas na mesma escala de Y."),
                avisos = avisos, apa = character(), extra = list(sugestoes = sug))
}

# ================================================================ GLM
.ee_mod_checa_y <- function(y, familia) {
  if (familia == "binomial" && !all(y %in% c(0, 1))) stop("binomial: o desfecho precisa ser 0/1 (ou fator com 2 níveis)")
  if (familia %in% c("poisson", "binomial_negativa") && (any(y < 0) || any(abs(y - round(y)) > 1e-8)))
    stop(familia, ": o desfecho precisa ser contagem (inteiro >= 0)")
  if (familia == "gama" && any(y <= 0))
    stop("gama: o desfecho precisa ser > 0. Somar uma constante muda o modelo (ver references/modelos.md); considere modelo em duas partes ou Tweedie")
}

.ee_mod_y01 <- function(d, vd, idioma) {
  nota <- character()
  if (is.factor(d[[vd]]) || is.character(d[[vd]]) || is.logical(d[[vd]])) {
    f <- droplevels(as.factor(d[[vd]]))
    if (nlevels(f) != 2) stop("o desfecho binário precisa ter 2 categorias")
    d[[vd]] <- as.integer(f == levels(f)[2])
    nota <- paste0("Desfecho binário codificado como 1 = '", levels(f)[2], "' (segundo nível do fator).")
  }
  list(d = d, nota = nota)
}

# Ajuste do GLM. Devolve o modelo com contrastes soma-zero (Tipo III) e o modelo com
# contrastes de tratamento (coeficientes contra a referência); os dois têm o mesmo ajuste.
.ee_mod_glm_fit <- function(formula, d, familia, ligacao) {
  ajusta <- function() {
    if (familia == "binomial_negativa") return(MASS::glm.nb(formula, data = d, link = log))
    fam <- switch(familia, normal = stats::gaussian(ligacao), gama = stats::Gamma(ligacao),
                  binomial = stats::binomial(ligacao), poisson = stats::poisson(ligacao))
    st <- NULL
    if (familia == "gama" && ligacao == "identity") {
      st <- tryCatch(stats::coef(stats::lm(formula, data = d)), error = function(e) NULL)
    }
    suppressWarnings(stats::glm(formula, family = fam, data = d, start = st, control = stats::glm.control(maxit = 100)))
  }
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op), add = TRUE)
  fs <- ajusta()
  options(contrasts = c("contr.treatment", "contr.poly"))
  ft <- ajusta()
  list(soma = fs, tratamento = ft)
}

# Escala (dispersão) do GLM. "ml": máxima verossimilhança, padrão do SPSS (GENLIN) e do
# curso; "pearson": χ² de Pearson / (n - p), padrão do glm() do R. Só normal e gama têm
# dispersão livre; Poisson, binomial e binomial negativa ficam com 1.
.ee_mod_glm_escala <- function(fit, familia, escala) {
  if (!familia %in% c("normal", "gama")) return(1)
  if (escala == "pearson") return(sum(stats::residuals(fit, "pearson")^2) / fit$df.residual)
  if (familia == "normal") return(sum(stats::residuals(fit, "response")^2) / length(fit$y))
  1 / MASS::gamma.shape(fit)$alpha
}

# Log-verossimilhança completa. Na gama, com a forma estimada por ML (o AIC do glm() do R
# usa outra aproximação da dispersão e não bate com o SPSS).
.ee_mod_glm_loglik <- function(fit, familia) {
  if (familia == "gama") {
    a <- MASS::gamma.shape(fit)$alpha
    ll <- sum(stats::dgamma(fit$y, shape = a, scale = stats::fitted(fit) / a, log = TRUE))
    return(c(ll = ll, k = length(stats::coef(fit)) + 1))
  }
  l <- stats::logLik(fit)
  c(ll = as.numeric(l), k = attr(l, "df"))
}

modelo_glm <- function(dados, vd, preditores, familia = c("normal", "gama", "binomial", "poisson", "binomial_negativa"),
                       ligacao = NULL, referencia = NULL, fatores = NULL, escala = c("ml", "pearson"), exposicao = NULL,
                       comparar = NULL, ajuste = "bonferroni", nivel = .95, decimal = ".", idioma = "pt") {
  familia <- match.arg(familia); escala <- match.arg(escala)
  .pr_exige("MASS", "emmeans")
  vars <- .ee_mod_vars(preditores)
  .ee_mod_nome_ok(c(vd, vars, exposicao))
  dd <- .ee_dados(dados, unique(c(vd, vars, exposicao)))
  d <- .ee_mod_fatores(dd$d, vars, referencia, fatores)
  decisoes <- character(); avisos <- .ee_mod_aviso_numericos(d, vars)
  if (familia == "binomial") { b01 <- .ee_mod_y01(d, vd, idioma); d <- b01$d; decisoes <- c(decisoes, b01$nota) }
  if (!is.numeric(d[[vd]])) stop("o desfecho precisa ser numérico (ou binário)")
  .ee_mod_checa_y(d[[vd]], familia)
  continuos <- vars[vapply(vars, function(v) is.numeric(d[[v]]), logical(1))]
  if (is.null(ligacao)) {
    ligacao <- switch(familia, normal = "identity", binomial = "logit", poisson = "log", binomial_negativa = "log",
                      gama = if (length(continuos)) "log" else "identity")
    if (familia == "gama") decisoes <- c(decisoes, if (length(continuos))
      "Gama com ligação log (há preditor contínuo: a identidade pode gerar médias negativas e não convergir); exp(B) = razão de médias." else
      "Gama com ligação identidade (opção do curso para ler diferença média; com só fatores, identidade e log dão as mesmas médias ajustadas).")
  }
  if (familia == "binomial_negativa" && ligacao != "log") stop("binomial negativa: só a ligação log está implementada")
  if (!is.null(exposicao)) {
    if (!familia %in% c("poisson", "binomial_negativa")) stop("exposicao só se aplica a contagens (Poisson ou binomial negativa)")
    if (any(d[[exposicao]] <= 0)) stop("a exposição precisa ser > 0")
    extra_f <- paste0("offset(log(", exposicao, "))")
    decisoes <- c(decisoes, paste0("Janela de observação variável: offset(log(", exposicao, ")); exp(B) é razão de taxas."))
  } else extra_f <- NULL
  f <- .ee_mod_formula(vd, preditores, extra_f)
  fits <- .ee_mod_glm_fit(f, d, familia, ligacao)
  fs <- fits$soma; ft <- fits$tratamento
  if (!isTRUE(fs$converged)) avisos <- c(avisos, "O modelo não convergiu: resultados não confiáveis (com gama-identidade, tentar ligação log).")
  if (familia == "gama" && ligacao == "identity" && any(stats::fitted(fs) <= 0))
    avisos <- c(avisos, "Médias ajustadas <= 0 com gama-identidade: usar ligação log.")
  phi <- .ee_mod_glm_escala(fs, familia, escala)
  Vs <- if (familia %in% c("normal", "gama")) summary(fs, dispersion = phi)$cov.scaled else stats::vcov(fs)
  Vt <- if (familia %in% c("normal", "gama")) summary(ft, dispersion = phi)$cov.scaled else stats::vcov(ft)
  mm <- stats::model.matrix(fs)
  rotulos <- attr(stats::terms(fs), "term.labels")
  wald <- .ee_mod_wald3(stats::coef(fs), Vs, attr(mm, "assign"), rotulos)
  n <- nrow(d); p <- sum(!is.na(stats::coef(fs)))
  if (familia == "normal") {
    # F clássico (escala n - p): a ANOVA/regressão equivalente, como sensibilidade.
    fator <- if (escala == "ml") (n - p) / n else 1
    wald$F <- wald$wald_qui2 * fator / wald$gl
    wald$gl_res <- n - p
    wald$p_F <- stats::pf(wald$F, wald$gl, n - p, lower.tail = FALSE)
  }
  rot_exp <- .ee_mod_rot_exp(familia, ligacao, idioma)
  coefs <- .ee_mod_coefs(stats::coef(ft), Vt, nivel, rot_exp)
  ll <- .ee_mod_glm_loglik(fs, familia)
  aic <- -2 * ll[["ll"]] + 2 * ll[["k"]]
  # efeito: B (identidade) ou exp(B) (log/logit) de cada coeficiente fora o intercepto
  sem_int <- coefs[coefs$termo != "(Intercept)", , drop = FALSE]
  efeito <- if (is.na(rot_exp)) {
    data.frame(indice = paste0("B: ", sem_int$termo), valor = sem_int$B, ic_inf = sem_int$ic_inf,
               ic_sup = sem_int$ic_sup, principal = rep(TRUE, nrow(sem_int)))
  } else {
    data.frame(indice = paste0(rot_exp, ": ", sem_int$termo), valor = sem_int$exp_B, ic_inf = sem_int$exp_ic_inf,
               ic_sup = sem_int$exp_ic_sup, principal = rep(TRUE, nrow(sem_int)))
  }
  # pós-hoc
  termos_ph <- if (is.null(comparar)) .ee_mod_termos_fator(rotulos, d) else comparar
  ph <- if (length(termos_ph)) .ee_mod_emm(fs, termos_ph, Vs, d, ajuste, nivel) else NULL
  avisos <- c(avisos, .ee_mod_aviso_pequenos(ph, "casos"))
  # pressupostos
  pres <- data.frame(checagem = "Convergência", estatistica = fs$iter, gl = NA, p = NA,
                     leitura = if (isTRUE(fs$converged)) "convergiu (número de iterações)" else "NÃO convergiu")
  if (familia == "normal") {
    r <- stats::residuals(fs, "response")
    if (length(r) >= 3 && length(r) <= 5000) {
      sw <- stats::shapiro.test(r)
      pres <- rbind(pres, data.frame(checagem = "Shapiro-Wilk dos resíduos", estatistica = unname(sw$statistic), gl = NA,
                                     p = sw$p.value, leitura = "indício entre vários; com n grande acusa desvios irrelevantes"))
    }
  }
  if (familia == "poisson") {
    x2 <- sum(stats::residuals(fs, "pearson")^2); gl <- fs$df.residual
    pdisp <- stats::pchisq(x2, gl, lower.tail = FALSE)
    nb <- tryCatch(suppressWarnings(MASS::glm.nb(f, data = d)), error = function(e) NULL)
    aic_nb <- if (!is.null(nb)) stats::AIC(nb) else NA_real_
    pres <- rbind(pres, data.frame(checagem = "Dispersão de Pearson (χ²/gl)", estatistica = x2 / gl, gl = gl, p = pdisp,
                                   leitura = if (pdisp < .05) "superdispersão: os EPs da Poisson ficam pequenos demais" else "compatível com a Poisson (razão perto de 1)"))
    pres <- rbind(pres, data.frame(checagem = "AIC Poisson vs binomial negativa", estatistica = aic - aic_nb, gl = NA, p = NA,
                                   leitura = paste0("AIC Poisson = ", round(aic, 2), "; binomial negativa = ", round(aic_nb, 2),
                                                    " (diferença positiva favorece a binomial negativa)")))
    if (pdisp < .05) avisos <- c(avisos, paste0("Superdispersão (χ²/gl = ", round(x2 / gl, 2), ", p ", sub("^p ", "", .pr_p(pdisp)),
                                                "): refazer com familia = \"binomial_negativa\"."))
  }
  if (familia == "binomial_negativa") {
    pres <- rbind(pres, data.frame(checagem = "θ da binomial negativa", estatistica = fs$theta, gl = NA, p = NA,
                                   leitura = paste0("variância = μ + μ²/θ; EP de θ = ", round(fs$SE.theta, 3))))
  }
  # textos
  nome <- paste0("Generalized Linear Model (", .ee_mod_t(idioma, "distribuição ", "distribution "),
                 .ee_mod_rot_familia(familia, idioma), ", ", .ee_mod_t(idioma, "ligação ", "link "),
                 .ee_mod_rot_ligacao(ligacao, idioma), ")")
  apa <- paste0(nome, "; AIC = ", .pr_num(aic, 2, TRUE, decimal), ": ",
                paste(.ee_mod_frase_wald(wald, decimal, idioma), collapse = "; "), ".")
  rot_dif <- if (is.na(rot_exp)) .ee_mod_t(idioma, "diferença", "difference") else
    if (ligacao == "logit") "OR" else .ee_mod_t(idioma, "razão", "ratio")
  apa <- c(apa, .ee_mod_frase_pares(ph, decimal, idioma, rot_dif, nivel))
  decisoes <- c(decisoes,
                "Testes dos termos: Wald χ² Tipo III com contrastes soma-zero (\"Tests of Model Effects\" do SPSS). O χ² de Wald tem um só gl; a forma \"Wald (2, 23)\" das notas do curso não é usada.",
                if (familia %in% c("normal", "gama")) paste0("Escala (dispersão) por ", if (escala == "ml") "máxima verossimilhança, padrão do SPSS e do curso; com n pequeno o Wald χ² fica liberal (escala = \"pearson\" usa n - p, padrão do R)." else "χ² de Pearson / (n - p), padrão do glm() do R."),
                if (familia == "normal") "Normal-identidade é o mesmo modelo da ANOVA/regressão: o F clássico (escala n - p) sai na tabela como sensibilidade.",
                if (any(vapply(vars, function(v) is.factor(d[[v]]), logical(1)))) "Coeficientes de fatores contra a categoria de referência (primeiro nível; o SPSS usa o último, ver argumento referencia), IC de Wald." else "Coeficientes com IC de Wald; preditores numéricos entram como contínuos (um coeficiente por unidade).",
                if (length(termos_ph)) paste0("Médias marginais estimadas (emmeans, escala da resposta) e pares com correção de ", ajuste, "; IC dos pares ajustados pela mesma correção."))
  .ee_mod_classe(.ee_resultado("Modelo linear generalizado (GLM)", n = n, n_excluidos = dd$n_excluidos, tabela = wald,
                efeito = efeito, pressupostos = pres, posthoc = ph, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(modelo = ft, modelo_soma = fs, coeficientes = coefs, aic = aic, loglik = ll[["ll"]],
                             k_parametros = ll[["k"]], escala_phi = phi, familia = familia, ligacao = ligacao)))
}

# AIC entre distribuições candidatas para o MESMO Y, mesmos casos e mesmos preditores.
comparar_familias <- function(dados, vd, preditores, familias = c("normal", "gama"), ligacoes = NULL,
                              fatores = NULL, decimal = ".", idioma = "pt") {
  .pr_exige("MASS")
  vars <- .ee_mod_vars(preditores)
  .ee_mod_nome_ok(c(vd, vars))
  dd <- .ee_dados(dados, unique(c(vd, vars)))
  d <- .ee_mod_fatores(dd$d, vars, fatores = fatores)
  if (is.null(ligacoes)) ligacoes <- vapply(familias, function(f) switch(f, normal = "identity", gama = "identity",
                                                                         binomial = "logit", "log"), character(1))
  if (length(ligacoes) != length(familias)) stop("ligacoes precisa ter o mesmo comprimento de familias")
  f <- .ee_mod_formula(vd, preditores)
  linhas <- list(); avisos <- .ee_mod_aviso_numericos(d, vars)
  for (i in seq_along(familias)) {
    fam <- familias[i]; lig <- ligacoes[i]
    ok <- tryCatch({ .ee_mod_checa_y(d[[vd]], fam); TRUE }, error = function(e) { avisos <<- c(avisos, conditionMessage(e)); FALSE })
    if (!ok) next
    fit <- tryCatch(.ee_mod_glm_fit(f, d, fam, lig)$soma, error = function(e) NULL)
    if (is.null(fit)) { avisos <- c(avisos, paste0(fam, "-", lig, ": falhou no ajuste")); next }
    ll <- .ee_mod_glm_loglik(fit, fam)
    linhas[[length(linhas) + 1]] <- data.frame(familia = fam, ligacao = lig, loglik = ll[["ll"]], k = ll[["k"]],
                                               AIC = -2 * ll[["ll"]] + 2 * ll[["k"]], convergiu = isTRUE(fit$converged))
  }
  if (!length(linhas)) stop("nenhuma família pôde ser ajustada")
  tab <- do.call(rbind, linhas)
  tab$delta_AIC <- tab$AIC - min(tab$AIC)
  tab$peso_akaike <- exp(-tab$delta_AIC / 2) / sum(exp(-tab$delta_AIC / 2))
  cont <- c("normal", "gama"); disc <- c("poisson", "binomial_negativa", "binomial")
  if (any(tab$familia %in% cont) && any(tab$familia %in% disc))
    avisos <- c(avisos, "AIC entre distribuição contínua (densidade) e discreta (probabilidade) não é comparável: compare só dentro de cada tipo.")
  melhor <- tab[which.min(tab$AIC), ]
  apa <- paste0(.ee_mod_t(idioma, "Menor AIC: ", "Lowest AIC: "), .ee_mod_rot_familia(melhor$familia, idioma), "-",
                .ee_mod_rot_ligacao(melhor$ligacao, idioma), " (AIC = ", .pr_num(melhor$AIC, 2, TRUE, decimal), ")",
                if (nrow(tab) > 1) paste0(.ee_mod_t(idioma, "; demais: ", "; others: "),
                                          paste(paste0(vapply(tab$familia[-which.min(tab$AIC)], .ee_mod_rot_familia, character(1), idioma = idioma), "-",
                                                       vapply(tab$ligacao[-which.min(tab$AIC)], .ee_mod_rot_ligacao, character(1), idioma = idioma),
                                                       " = ", vapply(tab$AIC[-which.min(tab$AIC)], .pr_num, character(1), casas = 2, zero = TRUE, decimal = decimal)),
                                                collapse = .ee_sep(decimal))) else "", ".")
  .ee_resultado("Comparação de distribuições por AIC (GLM)", n = nrow(d), n_excluidos = dd$n_excluidos, tabela = tab,
                decisoes = c("AIC = -2 logLik + 2k com verossimilhança completa (na gama, forma estimada por ML, como no SPSS); menor é melhor, sem ponto de corte; peso de Akaike ao lado.",
                             "Só vale para o mesmo Y (sem transformação), os mesmos casos e os mesmos preditores. No GEE não há verossimilhança: o QIC entre famílias diferentes compara quasi-verossimilhanças em escalas diferentes (o caso QIC 31,53 vs 16308,22 das notas do curso é artefato de escala) e não é usado para escolher a distribuição.",
                             "A distribuição se escolhe primeiro pela teoria e pelo tipo de Y; o AIC confirma entre candidatas plausíveis."),
                avisos = avisos, apa = apa)
}

# ================================================================ GEE
.ee_mod_gee_prepara <- function(dados, vd, preditores, id, tempo, referencia, familia, idioma, fatores = NULL) {
  vars <- .ee_mod_vars(preditores)
  .ee_mod_nome_ok(c(vd, vars, id, tempo))
  dd <- .ee_dados(dados, unique(c(vd, vars, id, tempo)))
  d <- dd$d; notas <- character()
  if (is.character(d[[tempo]])) {
    notas <- c(notas, "tempo veio como texto e foi ordenado alfabeticamente; passe um fator com os níveis na ordem cronológica.")
    d[[tempo]] <- factor(d[[tempo]])
  }
  if (tempo %in% vars && !is.factor(d[[tempo]])) d[[tempo]] <- factor(d[[tempo]])
  d <- .ee_mod_fatores(d, setdiff(vars, tempo), referencia, setdiff(fatores, tempo))
  notas <- c(notas, .ee_mod_aviso_numericos(d, vars, isentos = tempo))
  if (is.factor(d[[tempo]])) d[[tempo]] <- droplevels(d[[tempo]])
  if (familia == "binomial") { b01 <- .ee_mod_y01(d, vd, idioma); d <- b01$d; notas <- c(notas, b01$nota) }
  .ee_mod_checa_y(d[[vd]], familia)
  d <- d[order(d[[id]], d[[tempo]]), , drop = FALSE]
  d$.id <- as.integer(factor(as.character(d[[id]]), levels = unique(as.character(d[[id]]))))
  d$.onda <- as.integer(if (is.factor(d[[tempo]])) d[[tempo]] else factor(d[[tempo]]))
  if (anyDuplicated(d[, c(".id", ".onda")])) stop("há mais de uma linha por sujeito no mesmo tempo: o GEE pede formato longo com uma linha por sujeito e tempo")
  rownames(d) <- NULL
  list(d = d, n_excluidos = dd$n_excluidos, notas = notas, vars = vars)
}

.ee_mod_gee_familia <- function(familia, ligacao) {
  switch(familia, normal = stats::gaussian(ligacao), gama = stats::Gamma(ligacao),
         binomial = stats::binomial(ligacao), poisson = stats::poisson(ligacao),
         binomial_negativa = stop("binomial negativa não é nativa no geepack (só com θ fixo); use Poisson com EP robusto ou um GLMM"))
}

# geeglm por estrutura. A não estruturada usa corstr = "userdefined" com zcor de
# geepack::genZcor e SEM o argumento waves: no geepack 1.3.13, "unstructured" com tempos
# faltando no meio derruba a sessão do R (segfault) e, sem waves, pareia as correlações
# pela posição da linha, não pelo tempo. A AR(1) usa waves = índice do tempo.
.ee_mod_gee_fit <- function(formula, d, fam, corstr) {
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op), add = TRUE)
  ajusta <- function(contr) {
    options(contrasts = contr)
    if (corstr == "unstructured") {
      z <- geepack::genZcor(clusz = rle(d$.id)$lengths, waves = d$.onda, corstrv = 4)
      geepack::geeglm(formula, family = fam, data = d, id = .id, corstr = "userdefined", zcor = z)
    } else if (corstr == "ar1") {
      geepack::geeglm(formula, family = fam, data = d, id = .id, waves = .onda, corstr = "ar1")
    } else geepack::geeglm(formula, family = fam, data = d, id = .id, corstr = corstr)
  }
  fs <- ajusta(c("contr.sum", "contr.poly"))
  ft <- ajusta(c("contr.treatment", "contr.poly"))
  list(soma = fs, tratamento = ft)
}

# Matriz de correlação de trabalho completa (T x T) a partir dos alfas do geeglm.
.ee_mod_gee_R <- function(fit, corstr, T) {
  a <- fit$geese$alpha
  switch(corstr,
         independence = diag(T),
         exchangeable = { R <- matrix(a[1], T, T); diag(R) <- 1; R },
         ar1 = a[1]^abs(outer(seq_len(T), seq_len(T), "-")),
         unstructured = {
           R <- diag(T)
           pares <- do.call(rbind, strsplit(sub("^alpha\\.", "", names(a)), ":", fixed = TRUE))
           for (k in seq_along(a)) { i <- as.integer(pares[k, 1]); j <- as.integer(pares[k, 2]); R[i, j] <- R[j, i] <- a[k] }
           R
         })
}

# Covariância sanduíche (robusta) e com correção de viés de Mancl e DeRouen (2001):
# (I - H_ii)^-1 aplicado aos resíduos de cada cluster, H_ii = D_i B^-1 D_i' V_i^-1.
# A escala φ se cancela nas duas. O geepack não traz a correção.
.ee_mod_gee_vcov <- function(fit, d, corstr) {
  X <- stats::model.matrix(fit)
  ok <- !is.na(stats::coef(fit)); X <- X[, ok, drop = FALSE]
  eta <- fit$linear.predictors; mu <- fit$fitted.values; y <- fit$y
  fam <- fit$family
  dmu <- fam$mu.eta(eta); va <- fam$variance(mu)
  T <- max(d$.onda); R <- .ee_mod_gee_R(fit, corstr, T)
  blocos <- split(seq_len(nrow(d)), d$.id)
  p <- ncol(X); B <- matrix(0, p, p); partes <- vector("list", length(blocos))
  for (b in seq_along(blocos)) {
    i <- blocos[[b]]
    Di <- dmu[i] * X[i, , drop = FALSE]
    a <- sqrt(va[i]); Ri <- R[d$.onda[i], d$.onda[i], drop = FALSE]
    Vinv <- solve(outer(a, a) * Ri)
    B <- B + crossprod(Di, Vinv %*% Di)
    partes[[b]] <- list(Di = Di, Vinv = Vinv, r = y[i] - mu[i])
  }
  Binv <- solve(B)
  M0 <- matrix(0, p, p); M1 <- matrix(0, p, p); md_ok <- TRUE
  for (pt in partes) {
    u0 <- crossprod(pt$Di, pt$Vinv %*% pt$r)
    H <- pt$Di %*% Binv %*% crossprod(pt$Di, pt$Vinv)
    # Cluster que determina sozinho um parâmetro (por exemplo, nível de fator com um só
    # sujeito) tem alavanca 1: I - H é singular e a correção não se define.
    IH <- diag(nrow(H)) - H
    s1 <- if (md_ok && rcond(IH) > 1e-10) tryCatch(solve(IH, pt$r), error = function(e) NULL) else NULL
    if (is.null(s1)) md_ok <- FALSE else M1 <- M1 + tcrossprod(crossprod(pt$Di, pt$Vinv %*% s1))
    M0 <- M0 + tcrossprod(u0)
  }
  nm <- colnames(X)
  V0 <- Binv %*% M0 %*% Binv
  V1 <- if (md_ok) Binv %*% M1 %*% Binv else NULL
  dimnames(V0) <- list(nm, nm); if (md_ok) dimnames(V1) <- list(nm, nm)
  list(robusto = V0, mancl_derouen = V1)
}

# QIC de Pan (2001) com a quasi-verossimilhança dividida pela escala do modelo de
# independência (a mesma usada em Ω_I). O geepack::QIC não divide pela escala (na normal
# o termo vira soma de quadrados em unidades de Y²) e, na versão 1.3.13, escreve a
# quasi-verossimilhança da gama como -y/(μ - log μ) em vez de -y/μ - log μ.
.ee_mod_qic <- function(fit, fit_ind) {
  y <- fit$y; mu <- fit$fitted.values
  q <- switch(fit$family$family,
              gaussian = -(y - mu)^2 / 2,
              poisson = y * log(mu) - mu,
              binomial = ifelse(y == 1, log(mu), log(1 - mu)),
              Gamma = -y / mu - log(mu),
              stop("família sem quasi-verossimilhança implementada"))
  phi <- if (isTRUE(fit_ind$geese$model$scale.fix)) 1 else as.numeric(fit_ind$geese$gamma)[1]
  Q <- sum(q) / phi
  Om <- tryCatch(solve(fit_ind$geese$vbeta.naiv), error = function(e) MASS::ginv(fit_ind$geese$vbeta.naiv))
  tr <- as.numeric(sum(diag(Om %*% fit$geese$vbeta)))
  pp <- sum(!is.na(stats::coef(fit)))
  c(QIC = -2 * Q + 2 * tr, QICu = -2 * Q + 2 * pp, CIC = tr, quasi_verossimilhanca = Q)
}

.ee_mod_rot_corr <- function(cs, idioma) {
  switch(cs, independence = .ee_mod_t(idioma, "independente", "independence"),
         exchangeable = .ee_mod_t(idioma, "intercambiável", "exchangeable"),
         ar1 = "AR(1)", unstructured = .ee_mod_t(idioma, "não estruturada", "unstructured"))
}

modelo_gee <- function(dados, vd, preditores, id, tempo, familia = c("normal", "gama", "binomial", "poisson"),
                       ligacao = NULL, correlacao = c("unstructured", "exchangeable", "ar1", "independence"),
                       ep = c("robusto", "mancl_derouen"), referencia = NULL, fatores = NULL, comparar = NULL,
                       ajuste = "bonferroni", nivel = .95, decimal = ".", idioma = "pt") {
  familia <- match.arg(familia); correlacao <- match.arg(correlacao); ep <- match.arg(ep)
  .pr_exige("geepack", "emmeans")
  prep <- .ee_mod_gee_prepara(dados, vd, preditores, id, tempo, referencia, familia, idioma, fatores)
  d <- prep$d; vars <- prep$vars
  num <- grepl("e entra como contínua", prep$notas, fixed = TRUE)
  decisoes <- prep$notas[!num]; avisos <- prep$notas[num]
  if (is.null(ligacao)) {
    continuos <- vars[vapply(vars, function(v) is.numeric(d[[v]]), logical(1))]
    ligacao <- switch(familia, normal = "identity", binomial = "logit", poisson = "log",
                      gama = if (length(continuos)) "log" else "identity")
  }
  fam <- .ee_mod_gee_familia(familia, ligacao)
  f <- .ee_mod_formula(vd, preditores)
  K <- max(d$.id); T <- max(d$.onda); nobs <- nrow(d)
  por_sujeito <- tabulate(d$.id)
  # ajuste com retorno pela ordem do curso se a estrutura pedida falhar
  tentativa <- unique(c(correlacao, if (correlacao == "unstructured") c("ar1", "exchangeable")))
  fits <- NULL; usada <- NA_character_
  for (cs in tentativa) {
    fits <- tryCatch(.ee_mod_gee_fit(f, d, fam, cs), error = function(e) NULL)
    if (!is.null(fits) && isTRUE(fits$soma$geese$error == 0) && all(is.finite(fits$soma$geese$alpha))) {
      # a matriz de correlação de trabalho precisa ser positiva definida
      Rw <- .ee_mod_gee_R(fits$soma, cs, T)
      if (min(eigen(Rw, symmetric = TRUE, only.values = TRUE)$values) > 1e-8) { usada <- cs; break }
    }
    fits <- NULL
  }
  if (is.null(fits)) stop("o GEE não convergiu com nenhuma estrutura tentada (", paste(tentativa, collapse = ", "), ")")
  if (usada != correlacao) avisos <- c(avisos, paste0("A estrutura ", .ee_mod_rot_corr(correlacao, idioma), " não convergiu; usada ", .ee_mod_rot_corr(usada, idioma), "."))
  fs <- fits$soma; ft <- fits$tratamento
  Vs <- .ee_mod_gee_vcov(fs, d, usada); Vt <- .ee_mod_gee_vcov(ft, d, usada)
  dif <- max(abs(Vs$robusto - fs$geese$vbeta[!is.na(stats::coef(fs)), !is.na(stats::coef(fs))]))
  if (dif > 1e-6 * max(abs(fs$geese$vbeta))) avisos <- c(avisos, "A covariância robusta recalculada difere da do geepack; conferir.")
  md_ok <- !is.null(Vs$mancl_derouen) && !is.null(Vt$mancl_derouen)
  if (!md_ok) {
    avisos <- c(avisos, .ee_mod_t(idioma, "A correção de Mancl e DeRouen não se define aqui: algum sujeito determina sozinho um parâmetro (alavanca 1, como num nível de fator com um só sujeito). Saem só o EP robusto e o p dele; agrupar a categoria pequena.",
                                 "The Mancl-DeRouen correction is undefined here: some subject alone determines a parameter (leverage 1, as in a factor level with a single subject). Only the robust SE and its p are reported; merge the small category."))
    if (ep == "mancl_derouen") ep <- "robusto"
  }
  V_s <- Vs[[ep]]; V_t <- Vt[[ep]]
  mm <- stats::model.matrix(fs); rotulos <- attr(stats::terms(fs), "term.labels")
  wald <- .ee_mod_wald3(stats::coef(fs), V_s, attr(mm, "assign"), rotulos)
  poucos <- K < 40
  if (md_ok && (poucos || ep == "mancl_derouen")) {
    outro <- if (ep == "robusto") "mancl_derouen" else "robusto"
    w2 <- .ee_mod_wald3(stats::coef(fs), Vs[[outro]], attr(mm, "assign"), rotulos)
    wald[[paste0("wald_qui2_", outro)]] <- w2$wald_qui2; wald[[paste0("p_", outro)]] <- w2$p
  }
  rot_exp <- .ee_mod_rot_exp(familia, ligacao, idioma)
  coefs <- .ee_mod_coefs(stats::coef(ft), V_t, nivel, rot_exp)
  if (md_ok && (poucos || ep == "mancl_derouen")) {
    outro <- if (ep == "robusto") "mancl_derouen" else "robusto"
    c2 <- .ee_mod_coefs(stats::coef(ft), Vt[[outro]], nivel)
    coefs[[paste0("ep_", outro)]] <- c2$ep; coefs[[paste0("p_", outro)]] <- c2$p
  }
  sem_int <- coefs[coefs$termo != "(Intercept)", , drop = FALSE]
  efeito <- if (is.na(rot_exp)) {
    data.frame(indice = paste0("B: ", sem_int$termo), valor = sem_int$B, ic_inf = sem_int$ic_inf,
               ic_sup = sem_int$ic_sup, principal = rep(TRUE, nrow(sem_int)))
  } else {
    data.frame(indice = paste0(rot_exp, ": ", sem_int$termo), valor = sem_int$exp_B, ic_inf = sem_int$exp_ic_inf,
               ic_sup = sem_int$exp_ic_sup, principal = rep(TRUE, nrow(sem_int)))
  }
  termos_ph <- if (is.null(comparar)) .ee_mod_termos_fator(rotulos, d) else comparar
  ph <- if (length(termos_ph)) .ee_mod_emm(fs, termos_ph, V_s, d, ajuste, nivel, id = ".id") else NULL
  avisos <- c(avisos, .ee_mod_aviso_pequenos(ph, "sujeitos"))
  # qualidade do ajuste: QIC da estrutura usada (Pan, 2001), mesma família
  fi <- if (usada == "independence") fs else tryCatch(.ee_mod_gee_fit(f, d, fam, "independence")$soma, error = function(e) NULL)
  qic <- if (!is.null(fi)) .ee_mod_qic(fs, fi) else rep(NA_real_, 4)
  alfa <- fs$geese$alpha
  pres <- data.frame(checagem = c("Sujeitos (clusters)", "Observações por sujeito", "Tempos (ondas)", "Correlação de trabalho estimada"),
                     estatistica = c(K, mean(por_sujeito), T, if (length(alfa)) mean(alfa) else 0),
                     gl = NA, p = NA,
                     leitura = c(if (poucos && md_ok) "menos de 40 clusters: EP sanduíche viesado para baixo; correção de Mancl e DeRouen ao lado" else if (poucos) "menos de 40 clusters: EP sanduíche viesado para baixo; correção de Mancl e DeRouen indefinida (ver avisos)" else "clusters suficientes para o EP sanduíche",
                                 paste0("mín. ", min(por_sujeito), ", máx. ", max(por_sujeito), "; ", nobs, " observações; ", prep$n_excluidos, " linhas excluídas por ausência"),
                                 if (T == 2) "com 2 tempos as estruturas não estruturada, intercambiável e AR(1) são a mesma (uma só correlação)" else "estrutura importa com 3 ou mais tempos",
                                 if (length(alfa)) paste0("alfa: ", paste(round(alfa, 3), collapse = ", ")) else "independência (alfa = 0)"))
  if (poucos && md_ok) avisos <- c(avisos, paste0("Só ", K, " sujeitos: o EP sanduíche infla o erro tipo I com poucos clusters (Mancl e DeRouen, 2001). A tabela traz o EP e o p com a correção de Mancl e DeRouen ao lado; com menos de cerca de 40 clusters, preferir ep = \"mancl_derouen\"."))
  if (usada == "ar1" && is.numeric(dados[[tempo]])) {
    tt <- sort(unique(dados[[tempo]][!is.na(dados[[tempo]])]))
    if (length(tt) > 2 && length(unique(round(diff(tt), 8))) > 1)
      avisos <- c(avisos, "AR(1) com intervalos desiguais entre os tempos: o geepack usa a ordem dos tempos (1, 2, 3...), não a distância real.")
  }
  if (prep$n_excluidos > 0) avisos <- c(avisos, paste0(prep$n_excluidos, " linhas com ausência saíram (o sujeito fica com os tempos que tem). O GEE padrão só é válido sob ausência completamente aleatória (MCAR); sob MAR, usar GEE ponderado (IPW) ou imputação múltipla."))
  nome <- paste0("GEE (", .ee_mod_t(idioma, "distribuição ", "distribution "), .ee_mod_rot_familia(familia, idioma), ", ",
                 .ee_mod_t(idioma, "ligação ", "link "), .ee_mod_rot_ligacao(ligacao, idioma), ", ",
                 .ee_mod_t(idioma, "correlação de trabalho ", "working correlation "), .ee_mod_rot_corr(usada, idioma), ", ",
                 if (ep == "robusto") .ee_mod_t(idioma, "EP robusto", "robust SE") else .ee_mod_t(idioma, "EP com correção de Mancl e DeRouen", "Mancl-DeRouen corrected SE"),
                 "; ", K, .ee_mod_t(idioma, " sujeitos, ", " subjects, "), nobs, .ee_mod_t(idioma, " observações)", " observations)"))
  apa <- paste0(nome, ": ", paste(.ee_mod_frase_wald(wald, decimal, idioma), collapse = "; "), ".")
  rot_dif <- if (is.na(rot_exp)) .ee_mod_t(idioma, "diferença", "difference") else if (ligacao == "logit") "OR" else .ee_mod_t(idioma, "razão", "ratio")
  apa <- c(apa, .ee_mod_frase_pares(ph, decimal, idioma, rot_dif, nivel))
  decisoes <- c(decisoes,
                paste0("Formato longo ordenado por sujeito e tempo; cada linha ausente sai sozinha (", nobs, " observações de ", K, " sujeitos)."),
                paste0("Correlação de trabalho ", .ee_mod_rot_corr(usada, idioma), if (correlacao == "unstructured") " (padrão do curso sem hipótese forte sobre o padrão temporal; com falha, AR(1) e depois intercambiável)." else "."),
                "EP sanduíche (robusto) por padrão: as estimativas são consistentes mesmo com a correlação de trabalho errada; a estrutura afeta sobretudo a eficiência.",
                "Testes dos termos: Wald χ² Tipo III (contrastes soma-zero), um gl por coluna do termo; a forma \"Wald (2, 24)\" das notas do curso não é usada.",
                "Interação sempre tabelada e relatada; pós-hoc da interação só se interpreta como confirmatório se o termo for significativo.",
                "QIC desta estrutura em extra$qic; para escolher a estrutura use comparar_correlacoes_gee() (mesma família). QIC entre famílias não é usado.")
  .ee_mod_classe(.ee_resultado("Equações de estimação generalizadas (GEE)", n = K, n_excluidos = prep$n_excluidos, tabela = wald,
                efeito = efeito, pressupostos = pres, posthoc = ph, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(modelo = ft, modelo_soma = fs, coeficientes = coefs, vcov = Vt, qic = qic,
                             correlacao = usada, n_observacoes = nobs, familia = familia, ligacao = ligacao)))
}

# QIC, QICu e CIC por estrutura de correlação, mesma família, ligação e preditores.
comparar_correlacoes_gee <- function(dados, vd, preditores, id, tempo, familia = c("normal", "gama", "binomial", "poisson"),
                                     ligacao = NULL, correlacoes = c("independence", "exchangeable", "ar1", "unstructured"),
                                     fatores = NULL, decimal = ".", idioma = "pt") {
  familia <- match.arg(familia)
  .pr_exige("geepack")
  prep <- .ee_mod_gee_prepara(dados, vd, preditores, id, tempo, NULL, familia, idioma, fatores)
  d <- prep$d; vars <- prep$vars
  if (is.null(ligacao)) {
    continuos <- vars[vapply(vars, function(v) is.numeric(d[[v]]), logical(1))]
    ligacao <- switch(familia, normal = "identity", binomial = "logit", poisson = "log",
                      gama = if (length(continuos)) "log" else "identity")
  }
  fam <- .ee_mod_gee_familia(familia, ligacao)
  f <- .ee_mod_formula(vd, preditores)
  fi <- .ee_mod_gee_fit(f, d, fam, "independence")$soma
  linhas <- list(); avisos <- prep$notas[grepl("e entra como contínua", prep$notas, fixed = TRUE)]
  for (cs in correlacoes) {
    fit <- tryCatch(if (cs == "independence") fi else .ee_mod_gee_fit(f, d, fam, cs)$soma, error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$geese$error == 0)) { avisos <- c(avisos, paste0(cs, ": não convergiu")); next }
    q <- .ee_mod_qic(fit, fi)
    a <- fit$geese$alpha
    linhas[[length(linhas) + 1]] <- data.frame(correlacao = cs, QIC = q[["QIC"]], QICu = q[["QICu"]], CIC = q[["CIC"]],
                                               alfa_medio = if (length(a)) mean(a) else 0, n_alfa = length(a))
  }
  tab <- do.call(rbind, linhas)
  tab$delta_QIC <- tab$QIC - min(tab$QIC)
  T <- max(d$.onda)
  if (T == 2) avisos <- c(avisos, "Com 2 tempos, não estruturada, intercambiável e AR(1) são a mesma estrutura: diferenças de QIC entre elas são numéricas.")
  mq <- tab$correlacao[which.min(tab$QIC)]; mc <- tab$correlacao[which.min(tab$CIC)]
  apa <- paste0(.ee_mod_t(idioma, "Menor QIC: correlação ", "Lowest QIC: "), .ee_mod_rot_corr(mq, idioma), " (QIC = ",
                .pr_num(min(tab$QIC), 2, TRUE, decimal), ")",
                if (mq != mc) paste0(.ee_mod_t(idioma, "; o CIC prefere ", "; CIC prefers "), .ee_mod_rot_corr(mc, idioma)) else "", ".")
  .ee_resultado("Escolha da correlação de trabalho no GEE (QIC)", n = max(d$.id), n_excluidos = prep$n_excluidos, tabela = tab,
                decisoes = c("QIC (Pan, 2001) para a correlação de trabalho, menor é melhor, sempre dentro da mesma família e ligação; o CIC (Hin e Wang, 2009), que é o termo de penalidade do QIC, sai ao lado porque escolhe a estrutura com menos erro.",
                             "QICu (= -2Q + 2p) serve para comparar preditores com a mesma correlação; corresponde ao QICC impresso pelo SPSS segundo a documentação do SPSS (equivalência não conferida aqui). O \"QICC\" do geepack é outra coisa (QIC com correção para amostra pequena).",
                             "A quasi-verossimilhança é dividida pela escala do modelo de independência, como em Pan (2001); por isso os valores diferem do geepack::QIC para normal e gama."),
                avisos = avisos, apa = apa)
}

# ================================================================ multinível: auxiliares
# Centralização para modelos multinível.
#   tipo = "geral" (padrão do curso): X_cgm = X - média geral; para nível 2 a média é
#     tomada entre clusters (um valor por cluster). Para preditor de nível 1 cria também,
#     ao lado, a alternativa de Enders e Tofighi (2007): X_cwc e X_media_cgm.
#   tipo = "grupo": só X_cwc = X - média do cluster (efeito intragrupo puro) e
#     X_media_cgm = média do cluster - média geral (efeito entre grupos). Só nível 1.
# Devolve o banco de forma invisível, com classe cujo print mostra só agregados.
centralizar_multinivel <- function(dados, vars, cluster, tipo = c("geral", "grupo")) {
  tipo <- match.arg(tipo)
  if (!is.data.frame(dados)) stop("dados precisa ser um data.frame")
  falta <- setdiff(c(vars, cluster), names(dados))
  if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "))
  dados <- as.data.frame(dados)
  g <- dados[[cluster]]
  novas <- character(); niveis <- character()
  for (v in vars) {
    x <- dados[[v]]
    if (!is.numeric(x)) stop(v, " precisa ser numérica para centralizar")
    ok <- !is.na(x) & !is.na(g)
    distintos <- tapply(x[ok], droplevels(as.factor(g[ok])), function(z) length(unique(z)))
    nivel2 <- all(distintos <= 1)
    niveis[v] <- if (nivel2) "nível 2" else "nível 1"
    m_g <- stats::ave(ifelse(ok, x, NA), g, FUN = function(z) mean(z, na.rm = TRUE))
    if (tipo == "grupo" && nivel2) stop(v, " é constante dentro de cada cluster (nível 2): centralização no grupo zera a variável; use tipo = \"geral\"")
    if (tipo == "geral") {
      centro <- if (nivel2) mean(tapply(x[ok], droplevels(as.factor(g[ok])), mean)) else mean(x, na.rm = TRUE)
      dados[[paste0(v, "_cgm")]] <- x - centro
      novas <- c(novas, paste0(v, "_cgm"))
    }
    if (!nivel2) {
      dados[[paste0(v, "_cwc")]] <- x - m_g
      dados[[paste0(v, "_media_cgm")]] <- m_g - mean(x, na.rm = TRUE)
      novas <- c(novas, paste0(v, "_cwc"), paste0(v, "_media_cgm"))
    }
  }
  attr(dados, "centralizacao") <- list(tipo = tipo, cluster = cluster, vars = vars, niveis = niveis, novas = novas, n = nrow(dados))
  class(dados) <- c("pr_centralizacao", "data.frame")
  invisible(dados)
}

# Print só com agregados: nenhuma linha de participante.
print.pr_centralizacao <- function(x, ...) {
  info <- attr(x, "centralizacao")
  if (is.null(info) || !identical(nrow(x), info$n)) {  # subconjunto (head, [ ]): só dimensões, nunca linhas
    cat("pr_centralizacao: ", nrow(x), " linhas; colunas: ", paste(names(x), collapse = ", "), "\n", sep = "")
    return(invisible(x))
  }
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  cat("== Centralização multinível (tipo = \"", info$tipo, "\"; ", nrow(d), " linhas, ",
      length(unique(stats::na.omit(d[[info$cluster]]))), " clusters)\n", sep = "")
  cat("Variáveis: ", paste(paste0(info$vars, " (", info$niveis[info$vars], ")"), collapse = ", "), "\n", sep = "")
  res <- do.call(rbind, lapply(info$novas, function(v) {
    z <- d[[v]]; mc <- tapply(z, d[[info$cluster]], mean, na.rm = TRUE)
    data.frame(coluna = v, media = mean(z, na.rm = TRUE), dp = stats::sd(z, na.rm = TRUE),
               minimo = min(z, na.rm = TRUE), maximo = max(z, na.rm = TRUE),
               dp_medias_cluster = stats::sd(mc, na.rm = TRUE))
  }))
  num <- vapply(res, is.numeric, logical(1)); res[num] <- lapply(res[num], function(z) round(z, 4))
  print(res, row.names = FALSE)
  cat("Colunas _cwc: desvio da média do cluster (efeito intragrupo); _media_cgm: média do cluster centrada (efeito entre grupos); _cgm: média geral.\n")
  cat("Os dados ficam no objeto; use-o como banco (ex.: dados <- centralizar_multinivel(...)).\n")
  invisible(x)
}

.ee_mod_mlm_prepara <- function(dados, vd, cluster, nivel1, nivel2, inclinacao, interacao, referencias, fatores = NULL) {
  vars <- unique(c(vd, cluster, nivel1, nivel2))
  .ee_mod_nome_ok(vars)
  if (length(inclinacao) > 1) stop("inclinacao aceita um preditor de nível 1")
  if (length(inclinacao) && !inclinacao %in% nivel1) stop("a inclinação aleatória precisa ser de um preditor de nível 1 (está em nivel1?)")
  if (length(interacao) && !length(inclinacao)) stop("a interação entre níveis (M4) pede a inclinação aleatória do preditor moderado (M3)")
  if (length(interacao) && !all(interacao %in% nivel2)) stop("o moderador da interação entre níveis precisa ser de nível 2 (está em nivel2?)")
  dd <- .ee_dados(dados, vars)
  d <- dd$d
  d[[cluster]] <- droplevels(as.factor(d[[cluster]]))
  d <- .ee_mod_fatores(d, setdiff(vars, c(vd, cluster)), referencias, fatores)
  if (length(inclinacao) && !is.numeric(d[[inclinacao]]))
    stop("inclinação aleatória de preditor categórico (", inclinacao, ") não está implementada: a variância e a faixa das ",
         "inclinações pedem um só coeficiente. Use uma variável 0/1 numérica (dummy) para dois grupos, ou ajuste à mão com lme4::lmer.")
  avisos <- .ee_mod_aviso_numericos(d, setdiff(vars, c(vd, cluster)))
  for (v in nivel2) {
    distintos <- tapply(as.character(d[[v]]), d[[cluster]], function(z) length(unique(z)))
    if (any(distintos > 1)) stop(v, " varia dentro dos clusters: é de nível 1, não de nível 2")
  }
  for (v in nivel1) {
    distintos <- tapply(as.character(d[[v]]), d[[cluster]], function(z) length(unique(z)))
    if (all(distintos <= 1)) avisos <- c(avisos, paste0(v, " não varia dentro dos clusters: comporta-se como preditor de nível 2."))
    if (is.numeric(d[[v]]) && abs(mean(d[[v]])) > .05 * stats::sd(d[[v]]))
      avisos <- c(avisos, paste0(v, " não está centralizado (média ", round(mean(d[[v]]), 2), "): o intercepto se refere a ", v, " = 0. Ver centralizar_multinivel()."))
  }
  list(d = d, n_excluidos = dd$n_excluidos, avisos = avisos)
}

# Sequência do curso (páginas 0-1-2-3-4, prática e reporte): M0 nulo, M1 nível 1,
# M2 nível 2, M3 inclinação aleatória, M4 interação entre níveis.
.ee_mod_mlm_sequencia <- function(vd, cluster, nivel1, nivel2, inclinacao, interacao, idioma) {
  ri <- paste0("(1 | ", cluster, ")")
  rs <- if (length(inclinacao)) paste0("(1 + ", inclinacao, " | ", cluster, ")") else NULL
  seq <- list(M0 = list(fixos = character(), aleat = ri, adiciona = .ee_mod_t(idioma, "nulo (só intercepto)", "null (intercept only)")))
  fx <- character()
  if (length(nivel1)) { fx <- c(fx, nivel1); seq$M1 <- list(fixos = fx, aleat = ri, adiciona = .ee_mod_t(idioma, "preditores de nível 1", "level-1 predictors")) }
  if (length(nivel2)) { fx <- c(fx, nivel2); seq$M2 <- list(fixos = fx, aleat = ri, adiciona = .ee_mod_t(idioma, "preditores de nível 2", "level-2 predictors")) }
  if (length(inclinacao)) seq$M3 <- list(fixos = fx, aleat = rs, adiciona = paste0(.ee_mod_t(idioma, "inclinação aleatória de ", "random slope of "), inclinacao))
  if (length(interacao)) seq$M4 <- list(fixos = c(fx, paste0(inclinacao, ":", interacao)), aleat = rs,
                                        adiciona = .ee_mod_t(idioma, "interação entre níveis", "cross-level interaction"))
  lapply(seq, function(m) { m$formula <- stats::as.formula(paste(vd, "~", paste(c(if (length(m$fixos)) m$fixos else "1", m$aleat), collapse = " + "))); m })
}

# Variâncias de um modelo lme4: componentes, R² de Nakagawa e Schielzeth (2013) e
# Nakagawa et al. (2017) (variância dos aleatórios = média de z_i' Σ z_i).
.ee_mod_mlm_var <- function(fit, residual = NULL) {
  X <- lme4::getME(fit, "X"); b <- lme4::fixef(fit)
  vf <- stats::var(as.vector(X %*% b))
  vc <- lme4::VarCorr(fit)
  vr <- 0
  for (g in names(vc)) {
    S <- vc[[g]]
    nm <- colnames(S)
    Z <- matrix(1, nrow(X), length(nm))
    for (j in seq_along(nm)) if (nm[j] != "(Intercept)") Z[, j] <- X[, nm[j]]
    vr <- vr + mean(rowSums((Z %*% S) * Z))
  }
  ve <- if (is.null(residual)) stats::sigma(fit)^2 else residual
  c(fixos = vf, aleatorios = vr, residuo = ve, r2_marginal = vf / (vf + vr + ve), r2_condicional = (vf + vr) / (vf + vr + ve))
}

.ee_mod_mlm_comp <- function(fits, seq, residual = NULL) {
  linhas <- list()
  for (i in seq_along(fits)) {
    m <- fits[[i]]; l <- stats::logLik(m)
    v <- .ee_mod_mlm_var(m, residual)
    linha <- data.frame(modelo = names(fits)[i], adiciona = seq[[i]]$adiciona, n_parametros = attr(l, "df"),
                        loglik = as.numeric(l), desvio = -2 * as.numeric(l), AIC = stats::AIC(m), BIC = stats::BIC(m),
                        delta_desvio = NA_real_, delta_gl = NA_real_, p_lrt = NA_real_, p_mistura = NA_real_,
                        r2_marginal = v[["r2_marginal"]], r2_condicional = v[["r2_condicional"]])
    if (i > 1) {
      ant <- linhas[[i - 1]]
      dd <- ant$desvio - linha$desvio; dg <- linha$n_parametros - ant$n_parametros
      linha$delta_desvio <- dd; linha$delta_gl <- dg
      linha$p_lrt <- stats::pchisq(dd, dg, lower.tail = FALSE)
      if (grepl("^M3", names(fits)[i])) {
        # variância na fronteira (Self e Liang, 1987; Stram e Lee, 1994): mistura 50:50 de χ²(gl-1) e χ²(gl)
        p_baixo <- if (dg - 1 > 0) stats::pchisq(dd, dg - 1, lower.tail = FALSE) else as.numeric(dd <= 0)
        linha$p_mistura <- .5 * p_baixo + .5 * linha$p_lrt
      }
    }
    linhas[[i]] <- linha
  }
  do.call(rbind, linhas)
}

.ee_mod_mlm_diag <- function(fits) {
  sing <- names(fits)[vapply(fits, lme4::isSingular, logical(1))]
  conv <- names(fits)[vapply(fits, function(m) length(m@optinfo$conv$lme4$messages) > 0, logical(1))]
  list(singulares = sing, sem_convergir = conv)
}

.ee_mod_mlm_cluster_pres <- function(d, cluster, cross) {
  tam <- as.vector(table(d[[cluster]])); J <- length(tam)
  data.frame(checagem = c("Clusters", "Tamanho dos clusters", "Clusters com menos de 10 casos"),
             estatistica = c(J, mean(tam), sum(tam < 10)), gl = NA, p = NA,
             leitura = c(if (J < 30) "poucos clusters: EPs de nível 2 viesados (Maas e Hox, 2005, com 50 grupos ou menos); simular o desenho" else
                           if (J <= 50) "até 50 clusters: EPs de nível 2 podem ficar viesados (Maas e Hox, 2005)" else "clusters suficientes para efeitos fixos e EPs de nível 2",
                         paste0("mín. ", min(tam), ", máx. ", max(tam), if (cross && J < 50) "; interação entre níveis pede cerca de 50 clusters" else ""),
                         if (any(tam < 10)) "clusters pequenos pesam pouco nas inclinações aleatórias" else "sem clusters pequenos"))
}

.ee_mod_mlm_escolhe_final <- function(comp, modelo_final) {
  if (!is.null(modelo_final)) {
    if (!modelo_final %in% comp$modelo) stop("modelo_final precisa ser um destes: ", paste(comp$modelo, collapse = ", "))
    return(modelo_final)
  }
  pp <- ifelse(is.na(comp$p_mistura), comp$p_lrt, comp$p_mistura)
  melhora <- which(!is.na(pp) & pp < .05)
  if (!length(melhora)) "M0" else comp$modelo[max(melhora)]
}

# ================================================================ multinível linear
multinivel_linear <- function(dados, vd, cluster, nivel1 = NULL, nivel2 = NULL, inclinacao = NULL, interacao = NULL,
                              referencias = NULL, fatores = NULL, modelo_final = NULL, estimador_final = c("REML", "ML"),
                              nivel = .95, decimal = ".", idioma = "pt") {
  estimador_final <- match.arg(estimador_final)
  .pr_exige("lme4", "lmerTest")
  prep <- .ee_mod_mlm_prepara(dados, vd, cluster, nivel1, nivel2, inclinacao, interacao, referencias, fatores)
  d <- prep$d; avisos <- prep$avisos
  if (!is.numeric(d[[vd]])) stop("o desfecho precisa ser numérico (para 0/1 use multinivel_logistico)")
  seqm <- .ee_mod_mlm_sequencia(vd, cluster, nivel1, nivel2, inclinacao, interacao, idioma)
  op <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(op), add = TRUE)
  ajusta <- function(fm, reml) suppressMessages(lmerTest::lmer(fm, data = d, REML = reml))
  fits <- lapply(seqm, function(m) ajusta(m$formula, FALSE))
  comp <- .ee_mod_mlm_comp(fits, seqm)
  dg <- .ee_mod_mlm_diag(fits)
  # ICC do nulo em ML e REML; efeito de desenho
  m0r <- ajusta(seqm$M0$formula, TRUE)
  icc <- function(m) { v <- as.data.frame(lme4::VarCorr(m)); v$vcov[1] / sum(v$vcov) }
  icc_reml <- icc(m0r); icc_ml <- icc(fits$M0)
  J <- nlevels(d[[cluster]]); mbar <- nrow(d) / J
  deff <- 1 + (mbar - 1) * icc_reml
  # modelo final
  fin <- .ee_mod_mlm_escolhe_final(comp, modelo_final)
  mf <- if (estimador_final == "REML") ajusta(seqm[[fin]]$formula, TRUE) else fits[[fin]]
  sm <- summary(mf)$coefficients
  tq <- stats::qt(1 - (1 - nivel) / 2, sm[, "df"])
  fixos <- data.frame(termo = rownames(sm), b = sm[, "Estimate"], ep = sm[, "Std. Error"], gl = sm[, "df"],
                      t = sm[, "t value"], p = sm[, "Pr(>|t|)"], ic_inf = sm[, "Estimate"] - tq * sm[, "Std. Error"],
                      ic_sup = sm[, "Estimate"] + tq * sm[, "Std. Error"], row.names = NULL)
  vcf <- as.data.frame(lme4::VarCorr(mf))
  vfin <- .ee_mod_mlm_var(mf)
  tau00 <- vcf$vcov[vcf$grp == cluster & vcf$var1 == "(Intercept)" & is.na(vcf$var2)]
  sig2 <- vcf$vcov[vcf$grp == "Residual"]
  efeito <- data.frame(indice = c("ICC (M0, REML)", "ICC (M0, ML)", "Efeito de desenho (M0)",
                                  paste0("R² marginal (", fin, ")"), paste0("R² condicional (", fin, ")"),
                                  paste0("τ00 variância do intercepto (", fin, ")"), paste0("σ² residual (", fin, ")")),
                       valor = c(icc_reml, icc_ml, deff, vfin[["r2_marginal"]], vfin[["r2_condicional"]], tau00, sig2),
                       ic_inf = NA_real_, ic_sup = NA_real_, principal = c(TRUE, rep(FALSE, 6)))
  tem_inc <- length(inclinacao) && fin %in% c("M3", "M4")
  if (tem_inc) {
    tau11 <- vcf$vcov[vcf$grp == cluster & vcf$var1 == inclinacao & is.na(vcf$var2)]
    rho <- vcf$sdcor[vcf$grp == cluster & !is.na(vcf$var2)]
    g10 <- lme4::fixef(mf)[[inclinacao]]
    efeito <- rbind(efeito, data.frame(
      indice = c(paste0("τ11 variância da inclinação de ", inclinacao), "ρ01 correlação intercepto-inclinação",
                 paste0("Faixa de 95% das inclinações de ", inclinacao), paste0("Faixa de 68% das inclinações de ", inclinacao, " (curso)")),
      valor = c(tau11, rho, g10, g10), ic_inf = c(NA, NA, g10 - 1.96 * sqrt(tau11), g10 - sqrt(tau11)),
      ic_sup = c(NA, NA, g10 + 1.96 * sqrt(tau11), g10 + sqrt(tau11)), principal = FALSE))
  }
  pres <- .ee_mod_mlm_cluster_pres(d, cluster, length(interacao) > 0)
  r1 <- stats::residuals(mf)
  pres <- rbind(pres, data.frame(
    checagem = c("Ajuste singular", "Convergência", "Resíduos de nível 1 (assimetria)", "Resíduos de nível 1 (curtose)"),
    estatistica = c(length(dg$singulares), length(dg$sem_convergir), .ee_assimetria(r1), .ee_curtose(r1)), gl = NA, p = NA,
    leitura = c(if (length(dg$singulares)) paste0("singular em ", paste(dg$singulares, collapse = ", "), ": variância de algum efeito aleatório perto de 0 ou correlação em ±1") else "nenhum modelo singular",
                if (length(dg$sem_convergir)) paste0("aviso do otimizador em ", paste(dg$sem_convergir, collapse = ", ")) else "todos convergiram",
                vapply(c(.ee_assimetria(r1), .ee_curtose(r1)), function(x) if (abs(x) < 1) "abaixo de 1 em módulo: afastamento leve" else if (abs(x) <= 2) "entre 1 e 2: afastamento moderado" else "acima de 2: afastamento que preocupa", character(1)))))
  if (length(dg$singulares)) avisos <- c(avisos, paste0("Ajuste singular em ", paste(dg$singulares, collapse = ", "), ": simplificar a parte aleatória desse modelo."))
  if (icc_reml < .05) avisos <- c(avisos, "ICC abaixo de .05: o multinível continua apropriado com dados aninhados (o curso não usa o ICC como porta de entrada), mas τ00 perto de 0 pode gerar ajuste singular.")
  # textos
  pct <- function(x) paste0(.pr_num(100 * x, 1, TRUE, decimal), "%")
  apa <- paste0(.ee_mod_t(idioma, "Modelo nulo: ICC = ", "Null model: ICC = "), .pr_num(icc_reml, 3, FALSE, decimal),
                .ee_mod_t(idioma, " (REML; ", " (REML; "), .pr_num(icc_ml, 3, FALSE, decimal), .ee_mod_t(idioma, " em ML), ou seja, ", " under ML), i.e., "),
                pct(icc_reml), .ee_mod_t(idioma, " da variância de ", " of the variance of "), vd,
                .ee_mod_t(idioma, " está entre os clusters (", " lies between clusters ("), cluster, "); ",
                .ee_mod_t(idioma, "efeito de desenho = ", "design effect = "), .pr_num(deff, 2, TRUE, decimal), ".")
  for (i in seq_len(nrow(comp))[-1]) {
    l <- comp[i, ]
    pmix <- if (!is.na(l$p_mistura)) paste0(.ee_mod_t(idioma, "; mistura de χ² ", "; χ² mixture "), .pr_p(l$p_mistura, decimal)) else ""
    apa <- c(apa, paste0(l$modelo, " (", l$adiciona, ") vs ", comp$modelo[i - 1], ": Δ-2LL = ", .pr_num(l$delta_desvio, 2, TRUE, decimal),
                         .ee_sep(decimal), .ee_mod_t(idioma, "Δgl = ", "Δdf = "), l$delta_gl, .ee_sep(decimal), .pr_p(l$p_lrt, decimal), pmix,
                         .ee_sep(decimal), "AIC = ", .pr_num(l$AIC, 2, TRUE, decimal), .ee_sep(decimal), "BIC = ", .pr_num(l$BIC, 2, TRUE, decimal), "."))
  }
  fx <- fixos[fixos$termo != "(Intercept)", , drop = FALSE]
  if (nrow(fx)) apa <- c(apa, paste0(.ee_mod_t(idioma, "Modelo final ", "Final model "), fin, " (", estimador_final, "; ",
                                      .ee_mod_t(idioma, "gl de Satterthwaite): ", "Satterthwaite df): "),
                                      paste(vapply(seq_len(nrow(fx)), function(i) paste0(gsub(":", " × ", fx$termo[i], fixed = TRUE), ", b = ", .pr_num(fx$b[i], 2, TRUE, decimal), .ee_sep(decimal),
                                                                                       .ee_estat("t", fx$t[i], fx$gl[i], decimal), .ee_sep(decimal), .pr_p(fx$p[i], decimal), .ee_sep(decimal),
                                                                                       .ee_ic(fx$ic_inf[i], fx$ic_sup[i], 2, TRUE, decimal, nivel, idioma)), character(1)), collapse = "; "), "."))
  if (tem_inc) apa <- c(apa, paste0(.ee_mod_t(idioma, "Variância da inclinação de ", "Variance of the slope of "), inclinacao, " = ",
                                    .pr_num(tau11, 3, TRUE, decimal), " (", .ee_mod_t(idioma, "DP", "SD"), " = ", .pr_num(sqrt(tau11), 3, TRUE, decimal), "); ",
                                    .ee_mod_t(idioma, "95% dos clusters com inclinação entre ", "95% of clusters with slopes between "),
                                    .pr_num(g10 - 1.96 * sqrt(tau11), 2, TRUE, decimal), .ee_mod_t(idioma, " e ", " and "), .pr_num(g10 + 1.96 * sqrt(tau11), 2, TRUE, decimal),
                                    "; ", .ee_mod_t(idioma, "correlação intercepto-inclinação = ", "intercept-slope correlation = "), .pr_num(rho, 2, FALSE, decimal), "."))
  decisoes <- c("Sequência do curso (páginas 0-1-2-3-4, prática no JASP e reporte): M0 nulo, M1 nível 1, M2 nível 2, M3 inclinação aleatória, M4 interação entre níveis. A página 'Ajuste de Modelos Multiníveis' inverte M2 e M3; segue-se a ordem das outras três.",
                "Comparações com todos os modelos em ML (REML = FALSE), porque as partes fixas mudam; desvio = -2 logLik (o checklist do curso escreve 2 × logLik, erro de sinal).",
                paste0("Estimativas finais do ", fin, " em ", estimador_final, if (estimador_final == "REML") " (menos viés nas variâncias; o curso não fala da escolha)" else "", "; p dos efeitos fixos com gl de Satterthwaite (lmerTest)."),
                if (is.null(modelo_final)) paste0("Modelo final = último passo com melhora significativa no LRT (", fin, "); a escolha teórica prevalece (argumento modelo_final).") else paste0("Modelo final escolhido pelo usuário: ", fin, "."),
                "LRT de variância (M3 vs anterior) fica na fronteira do espaço de parâmetros: o p ingênuo é conservador; a mistura 50:50 de χ²(gl-1) e χ²(gl) sai ao lado (Self e Liang, 1987; Stram e Lee, 1994).",
                "ICC do nulo sem limiar de entrada; ICC = τ00/(τ00 + σ²). Efeito de desenho = 1 + (m - 1) ICC com m = tamanho médio do cluster; a regra 'deff < 2 dispensa multinível' não vale com preditores de nível 2 ou clusters com menos de 10 casos (Lai e Kwok, 2015).",
                "R² marginal e condicional de Nakagawa e Schielzeth (2013), com inclinações aleatórias pela extensão de Nakagawa et al. (2017).",
                if (length(interacao)) "M4 mantém a inclinação aleatória do preditor moderado (como no reporte do curso); tirar a inclinação aleatória infla o erro tipo I do termo entre níveis." else NULL,
                "Centralização é passo anterior (centralizar_multinivel, padrão do curso = média geral); para efeito intragrupo puro e para a interação entre níveis, usar as colunas _cwc e _media_cgm que a função cria ao lado (Enders e Tofighi, 2007).")
  .ee_resultado("Regressão linear multinível", n = nrow(d), n_excluidos = prep$n_excluidos, tabela = fixos, efeito = efeito,
                pressupostos = pres, posthoc = comp, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(modelo = mf, modelos_ml = fits, modelo_final = fin, n_clusters = J, icc_reml = icc_reml, icc_ml = icc_ml,
                             efeitos_aleatorios = vcf))
}

# ================================================================ multinível logístico
multinivel_logistico <- function(dados, vd, cluster, nivel1 = NULL, nivel2 = NULL, inclinacao = NULL, interacao = NULL,
                                 referencias = NULL, fatores = NULL, evento = NULL, modelo_final = NULL, nAGQ = 1,
                                 nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("lme4")
  prep <- .ee_mod_mlm_prepara(dados, vd, cluster, nivel1, nivel2, inclinacao, interacao, referencias, fatores)
  d <- prep$d; avisos <- prep$avisos; decisoes <- character()
  y <- d[[vd]]
  if (!is.null(evento)) {
    d[[vd]] <- as.integer(as.character(y) == as.character(evento))
    decisoes <- c(decisoes, paste0("Evento (1) = '", evento, "'."))
  } else if (!is.numeric(y)) {
    b01 <- .ee_mod_y01(d, vd, idioma); d <- b01$d; decisoes <- c(decisoes, b01$nota)
  }
  if (!all(d[[vd]] %in% c(0, 1))) stop("o desfecho precisa ser 0/1 (ou fator com 2 níveis, ou indique evento)")
  seqm <- .ee_mod_mlm_sequencia(vd, cluster, nivel1, nivel2, inclinacao, interacao, idioma)
  op <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(op), add = TRUE)
  ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  ajusta <- function(m) {
    q <- if (grepl("^\\(1 \\|", m$aleat)) nAGQ else 1
    lme4::glmer(m$formula, data = d, family = stats::binomial("logit"), nAGQ = q, control = ctrl)
  }
  fits <- lapply(seqm, ajusta)
  comp <- .ee_mod_mlm_comp(fits, seqm, residual = pi^2 / 3)
  # com nAGQ > 1 os modelos só com intercepto aleatório saem por quadratura e os com
  # inclinação aleatória por Laplace: o LRT do passo que acrescenta a inclinação usa o
  # modelo anterior reajustado por Laplace (mesma aproximação nos dois lados).
  if (nAGQ > 1 && "M3" %in% names(fits)) {
    i3 <- which(names(fits) == "M3")
    ant_l <- lme4::glmer(seqm[[i3 - 1]]$formula, data = d, family = stats::binomial("logit"), nAGQ = 1, control = ctrl)
    dd3 <- -2 * as.numeric(stats::logLik(ant_l)) - comp$desvio[i3]
    dg3 <- comp$n_parametros[i3] - attr(stats::logLik(ant_l), "df")
    comp$delta_desvio[i3] <- dd3; comp$delta_gl[i3] <- dg3
    comp$p_lrt[i3] <- stats::pchisq(dd3, dg3, lower.tail = FALSE)
    comp$p_mistura[i3] <- .5 * (if (dg3 > 1) stats::pchisq(dd3, dg3 - 1, lower.tail = FALSE) else as.numeric(dd3 <= 0)) + .5 * comp$p_lrt[i3]
    decisoes <- c(decisoes, paste0("nAGQ = ", nAGQ, " nos modelos só com intercepto aleatório; os com inclinação aleatória usam Laplace. O LRT de M3 compara com ",
                                   names(fits)[i3 - 1], " reajustado por Laplace, para não misturar aproximações da verossimilhança; AIC e BIC da tabela seguem a aproximação de cada modelo."))
  }
  dg <- .ee_mod_mlm_diag(fits)
  v0 <- as.data.frame(lme4::VarCorr(fits$M0))
  tau00_0 <- v0$vcov[1]; icc <- tau00_0 / (tau00_0 + pi^2 / 3)
  fin <- .ee_mod_mlm_escolhe_final(comp, modelo_final)
  mf <- fits[[fin]]
  sm <- summary(mf)$coefficients
  z <- stats::qnorm(1 - (1 - nivel) / 2)
  fixos <- data.frame(termo = rownames(sm), b = sm[, 1], ep = sm[, 2], z = sm[, 3], p = sm[, 4],
                      OR = exp(sm[, 1]), OR_ic_inf = exp(sm[, 1] - z * sm[, 2]), OR_ic_sup = exp(sm[, 1] + z * sm[, 2]),
                      row.names = NULL)
  vfin <- .ee_mod_mlm_var(mf, residual = pi^2 / 3)
  vcf <- as.data.frame(lme4::VarCorr(mf))
  tau00 <- vcf$vcov[vcf$var1 == "(Intercept)" & is.na(vcf$var2)]
  efeito <- rbind(
    data.frame(indice = "ICC latente (M0) = τ00/(τ00 + π²/3)", valor = icc, ic_inf = NA, ic_sup = NA, principal = TRUE),
    if (nrow(fixos) > 1) data.frame(indice = paste0("OR: ", fixos$termo[-1]), valor = fixos$OR[-1], ic_inf = fixos$OR_ic_inf[-1],
                                    ic_sup = fixos$OR_ic_sup[-1], principal = FALSE),
    data.frame(indice = c(paste0("R² marginal (", fin, ")"), paste0("R² condicional (", fin, ")"), paste0("τ00 (", fin, ", escala logit)")),
               valor = c(vfin[["r2_marginal"]], vfin[["r2_condicional"]], tau00), ic_inf = NA, ic_sup = NA, principal = FALSE))
  pres <- .ee_mod_mlm_cluster_pres(d, cluster, length(interacao) > 0)
  pres <- rbind(pres, data.frame(checagem = c("Ajuste singular", "Convergência", "Eventos (proporção de 1)"),
                                 estatistica = c(length(dg$singulares), length(dg$sem_convergir), mean(d[[vd]])), gl = NA, p = NA,
                                 leitura = c(if (length(dg$singulares)) paste0("singular em ", paste(dg$singulares, collapse = ", ")) else "nenhum modelo singular",
                                             if (length(dg$sem_convergir)) paste0("aviso do otimizador em ", paste(dg$sem_convergir, collapse = ", ")) else "todos convergiram",
                                             if (min(mean(d[[vd]]), 1 - mean(d[[vd]])) < .05) "evento raro: estimativas instáveis" else "proporção adequada")))
  apa <- paste0(.ee_mod_t(idioma, "Modelo nulo: ICC na escala latente = ", "Null model: latent-scale ICC = "), .pr_num(icc, 3, FALSE, decimal),
                " (τ00 = ", .pr_num(tau00_0, 3, TRUE, decimal), .ee_mod_t(idioma, "; variância de nível 1 fixada em π²/3 = ", "; level-1 variance fixed at π²/3 = "),
                .pr_num(pi^2 / 3, 2, TRUE, decimal), ").")
  for (i in seq_len(nrow(comp))[-1]) {
    l <- comp[i, ]
    pmix <- if (!is.na(l$p_mistura)) paste0(.ee_mod_t(idioma, "; mistura de χ² ", "; χ² mixture "), .pr_p(l$p_mistura, decimal)) else ""
    apa <- c(apa, paste0(l$modelo, " (", l$adiciona, ") vs ", comp$modelo[i - 1], ": Δ-2LL = ", .pr_num(l$delta_desvio, 2, TRUE, decimal),
                         .ee_sep(decimal), .ee_mod_t(idioma, "Δgl = ", "Δdf = "), l$delta_gl, .ee_sep(decimal), .pr_p(l$p_lrt, decimal), pmix, "."))
  }
  fx <- fixos[-1, , drop = FALSE]
  if (nrow(fx)) apa <- c(apa, paste0(.ee_mod_t(idioma, "Modelo final ", "Final model "), fin, ": ",
                                      paste(vapply(seq_len(nrow(fx)), function(i) paste0(gsub(":", " × ", fx$termo[i], fixed = TRUE), ", OR = ", .pr_num(fx$OR[i], 2, TRUE, decimal), .ee_sep(decimal),
                                                                                       .ee_ic(fx$OR_ic_inf[i], fx$OR_ic_sup[i], 2, TRUE, decimal, nivel, idioma), .ee_sep(decimal),
                                                                                       .ee_estat("z", fx$z[i], NULL, decimal), .ee_sep(decimal), .pr_p(fx$p[i], decimal)), character(1)), collapse = "; "),
                                      .ee_mod_t(idioma, " (OR condicionais: comparam casos do mesmo cluster).", " (cluster-specific ORs).")))
  if (length(dg$singulares)) avisos <- c(avisos, paste0("Ajuste singular em ", paste(dg$singulares, collapse = ", "), "."))
  decisoes <- c(decisoes,
                "Mesma sequência M0 a M4 da linear; todos os modelos por máxima verossimilhança (aproximação de Laplace; nAGQ > 1 só com intercepto aleatório). GLMM não tem REML.",
                "OR = exp(b) com IC de Wald. Os OR são condicionais ao cluster (comparam dois casos do mesmo cluster) e ficam maiores em módulo que os OR populacionais de uma logística simples ou do GEE quando τ00 é grande.",
                "ICC na escala latente com variância de nível 1 = π²/3 (Goldstein et al., 2002). Como essa variância é fixa, coeficientes e τ00 mudam de escala ao entrar preditores: comparar τ00 entre modelos como na linear não vale sem reescala.",
                "R² marginal e condicional de Nakagawa e Schielzeth (2013) com variância residual π²/3.",
                "Comparações por LRT (Δ-2LL); para a inclinação aleatória, mistura de χ² ao lado (fronteira).")
  .ee_resultado("Regressão logística multinível", n = nrow(d), n_excluidos = prep$n_excluidos, tabela = fixos, efeito = efeito,
                pressupostos = pres, posthoc = comp, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(modelo = mf, modelos = fits, modelo_final = fin, icc_latente = icc))
}

# ================================================================ meta-análise
.ee_mod_meta_transf <- function(medida) {
  switch(medida, OR = exp, RR = exp, COR = tanh, NULL)
}
.ee_mod_meta_rot <- function(medida, idioma) {
  switch(medida, SMD = "g", MD = .ee_mod_t(idioma, "DM", "MD"), OR = "OR", RR = "RR", COR = "r",
         .ee_mod_t(idioma, "efeito", "effect"))
}

# Tamanhos de efeito por estudo com metafor::escalc. SMD = g de Hedges (correção de viés
# exata); OR e RR em log; COR em z de Fisher (volta para r na síntese). Para OR e RR,
# n1 e n2 são os totais de cada grupo e eventos1/eventos2 os eventos.
calcular_efeitos_meta <- function(dados, medida = c("SMD", "MD", "OR", "RR", "COR"), estudo = NULL,
                                  media1 = NULL, dp1 = NULL, n1 = NULL, media2 = NULL, dp2 = NULL, n2 = NULL,
                                  eventos1 = NULL, eventos2 = NULL, r = NULL, n = NULL) {
  medida <- match.arg(medida)
  .pr_exige("metafor")
  if (!is.data.frame(dados)) stop("dados precisa ser um data.frame")
  col <- function(x) { if (is.null(x)) stop("faltou indicar uma coluna para a medida ", medida); if (!x %in% names(dados)) stop("coluna ausente: ", x); dados[[x]] }
  rot <- if (is.null(estudo)) paste("Estudo", seq_len(nrow(dados))) else as.character(dados[[estudo]])
  notas <- character()
  if (medida %in% c("SMD", "MD")) {
    es <- metafor::escalc(measure = medida, m1i = col(media1), sd1i = col(dp1), n1i = col(n1),
                          m2i = col(media2), sd2i = col(dp2), n2i = col(n2))
    ntot <- col(n1) + col(n2)
    if (medida == "SMD") notas <- "SMD = g de Hedges (d com correção de viés exata), variância de grandes amostras."
  } else if (medida %in% c("OR", "RR")) {
    a <- col(eventos1); n_1 <- col(n1); cc <- col(eventos2); n_2 <- col(n2)
    if (any(a > n_1 | cc > n_2)) stop("eventos maiores que o total do grupo")
    es <- metafor::escalc(measure = medida, ai = a, bi = n_1 - a, ci = cc, di = n_2 - cc, drop00 = TRUE)
    ntot <- n_1 + n_2
    notas <- paste0(medida, " em escala log; 0,5 somado às células de tabelas com zero; estudos sem evento nos dois grupos saem da síntese (drop00).")
  } else {
    es <- metafor::escalc(measure = "ZCOR", ri = col(r), ni = col(n))
    ntot <- col(n)
    notas <- "Correlações em z de Fisher, variância 1/(n - 3); a síntese volta para r (tanh)."
  }
  out <- data.frame(estudo = rot, yi = as.numeric(es$yi), vi = as.numeric(es$vi), sei = sqrt(as.numeric(es$vi)), n_total = ntot,
                    stringsAsFactors = FALSE)
  out <- cbind(out, dados[, setdiff(names(dados), names(out)), drop = FALSE])
  attr(out, "medida") <- medida
  if (medida %in% c("OR", "RR")) attr(out, "contagens") <- c(eventos1 = eventos1, n1 = n1, eventos2 = eventos2, n2 = n2)
  attr(out, "notas") <- notas
  out
}

# Teste de Peters et al. (2006): regressão de ln(OR) (ou ln(RR)) em 1/N, ponderada por
# S·F/N (S = eventos, F = não eventos, N = total do estudo); testa a inclinação com
# t(k - 2). É o meta::metabias(method.bias = "Peters").
.ee_mod_peters <- function(yi, a, n1, cc, n2) {
  N <- as.numeric(n1) + as.numeric(n2); S <- as.numeric(a) + as.numeric(cc); F <- N - S
  ok <- S > 0 & F > 0
  fit <- stats::lm(yi[ok] ~ I(1 / N[ok]), weights = (S * F / N)[ok])
  s <- summary(fit)$coefficients
  list(vies = s[2, 1], ep = s[2, 2], t = s[2, 3], gl = fit$df.residual, p = s[2, 4])
}

# Teste de Harbord et al. (2006) para OR: escore Z = a - n1·m1/N e variância
# hipergeométrica V = n1·n2·m1·m2/(N²(N - 1)); regressão de Z/√V em √V (MQO); testa o
# intercepto com t(k - 2).
.ee_mod_harbord <- function(a, n1, cc, n2) {
  a <- as.numeric(a); n1 <- as.numeric(n1); cc <- as.numeric(cc); n2 <- as.numeric(n2)
  N <- n1 + n2; m1 <- a + cc; m2 <- N - m1
  Z <- a - n1 * m1 / N; V <- n1 * n2 * m1 * m2 / (N^2 * (N - 1))
  ok <- V > 0
  fit <- stats::lm(I(Z[ok] / sqrt(V[ok])) ~ sqrt(V[ok]))
  s <- summary(fit)$coefficients
  list(vies = s[1, 1], ep = s[1, 2], t = s[1, 3], gl = fit$df.residual, p = s[1, 4])
}

# Intervalo de predição de Riley et al. (2011) e IntHout et al. (2016):
# μ ± t(k - 2) · √(τ² + EP²), com o EP do modelo ajustado (de Hartung-Knapp, se pedido).
# Não definido com k < 3 nem no modelo de efeito comum. O padrão do metafor usa z.
.ee_mod_meta_ip <- function(res, nivel = .95) {
  if (res$method %in% c("EE", "FE", "CE") || res$k < 3) return(c(NA_real_, NA_real_))
  tq <- stats::qt(1 - (1 - nivel) / 2, res$k - 2)
  as.numeric(res$b) + c(-1, 1) * tq * sqrt(res$tau2 + res$se^2)
}

.ee_mod_meta_linha <- function(res, rotulo, transf) {
  ip <- .ee_mod_meta_ip(res, if (res$level < .5) 1 - res$level else res$level)
  l <- data.frame(modelo = rotulo, k = res$k, estimativa = as.numeric(res$b), ep = res$se, ic_inf = res$ci.lb, ic_sup = res$ci.ub,
                  estatistica = res$zval, gl = if (identical(res$test, "z")) NA_real_ else res$ddf[[1]], p = res$pval,
                  ip_inf = ip[1], ip_sup = ip[2], tau2 = res$tau2, stringsAsFactors = FALSE)
  if (!is.null(transf)) {
    l$efeito <- transf(l$estimativa); l$efeito_ic_inf <- transf(l$ic_inf); l$efeito_ic_sup <- transf(l$ic_sup)
    l$efeito_ip_inf <- transf(l$ip_inf); l$efeito_ip_sup <- transf(l$ip_sup)
  }
  l
}

.ee_mod_meta_dados <- function(dados, yi, vi, estudo, extras = character()) {
  vars <- unique(c(yi, vi, estudo, extras))
  dd <- .ee_dados(dados, vars)
  d <- dd$d
  if (any(d[[vi]] <= 0)) stop("variâncias amostrais (vi) precisam ser > 0")
  d$.estudo <- if (is.null(estudo)) paste("Estudo", seq_len(nrow(d))) else as.character(d[[estudo]])
  list(d = d, n_excluidos = dd$n_excluidos)
}

.ee_mod_i2q <- function(Q, gl) if (Q > 0) 100 * max(0, (Q - gl) / Q) else 0

.ee_mod_faixa_i2 <- function(i2, idioma) {
  .ee_mod_t(idioma,
            if (i2 <= 40) "0 a 40%: pode não ser importante" else if (i2 < 60) "30 a 60%: pode ser moderada" else if (i2 < 75) "50 a 90%: pode ser substancial" else "75 a 100%: considerável",
            if (i2 <= 40) "0-40%: might not be important" else if (i2 < 60) "30-60%: may be moderate" else if (i2 < 75) "50-90%: may be substantial" else "75-100%: considerable")
}

meta_analise <- function(dados, yi = "yi", vi = "vi", estudo = NULL, medida = NULL, metodo = "REML",
                         hartung_knapp = FALSE, sensibilidade = c("PM", "DL"), vies = TRUE, k_min_vies = 10,
                         nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("metafor")
  if (is.null(medida)) medida <- if (!is.null(attr(dados, "medida"))) attr(dados, "medida") else "generica"
  if (is.null(estudo) && "estudo" %in% names(dados)) estudo <- "estudo"
  cont <- attr(dados, "contagens")
  extras <- c(if ("n_total" %in% names(dados)) "n_total", if (!is.null(cont)) unname(cont))
  md <- .ee_mod_meta_dados(dados, yi, vi, estudo, extras)
  d <- md$d; k <- nrow(d)
  if (k < 2) stop("meta-análise precisa de pelo menos 2 estudos")
  transf <- .ee_mod_meta_transf(medida); rot <- .ee_mod_meta_rot(medida, idioma)
  teste <- if (hartung_knapp) "knha" else "z"
  ajusta <- function(m, tst = teste) metafor::rma(yi = d[[yi]], vi = d[[vi]], method = m, test = tst, level = 100 * nivel, slab = d$.estudo)
  res <- ajusta(metodo)
  comum <- metodo %in% c("EE", "FE", "CE")
  rot_mod <- if (comum) .ee_mod_t(idioma, "Efeito comum ", "Common effect ") else .ee_mod_t(idioma, "Efeitos aleatórios ", "Random effects ")
  rot_principal <- paste0(rot_mod, metodo, if (hartung_knapp) " + Hartung-Knapp" else "")
  tab <- .ee_mod_meta_linha(res, rot_principal, transf)
  alt <- ajusta(metodo, if (hartung_knapp) "z" else "knha")
  tab <- rbind(tab, .ee_mod_meta_linha(alt, paste0(rot_mod, metodo, if (!hartung_knapp) " + Hartung-Knapp" else " (z)"), transf))
  for (m in setdiff(sensibilidade, metodo)) {
    r2 <- tryCatch(ajusta(m), error = function(e) NULL)
    if (!is.null(r2)) tab <- rbind(tab, .ee_mod_meta_linha(r2, paste0(.ee_mod_t(idioma, "Sensibilidade: ", "Sensitivity: "), m), transf))
  }
  if (!comum) {
    ee <- metafor::rma(yi = d[[yi]], vi = d[[vi]], method = "EE", level = 100 * nivel)
    tab <- rbind(tab, .ee_mod_meta_linha(ee, .ee_mod_t(idioma, "Sensibilidade: efeito comum", "Sensitivity: common effect"), transf))
  }
  # heterogeneidade
  Q <- res$QE; glq <- k - 1; i2q <- .ee_mod_i2q(Q, glq)
  ci <- if (comum) NULL else tryCatch(stats::confint(res), error = function(e) NULL)
  tau_ic <- if (!is.null(ci)) ci$random["tau^2", c("ci.lb", "ci.ub")] else c(NA, NA)
  i2_ic <- if (!is.null(ci)) ci$random["I^2(%)", c("ci.lb", "ci.ub")] else c(NA, NA)
  pres <- data.frame(checagem = c("Q de Cochran", "I² (a partir do Q)", "I² (a partir do τ², metafor)", "τ²"),
                     estatistica = c(Q, i2q, res$I2, res$tau2), gl = c(glq, NA, NA, NA), p = c(res$QEp, NA, NA, NA),
                     leitura = c(if (res$QEp < .10) "heterogeneidade além do acaso (critério p < .10 do curso; Q tem pouco poder com poucos estudos)" else "sem evidência de heterogeneidade pelo Q (pouco poder com poucos estudos; não prova homogeneidade)",
                                 paste0(.ee_mod_faixa_i2(i2q, idioma), " (faixas sobrepostas da Cochrane; não usar como corte)"),
                                 paste0("IC ", round(100 * nivel), "%: ", round(i2_ic[1], 1), " a ", round(i2_ic[2], 1), "%"),
                                 paste0("IC ", round(100 * nivel), "% (Q-profile): ", signif(tau_ic[1], 3), " a ", signif(tau_ic[2], 3))))
  avisos <- character()
  if (k < 5) avisos <- c(avisos, paste0("Só ", k, " estudos: τ² muito impreciso; o intervalo de predição (t com ", k - 2, " gl", if (k < 3) ", não definido" else "", ") e o I² são instáveis. Hartung-Knapp e o modelo de efeito comum saem ao lado como sensibilidade."))
  # leave-one-out
  loo <- tryCatch(as.data.frame(metafor::leave1out(res)), error = function(e) NULL)
  if (!is.null(loo)) {
    loo$estudo <- d$.estudo
    sig0 <- res$pval < 1 - nivel
    muda <- loo$estudo[(loo$pval < 1 - nivel) != sig0]
    faixa <- range(loo$estimate)
    pres <- rbind(pres, data.frame(checagem = "Leave-one-out", estatistica = faixa[2] - faixa[1], gl = NA, p = NA,
                                   leitura = paste0("estimativas de ", signif(if (is.null(transf)) faixa[1] else transf(faixa[1]), 3), " a ",
                                                    signif(if (is.null(transf)) faixa[2] else transf(faixa[2]), 3),
                                                    if (length(muda)) paste0("; a significância muda ao tirar: ", paste(muda, collapse = ", ")) else "; nenhuma exclusão muda a significância")))
  }
  # viés de estudos pequenos
  vies_tab <- NULL; tf <- NULL
  if (vies && k >= k_min_vies) {
    eg <- metafor::regtest(d[[yi]], d[[vi]], model = "lm", predictor = "sei", ret.fit = TRUE)
    vies_tab <- data.frame(teste = "Egger, regressão no EP", vies = unname(stats::coef(eg$fit)[2]), estatistica = eg$zval, gl = eg$dfs, p = eg$pval)
    if (medida %in% c("OR", "RR") && !is.null(cont)) {
      pe <- .ee_mod_peters(d[[yi]], d[[cont[["eventos1"]]]], d[[cont[["n1"]]]], d[[cont[["eventos2"]]]], d[[cont[["n2"]]]])
      vies_tab <- rbind(vies_tab, data.frame(teste = "Peters, regressão em 1/N", vies = pe$vies, estatistica = pe$t, gl = pe$gl, p = pe$p))
    }
    if (identical(medida, "OR") && !is.null(cont)) {
      hb <- .ee_mod_harbord(d[[cont[["eventos1"]]]], d[[cont[["n1"]]]], d[[cont[["eventos2"]]]], d[[cont[["n2"]]]])
      vies_tab <- rbind(vies_tab, data.frame(teste = "Harbord, escore", vies = hb$vies, estatistica = hb$t, gl = hb$gl, p = hb$p))
    }
    if (identical(medida, "SMD") && "n_total" %in% names(d)) {
      pn <- metafor::regtest(d[[yi]], d[[vi]], ni = d$n_total, model = "lm", predictor = "sqrtninv", ret.fit = TRUE)
      vies_tab <- rbind(vies_tab, data.frame(teste = "regressão em 1/√n, alternativa para SMD", vies = unname(stats::coef(pn$fit)[2]), estatistica = pn$zval, gl = pn$dfs, p = pn$pval))
    }
    principal_vies <- if (nrow(vies_tab) > 1 && medida %in% c("OR", "RR", "SMD")) 2 else 1
    vb <- vies_tab[principal_vies, ]
    pres <- rbind(pres, data.frame(checagem = paste0("Assimetria: ", vb$teste), estatistica = vb$estatistica, gl = vb$gl, p = vb$p,
                                   leitura = if (vb$p < .05) "assimetria compatível com efeitos de estudos pequenos (viés de publicação é uma das causas possíveis)" else "sem evidência de assimetria (o teste tem pouco poder)"))
    if (medida %in% c("OR", "RR", "SMD")) avisos <- c(avisos, if (nrow(vies_tab) > 1) paste0("Para ", medida, " o Egger clássico tem erro tipo I inflado (o EP se correlaciona com o efeito); a leitura principal usa o teste baseado em n", if (identical(medida, "OR")) " (Peters; Harbord ao lado)" else "", ".") else
      paste0("Para ", medida, " o Egger clássico tem erro tipo I inflado e não há contagens/n para o teste baseado em n: gere os efeitos com calcular_efeitos_meta()."))
    tf <- tryCatch(metafor::trimfill(res), error = function(e) NULL)
    if (!is.null(tf)) {
      tab <- rbind(tab, .ee_mod_meta_linha(tf, paste0(.ee_mod_t(idioma, "Sensibilidade: trim-and-fill (+", "Sensitivity: trim-and-fill (+"), tf$k0, ")"), transf))
      pres <- rbind(pres, data.frame(checagem = "Trim-and-fill (sensibilidade)", estatistica = tf$k0, gl = NA, p = NA,
                                     leitura = paste0(if (tf$k0 == 0) "nenhum estudo imputado" else paste0(tf$k0, if (tf$k0 == 1) " estudo imputado" else " estudos imputados",
                                                                                                         " no lado ", if (identical(tf$side, "left")) "esquerdo" else "direito"),
                                                      "; não é o efeito verdadeiro, só mostra a sensibilidade da conclusão")))
    }
  } else if (vies) {
    avisos <- c(avisos, paste0("Com ", k, " estudos (menos de ", k_min_vies, ") os testes de assimetria do funil e o trim-and-fill não rodam: poder insuficiente (Sterne et al., 2011)."))
  }
  # efeito
  l1 <- tab[1, ]
  val <- if (is.null(transf)) l1$estimativa else l1$efeito
  efeito <- data.frame(indice = c(paste0(rot, .ee_mod_t(idioma, " combinado", " pooled")), .ee_mod_t(idioma, "Intervalo de predição", "Prediction interval")),
                       valor = c(val, val),
                       ic_inf = c(if (is.null(transf)) l1$ic_inf else l1$efeito_ic_inf, if (is.null(transf)) l1$ip_inf else l1$efeito_ip_inf),
                       ic_sup = c(if (is.null(transf)) l1$ic_sup else l1$efeito_ic_sup, if (is.null(transf)) l1$ip_sup else l1$efeito_ip_sup),
                       principal = c(TRUE, FALSE))
  # textos
  lim <- !is.null(transf) && medida == "COR"
  est_txt <- .ee_efeito_txt(rot, efeito$valor[1], efeito$ic_inf[1], efeito$ic_sup[1], limitado = lim, casas = 2, decimal = decimal,
                            idioma = idioma, nivel = nivel)
  estat <- if (hartung_knapp) .ee_estat("t", l1$estatistica, l1$gl, decimal) else .ee_estat("z", l1$estatistica, NULL, decimal)
  tem_ip <- !is.na(efeito$ic_inf[2])
  apa <- paste0(if (comum) .ee_mod_t(idioma, "Meta-análise de efeito comum (", "Common-effect meta-analysis (") else .ee_mod_t(idioma, "Meta-análise de efeitos aleatórios (", "Random-effects meta-analysis ("),
                metodo, if (hartung_knapp) " + Hartung-Knapp" else "", "; k = ", k, "): ", est_txt, .ee_sep(decimal), estat, .ee_sep(decimal), .pr_p(l1$p, decimal),
                if (tem_ip) paste0(.ee_sep(decimal), .ee_mod_t(idioma, "intervalo de predição ", "prediction interval "), round(100 * nivel), "% [",
                                   .pr_num(efeito$ic_inf[2], 2, !lim, decimal), .ee_sep(decimal), .pr_num(efeito$ic_sup[2], 2, !lim, decimal), "]") else "", ".")
  apa <- c(apa, paste0(.ee_mod_t(idioma, "Heterogeneidade: ", "Heterogeneity: "), .ee_estat("Q", Q, glq, decimal), .ee_sep(decimal), .pr_p(res$QEp, decimal),
                       .ee_sep(decimal), "I² = ", .pr_num(i2q, 1, TRUE, decimal), "%",
                       if (!comum) paste0(.ee_sep(decimal), "τ² = ", .pr_num(res$tau2, 3, TRUE, decimal),
                                          if (all(is.finite(tau_ic))) paste0(" [", .pr_num(tau_ic[1], 3, TRUE, decimal), .ee_sep(decimal), .pr_num(tau_ic[2], 3, TRUE, decimal), "]") else "") else "", "."))
  if (!is.null(vies_tab)) {
    vb <- vies_tab[if (nrow(vies_tab) > 1 && medida %in% c("OR", "RR", "SMD")) 2 else 1, ]
    apa <- c(apa, paste0(.ee_mod_t(idioma, "Assimetria do funil (", "Funnel asymmetry ("), sub(",.*$", "", vb$teste), "): ", .ee_estat("t", vb$estatistica, vb$gl, decimal), .ee_sep(decimal), .pr_p(vb$p, decimal),
                         if (!is.null(tf)) paste0(.ee_mod_t(idioma, "; trim-and-fill (sensibilidade) imputou ", "; trim-and-fill (sensitivity) imputed "), tf$k0,
                                                  .ee_mod_t(idioma, if (tf$k0 == 1) " estudo, efeito ajustado " else " estudos, efeito ajustado ",
                                                            if (tf$k0 == 1) " study, adjusted effect " else " studies, adjusted effect "),
                                                  if (rot %in% c("efeito", "effect")) "" else paste0(rot, " "), "= ",
                                                  .pr_num(if (is.null(transf)) as.numeric(tf$b) else transf(as.numeric(tf$b)), 2, !lim, decimal), .ee_sep(decimal),
                                                  .ee_ic(if (is.null(transf)) tf$ci.lb else transf(tf$ci.lb), if (is.null(transf)) tf$ci.ub else transf(tf$ci.ub), 2, !lim, decimal, nivel, idioma)) else "", "."))
  }
  decisoes <- c(attr(dados, "notas"),
                paste0("Efeitos aleatórios como padrão (curso, 'Fixo vs Aleatórios'): a escolha do modelo é pela hipótese sobre os efeitos verdadeiros, não pelo I². τ² por ", metodo, " (REML ou Paule-Mandel em vez de DerSimonian-Laird: Veroniki et al., 2016; Langan et al., 2019)."),
                "Hartung-Knapp calculado ao lado (IC mais largo, reflete a incerteza de τ²; IntHout et al., 2014); vira o principal com hartung_knapp = TRUE.",
                paste0("Sensibilidade com outros estimadores de τ² (", paste(setdiff(sensibilidade, metodo), collapse = ", "), ") e com o modelo de efeito comum, como no curso."),
                "I² principal pela fórmula do curso, (Q - gl)/Q, truncado em 0 (o mesmo do pacote meta); o I² do metafor (baseado em τ²) sai ao lado.",
                "Intervalo de predição = μ ± t(k - 2) · √(τ² + EP²) (Riley et al., 2011; IntHout et al., 2016), relatado sempre que k >= 3; o padrão do metafor (normal) seria mais estreito.",
                if (comum) "Modelo de efeito comum pedido (metodo = \"EE\"): supõe um único efeito verdadeiro; sem τ² nem intervalo de predição." else NULL,
                if (vies) paste0("Assimetria do funil só com k >= ", k_min_vies, " (Sterne et al., 2011); para OR e RR o teste principal é o de Peters et al. (2006), regressão em 1/N com pesos S·F/N, com Harbord ao lado para OR (Harbord et al., 2006); trim-and-fill só como sensibilidade (Duval e Tweedie, 2000; Peters et al., 2007)."))
  zc <- stats::qnorm(1 - (1 - nivel) / 2)
  estudos <- data.frame(estudo = d$.estudo, yi = d[[yi]], ic_inf = d[[yi]] - zc * sqrt(d[[vi]]), ic_sup = d[[yi]] + zc * sqrt(d[[vi]]),
                        peso_aleatorio = 100 * (1 / (d[[vi]] + res$tau2)) / sum(1 / (d[[vi]] + res$tau2)),
                        peso_comum = 100 * (1 / d[[vi]]) / sum(1 / d[[vi]]), stringsAsFactors = FALSE)
  if (!is.null(transf)) { estudos$efeito <- transf(estudos$yi); estudos$efeito_ic_inf <- transf(estudos$ic_inf); estudos$efeito_ic_sup <- transf(estudos$ic_sup) }
  .ee_resultado("Meta-análise", n = k, n_excluidos = md$n_excluidos, tabela = tab, efeito = efeito, pressupostos = pres,
                posthoc = NULL, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(modelo = res, estudos = estudos, leave_one_out = loo, vies = vies_tab, trimfill = tf, medida = medida, i2_q = i2q))
}

meta_subgrupos <- function(dados, subgrupo, yi = "yi", vi = "vi", estudo = NULL, medida = NULL, metodo = "REML",
                           hartung_knapp = FALSE, nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("metafor")
  if (is.null(medida)) medida <- if (!is.null(attr(dados, "medida"))) attr(dados, "medida") else "generica"
  if (is.null(estudo) && "estudo" %in% names(dados)) estudo <- "estudo"
  md <- .ee_mod_meta_dados(dados, yi, vi, estudo, subgrupo)
  d <- md$d; d[[subgrupo]] <- droplevels(as.factor(d[[subgrupo]]))
  transf <- .ee_mod_meta_transf(medida); rot <- .ee_mod_meta_rot(medida, idioma)
  teste <- if (hartung_knapp) "knha" else "z"
  ajusta <- function(dd) suppressWarnings(metafor::rma(yi = dd[[yi]], vi = dd[[vi]], method = metodo, test = teste, level = 100 * nivel))
  glob <- ajusta(d)
  linhas <- list(); est <- c(); ep <- c()
  for (g in levels(d[[subgrupo]])) {
    dg <- d[d[[subgrupo]] == g, , drop = FALSE]
    r <- tryCatch(ajusta(dg), error = function(e) NULL)
    if (is.null(r)) next
    l <- .ee_mod_meta_linha(r, g, transf)
    l$Q <- r$QE; l$i2_q <- .ee_mod_i2q(r$QE, r$k - 1)
    linhas[[g]] <- l; est[g] <- as.numeric(r$b); ep[g] <- r$se
  }
  lg <- .ee_mod_meta_linha(glob, .ee_mod_t(idioma, "Global", "Overall"), transf)
  lg$Q <- glob$QE; lg$i2_q <- .ee_mod_i2q(glob$QE, glob$k - 1)
  tab <- rbind(do.call(rbind, linhas), lg); rownames(tab) <- NULL
  # teste entre subgrupos com τ² separado por subgrupo (padrão do pacote meta e do curso)
  entre <- metafor::rma(yi = est, sei = ep, method = "EE")
  # alternativa com τ² comum: meta-regressão com o subgrupo como moderador
  mr <- metafor::rma(yi = d[[yi]], vi = d[[vi]], mods = ~ d[[subgrupo]], method = metodo, test = teste)
  pres <- data.frame(checagem = c("Diferença entre subgrupos (τ² separado)", "Diferença entre subgrupos (τ² comum, meta-regressão)"),
                     estatistica = c(entre$QE, mr$QM), gl = c(length(est) - 1, mr$QMdf[1]), p = c(entre$QEp, mr$QMp),
                     leitura = c("Q entre subgrupos sobre as estimativas de cada subgrupo (teste principal)",
                                 if (hartung_knapp) "F de Hartung-Knapp com τ² comum" else "QM com τ² comum (sensibilidade)"))
  avisos <- character()
  ks <- vapply(linhas, function(l) l$k, numeric(1))
  if (any(ks < 5)) avisos <- c(avisos, paste0("Subgrupos com menos de 5 estudos (", paste(names(ks)[ks < 5], collapse = ", "), "): τ² do subgrupo muito impreciso; o teste entre subgrupos tem pouco poder."))
  if (hartung_knapp) avisos <- c(avisos, "Com Hartung-Knapp o EP de um subgrupo pode ficar menor que o usual (fator de ajuste < 1, sem o ajuste ad hoc), e o Q entre subgrupos fica mais liberal, não mais conservador (o pacote meta faz o mesmo). Conferir com hartung_knapp = FALSE.")
  avisos <- c(avisos, "Subgrupos são comparações observacionais entre estudos (características dos estudos podem se confundir): não provam causalidade; hipóteses de subgrupo devem estar no protocolo.")
  lim <- medida == "COR"
  frases <- vapply(linhas, function(l) paste0(l$modelo, " (k = ", l$k, "): ",
                                              .ee_efeito_txt(rot, if (is.null(transf)) l$estimativa else l$efeito, if (is.null(transf)) l$ic_inf else l$efeito_ic_inf,
                                                             if (is.null(transf)) l$ic_sup else l$efeito_ic_sup, limitado = lim, casas = 2, decimal = decimal,
                                                             idioma = idioma, nivel = nivel),
                                              .ee_sep(decimal), "I² = ", .pr_num(l$i2_q, 1, TRUE, decimal), "%"), character(1))
  apa <- c(paste0(.ee_mod_t(idioma, "Subgrupos por ", "Subgroups by "), subgrupo, ": ", paste(frases, collapse = "; "), "."),
           paste0(.ee_mod_t(idioma, "Teste de diferença entre subgrupos: ", "Test for subgroup differences: "), .ee_estat("Q", entre$QE, length(est) - 1, decimal),
                  .ee_sep(decimal), .pr_p(entre$QEp, decimal), "."))
  efeito <- data.frame(indice = paste0(rot, ": ", tab$modelo),
                       valor = if (is.null(transf)) tab$estimativa else tab$efeito,
                       ic_inf = if (is.null(transf)) tab$ic_inf else tab$efeito_ic_inf,
                       ic_sup = if (is.null(transf)) tab$ic_sup else tab$efeito_ic_sup,
                       principal = tab$modelo != .ee_mod_t(idioma, "Global", "Overall"))
  .ee_resultado("Meta-análise por subgrupos", n = nrow(d), n_excluidos = md$n_excluidos, tabela = tab, efeito = efeito,
                pressupostos = pres, decisoes = c(
                  paste0("Efeitos aleatórios (", metodo, ") dentro de cada subgrupo, com τ² próprio de cada um (padrão do pacote meta usado no curso); I² de cada subgrupo por (Q - gl)/Q."),
                  "A evidência de modificação de efeito vem do teste entre subgrupos, não de um subgrupo 'dar significativo' e o outro não.",
                  "Com τ² comum (meta-regressão com o subgrupo como moderador) o teste sai ao lado como sensibilidade."),
                avisos = avisos, apa = apa, extra = list(modelos = linhas, teste_entre = entre, meta_regressao = mr))
}

meta_regressao <- function(dados, moderadores, yi = "yi", vi = "vi", medida = NULL, metodo = "REML", hartung_knapp = FALSE,
                           fatores = NULL, permutacao = FALSE, semente = 2026, reamostras = 1000, nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("metafor")
  if (is.null(medida)) medida <- if (!is.null(attr(dados, "medida"))) attr(dados, "medida") else "generica"
  vars <- .ee_mod_vars(moderadores)
  .ee_mod_nome_ok(vars)
  md <- .ee_mod_meta_dados(dados, yi, vi, NULL, vars)
  d <- .ee_mod_fatores(md$d, vars, fatores = fatores)
  fm <- stats::as.formula(paste("~", paste(moderadores, collapse = " + ")))
  teste <- if (hartung_knapp) "knha" else "z"
  res <- metafor::rma(yi = d[[yi]], vi = d[[vi]], mods = fm, data = d, method = metodo, test = teste, level = 100 * nivel)
  cf <- data.frame(termo = rownames(res$beta), B = as.numeric(res$beta), ep = res$se, estatistica = res$zval, p = res$pval,
                   ic_inf = res$ci.lb, ic_sup = res$ci.ub, row.names = NULL)
  razao <- medida %in% c("OR", "RR")
  if (razao) { cf$exp_B <- exp(cf$B); cf$mudanca_pct <- 100 * (exp(cf$B) - 1) }
  perm <- NULL
  if (permutacao) {
    set.seed(semente)
    perm <- metafor::permutest(res, iter = reamostras, progbar = FALSE)
    cf$p_permutacao <- perm$pval
  }
  k <- res$k; ncoef <- nrow(cf) - 1
  pres <- data.frame(checagem = c("QM (moderadores)", "QE (heterogeneidade residual)", "τ² residual", "I² residual (metafor)", "R² (heterogeneidade explicada)"),
                     estatistica = c(res$QM, res$QE, res$tau2, res$I2, res$R2), gl = c(res$QMdf[1], res$k - res$p, NA, NA, NA),
                     p = c(res$QMp, res$QEp, NA, NA, NA),
                     leitura = c(if (hartung_knapp) "F de Hartung-Knapp" else "teste conjunto dos moderadores",
                                 "heterogeneidade que sobra depois dos moderadores", "", "",
                                 "pseudo-R² baseado em τ²: instável com poucos estudos; pode cair ao acrescentar moderador (comparar modelos por LRT em ML ou AIC)"))
  avisos <- .ee_mod_aviso_numericos(d, vars)
  if (k / max(ncoef, 1) < 10) avisos <- c(avisos, paste0(k, " estudos para ", ncoef, " coeficiente(s) de moderador: menos de 10 estudos por covariável (regra da Cochrane), poder baixo e risco de achado espúrio."))
  avisos <- c(avisos, "Associação entre estudos (ecológica): não autoriza conclusão causal nem sobre indivíduos.")
  cft <- cf[-1, , drop = FALSE]
  cs <- if (nrow(cft) && max(abs(cft$B)) < .1) 4 else 2  # moderador em unidade pequena (ex.: latitude em graus)
  frase_coef <- vapply(seq_len(nrow(cft)), function(i) paste0(cft$termo[i], ": B = ", .pr_num(cft$B[i], cs, TRUE, decimal), .ee_sep(decimal),
                                                            "EP = ", .pr_num(cft$ep[i], cs, TRUE, decimal), .ee_sep(decimal),
                                                            .ee_ic(cft$ic_inf[i], cft$ic_sup[i], cs, TRUE, decimal, nivel, idioma), .ee_sep(decimal), .pr_p(cft$p[i], decimal),
                                                            if (razao) paste0(" (", .pr_num(cft$mudanca_pct[i], 1, TRUE, decimal), .ee_mod_t(idioma, "% na razão por unidade)", "% change in the ratio per unit)")) else ""), character(1))
  apa <- c(paste0(.ee_mod_t(idioma, "Meta-regressão (", "Meta-regression ("), metodo, if (hartung_knapp) " + Hartung-Knapp" else "", "; k = ", k, "): ",
                  if (hartung_knapp) .ee_estat("F", res$QM, c(res$QMdf[1], res$QMdf[2]), decimal) else .ee_estat("QM", res$QM, res$QMdf[1], decimal),
                  .ee_sep(decimal), .pr_p(res$QMp, decimal), .ee_sep(decimal), "R² = ", .pr_num(res$R2 / 100, 3, FALSE, decimal), "; ",
                  paste(frase_coef, collapse = "; "), "."),
           paste0(.ee_mod_t(idioma, "Heterogeneidade residual: ", "Residual heterogeneity: "), .ee_estat("QE", res$QE, res$k - res$p, decimal), .ee_sep(decimal),
                  .pr_p(res$QEp, decimal), .ee_sep(decimal), "τ² = ", .pr_num(res$tau2, 3, TRUE, decimal), .ee_sep(decimal), "I² = ", .pr_num(res$I2, 1, TRUE, decimal), "%."))
  efeito <- data.frame(indice = if (razao) paste0("exp(B): ", cft$termo) else paste0("B: ", cft$termo),
                       valor = if (razao) cft$exp_B else cft$B, ic_inf = if (razao) exp(cft$ic_inf) else cft$ic_inf,
                       ic_sup = if (razao) exp(cft$ic_sup) else cft$ic_sup, principal = TRUE)
  .ee_resultado("Meta-regressão", n = k, n_excluidos = md$n_excluidos, tabela = cf, efeito = efeito, pressupostos = pres,
                decisoes = c(paste0("Meta-regressão de efeitos mistos (", metodo, "); moderador contínuo entra contínuo (categorizar perde informação); fator usa o primeiro nível como referência."),
                             "Hartung-Knapp recomendado em meta-regressão (hartung_knapp = TRUE); teste de permutação opcional (permutacao = TRUE).",
                             "I² residual e R² como no metafor (os números do curso vêm do metafor via meta::metareg); o R² é descritivo e não serve para escolher modelo."),
                avisos = avisos, apa = apa, extra = list(modelo = res, permutacao = perm))
}

# ================================================================ meta-análise de proporção (v1.4)
# Curso (meta-análise de prevalência): prevalência = eventos / n, IC de cada estudo pelo
# método exato de Clopper-Pearson, efeito comum ou aleatório pelo inverso da variância,
# com transformação antes de combinar. transformacao: "PLO" (logito, padrão), "PFT"
# (duplo arco-seno de Freeman-Tukey, volta pela média harmônica dos n) ou "PR" (bruta).
# O logito é o padrão porque a volta do Freeman-Tukey pode distorcer a estimativa
# combinada (Schwarzer et al., 2019); com muitos estudos de prevalência extrema, o GLMM
# (meta::metaprop, method = "GLMM") é a alternativa. Os evals conferem com meta::metaprop.
meta_proporcao <- function(dados, eventos, n, estudo = NULL, transformacao = c("PLO", "PFT", "PR"), metodo = "REML",
                           nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("metafor")
  transformacao <- match.arg(transformacao)
  d <- as.data.frame(dados)
  for (v in c(eventos, n)) if (!v %in% names(d)) stop("coluna ausente: ", v)
  d <- d[stats::complete.cases(d[, c(eventos, n)]), , drop = FALSE]
  x <- d[[eventos]]; ni <- d[[n]]
  if (any(x < 0 | x > ni)) stop("eventos fora de 0 a n em algum estudo")
  k <- nrow(d); if (k < 2) stop("meta-análise precisa de pelo menos 2 estudos")
  lab <- if (!is.null(estudo)) as.character(d[[estudo]]) else paste0("estudo ", seq_len(k))
  es <- metafor::escalc(measure = transformacao, xi = x, ni = ni)
  res <- metafor::rma(yi = es$yi, vi = es$vi, method = metodo, level = 100 * nivel)
  volta <- switch(transformacao, PLO = function(v) metafor::transf.ilogit(v),
                  PFT = function(v) metafor::transf.ipft.hm(v, targs = list(ni = ni)), PR = function(v) v)
  pr <- stats::predict(res, level = 100 * nivel)
  comb <- data.frame(modelo = if (metodo %in% c("EE", "FE", "CE")) "efeito comum" else paste0("efeitos aleatórios (", metodo, ")"),
                     prevalencia = volta(pr$pred), ic_inf = volta(pr$ci.lb), ic_sup = volta(pr$ci.ub),
                     pi_inf = if (!is.null(pr$pi.lb)) volta(pr$pi.lb) else NA_real_, pi_sup = if (!is.null(pr$pi.ub)) volta(pr$pi.ub) else NA_real_)
  cp <- t(mapply(function(a, b) stats::binom.test(a, b, conf.level = nivel)$conf.int, x, ni))
  estudos <- data.frame(estudo = lab, eventos = x, n = ni, prevalencia = x / ni, ic_inf = cp[, 1], ic_sup = cp[, 2],
                        peso = stats::weights(res), stringsAsFactors = FALSE)
  het <- data.frame(Q = res$QE, gl = k - 1, p = res$QEp, I2 = res$I2, tau2 = res$tau2, tau = sqrt(res$tau2))
  avisos <- character()
  if (transformacao == "PFT") avisos <- c(avisos, "Freeman-Tukey: a volta para proporção pode distorcer a estimativa combinada (Schwarzer et al., 2019); conferir com o logito.")
  if (any(x == 0 | x == ni) && transformacao == "PLO") avisos <- c(avisos, "Estudo com 0 ou n eventos: o logito usa correção de 0,5 nesses estudos (escalc); conferir com GLMM (meta::metaprop).")
  if (k < 10) avisos <- c(avisos, "Menos de 10 estudos: τ² e o intervalo de predição ficam imprecisos.")
  f <- function(v, kk = 1) .pr_num(100 * v, kk, TRUE, decimal)
  apa <- if (identical(idioma, "en"))
    paste0("The pooled prevalence (", comb$modelo[1], ", ", c(PLO = "logit", PFT = "Freeman-Tukey double arcsine", PR = "untransformed")[[transformacao]],
           " scale; k = ", k, ") was ", f(comb$prevalencia), "% (", round(100 * nivel), "% CI [", f(comb$ic_inf), ", ", f(comb$ic_sup), "]",
           if (!is.na(comb$pi_inf)) paste0("; prediction interval [", f(comb$pi_inf), ", ", f(comb$pi_sup), "]") else "", "), I² = ", f(res$I2 / 100), "%.") else
    paste0("A prevalência combinada (", comb$modelo[1], ", escala ", c(PLO = "logito", PFT = "duplo arco-seno de Freeman-Tukey", PR = "bruta")[[transformacao]],
           "; k = ", k, ") foi de ", f(comb$prevalencia), "% (IC ", round(100 * nivel), "% [", f(comb$ic_inf), "; ", f(comb$ic_sup), "]",
           if (!is.na(comb$pi_inf)) paste0("; intervalo de predição [", f(comb$pi_inf), "; ", f(comb$pi_sup), "]") else "", "), I² = ", f(res$I2 / 100), "%.")
  .ee_resultado("Meta-análise de prevalência", n = sum(ni), tabela = comb,
                efeito = data.frame(indice = "prevalência combinada", valor = comb$prevalencia, ic_inf = comb$ic_inf, ic_sup = comb$ic_sup, principal = TRUE),
                decisoes = c(paste0("Transformação ", transformacao, " (metafor::escalc), ", metodo, ", volta para proporção",
                                    if (transformacao == "PFT") " pela média harmônica dos n (transf.ipft.hm)" else "", "."),
                             "IC de cada estudo pelo método exato de Clopper-Pearson (binom.test), como no curso."),
                avisos = avisos, apa = apa, extra = list(modelo = res, tabelas_extra = list(`Estudos` = estudos, `Heterogeneidade` = het)))
}

# ================================================================ meta-análise de acurácia diagnóstica (v1.4)
# Modelo bivariado de efeitos aleatórios de Reitsma et al. (2005): logito da sensibilidade
# e logito da especificidade de cada estudo, variâncias dentro do estudo pela binomial
# (1/VP + 1/FN e 1/VN + 1/FP) e matriz de efeitos aleatórios não estruturada (REML,
# metafor::rma.mv). Mesmo modelo do mada::reitsma, que parametriza pela taxa de falso
# positivo (a correlação troca de sinal). Correção de 0,5 em todas as células quando algum
# estudo tem célula zero (padrão do mada). RV+, RV- e RDO saem das estimativas combinadas.
# Viés de publicação: teste de Deeks et al. (2005), regressão do ln RDO em 1/raiz(ESS)
# ponderada por ESS; p < .10 indica assimetria.
meta_acuracia <- function(dados, vp, fp, fn, vn, estudo = NULL, correcao = 0.5, nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("metafor")
  d <- as.data.frame(dados)
  for (v in c(vp, fp, fn, vn)) if (!v %in% names(d)) stop("coluna ausente: ", v)
  d <- d[stats::complete.cases(d[, c(vp, fp, fn, vn)]), , drop = FALSE]
  k <- nrow(d); if (k < 4) stop("o modelo bivariado pede pelo menos 4 estudos")
  tp <- d[[vp]]; fpos <- d[[fp]]; fneg <- d[[fn]]; tn <- d[[vn]]
  lab <- if (!is.null(estudo)) as.character(d[[estudo]]) else paste0("estudo ", seq_len(k))
  corrigido <- any(c(tp, fpos, fneg, tn) == 0)
  if (corrigido) { tp <- tp + correcao; fpos <- fpos + correcao; fneg <- fneg + correcao; tn <- tn + correcao }
  lse <- stats::qlogis(tp / (tp + fneg)); lsp <- stats::qlogis(tn / (tn + fpos))
  v_se <- 1 / tp + 1 / fneg; v_sp <- 1 / tn + 1 / fpos
  long <- data.frame(estudo = rep(seq_len(k), each = 2), medida = factor(rep(c("sens", "espec"), k), levels = c("sens", "espec")),
                     yi = as.vector(rbind(lse, lsp)), vi = as.vector(rbind(v_se, v_sp)))
  res <- metafor::rma.mv(yi, V = diag(long$vi), mods = ~ medida - 1, random = ~ medida | estudo, struct = "UN", data = long,
                         method = "REML", level = 100 * nivel)
  b <- stats::coef(res); V <- stats::vcov(res); z <- stats::qnorm(1 - (1 - nivel) / 2)
  il <- stats::plogis
  se_c <- il(b[1]); sp_c <- il(b[2])
  ic <- function(i) il(b[i] + c(-1, 1) * z * sqrt(V[i, i]))
  # RV e RDO pelo método delta na escala log
  lrp <- log(se_c / (1 - sp_c)); lrn <- log((1 - se_c) / sp_c); ldor <- b[1] + b[2]
  g_lrp <- c(1 - se_c, sp_c); g_lrn <- c(-se_c, -(1 - sp_c)); g_dor <- c(1, 1)
  ep <- function(g) sqrt(as.numeric(t(g) %*% V %*% g))
  comb <- data.frame(indice = c("sensibilidade", "especificidade", "RV+", "RV-", "razão de chances diagnóstica"),
                     valor = c(se_c, sp_c, exp(lrp), exp(lrn), exp(ldor)),
                     ic_inf = c(ic(1)[1], ic(2)[1], exp(lrp - z * ep(g_lrp)), exp(lrn - z * ep(g_lrn)), exp(ldor - z * ep(g_dor))),
                     ic_sup = c(ic(1)[2], ic(2)[2], exp(lrp + z * ep(g_lrp)), exp(lrn + z * ep(g_lrn)), exp(ldor + z * ep(g_dor))),
                     principal = c(TRUE, TRUE, FALSE, FALSE, FALSE), stringsAsFactors = FALSE)
  het <- data.frame(tau2_logito_sens = res$tau2[1], tau2_logito_espec = res$tau2[2], correlacao = res$rho)
  ess <- 4 * (tp + fneg) * (tn + fpos) / (tp + fneg + tn + fpos)
  ldor_i <- log((tp * tn) / (fpos * fneg))
  dk <- summary(stats::lm(ldor_i ~ I(1 / sqrt(ess)), weights = ess))$coefficients
  deeks <- data.frame(inclinacao = dk[2, 1], ep = dk[2, 2], t = dk[2, 3], gl = k - 2, p = dk[2, 4])
  estudos <- data.frame(estudo = lab, sensibilidade = il(lse), especificidade = il(lsp), n_doentes = tp + fneg, n_nao_doentes = tn + fpos,
                        stringsAsFactors = FALSE)
  avisos <- character()
  if (corrigido) avisos <- c(avisos, paste0("Célula zero em algum estudo: ", correcao, " somado a todas as células de todos os estudos (padrão do mada::reitsma)."))
  if (deeks$p < .10) avisos <- c(avisos, "Teste de Deeks com p < .10: assimetria do funil (possível viés de publicação ou efeito de estudos pequenos).")
  if (k < 10) avisos <- c(avisos, "Menos de 10 estudos: o teste de Deeks tem poder baixo e o τ² fica impreciso.")
  f <- function(x) .pr_num(x, 2, FALSE, decimal)
  apa <- if (identical(idioma, "en"))
    paste0("A bivariate random-effects model (Reitsma et al., 2005; k = ", k, ") yielded a pooled sensitivity of ", f(se_c), " (", round(100 * nivel), "% CI [",
           f(ic(1)[1]), ", ", f(ic(1)[2]), "]) and specificity of ", f(sp_c), " (", round(100 * nivel), "% CI [", f(ic(2)[1]), ", ", f(ic(2)[2]), "]).") else
    paste0("O modelo bivariado de efeitos aleatórios (Reitsma et al., 2005; k = ", k, ") estimou sensibilidade combinada de ", f(se_c), " (IC ", round(100 * nivel),
           "% [", f(ic(1)[1]), "; ", f(ic(1)[2]), "]) e especificidade de ", f(sp_c), " (IC ", round(100 * nivel), "% [", f(ic(2)[1]), "; ", f(ic(2)[2]), "]).")
  .ee_resultado("Meta-análise de acurácia diagnóstica (bivariada)", n = sum(tp + fpos + fneg + tn), tabela = comb, efeito = comb[1:2, ],
                decisoes = c("Bivariado de Reitsma (REML, metafor::rma.mv): logitos de sensibilidade e especificidade com efeitos aleatórios correlacionados.",
                             "RV+, RV- e razão de chances diagnóstica derivadas das estimativas combinadas (IC pelo método delta).",
                             "Deeks et al. (2005): ln RDO em 1/raiz(ESS), ponderado por ESS; p < .10 sugere assimetria."),
                avisos = avisos, apa = apa, extra = list(modelo = res, tabelas_extra = list(`Estudos` = estudos, `Heterogeneidade` = het, `Teste de Deeks` = deeks)))
}

# ---------------------------------------------------------------- v1.5: meta-análise em rede
# meta_rede(): meta-análise em rede frequentista (netmeta; Rücker, 2012), a do curso
# "Metanálise de Rede" (Franchini 2012). Entrada por braço, no formato longo: uma linha
# por braço de estudo, com estudo, tratamento e n + média + dp (contínuo) ou n + eventos
# (binário). Saem: efeitos contra a referência no modelo comum e no aleatório, τ² e I²,
# decomposição do Q (dentro dos desenhos = heterogeneidade; entre desenhos =
# inconsistência), separação de nós (evidência direta contra indireta, netsplit) e
# P-scores (Rücker & Schwarzer, 2015). As aulas usam só o modelo comum; a skill põe os
# dois lado a lado e decide pela heterogeneidade declarada.
meta_rede <- function(dados, estudo, tratamento, n, media = NULL, dp = NULL, eventos = NULL, medida = NULL,
                      referencia = NULL, modelo = c("aleatorio", "comum"), menor_melhor = TRUE, nivel = .95,
                      decimal = ".", idioma = "pt") {
  .pr_exige("netmeta")
  modelo <- match.arg(modelo)
  continuo <- !is.null(media)
  if (continuo && is.null(dp)) stop("desfecho contínuo pede media e dp")
  if (!continuo && is.null(eventos)) stop("informe media e dp (contínuo) ou eventos (binário)")
  if (is.null(medida)) medida <- if (continuo) "MD" else "OR"
  vars <- c(estudo, tratamento, n, media, dp, eventos)
  falta <- setdiff(vars, names(dados)); if (length(falta)) stop("colunas ausentes: ", paste(falta, collapse = ", "))
  d <- dados[stats::complete.cases(dados[vars]), vars]
  args <- list(treat = d[[tratamento]], n = d[[n]], studlab = d[[estudo]], sm = medida)
  if (continuo) { args$mean <- d[[media]]; args$sd <- d[[dp]] } else args$event <- d[[eventos]]
  # pairwise() saiu do netmeta 3 para o pacote meta
  pw_f <- if ("pairwise" %in% getNamespaceExports("netmeta")) getExportedValue("netmeta", "pairwise") else getExportedValue("meta", "pairwise")
  pw <- suppressWarnings(do.call(pw_f, args))
  trts <- sort(unique(as.character(d[[tratamento]])))
  if (is.null(referencia)) referencia <- trts[1]
  net <- suppressWarnings(netmeta::netmeta(pw, reference.group = referencia, common = TRUE, random = TRUE, level.ma = nivel))
  bt <- !continuo && medida %in% c("OR", "RR")
  tr <- function(x) if (bt) exp(x) else x
  outros <- setdiff(net$trts, referencia)
  lin <- function(sfx) data.frame(modelo = if (sfx == "common") "comum" else "aleatório", tratamento = outros,
                                  efeito = tr(net[[paste0("TE.", sfx)]][outros, referencia]),
                                  ic_inf = tr(net[[paste0("lower.", sfx)]][outros, referencia]),
                                  ic_sup = tr(net[[paste0("upper.", sfx)]][outros, referencia]),
                                  p = net[[paste0("pval.", sfx)]][outros, referencia], stringsAsFactors = FALSE, row.names = NULL)
  tab <- rbind(lin("common"), lin("random"))
  dd <- tryCatch(netmeta::decomp.design(net)$Q.decomp, error = function(e) NULL)
  if (is.null(dd) || nrow(dd) != 3) dd <- data.frame(Q = c(net$Q, NA, NA), df = c(net$df.Q, NA, NA), pval = c(net$pval.Q, NA, NA))
  qd <- data.frame(componente = c("total", "dentro dos desenhos (heterogeneidade)", "entre desenhos (inconsistência)"),
                   Q = dd$Q, gl = dd$df, p = dd$pval, stringsAsFactors = FALSE)
  ns <- tryCatch(netmeta::netsplit(net), error = function(e) NULL)
  sfx <- if (modelo == "comum") "common" else "random"
  so <- FALSE
  if (!is.null(ns)) { dir <- ns[[paste0("direct.", sfx)]]; ind <- ns[[paste0("indirect.", sfx)]]; cmp <- ns[[paste0("compare.", sfx)]]
    so <- !is.na(dir$TE) & !is.na(ind$TE) }
  split <- if (any(so)) data.frame(comparacao = dir$comparison[so], direto = tr(dir$TE[so]), indireto = tr(ind$TE[so]),
                                   diferenca = cmp$TE[so], p_inconsistencia = cmp$p[so], stringsAsFactors = FALSE, row.names = NULL)
  rk <- netmeta::netrank(net, small.values = if (menor_melhor) "desirable" else "undesirable")
  ps <- rk[[if (modelo == "comum") "ranking.common" else "ranking.random"]]
  rank <- data.frame(tratamento = names(ps), p_score = unname(ps), stringsAsFactors = FALSE)
  rank <- rank[order(-rank$p_score), ]; rownames(rank) <- NULL
  het <- data.frame(tau2 = net$tau2, tau = net$tau, I2 = net$I2, I2_inf = net$lower.I2, I2_sup = net$upper.I2,
                    k_estudos = net$k, n_tratamentos = net$n, n_desenhos = net$d)
  ef <- tab[tab$modelo == if (modelo == "comum") "comum" else "aleatório", ]
  nome_med <- switch(medida, MD = "DM", SMD = "DMP", OR = "RC", RR = "RR", medida)
  efeito <- data.frame(indice = paste0(nome_med, " ", ef$tratamento, " vs ", referencia), valor = ef$efeito, ic_inf = ef$ic_inf, ic_sup = ef$ic_sup,
                       principal = TRUE, stringsAsFactors = FALSE)
  avisos <- character()
  if (isTRUE(qd$p[3] < .10)) avisos <- c(avisos, "Inconsistência entre desenhos (Q entre desenhos, p < .10): a transitividade está em dúvida; ler o netsplit antes do ranking.")
  if (modelo == "comum" && isTRUE(net$I2 > .25)) avisos <- c(avisos, "Modelo comum com I² acima de 25%: o modelo aleatório é a escolha defensável.")
  en <- identical(idioma, "en")
  fr <- paste(vapply(seq_len(nrow(ef)), function(i) paste0(ef$tratamento[i], ": ", nome_med, " = ", .pr_num(ef$efeito[i], 2, TRUE, decimal), ", ",
                                                          .ee_ic(ef$ic_inf[i], ef$ic_sup[i], 2, TRUE, decimal, nivel, idioma)), ""), collapse = .ee_sep(decimal))
  apa <- if (en) paste0("A frequentist network meta-analysis (netmeta, ", if (modelo == "comum") "common" else "random", "-effects model; ", net$k, " studies, ", net$n,
                        " treatments) compared each treatment with ", referencia, ": ", fr, ". Heterogeneity: τ² = ", .pr_num(net$tau2, 3, TRUE, decimal), ", I² = ",
                        .pr_num(100 * net$I2, 1, TRUE, decimal), "%; inconsistency between designs: Q(", qd$gl[3], ") = ", .pr_num(qd$Q[3], 2, TRUE, decimal), ", ",
                        .pr_p(qd$p[3], decimal), ".") else
    paste0("A meta-análise em rede frequentista (netmeta, modelo ", if (modelo == "comum") "de efeito comum" else "de efeitos aleatórios", "; ", net$k, " estudos, ",
           net$n, " tratamentos) comparou cada tratamento com ", referencia, ": ", fr, ". Heterogeneidade: τ² = ", .pr_num(net$tau2, 3, TRUE, decimal), ", I² = ",
           .pr_num(100 * net$I2, 1, TRUE, decimal), "%; inconsistência entre desenhos: Q(", qd$gl[3], ") = ", .pr_num(qd$Q[3], 2, TRUE, decimal), ", ",
           .pr_p(qd$p[3], decimal), ".")
  extras <- list(`Heterogeneidade` = het, `Decomposição do Q` = qd, `Ranking (P-score)` = rank)
  if (!is.null(split)) extras$`Direto contra indireto (netsplit)` <- split
  .ee_resultado(paste0("Meta-análise em rede (", nome_med, ", referência ", referencia, ")"), n = net$k, tabela = tab, efeito = efeito,
                decisoes = c(paste0("Modelo relatado: ", if (modelo == "comum") "efeito comum (o das aulas)" else "efeitos aleatórios", "; o outro fica na tabela."),
                             "Inconsistência pela decomposição do Q (entre desenhos) e pela separação de nós (netsplit) nas comparações com evidência direta e indireta.",
                             paste0("P-score com ", if (menor_melhor) "valores menores" else "valores maiores", " como melhores; ranking é descritivo e não substitui os IC."),
                             if (bt) paste0(nome_med, " exponenciada; diferenças do netsplit na escala log.") else character()),
                avisos = avisos, apa = apa, extra = list(tabelas_extra = extras, modelo = net))
}

# ---------------------------------------------------------------- v1.5: extração de dados
# extrair_media_dp(): média e DP para meta-análise a partir do que o estudo relata, como
# na aula de extração de dados (Cochrane Handbook, cap. 6):
#   "mediana_iqr"            mediana, q1, q3        média de Luo et al. (2018) ou Wan et al. (2014); DP de Wan et al. (2014)
#   "mediana_amplitude"      mediana, mínimo, máximo
#   "mediana_iqr_amplitude"  os cinco números
#   "ep"                     DP = EP × √n
#   "ic"                     DP = √n (ls − li) / (2 t) da média (t com n − 1 gl; n ≥ 60 usa z, como o Handbook)
#   "combinar"               dois grupos num só (fórmula do Handbook, tabela 6.5.a)
#   "mudanca"                DP da mudança pré-pós com a correlação r informada
#   "correlacao"             r implícita num estudo que relata DP pré, pós e da mudança
# As conversões de mediana chamam as rotinas do pacote meta (que usam a tabela 2 de Wan
# et al., 2014, para n até 201); a média de Luo é o padrão do meta e da Cochrane.
extrair_media_dp <- function(tipo = c("mediana_iqr", "mediana_amplitude", "mediana_iqr_amplitude", "ep", "ic", "combinar", "mudanca", "correlacao"),
                             n = NULL, mediana = NULL, q1 = NULL, q3 = NULL, minimo = NULL, maximo = NULL, media = NULL, ep = NULL,
                             ic_inf = NULL, ic_sup = NULL, nivel = .95, n2 = NULL, media2 = NULL, dp1 = NULL, dp2 = NULL,
                             dp_pre = NULL, dp_pos = NULL, dp_mudanca = NULL, r = NULL, metodo_media = c("Luo", "Wan")) {
  tipo <- match.arg(tipo); metodo_media <- match.arg(metodo_media)
  amb <- environment()
  pede <- function(...) { nm <- c(...); f <- nm[vapply(nm, function(x) is.null(get(x, envir = amb)), TRUE)]
    if (length(f)) stop("para tipo = '", tipo, "' faltou: ", paste(f, collapse = ", "), call. = FALSE) }
  if (tipo %in% c("mediana_iqr", "mediana_amplitude", "mediana_iqr_amplitude")) .pr_exige("meta")
  out <- switch(tipo,
    mediana_iqr = { pede("n", "mediana", "q1", "q3")
      x <- getFromNamespace("mean_sd_iqr", "meta")(n = n, median = mediana, q1 = q1, q3 = q3, method.mean = metodo_media)
      data.frame(media = x$mean, dp = x$sd, metodo = paste0("média de ", metodo_media, if (metodo_media == "Luo") " et al. (2018)" else " et al. (2014)",
                                                            "; DP de Wan et al. (2014) pelo IIQ")) },
    mediana_amplitude = { pede("n", "mediana", "minimo", "maximo")
      x <- getFromNamespace("mean_sd_range", "meta")(n = n, median = mediana, min = minimo, max = maximo, method.mean = metodo_media)
      data.frame(media = x$mean, dp = x$sd, metodo = paste0("média de ", metodo_media, "; DP de Wan et al. (2014) pela amplitude")) },
    mediana_iqr_amplitude = { pede("n", "mediana", "q1", "q3", "minimo", "maximo")
      x <- getFromNamespace("mean_sd_iqr_range", "meta")(n = n, median = mediana, q1 = q1, q3 = q3, min = minimo, max = maximo, method.mean = metodo_media)
      data.frame(media = x$mean, dp = x$sd, metodo = paste0("média de ", metodo_media, "; DP de Wan et al. (2014) pelos cinco números")) },
    ep = { pede("n", "ep"); data.frame(media = if (is.null(media)) NA_real_ else media, dp = ep * sqrt(n), metodo = "DP = EP × √n") },
    ic = { pede("n", "ic_inf", "ic_sup")
      q <- ifelse(n < 60, stats::qt(1 - (1 - nivel) / 2, n - 1), stats::qnorm(1 - (1 - nivel) / 2))
      data.frame(media = if (is.null(media)) (ic_inf + ic_sup) / 2 else media, dp = sqrt(n) * (ic_sup - ic_inf) / (2 * q),
                 metodo = ifelse(n < 60, "DP pelo IC da média com t (n < 60)", "DP pelo IC da média com z")) },
    combinar = { pede("n", "media", "dp1", "n2", "media2", "dp2")
      nt <- n + n2; m <- (n * media + n2 * media2) / nt
      s <- sqrt(((n - 1) * dp1^2 + (n2 - 1) * dp2^2 + n * n2 / nt * (media^2 + media2^2 - 2 * media * media2)) / (nt - 1))
      data.frame(n = nt, media = m, dp = s, metodo = "grupos combinados (Cochrane Handbook, tabela 6.5.a)") },
    mudanca = { pede("dp_pre", "dp_pos", "r")
      data.frame(dp = sqrt(dp_pre^2 + dp_pos^2 - 2 * r * dp_pre * dp_pos), r = r, metodo = "DP da mudança com r imputada (Handbook 6.5.2.8): declarar a origem de r e testar outros valores") },
    correlacao = { pede("dp_pre", "dp_pos", "dp_mudanca")
      data.frame(r = (dp_pre^2 + dp_pos^2 - dp_mudanca^2) / (2 * dp_pre * dp_pos), metodo = "r implícita no estudo que relata os três DP (Handbook 6.5.2.8)") })
  rownames(out) <- NULL
  out
}

# ================================================================ desfecho ordinal e nominal (v1.6)
# Martinková e Hladká (2023, cap. 6.7) tratam a regressão ordinal (logito cumulativo) e a
# multinomial (logito de categoria-base) como modelos de item; aqui servem a qualquer
# desfecho com 3 ou mais categorias. Parametrização do MASS::polr e do PLUM do SPSS:
# logit P(Y <= j) = limiar_j - (b1 X1 + ...), de modo que B > 0 empurra para categorias
# mais altas e OR = exp(B) é a razão das chances de estar acima de qualquer corte.
# Teste das linhas paralelas como no PLUM: razão de verossimilhança entre o modelo de
# chances proporcionais e o modelo em que cada preditor tem um coeficiente por corte
# (ordinal::clm com nominal = ~ preditores); por termo, ordinal::nominal_test.

.ee_mod_pseudo_r2 <- function(ll1, ll0, n) {
  cs <- 1 - exp(2 * (ll0 - ll1) / n)
  data.frame(indice = c("Cox e Snell", "Nagelkerke", "McFadden"),
             valor = c(cs, cs / (1 - exp(2 * ll0 / n)), 1 - ll1 / ll0))
}

.ee_mod_y_ordenado <- function(y, vd, max_cat = 10) {
  nota <- character()
  if (is.ordered(y)) return(list(y = droplevels(y), nota = nota))
  if (is.factor(y)) {
    nota <- paste0(vd, " é fator sem ordem declarada: a ordem usada é a dos níveis (", paste(levels(droplevels(y)), collapse = " < "),
                   "). Conferir; se estiver errada, declarar com factor(..., levels = ..., ordered = TRUE).")
    return(list(y = factor(y, levels = levels(droplevels(y)), ordered = TRUE), nota = nota))
  }
  if (is.character(y) || is.logical(y)) stop("desfecho ordinal em texto: declarar a ordem com factor(", vd, ", levels = c(...), ordered = TRUE)")
  u <- sort(unique(y))
  if (length(u) > max_cat) stop(vd, " tem ", length(u), " valores distintos: com tantas categorias, tratar como contínuo (modelo_glm) ou agrupar")
  list(y = factor(y, levels = u, ordered = TRUE), nota = nota)
}

modelo_ordinal <- function(dados, vd, preditores, ligacao = c("logit", "probit"), referencia = NULL, fatores = NULL,
                           nivel = .95, decimal = ".", idioma = "pt") {
  ligacao <- match.arg(ligacao)
  .pr_exige("MASS")
  vars <- .ee_mod_vars(preditores)
  .ee_mod_nome_ok(c(vd, vars))
  dd <- .ee_dados(dados, unique(c(vd, vars)))
  d <- .ee_mod_fatores(dd$d, vars, referencia, fatores)
  decisoes <- character(); avisos <- .ee_mod_aviso_numericos(d, vars)
  yo <- .ee_mod_y_ordenado(d[[vd]], vd); d[[vd]] <- yo$y; if (length(yo$nota)) avisos <- c(avisos, yo$nota)
  k <- nlevels(d[[vd]])
  if (k < 3) stop("desfecho com ", k, " categorias: com 2, usar regressao_logistica() ou modelo_glm(familia = \"binomial\")")
  freq <- table(d[[vd]])
  if (any(freq < .ee_n_min_print)) avisos <- c(avisos, paste0("Categoria do desfecho com menos de ", .ee_n_min_print, " casos (",
    paste(names(freq)[freq < .ee_n_min_print], collapse = ", "), "): limiar mal estimado; considerar agrupar categorias vizinhas."))
  metodo <- if (ligacao == "logit") "logistic" else "probit"
  f <- .ee_mod_formula(vd, preditores); environment(f) <- environment()
  ajusta <- function(fo) suppressWarnings(MASS::polr(fo, data = d, Hess = TRUE, method = metodo))
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op), add = TRUE)
  fs <- ajusta(f)
  options(contrasts = c("contr.treatment", "contr.poly"))
  ft <- ajusta(f)
  f0 <- stats::as.formula(paste(vd, "~ 1"), env = environment()); f00 <- ajusta(f0)
  n <- nrow(d)
  ll1 <- as.numeric(stats::logLik(ft)); ll0 <- as.numeric(stats::logLik(f00))
  p_beta <- length(stats::coef(ft))
  lr <- 2 * (ll1 - ll0); p_lr <- stats::pchisq(lr, p_beta, lower.tail = FALSE)
  r2 <- .ee_mod_pseudo_r2(ll1, ll0, n)
  # coeficientes (contraste de tratamento) e limiares
  Vt <- stats::vcov(ft); bt <- stats::coef(ft)
  coefs <- .ee_mod_coefs(bt, Vt[names(bt), names(bt), drop = FALSE], nivel, if (ligacao == "logit") "OR" else NA)
  zeta <- ft$zeta; ep_z <- sqrt(diag(Vt)[names(zeta)]); zc <- stats::qnorm(1 - (1 - nivel) / 2)
  limiares <- data.frame(corte = names(zeta), limiar = unname(zeta), ep = unname(ep_z),
                         ic_inf = unname(zeta - zc * ep_z), ic_sup = unname(zeta + zc * ep_z))
  # Wald tipo III por termo (contrastes soma-zero)
  mm <- stats::model.matrix(stats::terms(fs), d)
  asg <- attr(mm, "assign")[-1]
  rotulos <- attr(stats::terms(fs), "term.labels")
  bs <- stats::coef(fs); Vs <- stats::vcov(fs)[names(bs), names(bs), drop = FALSE]
  wald <- .ee_mod_wald3(bs, Vs, asg, rotulos)
  # linhas paralelas
  pres <- data.frame(checagem = "Convergência", estatistica = NA, gl = NA, p = NA,
                     leitura = if (isTRUE(ft$convergence == 0)) "convergiu" else "NÃO convergiu")
  paralelas <- NULL
  if (requireNamespace("ordinal", quietly = TRUE)) {
    lk <- if (ligacao == "logit") "logit" else "probit"
    c_po <- tryCatch(ordinal::clm(f, data = d, link = lk), error = function(e) NULL)
    c_np <- tryCatch(suppressWarnings(ordinal::clm(f0, nominal = stats::as.formula(paste("~", paste(preditores, collapse = " + "))),
                                                    data = d, link = lk)), error = function(e) NULL)
    if (!is.null(c_po) && !is.null(c_np) && is.finite(as.numeric(stats::logLik(c_np)))) {
      lrp <- 2 * (as.numeric(stats::logLik(c_np)) - as.numeric(stats::logLik(c_po)))
      glp <- attr(stats::logLik(c_np), "df") - attr(stats::logLik(c_po), "df")
      pp <- stats::pchisq(max(lrp, 0), glp, lower.tail = FALSE)
      pres <- rbind(pres, data.frame(checagem = "Teste das linhas paralelas (RV)", estatistica = lrp, gl = glp, p = pp,
                                     leitura = if (pp < .05) "chances não proporcionais em algum preditor" else "compatível com chances proporcionais"))
      ll_po <- as.numeric(stats::logLik(c_po)); gl_po <- attr(stats::logLik(c_po), "df")
      paralelas <- do.call(rbind, lapply(vars, function(v) {
        m <- tryCatch(suppressWarnings(ordinal::clm(f, nominal = stats::as.formula(paste("~", v)), data = d, link = lk)), error = function(e) NULL)
        if (is.null(m) || !is.finite(as.numeric(stats::logLik(m)))) return(data.frame(termo = v, qui2_rv = NA_real_, gl = NA_real_, p = NA_real_))
        x2 <- max(2 * (as.numeric(stats::logLik(m)) - ll_po), 0); g <- attr(stats::logLik(m), "df") - gl_po
        data.frame(termo = v, qui2_rv = x2, gl = g, p = stats::pchisq(x2, g, lower.tail = FALSE))
      }))
      if (pp < .05) avisos <- c(avisos, paste0("Linhas paralelas rejeitadas (χ²(", glp, ") = ", round(lrp, 2), ", p ", sub("^p ", "", .pr_p(pp)),
        "): o OR único resume mal o preditor", if (!is.null(paralelas) && any(paralelas$p < .05, na.rm = TRUE))
          paste0(" (", paste(paralelas$termo[paralelas$p < .05], collapse = ", "), ")") else "",
        ". Com n grande o teste acusa desvios pequenos: comparar os OR por corte (regressões logísticas binárias em cada corte) antes de trocar de modelo; alternativa: modelo_multinomial() ou chances proporcionais parciais (ordinal::clm com nominal =)."))
      decisoes <- c(decisoes, "Teste das linhas paralelas como no PLUM do SPSS: razão de verossimilhança contra o modelo com um coeficiente por corte (ordinal::clm, nominal =); por variável, a mesma razão com nominal = ~ variável (como ordinal::nominal_test).")
    } else avisos <- c(avisos, "O modelo de chances não proporcionais não convergiu (comum com células vazias): teste das linhas paralelas indisponível.")
  } else avisos <- c(avisos, "Pacote ordinal ausente: teste das linhas paralelas não rodou (install.packages(\"ordinal\")).")
  cel <- character()
  fat <- vars[vapply(vars, function(v) is.factor(d[[v]]), logical(1))]
  for (v in fat) if (any(table(d[[v]], d[[vd]]) == 0)) cel <- c(cel, v)
  if (length(cel)) avisos <- c(avisos, paste0("Cruzamento do desfecho com ", paste(cel, collapse = ", "), " tem célula vazia: EP e limiares podem explodir; agrupar categorias."))
  sem <- coefs
  efeito <- if (ligacao == "logit")
    data.frame(indice = paste0("OR: ", sem$termo), valor = sem$exp_B, ic_inf = sem$exp_ic_inf, ic_sup = sem$exp_ic_sup, principal = TRUE) else
    data.frame(indice = paste0("B (probit): ", sem$termo), valor = sem$B, ic_inf = sem$ic_inf, ic_sup = sem$ic_sup, principal = TRUE)
  efeito <- rbind(efeito, data.frame(indice = paste0("R² de ", r2$indice), valor = r2$valor, ic_inf = NA, ic_sup = NA, principal = FALSE))
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  nome <- if (en) paste0("Ordinal ", if (ligacao == "logit") "logistic" else "probit", " regression (proportional odds)") else
    paste0("Regressão ", if (ligacao == "logit") "logística" else "probit", " ordinal (chances proporcionais)")
  termos_txt <- vapply(seq_len(nrow(sem)), function(i) {
    if (ligacao == "logit") paste0(sem$termo[i], ": ", .ee_efeito_txt("OR", sem$exp_B[i], sem$exp_ic_inf[i], sem$exp_ic_sup[i], FALSE, 2, decimal, idioma, nivel), sep, .pr_p(sem$p[i], decimal)) else
      paste0(sem$termo[i], ": ", .ee_efeito_txt("B", sem$B[i], sem$ic_inf[i], sem$ic_sup[i], FALSE, 2, decimal, idioma, nivel), sep, .pr_p(sem$p[i], decimal))
  }, character(1))
  apa <- paste0(nome, if (en) " of " else " de ", vd, " (", k, if (en) " categories" else " categorias", "): ",
                .ee_estat("χ²", lr, p_beta, decimal), sep, .pr_p(p_lr, decimal), sep, "R² Nagelkerke = ", .pr_num(r2$valor[2], 2, FALSE, decimal),
                ". ", paste(termos_txt, collapse = "; "), ".")
  decisoes <- c(decisoes,
                paste0("MASS::polr (", metodo, "), parametrização do PLUM do SPSS: B > 0 desloca para categorias mais altas", if (ligacao == "logit") "; OR = exp(B), com IC de Wald." else "; B na escala probit (desvios-padrão da variável latente)."),
                "Teste global: razão de verossimilhança contra o modelo só com os limiares; R² de Cox e Snell, Nagelkerke e McFadden (pseudo-R², não é variância explicada).",
                "Testes dos termos: Wald χ² Tipo III com contrastes soma-zero.")
  .ee_mod_classe(.ee_resultado(nome, n = n, n_excluidos = dd$n_excluidos, tabela = coefs, efeito = efeito, pressupostos = pres,
                decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(tabelas_extra = c(list("Limiares (cortes)" = limiares, "Wald tipo III por termo" = wald),
                                               if (!is.null(paralelas)) list("Linhas paralelas por termo" = paralelas)),
                             modelo = ft, coeficientes = coefs, limiares = limiares, wald = wald, pseudo_r2 = r2, rv_global = c(qui2 = lr, gl = p_beta, p = p_lr),
                             paralelas = paralelas, ligacao = ligacao)))
}

# Multinomial (logito de categoria-base, nnet::multinom): um conjunto de coeficientes por
# categoria contra a referência; RRR = exp(B). Testes dos termos por razão de
# verossimilhança (modelo sem o termo), como os "Likelihood Ratio Tests" do NOMREG do SPSS.
modelo_multinomial <- function(dados, vd, preditores, referencia_vd = NULL, referencia = NULL, fatores = NULL,
                               nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("nnet")
  vars <- .ee_mod_vars(preditores)
  .ee_mod_nome_ok(c(vd, vars))
  dd <- .ee_dados(dados, unique(c(vd, vars)))
  d <- .ee_mod_fatores(dd$d, vars, referencia, fatores)
  decisoes <- character(); avisos <- .ee_mod_aviso_numericos(d, vars)
  y <- droplevels(factor(d[[vd]], ordered = FALSE))
  if (!is.null(referencia_vd)) {
    if (!referencia_vd %in% levels(y)) stop("a categoria '", referencia_vd, "' não existe em ", vd)
    y <- stats::relevel(y, ref = as.character(referencia_vd))
  }
  d[[vd]] <- y; k <- nlevels(y)
  if (k < 3) stop("desfecho com ", k, " categorias: com 2, usar regressao_logistica()")
  freq <- table(y)
  if (any(freq < .ee_n_min_print)) avisos <- c(avisos, paste0("Categoria do desfecho com menos de ", .ee_n_min_print, " casos (",
    paste(names(freq)[freq < .ee_n_min_print], collapse = ", "), "): coeficientes instáveis; agrupar ou excluir a categoria."))
  f <- .ee_mod_formula(vd, preditores); environment(f) <- environment()
  ajusta <- function(fo) nnet::multinom(fo, data = d, trace = FALSE, maxit = 1000, Hess = TRUE)
  op <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(op), add = TRUE)
  fit <- ajusta(f); fit0 <- ajusta(stats::as.formula(paste(vd, "~ 1"), env = environment()))
  if (isTRUE(fit$convergence != 0)) avisos <- c(avisos, "O modelo não convergiu em 1000 iterações: resultados não confiáveis.")
  n <- nrow(d)
  ll1 <- as.numeric(stats::logLik(fit)); ll0 <- as.numeric(stats::logLik(fit0))
  gl_g <- length(stats::coef(fit)) - length(stats::coef(fit0))
  lr <- 2 * (ll1 - ll0); p_lr <- stats::pchisq(lr, gl_g, lower.tail = FALSE)
  r2 <- .ee_mod_pseudo_r2(ll1, ll0, n)
  B <- stats::coef(fit); if (is.null(dim(B))) B <- matrix(B, nrow = 1, dimnames = list(levels(y)[2], names(B)))
  V <- stats::vcov(fit); zc <- stats::qnorm(1 - (1 - nivel) / 2)
  coefs <- do.call(rbind, lapply(rownames(B), function(cat) {
    b <- B[cat, ]; nm <- paste0(cat, ":", names(b)); ep <- sqrt(diag(V)[nm])
    data.frame(categoria = cat, termo = names(b), B = unname(b), ep = unname(ep), wald_qui2 = unname((b / ep)^2), gl = 1,
               p = unname(2 * stats::pnorm(-abs(b / ep))), RRR = unname(exp(b)), ic_inf = unname(exp(b - zc * ep)),
               ic_sup = unname(exp(b + zc * ep)), row.names = NULL)
  }))
  # RV por termo (modelo sem o termo)
  rotulos <- attr(stats::terms(f), "term.labels")
  rv <- do.call(rbind, lapply(rotulos, function(tt) {
    f_sem <- stats::update(f, paste(". ~ . -", tt))
    m_sem <- tryCatch(ajusta(f_sem), error = function(e) NULL)
    if (is.null(m_sem)) return(data.frame(termo = tt, qui2_rv = NA, gl = NA, p = NA))
    x2 <- 2 * (ll1 - as.numeric(stats::logLik(m_sem))); g <- length(stats::coef(fit)) - length(stats::coef(m_sem))
    data.frame(termo = tt, qui2_rv = x2, gl = g, p = stats::pchisq(x2, g, lower.tail = FALSE))
  }))
  pred <- stats::predict(fit, type = "class")
  acerto <- mean(as.character(pred) == as.character(y)); base <- max(freq) / n
  pres <- data.frame(checagem = c("Convergência", "Acerto de classificação"), estatistica = c(fit$convergence, acerto), gl = NA, p = NA,
                     leitura = c(if (fit$convergence == 0) "convergiu" else "NÃO convergiu",
                                 paste0("contra ", round(100 * base, 1), "% de acerto chutando a categoria mais comum")))
  fat <- vars[vapply(vars, function(v) is.factor(d[[v]]), logical(1))]
  cel <- character(); for (v in fat) if (any(table(d[[v]], y) == 0)) cel <- c(cel, v)
  if (length(cel)) avisos <- c(avisos, paste0("Cruzamento do desfecho com ", paste(cel, collapse = ", "), " tem célula vazia: separação, RRR infinito ou EP enorme; agrupar categorias."))
  if (any(coefs$ep > 10, na.rm = TRUE)) avisos <- c(avisos, "Há EP > 10 em algum coeficiente: sinal de separação; não interpretar esses RRR.")
  sem <- coefs[coefs$termo != "(Intercept)", , drop = FALSE]
  efeito <- rbind(data.frame(indice = paste0("RRR ", sem$categoria, " vs ", levels(y)[1], ": ", sem$termo), valor = sem$RRR,
                             ic_inf = sem$ic_inf, ic_sup = sem$ic_sup, principal = TRUE),
                  data.frame(indice = paste0("R² de ", r2$indice), valor = r2$valor, ic_inf = NA, ic_sup = NA, principal = FALSE))
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  nome <- if (en) "Multinomial logistic regression" else "Regressão logística multinomial"
  txt <- vapply(seq_len(nrow(sem)), function(i) paste0(sem$categoria[i], " vs ", levels(y)[1], ", ", sem$termo[i], ": ",
    .ee_efeito_txt("RRR", sem$RRR[i], sem$ic_inf[i], sem$ic_sup[i], FALSE, 2, decimal, idioma, nivel), sep, .pr_p(sem$p[i], decimal)), character(1))
  rv_txt <- vapply(seq_len(nrow(rv)), function(i) paste0(rv$termo[i], ", ", .ee_estat("χ²", rv$qui2_rv[i], rv$gl[i], decimal), sep, .pr_p(rv$p[i], decimal)), character(1))
  apa <- c(paste0(nome, if (en) " of " else " de ", vd, " (", if (en) "reference: " else "referência: ", levels(y)[1], "): ",
                  .ee_estat("χ²", lr, gl_g, decimal), sep, .pr_p(p_lr, decimal), sep, "R² Nagelkerke = ", .pr_num(r2$valor[2], 2, FALSE, decimal),
                  if (en) ". Likelihood-ratio tests: " else ". Testes de razão de verossimilhança: ", paste(rv_txt, collapse = "; "), "."),
           paste0(paste(txt, collapse = "; "), "."))
  decisoes <- c(decisoes, paste0("nnet::multinom, categoria de referência do desfecho: ", levels(y)[1], " (argumento referencia_vd); RRR = exp(B), razão de riscos relativos contra a referência, IC de Wald."),
                "Testes dos termos por razão de verossimilhança (modelo sem o termo), como no NOMREG do SPSS; R² de Cox e Snell, Nagelkerke e McFadden.",
                "Com desfecho ordenado, começar por modelo_ordinal(): a multinomial ignora a ordem e gasta k - 1 coeficientes por preditor.")
  .ee_mod_classe(.ee_resultado(nome, n = n, n_excluidos = dd$n_excluidos, tabela = rv, efeito = efeito, pressupostos = pres,
                decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(tabelas_extra = list("Coeficientes por categoria" = coefs), modelo = fit, coeficientes = coefs,
                             pseudo_r2 = r2, rv_global = c(qui2 = lr, gl = gl_g, p = p_lr), referencia_vd = levels(y)[1])))
}

# Efeito contextual (v1.6): com X_cwc (desvio da média do cluster) e X_media_cgm (média
# do cluster centrada) no mesmo modelo, b_cwc é o efeito intragrupo e b_media o efeito entre
# grupos; o efeito contextual, diferença entre os dois, diz quanto a média do grupo acrescenta
# além da posição da pessoa dentro dele (Enders & Tofighi, 2007; aula de regressão multinível,
# "efeito contextual"). Com X centrado na média geral e a média do cluster no modelo, o
# coeficiente da média do cluster já é o efeito contextual, sem conta.
efeito_contextual <- function(resultado, var, nivel = .95, decimal = ".", idioma = "pt") {
  fit <- if (inherits(resultado, "merMod")) resultado else resultado$modelo
  if (is.null(fit) || !inherits(fit, "merMod")) stop("passe o resultado de multinivel_linear()/multinivel_logistico() ou um modelo lme4", call. = FALSE)
  b <- lme4::fixef(fit); V <- as.matrix(stats::vcov(fit))
  w <- paste0(var, "_cwc"); m <- paste0(var, "_media_cgm")
  if (!all(c(w, m) %in% names(b))) stop("o modelo precisa ter ", w, " e ", m, " (criados por centralizar_multinivel(..., tipo = \"grupo\") ou \"geral\")", call. = FALSE)
  est <- b[[m]] - b[[w]]; ep <- sqrt(V[m, m] + V[w, w] - 2 * V[m, w]); z <- est / ep
  zc <- stats::qnorm(1 - (1 - nivel) / 2); p <- 2 * stats::pnorm(-abs(z))
  tab <- data.frame(efeito = c("intragrupo (b_cwc)", "entre grupos (b_media)", "contextual (b_media - b_cwc)"),
                    B = c(b[[w]], b[[m]], est), ep = c(sqrt(V[w, w]), sqrt(V[m, m]), ep))
  tab$ic_inf <- tab$B - zc * tab$ep; tab$ic_sup <- tab$B + zc * tab$ep; tab$z <- tab$B / tab$ep; tab$p <- 2 * stats::pnorm(-abs(tab$z))
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  apa <- if (en) paste0("The contextual effect of ", var, " (between-cluster minus within-cluster slope) was B = ", .pr_num(est, 2, TRUE, decimal), sep,
                        .ee_ic(tab$ic_inf[3], tab$ic_sup[3], 2, TRUE, decimal, nivel, idioma), sep, .pr_p(p, decimal), ".") else
    paste0("O efeito contextual de ", var, " (inclinação entre grupos menos intragrupo) foi B = ", .pr_num(est, 2, TRUE, decimal), sep,
           .ee_ic(tab$ic_inf[3], tab$ic_sup[3], 2, TRUE, decimal, nivel, idioma), sep, .pr_p(p, decimal), ".")
  .ee_resultado("Efeito contextual (multinível)", n = stats::nobs(fit), tabela = tab,
                efeito = data.frame(indice = "efeito contextual", valor = est, ic_inf = tab$ic_inf[3], ic_sup = tab$ic_sup[3], principal = TRUE),
                decisoes = c("Efeito contextual = b_media - b_cwc com EP pela matriz de covariância dos efeitos fixos (Wald, z).",
                             "No modelo logístico, a diferença está na escala do logito; exp(B) é a razão de chances associada a uma unidade da média do grupo com a posição individual fixa."),
                apa = apa)
}
