# est_naoparametricos.R: testes não paramétricos para delineamentos complexos (v1.5).
# Carregar por carregar.R.
#
#   np_fatorial()     dois ou mais fatores entre sujeitos, por postos: estatística do tipo
#                     ANOVA (ATS) e do tipo Wald (WTS) do rankFD, efeitos relativos com IC e
#                     comparações múltiplas do nparcomp (mctp)
#   np_longitudinal() medidas repetidas com ou sem fatores entre sujeitos (F1-LD-F1,
#                     F2-LD-F1, LD-F1): ATS e WTS do nparLD, efeitos relativos (RTE) e
#                     comparações entre tempos por grupo (mctp.rm)
#   manova_np()       MANOVA sem normalidade multivariada nem homogeneidade das matrizes:
#                     WTS e MATS com bootstrap paramétrico do MANOVA.RM (MANOVA.wide para
#                     fatores entre sujeitos, multRM para medidas repetidas), mais os
#                     testes univariados com correção de Holm
#
# Rota do curso "Testes não paramétricos para delineamentos complexos" (rankFD, nparLD,
# nparcomp, MANOVA.RM). O efeito relativo p = P(X < Y) + 0,5 P(X = Y) compara a
# distribuição de uma célula com a média das distribuições do delineamento: 0,5 é ausência
# de efeito (Brunner, Konietschke, Pauly & Puri, 2017; Noguchi et al., 2012). A ATS é a
# estatística recomendada para amostras pequenas e moderadas; a WTS pede amostras grandes
# e fica ao lado. O nparLD 2.3 trocou a interface (sem os argumentos description e
# plot.CI das aulas, que davam erro) e não devolve mais IC para o RTE: o IC sai do mctp.rm,
# por grupo.

.np_formula <- function(vd, termos, op = "*") stats::as.formula(paste(.ee_assoc_bq(vd), "~", paste(.ee_assoc_bq(termos), collapse = paste0(" ", op, " "))))

.np_tab_stat <- function(m, tipo) {
  m <- as.matrix(m)
  data.frame(efeito = rownames(m), estatistica = tipo, valor = unname(m[, 1]),
             gl1 = unname(m[, 2]), gl2 = if (ncol(m) >= 4) unname(m[, 3]) else NA_real_,
             p = unname(m[, ncol(m)]), stringsAsFactors = FALSE, row.names = NULL)
}

.np_frase_ats <- function(tab, decimal, idioma) {
  a <- tab[tab$estatistica == "ATS", ]
  en <- identical(idioma, "en")
  paste(vapply(seq_len(nrow(a)), function(i) paste0(a$efeito[i], ": ATS(", .pr_num(a$gl1[i], 2, TRUE, decimal),
    if (!is.na(a$gl2[i]) && is.finite(a$gl2[i])) paste0(", ", .pr_num(a$gl2[i], 2, TRUE, decimal)) else "", ") = ",
    .pr_num(a$valor[i], 2, TRUE, decimal), ", ", .pr_p(a$p[i], decimal)), ""), collapse = .ee_sep(decimal))
}

# Comparações múltiplas do nparcomp: diferenças de efeitos relativos com IC simultâneo e p
# ajustado pelo próprio procedimento (Konietschke et al., 2012). As aulas usam
# asy.method = "log.odds", cujo estimador sai numa escala transformada (3 - 1 = 1,01 quando
# a diferença dos efeitos relativos é 0,39); a skill usa "fisher", que devolve a diferença na
# escala original, com o IC preso a [-1, 1].
.np_mctp <- function(d, vd, fator, tipo = "Tukey", rm = FALSE, id = NULL, rotulo = NULL, nivel = .95) {
  f <- .np_formula(vd, fator)
  r <- if (!rm) suppressWarnings(nparcomp::mctp(f, data = d, type = tipo, asy.method = "fisher",
                                               conf.level = nivel, info = FALSE)) else
    suppressWarnings(nparcomp::mctp.rm(f, data = d, type = tipo, asy.method = "fisher", conf.level = nivel, info = FALSE))
  a <- r$Analysis
  out <- data.frame(subamostra = if (is.null(rotulo)) "todos" else rotulo, comparacao = rownames(a),
                    diferenca_efeito = a$Estimator, ic_inf = a$Lower, ic_sup = a$Upper, estatistica = a$Statistic,
                    p_ajustado = a$p.Value, stringsAsFactors = FALSE, row.names = NULL)
  di <- r$Data.Info
  ef <- data.frame(subamostra = out$subamostra[1], nivel = di$Sample, n = di$Size, efeito_relativo = di$Effect,
                   ic_inf = di$Lower, ic_sup = di$Upper, stringsAsFactors = FALSE, row.names = NULL)
  list(comparacoes = out, efeitos = ef)
}

np_fatorial <- function(dados, vd, fatores, posthoc = TRUE, tipo_posthoc = "Tukey", nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("rankFD", "nparcomp")
  if (length(fatores) < 1) stop("indique ao menos um fator entre sujeitos")
  dados <- .ee_assoc_vazios(dados, c(vd, fatores))
  dd <- .ee_dados(dados, c(vd, fatores), numericas = vd)
  d <- dd$d
  for (f in fatores) d[[f]] <- droplevels(factor(d[[f]]))
  cel <- table(interaction(d[fatores], drop = FALSE))
  if (any(cel < 2)) stop("célula com menos de 2 casos: ", paste(names(cel)[cel < 2], collapse = ", "))
  fit <- rankFD::rankFD(.np_formula(vd, fatores), data = d, hypothesis = "H0F", CI.method = "logit", alpha = 1 - nivel)
  tab <- rbind(.np_tab_stat(fit$ANOVA.Type.Statistic, "ATS"), .np_tab_stat(fit$Wald.Type.Statistic, "WTS"))
  de <- fit$Descriptive
  efeito <- data.frame(celula = apply(de[fatores], 1, paste, collapse = ":"), n = de$Size, indice = "efeito relativo",
                       valor = de$Rel.Effect, ic_inf = de$Lower, ic_sup = de$Upper, principal = TRUE,
                       stringsAsFactors = FALSE, row.names = NULL)
  ph <- NULL; efs_ph <- NULL; dec_ph <- character()
  ats <- tab[tab$estatistica == "ATS", ]
  if (posthoc) {
    alfa <- 1 - nivel
    inter <- if (length(fatores) == 2) ats$p[ats$efeito == paste(fatores, collapse = ":")] else NA
    if (length(fatores) == 1 && ats$p[1] < alfa) {
      r <- .np_mctp(d, vd, fatores, tipo_posthoc, nivel = nivel); ph <- r$comparacoes; efs_ph <- r$efeitos
      dec_ph <- paste0("Comparações múltiplas (mctp, ", tipo_posthoc, ", IC de Fisher) depois da ATS significativa.")
    } else if (length(fatores) == 2 && isTRUE(inter < alfa)) {
      rr <- lapply(levels(d[[fatores[2]]]), function(l) .np_mctp(d[d[[fatores[2]]] == l, ], vd, fatores[1], tipo_posthoc,
                                                                  rotulo = paste0(fatores[2], " = ", l), nivel = nivel))
      ph <- do.call(rbind, lapply(rr, `[[`, "comparacoes")); efs_ph <- do.call(rbind, lapply(rr, `[[`, "efeitos"))
      dec_ph <- paste0("Interação significativa: comparações de ", fatores[1], " dentro de cada nível de ", fatores[2],
                       " (mctp por subamostra, como nas aulas); cada subamostra controla o próprio erro tipo I.")
    } else if (length(fatores) >= 2) {
      sig <- ats$efeito[ats$p < alfa & ats$efeito %in% fatores]
      if (length(sig)) {
        rr <- lapply(sig, function(f) .np_mctp(d, vd, f, tipo_posthoc, rotulo = paste0("efeito principal de ", f), nivel = nivel))
        ph <- do.call(rbind, lapply(rr, `[[`, "comparacoes")); efs_ph <- do.call(rbind, lapply(rr, `[[`, "efeitos"))
        dec_ph <- "Sem interação significativa: comparações só nos efeitos principais significativos (mctp com os dados agregados no fator)."
      }
    }
  }
  en <- identical(idioma, "en")
  apa <- if (en) paste0("A rank-based nonparametric factorial analysis (ANOVA-type statistic, rankFD; Brunner et al., 2017) was run for ", vd, ": ",
                        .np_frase_ats(tab, decimal, idioma), ".") else
    paste0("Foi feita uma análise fatorial não paramétrica por postos (estatística do tipo ANOVA, rankFD; Brunner et al., 2017) para ", vd, ": ",
           .np_frase_ats(tab, decimal, idioma), ".")
  .ee_resultado(paste0("Fatorial não paramétrica (rankFD): ", vd, " ~ ", paste(fatores, collapse = " × ")), n = nrow(d),
                n_excluidos = dd$n_excluidos, tabela = tab, efeito = efeito, posthoc = ph,
                decisoes = c("ATS (gl de Box) como teste principal; WTS ao lado, que pede amostra grande e infla o erro tipo I com n pequeno (Brunner et al., 2017).",
                             "Efeito relativo: 0,5 = distribuição da célula igual à média das distribuições; IC logit.",
                             dec_ph),
                avisos = .ee_aviso_nao_sig(ats$p, 1 - nivel, o_que = "efeito", decimal = decimal, rotulos = ats$efeito),
                apa = apa, extra = list(tabelas_extra = if (!is.null(efs_ph)) list(`Efeitos relativos por subamostra` = efs_ph), modelo = fit))
}

np_longitudinal <- function(dados, vd, id, tempo, grupos = NULL, posthoc = TRUE, tipo_posthoc = "Tukey", nivel = .95,
                            decimal = ".", idioma = "pt") {
  .pr_exige("nparLD", "nparcomp")
  if (length(grupos) > 2) stop("o nparLD cobre até dois fatores entre sujeitos com um fator de medida repetida")
  vars <- c(vd, id, tempo, grupos)
  dados <- .ee_assoc_vazios(dados, vars)
  dd <- .ee_dados(dados, vars, numericas = vd)
  d <- dd$d
  d[[tempo]] <- droplevels(factor(d[[tempo]]))
  for (g in grupos) d[[g]] <- droplevels(factor(d[[g]]))
  # Formato longo: uma linha por pessoa e tempo. Pessoa sem todos os tempos sai (a análise
  # do nparLD pede o delineamento completo dentro de cada pessoa).
  cont <- tapply(d[[tempo]], d[[id]], function(x) length(unique(x)))
  incompletos <- names(cont)[cont < nlevels(d[[tempo]])]
  if (any(duplicated(d[c(id, tempo)]))) stop("mais de uma linha para a mesma pessoa no mesmo tempo")
  if (length(incompletos)) d <- d[!d[[id]] %in% incompletos, ]
  d[[id]] <- factor(d[[id]])
  nsuj <- nlevels(droplevels(d[[id]]))
  fit <- nparLD::nparLD(.np_formula(vd, c(tempo, grupos)), data = d, subject = id)
  tab <- rbind(.np_tab_stat(fit$ATS, "ATS"), .np_tab_stat(fit$WTS, "WTS"))
  tab$gl2 <- ifelse(tab$estatistica == "ATS", Inf, NA)
  ef <- fit$effects
  efeito <- data.frame(celula = apply(ef[c(grupos, tempo)], 1, paste, collapse = ":"), n = ef$Nsubj, indice = "RTE",
                       valor = ef$RTE, ic_inf = NA_real_, ic_sup = NA_real_, principal = TRUE,
                       stringsAsFactors = FALSE, row.names = NULL)
  ph <- NULL; efs_ph <- NULL; dec_ph <- character()
  ats <- tab[tab$estatistica == "ATS", ]
  alfa <- 1 - nivel
  if (posthoc) {
    termos_tempo <- ats$efeito[grepl(tempo, ats$efeito, fixed = TRUE)]
    gsig <- if (length(grupos)) grupos[vapply(grupos, function(g) any(ats$p[ats$efeito %in% c(paste(tempo, g, sep = ":"), paste(g, tempo, sep = ":"))] < alfa), TRUE)] else character()
    if (length(gsig)) {
      rr <- unlist(lapply(gsig, function(g) lapply(levels(d[[g]]), function(l)
        .np_mctp(d[d[[g]] == l, ], vd, tempo, tipo_posthoc, rm = TRUE, rotulo = paste0(g, " = ", l), nivel = nivel))), recursive = FALSE)
      ph <- do.call(rbind, lapply(rr, `[[`, "comparacoes")); efs_ph <- do.call(rbind, lapply(rr, `[[`, "efeitos"))
      dec_ph <- paste0("Interação com o tempo significativa: comparações entre tempos (mctp.rm, ", tipo_posthoc,
                       ") dentro de cada nível do grupo que interage, como nas aulas.")
    } else if (ats$p[ats$efeito == tempo] < alfa) {
      r <- .np_mctp(d, vd, tempo, tipo_posthoc, rm = TRUE, nivel = nivel); ph <- r$comparacoes; efs_ph <- r$efeitos
      dec_ph <- paste0("Efeito do tempo significativo: comparações entre tempos (mctp.rm, ", tipo_posthoc, ").")
    }
  }
  desenho <- trimws(gsub("-Design", "", fit$Design))
  en <- identical(idioma, "en")
  apa <- if (en) paste0("A nonparametric analysis of longitudinal data (", desenho, " design, nparLD; Noguchi et al., 2012) was run for ", vd,
                        " (n = ", nsuj, "): ", .np_frase_ats(tab, decimal, idioma), ".") else
    paste0("Foi feita uma análise não paramétrica de dados longitudinais (delineamento ", desenho, ", nparLD; Noguchi et al., 2012) para ", vd,
           " (n = ", nsuj, "): ", .np_frase_ats(tab, decimal, idioma), ".")
  .ee_resultado(paste0("Longitudinal não paramétrica (nparLD, ", desenho, "): ", vd), n = nsuj, n_excluidos = dd$n_excluidos + length(incompletos),
                tabela = tab, efeito = efeito, posthoc = ph,
                decisoes = c(paste0("Delineamento ", desenho, "; ATS (gl1 de Box, gl2 infinito) como teste principal e WTS ao lado."),
                             "RTE: efeito relativo de cada célula contra a média das distribuições do delineamento (0,5 = sem efeito). O nparLD 2.3 não dá IC do RTE; os IC estão nas comparações por grupo (mctp.rm).",
                             if (length(incompletos)) paste0(length(incompletos), " pessoa(s) sem todos os tempos excluída(s).") else character(),
                             dec_ph),
                avisos = .ee_aviso_nao_sig(ats$p, alfa, o_que = "efeito", decimal = decimal, rotulos = ats$efeito),
                apa = apa, extra = list(tabelas_extra = if (!is.null(efs_ph)) list(`Efeitos relativos por subamostra (IC)` = efs_ph), modelo = fit))
}

manova_np <- function(dados, vds, fatores, id = NULL, tempo = NULL, iter = 10000, semente = 2026, univariados = TRUE,
                      nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("MANOVA.RM")
  if (length(vds) < 2) stop("MANOVA pede ao menos duas variáveis dependentes")
  rm <- !is.null(tempo)
  if (rm && is.null(id)) stop("com tempo, indique id (banco no formato longo)")
  vars <- c(vds, fatores, id, tempo)
  dados <- .ee_assoc_vazios(dados, vars)
  dd <- .ee_dados(dados, vars, numericas = vds)
  d <- dd$d
  for (f in c(fatores, tempo)) d[[f]] <- droplevels(factor(d[[f]]))
  incompletos <- character()
  if (rm) {
    cont <- tapply(d[[tempo]], d[[id]], function(x) length(unique(x)))
    incompletos <- names(cont)[cont < nlevels(d[[tempo]])]
    if (length(incompletos)) d <- d[!d[[id]] %in% incompletos, ]
    d <- d[order(d[[id]], d[[tempo]]), ]
  }
  termos <- c(fatores, tempo)
  ajusta <- function(y) {
    f <- stats::as.formula(paste(if (length(y) > 1) paste0("cbind(", paste(.ee_assoc_bq(y), collapse = ", "), ")") else .ee_assoc_bq(y),
                                 "~", paste(.ee_assoc_bq(termos), collapse = " * ")))
    utils::capture.output(m <- if (rm) MANOVA.RM::multRM(f, data = d, subject = id, within = tempo, iter = iter, alpha = 1 - nivel,
                                                          resampling = "paramBS", seed = semente) else
      MANOVA.RM::MANOVA.wide(f, data = d, iter = iter, alpha = 1 - nivel, resampling = "paramBS", seed = semente))
    m
  }
  fit <- ajusta(vds)
  w <- as.matrix(fit$WTS); ma <- as.matrix(fit$MATS); rs <- as.matrix(fit$resampling)
  tab <- data.frame(efeito = rownames(w), WTS = w[, 1], gl = w[, 2], p_WTS_assintotico = w[, 3], MATS = ma[, 1],
                    p_WTS_bootstrap = rs[, 1], p_MATS_bootstrap = rs[, 2], stringsAsFactors = FALSE, row.names = NULL)
  uni <- NULL
  if (univariados) {
    uni <- do.call(rbind, lapply(vds, function(v) { m <- ajusta(v); r <- as.matrix(m$resampling)
      data.frame(vd = v, efeito = rownames(r), p_MATS_bootstrap = r[, 2], stringsAsFactors = FALSE, row.names = NULL) }))
    uni$p_holm <- stats::ave(uni$p_MATS_bootstrap, uni$efeito, FUN = function(p) stats::p.adjust(p, "holm"))
  }
  frase <- paste(vapply(seq_len(nrow(tab)), function(i) paste0(tab$efeito[i], ": MATS = ", .pr_num(tab$MATS[i], 2, TRUE, decimal),
                                                               ", p bootstrap ", sub("^p ", "", .pr_p(tab$p_MATS_bootstrap[i], decimal))), ""),
                 collapse = .ee_sep(decimal))
  en <- identical(idioma, "en")
  apa <- if (en) paste0("A nonparametric ", if (rm) "repeated-measures " else "", "MANOVA (MANOVA.RM, modified ANOVA-type statistic with parametric bootstrap, ",
                        iter, " resamples; Friedrich & Pauly, 2018) was run for ", paste(vds, collapse = ", "), ": ", frase, ".") else
    paste0("Foi feita uma MANOVA ", if (rm) "de medidas repetidas " else "", "sem pressuposto de normalidade (MANOVA.RM, estatística do tipo ANOVA modificada com bootstrap paramétrico, ",
           iter, " reamostras; Friedrich & Pauly, 2018) para ", paste(vds, collapse = ", "), ": ", frase, ".")
  .ee_resultado(paste0("MANOVA robusta (", if (rm) "multRM" else "MANOVA.wide", "): ", paste(vds, collapse = ", ")),
                n = if (rm) length(unique(d[[id]])) else nrow(d), n_excluidos = dd$n_excluidos + length(incompletos), tabela = tab,
                decisoes = c("MATS com bootstrap paramétrico como teste principal: invariante à escala das VDs e válida com matrizes de covariância diferentes entre grupos (Friedrich & Pauly, 2018); WTS ao lado, que pede amostra grande.",
                             paste0("Semente ", semente, " e ", iter, " reamostras."),
                             if (univariados) "Testes univariados (mesma rota, uma VD por vez) com correção de Holm por efeito, para ver qual VD contribui, como nas aulas." else character(),
                             if (length(incompletos)) paste0(length(incompletos), " pessoa(s) sem todos os tempos excluída(s).") else character()),
                avisos = .ee_aviso_nao_sig(tab$p_MATS_bootstrap, 1 - nivel, o_que = "efeito", decimal = decimal, rotulos = tab$efeito),
                apa = apa, extra = list(tabelas_extra = if (!is.null(uni)) list(`Testes univariados (p bootstrap da MATS, Holm)` = uni), modelo = fit))
}
