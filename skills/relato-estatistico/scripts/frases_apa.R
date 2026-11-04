# ==============================================================================
# frases_apa.R  (relato-estatistico v1.5.0)
# Frases de resultado geradas direto dos objetos do R (htest, lm, glm, lavaan,
# lmerMod, coxph, rma, roc, bootnet), no formato que o auditar_texto_statcheck() reconhece.
# Depende de formatador_apa.R. Nenhum numero e digitado: tudo sai do objeto.
# ==============================================================================

.it <- function(s, italico) if (italico) paste0("*", s, "*") else s
.df_txt <- function(df, perfil) {
  if (abs(df - round(df)) < 1e-8) as.character(round(df)) else fmt_apa(df, "livre", perfil, 2)
}

#' Frase APA a partir de objeto do R
#'
#' @param obj htest (t.test, cor.test, chisq.test), lm, glm, lmerMod/glmerMod, coxph, rma ou roc
#' @param perfil perfil do formatador_apa.R
#' @param rotulo_ic "IC" (portugues) ou "CI"
#' @param italico marcar simbolos estatisticos com *...* (markdown)
#' @return vetor de caracteres; para modelos multivariados, nomeado por componente
frases_apa <- function(obj, perfil = "apa7_relato", rotulo_ic = "IC", italico = TRUE, ...) {
  if (inherits(obj, "htest")) return(.frase_htest(obj, perfil, rotulo_ic, italico))
  if (inherits(obj, "glm")) return(.frase_glm(obj, perfil, rotulo_ic, italico))
  if (inherits(obj, "lm")) return(.frase_lm(obj, perfil, rotulo_ic, italico))
  if (inherits(obj, "lavaan")) return(.frase_lavaan(obj, perfil, rotulo_ic, italico))
  if (inherits(obj, c("lmerMod", "glmerMod"))) return(.frase_lmer(obj, perfil, rotulo_ic, italico))
  if (inherits(obj, "coxph")) return(.frase_cox(obj, perfil, rotulo_ic, italico))
  if (inherits(obj, c("rma.uni", "rma"))) return(.frase_metafor(obj, perfil, rotulo_ic, italico))
  if (inherits(obj, "roc")) return(.frase_roc(obj, perfil, rotulo_ic, italico))
  if (inherits(obj, "bootnetResult")) return(frases_rede(obj, perfil = perfil, italico = italico, ...))
  stop("Classe nao suportada: ", paste(class(obj), collapse = "/"), 
       ". Suportadas: htest, lm, glm, lavaan, lmerMod, glmerMod, coxph, rma, roc, bootnetResult.")
}

.frase_htest <- function(h, perfil, rot, it) {
  s <- sep_apa(perfil); m <- tolower(h$method)
  conf <- attr(h$conf.int, "conf.level")
  ic <- function(tipo) if (!is.null(h$conf.int))
    paste0(s, rot, " ", round(100 * conf), "% ", fmt_ic_apa(h$conf.int[1], h$conf.int[2], tipo, perfil)) else ""
  if (grepl("t-test", m)) {
    return(paste0(.it("t", it), "(", .df_txt(h$parameter, perfil), ") = ", fmt_apa(h$statistic, "livre", perfil),
                  s, fmt_p_apa(h$p.value, perfil, italico = it), ic("livre")))
  }
  if (grepl("pearson", m) && grepl("correlation", m)) {
    return(paste0(.it("r", it), "(", .df_txt(h$parameter, perfil), ") = ", fmt_apa(h$estimate, "r", perfil),
                  s, fmt_p_apa(h$p.value, perfil, italico = it), ic("r")))
  }
  if (grepl("chi-squared", m)) {
    n <- if (!is.null(h$observed)) sum(h$observed) else NA
    ntxt <- if (!is.na(n)) paste0(", ", .it("N", it), " = ", n) else ""
    return(paste0("χ²(", .df_txt(h$parameter, perfil), ntxt, ") = ", fmt_apa(h$statistic, "livre", perfil),
                  s, fmt_p_apa(h$p.value, perfil, italico = it)))
  }
  stop("htest nao suportado: ", h$method)
}

.frase_lm <- function(fit, perfil, rot, it) {
  s <- sep_apa(perfil); sm <- summary(fit); f <- sm$fstatistic
  pF <- stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
  modelo <- paste0(.it("F", it), "(", f[2], ", ", f[3], ") = ", fmt_apa(f[1], "livre", perfil), s,
                   fmt_p_apa(pF, perfil, italico = it), s, .it("R", it), "² = ", fmt_apa(sm$r.squared, "R2", perfil),
                   s, .it("R", it), "² ajustado = ", fmt_apa(sm$adj.r.squared, "R2", perfil))
  co <- sm$coefficients; ci <- stats::confint(fit); dfr <- fit$df.residual
  mf <- stats::model.frame(fit); y <- stats::model.response(mf)
  X <- stats::model.matrix(fit)
  out <- c(modelo = modelo)
  for (nm in rownames(co)[-1]) {
    beta_txt <- ""
    if (nm %in% colnames(X) && length(unique(X[, nm])) > 2) {
      b <- co[nm, 1] * stats::sd(X[, nm]) / stats::sd(y)
      beta_txt <- paste0(s, "β = ", fmt_apa(b, "beta", perfil))
    }
    out[nm] <- paste0(.it("B", it), " = ", fmt_apa(co[nm, 1], "livre", perfil), s, .it("EP", it), " = ",
                      fmt_apa(co[nm, 2], "livre", perfil), beta_txt, s, .it("t", it), "(", dfr, ") = ",
                      fmt_apa(co[nm, 3], "livre", perfil), s, fmt_p_apa(co[nm, 4], perfil, italico = it),
                      s, rot, " 95% ", fmt_ic_apa(ci[nm, 1], ci[nm, 2], "livre", perfil))
  }
  out
}

.frase_lavaan <- function(fit, perfil, rot, it) {
  s <- sep_apa(perfil)
  fm <- lavaan::fitMeasures(fit)
  pega <- function(...) { for (n in c(...)) if (n %in% names(fm) && !is.na(fm[[n]])) return(fm[[n]]); NA }
  chi <- pega("chisq.scaled", "chisq"); df <- pega("df.scaled", "df"); p <- pega("pvalue.scaled", "pvalue")
  cfi <- pega("cfi.robust", "cfi.scaled", "cfi"); tli <- pega("tli.robust", "tli.scaled", "tli")
  rm <- pega("rmsea.robust", "rmsea.scaled", "rmsea")
  li <- pega("rmsea.ci.lower.robust", "rmsea.ci.lower.scaled", "rmsea.ci.lower")
  ls <- pega("rmsea.ci.upper.robust", "rmsea.ci.upper.scaled", "rmsea.ci.upper")
  sr <- pega("srmr")
  if (is.na(df) || df == 0) stop("Modelo saturado (gl = 0): indices de ajuste nao informam nada; relatar que o modelo e just-identified.")
  versao <- if (!is.na(pega("cfi.robust"))) "robusta" else if (!is.na(pega("cfi.scaled"))) "escalonada" else "padrao"
  frase <- paste0("χ²(", .df_txt(df, perfil), ") = ", fmt_apa(chi, "livre", perfil), s,
                  fmt_p_apa(p, perfil, italico = it), s, "CFI = ", fmt_apa(cfi, "cfi", perfil), s,
                  "TLI = ", fmt_apa(tli, "tli", perfil), s, "RMSEA = ", fmt_apa(rm, "rmsea", perfil),
                  if (!is.na(li)) paste0(s, rot, " 90% ", fmt_ic_apa(li, ls, "rmsea", perfil)) else "",
                  s, "SRMR = ", fmt_apa(sr, "srmr", perfil))
  attr(frase, "versao_indices") <- versao
  frase
}

.frase_glm <- function(fit, perfil, rot, it) {
  s <- sep_apa(perfil); sm <- summary(fit); co <- sm$coefficients
  ci <- tryCatch(stats::confint.default(fit), error = function(e) NULL)
  fam <- stats::family(fit)$family
  is_bin <- fam == "binomial"
  metric_name <- if (is_bin) "OR" else if (fam == "poisson") "RP" else "exp(B)"
  out <- character()
  for (nm in rownames(co)[-1]) {
    b <- co[nm, 1]; ep <- co[nm, 2]; z_val <- co[nm, 3]; p_val <- co[nm, 4]
    or_val <- exp(b)
    ic_txt <- if (!is.null(ci)) {
      paste0(s, rot, " 95% ", fmt_ic_apa(exp(ci[nm, 1]), exp(ci[nm, 2]), "livre", perfil))
    } else ""
    out[nm] <- paste0(metric_name, " = ", fmt_apa(or_val, "livre", perfil), ic_txt,
                      s, .it("B", it), " = ", fmt_apa(b, "livre", perfil),
                      s, .it("EP", it), " = ", fmt_apa(ep, "livre", perfil),
                      s, .it("z", it), " = ", fmt_apa(z_val, "livre", perfil),
                      s, fmt_p_apa(p_val, perfil, italico = it))
  }
  out
}

.frase_lmer <- function(fit, perfil, rot, it) {
  s <- sep_apa(perfil); sm <- summary(fit); co <- sm$coefficients
  vc <- as.data.frame(lme4::VarCorr(fit))
  tau00 <- vc$vcov[is.na(vc$var2) & vc$grp != "Residual"][1]
  sigma2 <- vc$vcov[vc$grp == "Residual"][1]
  icc <- if (!is.na(tau00) && !is.na(sigma2) && (tau00 + sigma2) > 0) tau00 / (tau00 + sigma2) else NA
  modelo <- paste0("ICC = ", fmt_apa(icc, "r", perfil),
                   s, "τ₀₀ = ", fmt_apa(tau00, "livre", perfil),
                   s, "σ² = ", fmt_apa(sigma2, "livre", perfil))
  out <- c(modelo = modelo)
  dfr <- tryCatch(stats::df.residual(fit), error = function(e) NA)
  for (nm in rownames(co)[-1]) {
    b <- co[nm, "Estimate"]; ep <- co[nm, "Std. Error"]
    stat_col <- grep("value|t value|z value", colnames(co), value = TRUE)[1]
    stat_val <- co[nm, stat_col]
    stat_letter <- if (grepl("z", stat_col)) "z" else "t"
    p_col <- grep("Pr\\(>", colnames(co), value = TRUE)
    p_val <- if (length(p_col)) co[nm, p_col] else 2 * stats::pnorm(abs(stat_val), lower.tail = FALSE)
    ci_inf <- b - 1.96 * ep; ci_sup <- b + 1.96 * ep
    df_str <- if (!is.na(dfr) && stat_letter == "t") paste0("(", .df_txt(dfr, perfil), ")") else ""
    out[nm] <- paste0(.it("B", it), " = ", fmt_apa(b, "livre", perfil),
                      s, .it("EP", it), " = ", fmt_apa(ep, "livre", perfil),
                      s, .it(stat_letter, it), df_str, " = ", fmt_apa(stat_val, "livre", perfil),
                      s, fmt_p_apa(p_val, perfil, italico = it),
                      s, rot, " 95% ", fmt_ic_apa(ci_inf, ci_sup, "livre", perfil))
  }
  out
}

.frase_cox <- function(fit, perfil, rot, it) {
  s <- sep_apa(perfil); sm <- summary(fit); co <- sm$coefficients; ci <- sm$conf.int
  out <- character()
  for (nm in rownames(co)) {
    hr <- co[nm, "exp(coef)"]; b <- co[nm, "coef"]; ep <- co[nm, "se(coef)"]
    z_val <- co[nm, "z"]; p_val <- co[nm, "p"]
    ci_inf <- ci[nm, "lower .95"]; ci_sup <- ci[nm, "upper .95"]
    out[nm] <- paste0("HR = ", fmt_apa(hr, "livre", perfil),
                      s, rot, " 95% ", fmt_ic_apa(ci_inf, ci_sup, "livre", perfil),
                      s, .it("B", it), " = ", fmt_apa(b, "livre", perfil),
                      s, .it("EP", it), " = ", fmt_apa(ep, "livre", perfil),
                      s, .it("z", it), " = ", fmt_apa(z_val, "livre", perfil),
                      s, fmt_p_apa(p_val, perfil, italico = it))
  }
  out
}

.frase_metafor <- function(fit, perfil, rot, it) {
  s <- sep_apa(perfil)
  b <- fit$beta[1, 1]; ep <- fit$se[1]; z_val <- fit$zval[1]; p_val <- fit$pval[1]
  ci_inf <- fit$ci.lb[1]; ci_sup <- fit$ci.ub[1]
  hetero <- paste0("Q(", fit$k - fit$p, ") = ", fmt_apa(fit$QE, "livre", perfil),
                   s, fmt_p_apa(fit$QEp, perfil, italico = it),
                   s, .it("I", it), "² = ", fmt_apa(fit$I2 / 100, "R2", perfil),
                   s, "τ² = ", fmt_apa(fit$tau2, "livre", perfil))
  resumo <- paste0("Tamanho de efeito combinado = ", fmt_apa(b, "livre", perfil),
                   s, rot, " 95% ", fmt_ic_apa(ci_inf, ci_sup, "livre", perfil),
                   s, .it("EP", it), " = ", fmt_apa(ep, "livre", perfil),
                   s, .it("z", it), " = ", fmt_apa(z_val, "livre", perfil),
                   s, fmt_p_apa(p_val, perfil, italico = it),
                   "; Heterogeneidade: ", hetero)
  c(resumo = resumo, heterogeneidade = hetero)
}

.frase_roc <- function(roc_obj, perfil, rot, it) {
  s <- sep_apa(perfil); auc_val <- as.numeric(roc_obj$auc)
  ci_auc <- tryCatch(pROC::ci.auc(roc_obj), error = function(e) NULL)
  ic_txt <- if (!is.null(ci_auc)) {
    paste0(s, rot, " 95% ", fmt_ic_apa(ci_auc[1], ci_auc[3], "cfi", perfil))
  } else ""
  coords_yd <- tryCatch(pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE), error = function(e) NULL)
  corte_txt <- if (!is.null(coords_yd) && nrow(coords_yd) > 0) {
    paste0("; Ponto de corte ótimo (Youden) = ", fmt_apa(coords_yd[1, "threshold"], "livre", perfil),
           " (Sensibilidade = ", fmt_apa(coords_yd[1, "sensitivity"], "cfi", perfil),
           ", Especificidade = ", fmt_apa(coords_yd[1, "specificity"], "cfi", perfil), ")")
  } else ""
  paste0("AUC = ", fmt_apa(auc_val, "cfi", perfil), ic_txt, corte_txt)
}

#' Frases APA de rede psicometrica estimada no bootnet
#'
#' Relata nos, estimador, N, arestas diferentes de zero, densidade e sinais; com cs, a
#' estabilidade da centralidade pelo coeficiente CS (Epskamp, Borsboom & Fried, 2018):
#' nao interpretar abaixo de 0,25 e, de preferencia, acima de 0,50 (cortes que os autores
#' chamam de algo arbitrarios, vindos da simulacao deles).
#' @param rede objeto de bootnet::estimateNetwork()
#' @param cs vetor nomeado de bootnet::corStability() ou NULL
#' @param cor_cs correlacao alvo usada no corStability (padrao do bootnet: 0,7)
frases_rede <- function(rede, cs = NULL, perfil = "apa7_relato", italico = TRUE, cor_cs = 0.7) {
  if (!inherits(rede, "bootnetResult")) stop("frases_rede: use o objeto de bootnet::estimateNetwork().", call. = FALSE)
  W <- as.matrix(rede$graph); p <- ncol(W); w <- W[upper.tri(W)]; nz <- w[abs(w) > 0]
  tun <- rede$arguments$tuning
  est <- if (identical(rede$default, "EBICglasso") && !is.null(tun)) sprintf("EBICglasso (\u03b3 = %s)", fmt_apa(tun, "livre", perfil, 2)) else as.character(rede$default)
  f_rede <- sprintf("A rede de %d n\u00f3s, estimada com %s em %s = %d, teve %d de %d arestas poss\u00edveis diferentes de zero (densidade = %s; %d positivas e %d negativas; peso absoluto m\u00e9dio = %s).",
    p, est, .it("N", italico), as.integer(rede$nPerson), length(nz), length(w), fmt_apa(length(nz) / length(w), "prop", perfil),
    sum(nz > 0), sum(nz < 0), if (length(nz)) fmt_apa(mean(abs(nz)), "r", perfil) else "\u2014")
  if (is.null(cs)) return(c(rede = f_rede))
  cs <- unlist(cs)
  rot <- c(strength = "for\u00e7a", expectedInfluence = "influ\u00eancia esperada", betweenness = "intermedia\u00e7\u00e3o",
           closeness = "proximidade", edge = "arestas")
  nm <- ifelse(names(cs) %in% names(rot), rot[names(cs)], names(cs))
  q1 <- fmt_apa(.25, "livre", perfil); q2 <- fmt_apa(.50, "livre", perfil)
  cls <- ifelse(cs >= .5, "est\u00e1vel", ifelse(cs >= .25, sprintf("abaixo de %s, interpretar com cautela", q2),
                sprintf("abaixo de %s, n\u00e3o interpretar a ordem", q1)))
  partes <- sprintf("%s = %s (%s)", nm, fmt_apa(cs, "livre", perfil), cls)
  f_cs <- sprintf("O coeficiente de estabilidade CS(cor = %s) foi %s; Epskamp et al. (2018) recomendam n\u00e3o interpretar diferen\u00e7as de centralidade com CS abaixo de %s e, de prefer\u00eancia, exigir CS acima de %s.",
    fmt_apa(cor_cs, "livre", perfil, 1), paste(partes, collapse = sep_apa(perfil)), q1, q2)
  c(rede = f_rede, estabilidade = f_cs)
}
