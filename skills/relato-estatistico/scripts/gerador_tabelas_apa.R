# ==============================================================================
# gerador_tabelas_apa.R
# Modulo de Geracao e Formatacao de Tabelas no Padrao APA 7
# Skill: relato-estatistico (v1.3.0)
# ==============================================================================

#' Formata valor de p de acordo com a norma APA 7 (sem zero a esquerda)
#' 
#' @param p Vetor numerico de probabilidades.
#' @param digits Numero de casas decimais (padrao 3).
#' @return Vetor de strings formatadas.
format_apa_p <- function(p, digits = 3) {
  sapply(p, function(val) {
    if (is.na(val)) return("—")
    if (val < 0.001) return("< .001")
    txt <- sprintf(paste0("%.", digits, "f"), val)
    sub("^0\\.", ".", txt)
  })
}

#' Formata coeficientes limitados a 1 (r, beta, cargas, omega, alfa) sem zero a esquerda
#' 
#' @param val Vetor numerico.
#' @param digits Numero de casas decimais (padrao 2).
#' @return Vetor de strings formatadas.
format_apa_cor <- function(val, digits = 2) {
  sapply(val, function(x) {
    if (is.na(x)) return("—")
    sinal <- if (x < 0) "-" else ""
    abs_x <- abs(x)
    txt <- sprintf(paste0("%.", digits, "f"), abs_x)
    txt_sem_zero <- sub("^0\\.", ".", txt)
    paste0(sinal, txt_sem_zero)
  })
}

#' Formata numeros continuos nao limitados a 1 (M, SD, t, F, chi2, z) mantendo zero a esquerda
#' 
#' @param val Vetor numerico.
#' @param digits Numero de casas decimais (padrao 2).
#' @return Vetor de strings formatadas.
format_apa_num <- function(val, digits = 2) {
  sapply(val, function(x) {
    if (is.na(x)) return("—")
    sprintf(paste0("%.", digits, "f"), x)
  })
}

#' Renderizador universal de tabelas nos formatos Markdown ou APA 7 (linhas horizontais limpas)
#' 
#' @param titulo Titulo ou identificador da tabela.
#' @param subtitulo Subtitulo ou descricao em italico.
#' @param headers Vetor de rotulos de colunas.
#' @param linhas_dados Lista de vetores de caracteres com os valores de cada linha.
#' @param alinhamentos Vetor com alinhamentos por coluna ("esq", "dir", "centro").
#' @param nota String opcional com a nota de rodape da tabela.
#' @param estilo Formato de saida: "markdown" (pipes |) ou "linhas" (APA 7 classico com 3 linhas horizontais).
renderizar_tabela_apa <- function(titulo,
                                  subtitulo,
                                  headers,
                                  linhas_dados,
                                  alinhamentos = NULL,
                                  nota = NULL,
                                  estilo = c("markdown", "linhas")) {
  estilo <- match.arg(estilo)
  if (estilo == "markdown") {
    sep_map <- c(esq = ":---", dir = "---:", centro = ":---:")
    if (is.null(alinhamentos)) alinhamentos <- c("esq", rep("centro", length(headers) - 1))
    seps <- sep_map[alinhamentos]
    seps[is.na(seps)] <- ":---:"
    
    header_str <- paste0("| ", paste(headers, collapse = " | "), " |")
    sep_str <- paste0("| ", paste(seps, collapse = " | "), " |")
    linhas_str <- vapply(linhas_dados, function(l) paste0("| ", paste(l, collapse = " | "), " |"), character(1))
    
    out <- c(
      paste0("**", titulo, "**"),
      paste0("*", subtitulo, "*"),
      "",
      header_str,
      sep_str,
      linhas_str
    )
    if (!is.null(nota)) out <- c(out, "", nota)
    return(paste(out, collapse = "\n"))
  } else {
    # Estilo "linhas": APA 7 com tres linhas horizontais e alinhamento de espacos
    limpar_md <- function(s) gsub("\\*|_|\\$|\\\\text\\{|\\}", "", s)
    headers_limpos <- limpar_md(headers)
    linhas_limpas <- lapply(linhas_dados, limpar_md)
    todas_mat <- rbind(headers_limpos, do.call(rbind, linhas_limpas))
    larguras <- apply(todas_mat, 2, function(col) max(nchar(col, type = "width"), 2))
    largura_total <- sum(larguras) + (length(larguras) - 1) * 2
    linha_div <- paste(rep("─", max(largura_total, 45)), collapse = "")
    
    pad_cell <- function(txt, w, alig) {
      n <- nchar(txt, type = "width")
      diff <- max(0, w - n)
      if (alig == "dir") paste0(paste(rep(" ", diff), collapse = ""), txt)
      else if (alig == "centro") {
        esq <- floor(diff / 2)
        dir <- ceiling(diff / 2)
        paste0(paste(rep(" ", esq), collapse = ""), txt, paste(rep(" ", dir), collapse = ""))
      } else paste0(txt, paste(rep(" ", diff), collapse = ""))
    }
    
    if (is.null(alinhamentos)) alinhamentos <- c("esq", rep("dir", length(headers) - 1))
    
    formatar_linha <- function(vals) {
      vapply(seq_along(vals), function(i) pad_cell(vals[i], larguras[i], alinhamentos[i]), character(1)) |>
        paste(collapse = "  ")
    }
    
    header_fmt <- formatar_linha(headers_limpos)
    corpo_fmt <- vapply(linhas_limpas, formatar_linha, character(1))
    
    out <- c(
      titulo,
      subtitulo,
      "",
      linha_div,
      header_fmt,
      linha_div,
      corpo_fmt,
      linha_div
    )
    if (!is.null(nota)) out <- c(out, "", nota)
    return(paste(out, collapse = "\n"))
  }
}

#' Gera Tabela de Descritivas e Correlacoes no Padrao APA 7 (Markdown ou Linhas)
#' 
#' @param df Data.frame com as variaveis numericas.
#' @param vars Vetor de nomes das variaveis a incluir.
#' @param titulos_amigaveis Vetor opcional de rotulos legiveis.
#' @param estilo Formato de saida: "markdown" (pipes |) ou "linhas" (APA 7 classico).
#' @return String com a tabela formatada.
gerar_tabela_correlacoes_apa <- function(df, vars, titulos_amigaveis = NULL, estilo = c("markdown", "linhas")) {
  estilo <- match.arg(estilo)
  df_sub <- df[, vars, drop = FALSE]
  k <- length(vars)
  
  if (is.null(titulos_amigaveis)) {
    titulos_amigaveis <- vars
  }
  
  # Calcula medias e desvios
  medias <- sapply(df_sub, mean, na.rm = TRUE)
  desvios <- sapply(df_sub, stats::sd, na.rm = TRUE)
  
  # Matriz de correlacoes
  cor_mat <- stats::cor(df_sub, use = "pairwise.complete.obs")
  
  headers <- c("Variavel", "*M*", "*SD*", as.character(1:k))
  alinhamentos <- c("esq", "dir", "dir", rep("dir", k))
  
  linhas_dados <- vector("list", k)
  for (i in 1:k) {
    vals_cor <- character(k)
    for (j in 1:k) {
      if (j == i) {
        vals_cor[j] <- "—"
      } else if (j < i) {
        vals_cor[j] <- format_apa_cor(cor_mat[i, j])
      } else {
        vals_cor[j] <- ""
      }
    }
    linhas_dados[[i]] <- c(
      sprintf("%d. %s", i, titulos_amigaveis[i]),
      format_apa_num(medias[i]),
      format_apa_num(desvios[i]),
      vals_cor
    )
  }
  
  # n por par (pairwise): relata a faixa quando ha ausentes, nao o nrow do banco
  n_par <- crossprod(!is.na(as.matrix(df_sub)))
  n_inf <- n_par[lower.tri(n_par)]
  n_txt <- if (length(unique(n_inf)) <= 1) sprintf("$N = %d$", as.integer(n_inf[1])) else
    sprintf("$n$ por par entre %d e %d (ausentes tratados par a par)", as.integer(min(n_inf)), as.integer(max(n_inf)))
  nota <- sprintf("*Nota.* %s. *M* = media; *SD* = desvio-padrao. Coeficientes na diagonal inferior representam correlacoes de Pearson.", n_txt)
  
  renderizar_tabela_apa(
    titulo = "Tabela",
    subtitulo = "Medias, Desvios-Padrao e Correlacoes de Pearson para as Variaveis do Estudo",
    headers = headers,
    linhas_dados = linhas_dados,
    alinhamentos = alinhamentos,
    nota = nota,
    estilo = estilo
  )
}

#' Gera Tabela de Regressao Linear Multipla no Padrao APA 7 (Markdown ou Linhas)
#' 
#' @param fit Objeto de modelo ajustado pela funcao lm().
#' @param nomes_preditores Vetor opcional com nomes amigaveis dos preditores.
#' @param se_tipo Tipo de erro-padrao: "padrao" (OLS) ou "HC3" (robusto a heterocedasticidade).
#' @param estilo Formato de saida: "markdown" (pipes |) ou "linhas" (APA 7 classico).
#' @return String com a tabela formatada.
gerar_tabela_regressao_apa <- function(fit, nomes_preditores = NULL, se_tipo = c("padrao", "HC3"), estilo = c("markdown", "linhas")) {
  estilo <- match.arg(estilo)
  se_tipo <- match.arg(se_tipo)
  resumo <- summary(fit)
  
  if (se_tipo == "HC3" && requireNamespace("sandwich", quietly = TRUE)) {
    vcov_hc3 <- sandwich::vcovHC(fit, type = "HC3")
    coefs <- summary(fit)$coefficients
    se_hc3 <- sqrt(diag(vcov_hc3))
    t_hc3 <- coefs[, 1] / se_hc3
    p_hc3 <- 2 * (1 - stats::pt(abs(t_hc3), df = fit$df.residual))
    coefs <- cbind(coefs[, 1], se_hc3, t_hc3, p_hc3)
    colnames(coefs) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    ci <- cbind(coefs[, 1] - 1.96 * se_hc3, coefs[, 1] + 1.96 * se_hc3)
  } else {
    coefs <- resumo$coefficients
    ci <- stats::confint(fit)
  }
  
  # Calcula betas padronizados
  df_model <- fit$model
  sd_y <- stats::sd(df_model[[1]], na.rm = TRUE)
  
  p_vars <- nrow(coefs)
  betas_padrao <- rep(NA_real_, p_vars)
  
  for (i in 2:p_vars) {
    nome_var <- rownames(coefs)[i]
    if (nome_var %in% colnames(df_model)) {
      sd_x <- stats::sd(df_model[[nome_var]], na.rm = TRUE)
      betas_padrao[i] <- coefs[i, 1] * (sd_x / sd_y)
    }
  }
  
  if (is.null(nomes_preditores)) {
    nomes_preditores <- rownames(coefs)
  }
  
  se_rotulo <- if (se_tipo == "HC3") "*SE* (HC3)" else "*SE*"
  headers <- c("Preditor", "*B*", se_rotulo, "*beta*", "*t*", "*p*", "IC 95% para *B*")
  alinhamentos <- c("esq", "dir", "dir", "dir", "dir", "dir", "dir")
  
  linhas_dados <- vector("list", p_vars)
  for (i in 1:p_vars) {
    beta_txt <- if (is.na(betas_padrao[i])) "—" else format_apa_cor(betas_padrao[i])
    ci_txt <- sprintf("[%s, %s]", format_apa_num(ci[i, 1]), format_apa_num(ci[i, 2]))
    linhas_dados[[i]] <- c(
      nomes_preditores[i],
      format_apa_num(coefs[i, 1]),
      format_apa_num(coefs[i, 2]),
      beta_txt,
      format_apa_num(coefs[i, 3]),
      format_apa_p(coefs[i, 4]),
      ci_txt
    )
  }
  
  f_stat <- resumo$fstatistic
  f_p <- 1 - stats::pf(f_stat[1], f_stat[2], f_stat[3])
  nota_se <- if (se_tipo == "HC3") " Erros-padrao robustos estimados via HC3 (sandwich)." else ""
  
  nota <- sprintf("*Nota.* $R^2 = %s$, $R^2_{\\text{adj}} = %s$, $F(%d, %d) = %s$, $p %s$. *B* = coeficiente nao padronizado; *SE* = erro-padrao; *beta* = coeficiente padronizado; IC = intervalo de confianca.%s",
                  format_apa_cor(resumo$r.squared),
                  format_apa_cor(resumo$adj.r.squared),
                  round(f_stat[2]),
                  round(f_stat[3]),
                  format_apa_num(f_stat[1]),
                  if (f_p < .001) "< .001" else paste("=", format_apa_p(f_p)),
                  nota_se)
  
  renderizar_tabela_apa(
    titulo = "Tabela",
    subtitulo = "Resultados do Modelo de Regressao Linear Multipla",
    headers = headers,
    linhas_dados = linhas_dados,
    alinhamentos = alinhamentos,
    nota = nota,
    estilo = estilo
  )
}

#' Gera Tabela de Cargas Fatoriais e Ajuste da AFC no Padrao APA 7 (Markdown ou Linhas)
#' 
#' @param fit Objeto de modelo ajustado pelo lavaan.
#' @param titulos_fatores Vetor nomeado opcional para rotular fatores.
#' @param rotulos_itens Vetor nomeado opcional para rotular itens.
#' @param estilo Formato de saida: "markdown" (pipes |) ou "linhas" (APA 7 classico).
#' @return String com a tabela formatada.
gerar_tabela_afc_apa <- function(fit, titulos_fatores = NULL, rotulos_itens = NULL, estilo = c("markdown", "linhas")) {
  estilo <- match.arg(estilo)
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan e necessario.")
  
  std <- lavaan::standardizedSolution(fit)
  cargas <- subset(std, op == "=~")
  
  residuos <- subset(std, op == "~~" & lhs == rhs & lhs %in% cargas$rhs)
  res_map <- stats::setNames(residuos$est.std, residuos$lhs)
  
  headers <- c("Item", "Fator", "*lambda* Padrao", "*SE*", "*z*", "*p*", "Unicidade (*delta*)", "*R*^2")
  alinhamentos <- c("esq", "esq", "dir", "dir", "dir", "dir", "dir", "dir")
  
  linhas_dados <- vector("list", nrow(cargas))
  for (i in seq_len(nrow(cargas))) {
    it <- cargas$rhs[i]
    f <- cargas$lhs[i]
    it_nome <- if (!is.null(rotulos_itens) && it %in% names(rotulos_itens)) rotulos_itens[[it]] else it
    f_nome <- if (!is.null(titulos_fatores) && f %in% names(titulos_fatores)) titulos_fatores[[f]] else f
    
    lam_std <- cargas$est.std[i]
    se_val <- cargas$se[i]
    z_val <- cargas$z[i]
    p_val <- cargas$pvalue[i]
    delta_val <- if (it %in% names(res_map)) res_map[[it]] else (1 - lam_std^2)
    r2_val <- lam_std^2
    
    linhas_dados[[i]] <- c(
      it_nome,
      f_nome,
      format_apa_cor(lam_std),
      format_apa_num(se_val),
      format_apa_num(z_val),
      format_apa_p(p_val),
      format_apa_cor(delta_val),
      format_apa_cor(r2_val)
    )
  }
  
  fm <- lavaan::fitMeasures(fit)
  pega <- function(nm) if (nm %in% names(fm)) fm[[nm]] else NA_real_
  
  # Uma unica versao de indices (robusta, senao escalonada, senao padrao), declarada na nota
  versao <- if (!is.na(pega("cfi.robust"))) "robust" else if (!is.na(pega("cfi.scaled"))) "scaled" else ""
  suf <- if (nzchar(versao)) paste0(".", versao) else ""
  rot_versao <- c(robust = "robustos", scaled = "escalonados")[versao]
  if (is.na(rot_versao)) rot_versao <- "padrao"
  cfi <- pega(paste0("cfi", suf)); tli <- pega(paste0("tli", suf)); rmsea <- pega(paste0("rmsea", suf))
  rmsea_li <- pega(paste0("rmsea.ci.lower", suf)); rmsea_ls <- pega(paste0("rmsea.ci.upper", suf)); srmr <- pega("srmr")
  chisq <- if (!is.na(pega("chisq.scaled"))) pega("chisq.scaled") else pega("chisq")
  df_fit <- if (!is.na(pega("df.scaled"))) pega("df.scaled") else pega("df")
  p_chi <- if (!is.na(pega("pvalue.scaled"))) pega("pvalue.scaled") else pega("pvalue")
  n_obs <- lavaan::lavInspect(fit, "nobs")
  
  nota_ajuste <- sprintf(
    "*Nota.* $N = %d$. *lambda* = carga fatorial padronizada; *SE* = erro-padrao; *delta* = variancia residual unica (unicidade); *R*^2 = proporcao de variancia explicada pelo fator. Ajuste global (indices %s): $\\chi^2(%d) = %s$, $p %s$, $\\text{CFI} = %s$, $\\text{TLI} = %s$, $\\text{RMSEA} = %s$, IC 90%% [%s, %s], $\\text{SRMR} = %s$.",
    n_obs, rot_versao, round(df_fit), format_apa_num(chisq),
    if (p_chi < .001) "< .001" else paste("=", format_apa_p(p_chi)),
    format_apa_cor(cfi, 3), format_apa_cor(tli, 3), format_apa_num(rmsea, 3), format_apa_num(rmsea_li, 3), format_apa_num(rmsea_ls, 3),
    format_apa_num(srmr, 3)
  )
  
  renderizar_tabela_apa(
    titulo = "Tabela",
    subtitulo = "Cargas Fatoriais Padronizadas, Erros-Padrao e Unicidades da Analise Fatorial Confirmatoria",
    headers = headers,
    linhas_dados = linhas_dados,
    alinhamentos = alinhamentos,
    nota = nota_ajuste,
    estilo = estilo
  )
}

#' Gera Tabela de Confiabilidade e Fidedignidade das Escalas (APA 7 Markdown ou Linhas)
#' 
#' @param tab_conf Data.frame retornado pela funcao confiabilidade() do psicometria-r.
#' @param estilo Formato de saida: "markdown" (pipes |) ou "linhas" (APA 7 classico).
#' @return String com a tabela formatada.
gerar_tabela_confiabilidade_apa <- function(tab_conf, estilo = c("markdown", "linhas")) {
  estilo <- match.arg(estilo)
  tem_ic <- "omega_ic_inf" %in% colnames(tab_conf)
  tem_delta <- "omega_delta" %in% colnames(tab_conf)
  
  if (tem_ic && tem_delta) {
    headers <- c("Fator / Escala", "Itens", "*omega* (escore observado)", "*omega* delta", "*EP*", "IC 95% *omega* delta", "*alpha*", "IC 95% *alpha* (Feldt)")
    alinhamentos <- c("esq", "centro", "dir", "dir", "dir", "dir", "dir", "dir")
    linhas_dados <- vector("list", nrow(tab_conf))
    for (i in seq_len(nrow(tab_conf))) {
      linhas_dados[[i]] <- c(
        as.character(tab_conf$fator[i]),
        as.character(tab_conf$itens[i]),
        format_apa_cor(tab_conf$omega[i], 3),
        format_apa_cor(tab_conf$omega_delta[i], 3),
        format_apa_num(tab_conf$omega_ep[i], 3),
        sprintf("[%s, %s]", format_apa_cor(tab_conf$omega_ic_inf[i], 3), format_apa_cor(tab_conf$omega_ic_sup[i], 3)),
        format_apa_cor(tab_conf$alfa[i], 3),
        sprintf("[%s, %s]", format_apa_cor(tab_conf$alfa_ic_inf[i], 3), format_apa_cor(tab_conf$alfa_ic_sup[i], 3))
      )
    }
    nota <- "*Nota.* *omega* (escore observado) = semTools::compRelSEM; *omega* delta = omega da variavel-resposta latente nos itens ordinais (igual ao anterior em itens continuos), com *EP* e IC Wald pelo metodo delta; *alpha* = alfa de Cronbach com IC exato de Feldt."
  } else if (tem_ic) {
    headers <- c("Fator / Escala", "Itens", "*omega* de McDonald", "IC 95% *omega*", "*alpha* de Cronbach", "IC 95% *alpha*")
    alinhamentos <- c("esq", "centro", "dir", "dir", "dir", "dir")
    linhas_dados <- vector("list", nrow(tab_conf))
    for (i in seq_len(nrow(tab_conf))) {
      linhas_dados[[i]] <- c(
        as.character(tab_conf$fator[i]),
        as.character(tab_conf$itens[i]),
        format_apa_cor(tab_conf$omega[i], 3),
        sprintf("[%s, %s]", format_apa_cor(tab_conf$omega_ic_inf[i], 3), format_apa_cor(tab_conf$omega_ic_sup[i], 3)),
        format_apa_cor(tab_conf$alfa[i], 3),
        sprintf("[%s, %s]", format_apa_cor(tab_conf$alfa_ic_inf[i], 3), format_apa_cor(tab_conf$alfa_ic_sup[i], 3))
      )
    }
    nota <- "*Nota.* *omega* = coeficiente omega de McDonald; *alpha* = alfa de Cronbach; IC 95% percentilico por bootstrap."
  } else {
    headers <- c("Fator / Escala", "Itens", "*omega* de McDonald", "*alpha* de Cronbach")
    alinhamentos <- c("esq", "centro", "dir", "dir")
    linhas_dados <- vector("list", nrow(tab_conf))
    for (i in seq_len(nrow(tab_conf))) {
      linhas_dados[[i]] <- c(
        as.character(tab_conf$fator[i]),
        as.character(tab_conf$itens[i]),
        format_apa_cor(tab_conf$omega[i], 3),
        format_apa_cor(tab_conf$alfa[i], 3)
      )
    }
    nota <- "*Nota.* *omega* = coeficiente omega de McDonald; *alpha* = alfa de Cronbach."
  }
  
  renderizar_tabela_apa(
    titulo = "Tabela",
    subtitulo = "Indices de Fidedignidade e Consistencia Interna das Subescalas",
    headers = headers,
    linhas_dados = linhas_dados,
    alinhamentos = alinhamentos,
    nota = nota,
    estilo = estilo
  )
}
