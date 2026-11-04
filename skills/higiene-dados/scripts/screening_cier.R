# ==============================================================================
# screening_cier.R
# Modulo de Triagem de Respostas Descuidadas (C/IER) e Outliers Robustos
# Skill: higiene-dados (v1.1)
# ==============================================================================

#' Recodifica vetores ou colunas de fatores/strings de Likert para inteiros ordenados
#' 
#' Detecta padroes comuns como "0 - Nao se aplicou", "1 - Pouco", "(2) Moderado"
#' ou converte via mapa de niveis declarado.
#' 
#' @param x Vetor ou data.frame de itens.
#' @param mapa_niveis Vetor nomeado opcional (ex: c("Discordo" = 1, "Neutro" = 2, "Concordo" = 3))
#' @return Vetor numerico ou data.frame com valores inteiros limpos.
recodificar_likert_labels <- function(x, mapa_niveis = NULL) {
  # v1.0.1: sem mapa_niveis, rotulos sem numero param com erro. A versao 1.0 caia em as.numeric(as.factor()),
  # que ordena alfabeticamente (na DASS-21, "Nao aconteceu comigo" virava 4 em vez de 0).
  norm <- function(s) gsub("\\s+", " ", trimws(s))
  if (!is.null(mapa_niveis)) names(mapa_niveis) <- norm(names(mapa_niveis))
  recod_vetor <- function(v, nome = "vetor") {
    if (is.numeric(v)) return(v)
    v_char <- norm(as.character(v)); v_char[v_char == ""] <- NA; obs <- !is.na(v_char)
    if (!is.null(mapa_niveis)) {
      v_num <- unname(mapa_niveis[v_char]); fora <- obs & is.na(v_num)
      if (any(fora)) { u <- unique(v_char[fora])
        warning(sprintf("%s: %d valor(es) fora do mapa_niveis viraram NA (%s)", nome, sum(fora), paste(utils::head(u, 3), collapse = " | ")), call. = FALSE) }
      return(as.numeric(v_num))
    }
    tem_dig <- grepl("^[^0-9]*[0-9]+", v_char)
    if (any(obs) && all(tem_dig[obs])) return(as.numeric(sub("^[^0-9]*([0-9]+).*", "\\1", v_char)))
    if (!any(obs)) return(rep(NA_real_, length(v_char)))
    stop(sprintf("%s: rotulos sem numero e sem mapa_niveis. Passe mapa_niveis = c(\"rotulo\" = valor, ...). Rotulos: %s",
                 nome, paste(sort(unique(v_char[obs])), collapse = " | ")), call. = FALSE)
  }
  if (is.data.frame(x)) { out <- x; for (nm in names(x)) out[[nm]] <- recod_vetor(x[[nm]], nm); out } else recod_vetor(x)
}

#' Calcula o indice Longstring para cada respondente
#' 
#' @param x Matriz ou data.frame contendo apenas os itens da escala.
#' @param ignorar_todos_zeros Booleano; se TRUE, sequencias completas de zero em inventarios clinicos unipolares nao sao tratadas como straightlining descuidado.
#' @param val_zero Valor numerico que define a resposta assintomatica (padrao 0).
#' @return Vetor numerico com o comprimento da maior sequencia consecutiva identica.
calc_longstring <- function(x, ignorar_todos_zeros = FALSE, val_zero = 0) {
  x_mat <- as.matrix(x)
  apply(x_mat, 1, function(row) {
    clean_row <- row[!is.na(row)]
    if (length(clean_row) == 0) return(0)
    if (isTRUE(ignorar_todos_zeros) && all(clean_row == val_zero)) {
      return(0)
    }
    rle_res <- rle(clean_row)
    max(rle_res$lengths)
  })
}

#' Calcula a Variancia Intraindividual de Resposta (IRV)
#' 
#' Formula: IRV_i = sqrt( (1 / (K - 1)) * sum((X_ij - X_bar_i)^2) )
#' 
#' @param x Matriz ou data.frame contendo os itens da escala.
#' @return Vetor numerico com o desvio-padrao intraindividual por linha.
calc_irv <- function(x) {
  x_mat <- as.matrix(x)
  apply(x_mat, 1, function(row) {
    clean_row <- row[!is.na(row)]
    if (length(clean_row) < 2) return(NA_real_)
    stats::sd(clean_row)
  })
}

#' Calcula Pares de Sinonimos e Antonimos Psicometricos
#' 
#' @param x Data.frame com itens numericos.
#' @param critval Limiar absoluto de correlacao empirica para selecionar pares (|r| >= critval).
#' @return Data.frame com as colunas cor_sinonimos e cor_antonimos por respondente.
calc_psychsyn_psychant <- function(x, critval = 0.60) {
  x_mat <- as.matrix(x)
  cor_mat <- stats::cor(x_mat, use = "pairwise.complete.obs")
  diag(cor_mat) <- 0
  
  # Identifica pares de sinonimos (r >= critval)
  sin_pairs <- which(cor_mat >= critval & upper.tri(cor_mat), arr.ind = TRUE)
  # Identifica pares de antonimos (r <= -critval)
  ant_pairs <- which(cor_mat <= -critval & upper.tri(cor_mat), arr.ind = TRUE)
  
  n <- nrow(x_mat)
  cor_sin <- rep(NA_real_, n)
  cor_ant <- rep(NA_real_, n)
  
  if (nrow(sin_pairs) >= 3) {
    cor_sin <- apply(x_mat, 1, function(row) {
      vals_1 <- row[sin_pairs[, 1]]
      vals_2 <- row[sin_pairs[, 2]]
      valid <- !is.na(vals_1) & !is.na(vals_2)
      if (sum(valid) >= 3 && stats::sd(vals_1[valid]) > 0 && stats::sd(vals_2[valid]) > 0) {
        stats::cor(vals_1[valid], vals_2[valid])
      } else {
        NA_real_
      }
    })
  }
  
  if (nrow(ant_pairs) >= 3) {
    cor_ant <- apply(x_mat, 1, function(row) {
      vals_1 <- row[ant_pairs[, 1]]
      vals_2 <- row[ant_pairs[, 2]]
      valid <- !is.na(vals_1) & !is.na(vals_2)
      if (sum(valid) >= 3 && stats::sd(vals_1[valid]) > 0 && stats::sd(vals_2[valid]) > 0) {
        stats::cor(vals_1[valid], vals_2[valid])
      } else {
        NA_real_
      }
    })
  }
  
  data.frame(
    cor_sinonimos = cor_sin,
    cor_antonimos = cor_ant,
    n_pares_sinonimos = nrow(sin_pairs),
    n_pares_antonimos = nrow(ant_pairs)
  )
}

#' Deteccao de Outliers Univariados via MAD (Leys et al., 2013)
#' 
#' @param v Vetor numerico continuo.
#' @param threshold Corte de desvios absolutos medianos (padrao 3.0).
#' @return Lista contendo indices booleanos de outliers, mediana, valor do MAD e escores normalizados.
calc_mad_outliers <- function(v, threshold = 3.0) {
  clean_v <- stats::na.omit(v)
  med <- stats::median(clean_v)
  # Constante de consistencia para normalidade b = 1.4826
  mad_val <- 1.4826 * stats::median(abs(clean_v - med))
  
  if (mad_val == 0) {
    # Fallback se mais de 50% dos dados forem identicos
    mad_val <- stats::IQR(clean_v) / 1.349
    if (mad_val == 0) mad_val <- 1e-6
  }
  
  escores_mad <- abs(v - med) / mad_val
  is_outlier <- !is.na(v) & (escores_mad > threshold)
  
  list(
    is_outlier = is_outlier,
    mediana = med,
    mad = mad_val,
    escores_mad = escores_mad
  )
}

#' Deteccao de Outliers Multivariados via Mahalanobis / FastMCD
#' 
#' Implementa busca em grade sobre alpha em covMcd para evitar singularidade em hiperplano
#' causada por inflacao de zeros, com fallback regularizado e semente de reproduzibilidade.
#' 
#' @param x Data.frame ou matriz contendo variaveis continuas.
#' @param p_val Limiar de probabilidade sob qui-quadrado (padrao .001).
#' @param semente Semente numerica para reproducibilidade do algoritmo estocastico FastMCD (padrao 2026).
#' @return Lista com distancias quadradas, valores criticos e flag booleana.
calc_mahad_robust <- function(x, p_val = 0.001, semente = 2026, alpha_mcd = 0.75) {
  # v1.2.0: MCD75 por padrao (ponto de ruptura .25; Leys et al., 2019, que reservam o MCD50 para suspeita de
  # mais de 25% de atipicos). Com matriz singular a grade sobe (.85, .90): h maior tira o subconjunto do hiperplano.
  x_mat <- as.matrix(stats::na.omit(x))
  p <- ncol(x_mat)
  n <- nrow(x_mat)
  corte_chi2 <- stats::qchisq(1 - p_val, df = p)
  
  dist_sq <- rep(NA_real_, nrow(x))
  idx_completos <- which(stats::complete.cases(x))
  
  usou_mcd <- FALSE
  estimador_usado <- "Mahalanobis Padrao Regularizado"
  
  if (requireNamespace("robustbase", quietly = TRUE) && n > (p * 2)) {
    if (!is.null(semente)) {
      if (exists(".Random.seed", envir = globalenv())) {
        velha_seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", velha_seed, envir = globalenv()), add = TRUE)
      } else {
        on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
      }
      set.seed(semente)
    }
    
    # Grade de busca sobre alpha para contornar singularidade por inflacao de zeros
    for (a_cand in unique(c(alpha_mcd, c(0.75, 0.85, 0.90)[c(0.75, 0.85, 0.90) > alpha_mcd]))) {
      mcd_fit <- tryCatch(suppressWarnings(robustbase::covMcd(x_mat, alpha = a_cand)), error = function(e) NULL)
      sing <- is.null(mcd_fit) || !is.null(mcd_fit$singularity) || is.null(mcd_fit$cov) || (rcond(mcd_fit$cov) < 1e-12)
      if (!sing) {
        dist_sq[idx_completos] <- if (!is.null(mcd_fit$mah)) mcd_fit$mah else stats::mahalanobis(x_mat, center = mcd_fit$center, cov = mcd_fit$cov)
        usou_mcd <- TRUE
        estimador_usado <- paste0("FastMCD (alpha = ", a_cand, ")")
        break
      }
    }
  }
  
  if (!usou_mcd) {
    # Fallback robusto via centroide mediano e matriz de covariancia regularizada
    center <- apply(x_mat, 2, stats::median)
    cov_mat <- stats::cov(x_mat)
    if (rcond(cov_mat) < 1e-12) {
      cov_mat <- cov_mat + diag(1e-4, p)
    }
    dist_sq[idx_completos] <- stats::mahalanobis(x_mat, center = center, cov = cov_mat)
  }
  
  is_outlier <- !is.na(dist_sq) & (dist_sq > corte_chi2)
  
  list(
    dist_mahalanobis_sq = dist_sq,
    corte_chi2 = corte_chi2,
    is_outlier = is_outlier,
    estimador = estimador_usado
  )
}

#' Indice even-odd (Jackson, 1976; Meade & Craig, 2012), como careless::evenodd, para itens nao contiguos
#'
#' @param x Matriz ou data.frame de itens.
#' @param fatores Lista de vetores (nomes ou indices) por subescala, na ordem de apresentacao, ou vetor de tamanhos
#'   para subescalas contiguas (a forma do careless).
#' @return Vetor em que valor maior indica mais descuido: menos a correlacao intrapessoal entre as medias pares e
#'   impares das subescalas, com Spearman-Brown. NA quando ha menos de 2 subescalas validas ou variancia nula.
calc_evenodd <- function(x, fatores, diag = FALSE) {
  x <- as.matrix(x)
  if (!is.list(fatores)) { ini <- cumsum(c(1, utils::head(fatores, -1))); fatores <- mapply(function(a, k) a:(a + k - 1), ini, fatores, SIMPLIFY = FALSE) }
  if (length(fatores) < 2) stop("even-odd exige pelo menos 2 subescalas.", call. = FALSE)
  idx <- lapply(fatores, function(f) if (is.character(f)) match(f, colnames(x)) else as.integer(f))
  if (anyNA(unlist(idx))) stop("itens de 'fatores' ausentes em x.", call. = FALSE)
  um <- function(r) {
    m <- t(vapply(idx, function(ii) { v <- r[ii]; pos <- seq_along(v)
      c(mean(v[pos %% 2 == 0], na.rm = TRUE), mean(v[pos %% 2 == 1], na.rm = TRUE)) }, numeric(2)))
    ok <- stats::complete.cases(m); n_ok <- sum(ok)
    if (n_ok < 2) return(c(NA_real_, n_ok))
    val <- suppressWarnings(stats::cor(m[ok, 1], m[ok, 2]))
    val <- (2 * val) / (1 + val); if (!is.na(val) && val < -1) val <- -1
    c(-val, n_ok)
  }
  r <- t(apply(x, 1, um))
  if (isTRUE(diag)) data.frame(evenodd = r[, 1], n_subescalas = r[, 2]) else r[, 1]
}

#' Distancia de Mahalanobis classica nos itens como indice de C/IER (Meade & Craig, 2012; DeSimone et al., 2015)
#'
#' D2 de cada caso completo contra o vetor de medias dos itens, comparada com qui-quadrado(k). DeSimone et al.
#' (2015) sugerem os 5% superiores (p_corte = .05) e lembram que a distribuicao qui-quadrado supoe itens normais;
#' a skill usa .001 por padrao, mais conservador com itens ordinais. Meade e Craig (2012) calculam um D por fator
#' quando ha muitos itens e tiram a media: 'fatores' faz isso.
calc_mahad_cier <- function(x, p_corte = 0.001, fatores = NULL) {
  x <- as.matrix(x)
  um <- function(M) {
    d2 <- rep(NA_real_, nrow(M)); ok <- stats::complete.cases(M)
    if (sum(ok) <= ncol(M) + 1) return(list(d2 = d2, p = d2))
    S <- stats::cov(M[ok, , drop = FALSE]); if (rcond(S) < 1e-12) S <- S + diag(1e-6, ncol(M))
    d2[ok] <- stats::mahalanobis(M[ok, , drop = FALSE], colMeans(M[ok, , drop = FALSE]), S)
    list(d2 = d2, p = stats::pchisq(d2, ncol(M), lower.tail = FALSE))
  }
  if (is.null(fatores)) { r <- um(x); return(data.frame(mahad_d2 = r$d2, mahad_p = r$p, flag_mahad = !is.na(r$p) & r$p < p_corte)) }
  idx <- lapply(fatores, function(f) if (is.character(f)) match(f, colnames(x)) else as.integer(f))
  ps <- sapply(idx, function(ii) um(x[, ii, drop = FALSE])$p)
  d2s <- sapply(idx, function(ii) um(x[, ii, drop = FALSE])$d2)
  # p combinado conservador: o menor p entre fatores com correcao de Bonferroni
  pmin_b <- apply(ps, 1, function(v) if (all(is.na(v))) NA_real_ else min(1, min(v, na.rm = TRUE) * sum(!is.na(v))))
  data.frame(mahad_d2 = rowMeans(d2s, na.rm = TRUE), mahad_p = pmin_b, flag_mahad = !is.na(pmin_b) & pmin_b < p_corte)
}

#' Itens de atencao (instrucao direta ou infrequencia; Meade & Craig, 2012; Huang et al., 2015)
#'
#' @param x data.frame ou matriz so com os itens de atencao
#' @param esperado vetor com a resposta correta de cada item (mesma ordem das colunas); NA = qualquer resposta
#'   diferente de 'incorreto' e aceita
#' @param incorreto lista opcional de respostas que contam como falha em cada item (itens de infrequencia)
#' @return numero de falhas por caso (NA nao conta como falha)
calc_falhas_atencao <- function(x, esperado = NULL, incorreto = NULL) {
  x <- as.data.frame(x); k <- ncol(x); f <- matrix(FALSE, nrow(x), k)
  for (j in seq_len(k)) {
    v <- x[[j]]
    if (!is.null(esperado) && !is.na(esperado[j])) f[, j] <- !is.na(v) & as.character(v) != as.character(esperado[j])
    if (!is.null(incorreto) && length(incorreto) >= j && length(incorreto[[j]])) f[, j] <- f[, j] | (!is.na(v) & as.character(v) %in% as.character(incorreto[[j]]))
  }
  rowSums(f)
}

#' Triagem de respostas descuidadas por familias de indices (v1.3.0)
#'
#' Repeticao (longstring >= metade da escala ou IRV = 0) conta como uma familia so, porque IRV = 0 implica
#' longstring = k. As outras familias sao tempo (< 2 s por item; Huang et al., 2012, citado por Curran, 2016),
#' psychsyn (correlacao intrapessoal < 0 nos pares com r > critval; Curran, 2016) e even-odd (Spearman-Brown < 0,
#' quando 'fatores_evenodd' e dado). Exclusao exige 2 familias (Ward & Meade, 2023). Um indice que assinala mais
#' que 'limite_informativo' dos casos em que e calculavel fica fora do criterio, com aviso (heuristica da skill).
triagem_cier_completa <- function(df_itens,
                                  tempo_total_seg = NULL,
                                  n_itens_total = ncol(df_itens),
                                  corte_seg_por_item = 2,
                                  max_longstring_proporcao = 0.5,
                                  ignorar_todos_zeros = FALSE,
                                  val_zero = 0,
                                  usar_psychsyn = TRUE,
                                  critval_psychsyn = 0.60,
                                  corte_psychsyn = 0,
                                  fatores_evenodd = NULL,
                                  corte_evenodd = 0,
                                  limite_informativo = 0.20,
                                  usar_mahad = FALSE,
                                  p_mahad = 0.001,
                                  fatores_mahad = NULL,
                                  itens_atencao = NULL,
                                  esperado_atencao = NULL,
                                  incorreto_atencao = NULL,
                                  max_falhas_atencao = 0,
                                  flag_ajuste_pessoa = NULL) {
  M <- as.matrix(df_itens); n <- nrow(M); k <- ncol(M); avisos <- character()
  todos_zeros <- apply(M, 1, function(r) { r <- r[!is.na(r)]; length(r) > 0 && all(r == val_zero) })
  ls_vals <- calc_longstring(M, ignorar_todos_zeros = ignorar_todos_zeros, val_zero = val_zero)
  limiar_ls <- max(3, ceiling(k * max_longstring_proporcao))
  flag_longstring <- ls_vals >= limiar_ls
  irv_vals <- calc_irv(M)
  flag_irv_zero <- !is.na(irv_vals) & irv_vals == 0
  if (isTRUE(ignorar_todos_zeros)) flag_irv_zero <- flag_irv_zero & !todos_zeros
  fam_repeticao <- flag_longstring | flag_irv_zero
  tempo_por_item <- rep(NA_real_, n); fam_tempo <- rep(FALSE, n)
  if (!is.null(tempo_total_seg) && length(tempo_total_seg) == n) {
    tempo_por_item <- tempo_total_seg / n_itens_total
    fam_tempo <- !is.na(tempo_por_item) & tempo_por_item < corte_seg_por_item
  } else avisos <- c(avisos, "tempo de resposta nao informado: familia tempo nao calculada")
  informativo <- function(flag, calc, nome) {
    p <- if (sum(calc) > 0) sum(flag & calc) / sum(calc) else NA
    if (!is.na(p) && p > limite_informativo) {
      avisos <<- c(avisos, sprintf("%s assinalou %.0f%% dos casos calculaveis (> %.0f%%): indice nao discrimina nesta amostra e ficou fora do criterio", nome, 100 * p, 100 * limite_informativo))
      return(FALSE) }
    TRUE
  }
  syn <- rep(NA_real_, n); fam_psychsyn <- rep(FALSE, n); usou_syn <- FALSE
  if (isTRUE(usar_psychsyn)) {
    ps <- tryCatch(calc_psychsyn_psychant(M, critval = critval_psychsyn), error = function(e) NULL)
    if (!is.null(ps) && ps$n_pares_sinonimos[1] >= 3) {
      syn <- ps$cor_sinonimos; f <- !is.na(syn) & syn < corte_psychsyn
      if (informativo(f, !is.na(syn), "psychsyn")) { fam_psychsyn <- f; usou_syn <- TRUE }
    } else avisos <- c(avisos, "psychsyn: menos de 3 pares de sinonimos; nao calculado")
  }
  eo <- rep(NA_real_, n); fam_evenodd <- rep(FALSE, n); usou_eo <- FALSE
  if (!is.null(fatores_evenodd)) {
    if (length(fatores_evenodd) < 4) avisos <- c(avisos, sprintf("even-odd com %d subescalas: a correlacao intrapessoal usa %d pontos e e instavel", length(fatores_evenodd), length(fatores_evenodd)))
    eo <- calc_evenodd(M, fatores_evenodd); f <- !is.na(eo) & eo > corte_evenodd
    if (informativo(f, !is.na(eo), "even-odd")) { fam_evenodd <- f; usou_eo <- TRUE }
  }
  md <- rep(NA_real_, n); fam_mahad <- rep(FALSE, n); usou_md <- FALSE
  if (isTRUE(usar_mahad)) {
    mr <- calc_mahad_cier(M, p_corte = p_mahad, fatores = fatores_mahad); md <- mr$mahad_p
    if (informativo(mr$flag_mahad, !is.na(mr$mahad_p), "Mahalanobis")) { fam_mahad <- mr$flag_mahad; usou_md <- TRUE }
  }
  fa <- rep(NA_real_, n); fam_atencao <- rep(FALSE, n)
  if (!is.null(itens_atencao)) {
    if (nrow(as.data.frame(itens_atencao)) != n) stop("itens_atencao com numero de linhas diferente de df_itens.", call. = FALSE)
    fa <- calc_falhas_atencao(itens_atencao, esperado_atencao, incorreto_atencao); fam_atencao <- fa > max_falhas_atencao
  }
  fam_ajuste <- rep(FALSE, n)
  if (!is.null(flag_ajuste_pessoa)) {
    if (is.data.frame(flag_ajuste_pessoa)) flag_ajuste_pessoa <- flag_ajuste_pessoa$flag_ajuste_pessoa
    if (length(flag_ajuste_pessoa) != n) stop("flag_ajuste_pessoa com comprimento diferente de df_itens.", call. = FALSE)
    fam_ajuste <- !is.na(flag_ajuste_pessoa) & as.logical(flag_ajuste_pessoa)
  }
  n_infr <- as.integer(fam_repeticao) + as.integer(fam_tempo) + as.integer(fam_psychsyn) + as.integer(fam_evenodd) +
    as.integer(fam_mahad) + as.integer(fam_atencao) + as.integer(fam_ajuste)
  out <- data.frame(longstring = ls_vals, flag_longstring = flag_longstring, irv = round(irv_vals, 3),
    flag_irv_zero = flag_irv_zero, familia_repeticao = fam_repeticao, tempo_por_item_seg = round(tempo_por_item, 2),
    flag_tempo_rapido = fam_tempo, psychsyn = round(syn, 3), flag_psychsyn = fam_psychsyn, evenodd = round(eo, 3),
    flag_evenodd = fam_evenodd, mahad_p = signif(md, 3), flag_mahad = fam_mahad, falhas_atencao = fa,
    flag_atencao = fam_atencao, flag_ajuste_pessoa = fam_ajuste, total_infracoes_cier = n_infr, excluir_por_cier = n_infr >= 2)
  attr(out, "decisoes") <- c(sprintf("longstring: corte >= %d de %d itens (metade da escala; Curran, 2016)", limiar_ls, k),
    "repeticao = longstring ou IRV = 0, contada como uma familia (IRV = 0 implica longstring = k)",
    sprintf("psychsyn %s; even-odd %s; Mahalanobis %s; atencao %s; ajuste de pessoa %s; exclusao com 2 ou mais familias",
      if (usou_syn) "no criterio" else "fora do criterio", if (usou_eo) "no criterio" else "fora do criterio",
      if (usou_md) sprintf("no criterio (p < %s contra qui-quadrado com k gl)", sub("^0", "", format(p_mahad))) else "fora do criterio",
      if (!is.null(itens_atencao)) sprintf("no criterio (mais de %d falha[s])", max_falhas_atencao) else "nao aplicada",
      if (!is.null(flag_ajuste_pessoa)) "no criterio (corte empirico de calc_ajuste_pessoa(); Niessen et al., 2016)" else "nao aplicado"))
  attr(out, "avisos") <- avisos
  if (length(avisos)) message(paste(avisos, collapse = "\n"))
  out
}
