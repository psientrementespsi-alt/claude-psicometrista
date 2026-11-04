# ==============================================================================
# ajuste_pessoa.R  (higiene-dados v1.5.0)
# Ajuste de pessoa para triagem de respostas descuidadas em escalas politomicas
# de tamanho medio (Niessen, Meijer & Tendeiro, 2016, J. Res. Pers. 63, 1-11):
# - G^p, numero de erros de Guttman com passos de item (Apendice A do artigo);
# - lz^p, log-verossimilhanca padronizada sob o modelo de resposta graduada
#   (Drasgow, Levine & Williams, 1985, citado no artigo), via mirt;
# - extensao multiteste: soma de G^p e de lz^p nas subescalas (secao 3.3).
# O artigo usa cortes empiricos (5%, 10% ou 20% mais extremos da amostra); por
# isso o indicador entra como UMA familia no criterio de 2 ou mais familias da
# triagem_cier_completa(), nunca sozinho. O mesmo artigo pede relato com e sem
# os casos excluidos.
# Privacidade: devolve vetores por pessoa para uso interno e um resumo agregado.
# ==============================================================================

#' Erros de Guttman politomicos (G^p) por pessoa
#'
#' @param dados data.frame ou matriz com os itens de UMA escala (mesmas categorias)
#' @param itens nomes ou indices das colunas
#' @param categorias categorias ordenadas da escala (padrao: valores observados)
#' @param p_passos popularidade dos passos, na ordem item 1 passo 1, item 1 passo 2, ...;
#'   NULL calcula da amostra como no Apendice A: tomaram o passo / tomaram o passo anterior
#' @return vetor de G^p (NA para quem tem item ausente), com atributos ordem_passos e p_passos
calc_guttman_poly <- function(dados, itens = colnames(dados), categorias = NULL, p_passos = NULL) {
  M <- as.matrix(as.data.frame(dados)[, itens, drop = FALSE]); storage.mode(M) <- "double"
  k <- ncol(M); if (is.null(categorias)) categorias <- sort(unique(as.vector(M[!is.na(M)])))
  s <- length(categorias) - 1L
  if (s < 1) stop("calc_guttman_poly: a escala precisa de pelo menos 2 categorias.", call. = FALSE)
  if (any(!is.na(M) & !(M %in% categorias))) stop("calc_guttman_poly: ha valores fora de 'categorias'.", call. = FALSE)
  nm_it <- if (is.character(itens)) itens else colnames(M)
  if (is.null(nm_it)) nm_it <- paste0("i", seq_len(k))
  Y <- do.call(cbind, lapply(seq_len(k), function(j) vapply(categorias[-1], function(cc) as.numeric(M[, j] >= cc), numeric(nrow(M)))))
  Y <- matrix(Y, nrow = nrow(M)); nome <- paste0(rep(nm_it, each = s), "-", rep(seq_len(s), k)); colnames(Y) <- nome
  completo <- stats::complete.cases(M)
  if (is.null(p_passos)) {
    p_passos <- unlist(lapply(seq_len(k), function(j) {
      tom <- colSums(Y[completo, (j - 1) * s + seq_len(s), drop = FALSE])
      ant <- c(sum(completo), tom[-s]); ifelse(ant > 0, tom / ant, 0)
    }))
  } else if (length(p_passos) != ncol(Y)) stop("calc_guttman_poly: p_passos precisa de ", ncol(Y), " valores.", call. = FALSE)
  names(p_passos) <- nome
  ordem <- order(-p_passos, seq_along(p_passos))   # empate: ordem dos itens e dos passos
  Yo <- Y[, ordem, drop = FALSE]; Yo[is.na(Yo)] <- 0
  nao_antes <- t(apply(1 - Yo, 1, function(r) c(0, cumsum(r)[-length(r)])))
  if (ncol(Yo) == 1) nao_antes <- matrix(0, nrow(Yo), 1)
  G <- rowSums(Yo * nao_antes); G[!completo] <- NA
  attr(G, "ordem_passos") <- nome[ordem]; attr(G, "p_passos") <- p_passos
  G
}

#' lz politomico (lz^p) por pessoa sob o modelo de resposta graduada
#'
#' @param theta escores da pessoa; NULL estima por EAP no proprio ajuste
#' @return vetor de lz^p (valores baixos = padrao improvavel; NA para quem tem item ausente)
calc_lz_poly <- function(dados, itens = colnames(dados), theta = NULL, metodo_theta = "EAP") {
  if (!requireNamespace("mirt", quietly = TRUE)) stop("calc_lz_poly precisa do pacote mirt.", call. = FALSE)
  D <- as.data.frame(as.data.frame(dados)[, itens, drop = FALSE]); completo <- stats::complete.cases(D)
  fit <- mirt::mirt(D, 1, itemtype = "graded", verbose = FALSE)
  if (is.null(theta)) theta <- mirt::fscores(fit, method = metodo_theta, verbose = FALSE)[, 1]
  n <- nrow(D); l0 <- E <- V <- numeric(n)
  for (j in seq_len(ncol(D))) {
    P <- mirt::probtrace(mirt::extract.item(fit, j), matrix(theta)); lp <- log(P)
    x <- match(D[[j]], sort(unique(D[[j]][!is.na(D[[j]])])))
    l0 <- l0 + ifelse(is.na(x), 0, lp[cbind(seq_len(n), ifelse(is.na(x), 1L, x))])
    e <- rowSums(P * lp); E <- E + e; V <- V + rowSums(P * lp^2) - e^2
  }
  lz <- (l0 - E) / sqrt(V); lz[!completo] <- NA
  lz
}

#' Ajuste de pessoa multiteste (G^p_m e lz^p_m) e indicador para a triagem C/IER
#'
#' @param subescalas lista nomeada de vetores de itens; uma escala so = lista de um elemento
#' @param estatistica "Gp" (nao parametrico, padrao) ou "lzp" (exige mirt) para o indicador
#' @param quantil corte empirico: G^p_m acima do quantil, ou lz^p_m abaixo de 1 - quantil
#' @return data.frame por pessoa (Gp_m, lzp_m, flag_ajuste_pessoa), com atributos decisoes e resumo
calc_ajuste_pessoa <- function(dados, subescalas, estatistica = c("Gp", "lzp"), quantil = 0.95, categorias = NULL, calcular_lz = NULL) {
  estatistica <- match.arg(estatistica)
  if (!is.list(subescalas)) subescalas <- list(escala = subescalas)
  if (is.null(calcular_lz)) calcular_lz <- estatistica == "lzp" || requireNamespace("mirt", quietly = TRUE)
  avisos <- character()
  if (any(lengths(subescalas) < 5)) avisos <- c(avisos, "subescala com menos de 5 itens: o ajuste de pessoa perde sensibilidade em escalas curtas (Niessen et al., 2016, secao 3.3)")
  Gm <- rowSums(do.call(cbind, lapply(subescalas, function(it) calc_guttman_poly(dados, it, categorias))))
  lzm <- if (isTRUE(calcular_lz)) rowSums(do.call(cbind, lapply(subescalas, function(it) calc_lz_poly(dados, it)))) else rep(NA_real_, length(Gm))
  if (estatistica == "Gp") { corte <- stats::quantile(Gm, quantil, na.rm = TRUE, names = FALSE); flag <- !is.na(Gm) & Gm > corte }
  else { corte <- stats::quantile(lzm, 1 - quantil, na.rm = TRUE, names = FALSE); flag <- !is.na(lzm) & lzm < corte }
  out <- data.frame(Gp_m = Gm, lzp_m = round(lzm, 4), flag_ajuste_pessoa = flag)
  attr(out, "decisoes") <- sprintf("ajuste de pessoa: %s somado em %d subescala(s); corte empirico %s %s (quantil %.2f da amostra; Niessen et al., 2016); entra como uma familia na triagem, nunca sozinho",
    if (estatistica == "Gp") "G^p" else "lz^p", length(subescalas), if (estatistica == "Gp") ">" else "<", format(round(corte, 3)), if (estatistica == "Gp") quantil else 1 - quantil)
  attr(out, "resumo") <- data.frame(n_calculaveis = sum(!is.na(if (estatistica == "Gp") Gm else lzm)), n_assinalados = sum(flag), corte = corte,
    r_Gp_lzp = if (isTRUE(calcular_lz)) round(stats::cor(Gm, lzm, use = "complete.obs"), 3) else NA)
  attr(out, "avisos") <- avisos; if (length(avisos)) message(paste(avisos, collapse = "\n"))
  out
}
