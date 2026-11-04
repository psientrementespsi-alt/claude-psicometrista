# ==============================================================================
# est_padroes.R  (psicometria-r 2.1.0)
# Definição de ponto de corte por julgamento (standard setting): Angoff
# modificado e Bookmark (Cook & Pitoniak, 2026, cap. 12, Tabela 12.2, pp.
# 862-863). Angoff: recomendação = soma dos julgamentos de item de cada juiz,
# agregada entre juízes, com erro-padrão entre juízes. Bookmark: itens
# ordenados pela localização na escala da TRI sob um critério de probabilidade
# de resposta (RP); o marcador de cada juiz vira um theta de corte, convertido
# em escore esperado pela curva característica do teste.
# ==============================================================================

#' Ponto de corte pelo Angoff modificado
#'
#' @param julgamentos matriz ou data.frame juízes x itens com a probabilidade (0 a 1) de o examinando limítrofe
#'   acertar cada item; em item politômico, o escore esperado do limítrofe (entre 0 e o máximo do item)
#' @param max_item vetor com o escore máximo de cada item (padrão 1 = dicotômico)
#' @param nivel nível do IC entre juízes
corte_angoff <- function(julgamentos, max_item = NULL, nivel = .95, decimal = ".") {
  J <- as.matrix(julgamentos); k <- ncol(J); nj <- nrow(J)
  if (is.null(max_item)) max_item <- rep(1, k)
  if (length(max_item) != k) stop("max_item precisa ter um valor por item.", call. = FALSE)
  if (any(J < 0, na.rm = TRUE) || any(sweep(J, 2, max_item, ">"), na.rm = TRUE)) stop("julgamento fora de [0, máximo do item].", call. = FALSE)
  if (anyNA(J)) stop("Há julgamentos ausentes: o corte do juiz seria a soma de itens diferentes. Completar ou excluir o juiz.", call. = FALSE)
  if (nj < 2) stop("Angoff precisa de pelo menos 2 juízes para o erro-padrão.", call. = FALSE)
  por_juiz <- rowSums(J); m <- mean(por_juiz); dp <- stats::sd(por_juiz); ep <- dp / sqrt(nj)
  q <- stats::qt(1 - (1 - nivel) / 2, nj - 1)
  tab <- data.frame(estatistica = c("corte (média dos juízes)", "mediana dos juízes", "DP entre juízes", "EP do corte", "IC inferior", "IC superior", "mínimo", "máximo", "escore máximo do teste"),
                    valor = c(m, stats::median(por_juiz), dp, ep, m - q * ep, m + q * ep, min(por_juiz), max(por_juiz), sum(max_item)))
  itens <- data.frame(item = if (is.null(colnames(J))) paste0("item", seq_len(k)) else colnames(J),
                      media = colMeans(J), dp = apply(J, 2, stats::sd), stringsAsFactors = FALSE)
  apa <- sprintf("Pelo método de Angoff modificado, %d juízes recomendaram, em média, um corte de %s pontos em %s possíveis (DP = %s; EP = %s; %s).",
                 nj, .pr_num(m, 2, TRUE, decimal), .pr_num(sum(max_item), 0, TRUE, decimal), .pr_num(dp, 2, TRUE, decimal),
                 .pr_num(ep, 2, TRUE, decimal), .ee_ic(m - q * ep, m + q * ep, 2, TRUE, decimal, nivel))
  .ee_resultado("Ponto de corte (Angoff modificado)", n = nj, tabela = tab, posthoc = itens,
    decisoes = c("Recomendação de cada juiz = soma dos julgamentos de item (Cook & Pitoniak, 2026, Tabela 12.2).",
                 sprintf("Agregação pela média entre juízes, EP = DP / sqrt(%d) e IC t com %d gl; a mediana sai ao lado.", nj, nj - 1),
                 "A tabela de itens (posthoc) mostra a média e o DP de cada item: item com DP alto entre juízes é candidato à discussão da rodada seguinte."),
    apa = apa, extra = list(corte_por_juiz = por_juiz))
}

.bm_params <- function(itens) {
  if (inherits(itens, "SingleGroupClass")) {
    if (itens@Model$nfact != 1) stop("Bookmark exige modelo unidimensional.", call. = FALSE)
    co <- mirt::coef(itens, IRTpars = TRUE, simplify = TRUE)$items
    if (!all(c("a", "b") %in% colnames(co)) || ncol(co) > 4) stop("Bookmark aqui só aceita itens dicotômicos (2PL, 3PL ou 4PL).", call. = FALSE)
    data.frame(item = rownames(co), a = co[, "a"], b = co[, "b"], c = if ("g" %in% colnames(co)) co[, "g"] else 0,
               d = if ("u" %in% colnames(co)) co[, "u"] else 1, stringsAsFactors = FALSE)
  } else {
    p <- as.data.frame(itens); if (!all(c("item", "a", "b") %in% names(p))) stop("Passe um ajuste do mirt ou data.frame com item, a, b (e c, d opcionais).", call. = FALSE)
    if (is.null(p$c)) p$c <- 0; if (is.null(p$d)) p$d <- 1; p
  }
}

#' Ponto de corte pelo Bookmark
#'
#' @param itens ajuste unidimensional dicotômico do mirt, ou data.frame com item, a, b, c, d (parametrização IRT, D = 1)
#' @param marcadores vetor com a página do caderno ordenado (1 = item mais fácil) em que cada juiz pôs o marcador
#' @param rp critério de probabilidade de resposta (o capítulo exige um critério e não fixa o valor; .67 é uma
#'   convenção usual, a declarar)
#' @param regra "item_marcado" (theta do item da página marcada) ou "ponto_medio" (entre a página anterior e a marcada)
#' @param agregacao "mediana" ou "media" dos thetas de corte dos juízes
corte_bookmark <- function(itens, marcadores, rp = .67, regra = c("item_marcado", "ponto_medio"), agregacao = c("mediana", "media"),
                           decimal = ".") {
  regra <- match.arg(regra); agregacao <- match.arg(agregacao)
  p <- .bm_params(itens)
  if (any(rp <= p$c | rp >= p$d)) stop("rp precisa ficar entre a assíntota inferior e a superior de todos os itens.", call. = FALSE)
  p$local <- p$b + log((rp - p$c) / (p$d - rp)) / p$a
  oib <- p[order(p$local), ]; oib$pagina <- seq_len(nrow(oib)); rownames(oib) <- NULL
  if (any(marcadores < 1 | marcadores > nrow(oib) | marcadores != round(marcadores))) stop("marcador fora das páginas do caderno.", call. = FALSE)
  th <- vapply(marcadores, function(k) if (regra == "item_marcado" || k == 1) oib$local[k] else mean(oib$local[c(k - 1, k)]), numeric(1))
  corte <- if (agregacao == "mediana") stats::median(th) else mean(th)
  tcc <- function(t) sum(p$c + (p$d - p$c) * stats::plogis(p$a * (t - p$b)))
  tab <- data.frame(estatistica = c("theta de corte", "escore esperado no corte", "DP dos thetas dos juízes", "EP do theta (DP / sqrt(J))", "escore máximo"),
                    valor = c(corte, tcc(corte), stats::sd(th), stats::sd(th) / sqrt(length(th)), nrow(p)))
  apa <- sprintf("Pelo método Bookmark (RP = %s), %d juízes chegaram a um theta de corte de %s (%s dos juízes), que corresponde a um escore esperado de %s em %d.",
                 .pr_num(rp, 2, FALSE, decimal), length(th), .pr_num(corte, 2, TRUE, decimal), agregacao, .pr_num(tcc(corte), 1, TRUE, decimal), nrow(p))
  .ee_resultado("Ponto de corte (Bookmark)", n = length(th), tabela = tab, posthoc = oib[, c("pagina", "item", "a", "b", "c", "d", "local")],
    decisoes = c(sprintf("Localização de cada item = theta em que P = %s; caderno ordenado do mais fácil ao mais difícil (Cook & Pitoniak, 2026, Tabela 12.2).", rp),
                 sprintf("Theta de corte de cada juiz pela regra '%s'; agregação pela %s; escore pela curva característica do teste.", regra, agregacao),
                 "O RP muda a ordem e o corte: declarar o valor e, se possível, mostrar o corte com RP = .50 ao lado."),
    apa = apa, extra = list(theta_por_juiz = th))
}
