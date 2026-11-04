# ==============================================================================
# est_rasch_cml.R  (psicometria-r 2.2.0)
# Modelo de Rasch dicotômico por máxima verossimilhança condicional (CML), a
# rota do eRm usada no curso (Rasch dicotômico - eRm.R): dificuldades com soma
# zero e EP, habilidade por ML para escores não extremos, infit e outfit (MNSQ e
# t), teste da razão de verossimilhança de Andersen, teste de Wald por item e
# PCA dos resíduos padronizados. Em R base; o eRm não está instalado no Mac.
# Os evals reproduzem a Tabela 1 do relatório do curso (BLOT, N = 150).
# ==============================================================================

# funções simétricas elementares gamma_0..gamma_k de eps (recursão de soma)
.rc_gamma <- function(eps) {
  g <- c(1, numeric(length(eps)))
  for (e in eps) g <- c(g[1], g[-1] + e * g[-length(g)])
  g
}

.rc_cml <- function(X) {
  X <- X[rowSums(X) > 0 & rowSums(X) < ncol(X), , drop = FALSE]
  k <- ncol(X); r <- rowSums(X); nr <- tabulate(r, nbins = k - 1); s <- colSums(X)
  ll <- function(th) { d <- c(th, -sum(th)); g <- .rc_gamma(exp(-d)); -sum(s * d) - sum(nr * log(g[2:k])) }
  gr <- function(th) {
    d <- c(th, -sum(th)); eps <- exp(-d); g <- .rc_gamma(eps)
    gd <- vapply(seq_len(k), function(i) { gi <- .rc_gamma(eps[-i]); sum(nr * eps[i] * gi[1:(k - 1)] / g[2:k]) }, 0)
    full <- -s + gd          # d logL / d d_i
    -full[-k] + full[k]      # gradiente de -logL em th, com d_k = -sum(th)
  }
  o <- stats::optim(rep(0, k - 1), function(th) -ll(th), function(th) gr(th), method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-14))
  H <- stats::optimHess(o$par, function(th) -ll(th), function(th) gr(th))
  V <- solve(H); d <- c(o$par, -sum(o$par))
  A <- rbind(diag(k - 1), rep(-1, k - 1)); Vd <- A %*% V %*% t(A)
  list(dificuldade = d, ep = sqrt(diag(Vd)), V = Vd, logLik = -o$value, convergiu = o$convergence == 0, n = nrow(X))
}

.rc_theta <- function(d, k) {
  vapply(1:(k - 1), function(r) stats::uniroot(function(t) sum(stats::plogis(t - d)) - r, c(-15, 15), tol = 1e-12)$root, 0)
}

#' Rasch dicotômico por CML (rota do eRm)
#'
#' @param dados data.frame com itens 0/1
#' @param itens nomes dos itens
#' @param divisao "mediana" (padrão do LRtest do eRm: escore > mediana contra os demais), "media", ou vetor de grupos
#'   para o teste de Andersen e o de Wald
rasch_cml <- function(dados, itens, divisao = "mediana", decimal = ".") {
  X <- as.matrix(dados[stats::complete.cases(dados[, itens]), itens]); n0 <- nrow(dados)
  if (!all(X %in% c(0, 1))) stop("rasch_cml() aceita só itens 0/1.", call. = FALSE)
  k <- ncol(X); est <- .rc_cml(X); d <- est$dificuldade
  r <- rowSums(X); nao_ext <- r > 0 & r < k
  th_r <- .rc_theta(d, k); th <- rep(NA_real_, nrow(X)); th[nao_ext] <- th_r[r[nao_ext]]
  Xf <- X[nao_ext, , drop = FALSE]; P <- stats::plogis(outer(th[nao_ext], d, "-")); W <- P * (1 - P)
  C <- P * (1 - P)^4 + (1 - P) * P^4; Z2 <- (Xf - P)^2 / W; nf <- nrow(Xf)
  outfit <- colMeans(Z2); infit <- colSums((Xf - P)^2) / colSums(W)
  q_out <- sqrt(pmax(colSums(C / W^2) / nf^2 - 1 / nf, 1e-12)); q_in <- sqrt(pmax(colSums(C - W^2) / colSums(W)^2, 1e-12))
  t_out <- (outfit^(1 / 3) - 1) * (3 / q_out) + q_out / 3; t_in <- (infit^(1 / 3) - 1) * (3 / q_in) + q_in / 3
  itens_tab <- data.frame(item = itens, dificuldade = d, ep = est$ep, outfit_msq = outfit, infit_msq = infit,
                          outfit_t = t_out, infit_t = t_in, row.names = NULL)
  # Andersen LR e Wald
  grp <- if (identical(divisao, "mediana")) ifelse(r > stats::median(r), "alto", "baixo") else
    if (identical(divisao, "media")) ifelse(r > mean(r), "alto", "baixo") else as.character(divisao[stats::complete.cases(dados[, itens])])
  gidx <- split(seq_len(nrow(X)), grp)
  # item com só acertos ou só erros (entre quem não tem escore extremo) num subgrupo sai do teste, como no LRtest do eRm
  extremo_em <- function(Xg, cols) { Xg <- Xg[, cols, drop = FALSE]; rr <- rowSums(Xg); Xg <- Xg[rr > 0 & rr < ncol(Xg), , drop = FALSE]
    cs <- colSums(Xg); cols[cs == 0 | cs == nrow(Xg)] }
  usa <- itens; repeat { fora <- unique(unlist(lapply(gidx, function(ix) extremo_em(X[ix, , drop = FALSE], usa)))); if (!length(fora)) break; usa <- setdiff(usa, fora) }
  avisos <- character(); lr <- NA_real_; gl <- NA_integer_; wald <- NULL
  if (length(usa) < k) avisos <- c(avisos, sprintf("Teste de Andersen e Wald sem %d item(ns) com só acertos ou só erros num subgrupo: %s (o LRtest do eRm faz a mesma exclusão).",
                                                k - length(usa), paste(setdiff(itens, usa), collapse = ", ")))
  if (length(usa) >= 3) {
    tot <- .rc_cml(X[, usa, drop = FALSE]); fits <- lapply(gidx, function(ix) .rc_cml(X[ix, usa, drop = FALSE]))
    lr <- 2 * (sum(vapply(fits, `[[`, 0, "logLik")) - tot$logLik); gl <- (length(usa) - 1) * (length(fits) - 1)
    if (length(fits) == 2) {
      dz <- (fits[[1]]$dificuldade - fits[[2]]$dificuldade) / sqrt(fits[[1]]$ep^2 + fits[[2]]$ep^2)
      wald <- data.frame(item = usa, z = dz, p = 2 * stats::pnorm(-abs(dz)))
    }
  }
  # PCA dos resíduos padronizados
  Zs <- (Xf - P) / sqrt(W); ev <- eigen(stats::cor(Zs), symmetric = TRUE, only.values = TRUE)$values
  # confiabilidade de separação de pessoas (variância de theta menos erro médio)
  ep_th <- 1 / sqrt(rowSums(W)); vt <- stats::var(th[nao_ext]); sep_rel <- (vt - mean(ep_th^2)) / vt
  tab <- data.frame(estatistica = c("N (casos completos)", "N sem escore extremo", "Andersen LR", "gl", "p", "1º autovalor dos resíduos", "confiabilidade de separação de pessoas"),
                    valor = c(nrow(X), nf, lr, gl, if (is.na(lr)) NA else stats::pchisq(lr, gl, lower.tail = FALSE), ev[1], sep_rel))
  apa <- sprintf("O modelo de Rasch (CML) teve teste de Andersen LR(%s) = %s, p %s; %d itens com infit ou outfit MNSQ fora de %s e o primeiro autovalor dos resíduos foi %s.",
                 gl, .pr_num(lr, 2, TRUE, decimal), if (is.na(lr)) "não calculado" else { p <- stats::pchisq(lr, gl, lower.tail = FALSE); if (p < .001) "< .001" else paste("=", sub("^0", "", sprintf("%.3f", p))) },
                 sum(outfit < .7 | outfit > 1.3 | infit < .7 | infit > 1.3), if (decimal == ",") "[0,7; 1,3]" else "[0.7, 1.3]", .pr_num(ev[1], 2, TRUE, decimal))
  .ee_resultado("Modelo de Rasch dicotômico (CML)", n = nrow(X), n_excluidos = n0 - nrow(X), tabela = tab, posthoc = itens_tab,
    decisoes = c("Máxima verossimilhança condicional, dificuldades com soma zero, como RM() do eRm (sum0 = TRUE); dificuldade = -betapar.",
                 "Habilidade por ML para escores de 1 a k-1; quem acerta tudo ou nada fica fora do ajuste e dos índices de infit e outfit, como no eRm.",
                 "Infit e outfit MNSQ com t pela transformação cúbica de Wright e Masters, como itemfit() do eRm; faixa 0,7 a 1,3 é régua (Bond et al., 2020, citado no relatório do curso).",
                 sprintf("Andersen LR com divisão '%s' (padrão do LRtest: mediana do escore); Wald por item na mesma divisão.", if (is.character(divisao) && length(divisao) == 1) divisao else "grupos informados"),
                 "PCA dos resíduos padronizados: primeiro autovalor perto de 2 ou mais sugere uma segunda dimensão (régua usual, a declarar)."),
    avisos = avisos, apa = apa, extra = list(wald = wald, autovalores_residuos = ev, theta_por_escore = stats::setNames(th_r, 1:(k - 1))))
}
