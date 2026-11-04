# ==============================================================================
# est_bayes.R  (psicometria-r 2.1.0)
# Fator de Bayes para teste t (JZS; Rouder et al., 2009) e para correlação de
# Pearson (prior beta esticada; Ly, Verhagen & Wagenmakers, 2016), em R base por
# integração numérica, com os priors padrão do JASP (Cauchy r = sqrt(2)/2;
# beta esticada kappa = 1). Quando o pacote BayesFactor está instalado, os evals
# comparam com ttestBF() e correlationBF(); sem ele, conferem contra uma
# estimativa Monte Carlo independente da verossimilhança marginal.
# ==============================================================================

.bf_jzs <- function(t, N, nu, r) {
  lden <- -(nu + 1) / 2 * log1p(t^2 / nu)
  f <- function(g) {
    a <- 1 + N * g * r^2
    exp(-0.5 * log(a) - (nu + 1) / 2 * log1p(t^2 / (a * nu)) - lden -
          0.5 * log(2 * pi) - 1.5 * log(g) - 1 / (2 * g))
  }
  I <- stats::integrate(f, 0, Inf, rel.tol = 1e-10, subdivisions = 1000L)
  I$value
}

.bf_rotulo <- function(bf10, idioma = "pt") {
  b <- if (bf10 >= 1) bf10 else 1 / bf10
  forca <- if (b < 3) "anedótica" else if (b < 10) "moderada" else if (b < 30) "forte" else if (b < 100) "muito forte" else "extrema"
  lado <- if (bf10 >= 1) "a favor de H1" else "a favor de H0"
  paste("evidência", forca, lado)
}

.bf_num <- function(x, decimal = ".") {
  s <- if (x >= 1000 || x < 0.001) formatC(x, format = "e", digits = 2) else formatC(x, format = "f", digits = 2)
  if (decimal == ",") s <- chartr(".", ",", s)
  s
}

#' Fator de Bayes para teste t (uma amostra, pareado ou independente)
#'
#' @param dados data.frame
#' @param vd variável dependente
#' @param grupo fator de 2 níveis (teste independente) ou NULL
#' @param pareado_com segunda medida (teste pareado: vd - pareado_com) ou NULL
#' @param mu valor de referência no teste de uma amostra
#' @param r_escala escala da Cauchy do tamanho de efeito (padrão do JASP sqrt(2)/2)
#' @param alternativa "bilateral" (H1: delta != 0). As unilaterais não estão implementadas.
bf_teste_t <- function(dados, vd, grupo = NULL, pareado_com = NULL, mu = 0, r_escala = sqrt(2) / 2,
                       decimal = ".", idioma = "pt") {
  if (!is.null(grupo) && !is.null(pareado_com)) stop("Use grupo (independente) ou pareado_com (pareado), não os dois.", call. = FALSE)
  decisoes <- character()
  if (!is.null(grupo)) {
    d <- dados[stats::complete.cases(dados[, c(vd, grupo)]), c(vd, grupo)]
    g <- droplevels(factor(d[[grupo]])); if (nlevels(g) != 2) stop("grupo precisa ter 2 níveis.", call. = FALSE)
    x1 <- d[[vd]][g == levels(g)[1]]; x2 <- d[[vd]][g == levels(g)[2]]
    tt <- stats::t.test(x1, x2, var.equal = TRUE); n1 <- length(x1); n2 <- length(x2)
    N <- n1 * n2 / (n1 + n2); nu <- n1 + n2 - 2; n <- n1 + n2; analise <- "Teste t bayesiano de amostras independentes"
    decisoes <- c(decisoes, sprintf("Diferença = %s - %s; t de Student (variâncias iguais), como o ttestBF.", levels(g)[1], levels(g)[2]))
  } else {
    cols <- c(vd, pareado_com)
    d <- dados[stats::complete.cases(dados[, cols, drop = FALSE]), cols, drop = FALSE]
    x <- if (is.null(pareado_com)) d[[vd]] - mu else d[[vd]] - d[[pareado_com]]
    tt <- stats::t.test(x); n <- length(x); N <- n; nu <- n - 1
    analise <- if (is.null(pareado_com)) "Teste t bayesiano de uma amostra" else "Teste t bayesiano pareado"
  }
  t <- unname(tt$statistic)
  escalas <- c(r_escala, 1, sqrt(2)); nomes <- c("informado", "amplo (r = 1)", "ultra-amplo (r = sqrt(2))")
  if (isTRUE(all.equal(r_escala, 1))) { escalas <- c(1, sqrt(2) / 2, sqrt(2)); nomes <- c("informado", "médio (r = .707)", "ultra-amplo") }
  bfs <- vapply(escalas, function(r) .bf_jzs(t, N, nu, r), numeric(1))
  tab <- data.frame(prior = nomes, r_escala = round(escalas, 3), BF10 = bfs, BF01 = 1 / bfs, log_BF10 = log(bfs), stringsAsFactors = FALSE)
  decisoes <- c(decisoes,
    sprintf("Prior Cauchy(0, %.3f) sobre delta, padrão do JASP para r = sqrt(2)/2 (Rouder et al., 2009); H1 bilateral.", r_escala),
    "Robustez: o BF é recalculado com r = 1 e r = sqrt(2); a conclusão vale se o lado e a força não mudarem.",
    "Rótulos de força (anedótica < 3, moderada < 10, forte < 30, muito forte < 100, extrema) são régua de leitura, não teste.")
  apa <- sprintf("%s, %s, %s, %s (%s).", analise, .ee_estat("t", t, nu, decimal), paste0("BF₁₀ = ", .bf_num(bfs[1], decimal)),
                 paste0("BF₀₁ = ", .bf_num(1 / bfs[1], decimal)), .bf_rotulo(bfs[1]))
  .ee_resultado(analise, n = n, tabela = tab, efeito = data.frame(indice = "BF10", valor = bfs[1], principal = TRUE),
                decisoes = decisoes, apa = apa, extra = list(t = t, gl = nu, N_efetivo = N))
}

# Densidade exata de r dado rho (Hotelling, forma integral) até uma constante em rho
.bf_lik_rho <- function(r, rho, n) {
  a <- rho * r
  # log(cosh(w) - a) sem overflow: w - log 2 + log1p(exp(-2w) - 2a exp(-w))
  f <- function(w) exp(-(n - 1) * (w - log(2) + log1p(exp(-2 * w) - 2 * a * exp(-w))))
  I <- stats::integrate(f, 0, Inf, rel.tol = 1e-10)$value
  (1 - rho^2)^((n - 1) / 2) * I
}

#' Fator de Bayes para correlação de Pearson (prior beta esticada em [-1, 1])
#'
#' @param kappa largura da prior (1 = uniforme, padrão do JASP)
bf_correlacao <- function(dados, x, y, kappa = 1, decimal = ".", idioma = "pt") {
  d <- dados[stats::complete.cases(dados[, c(x, y)]), c(x, y)]; n <- nrow(d)
  if (n < 4) stop("Correlação bayesiana exige pelo menos 4 pares completos.", call. = FALSE)
  r <- stats::cor(d[[x]], d[[y]])
  l0 <- .bf_lik_rho(r, 0, n)
  um <- function(k) {
    a <- 1 / k
    pri <- function(rho) stats::dbeta((rho + 1) / 2, a, a) / 2
    num <- stats::integrate(function(rr) vapply(rr, function(p) pri(p) * .bf_lik_rho(r, p, n) / l0, numeric(1)),
                            -1, 1, rel.tol = 1e-8, subdivisions = 500L)$value
    num
  }
  ks <- unique(c(kappa, 1, 0.5)); bfs <- vapply(ks, um, numeric(1))
  tab <- data.frame(kappa = ks, BF10 = bfs, BF01 = 1 / bfs, log_BF10 = log(bfs))
  apa <- sprintf("Correlação bayesiana de Pearson, %s, %s, %s (%s).", paste0("r(", n - 2, ") = ", .pr_num(r, 2, FALSE, decimal)),
                 paste0("BF₁₀ = ", .bf_num(bfs[1], decimal)), paste0("BF₀₁ = ", .bf_num(1 / bfs[1], decimal)), .bf_rotulo(bfs[1]))
  .ee_resultado("Correlação bayesiana de Pearson", n = n, tabela = tab, efeito = data.frame(indice = "BF10", valor = bfs[1], principal = TRUE),
    decisoes = c(sprintf("Prior beta esticada com largura kappa = %s sobre rho (kappa = 1 é a uniforme, padrão do JASP; Ly et al., 2016); H1 bilateral.", kappa),
                 "Verossimilhança exata de r (forma integral de Hotelling), integrada numericamente.",
                 "Robustez com kappa = 1 e kappa = .5 na tabela; rótulos de força são régua de leitura."),
    apa = apa, extra = list(r = r))
}
