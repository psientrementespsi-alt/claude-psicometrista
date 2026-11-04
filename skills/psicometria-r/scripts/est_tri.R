# est_tri.R: Teoria de Resposta ao Item da skill psicometria-r (v1.4; multidimensional na v1.5).
# Carregar por scripts/carregar.R (usa .ee_resultado e a formatação de est_comum.R e
# .pr_num, .pr_p e .pr_exige de psicometria.R).
#
# Fontes: curso Psicometria Online, módulos "Introdução à TRI" (Winsteps) e "R: Teoria de
# Resposta ao Item" (scripts com mirt e eRm); Martinková e Hladká (2023, caps. 6 a 9),
# Petersen (2024, cap. 8) e Wiley Handbook (2018, caps. 15 a 21 e 27). Réguas, rotas e
# divergências entre as fontes: references/tri.md.
#
# Funções públicas
#   tri_ajustar(dados, itens, modelo, ausentes, regua_mnsq, q3_limite, corrigir, theta,
#               semente, decimal, idioma, ...)
#   tri_comparar(dados, itens, modelos, semente, decimal, idioma)
#   tri_dif(dados, itens, grupo, modelo, metodo, ancoras, purificar, corrigir, alfa,
#           semente, decimal, idioma)
#   tri_multidimensional(...) (v1.5); escala_mokken(dados, itens, limite, aisp,
#           monotonicidade, iio) e simular_cat(dados, itens, modelo, fit, ep_max,
#           max_itens, metodo, criterio, semente) (v1.7)
#
# Privacidade: nenhuma função devolve θ, resíduo ou ajuste de pessoa. Saem parâmetros de
# item, índices de ajuste, contagens e proporções.

.tri_t <- function(idioma, pt, en) if (identical(idioma, "en")) en else pt
.tri_modelos <- c("graded", "gpcm", "pcm", "rsm", "rasch", "2PL", "3PL")
.tri_itemtype <- function(modelo) switch(modelo, graded = "graded", gpcm = "gpcm", pcm = "Rasch",
                                         rsm = "rsm", rasch = "Rasch", "2PL" = "2PL", "3PL" = "3PL")
.tri_familia_rasch <- function(modelo) modelo %in% c("rasch", "pcm", "rsm")

.tri_nome <- function(modelo, idioma = "pt") {
  pt <- c(graded = "modelo de resposta graduada (Samejima, 1969)",
          gpcm = "modelo de crédito parcial generalizado (Muraki, 1992)",
          pcm = "modelo de crédito parcial (Masters, 1982)",
          rsm = "modelo de escala de avaliação (Andrich, 1978)",
          rasch = "modelo de Rasch para itens dicotômicos",
          "2PL" = "modelo logístico de dois parâmetros (2PL)",
          "3PL" = "modelo logístico de três parâmetros (3PL)")
  en <- c(graded = "graded response model (Samejima, 1969)",
          gpcm = "generalized partial credit model (Muraki, 1992)",
          pcm = "partial credit model (Masters, 1982)",
          rsm = "rating scale model (Andrich, 1978)",
          rasch = "dichotomous Rasch model",
          "2PL" = "two-parameter logistic model (2PL)",
          "3PL" = "three-parameter logistic model (3PL)")
  (if (identical(idioma, "en")) en else pt)[[modelo]]
}

# Número de categorias observadas por item. O mirt recodifica os valores distintos,
# em ordem, para 0..K-1; categoria não usada some do item.
.tri_categorias <- function(d) vapply(d, function(x) length(unique(x[!is.na(x)])), integer(1))
.tri_recodifica <- function(x) match(x, sort(unique(x[!is.na(x)]))) - 1L

.tri_prepara <- function(dados, itens, ausentes) {
  if (!is.data.frame(dados)) stop("dados precisa ser um data.frame", call. = FALSE)
  falta <- setdiff(itens, names(dados))
  if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "), call. = FALSE)
  if (length(itens) < 3) stop("a TRI daqui pede pelo menos 3 itens", call. = FALSE)
  d <- as.data.frame(dados)[, itens, drop = FALSE]
  for (v in itens) if (!is.numeric(d[[v]])) stop("o item ", v, " precisa ser numérico (categorias como números inteiros)", call. = FALSE)
  nao_int <- itens[vapply(d, function(x) any(!is.na(x) & abs(x - round(x)) > 1e-8), logical(1))]
  if (length(nao_int)) stop("itens com valores não inteiros: ", paste(nao_int, collapse = ", "),
                            ". A TRI daqui é para respostas em categorias.", call. = FALSE)
  ok <- if (ausentes == "listwise") stats::complete.cases(d) else rowSums(!is.na(d)) > 0
  list(d = d[ok, , drop = FALSE], n_excluidos = sum(!ok), tem_na = anyNA(d[ok, , drop = FALSE]))
}

.tri_mirt <- function(d, itemtype, semente, grupo = NULL, invariancia = NULL, ...) {
  set.seed(semente)
  av <- character()
  chama <- function() {
    if (is.null(grupo)) mirt::mirt(d, 1, itemtype = itemtype, SE = TRUE, verbose = FALSE, ...)
    else mirt::multipleGroup(d, 1, group = grupo, itemtype = itemtype, invariance = invariancia,
                             SE = TRUE, verbose = FALSE, ...)
  }
  r <- withCallingHandlers(tryCatch(chama(), error = function(e) e),
                           warning = function(w) { av <<- c(av, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (inherits(r, "error")) return(list(fit = NULL, erro = conditionMessage(r), avisos = unique(av)))
  list(fit = r, erro = NULL, avisos = unique(av))
}

# Ajuste global: M2 com itens dicotômicos; M2* com politômicos (Cai & Hansen, 2013);
# C2 quando o M2* fica sem graus de liberdade (Cai & Monroe, 2014), como manda a
# documentação do mirt.
.tri_ajuste_global <- function(fit, dicotomico) {
  tenta <- function(tipo) tryCatch(suppressWarnings(mirt::M2(fit, type = tipo)), error = function(e) NULL)
  tipo <- if (dicotomico) "M2" else "M2*"
  r <- tenta(tipo)
  if (is.null(r) || !is.finite(r$df[1]) || r$df[1] <= 0) { tipo <- "C2"; r <- tenta("C2") }
  if (is.null(r)) return(NULL)
  data.frame(estatistica = tipo, valor = r$M2[1], gl = r$df[1], p = r$p[1], rmsea = r$RMSEA[1],
             rmsea_li = r$RMSEA_5[1], rmsea_ls = r$RMSEA_95[1], srmsr = r$SRMSR[1], cfi = r$CFI[1],
             tli = r$TLI[1], row.names = NULL)
}

# Parâmetros na parametrização da TRI (a, b), com erro-padrão pelo método delta do mirt.
# No modelo de escala de avaliação o mirt não converte: saem os limiares comuns (b) e a
# localização de cada item (c) na parametrização dele.
.tri_parametros <- function(fit, itens, modelo) {
  cf <- if (modelo == "rsm") NULL else tryCatch(mirt::coef(fit, IRTpars = TRUE, printSE = TRUE), error = function(e) NULL)
  irt <- !is.null(cf)
  if (!irt) cf <- mirt::coef(fit, printSE = TRUE)
  fixos <- switch(modelo, rasch = c("a", "a1", "g", "u"), pcm = c("a", "a1"), rsm = c("a1"),
                  "2PL" = c("g", "u"), "3PL" = c("u"), character())
  linhas <- lapply(itens, function(i) {
    m <- cf[[i]]
    nomes <- setdiff(colnames(m), fixos)
    v <- m[1, nomes]; e <- if (nrow(m) > 1) m[2, nomes] else rep(NA_real_, length(nomes))
    l <- as.list(c(rbind(v, e))); names(l) <- c(rbind(nomes, paste0(nomes, "_ep")))
    data.frame(item = i, l, check.names = FALSE, stringsAsFactors = FALSE)
  })
  todas <- unique(unlist(lapply(linhas, names)))
  linhas <- lapply(linhas, function(l) { for (cn in setdiff(todas, names(l))) l[[cn]] <- NA_real_; l[todas] })
  list(tabela = do.call(rbind, linhas), irt = irt)
}

# Limiares desordenados nos modelos de categorias adjacentes (PCM, GPCM): o passo k+1
# abaixo do passo k indica categoria que nunca é a mais provável em faixa alguma do traço.
.tri_desordenados <- function(par, modelo) {
  if (!modelo %in% c("pcm", "gpcm")) return(character())
  bs <- grep("^b[0-9]+$", names(par), value = TRUE)
  if (length(bs) < 2) return(character())
  B <- as.matrix(par[, bs, drop = FALSE])
  par$item[apply(B, 1, function(b) { b <- b[!is.na(b)]; length(b) > 1 && any(diff(b) < 0) })]
}

# PCA dos resíduos padronizados (pessoas com escore extremo fora, como no Winsteps):
# o maior autovalor da matriz de correlação dos resíduos, em unidades de item.
.tri_pca_residuos <- function(fit, d) {
  th <- tryCatch(suppressWarnings(mirt::fscores(fit, method = "WLE", verbose = FALSE)), error = function(e) NULL)
  if (is.null(th)) return(NA_real_)
  ok <- is.finite(th[, 1])
  if (sum(ok) < 10) return(NA_real_)
  z <- matrix(NA_real_, nrow(d), ncol(d))
  for (j in seq_len(ncol(d))) {
    it <- mirt::extract.item(fit, j)
    P <- mirt::probtrace(it, th[ok, 1, drop = FALSE])
    cats <- seq_len(ncol(P)) - 1
    E <- as.vector(P %*% cats); V <- as.vector(P %*% cats^2) - E^2
    x <- .tri_recodifica(d[[j]])[ok]
    z[ok, j] <- (x - E) / sqrt(pmax(V, 1e-10))
  }
  R <- suppressWarnings(stats::cor(z, use = "pairwise.complete.obs"))
  if (anyNA(R)) return(NA_real_)
  eigen(R, symmetric = TRUE, only.values = TRUE)$values[1]
}

# Confiabilidade marginal (Green et al., 1984): σ²/(σ² + E[EP²(θ)]), com a esperança na
# distribuição latente estimada N(μ, σ²). O marginal_rxx() do mirt 1.48 usa I/(I + σ²) e
# integra em N(0, 1): coincide quando σ² = 1 (GRM, GPCM, 2PL), mas distorce a família
# Rasch, em que a variância latente é estimada: infla com σ² < 1 (,94 contra ,77 num teste
# da construção) e encolhe com σ² > 1 (evals_tri.R: ,46 contra ,79 com σ² = 2,39, com a
# empírica EAP em ,80).
.tri_var_latente <- function(fit) {
  cf <- mirt::coef(fit, simplify = TRUE)
  c(media = as.vector(cf$means)[1], var = as.vector(cf$cov)[1])
}
.tri_rxx_marginal <- function(fit, pontos = 2001) {
  ml <- .tri_var_latente(fit); s <- sqrt(ml[["var"]])
  th <- seq(ml[["media"]] - 6 * s, ml[["media"]] + 6 * s, length.out = pontos)
  I <- mirt::testinfo(fit, matrix(th)); w <- stats::dnorm(th, ml[["media"]], s); w <- w / sum(w)
  ml[["var"]] / (ml[["var"]] + sum(w / I))
}

.tri_regua_mnsq <- function(regua, n) {
  if (regua == "linacre") return(c(.5, 1.5))
  if (n <= 500) c(.7, 1.3) else if (n <= 1000) c(.8, 1.2) else c(.9, 1.1)
}

tri_ajustar <- function(dados, itens, modelo = "graded", ausentes = c("listwise", "fiml"),
                        regua_mnsq = c("linacre", "amostra"), q3_limite = .20, corrigir = "BH",
                        theta = c(-4, 4), semente = 2026, decimal = ".", idioma = "pt", ...) {
  .pr_exige("mirt")
  modelo <- match.arg(modelo, .tri_modelos)
  ausentes <- match.arg(ausentes); regua_mnsq <- match.arg(regua_mnsq)
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  pp <- .tri_prepara(dados, itens, ausentes)
  d <- pp$d; n <- nrow(d)
  k <- .tri_categorias(d)
  dic <- all(k == 2)
  if (any(k < 2)) stop("item sem variação (uma categoria só): ", paste(itens[k < 2], collapse = ", "), call. = FALSE)
  if (modelo %in% c("rasch", "2PL", "3PL") && !dic)
    stop("o modelo ", modelo, " é para itens dicotômicos; para itens com mais categorias use graded, gpcm, pcm ou rsm.", call. = FALSE)
  if (modelo == "rsm" && length(unique(k)) > 1)
    stop("o modelo de escala de avaliação pede o mesmo número de categorias em todos os itens (observado: ",
         paste(sort(unique(k)), collapse = ", "), "). Juntar categorias vazias ou usar pcm.", call. = FALSE)
  decisoes <- character(); avisos <- character()
  if (modelo == "pcm" && dic) { modelo <- "rasch"; decisoes <- c(decisoes, "Itens dicotômicos: o PCM se reduz ao modelo de Rasch.") }
  # Tamanho amostral, régua do curso (Şahin & Anıl, 2017, para 1PL e 2PL; exemplos com o ENEM para 3PL).
  n_min <- switch(modelo, rasch = 150, pcm = 150, rsm = 150, "3PL" = 500, 200)
  if (n < n_min) avisos <- c(avisos, paste0("n = ", n, ": abaixo da régua do curso para o ", .tri_nome(modelo),
                                           " (", n_min, "; Şahin & Anıl, 2017, para 1PL e 2PL, com 500 como ideal fora do 1PL). Parâmetros e ajuste ficam instáveis; relatar como limitação."))
  if (any(k > 2) && length(unique(k)) > 1)
    avisos <- c(avisos, paste0("Itens com número diferente de categorias observadas (", paste(sort(unique(k)), collapse = ", "),
                               "): categoria vazia some do item no mirt. Conferir a frequência por categoria (auditar_itens()) antes de interpretar os limiares."))
  aj <- .tri_mirt(d, .tri_itemtype(modelo), semente, ...)
  if (is.null(aj$fit)) stop("o mirt não ajustou o modelo: ", aj$erro, call. = FALSE)
  fit <- aj$fit
  if (length(aj$avisos)) avisos <- c(avisos, paste0("mirt: ", aj$avisos))
  if (!isTRUE(mirt::extract.mirt(fit, "converged")))
    avisos <- c(avisos, "O EM não convergiu: aumentar o limite de ciclos (technical = list(NCYCLES = 2000)) ou simplificar o modelo antes de relatar.")
  decisoes <- c(decisoes, paste0(.tri_nome(modelo), ", unidimensional, por máxima verossimilhança marginal (mirt ",
                                 as.character(utils::packageVersion("mirt")), "; semente ", semente, ")."),
                if (ausentes == "listwise") paste0("Casos incompletos fora (", pp$n_excluidos, "). Alternativa: ausentes = \"fiml\", que usa todos os padrões de resposta na estimação; o S-X² e o Q3 continuam pedindo casos completos.")
                else "Ausentes tratados por máxima verossimilhança com informação completa; S-X², Q3 e confiabilidade empírica só com casos completos.")
  # Ajuste global
  glob <- .tri_ajuste_global(fit, dic)
  if (is.null(glob)) avisos <- c(avisos, "M2, M2* e C2 não puderam ser calculados (poucos graus de liberdade ou ausentes): ajuste global não calculado.")
  # Parâmetros
  pr <- .tri_parametros(fit, itens, modelo)
  par <- pr$tabela
  if (!pr$irt && modelo != "rsm") avisos <- c(avisos, "Parâmetros na parametrização inclinação-intercepto do mirt (a conversão para a TRI clássica falhou).")
  if (modelo == "rsm") decisoes <- c(decisoes, "Escala de avaliação: limiares b comuns a todos os itens e localização c por item, na parametrização do mirt (a do primeiro item fixa em 0).")
  if ("a" %in% names(par)) {
    a <- par$a
    if (any(a > 4, na.rm = TRUE)) avisos <- c(avisos, paste0("Discriminação muito alta (a > 4) em ", paste(par$item[which(a > 4)], collapse = ", "),
                                                         ": suspeitar de dependência local, item redundante ou amostra pequena antes de interpretar (no curso, a de 39 e 41 num banco de 150 casos sinalizou problema de ajuste)."))
    if (any(a < .4, na.rm = TRUE)) avisos <- c(avisos, paste0("Discriminação baixa (a < ,40, régua da skill) em ", paste(par$item[which(a < .4)], collapse = ", "),
                                                          ": o item quase não separa níveis do traço; negativo indica item invertido sem recodificar."))
  }
  deso <- .tri_desordenados(par, modelo)
  if (length(deso)) avisos <- c(avisos, paste0("Limiares desordenados em ", paste(deso, collapse = ", "),
                                               ": alguma categoria nunca é a mais provável. Conferir a frequência das categorias e considerar juntar as vizinhas, com justificativa de conteúdo."))
  # Ajuste por item
  sx2 <- tryCatch(suppressWarnings(mirt::itemfit(fit, fit_stats = "S_X2", na.rm = pp$tem_na)), error = function(e) NULL)
  ajuste_itens <- NULL
  if (!is.null(sx2)) {
    ajuste_itens <- data.frame(item = sx2$item, s_x2 = sx2$S_X2, gl = sx2$df.S_X2, rmsea = sx2$RMSEA.S_X2,
                               p = sx2$p.S_X2, p_ajust = stats::p.adjust(sx2$p.S_X2, corrigir), row.names = NULL)
  } else avisos <- c(avisos, "S-X² não calculado (o mirt pede casos completos e células suficientes).")
  mnsq <- NULL; lim <- .tri_regua_mnsq(regua_mnsq, n)
  if (.tri_familia_rasch(modelo)) {
    inf <- tryCatch(suppressWarnings(mirt::itemfit(fit, fit_stats = "infit", na.rm = pp$tem_na)), error = function(e) NULL)
    if (!is.null(inf)) {
      mnsq <- data.frame(item = inf$item, infit = inf$infit, z_infit = inf$z.infit, outfit = inf$outfit,
                         z_outfit = inf$z.outfit, row.names = NULL)
      mnsq$fora_regua <- mnsq$infit < lim[1] | mnsq$infit > lim[2] | mnsq$outfit < lim[1] | mnsq$outfit > lim[2]
      ajuste_itens <- if (is.null(ajuste_itens)) mnsq else merge(ajuste_itens, mnsq, by = "item", sort = FALSE)
    }
  }
  # Dependência local (Q3 de Yen, relativo à média dos Q3; Christensen et al., 2017)
  q3 <- tryCatch(mirt::residuals(fit, type = "Q3", verbose = FALSE), error = function(e) NULL)
  pares_q3 <- NULL; q3_max <- NA_real_; q3_media <- NA_real_
  if (!is.null(q3)) {
    off <- q3[upper.tri(q3)]; q3_media <- mean(off)
    rel <- q3 - q3_media
    q3_max <- max(rel[upper.tri(rel)])
    w <- which(upper.tri(rel) & rel > q3_limite, arr.ind = TRUE)
    if (nrow(w)) pares_q3 <- data.frame(item1 = colnames(q3)[w[, 1]], item2 = colnames(q3)[w[, 2]],
                                        q3 = q3[w], q3_relativo = rel[w], row.names = NULL)
  }
  if (!is.null(pares_q3)) avisos <- c(avisos, paste0(nrow(pares_q3), " par(es) de itens com Q3 mais de ", .pr_num(q3_limite, 2, FALSE, decimal),
                                                     " acima da média dos Q3: possível dependência local (conteúdo repetido, item derivado de outro). O valor crítico varia com n e número de itens (Christensen et al., 2017)."))
  # Confiabilidade e informação
  rxx_m <- tryCatch(.tri_rxx_marginal(fit), error = function(e) NA_real_)
  var_lat <- tryCatch(.tri_var_latente(fit)[["var"]], error = function(e) 1)
  fs <- tryCatch(suppressWarnings(mirt::fscores(fit, method = "EAP", full.scores.SE = TRUE, verbose = FALSE)), error = function(e) NULL)
  rxx_e <- if (is.null(fs)) NA_real_ else tryCatch(unname(mirt::empirical_rxx(fs))[1], error = function(e) NA_real_)
  rm(fs)
  # Precisão condicional = σ²/(σ² + EP²(θ)): ,80 com EP ≤ σ/2; ,90 com EP ≤ σ/3.
  grade <- seq(theta[1], theta[2], length.out = 161) * sqrt(var_lat)
  info <- tryCatch(mirt::testinfo(fit, matrix(grade)), error = function(e) rep(NA_real_, length(grade)))
  ep <- 1 / sqrt(info)
  prec <- var_lat / (var_lat + ep^2)
  faixa80 <- if (any(prec >= .80, na.rm = TRUE)) range(grade[prec >= .80]) else c(NA_real_, NA_real_)
  faixa90 <- if (any(prec >= .90, na.rm = TRUE)) range(grade[prec >= .90]) else c(NA_real_, NA_real_)
  informacao <- data.frame(var_latente = var_lat, theta_pico = grade[which.max(info)], info_pico = max(info, na.rm = TRUE),
                           precisao80_de = faixa80[1], precisao80_ate = faixa80[2],
                           precisao90_de = faixa90[1], precisao90_ate = faixa90[2])
  # Dimensionalidade e ajuste de pessoas na família Rasch
  autovalor <- NA_real_; pessoas <- NULL
  if (.tri_familia_rasch(modelo)) {
    autovalor <- .tri_pca_residuos(fit, d)
    pf <- tryCatch(suppressWarnings(mirt::personfit(fit)), error = function(e) NULL)
    if (!is.null(pf)) pessoas <- data.frame(n = sum(is.finite(pf$Zh)),
                                            pct_zh_abaixo_m196 = round(100 * mean(pf$Zh < -1.96, na.rm = TRUE), 1),
                                            pct_infit_acima_15 = round(100 * mean(pf$infit > 1.5, na.rm = TRUE), 1))
    if (is.finite(autovalor) && autovalor >= 2)
      avisos <- c(avisos, paste0("Primeiro contraste dos resíduos com autovalor ", .pr_num(autovalor, 2, TRUE, decimal),
                                 " (régua de Linacre usada no curso: abaixo de 2): sinal de dimensão secundária. Conferir com AFC de um fator ou EGA."))
  }
  # Tabelas
  tabela <- if (is.null(glob)) NULL else glob
  ef <- data.frame(indice = c("confiabilidade marginal", "confiabilidade empírica (EAP)"),
                   valor = c(rxx_m, rxx_e), ic_inf = NA_real_, ic_sup = NA_real_, principal = c(TRUE, FALSE))
  if (.tri_familia_rasch(modelo)) ef <- rbind(ef, data.frame(indice = "autovalor do 1º contraste (resíduos)", valor = autovalor,
                                                             ic_inf = NA_real_, ic_sup = NA_real_, principal = FALSE))
  ef <- rbind(ef, data.frame(indice = c("maior Q3 relativo", "média dos Q3"), valor = c(q3_max, q3_media),
                             ic_inf = NA_real_, ic_sup = NA_real_, principal = FALSE))
  n_desaj <- if (!is.null(ajuste_itens) && "p_ajust" %in% names(ajuste_itens)) sum(ajuste_itens$p_ajust < .05, na.rm = TRUE) else NA
  n_mnsq <- if (!is.null(mnsq)) sum(mnsq$fora_regua, na.rm = TRUE) else NA
  if (.tri_familia_rasch(modelo) && !is.null(mnsq))
    decisoes <- c(decisoes, paste0("Infit e outfit (média dos quadrados) lidos pela régua ",
                                   if (regua_mnsq == "linacre") "de Linacre (0,5 a 1,5, faixa produtiva para a medida)" else
                                     paste0("do curso por tamanho de amostra (", .pr_num(lim[1], 1, TRUE, ","), " a ", .pr_num(lim[2], 1, TRUE, ","), " com n = ", n, ")"),
                                   ". O z padronizado acusa desvio irrelevante com n grande; não decide sozinho."))
  # Frase APA
  rotulo_ajuste <- if (!is.null(glob)) switch(glob$estatistica, M2 = "M2 (Maydeu-Olivares & Joe, 2006)",
                                              "M2*" = "M2* (Cai & Hansen, 2013)", C2 = "C2 (Cai & Monroe, 2014)") else NULL
  txt_glob <- if (!is.null(glob)) paste0(glob$estatistica, "(", .ee_gl(glob$gl, decimal), ") = ", .pr_num(glob$valor, 2, TRUE, decimal), sep,
                                         .pr_p(glob$p, decimal), sep, "RMSEA = ", .pr_num(glob$rmsea, 3, FALSE, decimal), sep,
                                         .tri_t(idioma, "IC 90% [", "90% CI ["), .pr_num(glob$rmsea_li, 3, FALSE, decimal), sep,
                                         .pr_num(glob$rmsea_ls, 3, FALSE, decimal), "]", sep, "SRMSR = ", .pr_num(glob$srmsr, 3, FALSE, decimal),
                                         sep, "CFI = ", .pr_num(glob$cfi, 3, FALSE, decimal), sep, "TLI = ", .pr_num(glob$tli, 3, FALSE, decimal)) else NULL
  amp_a <- if ("a" %in% names(par) && any(is.finite(par$a))) range(par$a, na.rm = TRUE) else NULL
  if (en) {
    apa <- paste0("A unidimensional ", .tri_nome(modelo, "en"), " was fitted to the ", length(itens), " items (N = ", n,
                  ") by marginal maximum likelihood (mirt; Chalmers, 2012).",
                  if (!is.null(txt_glob)) paste0(" Global fit was assessed with ", rotulo_ajuste, ": ", txt_glob, ".") else " Global fit was not computed.",
                  if (!is.null(amp_a)) paste0(" Discriminations ranged from ", .pr_num(amp_a[1], 2, TRUE, decimal), " to ", .pr_num(amp_a[2], 2, TRUE, decimal), ".") else "",
                  if (is.finite(n_desaj)) paste0(" The S-X2 statistic (Orlando & Thissen, 2000), with ", corrigir, "-adjusted p values, flagged ", n_desaj, " item", if (n_desaj != 1) "s" else "", ".") else "",
                  if (is.finite(n_mnsq)) paste0(" Infit or outfit mean squares fell outside ", .pr_num(lim[1], 1, TRUE, decimal), "-", .pr_num(lim[2], 1, TRUE, decimal), " for ", n_mnsq, " item", if (n_mnsq != 1) "s" else "", ".") else "",
                  " Marginal reliability was ", .pr_num(rxx_m, 2, FALSE, decimal), " and empirical reliability ", .pr_num(rxx_e, 2, FALSE, decimal), ".",
                  if (is.finite(faixa80[1])) paste0(" Precision reached .80 or more from θ = ", .pr_num(faixa80[1], 1, TRUE, decimal), " to θ = ", .pr_num(faixa80[2], 1, TRUE, decimal), ".") else " Precision did not reach .80 anywhere on the θ grid.")
  } else {
    apa <- paste0("Ajustou-se o ", .tri_nome(modelo), ", unidimensional, aos ", length(itens), " itens (N = ", n,
                  "), por máxima verossimilhança marginal (mirt; Chalmers, 2012).",
                  if (!is.null(txt_glob)) paste0(" O ajuste global foi avaliado pelo ", rotulo_ajuste, ": ", txt_glob, ".") else " O ajuste global não foi calculado.",
                  if (!is.null(amp_a)) paste0(" As discriminações variaram de ", .pr_num(amp_a[1], 2, TRUE, decimal), " a ", .pr_num(amp_a[2], 2, TRUE, decimal), ".") else "",
                  if (is.finite(n_desaj)) paste0(" O S-X² (Orlando & Thissen, 2000), com p corrigido (", corrigir, "), indicou desajuste em ", n_desaj, " ite", if (n_desaj == 1) "m" else "ns", ".") else "",
                  if (is.finite(n_mnsq)) paste0(" Infit ou outfit ficaram fora de ", .pr_num(lim[1], 1, TRUE, decimal), " a ", .pr_num(lim[2], 1, TRUE, decimal), " em ", n_mnsq, " ite", if (n_mnsq == 1) "m" else "ns", ".") else "",
                  " A confiabilidade marginal foi ", .pr_num(rxx_m, 2, FALSE, decimal), " e a empírica, ", .pr_num(rxx_e, 2, FALSE, decimal), ".",
                  if (is.finite(faixa80[1])) paste0(" A precisão ficou em ", .pr_num(.80, 2, FALSE, decimal), " ou mais de θ = ", .pr_num(faixa80[1], 1, TRUE, decimal), " a θ = ", .pr_num(faixa80[2], 1, TRUE, decimal), ".") else
                    paste0(" A precisão não chegou a ", .pr_num(.80, 2, FALSE, decimal), " em nenhum ponto da grade de θ."))
  }
  .ee_resultado(paste0("TRI: ", .tri_nome(modelo)), n = n, n_excluidos = pp$n_excluidos, tabela = tabela, efeito = ef,
                decisoes = decisoes, avisos = unique(avisos), apa = apa,
                extra = list(parametros = par, ajuste_itens = ajuste_itens, pares_q3 = pares_q3, informacao = informacao,
                             pessoas = pessoas, categorias = data.frame(item = itens, categorias = unname(k)),
                             modelo_mirt = fit, titulo_posthoc = NULL,
                             tabelas_extra = list("Parâmetros dos itens" = par, "Ajuste dos itens" = ajuste_itens,
                                                  "Pares com dependência local" = pares_q3, "Informação do teste" = informacao,
                                                  "Ajuste de pessoas (agregado)" = pessoas)))
}

# Comparação de modelos: log-verossimilhança, AIC, BIC e BIC ajustado; razão de
# verossimilhança só entre modelos aninhados (Rasch ⊂ 2PL ⊂ 3PL; RSM ⊂ PCM ⊂ GPCM). GRM e
# GPCM não são aninhados: só os critérios de informação.
tri_comparar <- function(dados, itens, modelos = c("pcm", "gpcm", "graded"), semente = 2026,
                         decimal = ".", idioma = "pt") {
  .pr_exige("mirt")
  modelos <- vapply(modelos, function(m) match.arg(m, .tri_modelos), "")
  pp <- .tri_prepara(dados, itens, "listwise")
  d <- pp$d; n <- nrow(d)
  ajustes <- lapply(modelos, function(m) .tri_mirt(d, .tri_itemtype(m), semente))
  names(ajustes) <- modelos
  falhou <- vapply(ajustes, function(a) is.null(a$fit), logical(1))
  if (all(falhou)) stop("nenhum modelo ajustou: ", paste(vapply(ajustes, function(a) if (is.null(a$erro)) "" else a$erro, ""), collapse = "; "), call. = FALSE)
  ok <- modelos[!falhou]
  tab <- do.call(rbind, lapply(ok, function(m) {
    f <- ajustes[[m]]$fit
    ll <- mirt::extract.mirt(f, "logLik"); np <- mirt::extract.mirt(f, "nest")
    data.frame(modelo = m, logLik = ll, parametros = np, AIC = -2 * ll + 2 * np, BIC = -2 * ll + log(n) * np,
               SABIC = -2 * ll + log((n + 2) / 24) * np, convergiu = isTRUE(mirt::extract.mirt(f, "converged")), row.names = NULL)
  }))
  aninhados <- list(c("rasch", "2PL"), c("2PL", "3PL"), c("rasch", "3PL"), c("rsm", "pcm"), c("pcm", "gpcm"), c("rsm", "gpcm"))
  lrt <- do.call(rbind, lapply(aninhados, function(par) {
    if (!all(par %in% ok)) return(NULL)
    a <- tab[tab$modelo == par[1], ]; b <- tab[tab$modelo == par[2], ]
    qui <- 2 * (b$logLik - a$logLik); gl <- b$parametros - a$parametros
    data.frame(restrito = par[1], livre = par[2], qui2 = qui, gl = gl,
               p = if (gl > 0) stats::pchisq(qui, gl, lower.tail = FALSE) else NA_real_, row.names = NULL)
  }))
  melhor_bic <- tab$modelo[which.min(tab$BIC)]; melhor_aic <- tab$modelo[which.min(tab$AIC)]
  avisos <- character()
  if (any(falhou)) avisos <- c(avisos, paste0("Não ajustaram: ", paste(modelos[falhou], collapse = ", ")))
  if (any(!tab$convergiu)) avisos <- c(avisos, paste0("Sem convergência: ", paste(tab$modelo[!tab$convergiu], collapse = ", "), ". Não comparar antes de resolver."))
  if (melhor_aic != melhor_bic) avisos <- c(avisos, paste0("AIC e BIC discordam (AIC: ", melhor_aic, "; BIC: ", melhor_bic, "). O BIC pesa mais a parcimônia; decidir também pelo uso pretendido do escore."))
  sep <- .ee_sep(decimal)
  apa <- if (identical(idioma, "en")) paste0("Competing IRT models were compared by information criteria", if (!is.null(lrt)) " and likelihood-ratio tests for nested pairs" else "",
                                             "; the lowest BIC was obtained by the ", .tri_nome(melhor_bic, "en"), ".") else
    paste0("Os modelos concorrentes foram comparados pelos critérios de informação", if (!is.null(lrt)) " e, entre os aninhados, pela razão de verossimilhança" else "",
           "; o menor BIC foi o do ", .tri_nome(melhor_bic), ".")
  .ee_resultado("TRI: comparação de modelos", n = n, n_excluidos = pp$n_excluidos, tabela = tab, posthoc = lrt,
                decisoes = c("Mesma amostra (casos completos), mesma semente e mesmo estimador em todos os modelos.",
                             "GRM e GPCM não são aninhados: compará-los só pelos critérios de informação."),
                avisos = avisos, apa = apa, extra = list(titulo_posthoc = "Razão de verossimilhança (aninhados):",
                                                         modelos_mirt = lapply(ajustes[ok], `[[`, "fit")))
}

# DIF. Duas rotas:
#  - metodo = "mirt": modelo de grupos múltiplos. Sem âncoras, parte de todos os itens
#    iguais entre os grupos (média e variância latentes livres no focal) e libera um item
#    por vez em sequência ("drop_sequential"), como o script do curso; com âncoras
#    escolhidas pelo usuário, testa os demais adicionando as restrições ("add"). Para a
#    família Rasch sai também o contraste de localização (régua do Winsteps usada no
#    curso: < ,43 desprezível, ,43 a ,64 moderado, > ,64 grande, em logits).
#  - metodo = "logistica": regressão logística (binária para itens dicotômicos, ordinal
#    para os demais) com o escore total como critério de pareamento, purificado por
#    padrão (Swaminathan & Rogers, 1990; Zumbo, 1999). DIF uniforme (grupo), não
#    uniforme (grupo × escore) e o teste conjunto de 2 gl; tamanho pelo ΔR² de
#    Nagelkerke, classificado pela régua de Jodoin e Gierl (2001): < ,035 desprezível (A),
#    ,035 a ,070 moderado (B), > ,070 grande (C).
#  - metodo = "mh" (v1.6): Mantel-Haenszel para itens dicotômicos (Martinková & Hladká,
#    2023, 9.3.2), estratos pelo escore total (com o item), purificado por padrão; αMH
#    (chance de acerto do grupo de referência, o primeiro nível, sobre a do focal no mesmo
#    escore), χ² de MH com correção de continuidade e ΔMH = -2,35 ln αMH na régua do ETS
#    como no difR: |Δ| < 1 A (desprezível), 1 a 1,5 B (moderado), >= 1,5 C (grande).
tri_dif <- function(dados, itens, grupo, modelo = "graded", metodo = c("mirt", "logistica", "mh"),
                    ancoras = NULL, purificar = TRUE, corrigir = "BH", alfa = .05, semente = 2026,
                    decimal = ".", idioma = "pt") {
  metodo <- match.arg(metodo)
  modelo <- match.arg(modelo, .tri_modelos)
  if (!grupo %in% names(dados)) stop("variável de grupo ausente no banco: ", grupo, call. = FALSE)
  base <- as.data.frame(dados)
  gch <- as.character(base[[grupo]])
  base <- base[!is.na(gch) & nzchar(trimws(gch)), , drop = FALSE]
  pp <- .tri_prepara(base, itens, "listwise")
  d <- pp$d
  g <- droplevels(as.factor(base[rownames(d), grupo]))
  if (nlevels(g) != 2) stop("o DIF daqui compara 2 grupos; o grupo tem ", nlevels(g), " níveis com dados.", call. = FALSE)
  n_g <- table(g)
  avisos <- character(); decisoes <- character()
  if (any(n_g < 200)) avisos <- c(avisos, paste0("Grupo com menos de 200 casos (", paste(names(n_g), n_g, sep = " = ", collapse = "; "),
                                                 "): o teste de DIF tem pouco poder; ausência de DIF detectado não mostra invariância."))
  k <- .tri_categorias(d); dic <- all(k == 2)
  if (metodo == "mirt") {
    .pr_exige("mirt")
    if (modelo %in% c("rasch", "2PL", "3PL") && !dic) stop("modelo ", modelo, " é para itens dicotômicos.", call. = FALSE)
    if (modelo == "rsm") { modelo <- "pcm"; decisoes <- c(decisoes, "DIF com o modelo de crédito parcial: na escala de avaliação os limiares são comuns aos itens e não se liberam item a item.") }
    it <- .tri_itemtype(if (modelo == "pcm" && dic) "rasch" else modelo)
    if (!is.null(ancoras) && length(setdiff(ancoras, itens))) stop("âncoras fora da lista de itens: ", paste(setdiff(ancoras, itens), collapse = ", "), call. = FALSE)
    roda_dif <- function(fit, esquema, testar) {
      v <- mirt::mod2values(fit)
      wp <- unique(v$name[v$group == levels(g)[1] & v$item == testar[1] & v$est])
      set.seed(semente)
      r <- tryCatch(suppressWarnings(suppressMessages(mirt::DIF(fit, which.par = wp, scheme = esquema, items2test = testar,
                                                               p.adjust = "none", verbose = FALSE))), error = function(e) e)
      if (inherits(r, "error")) stop("mirt::DIF falhou: ", conditionMessage(r), call. = FALSE)
      data.frame(item = rownames(r), qui2 = r$X2, gl = r$df, p = r$p, row.names = NULL)
    }
    if (is.null(ancoras)) {
      # Etapa 1: todos os itens iguais entre os grupos; cada item liberado por vez (drop).
      b1 <- .tri_mirt(d, it, semente, grupo = g, invariancia = c(itens, "free_means", "free_var"))
      if (is.null(b1$fit)) stop("o mirt não ajustou o modelo de grupos múltiplos: ", b1$erro, call. = FALSE)
      t1 <- roda_dif(b1$fit, "drop", itens); t1$etapa <- 1L
      cand <- t1$item[stats::p.adjust(t1$p, corrigir) < alfa]
      tab <- t1
      # Etapa 2: os itens sem DIF na etapa 1 viram âncoras; os candidatos são testados de novo (add).
      if (length(cand) && length(cand) <= length(itens) - 2) {
        b2 <- .tri_mirt(d, it, semente, grupo = g, invariancia = c(setdiff(itens, cand), "free_means", "free_var"))
        if (!is.null(b2$fit)) {
          t2 <- roda_dif(b2$fit, "add", cand); t2$etapa <- 2L
          tab <- rbind(t1[!t1$item %in% cand, ], t2)
          tab <- tab[match(itens, tab$item), ]; rownames(tab) <- NULL
        }
      }
      decisoes <- c(decisoes, "Sem âncoras declaradas, em duas etapas: (1) todos os itens iguais entre os grupos e cada um liberado por vez; (2) os itens sem DIF na etapa 1 viram âncoras e os candidatos são testados de novo. O curso usa a busca sequencial do mirt (drop_sequential), que segue a mesma lógica; a skill fixa as duas etapas para devolver todos os itens na tabela.")
    } else {
      b2 <- .tri_mirt(d, it, semente, grupo = g, invariancia = c(ancoras, "free_means", "free_var"))
      if (is.null(b2$fit)) stop("o mirt não ajustou o modelo de grupos múltiplos: ", b2$erro, call. = FALSE)
      tab <- roda_dif(b2$fit, "add", setdiff(itens, ancoras)); tab$etapa <- 2L
      decisoes <- c(decisoes, paste0("Âncoras declaradas (", paste(ancoras, collapse = ", "), "): média e variância latentes livres no grupo focal; cada item restante testado ao igualar seus parâmetros entre os grupos (add)."))
    }
    tab$p_ajust <- stats::p.adjust(tab$p, corrigir)
    tab$dif <- tab$p_ajust < alfa
    if (.tri_familia_rasch(modelo)) {
      livre <- .tri_mirt(d, it, semente, grupo = g, invariancia = c(setdiff(itens, tab$item[tab$dif]), "free_means", "free_var"))
      if (!is.null(livre$fit)) {
        cf <- mirt::coef(livre$fit, IRTpars = modelo != "rsm", simplify = TRUE)
        loc <- function(m) { bs <- grep(if (modelo == "rsm") "^c$" else "^b[0-9]*$", colnames(m), value = TRUE); rowMeans(m[, bs, drop = FALSE]) }
        dl <- loc(cf[[2]]$items) - loc(cf[[1]]$items)
        tab$contraste_logits <- ifelse(tab$dif, unname(dl[tab$item]), NA_real_)
        tab$classe <- ifelse(!tab$dif, "âncora (sem DIF)", ifelse(abs(tab$contraste_logits) < .43, "desprezível", ifelse(abs(tab$contraste_logits) <= .64, "moderado", "grande")))
        decisoes <- c(decisoes, "Contraste = localização média do item no grupo 2 menos no grupo 1, em logits, com os itens sem DIF como âncora. Régua do Winsteps usada no curso (Linacre): < ,43 desprezível; ,43 a ,64 moderado; > ,64 grande.")
      }
    }
    metodo_txt <- "modelo de grupos múltiplos no mirt (razão de verossimilhança)"
  } else if (metodo == "mh") {
    if (!dic) stop("Mantel-Haenszel aqui é para itens dicotômicos; com itens ordinais, metodo = \"logistica\" ou \"mirt\".", call. = FALSE)
    d01 <- as.data.frame(lapply(d, .tri_recodifica))
    ref <- g == levels(g)[1]
    mh_item <- function(j, s) {
      y <- d01[[j]]; A <- B <- C <- D <- numeric()
      num <- den <- soma_a <- soma_e <- soma_v <- 0
      for (k in unique(s)) {
        idx <- s == k; nk <- sum(idx); if (nk <= 1) next
        a <- sum(y[idx & ref] == 1); b <- sum(y[idx & ref] == 0); c <- sum(y[idx & !ref] == 1); dd <- sum(y[idx & !ref] == 0)
        nr <- a + b; nf <- c + dd; m1 <- a + c; m0 <- b + dd
        if (nr == 0 || nf == 0) next
        num <- num + a * dd / nk; den <- den + b * c / nk
        soma_a <- soma_a + a; soma_e <- soma_e + nr * m1 / nk; soma_v <- soma_v + nr * nf * m1 * m0 / (nk^2 * (nk - 1))
      }
      alfa_mh <- num / den
      x2 <- (abs(soma_a - soma_e) - .5)^2 / soma_v
      data.frame(item = itens[j], qui2_mh = x2, p = stats::pchisq(x2, 1, lower.tail = FALSE), alfa_mh = alfa_mh,
                 delta_mh = -2.35 * log(alfa_mh), row.names = NULL)
    }
    roda <- function(livres) {
      s_base <- if (length(livres)) rowSums(d01[, setdiff(itens, livres), drop = FALSE]) else rowSums(d01)
      do.call(rbind, lapply(seq_along(itens), function(j) mh_item(j, if (itens[j] %in% livres) s_base + d01[[j]] else s_base)))
    }
    tab <- roda(character()); iter <- 1; marcados <- character()
    if (purificar) {
      repeat {
        novos <- tab$item[stats::p.adjust(tab$p, corrigir) < alfa]
        if (setequal(novos, marcados) || iter >= 10 || length(novos) >= length(itens) - 1) break
        marcados <- novos; tab <- roda(marcados); iter <- iter + 1
      }
      decisoes <- c(decisoes, paste0("Escore de pareamento purificado (", iter, " rodada", if (iter != 1) "s" else "", "; o item testado sempre entra no próprio escore, como no difR)."))
    }
    tab$p_ajust <- stats::p.adjust(tab$p, corrigir)
    tab$dif <- tab$p_ajust < alfa
    tab$classe <- ifelse(abs(tab$delta_mh) < 1, "A (desprezível)", ifelse(abs(tab$delta_mh) < 1.5, "B (moderado)", "C (grande)"))
    tab$dif_relevante <- tab$dif & abs(tab$delta_mh) >= 1
    tab$favorece <- ifelse(tab$alfa_mh > 1, levels(g)[1], levels(g)[2])
    decisoes <- c(decisoes, paste0("Mantel-Haenszel com estratos pelo escore total; grupo de referência = ", levels(g)[1], ", focal = ", levels(g)[2],
                                   "; αMH > 1 favorece a referência. χ² com correção de continuidade, p corrigido (", corrigir, ")."),
                  "Tamanho: ΔMH = -2,35 ln αMH (escala delta do ETS) pela régua do difR: |Δ| < 1 A, 1 a 1,5 B, >= 1,5 C. O MH só vê DIF uniforme; DIF não uniforme pede metodo = \"logistica\".")
    metodo_txt <- "Mantel-Haenszel"
  } else {
    .pr_exige("MASS")
    nn <- table(g)
    escore_total <- rowSums(d)
    yfac <- lapply(d, function(x) factor(.tri_recodifica(x), ordered = TRUE))
    ajusta <- function(j, mat) {
      y <- yfac[[j]]; dd <- data.frame(y = y, s = mat, g = g)
      bin <- nlevels(y) == 2
      f <- function(fo) if (bin) stats::glm(fo, data = dd, family = stats::binomial()) else MASS::polr(fo, data = dd, Hess = TRUE)
      m0 <- suppressWarnings(f(y ~ s)); m1 <- suppressWarnings(f(y ~ s + g)); m2 <- suppressWarnings(f(y ~ s * g))
      nulo <- suppressWarnings(f(y ~ 1))
      ll <- vapply(list(nulo, m0, m1, m2), function(m) as.numeric(stats::logLik(m)), 0)
      nag <- function(l) { r2cs <- 1 - exp(2 * (ll[1] - l) / length(y)); r2cs / (1 - exp(2 * ll[1] / length(y))) }
      mcf <- function(l) 1 - l / ll[1]
      data.frame(item = itens[j], qui2_uniforme = 2 * (ll[3] - ll[2]), p_uniforme = stats::pchisq(2 * (ll[3] - ll[2]), 1, lower.tail = FALSE),
                 qui2_nao_uniforme = 2 * (ll[4] - ll[3]), p_nao_uniforme = stats::pchisq(2 * (ll[4] - ll[3]), 1, lower.tail = FALSE),
                 qui2_conjunto = 2 * (ll[4] - ll[2]), p = stats::pchisq(2 * (ll[4] - ll[2]), 2, lower.tail = FALSE),
                 delta_r2_nagelkerke = nag(ll[4]) - nag(ll[2]), delta_r2_mcfadden = mcf(ll[4]) - mcf(ll[2]), row.names = NULL)
    }
    roda <- function(livres) {
      s_base <- if (length(livres)) rowSums(d[, setdiff(itens, livres), drop = FALSE]) else escore_total
      do.call(rbind, lapply(seq_along(itens), function(j) {
        # O item testado entra no escore de pareamento (como no difR), mesmo purificado.
        s <- if (itens[j] %in% livres) s_base + d[[j]] else s_base
        ajusta(j, s)
      }))
    }
    tab <- roda(character()); iter <- 1; marcados <- character()
    if (purificar) {
      repeat {
        novos <- tab$item[stats::p.adjust(tab$p, corrigir) < alfa]
        if (setequal(novos, marcados) || iter >= 10 || length(novos) >= length(itens) - 1) break
        marcados <- novos; tab <- roda(marcados); iter <- iter + 1
      }
      decisoes <- c(decisoes, paste0("Escore de pareamento purificado: itens com DIF saem do escore e o teste se repete até estabilizar (", iter,
                                     " rodada", if (iter != 1) "s" else "", "; o item testado sempre entra no próprio escore, como no difR)."))
    }
    tab$p_ajust <- stats::p.adjust(tab$p, corrigir)
    tab$dif <- tab$p_ajust < alfa
    tab$classe <- ifelse(tab$delta_r2_nagelkerke < .035, "A (desprezível)", ifelse(tab$delta_r2_nagelkerke <= .070, "B (moderado)", "C (grande)"))
    tab$dif_relevante <- tab$dif & tab$delta_r2_nagelkerke >= .035
    decisoes <- c(decisoes, paste0("Regressão ", if (dic) "logística" else "logística ordinal (MASS::polr)", " com o escore total como pareamento; teste conjunto de 2 gl (uniforme + não uniforme), p corrigido (", corrigir, ")."),
                  paste0("Tamanho: ΔR² de Nagelkerke pela régua de Jodoin e Gierl (2001), derivada para itens dicotômicos",
                         if (!dic) paste0("; com itens ordinais, Choi et al. (2011) usam ΔR² de McFadden ≥ ", .pr_num(.02, 2, FALSE, decimal), " (coluna ao lado)") else "",
                         ". DIF relevante = significativo e ΔR² ≥ ", .pr_num(.035, 3, FALSE, decimal), "."))
    metodo_txt <- if (dic) "regressão logística" else "regressão logística ordinal"
  }
  n_dif <- sum(tab$dif, na.rm = TRUE)
  rel <- if ("dif_relevante" %in% names(tab)) sum(tab$dif_relevante, na.rm = TRUE) else NA
  if (n_dif > 0 && n_dif >= length(itens) / 2) avisos <- c(avisos, "Metade ou mais dos itens com DIF: suspeitar de diferença no traço mal modelada ou de âncoras contaminadas antes de concluir que o teste é enviesado.")
  if (n_dif == 0) avisos <- c(avisos, "Nenhum item com DIF depois da correção: isso não mostra que os itens funcionam igual (o teste pode ter tido pouco poder). Relatar o tamanho dos grupos e, se possível, as curvas de escore esperado por grupo.")
  sep <- .ee_sep(decimal)
  apa <- if (identical(idioma, "en")) paste0("Differential item functioning by ", grupo, " (", paste(names(n_g), n_g, sep = ": n = ", collapse = "; "), ") was tested with ",
                                             if (metodo == "mirt") "multiple-group IRT models in mirt (likelihood-ratio tests)" else if (metodo == "mh") "the Mantel-Haenszel procedure with the purified total score as matching variable" else "logistic regression with the purified total score as matching variable",
                                             ", with ", corrigir, "-adjusted p values; ", n_dif, " of ", length(itens), " items showed DIF",
                                             if (is.finite(rel)) paste0(" (", rel, " with non-negligible effect, ", if (metodo == "mh") "|ΔMH| ≥ 1" else "ΔR² ≥ .035", ")") else "", ".") else
    paste0("O funcionamento diferencial dos itens por ", grupo, " (", paste(names(n_g), n_g, sep = ": n = ", collapse = sep), ") foi testado por ", metodo_txt,
           ", com p corrigido (", corrigir, "); ", n_dif, " dos ", length(itens), " itens mostraram DIF",
           if (is.finite(rel)) paste0(" (", rel, " com efeito não desprezível, ", if (metodo == "mh") paste0("|ΔMH| ≥ 1") else paste0("ΔR² ≥ ", .pr_num(.035, 3, FALSE, decimal)), ")") else "", ".")
  .ee_resultado(paste0("DIF: ", metodo_txt), n = nrow(d), n_excluidos = nrow(dados) - nrow(d),
                tabela = tab, decisoes = decisoes, avisos = unique(avisos), apa = apa,
                extra = list(grupos = as.data.frame(n_g, responseName = "n")))
}

# ---------------------------------------------------------------- TRI multidimensional
# tri_multidimensional(): GRM, GPCM, PCM (Rasch) ou 2PL com mais de uma dimensão, como nas
# aulas "GRM multidimensional (exploratório/confirmatório)" e "PCM multidimensional".
#   dimensoes = k (inteiro): exploratório com k dimensões e rotação (oblimin, como o padrão
#               do mirt e das aulas)
#   dimensoes = list(DEP = c(...), ANS = c(...)): confirmatório, dimensões correlacionadas
# Método: EM com 1 ou 2 dimensões; MHRM com 3 ou mais (as aulas usam MHRM a partir de 3
# dimensões; a quadratura do EM cresce exponencialmente). Confiabilidade empírica por
# dimensão a partir dos θ EAP (calculados aqui e descartados: nada por pessoa sai).
tri_multidimensional <- function(dados, itens, dimensoes, modelo = c("graded", "gpcm", "pcm", "rasch", "2PL"), rotacao = "oblimin",
                                 metodo = NULL, comparar_unidimensional = TRUE, corrigir = "BH", semente = 2026,
                                 decimal = ".", idioma = "pt") {
  .pr_exige("mirt")
  modelo <- match.arg(modelo)
  pp <- .tri_prepara(dados, itens, "listwise")
  d <- as.data.frame(lapply(pp$d, .tri_recodifica)); n <- nrow(d)
  confirmatorio <- is.list(dimensoes)
  if (confirmatorio) {
    if (is.null(names(dimensoes)) || any(!nzchar(names(dimensoes)))) stop("dimensoes: lista nomeada, por exemplo list(DEP = c('i3', 'i5'), ...)")
    fora <- setdiff(unlist(dimensoes), itens); if (length(fora)) stop("itens das dimensões fora de 'itens': ", paste(fora, collapse = ", "))
    k <- length(dimensoes)
    linhas <- vapply(names(dimensoes), function(nm) paste0(nm, " = ", paste(match(dimensoes[[nm]], itens), collapse = ", ")), "")
    sintaxe <- paste(c(linhas, if (k > 1) paste0("COV = ", paste(names(dimensoes), collapse = "*"))), collapse = "\n")
    mod <- mirt::mirt.model(sintaxe, itemnames = itens)
  } else {
    k <- as.integer(dimensoes); if (k < 2) stop("dimensoes: inteiro >= 2 (exploratório) ou lista de itens por dimensão (confirmatório)")
    if (.tri_familia_rasch(modelo)) stop("exploratório pede discriminações livres: usar graded, gpcm ou 2PL; Rasch multidimensional só confirmatório")
    mod <- k
  }
  if (is.null(metodo)) metodo <- if (k >= 3) "MHRM" else "EM"
  it <- .tri_itemtype(modelo)
  set.seed(semente)
  av <- character()
  ajusta <- function(m, met) withCallingHandlers(mirt::mirt(d, m, itemtype = it, method = met, SE = TRUE, verbose = FALSE),
                                                  warning = function(w) { av <<- c(av, conditionMessage(w)); invokeRestart("muffleWarning") })
  fit <- ajusta(mod, metodo)
  if (!mirt::extract.mirt(fit, "converged")) av <- c(av, "o mirt não declarou convergência")
  sm <- if (confirmatorio) mirt::summary(fit, verbose = FALSE) else mirt::summary(fit, rotate = rotacao, verbose = FALSE)
  cargas <- data.frame(item = itens, round(unclass(sm$rotF), 6), h2 = sm$h2, row.names = NULL, check.names = FALSE)
  fcor <- as.data.frame(sm$fcor); fcor <- cbind(dimensao = rownames(sm$fcor), fcor); rownames(fcor) <- NULL
  cf <- mirt::coef(fit, simplify = TRUE)$items
  par <- data.frame(item = rownames(cf), cf[, colSums(!is.na(cf)) > 0 & !colnames(cf) %in% c("g", "u"), drop = FALSE], row.names = NULL, check.names = FALSE)
  ag <- .tri_ajuste_global(fit, all(.tri_categorias(d) == 2))
  if (.tri_familia_rasch(modelo)) {
    inf <- tryCatch(suppressWarnings(mirt::itemfit(fit, fit_stats = "infit")), error = function(e) NULL)
    ajuste_itens <- if (!is.null(inf)) data.frame(item = inf$item, infit = inf$infit, outfit = inf$outfit, row.names = NULL)
  } else {
    sx2 <- tryCatch(suppressWarnings(mirt::itemfit(fit, fit_stats = "S_X2")), error = function(e) NULL)
    ajuste_itens <- if (!is.null(sx2)) data.frame(item = sx2$item, s_x2 = sx2$S_X2, gl = sx2$df.S_X2, rmsea = sx2$RMSEA.S_X2,
                                                  p = sx2$p.S_X2, p_ajust = stats::p.adjust(sx2$p.S_X2, corrigir), row.names = NULL)
  }
  th <- tryCatch(suppressWarnings(mirt::fscores(fit, full.scores.SE = TRUE, rotate = if (confirmatorio) "none" else rotacao, verbose = FALSE)),
                 error = function(e) NULL)
  rxx <- if (!is.null(th)) { r <- mirt::empirical_rxx(th); data.frame(dimensao = names(r), confiabilidade_empirica = unname(r), stringsAsFactors = FALSE) }
  rm(th)
  comp <- NULL
  if (comparar_unidimensional) {
    f1 <- tryCatch(suppressWarnings(mirt::mirt(d, 1, itemtype = it, verbose = FALSE)), error = function(e) NULL)
    if (!is.null(f1)) {
      a <- suppressWarnings(mirt::anova(f1, fit))
      comp <- data.frame(modelo = c("unidimensional", paste0(k, " dimensões")), aic = a$AIC, sabic = a$SABIC, bic = a$BIC, ll = a$logLik,
                         qui2 = c(NA, a$X2[2]), gl = c(NA, a$df[2]), p = c(NA, a$p[2]), row.names = NULL)
    }
  }
  efeito <- if (!is.null(rxx)) data.frame(indice = paste0("confiabilidade empírica ", rxx$dimensao), valor = rxx$confiabilidade_empirica,
                                          ic_inf = NA_real_, ic_sup = NA_real_, principal = TRUE, stringsAsFactors = FALSE)
  en <- identical(idioma, "en")
  aj_txt <- if (!is.null(ag)) paste0(ag$estatistica, "(", round(ag$gl), ") = ", .pr_num(ag$valor, 2, TRUE, decimal), ", ", .pr_p(ag$p, decimal),
                                     ", RMSEA = ", .pr_num(ag$rmsea, 3, FALSE, decimal), ", CFI = ", .pr_num(ag$cfi, 3, FALSE, decimal)) else
    .tri_t(idioma, "ajuste global não calculado", "global fit not computed")
  rx_txt <- if (!is.null(rxx)) paste(paste0(rxx$dimensao, " = ", vapply(rxx$confiabilidade_empirica, .pr_num, "", 2, FALSE, decimal)), collapse = .ee_sep(decimal)) else ""
  apa <- if (en) paste0("A ", if (confirmatorio) "confirmatory" else "exploratory", " ", k, "-dimensional ", .tri_nome(modelo, "en"), " was fitted (mirt, ", metodo,
                        if (!confirmatorio) paste0(", ", rotacao, " rotation") else "", "; n = ", n, "): ", aj_txt, ". Empirical reliability: ", rx_txt, ".") else
    paste0("Foi ajustado o ", .tri_nome(modelo), " ", if (confirmatorio) "confirmatório" else "exploratório", " com ", k, " dimensões (mirt, ", metodo,
           if (!confirmatorio) paste0(", rotação ", rotacao) else "", "; n = ", n, "): ", aj_txt, ". Confiabilidade empírica: ", rx_txt, ".")
  extras <- list(`Cargas (rotacionadas no exploratório)` = cargas, `Correlações entre dimensões` = fcor, `Parâmetros (inclinação-intercepto)` = par)
  if (!is.null(ajuste_itens)) extras$`Ajuste dos itens` <- ajuste_itens
  if (!is.null(comp)) extras$`Comparação com o unidimensional` <- comp
  .ee_resultado(paste0("TRI multidimensional (", modelo, ", ", k, " dimensões, ", if (confirmatorio) "confirmatória" else "exploratória", ")"),
                n = n, n_excluidos = pp$n_excluidos, tabela = ag, efeito = efeito,
                decisoes = c(paste0("Estimação ", metodo, if (metodo == "MHRM") " (3 ou mais dimensões, como nas aulas; semente fixa porque o MHRM é estocástico)" else "", "."),
                             if (!confirmatorio) paste0("Rotação ", rotacao, "; cargas e θ na solução rotacionada.") else "Confirmatório com dimensões correlacionadas (COV).",
                             "Parâmetros na parametrização inclinação-intercepto (a1..ak, d): o mirt não converte para a, b quando há mais de uma dimensão.",
                             paste0("Ajuste dos itens com p corrigido (", corrigir, ")."),
                             "Confiabilidade empírica por dimensão (var(θ EAP) / [var(θ EAP) + média dos EP²]); θ calculado e descartado."),
                avisos = unique(av), apa = apa, extra = list(tabelas_extra = extras, modelo = fit))
}

# ---------------------------------------------------------------- Mokken (v1.7)
# escala_mokken(): TRI não paramétrica de Mokken (Meijer & Tendeiro, 2018, Wiley, cap. 15;
# pacote mokken, van der Ark, 2012). Coeficientes de escalabilidade H (escala), Hi (item) e
# Hij (pares) com EP; seleção automática de itens (AISP) com limite c; monotonicidade
# manifesta (check.monotonicity) e, se pedido, ordenação invariante dos itens (check.iio).
# Régua de Sijtsma e Molenaar (2002, p. 60), citada no capítulo: H < ,30 não escalável;
# ,30 a ,40 fraca; ,40 a ,50 média; >= ,50 forte; Hi abaixo de c (,30) sinaliza item.
# Ressalva do capítulo: itens que se ajustam ao Rasch com discriminação perto de 1 podem
# ter Hi < ,30; a AISP descarta itens que um modelo paramétrico aceitaria.
.tri_mokken_classe <- function(h) cut(h, c(-Inf, .3, .4, .5, Inf), right = FALSE, labels = c("não escalável", "fraca", "média", "forte"))

escala_mokken <- function(dados, itens, limite = .30, aisp = TRUE, monotonicidade = TRUE, iio = FALSE,
                          decimal = ".", idioma = "pt") {
  .pr_exige("mokken")
  pp <- .tri_prepara(as.data.frame(dados), itens, "listwise")
  X <- as.matrix(as.data.frame(lapply(pp$d, as.numeric)))
  if (any(abs(X - round(X)) > 1e-8)) stop("Mokken pede escores inteiros (categorias ordenadas)", call. = FALSE)
  X <- X - min(X)  # mesmo deslocamento para todos os itens: não muda H nem a ordem dos itens
  colnames(X) <- itens
  n <- nrow(X); k <- .tri_categorias(pp$d)
  avisos <- character(); decisoes <- character()
  if (length(unique(k)) > 1) avisos <- c(avisos, paste0("Itens com número de categorias observadas diferente (", paste(range(k), collapse = " a "),
    "): se for categoria não usada, o escore somado ainda vale; se for desenho, o escore somado não serve como medida ordinal (aviso do mokken)."))
  h <- suppressWarnings(mokken::coefH(X, se = TRUE, nice.output = FALSE, results = FALSE))
  Hi <- as.numeric(h$Hi); seHi <- as.numeric(h$se.Hi)
  Hij <- h$Hij; diag(Hij) <- NA
  tab <- data.frame(item = itens, Hi = Hi, ep_Hi = seHi, ic_inf = Hi - 1.96 * seHi, ic_sup = Hi + 1.96 * seHi,
                    abaixo_limite = Hi < limite, Hij_minimo = apply(Hij, 1, min, na.rm = TRUE), row.names = NULL)
  H <- as.numeric(h$H); seH <- as.numeric(h$se.H)
  if (aisp) {
    a <- suppressWarnings(mokken::aisp(X, lowerbound = limite, verbose = FALSE))
    tab$escala_aisp <- as.integer(a[, 1])
    decisoes <- c(decisoes, paste0("AISP com limite c = ", .pr_num(limite, 2, FALSE, decimal), ": 0 = item não selecionado; 1, 2... = escalas formadas, na ordem em que a AISP as monta."))
    if (any(tab$escala_aisp == 0)) avisos <- c(avisos, paste0("A AISP deixou fora ", sum(tab$escala_aisp == 0), " item(ns): ", paste(itens[tab$escala_aisp == 0], collapse = ", "), "."))
    if (length(unique(tab$escala_aisp[tab$escala_aisp > 0])) > 1) avisos <- c(avisos, "A AISP formou mais de uma escala: conferir contra a estrutura teórica e a AFC antes de aceitar a partição.")
  }
  mono <- NULL
  if (monotonicidade) {
    sm <- suppressWarnings(summary(mokken::check.monotonicity(X)))
    mono <- data.frame(item = itens, violacoes = sm[, "#vi"], violacoes_sig = sm[, "#zsig"], crit = sm[, "crit"], row.names = NULL)
    tab$crit_monotonicidade <- mono$crit
    decisoes <- c(decisoes, "Monotonicidade manifesta (Junker, 1993, pelo check.monotonicity): crit > 80 indica violação séria, 40 a 80 pede inspeção (convenção de Molenaar e Sijtsma, do manual do MSP).")
    if (any(mono$crit > 80)) avisos <- c(avisos, paste0("Violação séria de monotonicidade (crit > 80): ", paste(mono$item[mono$crit > 80], collapse = ", "), "."))
  }
  io <- NULL
  if (iio) {
    si <- suppressWarnings(summary(mokken::check.iio(X)))
    io <- data.frame(item = rownames(si$item.summary), crit_iio = si$item.summary[, "crit"], row.names = NULL)
    tab$crit_iio <- io$crit_iio[match(itens, io$item)]
    decisoes <- c(decisoes, paste0("Ordenação invariante dos itens (check.iio, método ", si$method, "): HT = ", .pr_num(si$HT, 2, FALSE, decimal),
                                   "; HT >= ,30 sustenta usar a ordem dos itens (Ligtvoet et al., pelo manual do mokken)."))
    if (is.finite(si$HT) && si$HT < .3) avisos <- c(avisos, "HT < ,30: a ordem de dificuldade dos itens muda ao longo do traço; não interpretar a ordem dos itens como igual para todos.")
  }
  if (any(Hij < 0, na.rm = TRUE)) avisos <- c(avisos, "Há pares de itens com Hij negativo: item invertido sem recodificar ou item fora da escala.")
  if (any(tab$abaixo_limite)) avisos <- c(avisos, paste0("Hi < ", .pr_num(limite, 2, FALSE, decimal), " em: ", paste(itens[tab$abaixo_limite], collapse = ", "),
                                                         ". Ressalva do capítulo 15 do Wiley: item que se ajusta ao Rasch com discriminação perto de 1 pode ficar abaixo de ,30."))
  if (n < 250) avisos <- c(avisos, paste0("n = ", n, ": os EP de H ficam largos; ler o IC, não só o valor."))
  classe <- as.character(.tri_mokken_classe(H))
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  ic_H <- c(H - 1.96 * seH, H + 1.96 * seH)
  apa <- if (en) paste0("Mokken scale analysis (mokken) of ", length(itens), " items (n = ", n, "): H = ", .pr_num(H, 2, FALSE, decimal), ", 95% CI [",
                        .pr_num(ic_H[1], 2, FALSE, decimal), ", ", .pr_num(ic_H[2], 2, FALSE, decimal), "], a ", c(`não escalável` = "non-scalable", fraca = "weak", `média` = "medium", forte = "strong")[[classe]],
                        " scale by Sijtsma and Molenaar's (2002) guidelines; item Hi ranged from ", .pr_num(min(Hi), 2, FALSE, decimal), " to ", .pr_num(max(Hi), 2, FALSE, decimal), ".") else
    paste0("A análise de escalas de Mokken (mokken) dos ", length(itens), " itens (n = ", n, ") deu H = ", .pr_num(H, 2, FALSE, decimal), sep,
           .ee_ic(ic_H[1], ic_H[2], 2, FALSE, decimal, idioma = idioma), ", escala ", classe, " pela régua de Sijtsma e Molenaar (2002); Hi de ",
           .pr_num(min(Hi), 2, FALSE, decimal), " a ", .pr_num(max(Hi), 2, FALSE, decimal), ".")
  .ee_resultado("Escala de Mokken (TRI não paramétrica)", n = n, n_excluidos = pp$n_excluidos, tabela = tab,
                efeito = data.frame(indice = "H da escala", valor = H, ic_inf = ic_H[1], ic_sup = ic_H[2], principal = TRUE),
                decisoes = c("Coeficientes de escalabilidade H, Hi e Hij com EP delta (mokken::coefH); IC = ±1,96 EP.",
                             "Régua de Sijtsma e Molenaar (2002, p. 60), pelo capítulo 15 do Wiley: < ,30 não escalável; ,30 fraca; ,40 média; ,50 forte.", decisoes),
                avisos = avisos, apa = apa,
                extra = list(tabelas_extra = c(if (!is.null(mono)) list(`Monotonicidade manifesta` = mono)), H = H, ep_H = seH, Hij = h$Hij, monotonicidade = mono, iio = io))
}

# ---------------------------------------------------------------- CAT pós-hoc (v1.7)
# simular_cat(): simulação pós-hoc de testagem adaptativa (Martinková & Hladká, 2023, 10.2.6;
# Petersen, 2024, cap. 20) com o mirtCAT: calibra o modelo unidimensional nos dados (ou usa
# o ajuste dado), e para cada pessoa refaz o teste como se fosse adaptativo com as respostas
# que ela deu: primeiro item de máxima informação em θ = 0, θ por MAP, próximo item por
# máxima informação (MI) e parada quando o EP de θ cai abaixo de `ep_max` (padrão ,30, o do
# livro, confiabilidade perto de ,91) ou ao atingir `max_itens`. Sai só agregado: itens
# usados, correlação e RMSE entre o θ do CAT e o do teste completo, % que atingiu o EP e a
# exposição de cada item. Nenhum θ por pessoa sai no objeto (regra de privacidade do módulo).
simular_cat <- function(dados, itens, modelo = "2PL", fit = NULL, ep_max = .30, max_itens = NULL, metodo = "MAP",
                        criterio = "MI", semente = 2026, decimal = ".", idioma = "pt") {
  .pr_exige("mirt", "mirtCAT")
  pp <- .tri_prepara(as.data.frame(dados), itens, "listwise")
  d <- as.data.frame(lapply(pp$d, .tri_recodifica))
  if (is.null(fit)) {
    modelo <- match.arg(modelo, .tri_modelos)
    a <- .tri_mirt(d, .tri_itemtype(modelo), semente)
    if (is.null(a$fit)) stop("o mirt não ajustou o modelo: ", a$erro, call. = FALSE)
    fit <- a$fit
  }
  if (mirt::extract.mirt(fit, "nfact") != 1) stop("simular_cat() é para modelo unidimensional", call. = FALSE)
  desenho <- list(min_SEM = ep_max); if (!is.null(max_itens)) desenho$max_items <- max_itens
  set.seed(semente)
  sim <- suppressMessages(mirtCAT::mirtCAT(mo = fit, local_pattern = as.matrix(d), start_item = "MI", method = metodo,
                                           criteria = criterio, design = desenho))
  if (!is.null(sim$items_answered)) sim <- list(sim)
  th_cat <- vapply(sim, function(s) as.numeric(utils::tail(s$thetas, 1)), 0)
  ep_cat <- vapply(sim, function(s) as.numeric(utils::tail(s$SE_thetas, 1)), 0)
  n_it <- vapply(sim, function(s) length(s$items_answered), 0L)
  th_full <- mirt::fscores(fit, method = metodo, full.scores.SE = TRUE, response.pattern = as.matrix(d))
  th_c <- th_full[, "F1"]
  expo <- tabulate(unlist(lapply(sim, `[[`, "items_answered")), nbins = length(itens)) / length(sim)
  k <- length(itens); n <- length(sim)
  resumo <- data.frame(pessoas = n, itens_banco = k, itens_media = mean(n_it), itens_mediana = stats::median(n_it), itens_min = min(n_it), itens_max = max(n_it),
                       pct_itens_economizados = 100 * (1 - mean(n_it) / k), pct_atingiu_ep = 100 * mean(ep_cat <= ep_max + 1e-8),
                       r_theta_completo = stats::cor(th_cat, th_c), rmse_theta = sqrt(mean((th_cat - th_c)^2)), vies_theta = mean(th_cat - th_c))
  expos <- data.frame(item = itens, exposicao = expo, row.names = NULL)
  avisos <- character()
  if (resumo$pct_atingiu_ep < 100) avisos <- c(avisos, paste0(round(100 - resumo$pct_atingiu_ep, 1), "% das pessoas terminaram sem atingir EP <= ", .pr_num(ep_max, 2, FALSE, decimal),
    ": o banco não tem informação suficiente nos extremos do traço (conferir a curva de informação no tri_ajustar())."))
  if (max(expo) > .8) avisos <- c(avisos, paste0("Exposição acima de 80% em ", paste(itens[expo > .8], collapse = ", "), ": num CAT real, controlar a exposição (por exemplo, randomesque) para não gastar esses itens."))
  if (k < 30) avisos <- c(avisos, paste0("Banco de ", k, " itens: CAT costuma valer a pena com bancos maiores; com poucos itens, a forma curta fixa pode servir igual."))
  avisos <- c(avisos, "Simulação pós-hoc com as respostas que as pessoas deram ao teste completo: superestima a concordância com o teste inteiro (mesmas respostas) e não mede efeitos de ordem ou de contexto de um CAT real.")
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  apa <- if (en) paste0("A post-hoc CAT simulation (mirtCAT; MAP, maximum-information selection, stop at SE ≤ ", .pr_num(ep_max, 2, FALSE, decimal), ") used on average ",
                        .pr_num(resumo$itens_media, 1, TRUE, decimal), " of ", k, " items; CAT and full-test θ correlated r = ", .pr_num(resumo$r_theta_completo, 2, FALSE, decimal),
                        " (RMSE = ", .pr_num(resumo$rmse_theta, 2, TRUE, decimal), ").") else
    paste0("A simulação pós-hoc de CAT (mirtCAT; MAP, seleção por máxima informação, parada com EP ≤ ", .pr_num(ep_max, 2, FALSE, decimal), ") usou em média ",
           .pr_num(resumo$itens_media, 1, TRUE, decimal), " dos ", k, " itens; o θ do CAT e o do teste completo correlacionaram r = ", .pr_num(resumo$r_theta_completo, 2, FALSE, decimal),
           sep, "RMSE = ", .pr_num(resumo$rmse_theta, 2, TRUE, decimal), ".")
  .ee_resultado("Simulação pós-hoc de testagem adaptativa (CAT)", n = n, n_excluidos = pp$n_excluidos, tabela = resumo,
                decisoes = c(paste0("mirtCAT pós-hoc: primeiro item por máxima informação em θ = 0, θ por ", metodo, ", seleção ", criterio, ", parada com EP <= ", .pr_num(ep_max, 2, FALSE, decimal),
                                    if (!is.null(max_itens)) paste0(" ou ", max_itens, " itens") else "", " (Martinková & Hladká, 2023, 10.2.6)."),
                             paste0("Comparação com o θ por ", metodo, " do teste completo no mesmo modelo; exposição = proporção de pessoas que receberam o item.")),
                avisos = avisos, apa = apa,
                extra = list(tabelas_extra = list(`Exposição dos itens` = expos), modelo = fit,
                             distribuicao_itens = table(n_it)))
}
