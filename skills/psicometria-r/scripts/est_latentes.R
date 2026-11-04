# est_latentes.R: análise de classes latentes (LCA, itens categóricos) e de perfis
# latentes (LPA, indicadores contínuos), v1.4; três passos BCH na v1.5. Carregar por scripts/carregar.R.
#
#   classes_latentes()  poLCA, de 1 a k_max classes, com várias partidas aleatórias
#   perfis_latentes()   mclust (modelo EEI por padrão, o "modelo 1" do tidyLPA, o do
#                       curso), de 1 a k_max perfis
#
# Critérios do curso (módulo de classes e perfis latentes): menores AIC, BIC e SABIC
# (comparativos, tendem a cair sempre com mais classes); entropia relativa >= .80,
# probabilidade posterior média da classe mais provável e tamanho da menor classe
# (>= 5% da amostra) como critérios de adequação, não de comparação; testes de razão de
# verossimilhança entre k - 1 e k classes (LMR ajustado; BLRT, mais robusto e quase
# sempre significativo na prática, segundo o curso). A decisão final junta os índices,
# a interpretabilidade e a teoria; a skill não escolhe sozinha.
# Fórmulas (as do tidyLPA): AIC = -2LL + 2p; BIC = -2LL + p ln(n); SABIC = -2LL + p
# ln((n + 2)/24); entropia = 1 + Σ Σ z ln z / (n ln K); LMR = 2ΔLL / (1 + 1/((3K - 1 -
# (3(K - 1) - 1)) ln n)), p pela qui-quadrado com Δp gl (Lo, Mendell & Rubin, 2001,
# na implementação do tidyLPA::calc_lrt).
# prob_min e prob_max: as do tidyLPA (soma da posterior de k entre os alocados em k sobre
# a soma total da posterior de k); avepp_min: menor média da posterior entre os alocados
# em cada classe (a tabela "Average Latent Class Probabilities" do Mplus).
# Privacidade: a classe de cada pessoa fica em $classe (não impressa) para análise
# seguinte no computador do usuário; a impressão traz só agregados.

.lt_entropia <- function(z) {
  if (is.null(dim(z)) || ncol(z) == 1) return(1)
  1 + sum(z * log(z + 1e-12)) / (nrow(z) * log(ncol(z)))
}

.lt_avepp <- function(z, classe) {
  if (is.null(dim(z)) || ncol(z) == 1) return(c(1, 1))
  d <- vapply(seq_len(ncol(z)), function(k) if (any(classe == k)) mean(z[classe == k, k]) else NA_real_, 0)
  range(d, na.rm = TRUE)
}

# Probabilidade de classificação correta por classe, como tidyLPA (prob_min, prob_max):
# soma da posterior da classe k entre os alocados em k sobre a soma da posterior de k.
.lt_prob_tidylpa <- function(z, classe) {
  if (is.null(dim(z)) || ncol(z) == 1) return(c(1, 1))
  d <- vapply(seq_len(ncol(z)), function(k) sum(z[classe == k, k]) / sum(z[, k]), 0)
  range(d, na.rm = TRUE)
}

.lt_lmr <- function(ll0, p0, k0, ll1, p1, k1, n) {
  lr <- 2 * (ll1 - ll0)
  mod <- lr / (1 + (((3 * k1 - 1) - (3 * k0 - 1)) * log(n))^-1)
  c(lmr = mod, gl = p1 - p0, p = stats::pchisq(mod, p1 - p0, lower.tail = FALSE))
}

.lt_tabela <- function(fits, n) {
  tab <- do.call(rbind, lapply(fits, function(f) {
    pp <- .lt_avepp(f$z, f$classe); pt <- .lt_prob_tidylpa(f$z, f$classe); prop <- tabulate(f$classe, f$k) / n
    data.frame(classes = f$k, ll = f$ll, parametros = f$p, aic = -2 * f$ll + 2 * f$p, bic = -2 * f$ll + f$p * log(n),
               sabic = -2 * f$ll + f$p * log((n + 2) / 24), entropia = .lt_entropia(f$z), prob_min = pt[1], prob_max = pt[2], avepp_min = pp[1],
               menor_classe = min(prop), maior_classe = max(prop), convergiu = f$convergiu, stringsAsFactors = FALSE)
  }))
  tab$lmr <- NA_real_; tab$lmr_p <- NA_real_
  for (i in seq_len(nrow(tab))[-1]) {
    t <- .lt_lmr(tab$ll[i - 1], tab$parametros[i - 1], tab$classes[i - 1], tab$ll[i], tab$parametros[i], tab$classes[i], n)
    tab$lmr[i] <- t[["lmr"]]; tab$lmr_p[i] <- t[["p"]]
  }
  tab
}

.lt_leitura <- function(tab, rotulo, decimal, idioma) {
  f2 <- function(x) .pr_num(x, 2, FALSE, decimal)
  ok <- tab$entropia >= .80 & tab$menor_classe >= .05
  ok[tab$classes == 1] <- TRUE
  cand <- unique(c(tab$classes[which.min(tab$bic)], tab$classes[which.min(tab$aic)], tab$classes[which.min(tab$sabic)]))
  avisos <- character()
  if (which.min(tab$bic) == nrow(tab)) avisos <- c(avisos, paste0("O BIC ainda cai na maior solução testada (", max(tab$classes), " ", rotulo,
                                                                   "): ampliar k_max antes de decidir."))
  if (any(!tab$convergiu)) avisos <- c(avisos, paste0("Solução sem convergência ou sem replicar o máximo da verossimilhança em: ",
                                                      paste(tab$classes[!tab$convergiu], collapse = ", "), " ", rotulo, "; aumentar nrep ou maxiter."))
  pequenas <- tab$classes[tab$menor_classe < .05 & tab$classes > 1]
  if (length(pequenas)) avisos <- c(avisos, paste0("Menor classe abaixo de 5% da amostra em ", paste(pequenas, collapse = ", "), " ", rotulo,
                                                   " (critério do curso): classe instável, difícil de replicar."))
  list(candidatos = cand, adequadas = tab$classes[ok], avisos = avisos)
}

classes_latentes <- function(dados, itens, k_max = 5, nrep = 10, maxiter = 3000, semente = 2026, decimal = ".", idioma = "pt") {
  .pr_exige("poLCA")
  dd <- .ee_dados(dados, itens)
  d <- dd$d; n <- nrow(d)
  # poLCA pede categorias inteiras 1, 2, ...: recodifica pela ordem dos valores
  niveis <- lapply(d, function(x) sort(unique(x)))
  d1 <- as.data.frame(Map(function(x, lv) match(x, lv), d, niveis))
  if (any(lengths(niveis) < 2)) stop("item constante: ", paste(itens[lengths(niveis) < 2], collapse = ", "))
  if (any(lengths(niveis) > 10)) stop("item com mais de 10 categorias: LCA é para itens categóricos; para indicadores contínuos, perfis_latentes()")
  f <- stats::as.formula(paste0("cbind(", paste(sprintf("`%s`", itens), collapse = ", "), ") ~ 1"))
  restaura <- .ee_construcao_semente(semente); on.exit(restaura(), add = TRUE)
  fits <- lapply(seq_len(k_max), function(k) {
    av <- character()
    m <- withCallingHandlers(poLCA::poLCA(f, d1, nclass = k, nrep = if (k == 1) 1 else nrep, maxiter = maxiter, verbose = FALSE, calc.se = FALSE),
                             warning = function(w) { av <<- c(av, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(k = k, ll = m$llik, p = m$npar, z = if (k == 1) matrix(1, n, 1) else m$posterior, classe = if (k == 1) rep(1L, n) else m$predclass,
         convergiu = m$numiter < maxiter && !any(grepl("MAXIMUM LIKELIHOOD NOT FOUND", av, ignore.case = TRUE)), modelo = m)
  })
  tab <- .lt_tabela(fits, n)
  lt <- .lt_leitura(tab, "classes", decimal, idioma)
  perfis <- lapply(fits[-1], function(f) {
    pr <- f$modelo$probs
    do.call(rbind, lapply(names(pr), function(it) {
      m <- pr[[it]]; colnames(m) <- as.character(niveis[[it]])[seq_len(ncol(m))]
      data.frame(item = it, classe = seq_len(nrow(m)), m, check.names = FALSE, row.names = NULL)
    }))
  })
  names(perfis) <- paste0(seq_len(k_max)[-1], " classes")
  decisoes <- c(paste0("poLCA de 1 a ", k_max, " classes, ", nrep, " partidas aleatórias por solução (o maior log-verossimilhança fica), maxiter ", maxiter,
                       ", semente ", semente, "; categorias recodificadas para 1, 2, ... na ordem dos valores."),
                "Critérios do curso: menores AIC, BIC e SABIC; entropia >= .80, menor classe >= 5% e probabilidade posterior média alta como adequação; LMR entre k - 1 e k.",
                paste0("Candidatas pelos critérios de informação: ", paste(lt$candidatos, collapse = ", "), "; soluções que passam em entropia e menor classe: ",
                       paste(lt$adequadas, collapse = ", "), ". Decidir pela interpretabilidade das probabilidades de resposta ($perfis) e pela teoria."),
                "BLRT não calculado aqui (custo alto); o curso o considera o mais robusto, mas quase sempre significativo.")
  apa <- paste0(if (identical(idioma, "en")) "Latent class models with 1 to " else "Modelos de classes latentes de 1 a ", k_max,
                if (identical(idioma, "en")) " classes were compared (poLCA; " else " classes foram comparados (poLCA; ", nrep,
                if (identical(idioma, "en")) " random starts). The lowest BIC was for the " else " partidas aleatórias). O menor BIC foi o da solução de ",
                tab$classes[which.min(tab$bic)], if (identical(idioma, "en")) "-class solution (entropy = " else " classes (entropia = ",
                .pr_num(tab$entropia[which.min(tab$bic)], 2, FALSE, decimal),
                if (identical(idioma, "en")) ", smallest class = " else ", menor classe = ",
                .pr_num(100 * tab$menor_classe[which.min(tab$bic)], 1, TRUE, decimal), "%).")
  .ee_resultado("Análise de classes latentes", n = n, n_excluidos = dd$n_excluidos, tabela = tab, decisoes = decisoes, avisos = lt$avisos, apa = apa,
                extra = list(perfis = perfis, classe = lapply(fits, `[[`, "classe"), modelos = lapply(fits, `[[`, "modelo"),
                             posterior = lapply(fits, `[[`, "z"), linhas = rownames(d), candidatos = lt$candidatos))
}

# Parâmetros do LPA contados como no tidyLPA (calc_fitindices.Mclust).
.lt_npar_mclust <- function(modelo, k, v) {
  (k - 1) + v * k + switch(modelo, EEI = v, VVI = v * k, EEE = v + v * (v - 1), VVV = v * k + v * (v - 1) * k, v)
}

# mclust::Mclust e mclustBootstrapLRT avaliam a chamada de mclustBIC no ambiente de quem
# chamou: sem o mclust anexado, "could not find function mclustBIC". O invólucro põe a
# função no ambiente da chamada, sem anexar o pacote.
.lt_mclust <- function(x, G, modelNames) {
  mclustBIC <- mclust::mclustBIC
  mclust::Mclust(x, G = G, modelNames = modelNames, verbose = FALSE)
}
.lt_mclust_blrt <- function(x, modelName, maxG, nboot) {
  mclustBIC <- mclust::mclustBIC
  mclust::mclustBootstrapLRT(x, modelName = modelName, maxG = maxG, nboot = nboot, verbose = FALSE)
}

perfis_latentes <- function(dados, vars, k_max = 5, modelo = c("EEI", "VVI", "EEE", "VVV"), padronizar = FALSE, blrt = 0,
                            semente = 2026, decimal = ".", idioma = "pt") {
  modelo <- match.arg(modelo)
  .pr_exige("mclust")
  dd <- .ee_dados(dados, vars, numericas = vars)
  d <- dd$d; n <- nrow(d); v <- length(vars)
  x <- if (padronizar) scale(d) else as.matrix(d)
  restaura <- .ee_construcao_semente(semente); on.exit(restaura(), add = TRUE)
  fits <- lapply(seq_len(k_max), function(k) {
    m <- suppressWarnings(.lt_mclust(x, G = k, modelNames = modelo))
    if (is.null(m) && k == 1) m <- suppressWarnings(.lt_mclust(x, G = 1, modelNames = NULL))
    if (is.null(m)) return(list(k = k, ll = NA_real_, p = NA_real_, z = matrix(NA, n, k), classe = rep(NA_integer_, n), convergiu = FALSE, modelo = NULL))
    list(k = k, ll = m$loglik, p = .lt_npar_mclust(modelo, k, v), z = m$z, classe = m$classification, convergiu = TRUE, modelo = m)
  })
  # com k = 1 o mclust ajusta XII/XXX; a contagem de parâmetros segue o modelo pedido (tidyLPA)
  ok <- vapply(fits, function(f) !is.na(f$ll), TRUE)
  tab <- .lt_tabela(fits[ok], n)
  if (blrt > 0) {
    b <- .lt_mclust_blrt(x, modelName = modelo, maxG = max(tab$classes) - 1, nboot = blrt)
    tab$blrt_p <- c(NA, b$p.value)[seq_len(nrow(tab))]
  }
  lt <- .lt_leitura(tab, "perfis", decimal, idioma)
  medias <- lapply(fits[ok][-1], function(f) {
    mu <- f$modelo$parameters$mean
    data.frame(perfil = seq_len(ncol(mu)), proporcao = tabulate(f$classe, f$k) / n, t(mu), check.names = FALSE, row.names = NULL)
  })
  names(medias) <- paste0(tab$classes[-1], " perfis")
  decisoes <- c(paste0("mclust, modelo ", modelo, c(EEI = " (variâncias iguais entre perfis, covariâncias zero: o modelo 1 do tidyLPA, usado no curso)",
                                                    VVI = " (variâncias livres entre perfis, covariâncias zero)", EEE = " (variâncias e covariâncias iguais)",
                                                    VVV = " (variâncias e covariâncias livres)")[[modelo]],
                       ", de 1 a ", k_max, " perfis", if (padronizar) ", indicadores padronizados" else "", "."),
                "Critérios do curso: menores AIC, BIC e SABIC; entropia >= .80 e menor perfil >= 5% como adequação; LMR e, se pedido, BLRT entre k - 1 e k.",
                paste0("Candidatos pelos critérios de informação: ", paste(lt$candidatos, collapse = ", "), "; soluções que passam em entropia e menor perfil: ",
                       paste(lt$adequadas, collapse = ", "), ". Médias por perfil em $medias."))
  apa <- paste0(if (identical(idioma, "en")) "Latent profile models (mclust, " else "Modelos de perfis latentes (mclust, ", modelo,
                if (identical(idioma, "en")) ") with 1 to " else ") de 1 a ", k_max,
                if (identical(idioma, "en")) " profiles were compared; the lowest BIC was for the " else " perfis foram comparados; o menor BIC foi o da solução de ",
                tab$classes[which.min(tab$bic)], if (identical(idioma, "en")) "-profile solution (entropy = " else " perfis (entropia = ",
                .pr_num(tab$entropia[which.min(tab$bic)], 2, FALSE, decimal), if (identical(idioma, "en")) ", smallest profile = " else ", menor perfil = ",
                .pr_num(100 * tab$menor_classe[which.min(tab$bic)], 1, TRUE, decimal), "%).")
  names(tab)[names(tab) == "classes"] <- "perfis"
  .ee_resultado("Análise de perfis latentes", n = n, n_excluidos = dd$n_excluidos, tabela = tab, decisoes = decisoes, avisos = lt$avisos, apa = apa,
                extra = list(medias = medias, classe = lapply(fits, `[[`, "classe"), modelos = lapply(fits, `[[`, "modelo"),
                             posterior = lapply(fits, `[[`, "z"), linhas = rownames(d), candidatos = lt$candidatos))
}

# ---------------------------------------------------------------- v1.5: três passos (BCH)
# classes_distal(): relação entre as classes (ou perfis) e uma variável externa sem o viés
# da atribuição modal, pela correção de Bolck, Croon e Hagenaars (2004) na versão de
# Vermunt (2010), a abordagem de três passos que o curso roda no tidySEM. Passo 1: o
# modelo de classes de classes_latentes() ou perfis_latentes(). Passo 2: atribuição modal
# e matriz D de erro de classificação, D[t, s] = P(atribuído a s | classe verdadeira t),
# estimada pelas probabilidades posteriores. Passo 3: pesos BCH w*[i, t] = Σs W[i, s]
# (D⁻¹)[s, t] e médias (ou
# proporções) por classe ponderadas por eles, com EP sanduíche por pessoa (a incerteza do
# passo 1 fica de fora, como no tidySEM e no Latent GOLD). Pesos negativos são esperados.
classes_distal <- function(resultado, dados, distal, k, tipo = c("auto", "continua", "categorica"), nivel = .95, decimal = ".", idioma = "pt") {
  tipo <- match.arg(tipo)
  if (!inherits(resultado, "pr_resultado") || is.null(resultado$posterior)) stop("resultado: saída de classes_latentes() ou perfis_latentes() da v1.5")
  P <- resultado$posterior[[k]]
  if (is.null(P) || k < 2) stop("k: número de classes da solução escolhida (2 ou mais)")
  linhas <- resultado$linhas
  z <- dados[linhas, distal]
  ok <- !is.na(z); P <- P[ok, , drop = FALSE]; z <- z[ok]; n <- length(z)
  if (tipo == "auto") tipo <- if (is.numeric(z) && length(unique(z)) > 7) "continua" else "categorica"
  W <- diag(k)[max.col(P, ties.method = "first"), , drop = FALSE]
  if (any(colSums(W) == 0)) stop("classe sem nenhum caso na atribuição modal: a correção BCH não se define")
  # D[t, s] = P(atribuído a s | classe verdadeira t) = Σi P[i, t] W[i, s] / Σi P[i, t]
  D <- t(P) %*% W / colSums(P)
  wb <- W %*% solve(D)
  ent <- 1 + sum(P * log(P + 1e-300)) / (n * log(k))
  medias_ponderadas <- function(y) {
    m <- colSums(wb * y) / colSums(wb)
    U <- sweep(wb, 2, colSums(wb), "/") * outer(y, m, "-")
    list(est = m, V = crossprod(U))
  }
  wald <- function(est, V) {
    C <- cbind(matrix(1, k - 1, 1), -diag(k - 1))
    q <- drop(t(C %*% est) %*% solve(C %*% V %*% t(C)) %*% (C %*% est))
    c(qui2 = q, gl = k - 1, p = stats::pchisq(q, k - 1, lower.tail = FALSE))
  }
  zc <- stats::qnorm(1 - (1 - nivel) / 2)
  ingenuo <- tapply(if (tipo == "continua") z else rep(NA, n), max.col(P, ties.method = "first"), mean)
  if (tipo == "continua") {
    r <- medias_ponderadas(z); ep <- sqrt(diag(r$V)); wt <- wald(r$est, r$V)
    tab <- data.frame(classe = seq_len(k), n_modal = colSums(W), media_bch = r$est, ep = ep, ic_inf = r$est - zc * ep, ic_sup = r$est + zc * ep,
                      media_modal = as.numeric(ingenuo[as.character(seq_len(k))]), row.names = NULL)
    pares <- utils::combn(k, 2)
    ph <- data.frame(comparacao = paste(pares[1, ], "vs", pares[2, ]), diferenca = r$est[pares[1, ]] - r$est[pares[2, ]],
                     p = apply(pares, 2, function(pp) { cc <- numeric(k); cc[pp[1]] <- 1; cc[pp[2]] <- -1
                       stats::pchisq(sum(cc * r$est)^2 / drop(t(cc) %*% r$V %*% cc), 1, lower.tail = FALSE) }), row.names = NULL)
    ph$p_holm <- stats::p.adjust(ph$p, "holm")
  } else {
    niv <- sort(unique(z))
    tabs <- lapply(niv, function(v) { r <- medias_ponderadas(as.numeric(z == v)); list(est = r$est, V = r$V) })
    tab <- do.call(rbind, lapply(seq_along(niv), function(j) data.frame(categoria = niv[j], classe = seq_len(k), proporcao_bch = tabs[[j]]$est,
                                                                       ep = sqrt(diag(tabs[[j]]$V)), row.names = NULL)))
    wt <- wald(tabs[[1]]$est, tabs[[1]]$V)
    if (length(niv) > 2) {
      ws <- t(vapply(tabs[-length(tabs)], function(x) wald(x$est, x$V), numeric(3)))
      wt <- c(qui2 = NA, gl = NA, p = NA)
      ests <- unlist(lapply(tabs[-length(tabs)], `[[`, "est"))
      U <- do.call(cbind, lapply(niv[-length(niv)], function(v) { y <- as.numeric(z == v); m <- colSums(wb * y) / colSums(wb)
        sweep(wb, 2, colSums(wb), "/") * outer(y, m, "-") }))
      Vb <- crossprod(U); C <- kronecker(diag(length(niv) - 1), cbind(matrix(1, k - 1, 1), -diag(k - 1)))
      q <- drop(t(C %*% ests) %*% MASS::ginv(C %*% Vb %*% t(C)) %*% (C %*% ests)); g <- (length(niv) - 1) * (k - 1)
      wt <- c(qui2 = q, gl = g, p = stats::pchisq(q, g, lower.tail = FALSE))
    }
    ph <- NULL
  }
  en <- identical(idioma, "en")
  apa <- if (en) paste0("Using the three-step BCH approach (Vermunt, 2010), ", distal, " differed across the ", k, " classes: Wald χ²(", wt[["gl"]], ") = ",
                        .pr_num(wt[["qui2"]], 2, TRUE, decimal), ", ", .pr_p(wt[["p"]], decimal), ".") else
    paste0("Pela abordagem de três passos com correção BCH (Vermunt, 2010), a comparação de ", distal, " entre as ", k, " classes deu Wald χ²(", wt[["gl"]], ") = ",
           .pr_num(wt[["qui2"]], 2, TRUE, decimal), ", ", .pr_p(wt[["p"]], decimal), ".")
  avisos <- character()
  if (ent < .6) avisos <- c(avisos, paste0("Entropia ", round(ent, 2), ": com separação baixa a correção BCH fica instável e pode continuar enviesada (o curso avisa)."))
  if (tipo == "categorica" && any(tab$proporcao_bch < 0 | tab$proporcao_bch > 1)) avisos <- c(avisos, "Proporção BCH fora de [0, 1] (pesos negativos com classe pequena): ler com cautela.")
  .ee_resultado(paste0("Classes e variável externa (BCH, ", k, " classes): ", distal), n = n, n_excluidos = sum(!ok), tabela = tab, posthoc = ph,
                efeito = data.frame(indice = "Wald χ² entre classes", valor = wt[["qui2"]], ic_inf = NA_real_, ic_sup = NA_real_, principal = TRUE),
                decisoes = c("Três passos com correção BCH (Bolck et al., 2004; Vermunt, 2010): variável externa fora do modelo de medida, sem mudar as classes.",
                             "EP sanduíche por pessoa, sem a incerteza do passo 1; comparações par a par por Wald com Holm.",
                             if (tipo == "continua") "A coluna media_modal (atribuição modal ingênua) fica ao lado só para mostrar o viés que a correção remove; não relatar." else character()),
                avisos = avisos, apa = apa, extra = list(matriz_D = D, entropia = ent))
}
