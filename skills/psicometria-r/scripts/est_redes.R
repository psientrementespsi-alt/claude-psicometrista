# est_redes.R: redes psicométricas (estimação, precisão e estabilidade, ponte e
# previsibilidade). Carregar por scripts/carregar.R (usa .ee_resultado, .ee_dados,
# .ee_junta e .ee_ic de est_comum.R e .pr_num, .pr_exige de psicometria.R).
#
# Funções públicas:
#   rede_psicometrica(dados, itens, estimador = "EBICglasso", tuning = 0.5, cor = "cor_auto",
#                     ausentes = "pairwise", tamanho_amostra = "maximum", tipos = NULL,
#                     redundancia = TRUE, n_fortes = 10, decimal = ".", idioma = "pt", ...)
#   estabilidade_rede(rede, n_boot = 1000, indices = c("strength", "expectedInfluence"),
#                     semente = 2026, n_boot_casos = n_boot, caso_min = 0.05, caso_max = 0.75,
#                     caso_n = 10, executar = c("arestas", "casos"), comunidades = NULL,
#                     n_fortes = 10, n_nucleos = 1, decimal = ".", idioma = "pt")
#   ponte(rede, comunidades, decimal = ".", idioma = "pt")
#   previsibilidade(rede, tuning = NULL, decimal = ".", idioma = "pt")
#   rede_temporal(dados, vars, id, tempo, dia, bip, destendenciar, gamma, aleatorios, alfa)   (v1.7)
#
# Padrões: os do módulo Network do JASP 0.95.1, que o livro do usuário usa:
# EBICglasso, correlação "Auto" (cor_auto), γ = 0,5, ausentes par a par, tamanho de
# amostra "Maximum", centralidade "Normalized" (escores z). Os padrões do bootnet são
# outros (corMethod = "cor", sampleSize = "pairwise_average"), por isso tudo vai
# explícito na chamada. Bootstrap não paramétrico de arestas com 1000 reamostras
# e bootstrap de descarte de casos com CS por índice.

# ---------------------------------------------------------------- auxiliares
.ee_redes_rotulos_indice <- list(
  strength = c(pt = "força", en = "strength"),
  expectedInfluence = c(pt = "influência esperada", en = "expected influence"),
  closeness = c(pt = "proximidade", en = "closeness"),
  betweenness = c(pt = "intermediação", en = "betweenness"),
  bridgeStrength = c(pt = "força de ponte", en = "bridge strength"),
  bridgeExpectedInfluence = c(pt = "influência esperada de ponte", en = "bridge expected influence"),
  bridgeCloseness = c(pt = "proximidade de ponte", en = "bridge closeness"),
  bridgeBetweenness = c(pt = "intermediação de ponte", en = "bridge betweenness"))

.ee_redes_rotulo <- function(ind, idioma = "pt") {
  vapply(ind, function(i) {
    r <- .ee_redes_rotulos_indice[[i]]
    if (is.null(r)) i else unname(r[[if (idioma == "en") "en" else "pt"]])
  }, character(1), USE.NAMES = FALSE)
}

# Aresta como no livro (N1–N2, traço curto).
.ee_redes_aresta <- function(a, b) paste0(a, "–", b)

# Normalização "Normalized" do JASP: (x - média) / DP amostral, por coluna; com
# variância zero só centra (jaspNetwork, .networkAnalysisComputeCentrality).
.ee_redes_z <- function(x) {
  ok <- !is.na(x)
  if (!any(ok)) return(x)
  o <- x[ok] - mean(x[ok])
  s <- sqrt(sum(o^2) / max(1, length(o) - 1))
  x[ok] <- if (s == 0) o else o / s
  x
}

# Centralidade como o JASP: qgraph::centrality com pesos e sinal, sem todos os
# caminhos mínimos; intermediação dividida por 2 em rede não direcionada.
.ee_redes_centralidade <- function(W) {
  ct <- qgraph::centrality(W, weighted = TRUE, signed = TRUE, all.shortest.paths = FALSE)
  nos <- colnames(W)
  bruta <- data.frame(no = nos,
                      forca = unname(ct$InDegree),
                      influencia_esperada = unname(ct$InExpectedInfluence),
                      proximidade = unname(ct$Closeness),
                      intermediacao = unname(ct$Betweenness) / 2,
                      stringsAsFactors = FALSE)
  z <- bruta
  for (v in c("forca", "influencia_esperada", "proximidade", "intermediacao")) z[[v]] <- .ee_redes_z(bruta[[v]])
  list(bruta = bruta, z = z)
}

.ee_redes_arestas <- function(W) {
  ij <- which(upper.tri(W), arr.ind = TRUE)
  ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
  nos <- colnames(W)
  data.frame(aresta = .ee_redes_aresta(nos[ij[, 1]], nos[ij[, 2]]),
             no1 = nos[ij[, 1]], no2 = nos[ij[, 2]], peso = W[ij],
             stringsAsFactors = FALSE)
}

# Roda expr guardando os avisos (o print do resultado os mostra em português).
.ee_redes_capturar <- function(expr) {
  av <- character()
  valor <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    av <<- c(av, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(valor = valor, avisos = unique(av))
}

.ee_redes_traduz_aviso <- function(msg) {
  vapply(msg, function(m) {
    if (grepl("dense regularized network", m, fixed = TRUE))
      return(paste0("O EBICglasso escolheu uma rede densa (lambda < 0,1 vez o lambda máximo); o qgraph avisa ",
                    "possível perda de especificidade nas arestas menores e sugere threshold = TRUE, que ",
                    "ganha especificidade à custa de sensibilidade. Há sempre troca entre descoberta e cautela ",
                    "(Isvoranu & Epskamp, 2023): decidir pela pergunta de pesquisa, antes de ver o resultado, e declarar. ",
                    "As arestas de menor peso pedem leitura cautelosa."))
    paste0("bootnet/qgraph: ", m)
  }, character(1), USE.NAMES = FALSE)
}

# Mesmo cálculo de bootnet:::sampleSize_pairwise (R/defaultFunctions.R, bootnet 1.9.1).
.ee_redes_n_amostra <- function(d, tipo) {
  if (tipo == "maximum") return(sum(rowSums(!is.na(d)) > 0))
  if (tipo == "minimum") return(sum(stats::complete.cases(d)))
  m <- crossprod(!is.na(as.matrix(d)))
  v <- m[lower.tri(m)]
  switch(tipo, pairwise_maximum = max(v), pairwise_minimum = min(v), pairwise_average = mean(v),
         stop("tamanho_amostra não reconhecido: ", tipo))
}

# Descreve a correlação que entra no GGM.
.ee_redes_cor_txt <- function(d, cor, idioma) {
  if (cor == "cor_auto") {
    ncat <- vapply(d, function(x) length(unique(stats::na.omit(x))), numeric(1))
    nord <- sum(ncat <= 7)  # regra do qgraph::cor_auto (ordinalLevelMax = 7)
    if (idioma == "en") {
      if (nord == length(ncat)) return("polychoric correlations (cor_auto)")
      if (nord == 0) return("Pearson correlations (cor_auto; all items with more than 7 categories)")
      return("polychoric, polyserial and Pearson correlations according to the number of categories (cor_auto)")
    }
    if (nord == length(ncat)) return("correlações policóricas (cor_auto)")
    if (nord == 0) return("correlações de Pearson (cor_auto; todos os itens com mais de 7 categorias)")
    return("correlações policóricas, polisseriais e de Pearson conforme o número de categorias (cor_auto)")
  }
  txt <- list(cor = c("correlações de Pearson", "Pearson correlations"),
              spearman = c("correlações de Spearman", "Spearman correlations"),
              npn = c("correlações após transformação nonparanormal", "correlations after the nonparanormal transformation"),
              cov = c("covariâncias", "covariances"),
              cor_mantar = c("correlações pelo método mantar", "correlations by the mantar method"))[[cor]]
  if (is.null(txt)) return(cor)
  txt[[if (idioma == "en") 2 else 1]]
}

# Aceita o resultado de rede_psicometrica(), um objeto do bootnet::estimateNetwork()
# ou, em ponte(), uma matriz de pesos simétrica.
.ee_redes_obter <- function(rede, aceita_matriz = FALSE) {
  if (inherits(rede, "pr_resultado") && !is.null(rede$modelo)) {
    W <- rede$pesos
    return(list(modelo = rede$modelo, W = W, estimador = rede$estimador, tuning = rede$tuning,
                tipos = rede$tipos, n = rede$n))
  }
  if (inherits(rede, "bootnetResult")) {
    W <- rede$graph
    if (is.null(colnames(W))) dimnames(W) <- list(rede$labels, rede$labels)
    tun <- rede$arguments$tuning
    return(list(modelo = rede, W = W, estimador = rede$default, tuning = if (is.null(tun)) NA else tun,
                tipos = rede$arguments$type, n = rede$nPerson))
  }
  if (aceita_matriz && is.matrix(rede) && is.numeric(rede)) {
    if (nrow(rede) != ncol(rede) || !isSymmetric(unname(rede), tol = 1e-10)) stop("a matriz de pesos precisa ser quadrada e simétrica")
    if (is.null(colnames(rede))) dimnames(rede) <- list(paste0("V", seq_len(ncol(rede))), paste0("V", seq_len(ncol(rede))))
    return(list(modelo = NULL, W = rede, estimador = NA, tuning = NA, tipos = NULL, n = NA))
  }
  stop("rede precisa ser o resultado de rede_psicometrica() ou um objeto de bootnet::estimateNetwork()",
       if (aceita_matriz) " (ou uma matriz de pesos simétrica)" else "")
}

# Comunidades: vetor na ordem dos nós (ou nomeado pelos nós) ou lista nomeada
# comunidade -> nós. Devolve um vetor de texto na ordem dos nós.
.ee_redes_comunidades <- function(comunidades, nos) {
  if (is.list(comunidades)) {
    v <- rep(NA_character_, length(nos)); names(v) <- nos
    for (nm in names(comunidades)) {
      fora <- setdiff(comunidades[[nm]], nos)
      if (length(fora)) stop("nós ausentes da rede nas comunidades: ", paste(fora, collapse = ", "))
      v[comunidades[[nm]]] <- nm
    }
  } else {
    v <- as.character(comunidades)
    if (!is.null(names(comunidades)) && all(nos %in% names(comunidades))) v <- as.character(comunidades[nos])
    if (length(v) != length(nos)) stop("comunidades precisa ter um valor por nó (", length(nos), " nós)")
    names(v) <- nos
  }
  if (anyNA(v)) stop("há nós sem comunidade: ", paste(nos[is.na(v)], collapse = ", "))
  if (length(unique(v)) < 2) stop("ponte exige pelo menos 2 comunidades")
  v
}

.ee_redes_frase_alcance <- function(idioma) {
  if (idioma == "en") "The results describe conditional associations in cross-sectional data and do not support causal or intervention inferences."
  else "Os resultados descrevem associações condicionais em dados transversais e não autorizam inferências causais ou de intervenção."
}

.ee_redes_leitura_cs <- function(cs, caso_max, idioma = "pt") {
  vapply(cs, function(x) {
    if (is.na(x)) return(if (idioma == "en") "not computable" else "não calculável")
    teto <- abs(x - caso_max) < .005
    if (idioma == "en") {
      l <- if (x < .25) "below the minimum of .25: do not interpret the ordering of this index" else
        if (x <= .50) "at or above the minimum of .25 but not above the preferred .50: interpret with caution" else
          "above the preferred .50"
      if (teto) l <- paste0(l, " (ceiling of the tested grid)")
      return(l)
    }
    l <- if (x < .25) "abaixo do mínimo de .25: não interpretar a ordem deste índice" else
      if (x <= .50) "no mínimo de .25 ou acima, mas não acima do preferível de .50: interpretar com cautela" else
        "acima do preferível de .50"
    if (teto) l <- paste0(l, " (teto da grade de descarte testada)")
    l
  }, character(1), USE.NAMES = FALSE)
}

# ---------------------------------------------------------------- 1. estimação
rede_psicometrica <- function(dados, itens, estimador = "EBICglasso", tuning = 0.5, cor = "cor_auto",
                              ausentes = "pairwise", tamanho_amostra = "maximum", tipos = NULL,
                              redundancia = TRUE, n_fortes = 10, decimal = ".", idioma = "pt", ...) {
  .pr_exige("bootnet", "qgraph")
  estimador <- match.arg(estimador, c("EBICglasso", "IsingFit", "mgm"))
  ausentes <- match.arg(ausentes, c("pairwise", "listwise"))
  tamanho_amostra <- match.arg(tamanho_amostra, c("maximum", "pairwise_average", "minimum",
                                                  "pairwise_maximum", "pairwise_minimum"))
  if (!cor %in% c("cor_auto", "cor", "spearman", "npn", "cov", "cor_mantar"))
    stop("cor precisa ser um de: cor_auto, cor, spearman, npn, cov, cor_mantar")
  if (!is.data.frame(dados)) stop("dados precisa ser um data.frame")
  if (length(itens) < 3) stop("a rede precisa de pelo menos 3 itens")
  falta <- setdiff(itens, names(dados))
  if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "))
  if (!is.numeric(tuning) || length(tuning) != 1 || tuning < 0 || tuning > 1) stop("tuning (γ do EBIC) precisa estar entre 0 e 1")

  d <- as.data.frame(dados)[, itens, drop = FALSE]
  for (v in itens) if (!is.numeric(d[[v]])) stop("o item ", v, " precisa ser numérico (fator deve ser recodificado antes)")
  decisoes <- character(); avisos <- character()

  # Ausentes: o EBICglasso aceita par a par (padrão do JASP); IsingFit e mgm só listwise
  # (no JASP também, jaspNetwork .networkAnalysisMakeBootnetArgs).
  listwise <- ausentes == "listwise" || estimador != "EBICglasso"
  if (listwise) {
    prep <- .ee_dados(d, itens, numericas = itens)
    d <- prep$d; n_excluidos <- prep$n_excluidos
    if (estimador != "EBICglasso" && ausentes == "pairwise")
      decisoes <- c(decisoes, paste0(estimador, " não aceita exclusão par a par: casos incompletos saíram por listwise (", n_excluidos, " casos)."))
    n_usado <- nrow(d); n_incompletos <- 0L
  } else {
    vazias <- rowSums(!is.na(d)) == 0
    d <- d[!vazias, , drop = FALSE]; n_excluidos <- sum(vazias)
    n_incompletos <- sum(!stats::complete.cases(d))
    n_usado <- .ee_redes_n_amostra(d, tamanho_amostra)
  }
  dp <- vapply(d, stats::sd, numeric(1), na.rm = TRUE)
  if (any(is.na(dp) | dp == 0)) stop("itens sem variação: ", paste(itens[is.na(dp) | dp == 0], collapse = ", "))

  niveis <- NULL; recodificados <- FALSE
  if (estimador == "IsingFit") {
    ncat <- vapply(d, function(x) length(unique(x)), numeric(1))
    if (any(ncat != 2)) stop("o modelo de Ising exige itens binários; com mais categorias use EBICglasso, com dados mistos use mgm")
    if (!all(unlist(d) %in% c(0, 1))) {
      d[] <- lapply(d, function(x) as.integer(x == max(x)))
      recodificados <- TRUE
    }
  }
  if (estimador == "mgm") {
    if (is.null(tipos)) {
      tipos <- vapply(d, function(x) if (length(unique(x)) <= 2) "c" else "g", character(1))
      decisoes <- c(decisoes, "mgm: tipos inferidos (itens com 2 valores = categóricos \"c\", demais = gaussianos \"g\"); passar tipos para declarar contagens (\"p\") ou categóricas com mais níveis.")
    }
    tipos <- unname(as.character(tipos))
    if (length(tipos) != length(itens) || !all(tipos %in% c("g", "c", "p"))) stop("tipos precisa ter um valor \"g\", \"c\" ou \"p\" por item")
    niveis <- ifelse(tipos == "c", vapply(d, function(x) length(unique(x)), numeric(1)), 1)
  } else if (estimador == "IsingFit") {
    tipos <- rep("c", length(itens))
  } else {
    tipos <- rep("g", length(itens))
  }

  args <- list(data = d, default = estimador, tuning = tuning, verbose = FALSE)
  if (estimador == "EBICglasso") args <- c(args, list(corMethod = cor, missing = ausentes, sampleSize = tamanho_amostra))
  if (estimador == "IsingFit") args <- c(args, list(rule = "AND", split = "median"))
  if (estimador == "mgm") args <- c(args, list(type = tipos, level = niveis, rule = "AND", criterion = "EBIC"))
  extras <- list(...)
  args[names(extras)] <- extras
  cap <- .ee_redes_capturar(do.call(bootnet::estimateNetwork, args))
  modelo <- cap$valor
  if (length(cap$avisos)) avisos <- c(avisos, .ee_redes_traduz_aviso(cap$avisos))

  W <- modelo$graph
  dimnames(W) <- list(itens, itens)
  p <- length(itens); possiveis <- p * (p - 1) / 2
  ar <- .ee_redes_arestas(W)
  nz <- sum(ar$peso != 0); pos <- sum(ar$peso > 0); neg <- sum(ar$peso < 0)
  densidade <- nz / possiveis
  cent <- .ee_redes_centralidade(W)
  fortes <- ar[ar$peso != 0, , drop = FALSE]
  fortes <- fortes[order(-abs(fortes$peso)), , drop = FALSE]
  fortes <- utils::head(fortes, n_fortes); rownames(fortes) <- NULL
  ggm <- estimador == "EBICglasso"
  zero_peso <- !ggm  # peso de GGM é correlação parcial (sem zero à esquerda); Ising e mgm não são

  # Redundância antes da rede de itens: UVA do EGAnet.
  redund <- NULL
  if (isTRUE(redundancia)) {
    if (requireNamespace("EGAnet", quietly = TRUE)) {
      u <- tryCatch(suppressWarnings(suppressMessages(EGAnet::UVA(data = d, reduce = FALSE, verbose = FALSE))),
                    error = function(e) NULL)
      if (!is.null(u) && !is.null(u$wto$pairwise)) {
        pw <- as.data.frame(u$wto$pairwise)
        redund <- pw[pw$wto >= .25, c("node_i", "node_j", "wto"), drop = FALSE]
        rownames(redund) <- NULL
        if (nrow(redund)) {
          avisos <- c(avisos, paste0("UVA (EGAnet) sinalizou ", nrow(redund), " par(es) com sobreposição de vizinhança (wTO ≥ .25): ",
                                     paste0(.ee_redes_aresta(redund$node_i, redund$node_j), " (wTO = ",
                                            vapply(redund$wto, .pr_num, "", casas = 2, zero = FALSE), ")", collapse = "; "),
                                     ". Aresta forte entre itens redundantes pode refletir redação parecida e não relação substantiva; decidir antes da rede se mantém, combina ou remove, e declarar."))
        } else decisoes <- c(decisoes, "UVA (EGAnet): nenhum par com wTO ≥ .25 (sem sinal de itens redundantes).")
      } else decisoes <- c(decisoes, "UVA (EGAnet) não rodou nestes dados: conferir redundância de itens à parte antes de interpretar arestas fortes entre itens vizinhos.")
    } else decisoes <- c(decisoes, "EGAnet ausente: rodar EGAnet::UVA() antes da rede de itens para checar itens redundantes.")
  }

  # ------------------------------------------------ decisões (com fonte)
  nome_est <- switch(estimador, EBICglasso = "EBICglasso (GGM: correlações parciais regularizadas por LASSO, modelo escolhido pelo EBIC)",
                     IsingFit = "IsingFit (modelo de Ising para itens binários: regressões logísticas com LASSO por nó, regra AND)",
                     mgm = "mgm (modelo gráfico misto: regressões nodais com LASSO, regra AND)")
  decisoes <- c(decisoes,
    paste0("Estimador: ", nome_est, ". Padrão do livro e do JASP para itens ordinais e contínuos; Ising para binários e mgm para dados mistos seguem a escolha por tipo de dado."),
    paste0("γ (tuning do EBIC) = ", .pr_num(tuning, 2, TRUE, decimal), " declarado: mais alto deixa a rede mais esparsa. ",
           if (estimador != "EBICglasso") "No IsingFit e no mgm o bootnet usa 0,25 por padrão; 0,5 é o valor do JASP." else "Coincide com o padrão do bootnet, mas o relato precisa do valor."),
    if (ggm) paste0("Correlação: ", cor, if (cor == "cor_auto") " (itens com até 7 categorias viram policóricas; acima disso, Pearson). Declarado porque o padrão do bootnet_EBICglasso é \"cor\" (Pearson)." else "."),
    if (ggm) paste0("Ausentes: ", ausentes, "; tamanho de amostra: ", tamanho_amostra, " (N = ", round(n_usado, 1), "), os padrões do JASP. O bootnet usaria \"pairwise_average\"; com dado ausente, as duas escolhas podem dar redes diferentes."),
    "Centralidade como no JASP: qgraph::centrality com pesos e sinal (weighted, signed), sem todos os caminhos mínimos; intermediação dividida por 2 (rede não direcionada); escores z = (x - média) / DP amostral, o \"Normalized\" do JASP. Tabela bruta em $centralidade_bruta.",
    "Força é o índice principal; influência esperada quando arestas negativas representam oposição. Proximidade e intermediação não medem importância do nó (Bringmann et al., 2019) e só entram com CS adequado.",
    "Ordem de centralidade só se interpreta depois de estabilidade_rede() (CS por índice). Centralidade alta não é alvo de intervenção comprovado.",
    paste0("Tamanho da amostra: N = ", round(n_usado, 1), " para ", p, " nós e ", possiveis, " arestas possíveis, p(p-1)/2. Não há número mágico de participantes; para planejar, simular com bootnet::netSimulator ou powerly (Constantin et al., 2026)."),
    "Posição dos nós na figura não é centralidade: o layout é algorítmico.",
    "Ponte e previsibilidade, quando a pergunta pedir: ponte(res, comunidades) e previsibilidade(res).")
  decisoes <- decisoes[nzchar(decisoes)]
  if (recodificados) decisoes <- c(decisoes, "Itens binários recodificados para 0/1 (maior valor = 1) antes do IsingFit.")
  if (!ggm) decisoes <- c(decisoes, "Pesos do Ising e do mgm não são correlações parciais: saem com zero à esquerda e não se comparam com pesos de GGM.")
  if (!listwise && n_incompletos > 0)
    avisos <- c(avisos, paste0(n_incompletos, " casos incompletos entraram por exclusão par a par; o N do EBIC é o \"", tamanho_amostra, "\" (", round(n_usado, 1), ")."))
  avisos <- c(avisos, "Interpretar a ordem de centralidade só depois de estabilidade_rede(); a frase APA descreve a força sem afirmar estabilidade.")

  # ------------------------------------------------ tabelas
  tabela <- data.frame(
    indicador = c("nós", "arestas possíveis", "arestas não nulas", "arestas positivas", "arestas negativas",
                  "densidade", "esparsidade", "N usado na estimação", "casos incompletos (par a par)", "casos excluídos",
                  "γ (tuning)"),
    valor = c(p, possiveis, nz, pos, neg, densidade, 1 - densidade, n_usado, n_incompletos, n_excluidos, tuning),
    stringsAsFactors = FALSE)
  efeito <- if (nrow(fortes)) data.frame(indice = c(paste0("peso da aresta mais forte (", fortes$aresta[1], ")"), "densidade"),
                                         valor = c(fortes$peso[1], densidade), ic_inf = NA_real_, ic_sup = NA_real_,
                                         principal = c(TRUE, FALSE), stringsAsFactors = FALSE) else NULL
  checklist <- data.frame(
    elemento = c("método de estimação nomeado (modelo, estimador, γ, correlação, ausentes, N)",
                 "número exato de arestas não nulas sobre as possíveis",
                 "aresta mais forte com IC bootstrap",
                 "arestas menores indistinguíveis: IC sobrepostos ou teste de diferença",
                 "CS-coefficient por índice de centralidade",
                 "tamanho da amostra em relação ao número de nós",
                 "frase de alcance causal"),
    situacao = c("na frase APA", paste0("na frase APA (", nz, " de ", possiveis, ")"),
                 "pendente: estabilidade_rede()", "pendente: estabilidade_rede()", "pendente: estabilidade_rede()",
                 paste0("N = ", round(n_usado), "; ", p, " nós; ", possiveis, " arestas possíveis"), "na frase APA"),
    stringsAsFactors = FALSE)

  # ------------------------------------------------ frase APA
  en <- idioma == "en"
  zp <- function(x) .pr_num(x, 2, zero_peso, decimal)
  gtxt <- .pr_num(tuning, 2, TRUE, decimal)
  ntxt <- format(round(n_usado), big.mark = "")
  modelo_txt <- switch(estimador,
    EBICglasso = if (en) paste0("a Gaussian graphical model (edges = partial correlations) with LASSO regularization and EBIC model selection (EBICglasso; γ = ", gtxt, "), based on ", .ee_redes_cor_txt(d, cor, idioma))
                 else paste0("modelo gráfico gaussiano (arestas = correlações parciais), com regularização LASSO e seleção do modelo pelo EBIC (EBICglasso; γ = ", gtxt, "), a partir de ", .ee_redes_cor_txt(d, cor, idioma)),
    IsingFit = if (en) paste0("an Ising model (IsingFit: nodewise LASSO logistic regressions, EBIC with γ = ", gtxt, ", AND rule)")
               else paste0("modelo de Ising (IsingFit: regressões logísticas com LASSO por nó, EBIC com γ = ", gtxt, ", regra AND)"),
    mgm = if (en) paste0("a mixed graphical model (mgm: nodewise LASSO regressions, EBIC with γ = ", gtxt, ", AND rule)")
          else paste0("modelo gráfico misto (mgm: regressões nodais com LASSO, EBIC com γ = ", gtxt, ", regra AND)"))
  aus_txt <- if (listwise) { if (en) "listwise deletion of incomplete cases" else "exclusão listwise de casos incompletos" } else
    if (tamanho_amostra == "maximum") { if (en) "pairwise exclusion of missing data and maximum sample size" else "exclusão par a par de ausentes e tamanho de amostra máximo" } else
      { if (en) paste0("pairwise exclusion of missing data (sample size: ", tamanho_amostra, ")") else paste0("exclusão par a par de ausentes (tamanho de amostra: ", tamanho_amostra, ")") }
  f1 <- if (en) paste0("The network of the ", p, " items was estimated as ", modelo_txt, ", with ", aus_txt, " (N = ", ntxt, ").")
        else paste0("A rede dos ", p, " itens foi estimada como ", modelo_txt, ", com ", aus_txt, " (N = ", ntxt, ").")
  dens <- .pr_num(densidade, 2, FALSE, decimal)
  f2 <- if (en) paste0("Of the ", possiveis, " possible edges, ", nz, " were estimated as non-zero (",
                       .ee_junta(paste0("density = ", dens), paste0(pos, " positive and ", neg, " negative"), decimal = decimal), ").")
        else paste0("Das ", possiveis, " arestas possíveis, ", nz, " foram estimadas como não nulas (",
                    .ee_junta(paste0("densidade = ", dens), paste0(pos, " positivas e ", neg, " negativas"), decimal = decimal), ").")
  f3 <- character()
  if (nrow(fortes) >= 1) {
    f3 <- if (en) paste0("The strongest edge connected ", fortes$no1[1], " and ", fortes$no2[1], " (weight = ", zp(fortes$peso[1]), ")")
          else paste0("A aresta de maior peso ligou ", fortes$no1[1], " e ", fortes$no2[1], " (peso = ", zp(fortes$peso[1]), ")")
    if (nrow(fortes) >= 2) f3 <- paste0(f3, if (en) ", followed by " else ", seguida de ", fortes$aresta[2],
                                        if (en) " (weight = " else " (peso = ", zp(fortes$peso[2]), ")")
    f3 <- paste0(f3, ".")
  }
  zf <- cent$z$forca; names(zf) <- cent$z$no
  ordem <- order(-zf)
  topo <- utils::head(ordem, min(3, p - 1))
  zt <- function(i) paste0(names(zf)[i], " (", .pr_num(zf[i], 2, TRUE, decimal), ")")
  f4 <- if (en) paste0("In strength (z scores), the highest values were those of ", paste(vapply(topo, zt, ""), collapse = ", "),
                       ", and the lowest was that of ", zt(utils::tail(ordem, 1)), ".")
        else paste0("Na força (escores z), os maiores valores foram de ", paste(vapply(topo, zt, ""), collapse = ", "),
                    ", e o menor foi de ", zt(utils::tail(ordem, 1)), ".")
  apa <- paste(c(f1, f2, f3, f4, .ee_redes_frase_alcance(idioma)), collapse = " ")

  rotulo_z <- "Centralidade em escores z (o \"Normalized\" do JASP)"
  .ee_resultado(paste0("Rede psicométrica (", estimador, ")"), n = round(n_usado, 1), n_excluidos = n_excluidos,
                tabela = tabela, efeito = efeito, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = c(list(tabelas_extra = stats::setNames(list(fortes, cent$z, checklist),
                                                               c(paste0("Arestas mais fortes (", nrow(fortes), " de ", nz, " não nulas)"),
                                                                 rotulo_z, "Checklist de relato (livro e 520)")),
                               modelo = modelo, pesos = W, arestas = ar, centralidade_bruta = cent$bruta,
                               centralidade_z = cent$z, redundancia = redund, checklist = checklist,
                               estimador = estimador, tuning = tuning, tipos = tipos,
                               configuracao = list(estimador = estimador, tuning = tuning, cor = if (ggm) cor else NA,
                                                   ausentes = if (listwise) "listwise" else ausentes,
                                                   tamanho_amostra = if (ggm && !listwise) tamanho_amostra else NA,
                                                   n_usado = n_usado))))
}

# ---------------------------------------------------------------- 2. precisão e estabilidade
estabilidade_rede <- function(rede, n_boot = 1000, indices = c("strength", "expectedInfluence"), semente = 2026,
                              n_boot_casos = n_boot, caso_min = 0.05, caso_max = 0.75, caso_n = 10,
                              executar = c("arestas", "casos"), comunidades = NULL, n_fortes = 10,
                              n_nucleos = 1, decimal = ".", idioma = "pt") {
  .pr_exige("bootnet")
  obj <- .ee_redes_obter(rede)
  validos <- names(.ee_redes_rotulos_indice)
  if (!length(indices) || !all(indices %in% validos)) stop("indices aceitos: ", paste(validos, collapse = ", "))
  executar <- match.arg(executar, several.ok = TRUE)
  ponte_ind <- grepl("^bridge", indices)
  if (any(ponte_ind) && is.null(comunidades)) stop("índices de ponte exigem o argumento comunidades")
  if (caso_min <= 0 || caso_max >= 1 || caso_min >= caso_max) stop("caso_min e caso_max precisam satisfazer 0 < caso_min < caso_max < 1")
  nos <- colnames(obj$W)
  en <- idioma == "en"
  ggm <- identical(obj$estimador, "EBICglasso")
  zp <- function(x) .pr_num(x, 2, !ggm, decimal)
  base <- list(data = obj$modelo, nCores = n_nucleos, verbose = FALSE)
  if (!is.null(comunidades)) base$communities <- unname(.ee_redes_comunidades(comunidades, nos))
  decisoes <- character(); avisos <- character(); apa <- character()
  tab_cs <- NULL; ar <- NULL; comp <- NULL; b_np <- NULL; b_c <- NULL; niveis <- NULL; cs <- NULL
  avisos_boot <- character()

  if ("arestas" %in% executar) {
    set.seed(semente)
    cap <- .ee_redes_capturar(do.call(bootnet::bootnet, c(base, list(nBoots = n_boot, type = "nonparametric",
                                                                       statistics = unique(c("edge", indices))))))
    b_np <- cap$valor
    avisos_boot <- c(avisos_boot, cap$avisos)
    if (length(cap$avisos)) avisos <- c(avisos, paste0("Bootstrap das arestas: ", length(cap$avisos), " tipo(s) de aviso do bootnet/qgraph nas reamostras (ver $avisos_bootstrap)."))
    s <- as.data.frame(summary(b_np, statistics = "edge"))
    ar <- data.frame(aresta = .ee_redes_aresta(s$node1, s$node2), id = s$id, no1 = s$node1, no2 = s$node2,
                     peso = s$sample, media_boot = s$mean, ic_inf = s$q2.5, ic_sup = s$q97.5,
                     prop_nao_nula = 1 - s$prop0, stringsAsFactors = FALSE)
    ar <- ar[order(-abs(ar$peso)), , drop = FALSE]; rownames(ar) <- NULL
    nz <- ar[ar$peso != 0, , drop = FALSE]
    # IC em valor absoluto, para comparar força de arestas de sinais diferentes.
    abs_ic <- function(lo, hi) cbind(ifelse(lo >= 0, lo, ifelse(hi <= 0, -hi, 0)),
                                     ifelse(lo >= 0, hi, ifelse(hi <= 0, -lo, pmax(-lo, hi))))
    if (nrow(nz) >= 2) {
      bt <- as.data.frame(b_np$bootTable); bt <- bt[bt$type == "edge", c("name", "id", "value")]
      v1 <- bt$value[bt$id == nz$id[1]]; names(v1) <- bt$name[bt$id == nz$id[1]]
      v2 <- bt$value[bt$id == nz$id[2]]; names(v2) <- bt$name[bt$id == nz$id[2]]
      dif <- abs(v1) - abs(v2[names(v1)])
      ic_dif <- stats::quantile(dif, c(.025, .975), type = 6, na.rm = TRUE, names = FALSE)
      dist <- ic_dif[1] > 0 || ic_dif[2] < 0
      a <- abs_ic(nz$ic_inf, nz$ic_sup)
      sobrep <- sum(a[-1, 2] >= a[1, 1] & a[-1, 1] <= a[1, 2])
      comp <- data.frame(aresta_1 = nz$aresta[1], aresta_2 = nz$aresta[2], peso_1 = nz$peso[1], peso_2 = nz$peso[2],
                         dif_abs = abs(nz$peso[1]) - abs(nz$peso[2]), ic_inf = ic_dif[1], ic_sup = ic_dif[2],
                         distinguiveis = dist, n_ic_sobrepostos = sobrep, n_outras_nao_nulas = nrow(nz) - 1,
                         stringsAsFactors = FALSE)
    }
    decisoes <- c(decisoes,
      paste0("Precisão das arestas: bootstrap não paramétrico (reamostragem de participantes com reposição), ", n_boot,
             " reamostras (o livro usa 1000), semente ", semente, " fixada antes do bootstrap; IC 95% por percentis (quantil tipo 6, como o summary do bootnet)."),
      "IC de aresta em rede regularizada é descritivo da variabilidade, não teste: o LASSO encolhe pesos para zero, e incluir o zero no IC não prova que a aresta é nula.",
      "Leitura: arestas com IC sobrepostos não podem ser declaradas diferentes em força. A comparação formal é o teste de diferença do bootnet: a função compara a aresta mais forte com a segunda (|peso|, IC 95% da diferença nas reamostras); para outros pares, bootnet::differenceTest(res$boot_arestas, x = \"A\", x2 = \"B\", y = \"C\", y2 = \"D\", measure = \"edge\"), sem correção para comparações múltiplas.")
    if (n_boot < 1000) avisos <- c(avisos, paste0("Bootstrap de arestas com ", n_boot, " reamostras: serve para conferência rápida; para relato, 1000 ou mais."))
    if (!ggm) avisos <- c(avisos, "Bootstrap não paramétrico em Ising ou mgm: o bootnet roda, mas a literatura de redes binárias discute outras opções; conferir antes de relatar.")
    if (!is.null(comp)) {
      ic1 <- .ee_ic(nz$ic_inf[1], nz$ic_sup[1], 2, !ggm, decimal, .95, idioma)
      ic2 <- .ee_ic(nz$ic_inf[2], nz$ic_sup[2], 2, !ggm, decimal, .95, idioma)
      icd <- .ee_ic(comp$ic_inf, comp$ic_sup, 2, !ggm, decimal, .95, idioma)
      a1 <- paste0(nz$aresta[1], " (", .ee_junta(paste0(if (en) "weight = " else "peso = ", zp(nz$peso[1])), ic1, decimal = decimal), ")")
      a2 <- paste0(nz$aresta[2], " (", .ee_junta(paste0(if (en) "weight = " else "peso = ", zp(nz$peso[2])), ic2, decimal = decimal), ")")
      dtxt <- paste0(if (en) "difference in |weight| = " else "diferença de |peso| = ", zp(comp$dif_abs), .ee_sep(decimal), icd)
      apa <- c(apa, if (en)
        paste0("Edge accuracy was assessed with a nonparametric bootstrap (", n_boot, " resamples, seed ", semente,
               ", percentile 95% CIs), and the CIs were read as descriptive of variability rather than as tests of whether an edge is zero. ",
               "The strongest edge, ", a1, if (dist) ", was distinguishable from the second strongest, " else ", was not distinguishable from the second strongest, ",
               a2, "; ", dtxt, ". ", if (sobrep == 0) "None" else sobrep, " of the other ", nrow(nz) - 1, " non-zero edges had CIs overlapping that of the strongest edge.")
      else
        paste0("A precisão das arestas foi avaliada por bootstrap não paramétrico (", n_boot, " reamostras, semente ", semente,
               ", IC 95% por percentis), e os IC foram lidos como descritivos da variabilidade, não como teste de que a aresta é nula. ",
               "A aresta mais forte, ", a1, if (dist) ", distinguiu-se da segunda maior, " else ", não se distinguiu da segunda maior, ",
               a2, "; ", dtxt, ". ", if (sobrep == 0) "Nenhuma" else sobrep, " das outras ", nrow(nz) - 1,
               if (sobrep <= 1) " arestas não nulas teve IC sobreposto ao da aresta mais forte." else " arestas não nulas tiveram IC sobreposto ao da aresta mais forte."))
    }
  }

  if ("casos" %in% executar) {
    set.seed(semente)
    cap <- .ee_redes_capturar(do.call(bootnet::bootnet, c(base, list(nBoots = n_boot_casos, type = "case", statistics = indices,
                                                                       caseMin = caso_min, caseMax = caso_max, caseN = caso_n))))
    b_c <- cap$valor
    avisos_boot <- c(avisos_boot, cap$avisos)
    if (length(cap$avisos)) avisos <- c(avisos, paste0("Bootstrap de casos: ", length(cap$avisos), " tipo(s) de aviso do bootnet/qgraph nas reamostras (ver $avisos_bootstrap)."))
    cs_all <- suppressMessages(bootnet::corStability(b_c, cor = 0.7, verbose = FALSE))
    cs <- unname(cs_all[indices]); names(cs) <- indices
    niveis <- sort(unique(1 - b_c$bootTable$nPerson / b_c$sample$nPerson))
    tab_cs <- data.frame(indice = .ee_redes_rotulo(indices, idioma), indice_bootnet = indices, cs = cs,
                         leitura = .ee_redes_leitura_cs(cs, caso_max, idioma), stringsAsFactors = FALSE)
    rownames(tab_cs) <- NULL
    decisoes <- c(decisoes,
      paste0("Estabilidade da centralidade: bootstrap de descarte de casos (case-dropping), ", n_boot_casos, " reamostras, descarte de ",
             round(100 * caso_min), "% a ", round(100 * caso_max), "% dos casos em ", caso_n, " níveis (caseMin = ", caso_min,
             ", caseMax = ", caso_max, ", caseN = ", caso_n, "; padrões do bootnet, declarados porque o livro só ilustra 10%, 20%...); semente ", semente, "."),
      "CS-coefficient: maior proporção de casos que se pode descartar mantendo, com 95% de probabilidade, correlação de pelo menos .70 com a centralidade da amostra completa.",
      "Cortes do CS: não abaixo de .25 e de preferência acima de .50, convenções heurísticas. Sem patamar de .70: .70 é a correlação da definição, não um corte.",
      paste0("O CS só assume valores da grade testada (", paste(vapply(niveis, .pr_num, "", casas = 2, zero = FALSE), collapse = ", "), "); o teto é caseMax."))
    if (n_boot_casos < 1000) avisos <- c(avisos, paste0("Bootstrap de casos com ", n_boot_casos, " reamostras: o CS fica instável com poucas reamostras por nível; para relato, 1000 ou mais."))
    if (any(indices %in% c("closeness", "betweenness")))
      decisoes <- c(decisoes, "Proximidade e intermediação entraram porque foram pedidas; mesmo com CS alto, não medem importância do nó em redes psicológicas (Bringmann et al., 2019; Hallquist et al., 2021).")
    baixos <- indices[!is.na(cs) & cs < .25]
    if (length(baixos)) avisos <- c(avisos, paste0("CS abaixo de .25 em ", paste(.ee_redes_rotulo(baixos, "pt"), collapse = ", "),
                                                   ": não interpretar a ordem desse índice. CS baixo num índice não invalida a rede."))
    cs_txt <- vapply(seq_along(indices), function(i) {
      x <- cs[i]; rot <- .ee_redes_rotulo(indices[i], idioma)
      leitura <- if (is.na(x)) "" else if (x < .25) { if (en) ", below the minimum of .25" else ", abaixo do mínimo de .25" } else
        if (x <= .50) { if (en) ", above the minimum of .25 but not above the preferred .50" else ", acima do mínimo de .25, mas não acima do preferível de .50" } else
          { if (en) ", above the preferred .50" else ", acima do preferível de .50" }
      paste0("CS(", rot, ") = ", .pr_num(x, 2, FALSE, decimal), leitura)
    }, character(1))
    apa <- c(apa, if (en)
      paste0("Centrality stability was assessed with a case-dropping bootstrap (", n_boot_casos, " resamples, dropping ",
             round(100 * caso_min), "% to ", round(100 * caso_max), "% of cases in ", caso_n, " steps): ",
             paste(cs_txt, collapse = "; "), ". The CS-coefficient is the maximum proportion of cases that can be dropped while retaining, with 95% probability, a correlation of at least .70 with the full-sample centrality (Epskamp et al., 2018).")
    else
      paste0("A estabilidade da centralidade foi avaliada por bootstrap de descarte de casos (", n_boot_casos, " reamostras, descarte de ",
             round(100 * caso_min), "% a ", round(100 * caso_max), "% dos casos em ", caso_n, " níveis): ",
             paste(cs_txt, collapse = "; "), ". O CS é a maior proporção de casos que se pode descartar mantendo, com 95% de probabilidade, correlação de pelo menos .70 com a centralidade da amostra completa (Epskamp et al., 2018)."))
  }

  checklist <- data.frame(
    elemento = c("aresta mais forte com IC bootstrap",
                 "arestas menores indistinguíveis: IC sobrepostos ou teste de diferença",
                 "CS-coefficient por índice de centralidade"),
    situacao = c(if (!is.null(comp)) "na frase APA" else if (is.null(ar)) "não rodado (executar sem \"arestas\")" else "menos de 2 arestas não nulas",
                 if (!is.null(comp)) "na frase APA (sobreposição e diferença com a segunda maior)" else "não disponível",
                 if (!is.null(tab_cs)) "na frase APA" else "não rodado (executar sem \"casos\")"),
    stringsAsFactors = FALSE)
  efeito <- if (!is.null(tab_cs)) data.frame(indice = paste0("CS (", tab_cs$indice, ")"), valor = tab_cs$cs, ic_inf = NA_real_, ic_sup = NA_real_,
                                             principal = tab_cs$indice_bootnet == (if ("strength" %in% indices) "strength" else indices[1]),
                                             stringsAsFactors = FALSE) else NULL
  extras_tab <- list()
  if (!is.null(ar)) extras_tab[[paste0("IC bootstrap das arestas mais fortes (", min(n_fortes, sum(ar$peso != 0)), " de ", sum(ar$peso != 0), " não nulas)")]] <-
    utils::head(ar[ar$peso != 0, c("aresta", "peso", "media_boot", "ic_inf", "ic_sup", "prop_nao_nula")], n_fortes)
  extras_tab[["Checklist de relato (livro e 520)"]] <- checklist
  .ee_resultado("Precisão e estabilidade da rede (bootnet)", n = obj$n, tabela = tab_cs, efeito = efeito,
                posthoc = comp, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(titulo_posthoc = "Aresta mais forte contra a segunda (diferença de |peso|, IC 95% bootstrap):",
                             tabelas_extra = extras_tab, arestas_ic = ar, comparacao_arestas = comp, cs = cs,
                             niveis_descarte = niveis, boot_arestas = b_np, boot_casos = b_c, checklist = checklist,
                             avisos_bootstrap = unique(avisos_boot),
                             configuracao = list(n_boot = n_boot, n_boot_casos = n_boot_casos, semente = semente,
                                                 caseMin = caso_min, caseMax = caso_max, caseN = caso_n,
                                                 indices = indices, n_nucleos = n_nucleos)))
}

# ---------------------------------------------------------------- 3. ponte
ponte <- function(rede, comunidades, decimal = ".", idioma = "pt") {
  .pr_exige("networktools")
  obj <- .ee_redes_obter(rede, aceita_matriz = TRUE)
  W <- obj$W; nos <- colnames(W)
  com <- .ee_redes_comunidades(comunidades, nos)
  b <- networktools::bridge(W, communities = unname(com), useCommunities = "all", normalize = FALSE)
  tab <- data.frame(no = nos, comunidade = unname(com),
                    forca_ponte = unname(b[["Bridge Strength"]]),
                    influencia_esperada_ponte = unname(b[["Bridge Expected Influence (1-step)"]]),
                    influencia_esperada_ponte_2 = unname(b[["Bridge Expected Influence (2-step)"]]),
                    intermediacao_ponte = unname(b[["Bridge Betweenness"]]),
                    proximidade_ponte = unname(b[["Bridge Closeness"]]),
                    stringsAsFactors = FALSE)
  tab$forca_ponte_z <- .ee_redes_z(tab$forca_ponte)
  en <- idioma == "en"
  o <- order(-tab$forca_ponte); topo <- utils::head(o, min(3, nrow(tab)))
  k <- length(unique(com))
  ftxt <- paste(vapply(topo, function(i) paste0(tab$no[i], " (", .pr_num(tab$forca_ponte[i], 2, TRUE, decimal), ")"), ""), collapse = ", ")
  apa <- if (en) paste0("With ", k, " communities defined a priori, bridge strength (Jones et al., 2021) was highest for ", ftxt, ". ",
                        .ee_redes_frase_alcance(idioma))
         else paste0("Com ", k, " comunidades definidas a priori, a força de ponte (Jones et al., 2021) foi maior em ", ftxt, ". ",
                     .ee_redes_frase_alcance(idioma))
  decisoes <- c(
    "Índices de ponte de Jones, Ma e McNally (2021) pelo networktools::bridge: força de ponte (soma dos |pesos| com nós de outras comunidades), influência esperada de ponte de 1 e 2 passos, intermediação e proximidade de ponte; valores brutos, normalize = FALSE.",
    "Comunidades precisam vir antes do resultado: da teoria, das subescalas do instrumento ou de um EGA declarado; trocar a partição muda todos os índices.",
    "A função não fixa corte para \"nó ponte\": ordena e mostra o z da força de ponte. Conferir no artigo de Jones et al. (2021) antes de adotar um corte.",
    "Estabilidade dos índices de ponte: estabilidade_rede(res, indices = \"bridgeStrength\", comunidades = ...).",
    "Nó ponte numa rede transversal é posição na estrutura estimada, não mecanismo de comorbidade nem alvo de intervenção comprovado.")
  .ee_resultado("Centralidade de ponte (networktools)", n = obj$n, tabela = tab, decisoes = decisoes, apa = apa,
                extra = list(comunidades = com, bridge = b))
}

# ---------------------------------------------------------------- 4. previsibilidade
previsibilidade <- function(rede, tuning = NULL, decimal = ".", idioma = "pt") {
  .pr_exige("mgm")
  obj <- .ee_redes_obter(rede)
  if (is.null(obj$modelo$data)) stop("o objeto da rede não guarda os dados; reestimar com rede_psicometrica()")
  d <- as.data.frame(obj$modelo$data)
  nos <- colnames(obj$W)
  completos <- stats::complete.cases(d)
  n_exc <- sum(!completos)
  X <- as.matrix(d[completos, , drop = FALSE]); colnames(X) <- nos
  est <- obj$estimador
  tipos <- if (identical(est, "mgm") && !is.null(obj$tipos)) obj$tipos else if (identical(est, "IsingFit")) rep("c", ncol(X)) else rep("g", ncol(X))
  niveis <- ifelse(tipos == "c", apply(X, 2, function(x) length(unique(x))), 1)
  gam <- if (is.null(tuning)) { if (is.na(obj$tuning)) 0.25 else obj$tuning } else tuning
  fit <- suppressWarnings(suppressMessages(mgm::mgm(data = X, type = tipos, level = niveis, lambdaSel = "EBIC", lambdaGam = gam,
                                                     ruleReg = "AND", k = 2, pbar = FALSE, signInfo = FALSE)))
  pr <- stats::predict(fit, X, errorCon = "R2", errorCat = c("CC", "nCC"))
  er <- as.data.frame(pr$errors)
  tab <- data.frame(no = nos, tipo = tipos, stringsAsFactors = FALSE)
  tab$r2 <- if ("R2" %in% names(er)) suppressWarnings(as.numeric(er$R2)) else NA_real_
  tab$cc <- if ("CC" %in% names(er)) suppressWarnings(as.numeric(er$CC)) else NA_real_
  tab$ncc <- if ("nCC" %in% names(er)) suppressWarnings(as.numeric(er$nCC)) else NA_real_
  optwi <- obj$modelo$results$optwi
  if (identical(est, "EBICglasso") && is.matrix(optwi) && nrow(optwi) == length(nos))
    tab$r2_ggm_implicito <- 1 - 1 / diag(optwi)
  en <- idioma == "en"
  decisoes <- c(
    paste0("Previsibilidade (Haslbeck & Waldorp, 2018): proporção da variância de cada nó explicada pelos vizinhos, pelo mgm (EBIC com γ = ",
           .pr_num(gam, 2, TRUE, decimal), if (is.null(tuning)) ", o mesmo da rede" else "", ", regra AND), R² para nós gaussianos; para categóricos, acurácia (CC) e acurácia normalizada além da categoria modal (nCC)."),
    "O mgm reestima a rede por regressões nodais; os pesos dele não são os do EBICglasso. Itens ordinais entram como gaussianos (escala bruta padronizada).",
    if (!is.null(tab$r2_ggm_implicito)) "r2_ggm_implicito: R² implicado pela própria rede EBICglasso (1 - 1 / diagonal da matriz de precisão regularizada, na métrica da correlação usada, policórica com cor_auto); coerente com a figura, mas não é o R² de predição do mgm." else "",
    "Previsibilidade alta não é alvo de intervenção: indica quanto do nó a rede explica, numa amostra transversal.")
  decisoes <- decisoes[nzchar(decisoes)]
  avisos <- if (n_exc > 0) paste0(n_exc, " casos incompletos ficaram fora (o mgm exige casos completos).") else character()
  apa <- character()
  if (any(!is.na(tab$r2))) {
    g <- !is.na(tab$r2); i_min <- which(g)[which.min(tab$r2[g])]; i_max <- which(g)[which.max(tab$r2[g])]
    r <- function(x) .pr_num(x, 3, FALSE, decimal)
    apa <- if (en) paste0("Node predictability (Haslbeck & Waldorp, 2018), estimated with mgm, had a mean R² = ", r(mean(tab$r2[g])),
                          " (from ", r(tab$r2[i_min]), " for ", tab$no[i_min], " to ", r(tab$r2[i_max]), " for ", tab$no[i_max], ").")
           else paste0("A previsibilidade dos nós (Haslbeck & Waldorp, 2018), estimada com mgm, teve média de R² = ", r(mean(tab$r2[g])),
                       " (de ", r(tab$r2[i_min]), " em ", tab$no[i_min], " a ", r(tab$r2[i_max]), " em ", tab$no[i_max], ").")
  }
  if (any(!is.na(tab$ncc))) {
    g <- !is.na(tab$ncc)
    apa <- c(apa, if (en) paste0("For categorical nodes, mean normalized accuracy (nCC) = ", .pr_num(mean(tab$ncc[g]), 2, FALSE, decimal), ".")
                  else paste0("Nos nós categóricos, a acurácia normalizada média (nCC) foi ", .pr_num(mean(tab$ncc[g]), 2, FALSE, decimal), "."))
  }
  apa <- paste(apa, collapse = " ")
  .ee_resultado("Previsibilidade dos nós (mgm)", n = nrow(X), n_excluidos = n_exc, tabela = tab, decisoes = decisoes,
                avisos = avisos, apa = apa, extra = list(modelo_mgm = fit, tuning = gam))
}

# ---------------------------------------------------------------- v1.5: modelos gráficos probabilísticos
# rede_bayesiana(): rede bayesiana (grafo acíclico dirigido) aprendida dos dados com o
# bnlearn, a rota do módulo "Probabilistic Graphical Models" do curso (hc, tabu, pc.stable,
# mmhc). Para não ler o acaso de uma reamostra como estrutura, a skill estima a rede média
# de `reamostras` reamostras bootstrap (boot.strength; Friedman et al., 1999; Scutari &
# Nagarajan, 2013): cada aresta sai com a força (proporção de redes em que aparece) e a
# direção (proporção com aquele sentido, dado que aparece). O limiar de inclusão é o
# estimado pelo bnlearn (Scutari & Nagarajan, 2013), salvo outro informado.
# Variáveis com até 7 valores viram categorias (rede discreta, como no código das aulas);
# as demais, contínuas (rede gaussiana). Misturar os dois tipos pede discretizar.
# Direção de aresta em rede aprendida de dados transversais não é causalidade: arestas de
# uma mesma classe de equivalência têm sentido indeterminado (direção perto de 0,5).
rede_bayesiana <- function(dados, vars, algoritmo = c("hc", "tabu", "pc.stable", "mmhc"), reamostras = 500, limiar = NULL,
                           tipo = c("auto", "discreta", "gaussiana"), comparar_algoritmos = FALSE, semente = 2026, decimal = ".", idioma = "pt") {
  .pr_exige("bnlearn")
  algoritmo <- match.arg(algoritmo); tipo <- match.arg(tipo)
  dd <- .ee_dados(dados, vars); d <- dd$d
  if (tipo == "auto") tipo <- if (all(vapply(d, function(x) length(unique(x)) <= 7, TRUE))) "discreta" else
    if (all(vapply(d, function(x) is.numeric(x) && length(unique(x)) > 7, TRUE))) "gaussiana" else
      stop("variáveis com até 7 valores misturadas a contínuas: discretizar as contínuas ou informar tipo")
  d <- as.data.frame(lapply(d, function(x) if (tipo == "discreta") factor(x) else as.numeric(x)))
  names(d) <- make.names(vars)
  set.seed(semente)
  bs <- bnlearn::boot.strength(d, R = reamostras, algorithm = algoritmo)
  lim <- if (is.null(limiar)) attr(bs, "threshold") else limiar
  media <- bnlearn::averaged.network(bs, threshold = lim)
  arcos <- as.data.frame(bnlearn::arcs(media), stringsAsFactors = FALSE)
  tab <- merge(arcos, bs, by = c("from", "to"), all.x = TRUE, sort = FALSE)
  names(tab) <- c("de", "para", "forca", "direcao")
  tab <- tab[order(-tab$forca, -tab$direcao), ]; rownames(tab) <- NULL
  und <- bnlearn::undirected.arcs(media)
  if (nrow(und)) tab$direcao_indeterminada <- paste(tab$de, tab$para) %in% paste(und[, 1], und[, 2])
  ext <- list(`Todas as arestas candidatas (força >= 0,10)` = { b <- as.data.frame(bs); b <- b[b$strength >= .10, ]; b <- b[order(-b$strength), ]
    names(b) <- c("de", "para", "forca", "direcao"); rownames(b) <- NULL; b })
  if (comparar_algoritmos) {
    alg <- c("hc", "tabu", "pc.stable", "mmhc")
    cv <- lapply(alg, function(a) tryCatch({ set.seed(semente); r <- bnlearn::bn.cv(d, bn = a, k = 5, runs = 3); mean(vapply(r, function(x) mean(bnlearn::loss(x)), 0)) },
                                          error = function(e) NA_real_))
    ext$`Validação cruzada (perda esperada; menor é melhor)` <- data.frame(algoritmo = alg, perda = unlist(cv))
  }
  en <- identical(idioma, "en")
  apa <- if (en) paste0("A ", if (tipo == "discreta") "discrete" else "Gaussian", " Bayesian network was learned with the ", algoritmo, " algorithm (bnlearn) and averaged over ",
                        reamostras, " bootstrap samples (threshold = ", .pr_num(lim, 2, FALSE, decimal), "): ", nrow(tab), " arcs retained.") else
    paste0("Uma rede bayesiana ", if (tipo == "discreta") "discreta" else "gaussiana", " foi aprendida com o algoritmo ", algoritmo, " (bnlearn) e promediada em ",
           reamostras, " reamostras bootstrap (limiar = ", .pr_num(lim, 2, FALSE, decimal), "): ", nrow(tab), " arestas retidas.")
  .ee_resultado(paste0("Rede bayesiana (", tipo, ", ", algoritmo, ")"), n = nrow(d), n_excluidos = dd$n_excluidos, tabela = tab,
                decisoes = c(paste0("Rede média de ", reamostras, " reamostras bootstrap (boot.strength + averaged.network), semente ", semente, "."),
                             paste0("Limiar de força ", if (is.null(limiar)) "estimado pelo bnlearn" else "informado", ": ", round(lim, 3), "."),
                             "Força = proporção de redes com a aresta; direção = proporção com o sentido indicado, dado que a aresta aparece.",
                             "Sentido de aresta não é efeito causal: dados observacionais transversais só identificam a classe de equivalência."),
                avisos = if (nrow(dd$d) < 10 * length(vars)) "Amostra pequena para o número de variáveis: forças instáveis; ler só as arestas com força alta." else character(),
                apa = apa, extra = list(tabelas_extra = ext, rede = media, forcas = bs))
}

# ---------------------------------------------------------------- redes temporais (v1.7)
# rede_temporal(): rede de séries temporais intensivas (EMA, diário) pelo modelo GVAR de
# lag 1 (Isvoranu et al., 2022, cap. 9 e 10). Uma pessoa: graphicalVAR (regularizado,
# EBIC com gamma), rede temporal em correlações parciais dirigidas (PDC, de t - 1 para t) e
# contemporânea em correlações parciais (PCC). Várias pessoas: mlVAR (dois passos, lme4),
# efeitos fixos temporais com EP e p, rede contemporânea e entre pessoas, e o DP dos efeitos
# aleatórios (quanto a rede varia entre pessoas). Destendenciar (padrão): tira a tendência
# linear no tempo de cada variável (por pessoa), recomendação do capítulo 10 (Epskamp, van
# Borkulo et al., 2018: destendenciar tudo ou só o significativo deu o mesmo; não
# destendenciar baixou a especificidade). Redes por pessoa do mlVAR e os dados ficam no
# objeto ($modelo); não imprimir com dado de participante.
.ee_rede_destendencia <- function(d, vars, id, tempo) {
  grupos <- if (is.null(id)) list(seq_len(nrow(d))) else split(seq_len(nrow(d)), d[[id]])
  tend <- list()
  for (ix in grupos) {
    tt <- if (is.null(tempo)) seq_along(ix) else d[[tempo]][ix]
    for (v in vars) {
      y <- d[[v]][ix]; ok <- !is.na(y) & !is.na(tt)
      if (sum(ok) < 5 || stats::var(tt[ok]) == 0) next
      f <- stats::lm(y[ok] ~ tt[ok]); p <- summary(f)$coefficients
      if (nrow(p) >= 2) tend[[length(tend) + 1]] <- data.frame(variavel = v, p = p[2, 4])
      y[ok] <- stats::residuals(f) + mean(y[ok]); d[[v]][ix] <- y
    }
  }
  t <- if (length(tend)) do.call(rbind, tend) else data.frame(variavel = character(), p = numeric())
  list(d = d, tendencias = t)
}

rede_temporal <- function(dados, vars, id = NULL, tempo = NULL, dia = NULL, bip = NULL, destendenciar = TRUE, gamma = .5,
                          aleatorios = c("orthogonal", "correlated", "fixed"), alfa = .05, decimal = ".", idioma = "pt") {
  aleatorios <- match.arg(aleatorios)
  if (length(vars) < 2) stop("rede temporal pede pelo menos 2 variáveis", call. = FALSE)
  falta <- setdiff(c(vars, id, tempo, dia, bip), names(dados)); if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "), call. = FALSE)
  d <- as.data.frame(dados)[unique(c(vars, id, tempo, dia, bip))]
  for (v in vars) if (!is.numeric(d[[v]])) stop("as variáveis da rede precisam ser numéricas: ", v, call. = FALSE)
  ord_cols <- c(id, dia, bip, tempo)
  if (length(ord_cols)) d <- d[do.call(order, unname(as.list(d[ord_cols]))), , drop = FALSE]
  ids <- if (is.null(id)) "1" else unique(as.character(d[[id]]))
  n_p <- length(ids); n_t <- if (is.null(id)) nrow(d) else as.numeric(table(d[[id]]))
  avisos <- character(); decisoes <- character()
  tend <- NULL
  tempo_ord <- tempo
  if (destendenciar) {
    dt <- .ee_rede_destendencia(d, vars, id, tempo_ord); d <- dt$d; tend <- dt$tendencias
    sig <- unique(tend$variavel[tend$p < .05])
    decisoes <- c(decisoes, paste0("Tendência linear no tempo retirada de todas as variáveis", if (n_p > 1) ", pessoa por pessoa" else "", " (Isvoranu et al., 2022, cap. 10)",
                                   if (length(sig)) paste0("; tendência significativa em ", paste(sig, collapse = ", "), if (n_p > 1) " (em pelo menos uma pessoa)" else "") else "", "."))
  }
  if (n_p == 1) {
    .pr_exige("graphicalVAR")
    if (n_t < 100) avisos <- c(avisos, paste0("Só ", n_t, " pontos no tempo: Epskamp, Waldorp et al. (2018) recuperaram redes esparsas de 8 nós com 100 pontos; com menos, ler a rede como exploratória."))
    if (n_t < 3 * length(vars)^2 / 2) avisos <- c(avisos, "Muitos parâmetros para o número de pontos: a regularização vai zerar arestas fracas, e aresta ausente não é evidência de ausência.")
    args <- list(data = d, vars = vars, gamma = gamma, verbose = FALSE)
    if (!is.null(dia)) args$dayvar <- dia
    if (!is.null(bip)) args$beepvar <- bip
    avs_g <- character()
    g <- withCallingHandlers(suppressMessages(do.call(graphicalVAR::graphicalVAR, args)),
                             warning = function(w) { avs_g <<- c(avs_g, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (any(grepl("converge", avs_g, ignore.case = TRUE))) avisos <- c(avisos, "O graphicalVAR avisou falta de convergência em parte da grade de penalização: conferir se a solução escolhida pelo EBIC é estável (rodar com outro gamma).")
    PDC <- g$PDC; PCC <- g$PCC
    tem <- which(PDC != 0, arr.ind = TRUE)
    temporal <- data.frame(de = rownames(PDC)[tem[, 1]], para = colnames(PDC)[tem[, 2]], pdc = PDC[tem], row.names = NULL)
    temporal$tipo <- ifelse(temporal$de == temporal$para, "autorregressão", "cruzada")
    temporal <- temporal[order(-abs(temporal$pdc)), , drop = FALSE]
    up <- which(upper.tri(PCC) & PCC != 0, arr.ind = TRUE)
    contemp <- data.frame(v1 = rownames(PCC)[up[, 1]], v2 = colnames(PCC)[up[, 2]], pcc = PCC[up], row.names = NULL)
    contemp <- contemp[order(-abs(contemp$pcc)), , drop = FALSE]
    dens <- c(temporal = mean(PDC != 0), contemporanea = mean(PCC[upper.tri(PCC)] != 0))
    decisoes <- c(decisoes, paste0("graphicalVAR (Abegaz & Wit, 2013; Rothman et al., 2010, pelo capítulo 10), seleção por EBIC com gamma = ", gamma,
                                   if (!is.null(dia)) "; pares que cruzam a noite cortados por dia" else "", if (!is.null(bip)) "; bipe para respeitar medidas faltantes" else "", "."),
                  "Rede temporal em PDC (linha = t - 1, coluna = t; inclui autorregressões); contemporânea em PCC, o que sobra correlacionado no mesmo bipe depois do temporal.")
    tab <- temporal
    en <- identical(idioma, "en")
    apa <- if (en) paste0("A lag-1 graphical VAR network (graphicalVAR, EBIC γ = ", gamma, ") on ", n_t, " time points of ", length(vars), " variables retained ",
                          sum(temporal$tipo == "cruzada"), " cross-lagged and ", sum(temporal$tipo == "autorregressão"), " autoregressive temporal edges and ", nrow(contemp), " contemporaneous edges.") else
      paste0("A rede GVAR de lag 1 (graphicalVAR, EBIC γ = ", gamma, ") com ", n_t, " pontos no tempo de ", length(vars), " variáveis reteve ", sum(temporal$tipo == "cruzada"),
             " arestas temporais cruzadas, ", sum(temporal$tipo == "autorregressão"), " autorregressões e ", nrow(contemp), " arestas contemporâneas.")
    return(.ee_resultado("Rede temporal de uma pessoa (GVAR de lag 1)", n = n_t, tabela = tab, decisoes = decisoes, avisos = avisos, apa = apa,
                         extra = list(tabelas_extra = list(`Rede contemporânea (PCC)` = contemp, `Densidade` = data.frame(rede = names(dens), densidade = unname(dens))),
                                      temporal = PDC, contemporanea = PCC, tendencias = tend, modelo = g)))
  }
  .pr_exige("mlVAR")
  if (any(n_t < 20)) avisos <- c(avisos, paste0(sum(n_t < 20), " pessoa(s) com menos de 20 medidas: a centralização na média da pessoa fica imprecisa (o capítulo 10 pede pelo menos cerca de 20 por pessoa)."))
  if (n_p < 30) avisos <- c(avisos, paste0("Só ", n_p, " pessoas: efeitos fixos e DP dos aleatórios imprecisos; a rede entre pessoas, com ", n_p, " médias, é a mais frágil."))
  args <- list(data = d, vars = vars, idvar = id, lags = 1, temporal = aleatorios, contemporaneous = if (aleatorios == "fixed") "fixed" else aleatorios, verbose = FALSE)
  if (!is.null(dia)) args$dayvar <- dia
  if (!is.null(bip)) args$beepvar <- bip
  m <- suppressWarnings(suppressMessages(do.call(mlVAR::mlVAR, args)))
  B <- m$results$Beta
  ij <- expand.grid(para = seq_along(vars), de = seq_along(vars))
  temporal <- data.frame(de = vars[ij$de], para = vars[ij$para], b_fixo = B$mean[cbind(ij$para, ij$de, 1)], ep = B$SE[cbind(ij$para, ij$de, 1)],
                         p = B$P[cbind(ij$para, ij$de, 1)], dp_aleatorio = if (!is.null(B$SD)) B$SD[cbind(ij$para, ij$de, 1)] else NA_real_, row.names = NULL)
  temporal$tipo <- ifelse(temporal$de == temporal$para, "autorregressão", "cruzada")
  temporal$significativo <- temporal$p < alfa
  utils::capture.output(sm <- suppressWarnings(summary(m, show = c("contemporaneous", "between"), print.summary = FALSE)))
  contemp <- sm$contemporaneous; entre <- sm$between
  decisoes <- c(decisoes, paste0("mlVAR (Epskamp, Deserno et al., pelo capítulo 10): estimação em dois passos com lme4, efeitos aleatórios ", aleatorios,
                                 " (orthogonal é o mais rápido e o recomendado com muitas variáveis; correlated estima as covariâncias entre os aleatórios)."),
                "Três redes: temporal (efeitos fixos de t - 1 para t, com EP e p), contemporânea (dentro da pessoa, no mesmo bipe) e entre pessoas (médias das pessoas). O DP do efeito aleatório diz quanto a aresta temporal muda de pessoa para pessoa.",
                "Arestas pelo p dos efeitos fixos (sem regularização); com muitas variáveis, considerar correção para comparações múltiplas antes de ler aresta por aresta.")
  ncr <- sum(temporal$significativo & temporal$tipo == "cruzada")
  en <- identical(idioma, "en")
  apa <- if (en) paste0("A multilevel lag-1 VAR (mlVAR; ", n_p, " participants, median ", stats::median(n_t), " time points each) showed ", ncr, " significant cross-lagged fixed effects (p < ",
                        .pr_num(alfa, 2, FALSE, decimal), ") among ", length(vars), " variables; contemporaneous and between-person networks were also estimated.") else
    paste0("O VAR multinível de lag 1 (mlVAR; ", n_p, " participantes, mediana de ", stats::median(n_t), " pontos por pessoa) mostrou ", ncr, " efeitos fixos temporais cruzados significativos (p < ",
           .pr_num(alfa, 2, FALSE, decimal), ") entre ", length(vars), " variáveis; as redes contemporânea e entre pessoas saem ao lado.")
  .ee_resultado("Rede temporal multinível (mlVAR)", n = n_p, tabela = temporal, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(tabelas_extra = list(`Rede contemporânea` = contemp, `Rede entre pessoas` = entre,
                                                  `Pontos por pessoa` = data.frame(pessoas = n_p, minimo = min(n_t), mediana = stats::median(n_t), maximo = max(n_t))),
                             tendencias = tend, modelo = m))
}
