# est_vieses.R: controle de vieses de resposta (v1.4; desejabilidade social e escolha forçada na v1.5). Carregar por scripts/carregar.R.
#
#   controle_aquiescencia()   índice de aquiescência, recentralização, modelo de
#                             interceptos randômicos e MIMIC com o índice (curso
#                             Psicometria Online Academy, módulo de vieses)
#
# O curso ensina três controles para escalas com itens positivos e negativos
# (preferencialmente balanceadas, com pares semânticos opostos):
#   1. índice de aquiescência = média das respostas brutas (sem inverter os negativos) de
#      cada pessoa; numa escala balanceada o traço se cancela e sobra a tendência a
#      concordar. Aqui ele sai centrado no ponto médio da escala (0 = sem aquiescência);
#      `pares` restringe o índice aos itens com par semântico;
#   2. recentralização (ipsatização): resposta centrada no ponto médio menos o índice da
#      pessoa; o curso arredonda e trunca no intervalo da escala para seguir com análise
#      ordinal (arredondar = TRUE). Os escores recentralizados não comparam pessoas entre
#      si no nível médio (todas ficam com a mesma média);
#   3. modelo de interceptos randômicos (Maydeu-Olivares & Coffman, 2006): fator de
#      método com carga fixa 1 em todos os itens (na métrica original, com os negativos
#      sem inverter), ortogonal aos fatores de conteúdo; a variância dele é a variância
#      de aquiescência;
#   4. MIMIC: cada item regredido no índice, com os fatores de conteúdo ortogonais ao
#      índice (sintaxe do curso).
# Os modelos rodam com os itens sem inverter (a carga negativa do item negativo mostra o
# polo). Nenhuma linha de participante é impressa: o índice e os dados recentralizados
# ficam em $indice e $dados_recentralizados para a análise seguinte no computador do
# usuário, e a impressão traz só agregados.

.vi_polos <- function(polo, itens) {
  if (is.null(names(polo))) {
    if (length(polo) != length(itens)) stop("polo: vetor nomeado pelos itens ou na mesma ordem dos itens")
    names(polo) <- itens
  }
  falta <- setdiff(itens, names(polo))
  if (length(falta)) stop("polo sem valor para: ", paste(falta, collapse = ", "))
  p <- polo[itens]
  if (all(p %in% c(0, 1))) p <- ifelse(p == 1, 1, -1)
  if (!all(p %in% c(-1, 1))) stop("polo deve ser 1 (positivo) e 0 ou -1 (negativo)")
  p
}

.vi_sintaxe_conteudo <- function(fatores) {
  paste(vapply(names(fatores), function(f) paste0(f, " =~ ", paste(fatores[[f]], collapse = " + ")), ""), collapse = "\n")
}

controle_aquiescencia <- function(dados, fatores, polo, pares = NULL, minimo = NULL, maximo = NULL,
                                  metodos = c("interceptos", "mimic"), ordinal = NULL, estimador = NULL,
                                  arredondar = TRUE, decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  if (!is.list(fatores) || is.null(names(fatores))) stop("fatores: lista nomeada, por exemplo list(Abertura = c('y1', 'y2'), ...)")
  itens <- unlist(fatores, use.names = FALSE)
  if (anyDuplicated(itens)) stop("item em mais de um fator: o controle aqui pressupõe estrutura simples")
  if (length(metodos)) metodos <- match.arg(metodos, c("interceptos", "mimic"), several.ok = TRUE)
  dd <- .ee_dados(dados, itens, numericas = itens)
  d <- dd$d; n <- nrow(d)
  pl <- .vi_polos(polo, itens)
  if (is.null(minimo)) minimo <- min(unlist(d)); if (is.null(maximo)) maximo <- max(unlist(d))
  meio <- (minimo + maximo) / 2; amp <- (maximo - minimo) / 2
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(d, itens)
  avisos <- character(); decisoes <- character()
  # balanceamento
  bal <- do.call(rbind, lapply(names(fatores), function(f) {
    it <- fatores[[f]]; data.frame(fator = f, itens = length(it), positivos = sum(pl[it] == 1), negativos = sum(pl[it] == -1))
  }))
  bal <- rbind(bal, data.frame(fator = "total", itens = length(itens), positivos = sum(pl == 1), negativos = sum(pl == -1)))
  if (any(bal$negativos == 0)) stop("há fator sem item negativo: sem itens dos dois polos o índice confunde aquiescência com traço")
  if (any(bal$positivos != bal$negativos)) avisos <- c(avisos, paste0("Escala não balanceada (positivos e negativos em número diferente em ",
                                                                     paste(bal$fator[bal$positivos != bal$negativos], collapse = ", "),
                                                                     "): o índice mistura aquiescência com traço; o curso indica a recentralização também nesse caso, ",
                                                                     "mas o controle fica parcial."))
  # índice
  base_ind <- if (is.null(pares)) itens else {
    if (is.null(names(pares))) stop("pares: vetor nomeado item -> número do par semântico")
    intersect(itens, names(pares)[!is.na(pares)])
  }
  centr <- as.matrix(d[, itens]) - meio
  indice <- rowMeans(centr[, base_ind, drop = FALSE])
  rec <- sweep(centr, 1, indice, "-")
  if (arredondar) { rec <- round(rec); rec[rec > amp] <- floor(amp); rec[rec < -amp] <- -floor(amp) }
  rec <- as.data.frame(rec)
  desc_ind <- data.frame(n = n, media = mean(indice), dp = stats::sd(indice), minimo = min(indice), maximo = max(indice),
                         pct_acima = mean(indice > 0), pct_abaixo = mean(indice < 0))
  # escores de conteúdo (negativos invertidos) e correlação com o índice
  inv <- sweep(centr, 2, pl[itens], "*")
  cor_ind <- data.frame(fator = names(fatores), r_indice = vapply(names(fatores), function(f) stats::cor(rowMeans(inv[, fatores[[f]], drop = FALSE]), indice), 0))
  decisoes <- c(decisoes, paste0("Índice de aquiescência = média das respostas sem inverter, centradas no ponto médio (", meio, "; 0 = sem aquiescência), sobre ",
                                 length(base_ind), " itens", if (!is.null(pares)) " com par semântico" else "", "."),
                paste0("Recentralização: resposta centrada menos o índice da pessoa", if (arredondar) paste0(", arredondada e truncada em ±", floor(amp), " como no curso") else "",
                       "; dados em $dados_recentralizados (não impressos) para a AFE seguinte."))
  # modelos
  aj <- NULL; tab_cargas <- NULL; var_acq <- NULL; mimic_tab <- NULL; fits <- list()
  if (length(metodos)) {
    if (is.null(estimador)) estimador <- if (ordinal) "WLSMV" else "MLR"
    ajusta <- function(m, dd2) {
      args <- list(model = m, data = dd2, estimator = estimador)
      if (ordinal) args$ordered <- itens
      if (!ordinal) args$missing <- "fiml"
      suppressWarnings(do.call(lavaan::sem, args))
    }
    m_sem <- .vi_sintaxe_conteudo(fatores)
    fits$sem_controle <- ajusta(m_sem, d)
    if ("interceptos" %in% metodos) {
      m_ri <- paste0(m_sem, "\nAcq =~ ", paste0("1*", itens, collapse = " + "), "\n",
                     paste0(names(fatores), " ~~ 0*Acq", collapse = "\n"))
      fits$interceptos <- ajusta(m_ri, d)
      pe <- lavaan::parameterEstimates(fits$interceptos, standardized = TRUE)
      v <- pe[pe$lhs == "Acq" & pe$op == "~~" & pe$rhs == "Acq", ]
      var_acq <- data.frame(variancia = v$est, ep = v$se, z = v$z, p = v$pvalue,
                            carga_padronizada_media = mean(pe$std.all[pe$lhs == "Acq" & pe$op == "=~"]))
      if (isTRUE(v$pvalue < .05)) decisoes <- c(decisoes, "Variância do fator de aquiescência diferente de zero: há aquiescência a controlar.") else
        avisos <- c(avisos, "Variância do fator de aquiescência não significativa: pouca evidência de aquiescência nesses dados.")
    }
    if ("mimic" %in% metodos) {
      d2 <- d; d2$aq_indice_ <- indice
      m_mi <- paste0(m_sem, "\n", paste(itens, collapse = " + "), " ~ aq_indice_\n", paste0("aq_indice_ ~~ 0*", names(fatores), collapse = "\n"))
      fits$mimic <- ajusta(m_mi, d2)
      pm <- lavaan::parameterEstimates(fits$mimic, standardized = TRUE)
      pm <- pm[pm$op == "~" & pm$rhs == "aq_indice_", ]
      mimic_tab <- data.frame(item = pm$lhs, polo = unname(pl[pm$lhs]), b = pm$est, ep = pm$se, p = pm$pvalue, beta = pm$std.all)
    }
    ind <- c("chisq", "df", "cfi", "tli", "rmsea", "srmr")
    if (ordinal) ind <- c("chisq.scaled", "df.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr") else
      ind <- c("chisq.scaled", "df.scaled", "cfi.robust", "tli.robust", "rmsea.robust", "srmr")
    aj <- do.call(rbind, lapply(names(fits), function(k) {
      fm <- lavaan::fitMeasures(fits[[k]], ind)
      data.frame(modelo = c(sem_controle = "sem controle", interceptos = "interceptos randômicos", mimic = "MIMIC com o índice")[[k]],
                 qui2 = fm[[1]], gl = fm[[2]], cfi = fm[[3]], tli = fm[[4]], rmsea = fm[[5]], srmr = fm[[6]], stringsAsFactors = FALSE)
    }))
    carg <- function(f) { p <- lavaan::parameterEstimates(f, standardized = TRUE); p <- p[p$op == "=~" & p$lhs %in% names(fatores), ]
      stats::setNames(p$std.all, p$rhs)[itens] }
    tab_cargas <- data.frame(item = itens, fator = rep(names(fatores), lengths(fatores)), polo = unname(pl[itens]),
                             sem_controle = unname(carg(fits$sem_controle)))
    if (!is.null(fits$interceptos)) tab_cargas$interceptos <- unname(carg(fits$interceptos))
    if (!is.null(fits$mimic)) tab_cargas$mimic <- unname(carg(fits$mimic))
    decisoes <- c(decisoes, paste0("Modelos com os itens sem inverter, estimador ", estimador, if (ordinal) " (itens ordinais)" else " (FIML para ausentes)",
                                   ". Interceptos randômicos: fator Acq com carga fixa 1 em todos os itens e ortogonal aos de conteúdo (Maydeu-Olivares & Coffman, 2006). ",
                                   "MIMIC: itens regredidos no índice, fatores de conteúdo ortogonais ao índice (sintaxe do curso)."),
                  "Compare o ajuste e as cargas de conteúdo entre os modelos: com aquiescência, o modelo sem controle ajusta pior e as cargas dos itens negativos encolhem.")
  }
  f2 <- function(x, k = 2) .pr_num(x, k, TRUE, decimal)
  apa <- if (identical(idioma, "en")) {
    paste0("Acquiescence was indexed by each respondent's mean raw response centered at the scale midpoint (M = ", f2(desc_ind$media), ", SD = ", f2(desc_ind$dp), ")",
           if (!is.null(var_acq)) paste0("; a random-intercept factor model (Maydeu-Olivares & Coffman, 2006) estimated an acquiescence variance of ", f2(var_acq$variancia, 3),
                                         " (", .pr_p(var_acq$p, decimal), ")") else "", ".")
  } else {
    paste0("A aquiescência foi indexada pela média das respostas brutas de cada pessoa, centrada no ponto médio da escala (M = ", f2(desc_ind$media), ", DP = ", f2(desc_ind$dp), ")",
           if (!is.null(var_acq)) paste0("; o modelo de interceptos randômicos (Maydeu-Olivares & Coffman, 2006) estimou variância de aquiescência de ", f2(var_acq$variancia, 3),
                                         " (", .pr_p(var_acq$p, decimal), ")") else "", ".")
  }
  extras <- list(`Balanceamento por fator` = bal, `Índice de aquiescência (agregado)` = desc_ind, `Correlação do índice com os escores de conteúdo` = cor_ind)
  if (!is.null(tab_cargas)) extras$`Cargas padronizadas de conteúdo por modelo` <- tab_cargas
  if (!is.null(mimic_tab)) extras$`MIMIC: itens no índice` <- mimic_tab
  .ee_resultado("Controle de aquiescência", n = n, n_excluidos = dd$n_excluidos, tabela = aj,
                efeito = if (!is.null(var_acq)) data.frame(indice = "variância do fator de aquiescência", valor = var_acq$variancia,
                                                           ic_inf = var_acq$variancia - 1.96 * var_acq$ep, ic_sup = var_acq$variancia + 1.96 * var_acq$ep, principal = TRUE) else NULL,
                decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(indice = unname(indice), dados_recentralizados = rec, variancia_aquiescencia = var_acq, ajustes = fits,
                             tabelas_extra = extras))
}

# ---------------------------------------------------------------- v1.5: desejabilidade social
# controle_desejabilidade(): controle da desejabilidade social em escala com itens
# balanceados por desejabilidade (quádruplas: itens desejáveis e indesejáveis nos dois polos
# do traço; Bäckström et al., 2009) ou escala neutralizada, como nas aulas do curso de
# desejabilidade social. Compara, no mesmo banco: traço só; traço + interceptos randômicos
# (aquiescência); traço + fator de desejabilidade (DS) ortogonal, com o sinal de cada carga
# preso pelo polo de desejabilidade do item (as aulas: L > .01 nos desejáveis, L < -.01 nos
# indesejáveis). Com `criterio`, a validade do traço contra uma medida externa é estimada
# no modelo sem controle e no com controle (as aulas: ACME).
#   polo_ds: vetor nomeado com +1 (item desejável), -1 (indesejável) ou 0 (neutro, fora do DS)
controle_desejabilidade <- function(dados, itens, polo_ds, fator = "Traco", aquiescencia = TRUE, restringir_sinal = TRUE, criterio = NULL,
                                    ordinal = NULL, estimador = NULL, decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  if (is.null(names(polo_ds)) || !all(itens %in% names(polo_ds))) stop("polo_ds: vetor nomeado com +1, -1 ou 0 para cada item")
  falta <- setdiff(c(itens, criterio), names(dados)); if (length(falta)) stop("variáveis ausentes: ", paste(falta, collapse = ", "))
  pol <- polo_ds[itens]; ds_it <- itens[pol != 0]
  if (length(ds_it) < 3 || !all(c(1, -1) %in% pol)) stop("o fator de desejabilidade pede itens desejáveis e indesejáveis (pelo menos 3 no total)")
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(dados, itens)
  if (is.null(estimador)) estimador <- if (ordinal) "WLSMV" else "MLR"
  ord <- if (ordinal) itens else character()
  tr <- paste(fator, "=~", paste(itens, collapse = " + "))
  acq <- c(paste("Acq =~", paste(paste0("1*", itens), collapse = " + ")), paste(fator, "~~ 0*Acq"))
  lab <- paste0("ds_", seq_along(ds_it))
  ds <- c(paste("DS =~", paste(paste0("NA*", ds_it, " + ", lab, "*", ds_it), collapse = " + ")), "DS ~~ 1*DS", paste(fator, "~~ 0*DS"),
          if (aquiescencia) "DS ~~ 0*Acq",
          if (restringir_sinal) ifelse(pol[pol != 0] > 0, paste(lab, "> .01"), paste(lab, "< -.01")))
  modelos <- list(traco = tr)
  if (aquiescencia) modelos$traco_aquiescencia <- c(tr, acq)
  modelos$traco_desejabilidade <- c(tr, if (aquiescencia) acq, ds)
  reg <- if (!is.null(criterio)) paste(fator, "~", criterio)
  ajusta <- function(m) {
    mm <- paste(c(m, reg), collapse = "\n")
    suppressWarnings(if (length(ord)) lavaan::sem(mm, data = dados, ordered = ord, estimator = estimador, parameterization = "theta",
                                                  missing = .pr_ausentes(TRUE, estimador)) else
      lavaan::sem(mm, data = dados, estimator = estimador, missing = .pr_ausentes(FALSE, estimador)))
  }
  fits <- lapply(modelos, ajusta)
  tab <- .mee_compara(fits)
  cg <- do.call(cbind, lapply(fits, function(f) { s <- lavaan::standardizedSolution(f, se = FALSE); s$est.std[s$op == "=~" & s$lhs == fator][match(itens, s$rhs[s$op == "=~" & s$lhs == fator])] }))
  cargas <- data.frame(item = itens, polo_ds = unname(pol), cg, row.names = NULL, check.names = FALSE)
  sf <- lavaan::standardizedSolution(fits$traco_desejabilidade, se = FALSE)
  cds <- sf[sf$op == "=~" & sf$lhs == "DS", ]
  dsl <- data.frame(item = cds$rhs, polo_ds = unname(pol[cds$rhs]), carga_ds = cds$est.std, sinal_coerente = sign(cds$est.std) == sign(pol[cds$rhs]), row.names = NULL)
  val <- NULL
  if (!is.null(criterio)) {
    val <- do.call(rbind, lapply(names(fits), function(nm) { s <- lavaan::standardizedSolution(fits[[nm]]); s <- s[s$op == "~" & s$lhs == fator & s$rhs == criterio, ]
      data.frame(modelo = nm, beta = s$est.std, ic_inf = s$ci.lower, ic_sup = s$ci.upper, p = s$pvalue) }))
  }
  final <- tab[tab$modelo == "traco_desejabilidade", ]
  en <- identical(idioma, "en")
  apa <- if (en) paste0("Adding an orthogonal social desirability factor", if (aquiescencia) " (with a random-intercept acquiescence factor)" else "",
                        " gave ", formatar_ajuste_apa(final, decimal), "; ", sum(dsl$sinal_coerente), " of ", nrow(dsl), " desirability loadings had the expected sign.") else
    paste0("Com o fator ortogonal de desejabilidade social", if (aquiescencia) " (e o fator de interceptos randômicos para aquiescência)" else "",
           ", o modelo apresentou ", formatar_ajuste_apa(final, decimal), "; ", sum(dsl$sinal_coerente), " de ", nrow(dsl), " cargas de desejabilidade com o sinal esperado.")
  extras <- list(`Cargas padronizadas no traço por modelo` = cargas, `Cargas no fator de desejabilidade` = dsl)
  if (!is.null(val)) extras$`Validade contra o critério por modelo` <- val
  .ee_resultado(paste0("Controle de desejabilidade social (", estimador, ")"), n = lavaan::lavInspect(fits[[1]], "nobs"), tabela = tab,
                decisoes = c("Modelos aninhados: traço; + interceptos randômicos (aquiescência, cargas 1, ortogonal); + desejabilidade (DS ortogonal, variância 1).",
                             if (restringir_sinal) "Cargas do DS presas ao sinal do polo de desejabilidade (desejável > .01, indesejável < -.01), como nas aulas." else
                               "Cargas do DS livres (sem restrição de sinal): o sentido do fator fica indeterminado.",
                             "Itens com polo 0 (neutros) ficam fora do DS.",
                             "Com dado ordinal, o teste de diferença entre modelos com restrição de desigualdade é aproximado; ler o ΔCFI junto."),
                avisos = if (any(!dsl$sinal_coerente)) paste0("Carga de DS no limite da restrição ou com sinal oposto ao polo: ", paste(dsl$item[!dsl$sinal_coerente], collapse = ", "),
                                                              ". Rever o polo de desejabilidade dado pelos juízes.") else character(),
                apa = apa, extra = list(tabelas_extra = extras, modelos = fits))
}

# ---------------------------------------------------------------- v1.5: escolha forçada
# escolha_forcada(): modelo de TRI thurstoniana (Brown & Maydeu-Olivares, 2011) para blocos
# de escolha forçada, com o thurstonianIRT (Bürkner, 2019) e o lavaan (ULSMV), como na aula
# "Modelo de TRI Thurstoniano (R)". Entrada: postos dentro de cada bloco (formato "ranks";
# a recodificação para comparações binárias que as aulas fazem com o RecodeFCit é feita
# pelo make_TIRT_data) ou as comparações já binárias (formato "pairwise").
#   blocos: lista com um elemento por bloco: list(itens = c(...), tracos = c(...), sinais = c(1, -1, 1))
# Confiabilidade empírica por traço como na aula "Análises de precisão da escolha forçada
# (R)" (código das aulas): informação de cada comparação no θ estimado de cada
# pessoa, EP² = diag[(I + Φ⁻¹)⁻¹] e ρ = var(θ̂) / [var(θ̂) + média(EP²)]. Como nas aulas, a
# informação ignora a dependência local entre comparações do mesmo bloco. Os θ são
# calculados e descartados.
escolha_forcada <- function(dados, blocos, formato = c("ranks", "pairwise"), direcao = NULL, estimador = "ULSMV",
                            confiabilidade = TRUE, decimal = ".", idioma = "pt") {
  .pr_exige("thurstonianIRT", "lavaan")
  formato <- match.arg(formato)
  # Postos: 1 = mais preferido, então valor menor indica a escolha ("smaller"). Binário
  # das aulas (RecodeFCit): 1 = primeiro item do par escolhido ("larger").
  if (is.null(direcao)) direcao <- if (formato == "ranks") "smaller" else "larger"
  direcao <- match.arg(direcao, c("smaller", "larger"))
  if (!is.list(blocos) || !length(blocos)) stop("blocos: lista de list(itens, tracos, sinais)")
  bl <- Reduce(`+`, lapply(blocos, function(b) thurstonianIRT::set_block(b$itens, traits = b$tracos, signs = if (is.null(b$sinais)) rep(1, length(b$itens)) else b$sinais)))
  td <- if (inherits(dados, "TIRTdata")) dados else suppressMessages(thurstonianIRT::make_TIRT_data(dados, blocks = bl, direction = direcao, format = formato, range = c(0, 1)))
  fit <- suppressWarnings(thurstonianIRT::fit_TIRT_lavaan(td, estimator = estimador))
  lf <- fit$fit
  if (!lavaan::lavInspect(lf, "converged")) stop("o modelo não convergiu")
  aj <- .mee_linha_ajuste(lf)
  est <- lavaan::lavInspect(lf, "est")
  phi <- est$psi; tracos <- colnames(phi)
  # O modelo identifica cada traço só até o sinal: o lavaan pode convergir para a solução
  # refletida. Pelo sinal declarado de cada item (blocos ou dados do thurstonianIRT),
  # traço com a maioria das cargas contra o sinal declarado é refletido no relato.
  tdd <- as.data.frame(td)
  mapa <- unique(rbind(data.frame(item = tdd$item1, traco = tdd$trait1, sinal = tdd$sign1), data.frame(item = tdd$item2, traco = tdd$trait2, sinal = tdd$sign2)))
  ptab <- lavaan::parTable(lf)
  mapa$carga <- vapply(mapa$item, function(i) { v <- ptab$est[ptab$label == paste0("L", i)][1]
    if (is.na(v)) v <- -ptab$est[ptab$label == paste0("L", i, "n")][1]; v }, 0)
  refl <- vapply(seq_along(tracos), function(j) { m <- mapa[mapa$traco == j, ]; mean(sign(m$carga) == sign(m$sinal)) < .5 }, TRUE)
  sg <- ifelse(refl, -1, 1)
  Rt <- stats::cov2cor(phi) * outer(sg, sg)
  cor_t <- data.frame(traco = tracos, as.data.frame(unclass(Rt)), row.names = NULL, check.names = FALSE)
  pe <- lavaan::parameterEstimates(lf, standardized = TRUE)
  cg <- pe[pe$op == "=~", ]
  cargas <- data.frame(comparacao = cg$rhs, traco = cg$lhs, carga = cg$est, carga_padronizada = cg$std.all, ep = cg$se, row.names = NULL)
  cargas <- cargas[order(cargas$comparacao), ]
  rxx <- NULL
  if (confiabilidade) {
    L <- est$lambda; tau <- est$tau[, 1]; th <- diag(est$theta); Pinv <- solve(phi)
    sc <- as.data.frame(stats::predict(fit))
    Th <- stats::reshape(sc, idvar = "id", timevar = "trait", direction = "wide")
    Th <- as.matrix(Th[, paste0("estimate.", tracos)])
    a <- L / sqrt(th); cte <- tau / sqrt(th)
    ep2 <- t(apply(Th, 1, function(t0) {
      z <- drop(a %*% t0) - cte
      f <- stats::dnorm(z)^2 / pmax(stats::pnorm(z) * stats::pnorm(-z), 1e-300)
      diag(solve(crossprod(a * sqrt(f)) + Pinv))
    }))
    v <- apply(Th, 2, stats::var); m <- colMeans(ep2)
    rxx <- data.frame(traco = tracos, confiabilidade_empirica = v / (v + m), ep_medio = sqrt(m), row.names = NULL)
    rm(Th, sc)
  }
  efeito <- if (!is.null(rxx)) data.frame(indice = paste0("confiabilidade empírica ", rxx$traco), valor = rxx$confiabilidade_empirica,
                                          ic_inf = NA_real_, ic_sup = NA_real_, principal = TRUE)
  en <- identical(idioma, "en")
  rx <- if (!is.null(rxx)) paste(paste0(rxx$traco, " = ", vapply(rxx$confiabilidade_empirica, .pr_num, "", 2, FALSE, decimal)), collapse = .ee_sep(decimal)) else ""
  apa <- if (en) paste0("A Thurstonian IRT model (Brown & Maydeu-Olivares, 2011; thurstonianIRT with lavaan, ", estimador, ") was fitted to ", length(blocos),
                        " forced-choice blocks: ", formatar_ajuste_apa(aj, decimal), if (nzchar(rx)) paste0(". Empirical reliability: ", rx) else "", ".") else
    paste0("O modelo de TRI thurstoniana (Brown & Maydeu-Olivares, 2011; thurstonianIRT com lavaan, ", estimador, ") foi ajustado a ", length(blocos),
           " blocos de escolha forçada: ", formatar_ajuste_apa(aj, decimal), if (nzchar(rx)) paste0(". Confiabilidade empírica: ", rx) else "", ".")
  .ee_resultado(paste0("TRI thurstoniana (", length(blocos), " blocos, ", length(tracos), " traços)"), n = lavaan::lavInspect(lf, "nobs"), tabela = aj, efeito = efeito,
                decisoes = c(paste0("Comparações binárias montadas pelo make_TIRT_data (formato ", formato, ", direção ", direcao, "), ajuste ", estimador, " no lavaan."),
                             "Escores normativos só com itens de sinais mistos nos blocos; com todos os sinais iguais os escores ficam ipsativos (Brown & Maydeu-Olivares, 2011).",
                             if (confiabilidade) "Confiabilidade empírica pela informação de cada comparação no θ estimado (fórmula das aulas), sem a dependência local dentro do bloco." else character(),
                             if (any(refl)) paste0("Traço(s) refletido(s) para concordar com o sinal declarado dos itens: ", paste(tracos[refl], collapse = ", "),
                                                   " (correlações já no sentido declarado).") else "Sinais das cargas concordam com o sinal declarado dos itens em todos os traços.",
                             "Blocos com 3 ou mais itens têm redundâncias entre as comparações, e Brown e Maydeu-Olivares (2012) descrevem a correção dos graus de liberdade por elas; a skill não aplica essa correção, então o χ² e os índices que dependem dos gl ficam aproximados."),
                apa = apa, extra = list(tabelas_extra = list(`Correlações entre traços` = cor_t, `Cargas por comparação` = cargas,
                                                             `Confiabilidade` = rxx), refletidos = tracos[refl], modelo = fit))
}
