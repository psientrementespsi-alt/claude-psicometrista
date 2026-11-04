# est_mee.R: modelagem por equações estruturais além da AFC (v1.5). Carregar por carregar.R.
#
#   mee()                       modelo estrutural completo (mensuração + regressões entre
#                               latentes), efeitos indiretos definidos com := e IC por
#                               bootstrap percentil das estimativas padronizadas, como nas
#                               aulas de MEE no R
#   mee_multigrupo()            moderação por variável categórica: configural, métrica,
#                               escalar e igualdade das regressões, com teste de cada
#                               caminho liberado (group.partial)
#   interacao_latente()         interação entre fatores latentes por indicadores-produto
#                               com centralização residual (semTools::indProd) e inclinações
#                               simples (probe2WayRC)
#   curva_crescimento()         curva de crescimento latente (lavaan::growth), linear ou de
#                               base latente, com resíduos adjacentes correlacionados
#   invariancia_longitudinal()  invariância da medida entre ondas (measEq.syntax com
#                               longFacNames), passo prévio da curva de segunda ordem
#
# Mesma convenção da AFC: WLSMV com item ordinal (parametrização theta), MLR com FIML
# nos demais; a versão dos índices é a mesma em todos os modelos comparados.

.mee_vars <- function(modelo) {
  pt <- lavaan::lavaanify(modelo)
  lat <- unique(pt$lhs[pt$op == "=~"])
  def <- unique(pt$lhs[pt$op == ":="])
  ops <- pt$op %in% c("=~", "~", "~~", "~1")
  obs <- setdiff(unique(c(pt$lhs[ops], pt$rhs[ops & pt$op != "~1"])), c(lat, def, ""))
  ind <- unique(pt$rhs[pt$op == "=~" & !pt$rhs %in% lat])
  list(pt = pt, latentes = lat, indicadores = ind, observadas = obs)
}

.mee_ordinal <- function(dados, v, ordinal) {
  if (is.character(ordinal)) return(ordinal)
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(dados, v$indicadores)
  if (isTRUE(ordinal)) v$indicadores else character()
}

.mee_ajusta <- function(modelo, dados, ord, estimador, funcao = lavaan::sem, ...) {
  if (length(ord)) funcao(modelo, data = dados, ordered = ord, estimator = estimador, parameterization = "theta",
                          missing = .pr_ausentes(TRUE, estimador), ...) else
    funcao(modelo, data = dados, estimator = estimador, missing = .pr_ausentes(FALSE, estimador), ...)
}

# Tabela de ajuste comparativa com a mesma versão dos índices em todos os modelos
# (robusta se existir em todos, senão escalonada, senão padrão) e teste de diferença
# entre modelos vizinhos (lavTestLRT: diferença escalonada com WLSMV e MLR).
.mee_compara <- function(fits) {
  tabs <- Map(tabela_ajuste, fits, names(fits))
  versao <- intersect(c("robusta", "escalonada", "padrão"), Reduce(intersect, lapply(tabs, `[[`, "versao")))[1]
  tab <- do.call(rbind, lapply(tabs, function(t) t[t$versao == versao, ]))
  gl <- vapply(fits, function(f) unname(lavaan::fitMeasures(f, "df")), numeric(1))
  tab$delta_cfi <- c(NA, diff(tab$cfi)); tab$delta_rmsea <- c(NA, diff(tab$rmsea)); tab$delta_gl <- c(NA, diff(gl))
  tab$p_dif_qui2 <- c(NA, vapply(seq_along(fits)[-1], function(i) {
    if (gl[i] == gl[i - 1]) return(NA_real_)
    tryCatch(suppressWarnings(lavaan::lavTestLRT(fits[[i - 1]], fits[[i]])[2, "Pr(>Chisq)"]), error = function(e) NA_real_)
  }, numeric(1)))
  rownames(tab) <- NULL
  tab
}

.mee_linha_ajuste <- function(fit, rotulo = "modelo") {
  t <- tabela_ajuste(fit, rotulo)
  t[t$versao == intersect(c("robusta", "escalonada", "padrão"), t$versao)[1], ]
}

mee <- function(modelo, dados, ordinal = NULL, estimador = NULL, bootstrap = 0, semente = 2026, nivel = .95,
                decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  v <- .mee_vars(modelo)
  falta <- setdiff(v$observadas, names(dados))
  if (length(falta)) stop("variáveis do modelo ausentes no banco: ", paste(falta, collapse = ", "))
  if (!any(v$pt$op == "~")) stop("o modelo não tem regressão (~): para só mensuração, ajustar_afc()")
  ord <- .mee_ordinal(dados, v, ordinal)
  if (is.null(estimador)) estimador <- if (length(ord)) "WLSMV" else "MLR"
  fit <- suppressWarnings(.mee_ajusta(modelo, dados, ord, estimador))
  if (!lavaan::lavInspect(fit, "converged")) stop("o modelo não convergiu")
  ss <- lavaan::standardizedSolution(fit, level = nivel)
  pe <- lavaan::parameterEstimates(fit, level = nivel)
  sel <- function(x) x[x$op %in% c("~", ":="), ]
  pe <- sel(pe); ss <- sel(ss)
  tab <- data.frame(parametro = ifelse(pe$op == ":=", pe$lhs, paste(pe$lhs, "<-", pe$rhs)), tipo = ifelse(pe$op == ":=", "definido", "regressão"),
                    b = pe$est, ep = pe$se, z = pe$z, p = pe$pvalue, beta = ss$est.std, ic_inf = ss$ci.lower, ic_sup = ss$ci.upper,
                    stringsAsFactors = FALSE, row.names = NULL)
  rota_ic <- paste0("IC ", round(100 * nivel), "% das estimativas padronizadas pelo método delta.")
  if (bootstrap > 0) {
    set.seed(semente)
    alvo <- which(lavaan::standardizedSolution(fit, se = FALSE)$op %in% c("~", ":="))
    bt <- suppressWarnings(lavaan::bootstrapLavaan(fit, R = bootstrap, type = "ordinary",
                                                   FUN = function(x) lavaan::standardizedSolution(x, se = FALSE)$est.std[alvo]))
    ok <- stats::complete.cases(bt)
    q <- apply(bt[ok, , drop = FALSE], 2, stats::quantile, probs = c((1 - nivel) / 2, 1 - (1 - nivel) / 2), na.rm = TRUE)
    tab$ic_inf <- q[1, ]; tab$ic_sup <- q[2, ]
    rota_ic <- paste0("IC ", round(100 * nivel), "% percentil de ", sum(ok), " reamostras bootstrap das estimativas padronizadas (semente ", semente,
                      if (sum(!ok)) paste0("; ", sum(!ok), " reamostras sem convergência descartadas") else "", "), como nas aulas.")
  }
  aj <- .mee_linha_ajuste(fit)
  r2 <- lavaan::lavInspect(fit, "rsquare")
  endo <- unique(v$pt$lhs[v$pt$op == "~"])
  r2t <- data.frame(variavel = endo, r2 = unname(r2[endo]), stringsAsFactors = FALSE)
  cl <- lavaan::standardizedSolution(fit, se = FALSE); cl <- cl[cl$op == "=~", ]
  cargas <- data.frame(fator = cl$lhs, indicador = cl$rhs, carga = cl$est.std, stringsAsFactors = FALSE, row.names = NULL)
  def <- tab[tab$tipo == "definido", ]
  efe <- if (nrow(def)) def else tab
  efeito <- data.frame(indice = paste0("β ", efe$parametro), valor = efe$beta, ic_inf = efe$ic_inf, ic_sup = efe$ic_sup,
                       principal = seq_len(nrow(efe)) == 1, stringsAsFactors = FALSE)
  frase_par <- paste(vapply(seq_len(nrow(tab)), function(i) paste0(tab$parametro[i], ": β = ", .pr_num(tab$beta[i], 2, TRUE, decimal), ", ",
                                                                    .ee_ic(tab$ic_inf[i], tab$ic_sup[i], 2, TRUE, decimal, nivel, idioma)), ""),
                     collapse = .ee_sep(decimal))
  en <- identical(idioma, "en")
  apa <- if (en) paste0("The structural equation model (", estimador, ", n = ", lavaan::lavInspect(fit, "nobs"), ") showed ",
                        formatar_ajuste_apa(aj, decimal), " (", aj$versao, " indices). Standardized estimates: ", frase_par, ".") else
    paste0("O modelo de equações estruturais (", estimador, ", n = ", lavaan::lavInspect(fit, "nobs"), ") apresentou ",
           formatar_ajuste_apa(aj, decimal), " (índices na versão ", aj$versao, "). Estimativas padronizadas: ", frase_par, ".")
  avisos <- character()
  if (any(r2t$r2 > 1 | r2t$r2 < 0, na.rm = TRUE) || any(abs(cargas$carga) > 1, na.rm = TRUE)) avisos <- "Solução imprópria (carga padronizada > 1 ou R² fora de [0, 1]): revisar o modelo antes de relatar."
  if (nrow(def) && bootstrap == 0) avisos <- c(avisos, "Efeito indireto com IC pelo método delta: a distribuição do produto é assimétrica; preferir bootstrap (argumento bootstrap = 1000 ou mais).")
  .ee_resultado(paste0("MEE (", estimador, ")"), n = lavaan::lavInspect(fit, "nobs"), tabela = tab, efeito = efeito,
                decisoes = c(paste0("Estimador ", estimador, if (length(ord)) " com indicadores ordinais (parametrização theta)" else "", "."),
                             rota_ic, "Teste de cada efeito pelo z da estimativa não padronizada; o IC é da padronizada."),
                avisos = avisos, apa = apa,
                extra = list(tabelas_extra = list(Ajuste = aj, `R² das endógenas` = r2t, `Cargas padronizadas` = cargas), modelo = fit))
}

mee_multigrupo <- function(modelo, dados, grupo, ordinal = NULL, estimador = NULL, caminhos = NULL, igualar_latentes = FALSE,
                           nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  v <- .mee_vars(modelo)
  if (!grupo %in% names(dados)) stop("grupo ausente no banco")
  dados <- dados[!is.na(dados[[grupo]]), ]
  ng <- table(dados[[grupo]])
  if (any(ng < 5)) stop("grupo com menos de 5 casos")
  ord <- .mee_ordinal(dados, v, ordinal)
  if (is.null(estimador)) estimador <- if (length(ord)) "WLSMV" else "MLR"
  med <- if (length(ord)) c("loadings", "thresholds") else c("loadings", "intercepts")
  lat <- if (igualar_latentes) c("lv.variances", "lv.covariances") else character()
  passos <- list(configural = character(), metrica = "loadings", escalar = med)
  if (igualar_latentes) passos$latentes <- c(med, lat)
  passos$regressoes <- c(med, lat, "regressions")
  ajusta <- function(eq, parcial = NULL) {
    a <- list(group = grupo)
    if (length(eq)) a$group.equal <- eq
    if (!is.null(parcial)) a$group.partial <- parcial
    suppressWarnings(do.call(.mee_ajusta, c(list(modelo, dados, ord, estimador), a)))
  }
  fits <- lapply(passos, ajusta)
  tab <- .mee_compara(fits)
  ref <- if (igualar_latentes) fits$latentes else fits$escalar
  regs <- v$pt[v$pt$op == "~", ]
  if (is.null(caminhos)) caminhos <- paste(regs$lhs, "~", regs$rhs)
  # Cada caminho: o modelo com todas as regressões iguais contra o mesmo modelo com só
  # aquele caminho liberado (as aulas: group.partial).
  ptc <- do.call(rbind, lapply(caminhos, function(cm) {
    f1 <- ajusta(passos$regressoes, parcial = cm)
    p <- tryCatch(suppressWarnings(lavaan::lavTestLRT(f1, fits$regressoes)[2, "Pr(>Chisq)"]), error = function(e) NA_real_)
    data.frame(caminho = cm, p_dif_qui2 = p, stringsAsFactors = FALSE)
  }))
  ss <- lavaan::standardizedSolution(ref, level = nivel); ss <- ss[ss$op == "~", ]
  gl <- lavaan::lavInspect(ref, "group.label")
  betas <- data.frame(grupo = gl[ss$group], caminho = paste(ss$lhs, "~", ss$rhs), beta = ss$est.std, ic_inf = ss$ci.lower, ic_sup = ss$ci.upper,
                      stringsAsFactors = FALSE, row.names = NULL)
  p_reg <- tab$p_dif_qui2[nrow(tab)]
  en <- identical(idioma, "en")
  apa <- if (en) paste0("Moderation by ", grupo, " was tested by multigroup SEM (", estimador, "). Constraining all structural paths to equality ",
                        if (isTRUE(p_reg < 1 - nivel)) "worsened" else "did not significantly worsen", " the fit (Δχ² ", .pr_p(p_reg, decimal), ", ΔCFI = ",
                        .pr_num(tab$delta_cfi[nrow(tab)], 3, FALSE, decimal), ").") else
    paste0("A moderação por ", grupo, " foi testada por MEE multigrupo (", estimador, "). Igualar todos os caminhos estruturais entre os grupos ",
           if (isTRUE(p_reg < 1 - nivel)) "piorou" else "não piorou de forma significativa", " o ajuste (Δχ² ", .pr_p(p_reg, decimal), ", ΔCFI = ",
           .pr_num(tab$delta_cfi[nrow(tab)], 3, FALSE, decimal), ").")
  .ee_resultado(paste0("MEE multigrupo por ", grupo, " (", estimador, ")"), n = sum(ng), tabela = tab, posthoc = ptc,
                decisoes = c(paste0("Sequência: ", paste(names(passos), collapse = " → "), "; a igualdade das regressões só se interpreta com a medida invariante (métrica e escalar)."),
                             if (igualar_latentes) "Variâncias e covariâncias latentes igualadas antes das regressões, como nas aulas (pressuposto extra; se esse passo falhar, a comparação das regressões fica confundida)." else
                               "Regressões comparadas a partir do modelo escalar, sem igualar variâncias latentes (as aulas igualam antes; igualar_latentes = TRUE reproduz).",
                             "Teste por caminho: modelo com todas as regressões iguais contra o mesmo modelo com só aquele caminho livre (group.partial)."),
                apa = apa, extra = list(titulo_posthoc = "Teste de cada caminho (liberado contra igualado):",
                                        tabelas_extra = list(`β por grupo (modelo de referência)` = betas), modelos = fits))
}

interacao_latente <- function(dados, x, z, y, nome_x = "X", nome_z = "Z", nome_y = "Y", covariaveis = NULL,
                              valores_mod = c(-1, 0, 1), estimador = "MLR", nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("lavaan", "semTools")
  if (length(x) != length(z)) stop("indicadores-produto pareados pedem o mesmo número de itens em x e z (Marsh et al., 2004); ",
                                   "reduzir o fator maior a parcelas ou escolher os pares pelos maiores pesos")
  vars <- c(x, z, y, covariaveis)
  falta <- setdiff(vars, names(dados)); if (length(falta)) stop("variáveis ausentes: ", paste(falta, collapse = ", "))
  d <- dados[stats::complete.cases(dados[vars]), vars]
  n_exc <- nrow(dados) - nrow(d)
  d2 <- semTools::indProd(d, var1 = x, var2 = z, match = TRUE, meanC = FALSE, residualC = TRUE, doubleMC = FALSE)
  prods <- setdiff(names(d2), names(d))
  nm_xz <- paste0(nome_x, "x", nome_z)
  med <- function(nm, it) paste(nm, "=~", paste(it, collapse = " + "))
  y_lat <- length(y) > 1
  base <- c(med(nome_x, x), med(nome_z, z), if (y_lat) med(nome_y, y))
  yn <- if (y_lat) nome_y else y
  reg0 <- paste(yn, "~", paste(c(nome_x, nome_z, covariaveis), collapse = " + "))
  reg1 <- paste(yn, "~", paste(c(nome_x, nome_z, nm_xz, covariaveis), collapse = " + "))
  m0 <- paste(c(base, reg0), collapse = "\n")
  m1 <- paste(c(base, med(nm_xz, prods), reg1), collapse = "\n")
  f0 <- suppressWarnings(lavaan::sem(m0, data = d, estimator = estimador, meanstructure = TRUE, std.lv = TRUE))
  f1 <- suppressWarnings(lavaan::sem(m1, data = d2, estimator = estimador, meanstructure = TRUE, std.lv = TRUE))
  pe <- lavaan::parameterEstimates(f1, level = nivel); ss <- lavaan::standardizedSolution(f1, level = nivel)
  pe <- pe[pe$op == "~" & pe$lhs == yn, ]; ss <- ss[ss$op == "~" & ss$lhs == yn, ]
  tab <- data.frame(preditor = pe$rhs, b = pe$est, ep = pe$se, z = pe$z, p = pe$pvalue, beta = ss$est.std, ic_inf = ss$ci.lower, ic_sup = ss$ci.upper,
                    stringsAsFactors = FALSE, row.names = NULL)
  pr <- semTools::probe2WayRC(f1, nameX = c(nome_x, nome_z, nm_xz), nameY = yn, modVar = nome_z, valProbe = valores_mod)
  sl <- pr$SimpleSlope
  inc <- data.frame(moderador = sl[[1]], inclinacao = sl$est, ep = sl$se, z = sl$z, p = sl$pvalue, stringsAsFactors = FALSE, row.names = NULL)
  r2_0 <- unname(lavaan::lavInspect(f0, "rsquare")[yn]); r2_1 <- unname(lavaan::lavInspect(f1, "rsquare")[yn])
  li <- tab[tab$preditor == nm_xz, ]
  efeito <- data.frame(indice = c(paste0("β ", nm_xz), "ΔR²"), valor = c(li$beta, r2_1 - r2_0), ic_inf = c(li$ic_inf, NA), ic_sup = c(li$ic_sup, NA),
                       principal = c(TRUE, FALSE), stringsAsFactors = FALSE)
  en <- identical(idioma, "en")
  apa <- if (en) paste0("The latent interaction ", nm_xz, " (residual-centered product indicators; Little et al., 2006) predicting ", yn, ": b = ",
                        .pr_num(li$b, 2, TRUE, decimal), ", β = ", .pr_num(li$beta, 2, TRUE, decimal), ", ", .pr_p(li$p, decimal), "; ΔR² = ",
                        .pr_num(r2_1 - r2_0, 3, FALSE, decimal), ".") else
    paste0("A interação latente ", nm_xz, " (indicadores-produto com centralização residual; Little et al., 2006) na predição de ", yn, ": b = ",
           .pr_num(li$b, 2, TRUE, decimal), ", β = ", .pr_num(li$beta, 2, TRUE, decimal), ", ", .pr_p(li$p, decimal), "; ΔR² = ",
           .pr_num(r2_1 - r2_0, 3, FALSE, decimal), ".")
  .ee_resultado(paste0("Interação latente ", nome_x, " × ", nome_z, " → ", yn, " (", estimador, ")"), n = nrow(d), n_excluidos = n_exc,
                tabela = tab, efeito = efeito, posthoc = inc,
                decisoes = c("Indicadores-produto pareados pela ordem dos itens, centralizados pelos resíduos (semTools::indProd, residualC = TRUE), como nas aulas.",
                             paste0("Fatores com variância 1 (std.lv) e estrutura de médias; inclinações simples de ", yn, " em ", nome_x, " nos valores ",
                                    paste(valores_mod, collapse = ", "), " DP de ", nome_z, " (probe2WayRC)."),
                             "ΔR²: R² de Y no modelo com a interação menos o do modelo sem ela, ajustados à parte.",
                             "Rota com dado contínuo (ML/MLR); item ordinal com poucas categorias pede outra rota (LMS no Mplus) e sai do escopo testado."),
                apa = apa, extra = list(titulo_posthoc = paste0("Inclinações simples (", nome_x, " → ", yn, "):"),
                                        tabelas_extra = list(Ajuste = .mee_linha_ajuste(f1)), modelo = f1, modelo_sem_interacao = f0))
}

curva_crescimento <- function(dados, ondas, tempos = seq_along(ondas) - 1, forma = c("linear", "livre"), residuos_adjacentes = FALSE,
                              residuos_iguais = FALSE, preditores = NULL, estimador = "MLR", nivel = .95, decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  forma <- match.arg(forma)
  k <- length(ondas)
  if (k < 3) stop("curva de crescimento pede ao menos 3 ondas")
  if (forma == "livre" && k < 4) stop("a curva de base latente com 3 ondas fica sem grau de liberdade para testar o ajuste; usar forma = 'linear' ou mais ondas")
  falta <- setdiff(c(ondas, preditores), names(dados)); if (length(falta)) stop("variáveis ausentes: ", paste(falta, collapse = ", "))
  cs <- if (forma == "linear") paste0(tempos, "*", ondas) else c(paste0("0*", ondas[1]), paste0("1*", ondas[2]), ondas[-(1:2)])
  m <- c(paste("intercepto =~", paste(paste0("1*", ondas), collapse = " + ")),
         paste("inclinacao =~", paste(cs, collapse = " + ")))
  if (residuos_adjacentes) m <- c(m, paste(ondas[-k], "~~", ondas[-1]))
  if (residuos_iguais) m <- c(m, paste0(ondas, " ~~ th*", ondas))
  if (length(preditores)) m <- c(m, paste("intercepto + inclinacao ~", paste(preditores, collapse = " + ")))
  modelo <- paste(m, collapse = "\n")
  fit <- suppressWarnings(lavaan::growth(modelo, data = dados, estimator = estimador, missing = .pr_ausentes(FALSE, estimador)))
  pe <- lavaan::parameterEstimates(fit, level = nivel)
  pega <- function(lhs, op, rhs = "") pe[pe$lhs == lhs & pe$op == op & pe$rhs == rhs, ]
  linhas <- rbind(cbind(parametro = "média do intercepto", pega("intercepto", "~1")), cbind(parametro = "média da inclinação", pega("inclinacao", "~1")),
                  cbind(parametro = "variância do intercepto", pega("intercepto", "~~", "intercepto")),
                  cbind(parametro = "variância da inclinação", pega("inclinacao", "~~", "inclinacao")),
                  cbind(parametro = "covariância intercepto-inclinação", pega("intercepto", "~~", "inclinacao")))
  tab <- data.frame(parametro = linhas$parametro, estimativa = linhas$est, ep = linhas$se, z = linhas$z, p = linhas$pvalue,
                    ic_inf = linhas$ci.lower, ic_sup = linhas$ci.upper, stringsAsFactors = FALSE)
  ss <- lavaan::standardizedSolution(fit, se = FALSE)
  rcor <- ss$est.std[ss$lhs == "intercepto" & ss$op == "~~" & ss$rhs == "inclinacao"]
  extras <- list(Ajuste = .mee_linha_ajuste(fit))
  if (forma == "livre") {
    bl <- pe[pe$lhs == "inclinacao" & pe$op == "=~", ]
    extras$`Cargas da inclinação (base latente)` <- data.frame(onda = bl$rhs, carga = bl$est, ep = bl$se, stringsAsFactors = FALSE)
  }
  if (length(preditores)) {
    rp <- pe[pe$op == "~", ]
    extras$Preditores <- data.frame(desfecho = rp$lhs, preditor = rp$rhs, b = rp$est, ep = rp$se, p = rp$pvalue, stringsAsFactors = FALSE)
  }
  ms <- tab[2, ]
  efeito <- data.frame(indice = c("média da inclinação", "correlação intercepto-inclinação"), valor = c(ms$estimativa, rcor),
                       ic_inf = c(ms$ic_inf, NA), ic_sup = c(ms$ic_sup, NA), principal = c(TRUE, FALSE), stringsAsFactors = FALSE)
  aj <- extras$Ajuste
  en <- identical(idioma, "en")
  unid <- if (forma == "linear") (if (en) "per unit of time" else "por unidade de tempo") else (if (en) "from the first to the second wave (latent basis)" else "da primeira para a segunda onda (base latente)")
  apa <- if (en) paste0("A ", if (forma == "linear") "linear" else "latent-basis", " latent growth curve (", estimador, ", n = ", lavaan::lavInspect(fit, "nobs"),
                        ") showed ", formatar_ajuste_apa(aj, decimal), ". Mean slope = ", .pr_num(ms$estimativa, 2, TRUE, decimal), " ", unid, ", ",
                        .ee_ic(ms$ic_inf, ms$ic_sup, 2, TRUE, decimal, nivel, idioma), ", ", .pr_p(ms$p, decimal), ".") else
    paste0("A curva de crescimento latente ", if (forma == "linear") "linear" else "de base latente", " (", estimador, ", n = ", lavaan::lavInspect(fit, "nobs"),
           ") apresentou ", formatar_ajuste_apa(aj, decimal), ". Inclinação média = ", .pr_num(ms$estimativa, 2, TRUE, decimal), " ", unid, ", ",
           .ee_ic(ms$ic_inf, ms$ic_sup, 2, TRUE, decimal, nivel, idioma), ", ", .pr_p(ms$p, decimal), ".")
  avisos <- character()
  vv <- tab$estimativa[3:4]
  if (any(vv < 0, na.rm = TRUE)) avisos <- "Variância latente negativa (caso Heywood): rever a forma da curva ou fixar a variância da inclinação."
  if (isTRUE(aj$gl == 0)) avisos <- c(avisos, "Modelo saturado (gl = 0): o ajuste não é testável e os índices não informam nada; com poucas ondas, tirar os resíduos adjacentes ou usar a forma linear.")
  .ee_resultado(paste0("Curva de crescimento latente (", forma, ", ", estimador, ")"), n = lavaan::lavInspect(fit, "nobs"), tabela = tab, efeito = efeito,
                decisoes = c(paste0("Intercepto com cargas 1; inclinação ", if (forma == "linear") paste0("com cargas fixas nos tempos ", paste(tempos, collapse = ", ")) else
                                      "com as duas primeiras cargas fixas em 0 e 1 e as demais livres (base latente): a inclinação média é a mudança da onda 1 para a 2", "."),
                             if (residuos_adjacentes) "Resíduos de ondas adjacentes correlacionados, como na segunda curva das aulas." else character(),
                             if (residuos_iguais) "Variâncias residuais iguais entre ondas (equivale ao modelo misto com intercepto e inclinação aleatórios)." else character(),
                             paste0("Ausentes por FIML (", estimador, ").")),
                avisos = avisos, apa = apa, extra = list(tabelas_extra = extras, modelo = fit, sintaxe = modelo))
}

# ondas: lista nomeada com um vetor de itens por onda, na mesma ordem de itens em todas
# as ondas, por exemplo list(T1 = c("a1_t1", "a2_t1"), T2 = c("a1_t2", "a2_t2")).
invariancia_longitudinal <- function(dados, ondas, construto = "F", ordinal = NULL, estimador = NULL, residuos_correlacionados = TRUE,
                                     decimal = ".", idioma = "pt") {
  .pr_exige("lavaan", "semTools")
  if (!is.list(ondas) || length(ondas) < 2) stop("ondas: lista com os itens de cada onda (ao menos duas)")
  k <- unique(lengths(ondas)); if (length(k) != 1) stop("todas as ondas precisam do mesmo número de itens, na mesma ordem")
  itens <- unlist(ondas, use.names = FALSE)
  falta <- setdiff(itens, names(dados)); if (length(falta)) stop("itens ausentes: ", paste(falta, collapse = ", "))
  nf <- paste0(construto, "_", if (is.null(names(ondas))) seq_along(ondas) else names(ondas))
  modelo <- paste(vapply(seq_along(ondas), function(i) paste(nf[i], "=~", paste(ondas[[i]], collapse = " + ")), ""), collapse = "\n")
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(dados, itens)
  if (is.null(estimador)) estimador <- if (ordinal) "WLSMV" else "MLR"
  lf <- stats::setNames(list(nf), construto)
  li <- stats::setNames(lapply(seq_len(k), function(j) vapply(ondas, `[`, "", j)), paste0("item", seq_len(k)))
  passos <- if (ordinal) list(configural = character(), limiares = "thresholds", cargas = c("thresholds", "loadings"),
                              interceptos = c("thresholds", "loadings", "intercepts")) else
    list(configural = character(), metrica = "loadings", escalar = c("loadings", "intercepts"), estrita = c("loadings", "intercepts", "residuals"))
  fits <- lapply(passos, function(eq) {
    a <- list(configural.model = modelo, data = dados, longFacNames = lf, ID.fac = "std.lv", estimator = estimador, return.fit = TRUE,
              missing = .pr_ausentes(ordinal, estimador))
    if (residuos_correlacionados) a$longIndNames <- li
    if (!residuos_correlacionados) a$auto <- 0L
    if (ordinal) a <- c(a, list(ordered = itens, ID.cat = "Wu.Estabrook.2016", parameterization = "theta"))
    if (length(eq)) a$long.equal <- eq
    suppressWarnings(do.call(semTools::measEq.syntax, a))
  })
  tab <- .mee_compara(fits)
  falhou <- which(tab$delta_cfi < -.01 | (!is.na(tab$p_dif_qui2) & tab$p_dif_qui2 < .05))
  en <- identical(idioma, "en")
  apa <- if (en) paste0("Longitudinal measurement invariance of ", construto, " across ", length(ondas), " waves (", estimador, "): ",
                        paste(vapply(seq_len(nrow(tab))[-1], function(i) paste0(tab$modelo[i], " ΔCFI = ", .pr_num(tab$delta_cfi[i], 3, FALSE, decimal)), ""),
                              collapse = .ee_sep(decimal)), ".") else
    paste0("A invariância longitudinal da medida de ", construto, " em ", length(ondas), " ondas (", estimador, "): ",
           paste(vapply(seq_len(nrow(tab))[-1], function(i) paste0(tab$modelo[i], " ΔCFI = ", .pr_num(tab$delta_cfi[i], 3, FALSE, decimal),
                                                                    ", Δχ² ", .pr_p(tab$p_dif_qui2[i], decimal)), ""), collapse = .ee_sep(decimal)), ".")
  .ee_resultado(paste0("Invariância longitudinal de ", construto, " (", estimador, ")"), n = lavaan::lavInspect(fits[[1]], "nobs"), tabela = tab,
                decisoes = c(paste0("Passos: ", paste(names(passos), collapse = " → "), if (ordinal) " (Wu & Estabrook, 2016, parametrização theta)" else "", "."),
                             if (residuos_correlacionados) "Resíduos do mesmo item correlacionados entre ondas (longIndNames), padrão em invariância longitudinal." else character(),
                             "Critério: ΔCFI ≤ -0,01 (Cheung & Rensvold, 2002) lido junto com o teste de diferença, que com n grande rejeita diferenças pequenas."),
                avisos = if (length(falhou)) paste0("Passo(s) que pioraram o ajuste pelo critério: ", paste(tab$modelo[falhou], collapse = ", "),
                                                    ". Investigar invariância parcial antes de comparar médias ou montar a curva de segunda ordem.") else character(),
                apa = apa, extra = list(modelos = fits))
}
