# est_construcao.R: construção e validação de instrumentos (complementa a AFC de
# psicometria.R). Carregar por scripts/carregar.R.
#
# Validade de conteúdo (juízes nas linhas, itens nas colunas):
#   validade_conteudo_cvc()   CVC de Hernández-Nieto com correção do erro (Pe = (1/J)^J)
#   validade_conteudo_ivc()   I-IVC, S-IVC/Ave, S-IVC/UA e kappa modificado
#   validade_conteudo_rvc()   RVC de Lawshe com valor crítico binomial exato (Ayre & Scally)
#   rvc_critico()             tabela de valores críticos por número de juízes
#   validade_conteudo_aiken() V de Aiken com IC escore (Penfield & Giacobbi)
# Análise fatorial exploratória (participantes nas linhas, itens nas colunas):
#   afe_adequacao()           KMO geral e por item, Bartlett
#   afe_numero_fatores()      análise paralela por permutação, Hull (CFI e CAF), MAP, PA-MRFA (v1.6) e EGA opcionais
#   afe()                     extração (minres ou WLSMV/MLR), rotação oblíqua, cargas
# Fidedignidade e concordância:
#   confiabilidade_teste_reteste()  ICC(A,1) e (C,1), Pearson, EPM, mudança mínima detectável, Bland-Altman
#   concordancia_juizes()     kappa de Cohen (ponderado), Fleiss, ICC(A,1) e ICC(A,k); AC1/AC2 de Gwet e alfa de Krippendorff (v1.6)
#   spearman_brown(), corrigir_correlacao(), generalizabilidade()   profecia, atenuação e restrição de amplitude, estudos G e D (v1.6)
#   analise_distratores()     opções de múltipla escolha por grupo de escore (v1.6)
#   normas()                  tabela de normas: percentil, z, T, escalar, padrão, estanino, por estrato (v1.7)
#   confiabilidade_ema()      RkF, R1R, RkR e Rc em EMA (pessoa x tempo x item; Shrout & Lane, 2012)
# Análise de itens:
#   analise_classica_itens()  dificuldade, RIT, RIR, ULI e alfa sem o item (Martinková & Hladká, 2023)
# Itens com até 7 valores distintos são tratados como ordinais (correlação policórica),
# a mesma regra de .pr_e_ordinal() usada pela AFC. Nenhuma função imprime linha de
# participante: só agregados por item, por fator ou por matriz.

# ---------------------------------------------------------------- auxiliares
# psych::polychoric calcula os pares com parallel::mcmapply, que abre processos por fork.
# Com BLAS multithread (OpenBLAS com OpenMP, comum em Linux e em servidores), um fork
# pode travar o R sem erro; aconteceu nos evals da v1.3. A policórica não usa sorteio,
# então um núcleo só (options(mc.cores = 1)) dá o mesmo resultado, sem o risco.
.ee_policorica <- function(d, ...) {
  op <- options(mc.cores = 1L); on.exit(options(op), add = TRUE)
  psych::polychoric(d, ...)
}

.ee_construcao_semente <- function(semente) {
  tem <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  antiga <- if (tem) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(semente)
  function() {
    if (tem) assign(".Random.seed", antiga, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }
}

.ee_construcao_juizes <- function(dados, itens, minimo = NULL, maximo = NULL) {
  dd <- .ee_dados(dados, itens, numericas = itens)
  if (nrow(dd$d) < 1) stop("nenhum juiz com avaliação completa dos itens")
  v <- unlist(dd$d)
  if (!is.null(minimo) && any(v < minimo)) stop("há notas abaixo do mínimo da escala (", minimo, ")")
  if (!is.null(maximo) && any(v > maximo)) stop("há notas acima do máximo da escala (", maximo, ")")
  dd
}

# Dependência linear exata entre colunas (escore total junto dos itens, item duplicado,
# item e sua versão invertida, colunas derivadas como resíduos ou médias) deixa a matriz
# singular: não há KMO, Bartlett nem fatoração possíveis, e suavizar esconderia o erro.
# Conferida nos dados brutos (posto da matriz padronizada e menor autovalor da Pearson),
# antes de qualquer correlação policórica.
.ee_construcao_dependencia <- function(d) {
  p <- ncol(d); n <- nrow(d)
  dp <- vapply(d, stats::sd, 0)
  const <- colnames(d)[!is.finite(dp) | dp == 0]
  if (length(const)) stop("Item sem variação (constante) no banco: ", paste(const, collapse = ", "),
                          ". Retirar antes da análise fatorial.", call. = FALSE)
  if (n <= p) stop("Há ", n, " participantes com dados completos para ", p, " itens: com n menor ou igual ao número ",
                   "de itens a matriz de correlações é singular de saída (posto no máximo n - 1), sem relação com ",
                   "dependência entre itens. Análise fatorial pede bem mais participantes do que itens.", call. = FALSE)
  q <- qr(scale(as.matrix(d)), tol = 1e-7)
  ev <- min(eigen(stats::cor(d), symmetric = TRUE, only.values = TRUE)$values)
  if (q$rank < p || ev < 1e-8) {
    dep <- if (q$rank < p) colnames(d)[q$pivot[seq(q$rank + 1, p)]] else character()
    stop("A matriz de correlações é singular (posto ", q$rank, " de ", p, "; menor autovalor ", signif(ev, 2),
         "): há itens linearmente dependentes", if (length(dep)) paste0(", e ", paste(dep, collapse = ", "),
         if (length(dep) == 1) " é combinação exata" else " são combinações exatas", " das demais colunas") else "",
         ". Causa provável: escore total ou de subescala junto dos itens, item duplicado, item junto da sua versão ",
         "invertida, ou colunas derivadas (resíduos, médias). Retirar a coluna derivada e rodar de novo.", call. = FALSE)
  }
  invisible(TRUE)
}

# Correlação para a AFE: policórica com itens ordinais, Pearson nos demais. A policórica
# é estimada par a par e pode sair não positiva definida mesmo sem dependência exata nos
# dados; só nesse caso há suavização (psych::cor.smooth, com o eig.tol padrão, o mesmo
# que psych::polychoric faria), com aviso. A policórica é pedida sem a suavização interna
# do psych para a skill saber quando suavizou. A Pearson não chega aqui singular: a
# checagem de dependência acima para antes.
.ee_construcao_cor <- function(d, ordinal) {
  .ee_construcao_dependencia(d)
  avisos <- character()
  R <- withCallingHandlers(
    if (ordinal) { .pr_exige("psych"); .ee_policorica(d, smooth = FALSE)$rho } else stats::cor(d),
    warning = function(w) { avisos <<- c(avisos, paste0("Aviso da correlação: ", conditionMessage(w))); invokeRestart("muffleWarning") })
  if (anyNA(R)) stop("A matriz de correlações tem valores ausentes (par de itens sem estimativa): conferir categorias ",
                     "raras e itens quase constantes (auditar_itens).", call. = FALSE)
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
    .pr_exige("psych")
    R <- suppressWarnings(psych::cor.smooth(R))
    avisos <- c(avisos, paste0("Matriz de correlações ", if (ordinal) "policóricas " else "", "não positiva definida (sem ",
                               "dependência linear exata nos dados brutos): suavizada com psych::cor.smooth antes do KMO, do ",
                               "Bartlett e da extração. Comum com policórica em amostra pequena ou categoria rara; conferir a ",
                               "auditoria dos itens e relatar a suavização."))
  }
  dimnames(R) <- list(colnames(d), colnames(d))
  list(R = R, tipo = if (ordinal) "policórica" else "Pearson", avisos = unique(avisos))
}

.ee_construcao_kmo <- function(R) {
  Ri <- tryCatch(solve(R), error = function(e) stop("A matriz de correlações não pôde ser invertida (", conditionMessage(e),
                                                    "): itens linearmente dependentes ou quase redundantes. Retirar escore total, ",
                                                    "itens duplicados ou colunas derivadas e rodar de novo.", call. = FALSE))
  Q <- -Ri / sqrt(outer(diag(Ri), diag(Ri)))          # correlações parciais (anti-imagem)
  R0 <- R; diag(R0) <- 0; Q0 <- Q; diag(Q0) <- 0
  list(geral = sum(R0^2) / (sum(R0^2) + sum(Q0^2)),
       item = colSums(R0^2) / (colSums(R0^2) + colSums(Q0^2)))
}

.ee_construcao_bartlett <- function(R, n) {
  p <- ncol(R)
  qui2 <- -(n - 1 - (2 * p + 5) / 6) * as.numeric(determinant(R, logarithm = TRUE)$modulus)
  gl <- p * (p - 1) / 2
  list(qui2 = qui2, gl = gl, p = stats::pchisq(qui2, gl, lower.tail = FALSE))
}

# Rótulos de Kaiser (1974) para o KMO.
.ee_construcao_kmo_leitura <- function(x) {
  as.character(cut(x, c(-Inf, .50, .60, .70, .80, .90, Inf), right = FALSE,
                   labels = c("inaceitável", "miserável", "medíocre", "mediano", "meritório", "maravilhoso")))
}

# ---------------------------------------------------------------- CVC (Hernández-Nieto)
# Por item e critério: CVCi = média / máximo da escala; Pe = (1/J)^J (J = juízes); CVCc =
# CVCi - Pe. Total do critério: média dos CVCi menos a média dos Pe (= média dos CVCc).
# Corte .80 (curso, página "Coeficiente de Validade de Conteúdo (CVC)"). `criterio`:
# nome da coluna que diz a que critério cada linha se refere (formato longo: um juiz
# aparece uma vez por critério); sem ela, um critério só.
validade_conteudo_cvc <- function(dados, itens, criterio = NULL, escala_max = 5, corte = 0.80,
                                  decimal = ".", idioma = "pt") {
  grupos <- if (is.null(criterio)) list(geral = dados) else {
    if (!criterio %in% names(dados)) stop("coluna de critério ausente: ", criterio)
    split(dados, droplevels(as.factor(dados[[criterio]])))
  }
  tabs <- list(); efs <- list(); n_exc <- 0; avisos <- character(); juizes <- integer()
  for (g in names(grupos)) {
    dd <- .ee_construcao_juizes(grupos[[g]], itens, minimo = 1, maximo = escala_max)
    n_exc <- n_exc + dd$n_excluidos
    J <- nrow(dd$d); juizes <- c(juizes, J)
    media <- colMeans(dd$d)
    cvci <- media / escala_max
    pe <- (1 / J)^J
    cvcc <- cvci - pe
    tabs[[g]] <- data.frame(criterio = g, item = itens, juizes = J, media = media, cvc_i = cvci, pe = pe,
                            cvc_c = cvcc, adequado = cvcc >= corte, row.names = NULL, stringsAsFactors = FALSE)
    efs[[g]] <- data.frame(indice = paste0("CVCt (", g, ")"), valor = mean(cvci) - mean(pe), ic_inf = NA_real_,
                           ic_sup = NA_real_, principal = TRUE, stringsAsFactors = FALSE)
    if (J < 3) avisos <- c(avisos, paste0(g, ": ", J, " juízes. Com 2 juízes o Pe é .25 e nem notas máximas passam de .75; o curso pede no mínimo 3."))
  }
  tab <- do.call(rbind, tabs); ef <- do.call(rbind, efs)
  if (is.null(criterio)) tab$criterio <- NULL
  ok_item <- tapply(tab$adequado, tab$item, all)[itens]
  decisoes <- c(paste0("Corte CVCc >= ", corte, " por item e por critério (curso; parte da literatura brasileira aceita .70). ",
                       "A decisão é item a item: o CVC total resume, mas média alta pode esconder item ruim."),
                paste0("Erro Pe = (1/J)^J com J = ", paste(unique(juizes), collapse = "/"), " juízes: ",
                       paste(unique(signif((1 / unique(juizes))^unique(juizes), 3)), collapse = "/"), "."),
                "Diagnóstico do curso: clareza baixa aponta redação; pertinência baixa, utilidade ou redundância; relevância baixa, representação do construto. Os comentários dos juízes dizem o que mudar.")
  sep <- .ee_sep(decimal)
  valores <- vapply(ef$valor, .pr_num, "", casas = 2, zero = FALSE, decimal = decimal)
  totais <- if (is.null(criterio)) valores else
    paste(paste0(sub("^CVCt \\((.*)\\)$", "\\1", ef$indice), " = ", valores), collapse = sep)
  n_ok <- sum(ok_item); n_it <- length(itens)
  apa <- if (idioma == "en") {
    paste0("Content validity was assessed by ", paste(unique(juizes), collapse = "/"), " judges on a 1-", escala_max,
           " scale. Total content validity coefficients (CVC; Hernández-Nieto, 2002) were ", totais, "; ", n_ok, " of ", n_it,
           " items reached CVCc ≥ ", .pr_num(corte, 2, FALSE, decimal), if (!is.null(criterio)) " on all criteria" else "", ".")
  } else {
    paste0("A validade de conteúdo foi avaliada por ", paste(unique(juizes), collapse = "/"), " juízes numa escala de 1 a ", escala_max,
           ". O coeficiente de validade de conteúdo total (CVC; Hernández-Nieto, 2002) foi ", totais, "; ", n_ok, " de ", n_it,
           " itens ", if (n_ok == 1) "alcançou" else "alcançaram", " CVCc ≥ ", .pr_num(corte, 2, FALSE, decimal), if (!is.null(criterio)) " em todos os critérios" else "", ".")
  }
  .ee_resultado("Coeficiente de validade de conteúdo (CVC)", n = max(juizes), n_excluidos = n_exc, tabela = tab,
                efeito = ef, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(itens_adequados = names(ok_item)[ok_item]))
}

# ---------------------------------------------------------------- IVC (I-CVI, S-CVI)
# I-IVC = proporção de juízes que deram nota favorável (3 ou 4 na escala de 4 pontos).
# Corte do curso: .80 (página "Índice de Validade de Conteúdo (IVC)"). Ao lado, a
# recomendação de Polit e Beck (2006): I-IVC = 1 com 3 a 5 juízes e >= .78 com 6 ou mais.
# Kappa modificado (Polit, Beck & Owen, 2007): k* = (I-IVC - Pc) / (1 - Pc), Pc =
# C(J, A)·.5^J. Escala: S-IVC/Ave = média dos I-IVC; S-IVC/UA = proporção de itens com
# concordância universal (I-IVC = 1).
validade_conteudo_ivc <- function(dados, itens, favoravel = c(3, 4), corte = 0.80, corte_escala = 0.90,
                                  decimal = ".", idioma = "pt") {
  dd <- .ee_construcao_juizes(dados, itens)
  J <- nrow(dd$d)
  A <- colSums(sapply(dd$d, function(x) x %in% favoravel))
  ivc <- A / J
  pc <- choose(J, A) * 0.5^J
  kap <- (ivc - pc) / (1 - pc)
  leit_k <- as.character(cut(kap, c(-Inf, .40, .60, .75, Inf), right = FALSE,
                             labels = c("ruim", "regular", "bom", "excelente")))
  pb <- if (J <= 5) ivc >= 1 - 1e-12 else ivc >= .78
  tab <- data.frame(item = itens, juizes = J, favoraveis = A, i_ivc = ivc, kappa_mod = kap, leitura_kappa = leit_k,
                    adequado_curso = ivc >= corte - 1e-12, adequado_polit_beck = pb, row.names = NULL, stringsAsFactors = FALSE)
  ave <- mean(ivc); ua <- mean(ivc >= 1 - 1e-12)
  ef <- data.frame(indice = c("S-IVC/Ave", "S-IVC/UA"), valor = c(ave, ua), ic_inf = NA_real_, ic_sup = NA_real_,
                   principal = c(TRUE, FALSE), stringsAsFactors = FALSE)
  decisoes <- c(paste0("Nota favorável = ", paste(favoravel, collapse = " ou "), "; corte do curso I-IVC >= ", corte,
                       " e S-IVC/Ave >= ", corte_escala, " (S-IVC/Ave = ", round(ave, 3), ")."),
                paste0("Critério de Polit e Beck (2006) ao lado: com ", J, " juízes, ",
                       if (J <= 5) "I-IVC = 1 (todos favoráveis)" else "I-IVC >= .78", ". Os dois discordam com 3 a 5 juízes (4 de 5 = .80 passa no curso e não em Polit e Beck)."),
                "Kappa modificado: rótulos > .74 excelente, .60 a .74 bom, .40 a .59 regular (Polit et al., 2007); corrige a concordância ao acaso, que pesa com poucos juízes.",
                "S-IVC/UA cai quando o número de juízes cresce (basta um desfavorável para o item perder a concordância universal); o curso decide item a item e não usa o total como critério.")
  avisos <- if (J < 3) "Menos de 3 juízes: o I-IVC só assume os valores 0, .5 e 1." else character()
  sep <- .ee_sep(decimal)
  apa <- if (idioma == "en") {
    paste0("Content validity was rated by ", J, " judges (favorable = ", paste(favoravel, collapse = " or "), "). Item-level CVI ranged from ",
           .pr_num(min(ivc), 2, FALSE, decimal), " to ", .pr_num(max(ivc), 2, FALSE, decimal), "; ", sum(tab$adequado_curso), " of ", length(itens),
           " items reached I-CVI ≥ ", .pr_num(corte, 2, FALSE, decimal), sep, "S-CVI/Ave = ", .pr_num(ave, 2, FALSE, decimal), sep,
           "S-CVI/UA = ", .pr_num(ua, 2, FALSE, decimal), " (Polit & Beck, 2006).")
  } else {
    paste0("A relevância dos itens foi avaliada por ", J, " juízes (favorável = ", paste(favoravel, collapse = " ou "), "). O I-IVC variou de ",
           .pr_num(min(ivc), 2, FALSE, decimal), " a ", .pr_num(max(ivc), 2, FALSE, decimal), "; ", sum(tab$adequado_curso), " de ", length(itens),
           " itens ", if (sum(tab$adequado_curso) == 1) "alcançou" else "alcançaram", " I-IVC ≥ ", .pr_num(corte, 2, FALSE, decimal), sep, "S-IVC/Ave = ", .pr_num(ave, 2, FALSE, decimal), sep,
           "S-IVC/UA = ", .pr_num(ua, 2, FALSE, decimal), " (Polit & Beck, 2006).")
  }
  .ee_resultado("Índice de validade de conteúdo (IVC)", n = J, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef,
                decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- RVC (Lawshe)
# RVC = (Ne - N/2) / (N/2), Ne = juízes que marcaram "essencial". Valor crítico pelo método
# de Ayre e Scally (2014): binomial exata unilateral com p = .5; Ne crítico = menor Ne
# com P(X >= Ne) <= alfa. Reproduz a tabela da aula para 5 a 10 juízes.
rvc_critico <- function(n_juizes, alfa = 0.05) {
  do.call(rbind, lapply(n_juizes, function(N) {
    pr <- stats::pbinom((0:N) - 1, N, 0.5, lower.tail = FALSE)      # P(X >= ne), ne = 0..N
    ne <- (0:N)[which(pr <= alfa)[1]]
    if (is.na(ne)) return(data.frame(n_juizes = N, ne_critico = NA_integer_, rvc_critico = NA_real_, p_exato = NA_real_))
    data.frame(n_juizes = N, ne_critico = as.integer(ne), rvc_critico = (ne - N / 2) / (N / 2),
               p_exato = stats::pbinom(ne - 1, N, 0.5, lower.tail = FALSE))
  }))
}

validade_conteudo_rvc <- function(dados, itens, essencial = 3, alfa = 0.05, decimal = ".", idioma = "pt") {
  dd <- .ee_construcao_juizes(dados, itens)
  N <- nrow(dd$d)
  ne <- colSums(dd$d == essencial)
  rvc <- (ne - N / 2) / (N / 2)
  cr <- rvc_critico(N, alfa)
  atende <- if (is.na(cr$rvc_critico)) rep(NA, length(itens)) else rvc >= cr$rvc_critico - 1e-12
  tab <- data.frame(item = itens, juizes = N, essenciais = ne, rvc = rvc, rvc_critico = cr$rvc_critico, atende = atende,
                    row.names = NULL, stringsAsFactors = FALSE)
  mantidos <- if (all(is.na(atende))) logical(length(itens)) else atende
  ivc_lawshe <- if (any(mantidos)) mean(rvc[mantidos]) else NA_real_
  ef <- data.frame(indice = "IVC de Lawshe (média do RVC dos itens mantidos)", valor = ivc_lawshe, ic_inf = NA_real_,
                   ic_sup = NA_real_, principal = TRUE, stringsAsFactors = FALSE)
  avisos <- character()
  if (is.na(cr$rvc_critico)) avisos <- c(avisos, paste0("Com ", N, " juízes nenhum Ne alcança p <= ", alfa,
                                                         " (unilateral): o RVC não tem valor crítico. O mínimo é 5 juízes."))
  if (N > 10) avisos <- c(avisos, paste0("Valor crítico para ", N, " juízes calculado pelo método de Ayre e Scally (2014); ",
                                         "conferido contra a tabela da aula só de 5 a 10 juízes. Conferir na tabela do artigo antes de citar."))
  decisoes <- c(paste0("Essencial = nota ", essencial, ". Valor crítico com ", N, " juízes: Ne >= ", cr$ne_critico, ", RVC >= ",
                       round(cr$rvc_critico, 3), " (binomial exata unilateral, p = .5, α = ", alfa, "; Ayre & Scally, 2014)."),
                "Relevância (IVC) e essencialidade (RVC) são perguntas diferentes: item correto mas redundante passa no IVC e cai no RVC.")
  apa <- if (idioma == "en") {
    paste0("Essentiality was rated by ", N, " judges. Using the exact critical value for ", N, " judges (CVR = ",
           .pr_num(cr$rvc_critico, 2, FALSE, decimal), "; Ayre & Scally, 2014), ", sum(mantidos), " of ", length(itens),
           " items met Lawshe's content validity ratio criterion.")
  } else {
    paste0("A essencialidade dos itens foi avaliada por ", N, " juízes. Com o valor crítico exato para ", N, " juízes (RVC = ",
           .pr_num(cr$rvc_critico, 2, FALSE, decimal), "; Ayre & Scally, 2014), ", sum(mantidos), " de ", length(itens),
           " itens atenderam à razão de validade de conteúdo de Lawshe.")
  }
  .ee_resultado("Razão de validade de conteúdo (RVC de Lawshe)", n = N, n_excluidos = dd$n_excluidos, tabela = tab,
                efeito = ef, decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- V de Aiken
# V = soma(nota - mínimo) / (J·(c - 1)), c = número de categorias. IC escore de Penfield e
# Giacobbi (2004), que é o intervalo de Wilson para a proporção V com J·(c - 1) ensaios.
validade_conteudo_aiken <- function(dados, itens, minimo = 1, maximo = 5, nivel = 0.95, corte = NULL,
                                    decimal = ".", idioma = "pt") {
  dd <- .ee_construcao_juizes(dados, itens, minimo = minimo, maximo = maximo)
  J <- nrow(dd$d); k <- maximo - minimo
  S <- colSums(dd$d - minimo)
  v <- S / (J * k)
  z <- stats::qnorm(1 - (1 - nivel) / 2); nk <- J * k
  raiz <- z * sqrt(4 * nk * v * (1 - v) + z^2)
  lo <- (2 * nk * v + z^2 - raiz) / (2 * (nk + z^2)); hi <- (2 * nk * v + z^2 + raiz) / (2 * (nk + z^2))
  tab <- data.frame(item = itens, juizes = J, v = v, ic_inf = lo, ic_sup = hi, row.names = NULL)
  if (!is.null(corte)) tab$adequado <- lo > corte
  ef <- data.frame(indice = "V médio", valor = mean(v), ic_inf = NA_real_, ic_sup = NA_real_, principal = TRUE)
  decisoes <- c(paste0("Escala de ", minimo, " a ", maximo, " (", k + 1, " categorias); IC ", round(100 * nivel),
                       "% escore de Penfield e Giacobbi (2004)."),
                if (is.null(corte)) "Sem corte: as fontes do usuário não trazem corte para o V; se adotar um, declarar e comparar com o limite inferior do IC (argumento corte)." else
                  paste0("Item adequado quando o limite inferior do IC passa de ", corte, " (corte informado pelo usuário)."))
  apa <- if (idioma == "en") {
    paste0("Aiken's V (Aiken, 1985), rated by ", J, " judges on a ", minimo, "-", maximo, " scale, ranged from ",
           .pr_num(min(v), 2, FALSE, decimal), " to ", .pr_num(max(v), 2, FALSE, decimal), " (score confidence intervals; Penfield & Giacobbi, 2004).")
  } else {
    paste0("O V de Aiken (Aiken, 1985), com ", J, " juízes numa escala de ", minimo, " a ", maximo, ", variou de ",
           .pr_num(min(v), 2, FALSE, decimal), " a ", .pr_num(max(v), 2, FALSE, decimal), " (intervalos escore de Penfield & Giacobbi, 2004).")
  }
  .ee_resultado("V de Aiken", n = J, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef, decisoes = decisoes, apa = apa)
}

# ---------------------------------------------------------------- AFE: adequação
afe_adequacao <- function(dados, itens, ordinal = NULL, decimal = ".", idioma = "pt") {
  dd <- .ee_dados(dados, itens, numericas = itens)
  d <- dd$d; n <- nrow(d)
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(d, itens)
  cr <- .ee_construcao_cor(d, ordinal)
  k <- .ee_construcao_kmo(cr$R); b <- .ee_construcao_bartlett(cr$R, n)
  leit <- .ee_construcao_kmo_leitura(k$item)
  tab <- data.frame(item = itens, msa = unname(k$item), leitura = leit, row.names = NULL, stringsAsFactors = FALSE)
  pres <- data.frame(checagem = c("KMO geral", "Bartlett (esfericidade)"), estatistica = c(k$geral, b$qui2),
                     gl = c(NA, b$gl), p = c(NA, b$p),
                     leitura = c(.ee_construcao_kmo_leitura(k$geral),
                                 if (b$p < .05) "rejeita matriz identidade" else "não rejeita matriz identidade: correlações fracas demais para fatorar"),
                     stringsAsFactors = FALSE)
  avisos <- cr$avisos
  baixos <- itens[k$item < .50]
  if (length(baixos)) avisos <- c(avisos, paste0("MSA abaixo de .50 em: ", paste(baixos, collapse = ", "),
                                                 ". Considerar retirar um item por vez (o de menor MSA) e recalcular."))
  if (k$geral < .50) avisos <- c(avisos, "KMO geral abaixo de .50: a matriz não é adequada para análise fatorial.")
  if (k$geral >= .50 && k$geral < .60) avisos <- c(avisos, paste0(
    "KMO geral na faixa miserável de Kaiser (.50 a .60): fatorabilidade fraca. A frase APA não afirma fatorabilidade; ",
    "rever itens de MSA baixo antes de extrair fatores."))
  if (b$p >= .05) avisos <- c(avisos, paste0(
    "Bartlett não rejeitou a matriz identidade (", .pr_p(b$p, decimal), "): as correlações não se distinguem de zero ",
    "nesta amostra, e a frase APA não afirma fatorabilidade."))
  # Fatorabilidade indicada na frase só com KMO >= .60 (acima da faixa miserável de Kaiser)
  # e Bartlett com p < .05; nos demais casos a frase dá os números e diz que não foi indicada.
  fatoravel <- k$geral >= .60 && b$p < .05
  decisoes <- c(paste0("Correlação ", cr$tipo, ": ", if (ordinal) "itens com até 7 valores distintos (regra de .pr_e_ordinal, a mesma da AFC)" else
    "itens com mais de 7 valores distintos, tratados como contínuos", "."),
    "Rótulos do KMO de Kaiser (1974): .90 maravilhoso, .80 meritório, .70 mediano, .60 medíocre, .50 miserável, abaixo de .50 inaceitável.",
    if (ordinal) "Bartlett sobre a matriz policórica é aproximado (a estatística foi derivada para a correlação de Pearson com normalidade); com n grande quase sempre rejeita, então pesa pouco." else
      "Bartlett com n grande quase sempre rejeita; pesa pouco diante do KMO.",
    paste0("Frase APA: fatorabilidade indicada só com KMO >= .60 e Bartlett com p < .05 (convenção da skill: a faixa ",
           "miserável de Kaiser, .50 a .60, não basta); aqui ", if (fatoravel) "indicada." else "não indicada."))
  sep <- .ee_sep(decimal)
  kmo_txt <- .pr_num(k$geral, 2, FALSE, decimal)
  msa_min <- .pr_num(min(k$item), 2, FALSE, decimal); msa_max <- .pr_num(max(k$item), 2, FALSE, decimal)
  bart_txt <- paste0(.ee_estat("χ²", b$qui2, b$gl, decimal), sep, .pr_p(b$p, decimal))
  rot_en <- as.character(cut(k$geral, c(-Inf, .50, .60, .70, .80, .90, Inf), right = FALSE,
                             labels = c("unacceptable", "miserable", "mediocre", "middling", "meritorious", "marvelous")))
  apa <- if (idioma == "en") {
    if (fatoravel) paste0("Factorability of the ", if (ordinal) "polychoric" else "Pearson", " correlation matrix was supported by KMO = ",
                          kmo_txt, " (item MSA ", msa_min, " to ", msa_max, ") and Bartlett's test, ", bart_txt, ".") else
      paste0("Factorability of the ", if (ordinal) "polychoric" else "Pearson", " correlation matrix was not supported: KMO = ",
             kmo_txt, " (", rot_en, " by Kaiser's, 1974, labels; item MSA ", msa_min, " to ", msa_max, ") and Bartlett's test, ",
             bart_txt, if (b$p < .05) ", rejecting the identity matrix." else ", not rejecting the identity matrix.")
  } else {
    if (fatoravel) paste0("A fatorabilidade da matriz de correlações ", if (ordinal) "policóricas" else "de Pearson", " foi indicada pelo KMO = ",
                          kmo_txt, " (MSA dos itens de ", msa_min, " a ", msa_max, ") e pelo teste de Bartlett, ", bart_txt, ".") else
      paste0("A fatorabilidade da matriz de correlações ", if (ordinal) "policóricas" else "de Pearson", " não foi indicada: KMO = ",
             kmo_txt, " (", .ee_construcao_kmo_leitura(k$geral), " pelos rótulos de Kaiser, 1974; MSA dos itens de ", msa_min, " a ",
             msa_max, ") e teste de Bartlett, ", bart_txt, if (b$p < .05) ", que rejeitou a matriz identidade." else
               ", que não rejeitou a matriz identidade.")
  }
  .ee_resultado("Adequação da matriz para análise fatorial", n = n, n_excluidos = dd$n_excluidos, tabela = tab,
                efeito = data.frame(indice = "KMO", valor = k$geral, ic_inf = NA_real_, ic_sup = NA_real_, principal = TRUE),
                pressupostos = pres, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(correlacao = cr$R, tipo_correlacao = cr$tipo, ordinal = ordinal, kmo = k, bartlett = b,
                             fatoravel = fatoravel))
}

# ---------------------------------------------------------------- AFE: número de fatores
# Análise paralela (Horn, 1965) com dados permutados: cada coluna é embaralhada de forma
# independente (Buja & Eyuboglu, 1992), o que mantém a distribuição de cada item e
# desfaz as correlações; a correlação é a mesma da análise (policórica ou Pearson).
# Retém enquanto o autovalor observado passa do percentil (95, Glorfeld, 1995) dos
# autovalores aleatórios. base = "componentes": autovalores da matriz de correlação (o
# critério original de Horn); "fatores": autovalores da matriz reduzida com a correlação
# múltipla ao quadrado (SMC) na diagonal. As duas saem sempre; `base` escolhe a principal.
.ee_construcao_autovalores <- function(R) {
  Rr <- R; diag(Rr) <- 1 - 1 / diag(solve(R))
  list(comp = eigen(R, symmetric = TRUE, only.values = TRUE)$values,
       fat = eigen(Rr, symmetric = TRUE, only.values = TRUE)$values)
}

.ee_construcao_paralela <- function(d, R, ordinal, n_iter, percentil, semente) {
  restaura <- .ee_construcao_semente(semente); on.exit(restaura())
  obs <- .ee_construcao_autovalores(R)
  p <- ncol(d)
  sim_c <- matrix(NA_real_, n_iter, p); sim_f <- matrix(NA_real_, n_iter, p)
  for (i in seq_len(n_iter)) {
    for (tentativa in 1:5) {
      dp <- as.data.frame(lapply(d, sample))
      Rp <- tryCatch(suppressWarnings(if (ordinal) .ee_policorica(dp)$rho else stats::cor(dp)), error = function(e) NULL)
      if (!is.null(Rp) && !anyNA(Rp) && min(eigen(Rp, symmetric = TRUE, only.values = TRUE)$values) > 1e-8) break
      Rp <- NULL
    }
    if (is.null(Rp)) next
    a <- .ee_construcao_autovalores(Rp)
    sim_c[i, ] <- a$comp; sim_f[i, ] <- a$fat
  }
  lim_c <- apply(sim_c, 2, stats::quantile, probs = percentil, na.rm = TRUE)
  lim_f <- apply(sim_f, 2, stats::quantile, probs = percentil, na.rm = TRUE)
  conta <- function(o, l) { f <- which(!(o > l)); if (length(f)) f[1] - 1 else length(o) }
  list(obs = obs, lim_c = lim_c, lim_f = lim_f, n_comp = conta(obs$comp, lim_c), n_fat = conta(obs$fat, lim_f),
       validas = sum(stats::complete.cases(sim_c)))
}

# ---------------------------------------------------------------- AFE: Hull e MAP
# Hull (Lorenzo-Seva, Timmerman & Kiers, 2011; Timmerman, Lorenzo-Seva & Ceulemans, 2018,
# cap. 11 do Wiley Handbook): ajusta soluções de 0 a J fatores por máxima verossimilhança
# (stats::factanal sobre a mesma matriz da análise), registra o ajuste e os graus de
# liberdade de cada uma, descarta as que não melhoram o ajuste e as que ficam abaixo da
# reta entre vizinhas (contorno convexo superior) e escolhe a de maior razão
#   st_i = [(f_i - f_{i-1}) / (gl_i - gl_{i-1})] / [(f_{i+1} - f_i) / (gl_{i+1} - gl_i)],
# o cotovelo do scree formalizado. J = maior número sugerido pela análise paralela (ou o
# teórico, se maior) + 1, limitado pelo limite de Ledermann (maior número de fatores com
# gl >= 0) e com mínimo de 3: a paralela fixa o teto e protege contra superfatoração.
# Ajuste principal: CFI (Bentler, 1990) com F = tr(S^-1 R) + ln|S| - ln|R| - p,
# qui2 = (N - 1) F, modelo nulo qui2_0 = -(N - 1) ln|R| com gl_0 = p(p - 1)/2; a solução de
# 0 fatores entra com CFI = 0. Secundário: CAF (common part accounted for) = 1 - KMO da
# matriz residual (diagonal 1), o índice indicado quando não há qui-quadrado; a solução de
# 0 fatores entra com 1 - KMO de R. Conferido contra EFAtools::HULL (evals_construcao.R).
.ee_hull_ledermann <- function(p) max(0, floor((2 * p + 1 - sqrt(8 * p + 9)) / 2))

.ee_hull_fml <- function(R, L) {
  S <- tcrossprod(L); diag(S) <- 1
  ls <- determinant(S, logarithm = TRUE); lr <- determinant(R, logarithm = TRUE)
  if (ls$sign <= 0 || lr$sign <= 0) return(NA_real_)
  v <- sum(diag(solve(S, R))) + as.numeric(ls$modulus) - as.numeric(lr$modulus) - ncol(R)
  if (!is.finite(v)) NA_real_ else max(0, v)
}

.ee_hull_cfi <- function(Fml, q, R, n) {
  p <- ncol(R); gl <- ((p - q)^2 - (p + q)) / 2
  if (!is.finite(Fml)) return(NA_real_)
  q2 <- Fml * (n - 1); q20 <- -as.numeric(determinant(R, logarithm = TRUE)$modulus) * (n - 1)
  dm <- max(q2 - gl, 0); d0 <- max(q20 - p * (p - 1) / 2, 0)
  if (gl == 0 || max(dm, d0) == 0) 1 else 1 - dm / max(dm, d0)
}

.ee_hull_caf <- function(R, L) {
  D <- R - tcrossprod(L); diag(D) <- 1
  tryCatch(1 - .ee_construcao_kmo(D)$geral, error = function(e) NA_real_)
}

# Uma solução de q fatores por ML; factanal primeiro (R base), psych::fa (fm = "ml") se o
# otimizador do factanal falhar. Heywood: comunalidade >= .995 (unicidade no limite).
.ee_hull_ajusta <- function(R, q, n) {
  L <- tryCatch(unclass(stats::factanal(covmat = R, factors = q, n.obs = n, rotation = "none")$loadings),
                error = function(e) NULL)
  via <- "factanal"
  if (is.null(L) && requireNamespace("psych", quietly = TRUE)) {
    L <- tryCatch(unclass(suppressWarnings(psych::fa(R, nfactors = q, n.obs = n, fm = "ml", rotate = "none",
                                                     warnings = FALSE))$loadings), error = function(e) NULL)
    via <- "psych::fa"
  }
  if (is.null(L)) return(list(L = NULL, F = NA_real_, heywood = NA, via = NA_character_))
  list(L = L, F = .ee_hull_fml(R, L), heywood = any(rowSums(L^2) >= .995), via = via)
}

# Contorno convexo e razão st, na ordem de complexidade crescente (0, 1, ..., J fatores).
.ee_hull_calc <- function(nfat, f, gl) {
  s <- data.frame(nfat = nfat, f = f, gl = gl)
  h <- s[is.finite(s$f), , drop = FALSE]
  st <- rep(NA_real_, nrow(s))
  if (nrow(h) < 1) return(list(st = st, contorno = rep(FALSE, nrow(s)), escolha = NA_integer_, poucos = TRUE))
  d <- diff(h$f)
  while (length(d) && any(d < 0)) {
    h <- h[c(1, d) > 0, , drop = FALSE]
    if (nrow(h) == 1) break
    d <- diff(h$f)
  }
  i <- 2
  while (i <= nrow(h) - 1) {
    linha <- h$f[i - 1] + (h$f[i + 1] - h$f[i - 1]) / (h$gl[i + 1] - h$gl[i - 1]) * (h$gl[i] - h$gl[i - 1])
    if (h$f[i] <= linha) { h <- h[-i, , drop = FALSE]; i <- 1 }
    i <- i + 1
  }
  if (nrow(h) < 3) {
    escolha <- h$nfat[which.max(h$f)]
  } else {
    sth <- rep(NA_real_, nrow(h))
    for (i in 2:(nrow(h) - 1))
      sth[i] <- ((h$f[i] - h$f[i - 1]) / (h$gl[i] - h$gl[i - 1])) / ((h$f[i + 1] - h$f[i]) / (h$gl[i + 1] - h$gl[i]))
    st[match(h$nfat[2:(nrow(h) - 1)], s$nfat)] <- sth[2:(nrow(h) - 1)]
    escolha <- h$nfat[which.max(sth)]
  }
  list(st = st, contorno = s$nfat %in% h$nfat, escolha = as.integer(escolha), poucos = nrow(h) < 3)
}

.ee_hull <- function(R, n, J) {
  p <- ncol(R)
  aj <- lapply(seq_len(J), function(q) .ee_hull_ajusta(R, q, n))
  cfi <- c(0, vapply(seq_len(J), function(q) .ee_hull_cfi(aj[[q]]$F, q, R, n), 0))
  caf <- c(1 - .ee_construcao_kmo(R)$geral,
           vapply(seq_len(J), function(q) if (is.null(aj[[q]]$L)) NA_real_ else .ee_hull_caf(R, aj[[q]]$L), 0))
  gl <- c(p * (p - 1) / 2, ((p - seq_len(J))^2 - (p + seq_len(J))) / 2)
  hc <- .ee_hull_calc(0:J, cfi, gl); ha <- .ee_hull_calc(0:J, caf, gl)
  tab <- data.frame(fatores = 0:J, gl = gl, cfi = cfi, st_cfi = hc$st, contorno_cfi = hc$contorno,
                    caf = caf, st_caf = ha$st, contorno_caf = ha$contorno,
                    heywood = c(FALSE, vapply(aj, function(a) isTRUE(a$heywood), TRUE)), stringsAsFactors = FALSE)
  list(tabela = tab, n_cfi = hc$escolha, n_caf = ha$escolha, poucos_cfi = hc$poucos, poucos_caf = ha$poucos,
       falhas = sum(vapply(aj, function(a) is.null(a$L), TRUE)),
       fallback = sum(vapply(aj, function(a) identical(a$via, "psych::fa"), TRUE)))
}

# MAP de Velicer (1976): para m = 0, 1, ..., p - 1 componentes principais parcializados, a
# média dos quadrados das correlações parciais fora da diagonal; retém o m do mínimo. A
# versão revisada (Velicer, Eaton & Fava, 2000) usa a quarta potência. m = 0 é a própria R.
.ee_map <- function(R) {
  p <- ncol(R); e <- eigen(R, symmetric = TRUE)
  A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), p)
  m2 <- m4 <- rep(NA_real_, p)
  R0 <- R; diag(R0) <- 0
  m2[1] <- sum(R0^2) / (p * (p - 1)); m4[1] <- sum(R0^4) / (p * (p - 1))
  for (m in seq_len(p - 1)) {
    C <- R - tcrossprod(A[, seq_len(m), drop = FALSE]); dg <- diag(C)
    if (any(dg <= 1e-10)) next
    P <- C / sqrt(outer(dg, dg)); diag(P) <- 0
    m2[m + 1] <- sum(P^2) / (p * (p - 1)); m4[m + 1] <- sum(P^4) / (p * (p - 1))
  }
  list(tabela = data.frame(componentes = 0:(p - 1), map = m2, map4 = m4),
       n = which.min(m2) - 1L, n4 = which.min(m4) - 1L)
}

afe_numero_fatores <- function(dados, itens, ordinal = NULL, n_iter = 100, percentil = 0.95,
                               base = c("componentes", "fatores"), ega = FALSE, hull = TRUE, map = TRUE, mrfa = FALSE,
                               n_teorico = NULL, semente = 2026, decimal = ".", idioma = "pt") {
  base <- match.arg(base)
  dd <- .ee_dados(dados, itens, numericas = itens)
  d <- dd$d; n <- nrow(d)
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(d, itens)
  if (ordinal) .pr_exige("psych")
  if (!is.null(n_teorico) && (!is.numeric(n_teorico) || length(n_teorico) != 1 || n_teorico < 1 || n_teorico != round(n_teorico)))
    stop("n_teorico deve ser um inteiro positivo (número de fatores previsto pela teoria)")
  cr <- .ee_construcao_cor(d, ordinal)
  pa <- .ee_construcao_paralela(d, cr$R, ordinal, n_iter, percentil, semente)
  p <- length(itens)
  tab <- data.frame(posicao = seq_len(p), autovalor_comp = pa$obs$comp, limiar_comp = pa$lim_c,
                    retem_comp = seq_len(p) <= pa$n_comp, autovalor_fat = pa$obs$fat, limiar_fat = pa$lim_f,
                    retem_fat = seq_len(p) <= pa$n_fat)
  sugerido <- if (base == "componentes") pa$n_comp else pa$n_fat
  outro <- if (base == "componentes") pa$n_fat else pa$n_comp
  ef <- data.frame(indice = c(paste0("fatores (paralela, ", base, ")"),
                              paste0("fatores (paralela, ", setdiff(c("componentes", "fatores"), base), ")")),
                   valor = c(sugerido, outro), ic_inf = NA_real_, ic_sup = NA_real_, principal = c(TRUE, FALSE),
                   stringsAsFactors = FALSE)
  avisos <- cr$avisos; decisoes <- character(); extras <- list()
  if (pa$validas < n_iter) avisos <- c(avisos, paste0(n_iter - pa$validas, " reamostras descartadas (matriz sem solução)."))
  linha_ef <- function(indice, valor) data.frame(indice = indice, valor = valor, ic_inf = NA_real_, ic_sup = NA_real_,
                                                principal = FALSE, stringsAsFactors = FALSE)
  # Hull
  n_hull <- NA_integer_; n_hull_caf <- NA_integer_; tab_hull <- NULL
  if (hull) {
    lim <- .ee_hull_ledermann(p)
    if (p < 6 || lim < 1) {
      avisos <- c(avisos, paste0("Hull não calculado: pede pelo menos 6 itens (há ", p, ")."))
    } else {
      J <- max(c(sugerido, n_teorico), na.rm = TRUE) + 1
      if (J > lim) J <- lim
      J <- max(J, min(3, lim))
      hl <- .ee_hull(cr$R, n, J)
      n_hull <- hl$n_cfi; n_hull_caf <- hl$n_caf; tab_hull <- hl$tabela
      ef <- rbind(ef, linha_ef("fatores (Hull-CFI)", n_hull), linha_ef("fatores (Hull-CAF)", n_hull_caf))
      decisoes <- c(decisoes, paste0("Hull (Lorenzo-Seva et al., 2011) com soluções de 0 a ", J, " fatores por máxima verossimilhança ",
                                     "sobre a mesma matriz (teto = paralela", if (!is.null(n_teorico)) " ou número teórico" else "",
                                     " + 1, limite de Ledermann ", lim, "). Índice principal CFI; CAF como secundário. ",
                                     "Tabela das soluções em $tabelas_extra$hull."))
      if (ordinal) decisoes <- c(decisoes, paste0("Com correlação policórica, o qui-quadrado de máxima verossimilhança é aproximado; ",
                                                 "o Hull usa só a posição relativa das soluções, e o capítulo 11 do Wiley Handbook mostra ",
                                                 "o Hull-CFI no modelo ordinal concordando com o Hull-CFI por ULS sobre policóricas."))
      if (hl$falhas > 0) avisos <- c(avisos, paste0(hl$falhas, " solução(ões) do Hull sem convergência na máxima verossimilhança: ",
                                                   "excluída(s) do contorno."))
      if (hl$poucos_cfi) avisos <- c(avisos, "Menos de três soluções no contorno do Hull-CFI: a escolha vira o maior CFI, heurística fraca; pesar a paralela.")
      if (!is.na(n_hull) && n_hull >= 1 && isTRUE(tab_hull$heywood[tab_hull$fatores == n_hull]))
        avisos <- c(avisos, paste0("A solução de ", n_hull, " fator(es) escolhida pelo Hull-CFI tem caso Heywood (comunalidade perto de 1): ",
                                   "conferir as cargas antes de aceitar esse número."))
      if (!is.na(n_hull) && n_hull != sugerido)
        avisos <- c(avisos, paste0("Paralela (", base, ") sugere ", sugerido, " e Hull-CFI sugere ", n_hull, ". O capítulo 11 do Wiley Handbook ",
                                   "recomenda considerar os dois e decidir pela interpretabilidade e plausibilidade teórica das soluções ",
                                   "indicadas; relatar ambos."))
      itens_por_fator <- p / max(1, sugerido)
      if (n <= 100 || itens_por_fator <= 5)
        decisoes <- c(decisoes, paste0("Amostra ou razão itens/fator na faixa em que a paralela superou o Hull na simulação de Lorenzo-Seva ",
                                       "et al. (2011) (N = 100 ou 5 itens por fator; aqui N = ", n, ", ", .pr_num(itens_por_fator, 1, decimal = decimal),
                                       " itens por fator pela paralela): havendo discordância, pesar mais a paralela."))
    }
  }
  # MAP
  n_map <- NA_integer_; n_map4 <- NA_integer_; tab_map <- NULL
  if (map) {
    mp <- .ee_map(cr$R)
    n_map <- mp$n; n_map4 <- mp$n4; tab_map <- mp$tabela
    ef <- rbind(ef, linha_ef("componentes (MAP, Velicer 1976)", n_map), linha_ef("componentes (MAP revisado, 4ª potência)", n_map4))
    decisoes <- c(decisoes, paste0("MAP de Velicer (1976) e revisado (Velicer, Eaton & Fava, 2000), sobre componentes principais. ",
                                   "Critério auxiliar: Lorenzo-Seva et al. (2011) desaconselham o MAP para fatores comuns, porque o ",
                                   "desempenho piora com unicidades altas (Wiley Handbook, cap. 11)."))
  }
  # PA-MRFA (v1.6): análise paralela sobre a variância comum explicada pela análise fatorial
  # de posto mínimo (Timmerman & Lorenzo-Seva, 2011), a versão que o capítulo 11 do Wiley
  # Handbook prefere para itens ordinais; via EFA.MRFA::parallelMRFA, com policórica se os
  # itens forem ordinais e o percentil da paralela daqui.
  n_mrfa <- NA_integer_
  if (mrfa) {
    if (!requireNamespace("EFA.MRFA", quietly = TRUE)) {
      avisos <- c(avisos, "PA-MRFA pedida, mas o pacote EFA.MRFA não está instalado (install.packages(\"EFA.MRFA\")); seguem os demais critérios.")
    } else {
      rest <- .ee_construcao_semente(semente)
      mr <- tryCatch(suppressWarnings(EFA.MRFA::parallelMRFA(as.matrix(d), Ndatsets = n_iter, percent = round(100 * percentil),
                                                             corr = if (ordinal) "Polychoric" else "Pearson", display = FALSE, graph = FALSE)),
                     error = function(e) e)
      rest()
      if (inherits(mr, "error")) avisos <- c(avisos, paste0("PA-MRFA falhou: ", conditionMessage(mr))) else {
        n_mrfa <- as.integer(mr$N_factors_percentiles)
        ef <- rbind(ef, linha_ef("fatores (PA-MRFA)", n_mrfa))
        extras_mrfa <- data.frame(posicao = seq_along(mr$Real_Data), pct_variancia_comum = mr$Real_Data, limiar = mr$Percentile_random)
        decisoes <- c(decisoes, paste0("PA-MRFA (Timmerman & Lorenzo-Seva, 2011; EFA.MRFA), ", n_iter, " bancos aleatórios, percentil ", round(100 * percentil),
                                       ", correlação ", if (ordinal) "policórica" else "de Pearson", ": compara a % de variância comum de cada fator com a do acaso."))
        if (n_mrfa != sugerido) avisos <- c(avisos, paste0("PA-MRFA sugere ", n_mrfa, " e a paralela (", base, ") sugere ", sugerido,
                                                           ". Com itens ordinais e cargas moderadas, o capítulo 11 do Wiley Handbook pesa mais a PA-MRFA."))
      }
    }
  }
  # EGA
  n_ega <- NA
  if (ega) {
    .pr_exige("EGAnet")
    g <- EGAnet::EGA(d, plot.EGA = FALSE, verbose = FALSE)
    n_ega <- g$n.dim
    ef <- rbind(ef, linha_ef("dimensões (EGA)", n_ega))
    decisoes <- c(decisoes, "EGA com os padrões do EGAnet (correlação automática, GLASSO, Walktrap). Estabilidade com bootEGA: ver references/redes-ega.md.")
  } else {
    decisoes <- c(decisoes, "Complemento recomendado pelo curso: EGA e bootEGA (EGAnet; ega = TRUE aqui ou references/redes-ega.md).")
  }
  if (sugerido != outro) avisos <- c(avisos, paste0("As duas bases da análise paralela discordam (", base, ": ", sugerido, "; outra: ", outro,
                                                   "). Decidir pela teoria e pela interpretabilidade das soluções vizinhas, e relatar as duas."))
  if (sugerido == 0) avisos <- c(avisos, "Nenhum autovalor passou do aleatório: os itens não compartilham variância suficiente para extrair fator.")
  decisoes <- c(paste0("Análise paralela por permutação (", n_iter, " reamostras, percentil ", round(100 * percentil), ", correlação ", cr$tipo,
                       ", semente ", semente, "). Base principal: ", base, "."), decisoes)
  extra_en <- c(if (mrfa && !is.na(n_mrfa)) paste0("PA-MRFA (Timmerman & Lorenzo-Seva, 2011) suggested ", n_mrfa),
                if (hull && !is.na(n_hull)) paste0("the Hull method with CFI (Lorenzo-Seva et al., 2011) suggested ", n_hull),
                if (map) paste0("Velicer's MAP suggested ", n_map),
                if (ega) paste0("exploratory graph analysis suggested ", n_ega, " dimension", if (n_ega != 1) "s" else ""))
  extra_pt <- c(if (mrfa && !is.na(n_mrfa)) paste0("a PA-MRFA (Timmerman & Lorenzo-Seva, 2011) sugeriu ", n_mrfa),
                if (hull && !is.na(n_hull)) paste0("o método de Hull com CFI (Lorenzo-Seva et al., 2011) sugeriu ", n_hull),
                if (map) paste0("o MAP de Velicer sugeriu ", n_map),
                if (ega) paste0("a análise exploratória de grafos sugeriu ", n_ega, if (n_ega == 1) " dimensão" else " dimensões"))
  junta <- function(x, e) if (!length(x)) "" else if (length(x) == 1) paste0("; ", x) else
    paste0("; ", paste(x[-length(x)], collapse = ", "), if (e == "en") ", and " else " e ", x[length(x)])
  apa <- if (idioma == "en") {
    paste0("Parallel analysis (Horn, 1965) on ", if (ordinal) "polychoric" else "Pearson", " correlations, with ", n_iter,
           " permuted data sets and the ", round(100 * percentil), "th percentile of random eigenvalues, suggested ", sugerido,
           " factor", if (sugerido != 1) "s" else "", junta(extra_en, "en"), ".")
  } else {
    paste0("A análise paralela (Horn, 1965) sobre correlações ", if (ordinal) "policóricas" else "de Pearson", ", com ", n_iter,
           " bancos permutados e o percentil ", round(100 * percentil), " dos autovalores aleatórios, sugeriu ", sugerido,
           if (sugerido == 1) " fator" else " fatores", junta(extra_pt, "pt"), ".")
  }
  extras <- list()
  if (!is.null(tab_hull)) extras$hull <- tab_hull
  if (!is.null(tab_map)) extras$map <- tab_map
  if (mrfa && !is.na(n_mrfa)) extras$mrfa <- extras_mrfa
  .ee_resultado("Número de fatores (paralela, Hull, MAP)", n = n, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef,
                decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(n_sugerido = sugerido, n_outra_base = outro, n_hull = n_hull, n_hull_caf = n_hull_caf,
                             n_map = n_map, n_map4 = n_map4, n_ega = n_ega, n_mrfa = n_mrfa, tipo_correlacao = cr$tipo,
                             tabelas_extra = if (length(extras)) extras else NULL))
}

# ---------------------------------------------------------------- AFE: extração e rotação
# extracao = "minres": psych::fa sobre a matriz (policórica ou Pearson), o padrão do psych.
# extracao = "wlsmv": lavaan::efa com WLSMV para itens ordinais e MLR para contínuos (os
# mesmos estimadores da AFC), com índices de ajuste. Rotação oblíqua por padrão (oblimin):
# fatores psicológicos raramente são independentes, e com correlação zero a solução
# oblíqua devolve o mesmo que a ortogonal. Fatores ordenados pela variância explicada e
# com sinal tal que a soma das cargas seja positiva.
.ee_construcao_alinha <- function(L, Phi) {
  ss <- colSums(L * (L %*% Phi))
  o <- order(ss, decreasing = TRUE)
  L <- L[, o, drop = FALSE]; Phi <- Phi[o, o, drop = FALSE]
  s <- ifelse(colSums(L) < 0, -1, 1)
  L <- sweep(L, 2, s, `*`); Phi <- Phi * outer(s, s)
  colnames(L) <- rownames(Phi) <- colnames(Phi) <- paste0("F", seq_len(ncol(L)))
  list(L = L, Phi = Phi)
}

afe <- function(dados, itens, n_fatores = NULL, extracao = c("minres", "wlsmv"), rotacao = "oblimin",
                ordinal = NULL, carga_min = 0.40, cruzada_min = 0.30, n_iter = 100, semente = 2026,
                decimal = ".", idioma = "pt") {
  extracao <- match.arg(extracao)
  dd <- .ee_dados(dados, itens, numericas = itens)
  d <- dd$d; n <- nrow(d); p <- length(itens)
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(d, itens)
  .pr_exige("psych", "GPArotation")
  adeq <- afe_adequacao(d, itens, ordinal = ordinal, decimal = decimal, idioma = idioma)
  R <- adeq$correlacao
  decisoes <- character(); avisos <- adeq$avisos; apa_nf <- character()
  if (is.null(n_fatores)) {
    nf <- afe_numero_fatores(d, itens, ordinal = ordinal, n_iter = n_iter, semente = semente, decimal = decimal, idioma = idioma)
    n_fatores <- max(1, nf$n_sugerido)
    decisoes <- c(decisoes, paste0("Número de fatores pela análise paralela (componentes): ", nf$n_sugerido,
                                   " (base fatores: ", nf$n_outra_base, "; Hull-CFI: ", nf$n_hull, "; MAP: ", nf$n_map,
                                   "). Detalhe em afe_numero_fatores()."))
    avisos <- c(avisos, setdiff(nf$avisos, avisos)); apa_nf <- nf$apa
  }
  if (n_fatores >= p) stop("n_fatores precisa ser menor que o número de itens")
  rot <- if (n_fatores == 1) "none" else rotacao
  avisos_ext <- character()
  pega <- function(expr) withCallingHandlers(expr, warning = function(w) {
    avisos_ext <<- c(avisos_ext, conditionMessage(w)); invokeRestart("muffleWarning") })
  ajuste <- NULL; modelo <- NULL
  if (extracao == "minres") {
    fa <- pega(psych::fa(R, nfactors = n_fatores, n.obs = n, fm = "minres", rotate = rot))
    L <- unclass(fa$loadings)[, seq_len(n_fatores), drop = FALSE]
    Phi <- if (is.null(fa$Phi)) diag(n_fatores) else fa$Phi
    modelo <- fa
    est_txt <- paste0("minres sobre correlações ", if (ordinal) "policóricas" else "de Pearson")
    ajuste <- data.frame(indice = "RMSR (resíduos fora da diagonal)", valor = fa$rms, ic_inf = NA_real_, ic_sup = NA_real_, principal = FALSE)
  } else {
    .pr_exige("lavaan")
    est <- if (ordinal) "WLSMV" else "MLR"
    args <- list(data = d, nfactors = n_fatores, estimator = est, rotation = if (rot == "none") "geomin" else rot)
    if (ordinal) args$ordered <- itens
    e <- pega(do.call(lavaan::efa, args))
    modelo <- e[[1]]
    std <- lavaan::lavInspect(modelo, "std")
    L <- std$lambda[, seq_len(n_fatores), drop = FALSE]
    Phi <- if (n_fatores == 1) matrix(1) else std$psi
    est_txt <- paste0(est, if (ordinal) " com itens ordinais" else "")
    ta <- tabela_ajuste(modelo, paste(n_fatores, "fatores"))
    v <- if ("robusta" %in% ta$versao && !is.na(ta$cfi[ta$versao == "robusta"])) "robusta" else ta$versao[1]
    t1 <- ta[ta$versao == v, ]
    ajuste <- data.frame(indice = paste0(c("CFI", "TLI", "RMSEA", "SRMR"), " (", v, ")"),
                         valor = c(t1$cfi, t1$tli, t1$rmsea, t1$srmr), ic_inf = c(NA, NA, t1$rmsea_li, NA),
                         ic_sup = c(NA, NA, t1$rmsea_ls, NA), principal = FALSE)
    decisoes <- c(decisoes, paste0("Índices de ajuste da AFE na versão ", v, " (lavaan); com WLSMV, os cortes de Hu e Bentler (1999) tendem a não acusar desajuste (references/afc-ordinal.md)."))
  }
  rownames(L) <- itens
  al <- .ee_construcao_alinha(L, as.matrix(Phi)); L <- al$L; Phi <- al$Phi
  h2 <- rowSums(L * (L %*% Phi))
  absL <- abs(L)
  princ <- apply(absL, 1, which.max)
  carga_p <- L[cbind(seq_len(p), princ)]
  n_sal <- rowSums(absL >= cruzada_min)
  # Caso Heywood: carga padronizada ou comunalidade igual ou acima de 1 (tolerância
  # numérica 1e-6), ou aviso de Heywood/variância negativa vindo da extração (psych ou
  # lavaan). A solução é imprópria e a frase APA diz isso em vez de descrever as cargas.
  heywood <- h2 >= 1 - 1e-6 | apply(absL, 1, max) >= 1 - 1e-6
  heywood_ext <- any(grepl("Heywood|variances are negative|negative variance", avisos_ext, ignore.case = TRUE))
  improprio <- any(heywood) || heywood_ext
  alerta <- ifelse(heywood, "Heywood (h² ou carga >= 1)",
                   ifelse(n_sal >= 2, "carga cruzada",
                          ifelse(apply(absL, 1, max) < carga_min, "carga baixa", "")))
  tab <- data.frame(item = itens, round(L, 3), h2 = h2, u2 = 1 - h2, fator = colnames(L)[princ], alerta = alerta,
                    row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  ss <- colSums(L * (L %*% Phi))
  ef <- data.frame(indice = c(paste0(colnames(L), ": proporção da variância"), "Total: proporção da variância comum"),
                   valor = c(ss / p, sum(ss) / p), ic_inf = NA_real_, ic_sup = NA_real_,
                   principal = c(rep(FALSE, length(ss)), TRUE), stringsAsFactors = FALSE)
  if (n_fatores > 1) {
    pares <- which(upper.tri(Phi), arr.ind = TRUE)
    ef <- rbind(ef, data.frame(indice = paste0("r(", colnames(Phi)[pares[, 1]], ", ", colnames(Phi)[pares[, 2]], ")"),
                               valor = Phi[pares], ic_inf = NA_real_, ic_sup = NA_real_, principal = FALSE))
  }
  if (!is.null(ajuste)) ef <- rbind(ef, ajuste)
  if (length(avisos_ext)) avisos <- c(avisos, paste0("Aviso da extração: ", unique(avisos_ext)))
  if (improprio) avisos <- c(avisos, paste0(
    "Caso Heywood (carga ou comunalidade igual ou acima de 1", if (any(heywood)) paste0(" em ", paste(itens[heywood], collapse = ", ")) else
      " indicado pela extração", "): solução imprópria, que não deve ser interpretada; rever o número de fatores, itens quase ",
    "redundantes ou a amostra. A frase APA registra a solução imprópria."))
  if (any(alerta == "carga cruzada")) avisos <- c(avisos, paste0("Carga cruzada (>= ", cruzada_min, " em 2 ou mais fatores): ",
                                                                paste(itens[alerta == "carga cruzada"], collapse = ", "), "."))
  if (any(alerta == "carga baixa")) avisos <- c(avisos, paste0("Carga principal abaixo de ", carga_min, ": ",
                                                              paste(itens[alerta == "carga baixa"], collapse = ", "), "."))
  if (!rot %in% c("oblimin", "promax", "geomin", "quartimin", "none")) avisos <- c(avisos, paste0("Rotação '", rot, "': se for ortogonal, ela força fatores não correlacionados; a skill usa oblíqua por padrão."))
  decisoes <- c(decisoes, paste0("Extração ", est_txt, "; rotação ", rot, "; ", n_fatores, if (n_fatores == 1) " fator." else " fatores."),
                paste0("Alertas por item: carga cruzada se >= ", cruzada_min, " em 2 ou mais fatores; carga baixa se a maior carga < ", carga_min,
                       " (convenção da skill, a mesma leitura de .40 da AFC em references/afc-ordinal.md; o curso não dá corte)."),
                "Variância por fator = soma de λ·(λΦ) do padrão, que soma a comunalidade total mesmo com fatores correlacionados.",
                "Exploratória não é sem teoria: aceitar uma solução só se os agrupamentos fizerem sentido para o construto (curso, página 'Estrutura Interna'); depois, AFC em outra amostra.")
  faixa <- range(abs(carga_p))
  cruz <- sum(alerta == "carga cruzada")
  sep <- .ee_sep(decimal)
  ex_pt <- c(minres = "extração por resíduos mínimos", wlsmv = paste("estimação", if (ordinal) "WLSMV" else "MLR"))[[extracao]]
  if (improprio) {
    # Solução imprópria (Heywood): a frase registra o problema e não descreve as cargas
    # como estimativas; os números ficam entre parênteses para rastrear a solução.
    pct <- .pr_num(100 * sum(ss) / p, 1, TRUE, decimal)
    a1 <- .pr_num(faixa[1], 2, FALSE, decimal); a2 <- .pr_num(faixa[2], 2, FALSE, decimal)
    hw <- paste(itens[heywood], collapse = ", ")
    frase <- if (idioma == "en") paste0(
      "An exploratory factor analysis (", c(minres = "minimum residuals extraction", wlsmv = paste(if (ordinal) "WLSMV" else "MLR", "estimation"))[[extracao]],
      " on ", if (ordinal) "polychoric" else "Pearson", " correlations", if (rot != "none") paste0(", ", rot, " rotation") else "",
      ") with ", n_fatores, " factor", if (n_fatores != 1) "s" else "", " yielded an improper solution, with a Heywood case ",
      "(loading or communality at or above 1)", if (nzchar(hw)) paste0(" for ", hw) else " flagged by the extraction",
      "; values from this solution (", pct, "% of the variance; primary loadings from ", a1, " to ", a2,
      ") should not be interpreted as estimates.") else paste0(
      "A análise fatorial exploratória (", ex_pt, " sobre correlações ", if (ordinal) "policóricas" else "de Pearson",
      if (rot != "none") paste0(", rotação ", rot) else "", ") com ", n_fatores, if (n_fatores == 1) " fator" else " fatores",
      " resultou em solução imprópria, com caso Heywood (carga ou comunalidade igual ou acima de 1)",
      if (nzchar(hw)) paste0(" em ", hw) else " indicado pela extração", "; os valores dessa solução (", pct,
      "% da variância; cargas principais de ", a1, " a ", a2, ") não devem ser interpretados como estimativas.")
  } else if (idioma == "en") {
    frase <- paste0("An exploratory factor analysis (", c(minres = "minimum residuals extraction", wlsmv = paste(if (ordinal) "WLSMV" else "MLR", "estimation"))[[extracao]],
                    " on ", if (ordinal) "polychoric" else "Pearson", " correlations", if (rot != "none") paste0(", ", rot, " rotation") else "",
                    ") with ", n_fatores, " factor", if (n_fatores != 1) "s" else "", " accounted for ", .pr_num(100 * sum(ss) / p, 1, TRUE, decimal),
                    "% of the variance; primary loadings ranged from ", .pr_num(faixa[1], 2, FALSE, decimal), " to ", .pr_num(faixa[2], 2, FALSE, decimal),
                    if (cruz) paste0(", with ", cruz, if (cruz == 1) " cross-loading item" else " cross-loading items") else ", with no cross-loadings",
                    if (n_fatores == 2) paste0("; the factor correlation was ", .pr_num(Phi[1, 2], 2, FALSE, decimal)) else
                      if (n_fatores > 2) paste0("; factor correlations ranged from ", .pr_num(min(Phi[upper.tri(Phi)]), 2, FALSE, decimal), " to ",
                                                .pr_num(max(Phi[upper.tri(Phi)]), 2, FALSE, decimal)) else "", ".")
  } else {
    frase <- paste0("A análise fatorial exploratória (", ex_pt, " sobre correlações ", if (ordinal) "policóricas" else "de Pearson",
                    if (rot != "none") paste0(", rotação ", rot) else "", ") com ", n_fatores, if (n_fatores == 1) " fator" else " fatores",
                    " explicou ", .pr_num(100 * sum(ss) / p, 1, TRUE, decimal), "% da variância; as cargas principais variaram de ",
                    .pr_num(faixa[1], 2, FALSE, decimal), " a ", .pr_num(faixa[2], 2, FALSE, decimal),
                    if (cruz) paste0(", com ", cruz, if (cruz == 1) " item com carga cruzada" else " itens com carga cruzada") else ", sem cargas cruzadas",
                    if (n_fatores == 2) paste0("; a correlação entre os fatores foi ", .pr_num(Phi[1, 2], 2, FALSE, decimal)) else
                      if (n_fatores > 2) paste0("; as correlações entre fatores variaram de ", .pr_num(min(Phi[upper.tri(Phi)]), 2, FALSE, decimal),
                                                " a ", .pr_num(max(Phi[upper.tri(Phi)]), 2, FALSE, decimal)) else "", ".")
  }
  if (extracao == "wlsmv") {
    t1 <- ajuste
    frase <- paste0(sub("\\.$", "", frase), if (idioma == "en") "; fit: " else "; ajuste: ",
                    "CFI = ", .pr_num(t1$valor[1], 3, FALSE, decimal), sep, "TLI = ", .pr_num(t1$valor[2], 3, FALSE, decimal), sep,
                    "RMSEA = ", .pr_num(t1$valor[3], 3, FALSE, decimal), sep, "SRMR = ", .pr_num(t1$valor[4], 3, FALSE, decimal), ".")
  }
  .ee_resultado(paste0("Análise fatorial exploratória (", n_fatores, if (n_fatores == 1) " fator)" else " fatores)"),
                n = n, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef, pressupostos = adeq$pressupostos,
                decisoes = c(adeq$decisoes[1], decisoes), avisos = unique(avisos), apa = c(adeq$apa, apa_nf, frase),
                extra = list(cargas = L, correlacoes_fatores = Phi, comunalidades = h2, msa = adeq$tabela,
                             correlacao = R, tipo_correlacao = adeq$tipo_correlacao, modelo = modelo,
                             solucao_impropria = improprio))
}

# ================================================================ precisão: reteste e juízes (v1.4)
# O curso calcula o teste-reteste pela correlação de Pearson entre T1 e T2. A Pearson só
# mede ordem relativa: uma escala em que todo mundo sobe 5 pontos no reteste tem r = 1 e
# não é estável. A literatura de confiabilidade (Petersen, 2024, cap. 3; Martinková &
# Hladká, 2023, cap. 4; Koo & Li, 2016) pede o ICC de concordância absoluta, de duas
# vias, medida única: ICC(A,1) de McGraw e Wong (1996), o ICC(2,1) de Shrout e Fleiss
# (1979). A Pearson sai ao lado, com a mudança sistemática de média que ela não enxerga.
# Erro de medida: EPM = DP × √(1 − ICC), com o DP de todas as observações, e mudança
# mínima detectável MMD95 = 1,96 × √2 × EPM (Weir, 2005); limites de concordância de
# Bland e Altman (1986) = média da diferença ± 1,96 DP da diferença.
.ee_icc_tabela <- function(x, nivel) {
  .pr_exige("psych")
  r <- suppressWarnings(suppressMessages(psych::ICC(x, alpha = 1 - nivel, lmer = FALSE)))$results
  data.frame(tipo = c("ICC(1)", "ICC(A,1)", "ICC(C,1)", "ICC(1,k)", "ICC(A,k)", "ICC(C,k)"),
             mcgraw_wong = c("uma via, medida única", "duas vias, concordância absoluta, medida única",
                             "duas vias, consistência, medida única", "uma via, média", "duas vias, concordância absoluta, média",
                             "duas vias, consistência, média"),
             icc = r$ICC, ic_inf = r[["lower bound"]], ic_sup = r[["upper bound"]], F = r$F,
             gl1 = r$df1, gl2 = r$df2, p = r$p, row.names = NULL)
}

.ee_koo_li <- function(v, idioma = "pt") {
  if (is.na(v)) return(NA_character_)
  pt <- c("fraca", "moderada", "boa", "excelente"); en <- c("poor", "moderate", "good", "excellent")
  (if (identical(idioma, "en")) en else pt)[findInterval(v, c(.5, .75, .9)) + 1]
}

.ee_koo_li_regua <- function(decimal = ".") {
  f <- function(x) .pr_num(x, 2, FALSE, decimal)
  paste0("< ", f(.5), " fraca; ", f(.5), " a ", f(.75), " moderada; ", f(.75), " a ", f(.9), " boa; > ", f(.9), " excelente")
}

confiabilidade_teste_reteste <- function(dados, t1, t2, nivel = .95, decimal = ".", idioma = "pt") {
  dd <- .ee_dados(dados, c(t1, t2), numericas = c(t1, t2))
  d <- dd$d; n <- nrow(d)
  if (n < 5) stop("o reteste pede pelo menos 5 pessoas com as duas medidas", call. = FALSE)
  if (stats::sd(d[[t1]]) == 0 || stats::sd(d[[t2]]) == 0) stop("uma das ocasiões não tem variação", call. = FALSE)
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  icc <- .ee_icc_tabela(d, nivel)
  a1 <- icc[icc$tipo == "ICC(A,1)", ]; c1 <- icc[icc$tipo == "ICC(C,1)", ]
  ct <- stats::cor.test(d[[t1]], d[[t2]], conf.level = nivel)
  dif <- d[[t2]] - d[[t1]]
  tt <- stats::t.test(d[[t2]], d[[t1]], paired = TRUE, conf.level = nivel)
  dp_todas <- stats::sd(c(d[[t1]], d[[t2]]))
  epm <- dp_todas * sqrt(1 - max(a1$icc, 0)); mmd <- 1.96 * sqrt(2) * epm
  ba <- mean(dif) + c(-1.96, 1.96) * stats::sd(dif)
  ef <- data.frame(indice = c("ICC(A,1)", "ICC(C,1)", "r de Pearson"),
                   valor = c(a1$icc, c1$icc, unname(ct$estimate)),
                   ic_inf = c(a1$ic_inf, c1$ic_inf, ct$conf.int[1]), ic_sup = c(a1$ic_sup, c1$ic_sup, ct$conf.int[2]),
                   principal = c(TRUE, FALSE, FALSE))
  erro <- data.frame(media_t1 = mean(d[[t1]]), media_t2 = mean(d[[t2]]), diferenca_media = mean(dif),
                     dif_ic_inf = tt$conf.int[1], dif_ic_sup = tt$conf.int[2], p_mudanca = tt$p.value,
                     epm = epm, mmd95 = mmd, ba_inf = ba[1], ba_sup = ba[2])
  avisos <- character()
  if (abs(c1$icc - a1$icc) >= .05) avisos <- c(avisos, paste0("ICC de consistência e de concordância diferem (", .pr_num(c1$icc, 2, FALSE, decimal), " contra ",
                                                          .pr_num(a1$icc, 2, FALSE, decimal), "): há mudança sistemática entre as ocasiões (efeito de prática, de reteste ou mudança real). A Pearson não enxerga isso."))
  if (tt$p.value < .05) avisos <- c(avisos, paste0("A média mudou entre as ocasiões (diferença ", .pr_num(mean(dif), 2, TRUE, decimal), sep, .pr_p(tt$p.value, decimal),
                                                   "): o ICC(A,1) já desconta essa mudança, a Pearson não. Discutir efeito de prática, de reteste ou mudança real no intervalo."))
  if (n < 30) avisos <- c(avisos, paste0("n = ", n, ": o IC do ICC fica largo; Koo e Li (2016) pedem pelo menos 30 participantes."))
  faixa <- paste0(.ee_koo_li(a1$ic_inf, idioma), " a ", .ee_koo_li(a1$ic_sup, idioma))
  decisoes <- c("ICC(A,1), duas vias, concordância absoluta, medida única (McGraw & Wong, 1996), como índice principal; o curso usa a Pearson, que sai ao lado.",
                paste0("Leitura pela régua de Koo e Li (2016) aplicada ao IC, como eles recomendam: ", .ee_koo_li_regua(decimal), ". Aqui: ", faixa, "."),
                "Intervalo entre as ocasiões e o que se espera do construto nesse intervalo (traço ou estado) entram no relato: sem isso o ICC não se interpreta.")
  apa <- if (en) paste0("Test-retest reliability was assessed with a two-way, absolute-agreement, single-measure intraclass correlation (McGraw & Wong, 1996; N = ", n, "): ",
                        .ee_efeito_txt("ICC(A,1)", a1$icc, a1$ic_inf, a1$ic_sup, TRUE, 2, decimal, "en", nivel), ", ", faixa, " by Koo and Li's (2016) guidelines (Pearson's ",
                        .ee_efeito_txt("r", unname(ct$estimate), ct$conf.int[1], ct$conf.int[2], TRUE, 2, decimal, "en", nivel), "). The mean change between occasions was ",
                        .pr_num(mean(dif), 2, TRUE, decimal), " (", .ee_ic(tt$conf.int[1], tt$conf.int[2], 2, TRUE, decimal, nivel, "en"), "; ", .pr_p(tt$p.value, decimal), "); SEM = ",
                        .pr_num(epm, 2, TRUE, decimal), " and MDC95 = ", .pr_num(mmd, 2, TRUE, decimal), ".") else
    paste0("A estabilidade teste-reteste foi avaliada pelo coeficiente de correlação intraclasse de duas vias, concordância absoluta e medida única (McGraw & Wong, 1996; N = ", n, "): ",
           .ee_efeito_txt("ICC(A,1)", a1$icc, a1$ic_inf, a1$ic_sup, TRUE, 2, decimal, "pt", nivel), ", concordância ", faixa, " pela régua de Koo e Li (2016) (r de Pearson: ",
           .ee_efeito_txt("r", unname(ct$estimate), ct$conf.int[1], ct$conf.int[2], TRUE, 2, decimal, "pt", nivel), "). A mudança média entre as ocasiões foi de ",
           .pr_num(mean(dif), 2, TRUE, decimal), sep, .ee_ic(tt$conf.int[1], tt$conf.int[2], 2, TRUE, decimal, nivel, "pt"), sep, .pr_p(tt$p.value, decimal),
           "; EPM = ", .pr_num(epm, 2, TRUE, decimal), " e mudança mínima detectável (95%) = ", .pr_num(mmd, 2, TRUE, decimal), ".")
  .ee_resultado("Confiabilidade teste-reteste", n = n, n_excluidos = dd$n_excluidos, tabela = icc, efeito = ef,
                decisoes = decisoes, avisos = avisos, apa = apa, extra = list(tabelas_extra = list("Mudança e erro de medida" = erro)))
}

# Concordância entre juízes ou avaliadores (pessoas ou itens nas linhas, avaliadores nas
# colunas). Nominal: kappa de Cohen com 2 avaliadores e de Fleiss com 3 ou mais; ordinal:
# kappa ponderado (pesos quadráticos) com 2 avaliadores e Fleiss com 3 ou mais (que trata
# as categorias como nominais, dito na saída); intervalar: ICC(A,1) e ICC(A,k). Régua de
# Landis e Koch (1977) para o kappa, declarada como régua. O curso usa o kappa na triagem
# de revisões e o kappam.fleiss e o ICC na análise de juízes de desejabilidade social.
.ee_fleiss <- function(m) {
  cats <- sort(unique(as.vector(m[!is.na(m)])))
  r <- ncol(m)
  nij <- t(apply(m, 1, function(l) table(factor(l, levels = cats))))
  if (length(cats) == 1) nij <- t(nij)
  Pi <- (rowSums(nij^2) - r) / (r * (r - 1))
  pj <- colSums(nij) / (nrow(m) * r)
  Pe <- sum(pj^2)
  (mean(Pi) - Pe) / (1 - Pe)
}

# AC1 (nominal) e AC2 (pesos quadráticos) de Gwet (2008) para r avaliadores com dados
# completos: pa pelo acordo ponderado entre pares e pe = Tw/(q(q - 1)) Σ πk(1 - πk), que não
# despenca quando quase todos os casos caem numa categoria (o paradoxo do kappa:
# Feinstein & Cicchetti, 1990). `cats` fixa as categorias possíveis (os níveis do fator):
# com todos os casos em "sim", q = 2 e AC1 = 1, enquanto o kappa nem é calculável.
.ee_gwet <- function(m, cats, ponderado = FALSE) {
  q <- length(cats); r <- ncol(m)
  if (q < 2) return(NA_real_)
  w <- if (ponderado) { v <- suppressWarnings(as.numeric(cats)); if (anyNA(v)) v <- seq_len(q)
    1 - outer(v, v, "-")^2 / (max(v) - min(v))^2 } else diag(q)
  rik <- t(apply(m, 1, function(l) table(factor(l, levels = cats))))
  if (q == 1) rik <- t(rik)
  rstar <- rik %*% t(w)
  pa <- mean(rowSums(rik * (rstar - 1)) / (r * (r - 1)))
  pik <- colMeans(rik / r)
  pe <- sum(w) / (q * (q - 1)) * sum(pik * (1 - pik))
  (pa - pe) / (1 - pe)
}

# Alfa de Krippendorff (Hayes & Krippendorff, 2007) com dados completos: matriz de
# coincidências entre pares de avaliadores do mesmo caso; métricas nominal, ordinal
# (pelos postos acumulados) e intervalar.
.ee_krippendorff <- function(m, cats, metrica = c("nominal", "ordinal", "intervalar")) {
  metrica <- match.arg(metrica)
  r <- ncol(m); q <- length(cats)
  if (q < 2) return(NA_real_)
  idx <- matrix(match(as.character(m), as.character(cats)), nrow(m))
  o <- matrix(0, q, q)
  for (i in seq_len(nrow(idx))) {
    v <- idx[i, ]
    tab <- tabulate(v, q)
    o <- o + (outer(tab, tab) - diag(tab)) / (r - 1)
  }
  nc <- rowSums(o); N <- sum(nc)
  d2 <- switch(metrica,
    nominal = 1 - diag(q),
    intervalar = { v <- as.numeric(cats); outer(v, v, "-")^2 },
    ordinal = { cs <- cumsum(nc); outer(seq_len(q), seq_len(q), Vectorize(function(a, b) {
      lo <- min(a, b); hi <- max(a, b); (sum(nc[lo:hi]) - (nc[a] + nc[b]) / 2)^2 })) })
  Do <- sum(o * d2) / N
  De <- sum(outer(nc, nc) * d2) / (N * (N - 1))
  if (De == 0) return(NA_real_)
  1 - Do / De
}

.ee_landis_koch <- function(k, idioma = "pt") {
  if (is.na(k)) return(NA_character_)
  pt <- c("pobre", "leve", "razoável", "moderada", "substancial", "quase perfeita")
  en <- c("poor", "slight", "fair", "moderate", "substantial", "almost perfect")
  (if (identical(idioma, "en")) en else pt)[findInterval(k, c(0, .2, .4, .6, .8), left.open = TRUE) + 1]
}

concordancia_juizes <- function(dados, avaliadores, tipo = c("nominal", "ordinal", "intervalar"), reamostras = 2000,
                                semente = 2026, nivel = .95, categorias = NULL, decimal = ".", idioma = "pt") {
  tipo <- match.arg(tipo)
  if (length(avaliadores) < 2) stop("concordância pede pelo menos 2 avaliadores", call. = FALSE)
  dd <- .ee_dados(dados, avaliadores)
  d <- dd$d; n <- nrow(d); r <- length(avaliadores)
  if (n < 2) stop("poucos casos com todas as avaliações", call. = FALSE)
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  # categorias possíveis: argumento, níveis dos fatores ou valores observados
  if (is.null(categorias)) {
    niv <- unique(unlist(lapply(d, function(x) if (is.factor(x)) levels(x))))
    obs <- sort(unique(unlist(lapply(d, as.character))))
    categorias <- if (length(niv)) union(niv, obs) else {
      u <- unique(unlist(lapply(d, function(x) if (is.numeric(x)) x else as.character(x))))
      if (is.numeric(u)) sort(u) else sort(u) }
  }
  categorias <- as.character(categorias)
  m <- as.matrix(as.data.frame(lapply(d, as.character), stringsAsFactors = FALSE))
  fora <- setdiff(unique(as.vector(m)), categorias)
  if (length(fora)) stop("avaliação fora das categorias declaradas: ", paste(fora, collapse = ", "), call. = FALSE)
  acordo <- mean(apply(m, 1, function(l) length(unique(l)) == 1))
  decisoes <- character(); avisos <- character()
  rest <- .ee_construcao_semente(semente); on.exit(rest(), add = TRUE)
  bs_ic <- function(fun) {
    b <- replicate(reamostras, fun(m[sample.int(n, n, replace = TRUE), , drop = FALSE]))
    stats::quantile(b, c((1 - nivel) / 2, 1 - (1 - nivel) / 2), na.rm = TRUE, names = FALSE)
  }
  metr <- switch(tipo, nominal = "nominal", ordinal = "ordinal", intervalar = "intervalar")
  ka_fun <- function(x) .ee_krippendorff(x, categorias, metr)
  ka <- ka_fun(m); ka_ic <- bs_ic(ka_fun)
  ka_lin <- data.frame(indice = paste0("alfa de Krippendorff (", metr, ")"), valor = ka, ic_inf = ka_ic[1], ic_sup = ka_ic[2], principal = FALSE)
  if (tipo == "intervalar") {
    for (v in avaliadores) if (!is.numeric(d[[v]])) stop("com tipo intervalar as avaliações precisam ser numéricas", call. = FALSE)
    icc <- .ee_icc_tabela(d, nivel)
    a1 <- icc[icc$tipo == "ICC(A,1)", ]; ak <- icc[icc$tipo == "ICC(A,k)", ]
    ef <- rbind(data.frame(indice = c("ICC(A,1)", "ICC(A,k)"), valor = c(a1$icc, ak$icc), ic_inf = c(a1$ic_inf, ak$ic_inf),
                           ic_sup = c(a1$ic_sup, ak$ic_sup), principal = c(TRUE, FALSE)), ka_lin)
    tab <- icc
    decisoes <- c("ICC de duas vias e concordância absoluta: ICC(A,1) descreve um avaliador; ICC(A,k), a média dos k avaliadores (usar este se o escore final for a média).",
                  paste0("Régua de Koo e Li (2016) aplicada ao IC (", .ee_koo_li_regua(decimal), "): ", .ee_koo_li(a1$ic_inf, idioma), " a ",
                         .ee_koo_li(a1$ic_sup, idioma), "."))
    apa <- if (en) paste0("Inter-rater agreement among ", r, " raters on ", n, " targets: ", .ee_efeito_txt("ICC(A,1)", a1$icc, a1$ic_inf, a1$ic_sup, TRUE, 2, decimal, "en", nivel),
                          "; ", .ee_efeito_txt("ICC(A,k)", ak$icc, ak$ic_inf, ak$ic_sup, TRUE, 2, decimal, "en", nivel), ".") else
      paste0("A concordância entre ", r, " avaliadores em ", n, " casos foi ", .ee_efeito_txt("ICC(A,1)", a1$icc, a1$ic_inf, a1$ic_sup, TRUE, 2, decimal, "pt", nivel),
             sep, .ee_efeito_txt("ICC(A,k)", ak$icc, ak$ic_inf, ak$ic_sup, TRUE, 2, decimal, "pt", nivel), ".")
  } else {
    ponderado <- tipo == "ordinal"
    g_fun <- function(x) .ee_gwet(x, categorias, ponderado)
    g <- g_fun(m); g_ic <- bs_ic(g_fun)
    g_nome <- if (ponderado) "AC2 de Gwet (pesos quadráticos)" else "AC1 de Gwet"
    g_lin <- data.frame(indice = g_nome, valor = g, ic_inf = g_ic[1], ic_sup = g_ic[2], principal = FALSE)
    if (r == 2) {
      lin <- if (ponderado) "weighted kappa" else "unweighted kappa"
      k <- li <- ls <- NA_real_
      if (length(unique(as.vector(m))) > 1) {
        .pr_exige("psych")
        mc <- matrix(match(m, categorias), nrow(m))
        ck <- suppressWarnings(psych::cohen.kappa(mc, alpha = 1 - nivel, levels = seq_along(categorias)))
        k <- ck$confid[lin, "estimate"]; li <- ck$confid[lin, "lower"]; ls <- ck$confid[lin, "upper"]
      }
      ef <- rbind(data.frame(indice = c(if (ponderado) "kappa ponderado (quadrático)" else "kappa de Cohen", "concordância simples"),
                             valor = c(k, acordo), ic_inf = c(li, NA), ic_sup = c(ls, NA), principal = c(TRUE, FALSE)), g_lin, ka_lin)
      tab <- data.frame(kappa = k, concordancia = acordo, gwet = g, krippendorff = ka, n = n)
      decisoes <- c(paste0(if (ponderado) "Kappa ponderado de Cohen (1968), pesos quadráticos (psych::cohen.kappa)" else "Kappa de Cohen (1960)", "; IC assintótico."),
                    paste0("Régua de Landis e Koch (1977): ", .ee_landis_koch(k, idioma), "."))
      simb <- if (ponderado) "κw" else "κ"
    } else {
      k <- .ee_fleiss(m); kic <- bs_ic(.ee_fleiss); li <- kic[1]; ls <- kic[2]
      ef <- rbind(data.frame(indice = c("kappa de Fleiss", "concordância simples (todos iguais)"), valor = c(k, acordo),
                             ic_inf = c(li, NA), ic_sup = c(ls, NA), principal = c(TRUE, FALSE)), g_lin, ka_lin)
      tab <- data.frame(kappa_fleiss = k, concordancia = acordo, gwet = g, krippendorff = ka, avaliadores = r, n = n)
      decisoes <- c(paste0("Kappa de Fleiss (1971) para ", r, " avaliadores; IC por bootstrap percentílico dos casos."),
                    paste0("Régua de Landis e Koch (1977): ", .ee_landis_koch(k, idioma), "."))
      if (ponderado) avisos <- c(avisos, "O kappa de Fleiss trata as categorias como nominais: discordância de um ponto pesa como discordância extrema. Com categorias ordenadas, ler o AC2 de Gwet e o alfa de Krippendorff ordinal da tabela, ou o ICC (tipo = \"intervalar\").")
      simb <- "κ de Fleiss"
    }
    paradoxo <- acordo >= .80 && (is.na(k) || k < .40)
    if (paradoxo) avisos <- c(avisos, paste0("Paradoxo do kappa (Feinstein & Cicchetti, 1990; aula de limitações do kappa): concordância simples de ",
      .pr_num(100 * acordo, 1, TRUE, decimal), "% com κ ", if (is.na(k)) "incalculável (todos os casos na mesma categoria: a concordância esperada pelo acaso vale 1)" else paste0("= ", .pr_num(k, 2, FALSE, decimal)),
      ". Com prevalência muito desigual entre categorias, relatar o ", g_nome, " e a concordância simples ao lado do κ."))
    decisoes <- c(decisoes, paste0(g_nome, " (Gwet, 2008) e alfa de Krippendorff ", metr, " (Hayes & Krippendorff, 2007) ao lado, com IC por bootstrap percentílico dos casos (",
                                   reamostras, " reamostras, semente ", semente, "); categorias consideradas: ", paste(categorias, collapse = ", "), "."))
    k_txt <- if (is.na(k)) paste0(simb, if (en) " not computable" else " incalculável") else .ee_efeito_txt(simb, k, li, ls, TRUE, 2, decimal, idioma, nivel)
    g_simb <- if (ponderado) "AC2" else "AC1"
    apa <- if (en) paste0("Agreement ", if (r == 2) "between the two raters" else paste0("among ", r, " raters"), " on ", n, " cases: ", k_txt,
                          if (!is.na(k)) paste0(", ", .ee_landis_koch(k, "en"), " by Landis and Koch (1977)") else "",
                          "; Gwet's ", .ee_efeito_txt(g_simb, g, g_ic[1], g_ic[2], TRUE, 2, decimal, "en", nivel), "; raw agreement ", .pr_num(100 * acordo, 1, TRUE, decimal), "%.") else
      paste0("A concordância ", if (r == 2) "entre os dois avaliadores" else paste0("entre ", r, " avaliadores"), " em ", n, " casos foi ", k_txt,
             if (!is.na(k)) paste0(", ", .ee_landis_koch(k), " pela régua de Landis e Koch (1977)") else "",
             sep, .ee_efeito_txt(paste0(g_simb, " de Gwet"), g, g_ic[1], g_ic[2], TRUE, 2, decimal, "pt", nivel), sep, "concordância simples de ", .pr_num(100 * acordo, 1, TRUE, decimal), "%.")
  }
  .ee_resultado(paste0("Concordância entre avaliadores (", tipo, ")"), n = n, n_excluidos = dd$n_excluidos, tabela = tab, efeito = ef,
                decisoes = decisoes, avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- análise clássica de itens (v1.4)
# Martinková e Hladká (2023, cap. 5): dificuldade (dicotômico: proporção de acerto ou
# endosso; ordinal: (média - mínimo)/(máximo - mínimo), a dificuldade corrigida),
# discriminação por correlação item-total (RIT) e item-resto (RIR) e pelo índice
# superior-inferior (ULI): diferença da dificuldade entre o terço superior e o inferior do
# escore total, com os grupos pelos quantis do total (como ShinyItemAnalysis::gDiscrim;
# empates no limite mudam o tamanho dos grupos). Régua de Ebel (1954), pelo livro: ULI por
# terços >= .20, salvo item muito fácil ou muito difícil. Alfa sem o item ao lado. RIT e
# RIR conferidos com psych::alpha (raw.r e r.drop) nos evals.
analise_classica_itens <- function(dados, itens, minimo = NULL, maximo = NULL, grupos = 3, decimal = ".", idioma = "pt") {
  dd <- .ee_dados(dados, itens, numericas = itens)
  d <- dd$d; n <- nrow(d); k <- length(itens)
  if (k < 3) stop("análise de itens pede pelo menos 3 itens")
  mn <- if (is.null(minimo)) vapply(d, min, 0) else rep(minimo, k)
  mx <- if (is.null(maximo)) vapply(d, max, 0) else rep(maximo, k)
  if (any(mx <= mn)) stop("item sem variação ou máximo <= mínimo: ", paste(itens[mx <= mn], collapse = ", "))
  tot <- rowSums(d)
  q <- stats::quantile(tot, probs = seq(0, 1, length.out = grupos + 1))
  gr <- cut(tot, unique(q), include.lowest = TRUE, labels = FALSE)
  inf <- gr == 1; sup <- gr == max(gr)
  dific <- (colMeans(d) - mn) / (mx - mn)
  rit <- vapply(itens, function(i) stats::cor(d[[i]], tot), 0)
  rir <- vapply(itens, function(i) stats::cor(d[[i]], tot - d[[i]]), 0)
  uli <- vapply(seq_len(k), function(j) (mean(d[[j]][sup]) - mean(d[[j]][inf])) / (mx[j] - mn[j]), 0)
  alfa <- function(m) { kk <- ncol(m); (kk / (kk - 1)) * (1 - sum(apply(m, 2, stats::var)) / stats::var(rowSums(m))) }
  a_tot <- alfa(as.matrix(d))
  a_sem <- vapply(seq_len(k), function(j) alfa(as.matrix(d[, -j, drop = FALSE])), 0)
  tab <- data.frame(item = itens, media = colMeans(d), dp = vapply(d, stats::sd, 0), dificuldade = unname(dific),
                    rit = unname(rit), rir = unname(rir), uli = uli, alfa_sem_item = a_sem, row.names = NULL)
  avisos <- character()
  fracos <- tab$item[tab$uli < .20 & tab$dificuldade > .20 & tab$dificuldade < .80]
  if (length(fracos)) avisos <- c(avisos, paste0("ULI < ", .pr_num(.20, 2, FALSE, decimal), " em item de dificuldade intermediária (Ebel, 1954, pela régua de Martinková & Hladká, 2023): ",
                                                 paste(fracos, collapse = ", "), ". Rever conteúdo antes de excluir."))
  neg <- tab$item[tab$rir < 0]
  if (length(neg)) avisos <- c(avisos, paste0("RIR negativo em ", paste(neg, collapse = ", "), ": item invertido sem recodificar ou mal escrito."))
  sobe <- tab$item[tab$alfa_sem_item > a_tot + .01]
  if (length(sobe)) decisoes_extra <- paste0("Alfa sobe sem: ", paste(sobe, collapse = ", "), " (decidir pelo conteúdo, não pelo alfa).") else decisoes_extra <- character()
  apa <- if (identical(idioma, "en"))
    paste0("Item difficulty ranged from ", .pr_num(min(dific), 2, FALSE, decimal), " to ", .pr_num(max(dific), 2, FALSE, decimal), ", and corrected item-total correlations (RIR) from ",
           .pr_num(min(rir), 2, FALSE, decimal), " to ", .pr_num(max(rir), 2, FALSE, decimal), " (N = ", n, "; α = ", .pr_num(a_tot, 2, FALSE, decimal), ").") else
    paste0("A dificuldade dos itens variou de ", .pr_num(min(dific), 2, FALSE, decimal), " a ", .pr_num(max(dific), 2, FALSE, decimal), ", e as correlações item-resto (RIR) de ",
           .pr_num(min(rir), 2, FALSE, decimal), " a ", .pr_num(max(rir), 2, FALSE, decimal), " (N = ", n, "; α = ", .pr_num(a_tot, 2, FALSE, decimal), ").")
  .ee_resultado("Análise clássica de itens", n = n, n_excluidos = dd$n_excluidos, tabela = tab,
                efeito = data.frame(indice = "alfa de Cronbach (todos os itens)", valor = a_tot, ic_inf = NA_real_, ic_sup = NA_real_, principal = TRUE),
                decisoes = c(paste0("Dificuldade corrigida (média - mínimo)/(máximo - mínimo), com mínimo e máximo ", if (is.null(minimo)) "observados" else "informados",
                                    "; ULI entre o grupo ", grupos, " e o grupo 1 do escore total (quantis)."),
                             "Discriminação pela TCT depende da amostra; com estrutura multidimensional, calcular por subescala.", decisoes_extra),
                avisos = avisos, apa = apa)
}

# ---------------------------------------------------------------- confiabilidade em EMA (v1.4)
# Delineamento pessoa x tempo x item (EMA, diário), pela teoria da generalizabilidade de
# Shrout e Lane (2012), como psych::multilevel.reliability, que esta função envolve:
#   RkF  confiabilidade entre pessoas da média de todos os tempos (tempos fixos);
#   R1R  de uma medida num tempo sorteado; RkR, da média de k tempos sorteados;
#   Rc   confiabilidade da mudança dentro da pessoa (variância pessoa x tempo contra o
#        resíduo dividido pelo número de itens; Cranford et al., 2006).
# Dados em formato largo: uma linha por pessoa e tempo, uma coluna por item. Com dados
# desbalanceados (falta de resposta em EMA), os componentes saem do lme4 (lmer = TRUE).
confiabilidade_ema <- function(dados, id, tempo, itens, decimal = ".", idioma = "pt") {
  .pr_exige("psych")
  d <- as.data.frame(dados)
  for (v in c(id, tempo, itens)) if (!v %in% names(d)) stop("coluna ausente: ", v)
  d <- d[!is.na(d[[id]]) & !is.na(d[[tempo]]), c(id, tempo, itens)]
  completo <- stats::complete.cases(d[, itens])
  balanceado <- all(completo) && length(unique(table(d[[id]]))) == 1
  if (!balanceado) .pr_exige("lme4")
  r <- suppressWarnings(suppressMessages(psych::multilevel.reliability(d, grp = id, Time = tempo, items = itens, lme = FALSE,
                                                                       lmer = !balanceado, plot = FALSE, alpha = FALSE)))
  comp <- r$components
  comp <- data.frame(componente = rownames(comp), variancia = comp[, 1], proporcao = comp[, 2], row.names = NULL)
  comp <- comp[stats::complete.cases(comp), ]
  coef <- data.frame(indice = c("RkF", "R1R", "RkR", "Rc"),
                     descricao = c("entre pessoas, média de todos os tempos (tempos fixos)", "entre pessoas, uma medida num tempo sorteado",
                                   "entre pessoas, média de k tempos sorteados", "mudança dentro da pessoa"),
                     valor = c(r$RkF, r$R1R, r$RkR, r$Rc), stringsAsFactors = FALSE)
  avisos <- character()
  if (isTRUE(r$Rc < .70)) avisos <- c(avisos, paste0("Rc = ", .pr_num(r$Rc, 2, FALSE, decimal), ": a escala detecta mal a variação dentro da pessoa; análise intraindividual (redes, mudança) fica atenuada."))
  apa <- if (identical(idioma, "en"))
    paste0("Generalizability analysis of the person × time × item design (Shrout & Lane, 2012; ", r$n.obs, " persons, ", r$n.time, " occasions, ", r$n.items,
           " items) yielded between-person reliability RkF = ", .pr_num(r$RkF, 2, FALSE, decimal), " and within-person reliability of change Rc = ", .pr_num(r$Rc, 2, FALSE, decimal), ".") else
    paste0("A análise de generalizabilidade do delineamento pessoa × tempo × item (Shrout & Lane, 2012; ", r$n.obs, " pessoas, ", r$n.time, " ocasiões, ", r$n.items,
           " itens) estimou confiabilidade entre pessoas RkF = ", .pr_num(r$RkF, 2, FALSE, decimal), " e confiabilidade da mudança dentro da pessoa Rc = ",
           .pr_num(r$Rc, 2, FALSE, decimal), ".")
  .ee_resultado("Confiabilidade em medidas intensivas (EMA)", n = r$n.obs, tabela = coef,
                efeito = data.frame(indice = c("RkF", "Rc"), valor = c(r$RkF, r$Rc), ic_inf = NA_real_, ic_sup = NA_real_, principal = c(TRUE, TRUE)),
                decisoes = c(paste0("Componentes de variância por ", if (balanceado) "ANOVA (dados balanceados)" else "lme4 (dados desbalanceados)",
                                    ", via psych::multilevel.reliability (Shrout & Lane, 2012)."),
                             "Rc é a confiabilidade relevante para análise intraindividual (Cranford et al., 2006); RkF, para comparar pessoas pela média."),
                avisos = avisos, apa = apa, extra = list(tabelas_extra = list(`Componentes de variância` = comp)))
}

# ---------------------------------------------------------------- v1.5: duas metades
# confiabilidade_metades(): precisão por duas metades, a da aula "Fidedignidade Split-Half"
# do curso de análise fatorial. Para uma divisão (ímpares contra pares, ou primeira contra
# segunda metade): correlação entre as metades, Spearman-Brown 2r/(1 + r) e Guttman
# (Flanagan-Rulon) λ4 = 2[1 - (s²A + s²B)/s²X], que não supõe metades com a mesma
# variância. Com metodo = "todas", λ4 de todas as divisões ao meio (média, mínimo e
# máximo): com número par de itens, a média é o alfa (Cronbach, 1951); o máximo superestima
# a precisão com muitos itens e amostra pequena.
confiabilidade_metades <- function(dados, itens, metodo = c("impar_par", "metades", "todas"), invertidos = NULL, minimo = NULL, maximo = NULL,
                                   decimal = ".", idioma = "pt") {
  metodo <- match.arg(metodo)
  dd <- .ee_dados(dados, itens, numericas = itens); d <- dd$d; k <- length(itens)
  if (k < 4) stop("duas metades pedem pelo menos 4 itens")
  if (length(invertidos)) {
    if (is.null(minimo) || is.null(maximo)) stop("para inverter itens, informe minimo e maximo da escala")
    for (v in invertidos) d[[v]] <- minimo + maximo - d[[v]]
  }
  x <- rowSums(d); vx <- stats::var(x)
  if (metodo == "todas") {
    # Todas as divisões em duas metades (k par: metades iguais; k ímpar: metades de
    # floor(k/2) e ceiling(k/2)); com mais de 16 itens, 10.000 divisões sorteadas.
    h <- floor(k / 2); S <- stats::cov(d)
    divs <- if (k <= 16) utils::combn(k, h, simplify = FALSE) else { set.seed(2026); lapply(seq_len(10000), function(i) sort(sample(k, h))) }
    if (k %% 2 == 0 && k <= 16) divs <- divs[vapply(divs, function(a) a[1] == 1, TRUE)]
    l4 <- vapply(divs, function(a) { b <- setdiff(seq_len(k), a); 2 * (1 - (sum(S[a, a]) + sum(S[b, b])) / sum(S)) }, 0)
    alfa <- k / (k - 1) * (1 - sum(diag(S)) / sum(S))
    tab <- data.frame(divisao = c("média das divisões", "mínima", "máxima"), lambda4 = c(mean(l4), min(l4), max(l4)), alfa = alfa,
                      n_divisoes = length(l4), row.names = NULL)
    efeito <- data.frame(indice = "λ4 médio", valor = mean(l4), ic_inf = NA_real_, ic_sup = NA_real_, principal = TRUE)
    frase <- paste0("média = ", .pr_num(mean(l4), 2, FALSE, decimal), ", mínima = ", .pr_num(min(l4), 2, FALSE, decimal), ", máxima = ", .pr_num(max(l4), 2, FALSE, decimal))
  } else {
    a <- if (metodo == "impar_par") itens[seq(1, k, 2)] else itens[seq_len(floor(k / 2))]
    b <- setdiff(itens, a)
    xa <- rowSums(d[a]); xb <- rowSums(d[b]); r <- stats::cor(xa, xb)
    sb <- 2 * r / (1 + r); l4 <- 2 * (1 - (stats::var(xa) + stats::var(xb)) / vx)
    tab <- data.frame(metade_a = paste(a, collapse = " "), metade_b = paste(b, collapse = " "), r_metades = r, spearman_brown = sb, guttman_lambda4 = l4,
                      razao_dp_metades = stats::sd(xa) / stats::sd(xb), row.names = NULL)
    efeito <- data.frame(indice = c("Spearman-Brown", "Guttman λ4"), valor = c(sb, l4), ic_inf = NA_real_, ic_sup = NA_real_, principal = c(FALSE, TRUE))
    frase <- paste0("r = ", .pr_num(r, 2, FALSE, decimal), ", Spearman-Brown = ", .pr_num(sb, 2, FALSE, decimal), ", Guttman λ4 = ", .pr_num(l4, 2, FALSE, decimal))
  }
  en <- identical(idioma, "en")
  apa <- if (en) paste0("Split-half reliability (", switch(metodo, impar_par = "odd-even split", metades = "first-second half", todas = "all splits"), "; n = ", nrow(d), "): ", frase, ".") else
    paste0("A precisão por duas metades (", switch(metodo, impar_par = "itens ímpares contra pares", metades = "primeira contra segunda metade", todas = "todas as divisões"),
           "; n = ", nrow(d), "): ", frase, ".")
  .ee_resultado(paste0("Precisão por duas metades (", metodo, ")"), n = nrow(d), n_excluidos = dd$n_excluidos, tabela = tab, efeito = efeito,
                decisoes = c("Guttman λ4 como índice principal: Spearman-Brown supõe metades paralelas e superestima quando as variâncias diferem.",
                             if (length(invertidos)) paste0("Itens invertidos antes da soma: ", paste(invertidos, collapse = ", "), ".") else character(),
                             "Duas metades depende da divisão escolhida; o ômega da AFC segue como estimativa principal de consistência interna."),
                apa = apa)
}

# ---------------------------------------------------------------- v1.6: profecia, correções e teoria G
# Profecia de Spearman-Brown (Martinková & Hladká, 2023, 4.2.5.1): confiabilidade de um
# teste com k_novo itens a partir de um com k_atual itens e confiabilidade rxx; ou, com
# alvo, o número de itens paralelos que leva a ela (arredondado para cima). Supõe itens
# novos paralelos (ou essencialmente tau-equivalentes) aos atuais.
spearman_brown <- function(confiabilidade, k_atual, k_novo = NULL, alvo = NULL, decimal = ".", idioma = "pt") {
  if (confiabilidade <= 0 || confiabilidade >= 1) stop("confiabilidade precisa estar entre 0 e 1", call. = FALSE)
  if (is.null(k_novo) == is.null(alvo)) stop("informe k_novo (confiabilidade projetada) ou alvo (itens necessários), um dos dois", call. = FALSE)
  en <- identical(idioma, "en")
  if (!is.null(k_novo)) {
    m <- k_novo / k_atual; rn <- m * confiabilidade / (1 + (m - 1) * confiabilidade)
    tab <- data.frame(k_atual = k_atual, confiabilidade_atual = confiabilidade, k_novo = k_novo, razao = m, confiabilidade_projetada = rn)
    apa <- if (en) paste0("By the Spearman-Brown prophecy formula, a ", k_novo, "-item version of the ", k_atual, "-item test (reliability ",
                          .pr_num(confiabilidade, 2, FALSE, decimal), ") would reach ", .pr_num(rn, 2, FALSE, decimal), ", if the added items are parallel to the current ones.") else
      paste0("Pela profecia de Spearman-Brown, a versão com ", k_novo, " itens do teste de ", k_atual, " itens (confiabilidade ", .pr_num(confiabilidade, 2, FALSE, decimal),
             ") chegaria a ", .pr_num(rn, 2, FALSE, decimal), ", se os itens acrescentados forem paralelos aos atuais.")
    ef <- data.frame(indice = "confiabilidade projetada", valor = rn, ic_inf = NA, ic_sup = NA, principal = TRUE)
  } else {
    if (alvo <= 0 || alvo >= 1) stop("alvo precisa estar entre 0 e 1", call. = FALSE)
    m <- alvo * (1 - confiabilidade) / (confiabilidade * (1 - alvo)); k <- ceiling(m * k_atual - 1e-9)
    tab <- data.frame(k_atual = k_atual, confiabilidade_atual = confiabilidade, alvo = alvo, razao = m, itens_necessarios = k)
    apa <- if (en) paste0("Reaching a reliability of ", .pr_num(alvo, 2, FALSE, decimal), " would require ", k, " parallel items (Spearman-Brown; currently ", k_atual,
                          " items with ", .pr_num(confiabilidade, 2, FALSE, decimal), ").") else
      paste0("Chegar à confiabilidade de ", .pr_num(alvo, 2, FALSE, decimal), " exigiria ", k, " itens paralelos (Spearman-Brown; hoje são ", k_atual, " itens com ",
             .pr_num(confiabilidade, 2, FALSE, decimal), ").")
    ef <- data.frame(indice = "itens necessários", valor = k, ic_inf = NA, ic_sup = NA, principal = TRUE)
  }
  .ee_resultado("Profecia de Spearman-Brown", n = NA, tabela = tab, efeito = ef,
                decisoes = c("Fórmula de Spearman (1910) e Brown (1910), como em Martinková e Hladká (2023, 4.2.5.1).",
                             "Vale para itens novos paralelos aos atuais: itens mais fracos rendem menos; com itens congenéricos, a teoria G (generalizabilidade()) projeta melhor."),
                apa = apa)
}

# Correções de correlação. Atenuação (Spearman; Martinková & Hladká, 2023, 4.2.3; Revelle &
# Condon, 2018): r corrigido = r / sqrt(rxx ryy); valor acima de 1 denuncia confiabilidade
# subestimada ou erros correlacionados. Restrição de amplitude por seleção direta em X
# (Thorndike, 1949, caso 2; Martinková & Hladká, 2023, 2.4.3; como psych::rangeCorrection):
# r* = (S/s) r / sqrt(1 + r² (S²/s² - 1)), com S o DP de X na população sem seleção e s o
# DP de X na amostra selecionada. Quando as duas se aplicam, a ordem segue Hunter e Schmidt
# (corrigir a restrição com a correlação observada e depois a atenuação) e fica declarada.
corrigir_correlacao <- function(r, confiabilidade_x = NULL, confiabilidade_y = NULL, dp_irrestrito = NULL, dp_restrito = NULL,
                                n = NULL, nivel = .95, decimal = ".", idioma = "pt") {
  if (abs(r) >= 1) stop("r precisa estar entre -1 e 1", call. = FALSE)
  tem_at <- !is.null(confiabilidade_x) || !is.null(confiabilidade_y)
  tem_rr <- !is.null(dp_irrestrito) || !is.null(dp_restrito)
  if (!tem_at && !tem_rr) stop("informe as confiabilidades (atenuação) e/ou os DPs irrestrito e restrito (restrição de amplitude)", call. = FALSE)
  if (tem_rr && (is.null(dp_irrestrito) || is.null(dp_restrito))) stop("restrição de amplitude pede dp_irrestrito e dp_restrito", call. = FALSE)
  rxx <- if (is.null(confiabilidade_x)) 1 else confiabilidade_x; ryy <- if (is.null(confiabilidade_y)) 1 else confiabilidade_y
  corr_rr <- function(x) { u <- dp_irrestrito / dp_restrito; u * x / sqrt(1 + x^2 * (u^2 - 1)) }
  corr_at <- function(x) x / sqrt(rxx * ryy)
  passos <- data.frame(etapa = "observada", r = r)
  atual <- r; f <- function(x) x
  if (tem_rr) { atual <- corr_rr(atual); passos <- rbind(passos, data.frame(etapa = "restrição de amplitude (Thorndike caso 2)", r = atual)); f0 <- f; f <- function(x) corr_rr(f0(x)) }
  if (tem_at) { atual <- corr_at(atual); passos <- rbind(passos, data.frame(etapa = paste0("atenuação (rxx = ", rxx, ", ryy = ", ryy, ")"), r = atual)); f1 <- f; f <- function(x) corr_at(f1(x)) }
  ic <- c(NA_real_, NA_real_)
  if (!is.null(n)) {
    z <- atanh(r); se <- 1 / sqrt(n - 3); q <- stats::qnorm(1 - (1 - nivel) / 2)
    ic <- f(tanh(c(z - q * se, z + q * se)))
  }
  avisos <- character()
  if (abs(atual) > 1) avisos <- c(avisos, "Correlação corrigida acima de 1: confiabilidade subestimada (alfa como limite inferior), erros correlacionados ou amostra pequena. Não relatar como está.")
  if (tem_rr && dp_irrestrito < dp_restrito) avisos <- c(avisos, "DP irrestrito menor que o restrito: a amostra tem mais variação que a população de referência e a correção vai reduzir r (Martinková & Hladká, 2023, notam que tirar o meio da distribuição infla r).")
  if (tem_at && (rxx < .5 || ryy < .5)) avisos <- c(avisos, "Confiabilidade abaixo de ,50: a correção para atenuação fica instável e o valor corrigido é especulativo.")
  en <- identical(idioma, "en")
  ic_txt <- if (!is.null(n)) paste0(.ee_sep(decimal), .ee_ic(ic[1], ic[2], 2, FALSE, decimal, nivel, idioma)) else ""
  apa <- if (en) paste0("The observed correlation (r = ", .pr_num(r, 2, FALSE, decimal), ") corrected for ", paste(c(if (tem_rr) "direct range restriction (Thorndike case 2)", if (tem_at) "attenuation"), collapse = " and then "),
                        " was ", .pr_num(atual, 2, FALSE, decimal), ic_txt, ".") else
    paste0("A correlação observada (r = ", .pr_num(r, 2, FALSE, decimal), "), corrigida para ", paste(c(if (tem_rr) "restrição de amplitude por seleção direta (Thorndike, caso 2)", if (tem_at) "atenuação"), collapse = " e depois para "),
           ", foi ", .pr_num(atual, 2, FALSE, decimal), ic_txt, ".")
  .ee_resultado("Correção de correlação", n = if (is.null(n)) NA else n, tabela = passos,
                efeito = data.frame(indice = "r corrigido", valor = atual, ic_inf = ic[1], ic_sup = ic[2], principal = TRUE),
                decisoes = c(if (tem_rr) "Thorndike (1949), caso 2: a seleção acontece direto em X (corte no próprio escore); seleção por outra variável (caso 3) pede a correlação dela com X e Y.",
                             if (tem_at) "Atenuação de Spearman: r / sqrt(rxx ryy); com só uma confiabilidade informada, corrige só aquele lado (a outra vale 1).",
                             if (tem_rr && tem_at) "Ordem: restrição primeiro, atenuação depois, com as confiabilidades da amostra restrita (Hunter & Schmidt).",
                             if (!is.null(n)) "IC: Fisher z da correlação observada, com os limites passados pelas mesmas correções (não inclui a incerteza das confiabilidades)."),
                avisos = avisos, apa = apa)
}

# Teoria da generalizabilidade, estudos G e D em delineamento cruzado (Martinková &
# Hladká, 2023, 4.5; Petersen, 2024, cap. 5): componentes de variância por REML
# (lme4::lmer) com uma faceta (p × i) ou duas (p × i × o); o resíduo absorve a interação
# de maior ordem. Estudo D: coeficiente de generalizabilidade Eρ² = σ²p/(σ²p + σ²rel)
# (decisões relativas) e de dependabilidade Φ = σ²p/(σ²p + σ²abs) (decisões absolutas),
# para o número observado de condições e para os tamanhos em `estudo_d`.
#   dados longos: pessoa, facetas (1 ou 2 colunas) e escore; ou largos com `itens`
#   (uma faceta, itens nas colunas).
generalizabilidade <- function(dados, pessoa, escore = NULL, facetas = NULL, itens = NULL, estudo_d = NULL,
                               decimal = ".", idioma = "pt") {
  .pr_exige("lme4")
  if (!is.null(itens)) {
    dd <- as.data.frame(dados)[, c(pessoa, itens)]
    longo <- data.frame(p = rep(dd[[pessoa]], length(itens)), item = rep(itens, each = nrow(dd)),
                        y = unlist(dd[itens], use.names = FALSE))
    facetas <- "item"; names(longo)[1] <- pessoa; escore <- "y"; dados <- longo
  }
  if (is.null(escore) || !length(facetas) || length(facetas) > 2) stop("informe escore e 1 ou 2 facetas (ou itens = colunas, para uma faceta)", call. = FALSE)
  base <- .ee_dados(dados, c(pessoa, facetas, escore), numericas = escore)
  d <- base$d
  for (v in c(pessoa, facetas)) d[[v]] <- factor(d[[v]])
  cel <- table(do.call(paste, c(d[c(pessoa, facetas)], sep = "|")))
  if (any(cel > 1)) stop("há mais de uma observação por célula (pessoa × condição): o delineamento daqui é cruzado com uma observação por célula", call. = FALSE)
  P <- pessoa; F1 <- facetas[1]; F2 <- if (length(facetas) == 2) facetas[2] else NULL
  termos <- c(P, F1, if (!is.null(F2)) c(F2, paste0(P, ":", F1), paste0(P, ":", F2), paste0(F1, ":", F2)))
  fo <- stats::as.formula(paste(escore, "~ 1 +", paste0("(1|", termos, ")", collapse = " + ")))
  fit <- suppressMessages(suppressWarnings(lme4::lmer(fo, data = d, REML = TRUE)))
  vc <- as.data.frame(lme4::VarCorr(fit))
  comp <- stats::setNames(vc$vcov, ifelse(vc$grp == "Residual", "residuo", vc$grp))
  rot_res <- if (is.null(F2)) paste0(P, ":", F1, ", e") else paste0(P, ":", F1, ":", F2, ", e")
  g_tab <- data.frame(componente = c(termos, rot_res), variancia = unname(comp[c(termos, "residuo")]))
  g_tab$percentual <- 100 * g_tab$variancia / sum(g_tab$variancia)
  n_obs <- c(nlevels(d[[F1]]), if (!is.null(F2)) nlevels(d[[F2]]))
  names(n_obs) <- facetas
  vp <- comp[[P]]
  d_coef <- function(n1, n2 = 1) {
    if (is.null(F2)) {
      rel <- comp[["residuo"]] / n1; abs_ <- rel + comp[[F1]] / n1
    } else {
      rel <- comp[[paste0(P, ":", F1)]] / n1 + comp[[paste0(P, ":", F2)]] / n2 + comp[["residuo"]] / (n1 * n2)
      abs_ <- rel + comp[[F1]] / n1 + comp[[F2]] / n2 + comp[[paste0(F1, ":", F2)]] / (n1 * n2)
    }
    c(var_rel = rel, var_abs = abs_, ep_rel = sqrt(rel), ep_abs = sqrt(abs_), g = vp / (vp + rel), phi = vp / (vp + abs_))
  }
  grade <- if (is.null(estudo_d)) as.data.frame(as.list(n_obs)) else {
    if (is.null(names(estudo_d))) stop("estudo_d precisa ser lista nomeada pelas facetas, ex.: list(item = c(5, 10, 20))", call. = FALSE)
    e <- lapply(facetas, function(f) if (!is.null(estudo_d[[f]])) estudo_d[[f]] else n_obs[[f]]); names(e) <- facetas
    unique(rbind(as.data.frame(as.list(n_obs)), expand.grid(e)))
  }
  d_tab <- cbind(grade, t(apply(grade, 1, function(l) d_coef(l[[1]], if (!is.null(F2)) l[[2]] else 1))))
  d_tab$observado <- apply(grade, 1, function(l) all(l == n_obs))
  obs <- d_tab[d_tab$observado, ][1, ]
  avisos <- character()
  if (lme4::isSingular(fit)) avisos <- c(avisos, paste0("Componente de variância estimado em zero (ajuste singular: ", paste(g_tab$componente[g_tab$variancia < 1e-10], collapse = ", "),
                                                       "): comum com poucos níveis da faceta; o estudo D trata como zero."))
  if (nlevels(d[[P]]) < 20) avisos <- c(avisos, paste0("Só ", nlevels(d[[P]]), " pessoas: componentes de variância imprecisos; o estudo D é ilustrativo."))
  if (base$n_excluidos > 0) avisos <- c(avisos, paste0(base$n_excluidos, " observação(ões) ausente(s): o REML usa o resto, mas o delineamento fica desbalanceado."))
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  n_txt <- paste(paste(n_obs, facetas), collapse = if (en) " and " else " e ")
  apa <- if (en) paste0("A G-study (", paste(c(P, facetas), collapse = " × "), ", crossed, REML) attributed ", .pr_num(g_tab$percentual[1], 1, TRUE, decimal),
                        "% of the variance to persons. For ", n_txt, ", Eρ² = ", .pr_num(obs$g, 2, FALSE, decimal), " and Φ = ", .pr_num(obs$phi, 2, FALSE, decimal), ".") else
    paste0("O estudo G (", paste(c(P, facetas), collapse = " × "), ", cruzado, REML) atribuiu ", .pr_num(g_tab$percentual[1], 1, TRUE, decimal),
           "% da variância às pessoas. Com ", n_txt, ", Eρ² = ", .pr_num(obs$g, 2, FALSE, decimal), sep, "Φ = ", .pr_num(obs$phi, 2, FALSE, decimal), ".")
  .ee_resultado(paste0("Teoria da generalizabilidade (", if (is.null(F2)) "uma faceta" else "duas facetas", ")"), n = nlevels(d[[P]]), n_excluidos = base$n_excluidos,
                tabela = g_tab,
                efeito = data.frame(indice = c("Eρ² (generalizabilidade)", "Φ (dependabilidade)"), valor = c(obs$g, obs$phi), ic_inf = NA, ic_sup = NA, principal = c(TRUE, TRUE)),
                decisoes = c("Componentes por REML (lme4), como Martinková e Hladká (2023, 4.4.3 e 4.5) e Petersen (2024, cap. 5); com dados balanceados coincidem com os da ANOVA.",
                             "Eρ² usa o erro relativo (interações com a pessoa): ordenar pessoas. Φ soma o erro absoluto (efeitos principais das facetas): decidir por ponto de corte.",
                             "Estudo D na tabela extra; a projeção supõe as mesmas condições amostradas do mesmo universo, a generalização da profecia de Spearman-Brown."),
                avisos = avisos, apa = apa,
                extra = list(tabelas_extra = list("Estudo D" = d_tab), componentes = comp, estudo_d = d_tab, modelo = fit))
}

# Análise de distratores (v1.6; Martinková & Hladká, 2023, 5.5): itens de múltipla escolha
# com as respostas não corrigidas (letras) e o gabarito. Escore total = número de acertos;
# grupos pelos quantis do total, como em analise_classica_itens() e no
# ShinyItemAnalysis::DistractorAnalysis. Para cada opção: proporção que a escolhe em cada
# grupo e discriminação (proporção no grupo superior menos no inferior). O gabarito deve
# subir com o escore; distrator que sobe com o escore ou que ninguém escolhe é sinal de
# item a revisar.
analise_distratores <- function(dados, itens, gabarito, grupos = 3, minimo_escolha = .02, decimal = ".", idioma = "pt") {
  if (length(gabarito) != length(itens)) stop("gabarito precisa ter uma resposta correta por item", call. = FALSE)
  dd <- .ee_dados(dados, itens)
  d <- as.data.frame(lapply(dd$d, as.character), stringsAsFactors = FALSE); n <- nrow(d)
  gab <- as.character(unlist(gabarito))
  acertos <- as.data.frame(Map(function(x, g) as.integer(x == g), d, gab))
  tot <- rowSums(acertos)
  q <- stats::quantile(tot, probs = seq(0, 1, length.out = grupos + 1))
  gr <- cut(tot, unique(q), include.lowest = TRUE, labels = FALSE)
  ng <- max(gr)
  linhas <- list()
  for (j in seq_along(itens)) {
    niv <- if (is.factor(dados[[itens[j]]])) levels(dados[[itens[j]]]) else sort(unique(d[[j]]))
    niv <- union(niv, gab[j])
    for (o in niv) {
      props <- vapply(seq_len(ng), function(g) mean(d[[j]][gr == g] == o), 0)
      linhas[[length(linhas) + 1]] <- data.frame(item = itens[j], opcao = o, gabarito = o == gab[j], proporcao_total = mean(d[[j]] == o),
                                                 t(stats::setNames(props, paste0("grupo", seq_len(ng)))), discriminacao = props[ng] - props[1],
                                                 row.names = NULL, check.names = FALSE)
    }
  }
  tab <- do.call(rbind, linhas)
  avisos <- character()
  gab_ruim <- tab$item[tab$gabarito & tab$discriminacao <= 0]
  dis_pos <- tab[!tab$gabarito & tab$discriminacao > 0, c("item", "opcao")]
  vazios <- tab[!tab$gabarito & tab$proporcao_total < minimo_escolha, c("item", "opcao")]
  if (length(gab_ruim)) avisos <- c(avisos, paste0("Gabarito que não cresce com o escore (discriminação <= 0): ", paste(unique(gab_ruim), collapse = ", "), ". Conferir a chave e a redação."))
  if (nrow(dis_pos)) avisos <- c(avisos, paste0("Distrator escolhido mais pelo grupo superior que pelo inferior: ", paste(dis_pos$item, dis_pos$opcao, sep = " ", collapse = "; "),
                                                ". Pode haver segunda resposta defensável ou ambiguidade."))
  if (nrow(vazios)) avisos <- c(avisos, paste0("Distrator escolhido por menos de ", .pr_num(100 * minimo_escolha, 0, TRUE, decimal), "% (não funciona como alternativa): ",
                                               paste(vazios$item, vazios$opcao, sep = " ", collapse = "; "), "."))
  en <- identical(idioma, "en")
  apa <- if (en) paste0("Distractor analysis (", ng, " total-score groups; Martinková & Hladká, 2023): ", length(unique(gab_ruim)), " keys did not increase with total score, ",
                        nrow(dis_pos), " distractors were chosen more often by the upper group, and ", nrow(vazios), " were rarely chosen.") else
    paste0("Na análise de distratores (", ng, " grupos pelo escore total; Martinková & Hladká, 2023), ", length(unique(gab_ruim)), " gabaritos não cresceram com o escore, ",
           nrow(dis_pos), " distratores foram mais escolhidos pelo grupo superior e ", nrow(vazios), " quase não foram escolhidos.")
  .ee_resultado("Análise de distratores", n = n, n_excluidos = dd$n_excluidos, tabela = tab,
                decisoes = c(paste0("Escore = número de acertos; ", ng, " grupos pelos quantis do escore (empates no limite mudam o tamanho dos grupos, como no ShinyItemAnalysis)."),
                             "Discriminação da opção = proporção no grupo superior menos no inferior: positiva no gabarito, negativa nos distratores."),
                avisos = avisos, apa = apa, extra = list(grupo = gr, escore = tot))
}

# ---------------------------------------------------------------- normas (v1.7)
# normas(): tabela de normas a partir de uma amostra normativa (Petersen, 2024, cap. 1,
# escores referidos à norma): para cada escore bruto possível, percentil, z (com a média e o
# DP da amostra normativa, por estrato se houver `grupos`), T = 50 + 10z, escore escalar
# 10 + 3z, padrão 100 + 15z e estanino pelas faixas de Petersen (4, 11, 23, 40, 60, 77, 89,
# 96%). Percentil: "medio" (padrão, posto médio: abaixo + metade dos iguais, sobre n) ou
# "petersen" (dplyr::percent_rank do livro: (posto mínimo - 1)/(n - 1)). normalizar = TRUE
# acrescenta o z normalizado (qnorm do percentil médio), para escores assimétricos. A
# tabela é agregada; o print não mostra frequência por escore.
normas <- function(dados, escore, grupos = NULL, percentil = c("medio", "petersen"), normalizar = FALSE,
                   valores = NULL, minimo_n = 30, decimal = ".", idioma = "pt") {
  percentil <- match.arg(percentil)
  dd <- .ee_dados(dados, c(escore, grupos), numericas = escore)
  d <- dd$d
  est <- if (is.null(grupos)) rep("(todos)", nrow(d)) else interaction(d[grupos], drop = TRUE, sep = " | ")
  est <- factor(est)
  linhas <- list(); resumo <- list(); avisos <- character()
  for (g in levels(est)) {
    x <- d[[escore]][est == g]; n <- length(x); m <- mean(x); s <- stats::sd(x)
    vv <- if (is.null(valores)) sort(unique(x)) else sort(unique(valores))
    pr <- vapply(vv, function(v) if (percentil == "medio") 100 * (sum(x < v) + .5 * sum(x == v)) / n else
      100 * sum(x < v) / (n - 1), 0)
    z <- (vv - m) / s
    est9 <- cut(pr, c(-Inf, 4, 11, 23, 40, 60, 77, 89, 96, Inf), labels = 1:9, right = TRUE)
    tab <- data.frame(grupo = g, bruto = vv, percentil = pr, z = z, T = 50 + 10 * z, escalar = 10 + 3 * z, padrao = 100 + 15 * z,
                      estanino = as.integer(as.character(est9)))
    if (normalizar) {
      pm <- vapply(vv, function(v) (sum(x < v) + .5 * sum(x == v)) / n, 0)
      tab$z_normalizado <- stats::qnorm(pmin(pmax(pm, .5 / n), 1 - .5 / n)); tab$T_normalizado <- 50 + 10 * tab$z_normalizado
    }
    linhas[[g]] <- tab
    resumo[[g]] <- data.frame(grupo = g, n = n, media = m, dp = s, assimetria = .ee_assimetria(x), minimo = min(x), maximo = max(x))
    if (n < minimo_n) avisos <- c(avisos, paste0("Estrato ", g, " com n = ", n, " (< ", minimo_n, "): percentis instáveis nas pontas; juntar estratos ou ampliar a amostra normativa."))
  }
  tab <- do.call(rbind, linhas); rownames(tab) <- NULL
  res <- do.call(rbind, resumo); rownames(res) <- NULL
  if (any(abs(res$assimetria) > 1, na.rm = TRUE) && !normalizar)
    avisos <- c(avisos, "Assimetria acima de 1 em algum estrato: z, T e padrão lineares herdam a assimetria e o percentil deixa de bater com a curva normal; considerar normalizar = TRUE.")
  en <- identical(idioma, "en")
  apa <- if (en) paste0("Norms for ", escore, " were derived from the normative sample (", paste0(res$grupo, ": n = ", res$n, collapse = "; "),
                        "), with ", if (percentil == "medio") "mid-point percentile ranks" else "percentile ranks as in Petersen (2024)", ", z, T, and stanine scores.") else
    paste0("As normas de ", escore, " foram derivadas da amostra normativa (", paste0(res$grupo, ": n = ", res$n, collapse = "; "), "), com percentil ",
           if (percentil == "medio") "pelo posto médio" else "como em Petersen (2024)", ", z, T e estanino.")
  .ee_resultado("Tabela de normas", n = nrow(d), n_excluidos = dd$n_excluidos, tabela = res,
                decisoes = c(paste0("Percentil ", if (percentil == "medio") "pelo posto médio: (abaixo + metade dos iguais)/n, que não dá 0 nem 100 a quem está no extremo" else "do livro de Petersen (percent_rank: (posto mínimo - 1)/(n - 1))", "."),
                             "z com a média e o DP do estrato; T = 50 + 10z; escalar = 10 + 3z; padrão = 100 + 15z; estanino pelas faixas cumulativas 4, 11, 23, 40, 60, 77, 89 e 96% (Petersen, 2024, tabela 1.1).",
                             "Estratos por grupos (idade, escolaridade): cada um com a própria média, DP e percentis. A norma vale para a população que a amostra normativa representa; descrever a amostragem no relato."),
                avisos = avisos, apa = apa, extra = list(tabelas_extra = list(`Normas por escore bruto` = tab), normas = tab, resumo = res))
}
