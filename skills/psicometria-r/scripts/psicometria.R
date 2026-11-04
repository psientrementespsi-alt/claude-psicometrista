# psicometria.R, funções de psicometria da skill psicometria-r (v1.5.0)
#
# Uso, no R do computador onde estão os dados:
#   source("~/.local/share/psicometria-r/carregar.R")  (traz este arquivo e os módulos est_*.R)
#
# Regra de privacidade: nenhuma função da skill imprime linha de participante. Tudo o
# que a impressão mostra é agregado (contagens por categoria, médias, cargas,
# índices). As poucas funções que criam variável por pessoa (índice de aquiescência e
# dados recentralizados, classe latente) guardam o vetor no objeto, para a análise
# seguinte no computador do usuário, e não o imprimem. Dado bruto de paciente não sai
# do computador do usuário.
#
# Testado com R 4.5.3 (evals da v1.3, 26/10/2026) e R 4.6.1 (v1.1), lavaan 0.7.3,
# semTools 0.5.10, mirt 1.48, bootnet 1.9.1, EGAnet 2.4.1 e NetworkComparisonTest 2.2.4
# (scripts/evals_simulados.R).

.pr_exige <- function(...) {
  p <- c(...)
  falta <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
  if (length(falta)) {
    stop("Pacotes ausentes: ", paste(falta, collapse = ", "),
         ". Instale no R do computador onde estão os dados.", call. = FALSE)
  }
}

# Itens com até 7 valores distintos contam como ordinais (a mesma regra do
# cor_auto do qgraph, usado pelo bootnet).
.pr_e_ordinal <- function(dados, itens, max_cat = 7) {
  all(vapply(itens, function(i) length(unique(stats::na.omit(dados[[i]]))) <= max_cat, logical(1)))
}

# Fatores com indicadores observados (primeira ordem) e se há fator de ordem superior.
.pr_estrutura <- function(fit, dados) {
  pt <- lavaan::parTable(fit)
  carga <- unique(pt[pt$op == "=~" & pt$rhs %in% names(dados), c("lhs", "rhs")])
  latentes <- unique(pt$lhs[pt$op == "=~"])
  list(carga = carga, fatores = unique(carga$lhs),
       superior = any(pt$op == "=~" & pt$rhs %in% latentes))
}

# 1. Auditoria dos itens -------------------------------------------------------
# Por item: n válido, % ausente, número de valores distintos, frequência da
# categoria mais rara (itens com até 10 valores), assimetria e curtose em excesso.
# Com `grupo`, também a menor frequência de categoria dentro de cada grupo,
# contando zero para categoria que existe no total e falta num grupo (a AFC
# ordinal multigrupo falha com categoria vazia num grupo).
auditar_itens <- function(dados, itens, grupo = NULL, min_celula = 10) {
  linhas <- lapply(itens, function(i) {
    x <- dados[[i]]
    ok <- !is.na(x)
    v <- as.numeric(x[ok])
    k <- length(unique(v))
    m <- mean(v)
    s <- stats::sd(v)
    lin <- data.frame(
      item = i, n = length(v), pct_ausente = round(100 * mean(!ok), 1),
      categorias = k,
      menor_categoria = if (k <= 10) as.integer(min(table(v))) else NA_integer_,
      assimetria = round(mean((v - m)^3) / s^3, 2),
      curtose_excesso = round(mean((v - m)^4) / s^4 - 3, 2),
      stringsAsFactors = FALSE
    )
    if (!is.null(grupo)) {
      if (k <= 10) {
        tab <- table(factor(v, levels = sort(unique(v))), droplevels(as.factor(dados[[grupo]][ok])))
        lin$menor_categoria_no_grupo <- as.integer(min(tab))
        lin$categoria_vazia_em_grupo <- any(tab == 0)
      } else {
        lin$menor_categoria_no_grupo <- NA_integer_
        lin$categoria_vazia_em_grupo <- NA
      }
    }
    lin
  })
  out <- do.call(rbind, linhas)
  alerta <- ifelse(!is.na(out$menor_categoria) & out$menor_categoria < min_celula,
                   paste0("categoria com menos de ", min_celula,
                          " casos: avaliar colapsar categorias adjacentes"), "")
  if (!is.null(grupo)) {
    vazia <- out$categoria_vazia_em_grupo %in% TRUE
    alerta[vazia] <- trimws(paste(alerta[vazia],
      "| categoria vazia em algum grupo: colapsar do mesmo jeito em todos os grupos antes da invariância"))
    alerta <- sub("^\\| ", "", alerta)
  }
  out$alerta <- alerta
  rownames(out) <- NULL
  out
}

# 2. Sugestão de estimador -----------------------------------------------------
# Rhemtulla et al. (2012): cat-LS (família do WLSMV) para menos de 5 categorias;
# ML robusto com 5 ou mais categorias quando a amostra é pequena e os limiares
# são aproximadamente simétricos; com 6 ou 7 categorias, qualquer um dos dois em
# muitas condições. ML fica enviesado com limiares assimétricos. A assimetria do
# item entra como aproximação da assimetria dos limiares, e "amostra pequena"
# como n até 200 (Li, 2016, achou os problemas do WLSMV com N = 200): as duas
# operacionalizações são da skill, não do artigo.
sugerir_estimador <- function(auditoria, assimetria_max = 1, n_pequeno = 200) {
  k <- min(auditoria$categorias)
  a <- max(abs(auditoria$assimetria), na.rm = TRUE)
  n <- min(auditoria$n)
  wlsmv <- function(motivo, sens = "MLR") list(estimador = "WLSMV", ordinal = TRUE, motivo = motivo, sensibilidade = sens)
  if (k < 5) {
    est <- wlsmv(sprintf("menor número de categorias = %d (< 5): cat-LS recomendado (Rhemtulla et al., 2012)", k), NA_character_)
  } else if (k <= 7) {
    simetrico <- a <= assimetria_max
    if (!simetrico) {
      est <- wlsmv(sprintf("%d categorias com assimetria de até %.2f: limiares provavelmente assimétricos, cenário em que ML fica enviesado (Rhemtulla et al., 2012); MLR só como sensibilidade", k, a))
    } else if (n <= n_pequeno) {
      est <- list(estimador = "MLR", ordinal = FALSE, sensibilidade = "WLSMV",
                  motivo = sprintf("%d categorias, limiares aproximadamente simétricos (assimetria até %.2f) e n = %d: com 5 ou mais categorias, amostra pequena e limiares simétricos, Rhemtulla et al. (2012) recomendam ML robusto", k, a, n))
    } else if (k == 5) {
      est <- wlsmv(sprintf("5 categorias, limiares aproximadamente simétricos e n = %d: a recomendação de ML robusto de Rhemtulla et al. (2012) é para amostra pequena; cat-LS estima as cargas com mais precisão (Rhemtulla et al., 2012; Li, 2016); MLR como sensibilidade", n))
    } else {
      est <- wlsmv(sprintf("%d categorias, limiares aproximadamente simétricos: com 6 ou 7 categorias os dois métodos são aceitáveis em muitas condições (Rhemtulla et al., 2012); WLSMV como principal por escolha da skill, MLR como sensibilidade", k))
    }
  } else {
    est <- list(estimador = "MLR", ordinal = FALSE, sensibilidade = NA_character_,
                motivo = sprintf("%d ou mais valores distintos (escala quase contínua, como 0-100): MLR com FIML para ausentes", k))
  }
  if (est$estimador == "WLSMV" && n <= n_pequeno) {
    est$motivo <- paste0(est$motivo, sprintf(
      "; atenção: n = %d, e cat-LS é mais sensível a amostra pequena (Rhemtulla et al., 2012); com N = 200, Li (2016) achou superestimação moderada das correlações entre fatores", n))
  }
  est
}

# 3. AFC -----------------------------------------------------------------------
# Ausentes: com WLSMV o lavaan só oferece listwise ou pairwise, e a deleção
# pairwise pressupõe ausência completamente ao acaso (MCAR); se não for o caso,
# a alternativa é imputação múltipla. Com ML ou MLR, FIML.
.pr_ausentes <- function(ordinal, estimador) {
  if (ordinal) "pairwise" else if (estimador %in% c("ML", "MLR")) "fiml" else "listwise"
}

ajustar_afc <- function(modelo, dados, itens, estimador = NULL, ordinal = NULL, ...) {
  .pr_exige("lavaan")
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(dados, itens)
  if (is.null(estimador)) estimador <- if (ordinal) "WLSMV" else "MLR"
  
  args <- list(
    model = modelo,
    data = quote(dados),
    estimator = estimador,
    missing = .pr_ausentes(ordinal, estimador)
  )
  if (ordinal) {
    args$ordered <- itens
    args$parameterization <- "theta"
  }
  extra <- list(...)
  if (length(extra) > 0) args <- c(args, extra)
  
  cl <- as.call(c(list(quote(lavaan::cfa)), args))
  eval(cl)
}

# Índices de ajuste. Para WLSMV e MLR o lavaan dá duas versões dos índices
# incrementais e do RMSEA: "escalonada" (a tradicional) e "robusta" (a versão
# mais nova, que a documentação do lavaan 0.7 recomenda usar sempre que
# disponível). A tabela traz as duas, rotuladas; o relato diz qual foi usada.
# Com dado ordinal, a robusta sai NA quando a matriz policórica não é positiva
# definida (comum em amostra pequena, segundo a mesma documentação).
tabela_ajuste <- function(fit, rotulo = "modelo") {
  .pr_exige("lavaan")
  fm <- lavaan::fitMeasures(fit)
  pega <- function(nome) if (nome %in% names(fm)) unname(fm[[nome]]) else NA_real_
  versoes <- c(if (!is.na(pega("cfi.scaled"))) "escalonada",
               if (!is.na(pega("cfi.robust"))) "robusta")
  if (!length(versoes)) versoes <- "padrão"
  if (!is.na(pega("cfi.scaled")) && is.na(pega("cfi.robust"))) {
    message("Índices robustos indisponíveis (NA) para ", rotulo,
            ": relatar a versão escalonada e dizer por quê.")
  }
  linhas <- lapply(versoes, function(v) {
    suf <- switch(v, escalonada = ".scaled", robusta = ".robust", "")
    # Com ML e FIML o lavaan dá índices robustos sem qui-quadrado escalonado: o teste
    # do modelo é o padrão.
    qui_suf <- if (v == "padrão" || is.na(pega("chisq.scaled"))) "" else ".scaled"
    data.frame(
      modelo = rotulo, versao = v,
      qui2 = pega(paste0("chisq", qui_suf)), gl = pega(paste0("df", qui_suf)),
      p = pega(paste0("pvalue", qui_suf)),
      cfi = pega(paste0("cfi", suf)), tli = pega(paste0("tli", suf)),
      rmsea = pega(paste0("rmsea", suf)),
      rmsea_li = pega(paste0("rmsea.ci.lower", suf)), rmsea_ls = pega(paste0("rmsea.ci.upper", suf)),
      srmr = pega("srmr"), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, linhas)
}

# 4. Confiabilidade -----------------------------------------------------------
# Ômega por fator via semTools::compRelSEM. Com itens ordinais (ord.scale = TRUE,
# o padrão), é a confiabilidade do escore somado observado (Green & Yang, 2009).
# Alfa de Cronbach calculado nos escores observados, só para relato conjunto
# (a fabrica-de-artigo exige alfa e ômega; McNeish, 2018).
# Só fatores com indicadores observados (primeira ordem). Em modelo de segunda
# ordem o semTools trata o fator geral como único escore verdadeiro, e o ômega
# de cada subescala sai como a parte da variância dela devida ao fator geral,
# não como a confiabilidade da subescala; para isso, ajustar o modelo de
# fatores correlacionados (os evals mostram a diferença). Modelo de um grupo.
# ic = TRUE (v1.6): intervalo por bootstrap percentílico dos casos, reajustando o mesmo
# modelo (mesmo estimador e tratamento de ausentes) em cada reamostra; é a rota que
# Kelley e Pornprasertmanit (2016) recomendam para o ômega, e vale também para o alfa
# (Martinková & Hladká, 2023, 4.4.4). Reamostra que não converge é descartada e contada.
# Confiabilidade de cada subescala (ômega de McDonald e alfa de Cronbach).
# Suporta metodo = "delta" (analitico instantaneo via Raykov, 2001 e Feldt, 1965)
# e metodo = "bootstrap" (com correcao de ambiente e descarte de falhas).
confiabilidade <- function(fit, dados, ic = FALSE, metodo = c("delta", "bootstrap"), reamostras = 1000, nivel = .95, semente = 2026) {
  metodo <- match.arg(metodo)
  .pr_exige("lavaan", "semTools")
  if (lavaan::lavInspect(fit, "ngroups") > 1) {
    stop("confiabilidade() é para modelo de um grupo: ajustar e rodar por grupo.", call. = FALSE)
  }
  est <- .pr_estrutura(fit, dados)
  carga <- est$carga
  fatores <- est$fatores
  if (est$superior) {
    message("Modelo com fator de ordem superior: o ômega por subescala é relativo ao fator geral. ",
            "Para a confiabilidade de cada subescala, usar o modelo de fatores correlacionados.")
  }
  std <- lavaan::standardizedSolution(fit)
  neg <- unique(std$rhs[std$op == "=~" & std$rhs %in% carga$rhs & std$est.std < 0])
  if (length(neg)) {
    message("Carga padronizada negativa em: ", paste(neg, collapse = ", "),
            ". Item invertido sem recodificar distorce alfa e ômega de pesos unitários: recodificar antes.")
  }
  omega_raw <- semTools::compRelSEM(fit)
  omega_vals <- vapply(fatores, function(f) {
    if (f %in% names(omega_raw)) as.numeric(omega_raw[[f]]) else NA_real_
  }, numeric(1))
  
  alfa <- vapply(fatores, function(f) {
    x <- stats::na.omit(dados[, carga$rhs[carga$lhs == f], drop = FALSE])
    x <- as.data.frame(lapply(x, as.numeric))
    k <- ncol(x)
    if (k < 2) return(NA_real_)
    (k / (k - 1)) * (1 - sum(vapply(x, stats::var, 0)) / stats::var(rowSums(x)))
  }, numeric(1))
  
  itens_por_fator <- as.vector(table(carga$lhs)[fatores])
  out <- data.frame(fator = fatores, omega = round(omega_vals, 3),
                    alfa = round(alfa, 3), itens = itens_por_fator,
                    row.names = NULL)
  if (!isTRUE(ic)) return(out)
  
  n <- nrow(dados)
  
  if (metodo == "delta") {
    # v1.8.1: o IC do omega sai de um refit com std.lv = TRUE em que o omega e um parametro definido (delta
    # method do lavaan). Itens ordinais: theta = 1 - lambda^2 (parametrizacao delta); continuos: theta rotulado.
    # O IC fica centrado em omega_delta. A v1.8.0 aplicava o gradiente de L^2/(L^2 + E) sem a variancia do fator e
    # sem os parametros fixos, e centrava o resultado no omega de compRelSEM, que e outra quantidade.
    # IC do alfa: Feldt (1965), com o n listwise de cada fator.
    k_f <- length(fatores)
    alfa_ic_inf <- alfa_ic_sup <- om_d <- om_ep <- om_inf <- om_sup <- rep(NA_real_, k_f)
    for (idx_f in seq_len(k_f)) {
      itf <- carga$rhs[carga$lhs == fatores[idx_f]]
      n_a <- nrow(stats::na.omit(dados[, itf, drop = FALSE])); k <- length(itf); a_val <- alfa[idx_f]
      if (!is.na(a_val) && k >= 2 && n_a > 2) {
        df1 <- n_a - 1; df2 <- (n_a - 1) * (k - 1)
        alfa_ic_inf[idx_f] <- round(max(0, 1 - (1 - a_val) * stats::qf(1 - (1 - nivel) / 2, df1, df2)), 3)
        alfa_ic_sup[idx_f] <- round(min(1, 1 - (1 - a_val) * stats::qf((1 - nivel) / 2, df1, df2)), 3)
      }
    }
    pt <- lavaan::parTable(fit); ov <- lavaan::lavNames(fit, "ov"); ord <- lavaan::lavNames(fit, "ov.ord")
    cruzadas <- any(duplicated(carga$rhs))
    cov_res <- any(pt$op == "~~" & pt$lhs != pt$rhs & pt$lhs %in% ov & pt$rhs %in% ov & pt$free > 0)
    if (isTRUE(est$superior) || cruzadas || cov_res) {
      message("omega_delta exige fatores de primeira ordem sem carga cruzada nem covariancia residual: IC do omega nao calculado.")
    } else {
      linhas <- character(); defs <- character()
      for (idx_f in seq_len(k_f)) {
        f <- fatores[idx_f]; itf <- carga$rhs[carga$lhs == f]; lb <- sprintf("prl%d_%d", idx_f, seq_along(itf))
        linhas <- c(linhas, paste0(f, " =~ ", paste0(lb, "*", itf, collapse = " + ")))
        E <- if (all(itf %in% ord)) paste0("(1 - ", lb, "^2)", collapse = " + ") else if (!any(itf %in% ord)) {
          th <- sprintf("prt%d_%d", idx_f, seq_along(itf)); linhas <- c(linhas, paste0(itf, " ~~ ", th, "*", itf)); paste(th, collapse = " + ")
        } else NA_character_
        if (!is.na(E)) { L <- paste(lb, collapse = " + "); defs <- c(defs, sprintf("prom%d := (%s)^2 / ((%s)^2 + %s)", idx_f, L, L, E)) }
      }
      opt <- lavaan::lavInspect(fit, "options")
      fit_d <- tryCatch(suppressWarnings(lavaan::cfa(paste(c(linhas, defs), collapse = "\n"), data = dados,
        ordered = if (length(ord)) ord else NULL, estimator = opt$estimator, se = opt$se, test = opt$test,
        missing = opt$missing, std.lv = TRUE, parameterization = "delta")), error = function(e) NULL)
      if (is.null(fit_d) || !isTRUE(lavaan::lavInspect(fit_d, "converged"))) {
        message("Refit do omega delta nao convergiu: IC do omega nao calculado.")
      } else {
        d_chi <- abs(unname(lavaan::fitMeasures(fit_d, "chisq") - lavaan::fitMeasures(fit, "chisq")))
        if (is.finite(d_chi) && d_chi > 1e-3) message("Refit do omega delta difere do modelo original (diferenca de qui2 = ", round(d_chi, 4), "): conferir a especificacao.")
        pe <- lavaan::parameterEstimates(fit_d, level = nivel)
        for (idx_f in seq_len(k_f)) {
          r <- pe[pe$op == ":=" & pe$lhs == sprintf("prom%d", idx_f), ]
          if (nrow(r)) { om_d[idx_f] <- r$est; om_ep[idx_f] <- r$se; om_inf[idx_f] <- r$ci.lower; om_sup[idx_f] <- r$ci.upper }
        }
        if (length(ord)) message("Itens ordinais: omega (compRelSEM) e do escore somado observado; omega_delta e da variavel-resposta latente, e o IC vale para omega_delta.")
      }
    }
    out$omega_delta <- round(om_d, 3); out$omega_ep <- round(om_ep, 4)
    out$omega_ic_inf <- round(pmax(0, om_inf), 3); out$omega_ic_sup <- round(pmin(1, om_sup), 3)
    out$alfa_ic_inf <- alfa_ic_inf; out$alfa_ic_sup <- alfa_ic_sup
    attr(out, "metodo_ic") <- "omega: delta method do lavaan sobre omega_delta (Raykov, 2001); alfa: Feldt (1965)"
    return(out)
  }
  
  # Metodo Bootstrap: sanitiza fit@call contra simbolos nao avaliados
  cl <- fit@call
  opt <- lavaan::lavInspect(fit, "options")
  cl$estimator <- opt$estimator
  cl$missing <- opt$missing
  ord_names <- lavaan::lavNames(fit, "ov.ord")
  if (length(ord_names) > 0) cl$ordered <- ord_names
  if (!is.null(opt$parameterization)) cl$parameterization <- opt$parameterization
  fit_boot <- fit
  fit_boot@call <- cl
  
  alfa_de <- function(dd) vapply(fatores, function(f) {
    x <- stats::na.omit(dd[, carga$rhs[carga$lhs == f], drop = FALSE])
    x <- as.data.frame(lapply(x, as.numeric))
    k <- ncol(x)
    if (k < 2) return(NA_real_)
    (k / (k - 1)) * (1 - sum(vapply(x, stats::var, 0)) / stats::var(rowSums(x)))
  }, numeric(1))
  
  if (exists(".Random.seed", envir = globalenv())) {
    velha <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", velha, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(semente)
  message("Bootstrap da confiabilidade: ", reamostras, " reajustes do modelo (pode levar alguns minutos com WLSMV).")
  
  B <- t(vapply(seq_len(reamostras), function(b) {
    db <- dados[sample.int(n, n, replace = TRUE), , drop = FALSE]
    fb <- tryCatch(suppressWarnings(lavaan::update(fit_boot, data = db)), error = function(e) NULL)
    if (is.null(fb) || !isTRUE(lavaan::lavInspect(fb, "converged"))) return(rep(NA_real_, 2 * length(fatores)))
    om_boot <- tryCatch({
      rel <- semTools::compRelSEM(fb)
      vapply(fatores, function(f) if (f %in% names(rel)) as.numeric(rel[[f]]) else NA_real_, numeric(1))
    }, error = function(e) rep(NA_real_, length(fatores)))
    c(om_boot, alfa_de(db))
  }, numeric(2 * length(fatores))))
  
  falhas <- sum(!stats::complete.cases(B))
  pr <- c((1 - nivel) / 2, 1 - (1 - nivel) / 2)
  q <- apply(B, 2, stats::quantile, probs = pr, na.rm = TRUE, names = FALSE)
  j <- seq_along(fatores)
  out$omega_ic_inf <- round(q[1, j], 3)
  out$omega_ic_sup <- round(q[2, j], 3)
  out$alfa_ic_inf <- round(q[1, length(fatores) + j], 3)
  out$alfa_ic_sup <- round(q[2, length(fatores) + j], 3)
  attr(out, "bootstrap") <- list(reamostras = reamostras, falhas = falhas, nivel = nivel, semente = semente)
  attr(out, "metodo_ic") <- "bootstrap"
  if (falhas > 0) message(falhas, " reamostra(s) sem convergência descartada(s); se passar de 5%, o IC fica otimista.")
  out
}

# Bootstrap da confiabilidade em blocos (v1.9.0) ---------------------------
# O bootstrap de confiabilidade() reajusta o modelo B vezes (cerca de 0,6 s por
# reamostra com WLSMV e 21 itens): 1.000 reamostras passam do limite de 60 s
# do r-btw. Esta funcao gera o plano de reamostras uma vez (mesma sequencia de
# sample.int que confiabilidade() usa com a mesma semente), roda blocos ate
# 'max_segundos', grava cada bloco em 'pasta' e, chamada de novo, continua de
# onde parou. Com todos os blocos feitos, devolve a mesma tabela que
# confiabilidade(ic = TRUE, metodo = "bootstrap"): os evals conferem a igualdade.
confiabilidade_bootstrap_blocos <- function(fit, dados, pasta, reamostras = 1000, bloco = 50, nivel = .95,
                                            semente = 2026, max_segundos = 45) {
  .pr_exige("lavaan", "semTools")
  t0 <- Sys.time()
  base <- confiabilidade(fit, dados, ic = FALSE)
  est <- .pr_estrutura(fit, dados); carga <- est$carga; fatores <- est$fatores; n <- nrow(dados)
  dir.create(pasta, showWarnings = FALSE, recursive = TRUE)
  arq_plano <- file.path(pasta, "plano.rds")
  impressao <- c(n = n, reamostras = reamostras, semente = semente,
                 qui2 = round(unname(lavaan::fitMeasures(fit, "chisq")), 6), itens = nrow(carga))
  if (file.exists(arq_plano)) {
    plano <- readRDS(arq_plano)
    if (!isTRUE(all.equal(plano$impressao, impressao)))
      stop("A pasta tem um plano de outro modelo, banco, semente ou numero de reamostras: usar outra pasta.", call. = FALSE)
  } else {
    if (exists(".Random.seed", envir = globalenv())) {
      velha <- get(".Random.seed", envir = globalenv()); on.exit(assign(".Random.seed", velha, envir = globalenv()), add = TRUE)
    } else on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(semente)
    idx <- vapply(seq_len(reamostras), function(b) sample.int(n, n, replace = TRUE), integer(n))
    plano <- list(impressao = impressao, idx = idx, bloco = bloco)
    saveRDS(plano, arq_plano)
  }
  cl <- fit@call; opt <- lavaan::lavInspect(fit, "options")
  cl$estimator <- opt$estimator; cl$missing <- opt$missing
  ord_names <- lavaan::lavNames(fit, "ov.ord"); if (length(ord_names)) cl$ordered <- ord_names
  if (!is.null(opt$parameterization)) cl$parameterization <- opt$parameterization
  fit_boot <- fit; fit_boot@call <- cl
  alfa_de <- function(dd) vapply(fatores, function(f) {
    x <- stats::na.omit(dd[, carga$rhs[carga$lhs == f], drop = FALSE]); x <- as.data.frame(lapply(x, as.numeric)); k <- ncol(x)
    if (k < 2) return(NA_real_); (k / (k - 1)) * (1 - sum(vapply(x, stats::var, 0)) / stats::var(rowSums(x)))
  }, numeric(1))
  uma <- function(b) {
    db <- dados[plano$idx[, b], , drop = FALSE]
    fb <- tryCatch(suppressWarnings(lavaan::update(fit_boot, data = db)), error = function(e) NULL)
    if (is.null(fb) || !isTRUE(lavaan::lavInspect(fb, "converged"))) return(rep(NA_real_, 2 * length(fatores)))
    om <- tryCatch({ rel <- semTools::compRelSEM(fb)
      vapply(fatores, function(f) if (f %in% names(rel)) as.numeric(rel[[f]]) else NA_real_, numeric(1)) },
      error = function(e) rep(NA_real_, length(fatores)))
    c(om, alfa_de(db))
  }
  inicios <- seq(1, reamostras, by = plano$bloco)
  arq_bloco <- function(i) file.path(pasta, sprintf("bloco_%04d.rds", i))
  for (i in inicios) {   # roda pelo menos um bloco por chamada; depois para quando passar de max_segundos
    if (file.exists(arq_bloco(i))) next
    bs <- i:min(reamostras, i + plano$bloco - 1)
    saveRDS(t(vapply(bs, uma, numeric(2 * length(fatores)))), arq_bloco(i))
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > max_segundos) break
  }
  feitos <- sum(file.exists(arq_bloco(inicios)))
  if (feitos < length(inicios)) {
    msg <- sprintf("Bootstrap em blocos: %d de %d blocos prontos (%d de %d reamostras). Chamar de novo com a mesma pasta para continuar.",
                   feitos, length(inicios), min(reamostras, feitos * plano$bloco), reamostras)
    message(msg)
    return(invisible(structure(list(completo = FALSE, blocos_feitos = feitos, blocos_total = length(inicios), pasta = pasta), class = "pr_boot_parcial")))
  }
  B <- do.call(rbind, lapply(inicios, function(i) readRDS(arq_bloco(i))))
  falhas <- sum(!stats::complete.cases(B)); pr <- c((1 - nivel) / 2, 1 - (1 - nivel) / 2)
  q <- apply(B, 2, stats::quantile, probs = pr, na.rm = TRUE, names = FALSE); j <- seq_along(fatores)
  out <- base
  out$omega_ic_inf <- round(q[1, j], 3); out$omega_ic_sup <- round(q[2, j], 3)
  out$alfa_ic_inf <- round(q[1, length(fatores) + j], 3); out$alfa_ic_sup <- round(q[2, length(fatores) + j], 3)
  attr(out, "bootstrap") <- list(reamostras = reamostras, falhas = falhas, nivel = nivel, semente = semente, blocos = length(inicios), pasta = pasta)
  attr(out, "metodo_ic") <- "bootstrap percentilico em blocos"
  if (falhas > 0) message(falhas, " reamostra(s) sem convergência descartada(s); se passar de 5%, o IC fica otimista.")
  out
}

# Confiabilidade de compostos definidos pelo usuário (por exemplo, escore de
# flexibilidade = soma de 12 itens). Usa o argumento W= do compRelSEM; no
# semTools 0.5.10 os argumentos higher= e omit.factors= estão obsoletos.
# A leitura depende do modelo ajustado:
#   - modelo de segunda ordem: o semTools trata só o fator de segunda ordem como
#     escore verdadeiro (os resíduos dos fatores de primeira ordem viram erro),
#     então o valor é um ômega hierárquico do composto;
#   - modelo de fatores correlacionados: todos os fatores comuns contam como
#     escore verdadeiro, então o valor é o ômega total do composto.
# Os evals conferem as duas leituras contra os valores teóricos.
confiabilidade_compostos <- function(fit, compostos, ord.scale = TRUE) {
  .pr_exige("semTools")
  stopifnot(is.list(compostos), length(names(compostos)) == length(compostos))
  W <- paste(sprintf("%s <~ %s", names(compostos),
                     vapply(compostos, paste, "", collapse = " + ")), collapse = "\n")
  r <- unlist(semTools::compRelSEM(fit, W = W, ord.scale = ord.scale))
  data.frame(composto = names(compostos), omega = round(unname(r[names(compostos)]), 3),
             itens = lengths(compostos), row.names = NULL)
}

# 4b. Escore somado: o modelo que a soma pressupõe ---------------------------
# Somar (ou tirar a média de) os itens de uma escala é pontuar com o modelo
# paralelo: cargas não padronizadas iguais (não precisam valer 1, só ser iguais)
# e variâncias residuais iguais. A AFC de validação costuma ser o modelo
# congenérico, com cargas e variâncias residuais livres; validar com um e
# pontuar com o outro deixa a soma sem evidência própria (McNeish & Wolf, 2020,
# que testam o paralelo contra o congenérico pela razão de verossimilhança).
# checar_escore_somado() ajusta, por escala ou subescala, um modelo de um fator
# (std.lv = TRUE, mesmo estimador e mesmo tratamento de ausentes de ajustar_afc)
# e compara os modelos aninhados com lavaan::lavTestLRT:
#   itens contínuos, 3 ou mais: congenérico, cargas iguais (tau-equivalente) e
#     paralelo; a comparação decisiva é paralelo contra congenérico, a do artigo;
#     as duas intermediárias mostram qual restrição falha (sequência da skill);
#   itens contínuos, 2: o congenérico não é identificado e o de cargas iguais é
#     saturado, então só a igualdade das variâncias residuais tem teste;
#   itens ordinais (WLSMV, parametrização theta), 3 ou mais: cargas iguais contra
#     livres. Na theta as variâncias residuais das respostas latentes ficam fixas
#     em 1, então cargas iguais já é o modelo paralelo dessas respostas; a soma
#     das categorias observadas ainda as trata como igualmente espaçadas, o que o
#     teste não alcança. O artigo trata os itens como contínuos (ML): a extensão
#     ao caso ordinal é da skill. Com 2 itens ordinais não há teste.
# O método do teste de diferença sai do cabeçalho do lavTestLRT: padrão com ML,
# satorra.bentler.2001 com MLR e satorra.2000 (escalonado e deslocado) com WLSMV
# no lavaan 0.7.3 (os evals conferem). Se o de 2001 der diferença negativa, a
# função refaz com satorra.bentler.2010, que a documentação do lavTestLRT indica
# para garantir valor positivo.
# Subescalas: cada uma no seu modelo de um fator. No exemplo de McNeish e Wolf
# (2020) os fatores entram juntos e sem covariância, para que o escore do
# paralelo corresponda exatamente à soma de cada subescala; a covariância entre
# subescalas fica fora desta checagem. Para um escore total, rodar de novo com
# todos os itens num vetor só. Não rejeitar o paralelo não mostra que as
# restrições valem (o teste pode ter pouco poder); a amplitude das cargas
# padronizadas do modelo menos restrito sai ao lado, sem corte (o artigo não dá
# um). Usa .ee_resultado de est_comum.R: carregar a skill por carregar.R.
# Em locale UTF-8 o lavaan põe o método entre aspas curvas (dQuote); a regra aceita
# qualquer aspa.
.pr_metodo_lrt <- function(cab) {
  m <- regmatches(cab, regexec("method = [^A-Za-z0-9._]*([A-Za-z0-9._]+)", cab))[[1]]
  if (length(m) >= 2) m[2] else "standard"
}

.pr_rotulo_lrt <- function(metodo, idioma = "pt") {
  pt <- c(standard = "teste de diferença de qui-quadrado",
          satorra.bentler.2001 = "teste de diferença de qui-quadrado escalonado",
          satorra.bentler.2010 = "teste de diferença de qui-quadrado escalonado",
          satorra.2000 = "teste de diferença de qui-quadrado escalonado e deslocado")
  en <- c(standard = "chi-square difference test",
          satorra.bentler.2001 = "scaled chi-square difference test",
          satorra.bentler.2010 = "scaled chi-square difference test",
          satorra.2000 = "scaled and shifted chi-square difference test")
  ref <- c(standard = "", satorra.bentler.2001 = "Satorra & Bentler, 2001",
           satorra.bentler.2010 = "Satorra & Bentler, 2010", satorra.2000 = "Satorra, 2000")
  rot <- if (identical(idioma, "en")) en else pt
  r <- if (metodo %in% names(rot)) rot[[metodo]] else paste0(rot[["standard"]], " (", metodo, ")")
  list(texto = r, ref = if (metodo %in% names(ref)) ref[[metodo]] else "")
}

# Ajusta capturando avisos e erro do lavaan (o aviso vai para o resultado).
.pr_ajuste_capturado <- function(f) {
  av <- character()
  r <- withCallingHandlers(tryCatch(f(), error = function(e) e),
                           warning = function(w) { av <<- c(av, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (inherits(r, "error")) return(list(fit = NULL, erro = conditionMessage(r), avisos = unique(av)))
  list(fit = r, erro = NULL, avisos = unique(av))
}

.pr_lrt_par <- function(restrito, livre) {
  roda <- function(metodo) .pr_ajuste_capturado(function() lavaan::lavTestLRT(livre, restrito, method = metodo))
  a <- roda("default")
  if (is.null(a$fit)) return(list(qui2 = NA_real_, gl = NA_real_, p = NA_real_, metodo = NA_character_,
                                  avisos = c(a$avisos, a$erro), refeito = FALSE))
  met <- .pr_metodo_lrt(attr(a$fit, "heading"))
  refeito <- FALSE
  if (identical(met, "satorra.bentler.2001") && is.finite(a$fit[2, "Chisq diff"]) && a$fit[2, "Chisq diff"] < 0) {
    b <- roda("satorra.bentler.2010")
    if (!is.null(b$fit)) { a <- b; met <- .pr_metodo_lrt(attr(b$fit, "heading")); refeito <- TRUE }
  }
  list(qui2 = unname(a$fit[2, "Chisq diff"]), gl = unname(a$fit[2, "Df diff"]), p = unname(a$fit[2, "Pr(>Chisq)"]),
       metodo = met, avisos = a$avisos, refeito = refeito)
}

checar_escore_somado <- function(dados, itens, ordinal = NULL, estimador = NULL, alfa = .05,
                                 decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  if (!exists(".ee_resultado", mode = "function")) {
    stop("checar_escore_somado() usa est_comum.R: carregar a skill por carregar.R.", call. = FALSE)
  }
  if (!is.data.frame(dados)) stop("dados precisa ser um data.frame", call. = FALSE)
  en <- identical(idioma, "en")
  nomeada <- is.list(itens)
  escalas <- if (nomeada) itens else list(escala = itens)
  if (nomeada && (is.null(names(escalas)) || any(!nzchar(names(escalas))))) {
    names(escalas) <- paste0("escala", seq_along(escalas))
  }
  escalas <- lapply(escalas, function(x) unique(as.character(x)))
  todos <- unique(unlist(escalas))
  falta <- setdiff(todos, names(dados))
  if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "), call. = FALSE)
  nao_num <- todos[!vapply(as.data.frame(dados)[todos], is.numeric, logical(1))]
  if (length(nao_num)) stop("itens precisam ser numéricos: ", paste(nao_num, collapse = ", "), call. = FALSE)
  constantes <- todos[vapply(todos, function(i) length(unique(stats::na.omit(dados[[i]]))) < 2, logical(1))]
  if (length(constantes)) stop("item sem variação (um valor só): ", paste(constantes, collapse = ", "), call. = FALSE)
  curtas <- names(escalas)[lengths(escalas) < 2]
  if (length(curtas)) stop("cada escala precisa de pelo menos 2 itens: ", paste(curtas, collapse = ", "), call. = FALSE)
  lat <- "F_soma"
  if (lat %in% todos) stop("renomear o item F_soma: o nome é usado pelo fator", call. = FALSE)
  familia_ml <- c("ML", "MLR", "MLM", "MLMV", "MLMVS", "MLF")
  if (is.null(ordinal)) {
    ordinal <- if (!is.null(estimador) && toupper(estimador) %in% familia_ml) FALSE else .pr_e_ordinal(dados, todos)
  }
  if (is.null(estimador)) estimador <- if (ordinal) "WLSMV" else "MLR"
  sep <- .ee_sep(decimal)
  rot_mod <- if (ordinal) c(congenerico = "cargas livres", tau = "cargas iguais") else
    c(congenerico = "congenérico", tau = "cargas iguais (tau-equivalente)", paralelo = "paralelo")
  mod_txt <- function(s) paste0("modelo ", if (startsWith(rot_mod[[s]], "cargas")) "de " else "", rot_mod[[s]])

  comparacoes <- list(); ajustes <- list(); cargas <- list(); amplitude <- list()
  fits_todos <- list(); avisos <- character(); apa <- character(); metodos <- character()
  decisivo <- stats::setNames(rep(NA, length(escalas)), names(escalas)); refeitos <- character()
  versoes <- character(); n_esc <- integer()

  for (nm in names(escalas)) {
    its <- escalas[[nm]]; k <- length(its)
    d <- as.data.frame(dados)[, its, drop = FALSE]
    d <- d[rowSums(!is.na(d)) > 0, , drop = FALSE]  # linha sem nenhum item respondido não entra
    livre <- paste0(lat, " =~ ", paste(its, collapse = " + "))
    iguais <- paste0(lat, " =~ ", paste0("l*", its, collapse = " + "))
    residuos <- paste0(its, " ~~ e*", its, collapse = "\n")
    specs <- list(congenerico = if (k >= 3) livre, tau = iguais,
                  paralelo = if (!ordinal) paste(iguais, residuos, sep = "\n"))
    specs <- Filter(Negate(is.null), specs)
    if (ordinal && k == 2) specs <- list()
    fits <- list()
    for (s in names(specs)) {
      a <- .pr_ajuste_capturado(function() ajustar_afc(specs[[s]], d, its, estimador = estimador,
                                                         ordinal = ordinal, std.lv = TRUE))
      rot <- paste0(mod_txt(s), if (nomeada) paste0(" (", nm, ")") else "")
      if (length(a$avisos)) avisos <- c(avisos, paste0("Aviso do lavaan (", rot, "): ", a$avisos))
      if (is.null(a$fit)) { avisos <- c(avisos, paste0("O ", rot, " não foi ajustado: ", a$erro)); next }
      if (!isTRUE(lavaan::lavInspect(a$fit, "converged"))) {
        avisos <- c(avisos, paste0("O ", rot, " não convergiu: comparações com ele ficam sem teste."))
        next
      }
      if (!isTRUE(suppressWarnings(lavaan::lavInspect(a$fit, "post.check")))) {
        avisos <- c(avisos, paste0("Solução imprópria no ", rot, " (variância negativa ou matriz não positiva definida): ler as comparações com ele com cautela."))
      }
      fits[[s]] <- a$fit
    }
    fits_todos[[nm]] <- fits
    n_esc[nm] <- if (length(fits)) as.integer(lavaan::lavInspect(fits[[1]], "nobs")) else nrow(d)

    # ajuste: a mesma versão dos índices em todos os modelos da escala com gl > 0
    # (o congenérico de 3 itens é saturado e não entra na escolha)
    if (length(fits)) {
      tabs <- lapply(names(fits), function(s) suppressWarnings(suppressMessages(tabela_ajuste(fits[[s]], rot_mod[[s]]))))
      gl_mod <- vapply(fits, function(f) unname(lavaan::fitMeasures(f, "df")), numeric(1))
      base_v <- if (any(gl_mod > 0)) tabs[gl_mod > 0] else tabs
      comum <- Reduce(intersect, lapply(base_v, `[[`, "versao"))
      v <- intersect(c("robusta", "escalonada", "padrão"), comum)[1]
      if (is.na(v)) v <- "padrão"
      versoes <- c(versoes, v)
      aj <- do.call(rbind, lapply(tabs, function(t) {
        l <- t[t$versao == v, ]
        if (!nrow(l)) { l <- t[1, ]; l[, c("cfi", "tli", "rmsea", "rmsea_li", "rmsea_ls")] <- NA; l$versao <- v }
        l
      }))
      aj <- cbind(escala = nm, aj, stringsAsFactors = FALSE)
      ajustes[[nm]] <- aj
      # cargas do modelo menos restrito que foi identificado
      base <- names(fits)[1]
      pe <- lavaan::parameterEstimates(fits[[base]], standardized = TRUE)
      ld <- pe[pe$op == "=~" & pe$lhs == lat, ]
      rs <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% its, ]
      cargas[[nm]] <- data.frame(escala = nm, modelo = rot_mod[[base]], item = ld$rhs, carga = ld$est,
                                 carga_padronizada = ld$std.all,
                                 var_residual = rs$est[match(ld$rhs, rs$lhs)], stringsAsFactors = FALSE)
      amplitude[[nm]] <- data.frame(escala = nm, modelo = rot_mod[[base]], itens = k,
                                    carga_pad_min = min(ld$std.all), carga_pad_max = max(ld$std.all),
                                    amplitude = max(ld$std.all) - min(ld$std.all), stringsAsFactors = FALSE)
    }

    # comparações
    pares <- if (ordinal) {
      if (k >= 3) list(c("tau", "congenerico", "cargas iguais", "sim")) else list()
    } else if (k >= 3) {
      list(c("tau", "congenerico", "cargas iguais", "não"),
           c("paralelo", "tau", "variâncias residuais iguais, dadas cargas iguais", "não"),
           c("paralelo", "congenerico", "modelo paralelo (cargas e variâncias residuais iguais)", "sim"))
    } else list(c("paralelo", "tau", "variâncias residuais iguais (com 2 itens, a única restrição com teste)", "sim"))
    for (pr in pares) {
      ok <- all(pr[1:2] %in% names(fits))
      t <- if (ok) .pr_lrt_par(fits[[pr[1]]], fits[[pr[2]]]) else
        list(qui2 = NA_real_, gl = NA_real_, p = NA_real_, metodo = NA_character_, avisos = character(), refeito = FALSE)
      if (length(t$avisos)) avisos <- c(avisos, paste0("Aviso do lavTestLRT (", nm, ", ", rot_mod[[pr[1]]], " contra ",
                                                      rot_mod[[pr[2]]], "): ", t$avisos))
      if (isTRUE(t$refeito)) refeitos <- c(refeitos, nm)
      if (!is.na(t$metodo)) metodos <- c(metodos, t$metodo)
      comparacoes[[length(comparacoes) + 1]] <- data.frame(
        escala = nm, n = if (nm %in% names(n_esc)) n_esc[[nm]] else NA_integer_,
        restrito = rot_mod[[pr[1]]], livre = rot_mod[[pr[2]]], restricao = pr[3], decisiva = pr[4] == "sim",
        qui2_dif = t$qui2, gl_dif = t$gl, p = t$p, metodo = t$metodo,
        resultado = if (is.na(t$p)) "sem teste" else if (t$p < alfa) "restrição rejeitada" else "restrição não rejeitada",
        stringsAsFactors = FALSE)
      if (pr[4] == "sim" && !is.na(t$p)) decisivo[[nm]] <- t$p < alfa
    }

    # frase APA da escala
    em <- if (nomeada) (if (en) paste0("For ", nm, ", ") else paste0("Em ", nm, ", ")) else ""
    cap <- function(s) if (nzchar(em)) s else paste0(toupper(substring(s, 1, 1)), substring(s, 2))
    if (ordinal && k == 2) {
      apa <- c(apa, if (en) paste0(em, cap("with two ordinal items the sum-score constraints cannot be tested: the equal-loadings model is already saturated.")) else
        paste0(em, cap("com dois itens ordinais as restrições da soma não têm teste: o modelo de cargas iguais já é saturado.")))
      next
    }
    cmp <- comparacoes[[length(comparacoes)]]
    if (is.na(cmp$p)) {
      apa <- c(apa, if (en) paste0(em, cap("the sum-score constraints could not be tested (see warnings).")) else
        paste0(em, cap("as restrições da soma ficaram sem teste (ver avisos).")))
      next
    }
    rl <- .pr_rotulo_lrt(cmp$metodo, idioma)
    teste_txt <- paste0("(", rl$texto, if (nzchar(rl$ref)) paste0("; ", rl$ref) else "", ")")
    dq <- paste0("Δχ²(", .ee_gl(cmp$gl_dif, decimal), ") = ", .pr_num(cmp$qui2_dif, 2, TRUE, decimal), sep, .pr_p(cmp$p, decimal))
    rest <- cmp$restrito
    linha <- ajustes[[nm]][ajustes[[nm]]$modelo == rest, ]
    aj_txt <- if (nrow(linha) == 1) formatar_ajuste_apa(linha, decimal) else ""
    if (en) aj_txt <- sub("IC 90%", "90% CI", aj_txt, fixed = TRUE)
    pior <- cmp$p < alfa
    frase <- if (ordinal) {
      if (en) paste0(em, cap("with ordinal items (WLSMV, theta parameterization), the equal-loadings model "),
                     if (pior) "fit worse than the free-loadings model " else "did not fit significantly worse than the free-loadings model ",
                     teste_txt, ", ", dq, if (nzchar(aj_txt)) paste0("; equal-loadings model fit: ", aj_txt) else "", ".") else
        paste0(em, cap("com itens ordinais (WLSMV, parametrização theta), o modelo de cargas iguais "),
               if (pior) "ajustou-se pior que o de cargas livres " else "não se ajustou significativamente pior que o de cargas livres ",
               teste_txt, ", ", dq, if (nzchar(aj_txt)) paste0("; ajuste do modelo de cargas iguais: ", aj_txt) else "", ".")
    } else if (k >= 3) {
      if (en) paste0(em, cap("the parallel model (equal loadings and residual variances), which corresponds to summing the items (McNeish & Wolf, 2020), "),
                     if (pior) "fit worse than the congeneric model " else "did not fit significantly worse than the congeneric model ",
                     teste_txt, ", ", dq, if (nzchar(aj_txt)) paste0("; parallel model fit: ", aj_txt) else "", ".") else
        paste0(em, cap("o modelo paralelo (cargas e variâncias residuais iguais), que corresponde à soma dos itens (McNeish & Wolf, 2020), "),
               if (pior) "ajustou-se pior que o congenérico " else "não se ajustou significativamente pior que o congenérico ",
               teste_txt, ", ", dq, if (nzchar(aj_txt)) paste0("; ajuste do paralelo: ", aj_txt) else "", ".")
    } else {
      if (en) paste0(em, cap("with two items, the parallel model "),
                     if (pior) "fit worse than the equal-loadings model with free residual variances " else
                       "did not fit significantly worse than the equal-loadings model with free residual variances ",
                     teste_txt, ", ", dq, "; with two items, equal loadings cannot be tested.") else
        paste0(em, cap("com dois itens, o modelo paralelo "),
               if (pior) "ajustou-se pior que o de cargas iguais com variâncias residuais livres " else
                 "não se ajustou significativamente pior que o de cargas iguais com variâncias residuais livres ",
               teste_txt, ", ", dq, "; com dois itens, a igualdade das cargas não tem teste.")
    }
    apa <- c(apa, frase)
  }

  tab <- if (length(comparacoes)) do.call(rbind, comparacoes) else NULL
  aj <- if (length(ajustes)) do.call(rbind, ajustes) else NULL
  cg <- if (length(cargas)) do.call(rbind, cargas) else NULL
  amp <- if (length(amplitude)) do.call(rbind, amplitude) else NULL
  if (!is.null(tab)) rownames(tab) <- NULL
  if (!is.null(aj)) rownames(aj) <- NULL
  if (!is.null(cg)) rownames(cg) <- NULL
  if (!is.null(amp)) rownames(amp) <- NULL

  metodos <- unique(metodos)
  aus <- .pr_ausentes(ordinal, estimador)
  decisoes <- c(
    "Somar ou tirar a média dos itens é pontuar com o modelo paralelo: cargas não padronizadas iguais (não precisam valer 1) e variâncias residuais iguais. A AFC usual é o congenérico, com cargas e variâncias residuais livres; validar com o congenérico e pontuar com a soma deixa a soma sem evidência própria (McNeish & Wolf, 2020).",
    paste0("Um modelo de um fator por escala (std.lv = TRUE), estimador ", estimador, if (ordinal) " com itens ordinais e parametrização theta" else "",
           ", ausentes por ", c(fiml = "FIML", pairwise = "deleção pairwise", listwise = "deleção listwise")[[aus]],
           "; teste de diferença pelo lavTestLRT, método ", if (length(metodos)) paste(metodos, collapse = " e ") else "nenhum", "."),
    if (ordinal) "Itens ordinais: só cargas iguais contra livres. Na parametrização theta as variâncias residuais das respostas latentes ficam fixas em 1, então cargas iguais já é o modelo paralelo dessas respostas; a soma das categorias observadas ainda as trata como igualmente espaçadas, o que o teste não alcança. McNeish e Wolf (2020) tratam os itens como contínuos: a extensão ao caso ordinal é da skill." else
      "Comparação decisiva: paralelo contra congenérico, a mesma de McNeish e Wolf (2020). As intermediárias (cargas iguais contra congenérico; paralelo contra cargas iguais) mostram qual restrição falha; essa sequência é da skill.",
    if (length(versoes)) paste0("Índices de ajuste na versão ", paste(unique(versoes), collapse = " e "), ", a mesma em todos os modelos de cada escala."),
    if (length(refeitos)) paste0("Em ", paste(unique(refeitos), collapse = ", "), " o método satorra.bentler.2001 deu diferença negativa; o teste foi refeito com satorra.bentler.2010, que a documentação do lavTestLRT indica para garantir valor positivo."),
    if (length(escalas) > 1) "Cada subescala foi testada no seu modelo de um fator. No exemplo de McNeish e Wolf (2020) os fatores entram juntos e sem covariância, para que o escore do paralelo corresponda à soma de cada subescala; a covariância entre subescalas fica fora desta checagem. Para um escore total, rodar de novo com todos os itens num vetor só.",
    "Correlação alta entre a soma e o escore fatorial não substitui o teste: no exemplo de McNeish e Wolf (2020), R² = .970 entre os dois e a conclusão sobre a diferença entre grupos mudou."
  )
  rejeitadas <- names(decisivo)[decisivo %in% TRUE]
  mantidas <- names(decisivo)[decisivo %in% FALSE]
  onde <- function(v) if (nomeada) paste0(" em ", paste(v, collapse = ", ")) else ""
  if (length(rejeitadas)) decisoes <- c(decisoes, paste0(
    "Restrições da soma rejeitadas", onde(rejeitadas), ": em vez da soma, usar o modelo de mensuração dentro de uma SEM ",
    "(abordagem simultânea) ou escore fatorial com a correção de Croon; para o escore fatorial, Bartlett quando ele for desfecho e ",
    "regressão (MAP) quando for preditor, como McNeish e Wolf (2020) resumem de Skrondal e Laake (2001)."))
  if (length(mantidas)) avisos <- c(avisos, paste0(
    "Restrições da soma não rejeitadas", onde(mantidas), ": isso não mostra que elas valem, porque o teste pode ter tido pouco poder ",
    "(Greenland et al., 2016). Conferir o ajuste do ", if (ordinal) "modelo de cargas iguais" else "paralelo",
    " e a amplitude das cargas padronizadas", if (ordinal) "." else
      "; um congenérico com cargas padronizadas quase iguais pode ser aproximado pelo paralelo (McNeish & Wolf, 2020)."))
  k2 <- names(escalas)[lengths(escalas) == 2]
  if (length(k2)) avisos <- c(avisos, paste0(
    if (nomeada) paste0("Em ", paste(k2, collapse = ", "), ", com 2 itens") else "Com 2 itens", ", ",
    if (ordinal) "nenhuma restrição da soma tem teste: o modelo de cargas iguais já é saturado." else
      "o congenérico não é identificado e o modelo de cargas iguais é saturado: a igualdade das cargas não tem teste, só a das variâncias residuais."))
  .ee_resultado(
    "Checagem do escore somado (modelo paralelo contra congenérico)",
    n = if (length(n_esc)) min(n_esc) else 0L, n_excluidos = if (length(n_esc)) nrow(dados) - min(n_esc) else nrow(dados),
    tabela = tab, decisoes = decisoes, avisos = unique(avisos), apa = apa,
    extra = list(tabelas_extra = list(`Ajuste dos modelos` = aj, `Cargas do modelo menos restrito` = cg,
                                      `Amplitude das cargas padronizadas` = amp),
                 comparacoes = tab, ajuste = aj, cargas = cg, amplitude = amp, fits = fits_todos,
                 soma_rejeitada = decisivo, metodo_lrt = metodos, ordinal = ordinal, estimador = estimador,
                 n_por_escala = n_esc))
}

# 4c. Modelo bifator: índices além do ajuste ---------------------------------
# Rodriguez, Reise e Haviland (2016) e Reise, Bonifay e Haviland (2018, cap. 22 do
# Wiley Handbook): no bifator o ajuste quase sempre favorece o modelo (ele absorve
# mais variância), então a decisão entre escore total e subescalas sai dos índices.
# Cargas padronizadas do lavaan (lavInspect "std"); fatores ortogonais obrigatórios.
#   ômega total (ωT) e hierárquico (ωH) do escore total, e por subescala ωS e ωHS
#   (Eq. 22.3 e 22.4 do capítulo), com a variância do escore somado implicada pelo
#   modelo (1'Σ1, inclui correlação residual especificada);
#   ECV (variância comum devida ao geral), ECV de cada específico dentro da
#   subescala (ECV_SS) e em relação ao total (ECV_SG), I-ECV por item;
#   PUC = 1 - soma das correlações dentro dos grupos / total de correlações;
#   H de Hancock e Mueller (2001) = 1 / (1 + 1 / Σ[λ²/(1 - λ²)]) por fator;
#   determinação fatorial FD = raiz da diagonal de Φ Λ' Σ^-1 Λ Φ (correlação entre o
#   fator e o escore fatorial por regressão);
#   ARPB (viés relativo médio das cargas do modelo de um fator em relação às cargas
#   gerais), quando o ajuste unidimensional entra em `unidimensional`.
# Itens ordinais (WLSMV): as cargas são das respostas latentes, então ω e ωH aqui são
# da soma das respostas latentes; o ω da soma observada das categorias (Green & Yang,
# 2009) sai do BifactorIndicesCalculator::bifactorIndices, que os evals usam para
# conferir o caso contínuo. Réguas declaradas como réguas, com a fonte na saída.
indices_bifator <- function(fit, geral = NULL, unidimensional = NULL, decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  if (lavaan::lavInspect(fit, "ngroups") > 1) stop("indices_bifator() é para modelo de um grupo.", call. = FALSE)
  st <- lavaan::lavInspect(fit, "std")
  L <- st$lambda; Psi <- st$psi; Th <- st$theta
  ov <- rownames(L); fat <- colnames(L)
  nz <- abs(L) > 1e-12
  if (is.null(geral)) {
    cand <- fat[colSums(nz) == nrow(L)]
    if (length(cand) != 1) stop("não achei um fator único com carga em todos os itens: informar `geral` (nome do fator geral).", call. = FALSE)
    geral <- cand
  }
  if (!geral %in% fat) stop("fator geral '", geral, "' não está no modelo", call. = FALSE)
  esp <- setdiff(fat, geral)
  if (!length(esp)) stop("o modelo não tem fatores específicos", call. = FALSE)
  if (any(!nz[, geral])) stop("há itens sem carga no fator geral: ", paste(ov[!nz[, geral]], collapse = ", "), call. = FALSE)
  fora <- Psi; diag(fora) <- 0
  if (max(abs(fora)) > 1e-6) stop("os fatores estão correlacionados: o bifator pede fatores ortogonais (orthogonal = TRUE no cfa).", call. = FALSE)
  if (any(rowSums(nz[, esp, drop = FALSE]) > 1)) stop("item com carga em mais de um fator específico: estrutura fora do bifator.", call. = FALSE)
  grupo <- apply(nz[, esp, drop = FALSE], 1, function(z) if (any(z)) esp[z] else NA_character_)
  lg <- L[, geral]; ls <- rowSums(L[, esp, drop = FALSE])
  Sig <- L %*% Psi %*% t(L) + Th
  den <- sum(Sig)
  cs <- c(sum(lg), vapply(esp, function(f) sum(L[, f]), 0))
  omega_t <- sum(cs^2) / den; omega_h <- cs[1]^2 / den
  ecv <- sum(lg^2) / (sum(lg^2) + sum(L[, esp]^2))
  p <- length(ov); ns <- vapply(esp, function(f) sum(grupo == f, na.rm = TRUE), 0)
  puc <- 1 - sum(ns * (ns - 1) / 2) / (p * (p - 1) / 2)
  hh <- function(l) { l <- l[abs(l) > 1e-12]; s <- sum(l^2 / (1 - l^2)); 1 / (1 + 1 / s) }
  Sinv <- solve(Sig)
  fd <- sqrt(diag(Psi %*% t(L) %*% Sinv %*% L %*% Psi)); names(fd) <- fat
  por_esp <- do.call(rbind, lapply(esp, function(f) {
    it <- which(grupo == f); S <- Sig[it, it, drop = FALSE]; d <- sum(S)
    data.frame(fator = f, itens = length(it), omega_s = (sum(lg[it])^2 + sum(L[it, f])^2) / d,
               omega_hs = sum(L[it, f])^2 / d, ecv_ss = sum(L[it, f]^2) / (sum(lg[it]^2) + sum(L[it, f]^2)),
               ecv_sg = sum(L[it, f]^2) / (sum(lg^2) + sum(L[, esp]^2)), h = hh(L[it, f]), fd = fd[[f]],
               stringsAsFactors = FALSE)
  }))
  tab <- rbind(data.frame(fator = geral, itens = p, omega_s = omega_t, omega_hs = omega_h, ecv_ss = ecv, ecv_sg = ecv,
                          h = hh(lg), fd = fd[[geral]], stringsAsFactors = FALSE), por_esp)
  names(tab)[names(tab) == "omega_s"] <- "omega"; names(tab)[names(tab) == "omega_hs"] <- "omega_h"
  names(tab)[names(tab) == "ecv_ss"] <- "ecv"; names(tab)[names(tab) == "ecv_sg"] <- "ecv_total"
  itens <- data.frame(item = ov, especifico = grupo, carga_geral = unname(lg), carga_especifica = unname(ls),
                      i_ecv = unname(lg^2 / (lg^2 + ls^2)), stringsAsFactors = FALSE)
  arpb <- NA_real_
  if (!is.null(unidimensional)) {
    lu <- lavaan::lavInspect(unidimensional, "std")$lambda
    if (ncol(lu) != 1) stop("`unidimensional` deve ser o ajuste do modelo de um fator", call. = FALSE)
    lu <- lu[ov, 1]
    itens$carga_unidimensional <- unname(lu)
    itens$vies_relativo <- unname((lu - lg) / lg)
    arpb <- mean(abs(itens$vies_relativo))
  }
  ordinal <- length(lavaan::lavNames(fit, "ov.ord")) > 0
  ef <- data.frame(indice = c("ômega hierárquico (ωH)", "ômega total (ωT)", "ECV", "PUC", "H do fator geral", "FD do fator geral",
                              if (!is.na(arpb)) "ARPB"),
                   valor = c(omega_h, omega_t, ecv, puc, hh(lg), fd[[geral]], if (!is.na(arpb)) arpb),
                   ic_inf = NA_real_, ic_sup = NA_real_, principal = c(TRUE, rep(FALSE, 5 + !is.na(arpb))), stringsAsFactors = FALSE)
  f2 <- function(x) .pr_num(x, 2, FALSE, decimal)
  decisoes <- c(paste0("Fator geral: ", geral, "; específicos: ", paste(esp, collapse = ", "), ". Cargas padronizadas; variância do escore total ",
                       "implicada pelo modelo."),
                paste0("Réguas declaradas: ωH ≥ ", f2(.80), " para interpretar o escore total como medida do fator geral (Rodriguez et al., 2016); ",
                       "ECV > ", f2(.70), " com PUC > ", f2(.70), " para tratar a variância comum como essencialmente unidimensional ",
                       "(Rodriguez et al., 2016); com PUC < ", f2(.80), ", ECV > ", f2(.60), " e ωH > ", f2(.70), " (Reise et al., 2013); ",
                       "H ≥ ", f2(.70), " para fator bem definido (Hancock & Mueller, 2001); FD ≥ ", f2(.90), " para usar escore fatorial ",
                       "(Gorsuch, 1983); ARPB < 10 a 15% (Muthén, Kaplan & Hollis, 1987). O curso usa ECV > ", f2(.85),
                       " no critério de unidimensionalidade do FACTOR (Ferrando & Lorenzo-Seva, 2018)."),
                "ωHS sem régua fixa: compare com o ωS da mesma subescala; a diferença é a confiabilidade que a subescala toma emprestada do fator geral (Reise, Bonifay & Haviland, 2018).")
  if (ordinal) decisoes <- c(decisoes, "Itens ordinais: ω e ωH aqui são da soma das respostas latentes. Para a soma observada das categorias (Green & Yang, 2009), BifactorIndicesCalculator::bifactorIndices(fit).")
  avisos <- character()
  if (omega_h < .5) avisos <- c(avisos, paste0("ωH = ", f2(omega_h), ": menos da metade da variância do escore total vem do fator geral; o escore total não representa bem o geral."))
  fracos <- por_esp$fator[por_esp$h < .70]
  if (length(fracos)) avisos <- c(avisos, paste0("Específico(s) com H < ", f2(.70), " (", paste(fracos, collapse = ", "),
                                                 "): fator mal definido pelos itens, pouco replicável; escore dessa subescala não se sustenta como medida própria."))
  decisoes <- c(decisoes, "I-ECV por item na tabela de itens: perto de 1, o item mede sobretudo o fator geral; perto de 0, sobretudo o específico.")
  essencial <- (ecv > .70 && puc > .70) || (puc < .80 && ecv > .60 && omega_h > .70)
  apa <- if (identical(idioma, "en")) {
    paste0("In the bifactor model, the general factor accounted for ", .pr_num(100 * ecv, 1, TRUE, decimal), "% of the common variance (ECV = ", f2(ecv),
           "; PUC = ", f2(puc), "). Omega hierarchical for the total score was ", f2(omega_h), " (omega total = ", f2(omega_t), "); H = ", f2(hh(lg)),
           " and factor determinacy = ", f2(fd[[geral]]), " for the general factor. Omega hierarchical subscale values ranged from ",
           f2(min(por_esp$omega_hs)), " to ", f2(max(por_esp$omega_hs)), " (Rodriguez et al., 2016).",
           if (essencial) " These values support treating the scale as essentially unidimensional." else "")
  } else {
    paste0("No modelo bifator, o fator geral respondeu por ", .pr_num(100 * ecv, 1, TRUE, decimal), "% da variância comum (ECV = ", f2(ecv),
           "; PUC = ", f2(puc), "). O ômega hierárquico do escore total foi ", f2(omega_h), " (ômega total = ", f2(omega_t), "); H = ", f2(hh(lg)),
           " e determinação fatorial = ", f2(fd[[geral]]), " para o fator geral. O ômega hierárquico das subescalas variou de ",
           f2(min(por_esp$omega_hs)), " a ", f2(max(por_esp$omega_hs)), " (Rodriguez et al., 2016).",
           if (essencial) " Esses valores sustentam tratar a escala como essencialmente unidimensional." else "")
  }
  .ee_resultado("Índices do modelo bifator", n = lavaan::lavInspect(fit, "nobs"), tabela = tab, efeito = ef,
                decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(itens = itens, arpb = arpb, essencialmente_unidimensional = essencial,
                             tabelas_extra = list(`Itens (cargas e I-ECV)` = itens)))
}

# 5. Validade convergente e discriminante -------------------------------------
# AVE por fator, correlação entre fatores, HTMT (Henseler et al., 2015) e o
# critério de Fornell-Larcker, só entre fatores de primeira ordem (com
# indicadores observados), alinhados por nome. Com menos de 2 fatores devolve
# NA no lugar do critério. Em modelo de segunda ordem, o HTMT é calculado sobre
# a parte de primeira ordem do modelo.
validade_discriminante <- function(fit, dados) {
  .pr_exige("lavaan", "semTools")
  if (lavaan::lavInspect(fit, "ngroups") > 1) {
    stop("validade_discriminante() é para modelo de um grupo: ajustar e rodar por grupo.", call. = FALSE)
  }
  est <- .pr_estrutura(fit, dados)
  fatores <- est$fatores
  ave <- semTools::AVE(fit)
  ave <- ave[intersect(fatores, names(ave))]
  if (length(fatores) < 2) {
    message("Validade discriminante pede pelo menos 2 fatores de primeira ordem.")
    return(list(ave = round(ave, 3), raiz_ave = round(sqrt(ave), 3), cor_fatores = NA,
                htmt = NA, fornell_larcker_ok = NA, segunda_ordem = est$superior))
  }
  cor_lv <- lavaan::lavInspect(fit, "cor.lv")[fatores, fatores]
  sintaxe <- paste(vapply(fatores, function(f)
    sprintf("%s =~ %s", f, paste(est$carga$rhs[est$carga$lhs == f], collapse = " + ")), ""),
    collapse = "\n")
  ordenados <- lavaan::lavNames(fit, "ov.ord")
  ht <- semTools::htmt(sintaxe, data = dados,
                       ordered = if (length(ordenados)) intersect(est$carga$rhs, ordenados) else NULL)
  ht <- ht[fatores, fatores]
  raiz <- sqrt(ave[fatores])
  lim <- outer(raiz, raiz, pmin)
  sup <- upper.tri(cor_lv)
  list(ave = round(ave[fatores], 3), raiz_ave = round(raiz, 3),
       cor_fatores = round(cor_lv, 3), htmt = round(ht, 3),
       fornell_larcker_ok = all(lim[sup] > abs(cor_lv[sup])),
       segunda_ordem = est$superior)
}

# 6. Invariância entre grupos -------------------------------------------------
# Ordinal: configural -> limiares -> cargas -> interceptos, com a identificação
# de Wu e Estabrook (2016) via semTools::measEq.syntax. Igualar os limiares
# primeiro libera restrições de identificação, a começar pelo intercepto da
# resposta latente de cada item nos grupos além do primeiro (documentação do
# measEq.syntax). Consequência prática, mostrada nos evals: deslocamento
# uniforme de todos os limiares de um item aparece no passo dos interceptos;
# o passo dos limiares acusa mudança no espaçamento entre eles. Com itens de 2
# ou 3 categorias, a documentação diz que não há teste independente dos
# limiares; o passo fica sem p quando não ganha grau de liberdade.
# Contínuo: configural -> métrica -> escalar -> estrita.
# Modelo com fator de ordem superior não é aceito: o measEq.syntax troca a
# identificação para ID.fac = "ul" e os passos daqui não alcançam o nível de
# segunda ordem (ver references/invariancia.md).
invariancia <- function(modelo, dados, itens, grupo, ordinal = NULL, estimador = NULL) {
  .pr_exige("lavaan", "semTools")
  pt0 <- lavaan::lavaanify(modelo)
  latentes <- unique(pt0$lhs[pt0$op == "=~"])
  if (any(pt0$op == "=~" & pt0$rhs %in% latentes)) {
    stop("invariancia() não cobre modelo com fator de ordem superior: o measEq.syntax troca a ",
         "identificação para ID.fac = 'ul' e os passos de cargas e interceptos não alcançam o ",
         "nível de segunda ordem. Testar a medida no modelo de fatores correlacionados de primeira ",
         "ordem; o nível de segunda ordem pede passos com 'regressions' e 'means' montados à mão ",
         "(references/invariancia.md).", call. = FALSE)
  }
  if (is.null(ordinal)) ordinal <- .pr_e_ordinal(dados, itens)
  if (is.null(estimador)) estimador <- if (ordinal) "WLSMV" else "MLR"
  if (ordinal) {
    kmin <- min(vapply(itens, function(i) length(unique(stats::na.omit(dados[[i]]))), numeric(1)))
    if (kmin <= 3) {
      message("Item com ", kmin, " categorias: pela documentação do measEq.syntax, com item binário ",
              "não há teste independente dos limiares e com 3 categorias a igualdade dos limiares é ",
              "pressuposta; ler o passo 'limiares' com isso em mente.")
    }
  }
  passos <- if (ordinal) {
    list(configural = "", limiares = "thresholds",
         cargas = c("thresholds", "loadings"),
         interceptos = c("thresholds", "loadings", "intercepts"))
  } else {
    list(configural = "", metrica = "loadings",
         escalar = c("loadings", "intercepts"),
         estrita = c("loadings", "intercepts", "residuals"))
  }
  fits <- lapply(passos, function(eq) {
    args <- list(configural.model = modelo, data = dados, group = grupo,
                 ID.fac = "std.lv", estimator = estimador, return.fit = TRUE,
                 missing = .pr_ausentes(ordinal, estimador))
    if (ordinal) {
      args <- c(args, list(ordered = itens, ID.cat = "Wu.Estabrook.2016",
                           parameterization = "theta"))
    }
    if (any(nzchar(eq))) args$group.equal <- eq
    do.call(semTools::measEq.syntax, args)
  })
  tabs <- Map(tabela_ajuste, fits, names(fits))
  # mesma versão dos índices em todos os passos: robusta se existir em todos,
  # senão escalonada, senão padrão
  comum <- Reduce(intersect, lapply(tabs, `[[`, "versao"))
  versao <- intersect(c("robusta", "escalonada", "padrão"), comum)[1]
  tab <- do.call(rbind, lapply(tabs, function(t) t[t$versao == versao, ]))
  tab$delta_cfi <- c(NA, diff(tab$cfi))
  tab$delta_rmsea <- c(NA, diff(tab$rmsea))
  tab$delta_srmr <- c(NA, diff(tab$srmr))
  gl <- vapply(fits, function(f) unname(lavaan::fitMeasures(f, "df")), numeric(1))
  tab$delta_gl <- c(NA, diff(gl))
  tab$p_dif_qui2 <- c(NA, vapply(seq_along(fits)[-1], function(i) {
    if (gl[i] == gl[i - 1]) return(NA_real_)  # passo sem grau de liberdade: sem teste
    tryCatch(lavaan::lavTestLRT(fits[[i - 1]], fits[[i]])[2, "Pr(>Chisq)"],
             error = function(e) NA_real_)
  }, numeric(1)))
  rownames(tab) <- NULL
  list(tabela = tab, fits = fits)
}

# 7. Relato no formato APA 7 --------------------------------------------------
# Sem zero à esquerda em p, CFI, TLI, RMSEA e SRMR, como na fabrica-de-artigo
# (o TLI pode passar de 1; a convenção de relato dele é a mesma do CFI). Com
# decimal = "," a vírgula vira separador decimal e os itens da frase passam a
# ser separados por ponto e vírgula; o zero à esquerda com vírgula fica a
# critério do periódico (zero_inicial = TRUE liga): conferir a norma.
.pr_num <- function(x, casas, zero = TRUE, decimal = ".") {
  if (is.na(x)) return("NA")
  # Meio para cima (afastando do zero), como o JASP e o SPSS imprimem: round() do R
  # arredonda 155.045 para 155.04 por causa da representação binária.
  x <- sign(x) * floor(abs(x) * 10^casas + 0.5 + 1e-9) / 10^casas
  if (x == 0) x <- 0  # tira o sinal de -0
  s <- formatC(x, format = "f", digits = casas)
  if (!zero) s <- sub("^(-?)0\\.", "\\1.", s)
  if (decimal != ".") s <- chartr(".", decimal, s)
  s
}

.pr_p <- function(p, decimal = ".", zero = FALSE) {
  if (is.na(p)) return("p = NA")
  if (p < .001) return(paste0("p < ", .pr_num(.001, 3, zero, decimal)))
  if (p > .999) return(paste0("p > ", .pr_num(.999, 3, zero, decimal)))
  paste0("p = ", .pr_num(p, 3, zero, decimal))
}

formatar_ajuste_apa <- function(linha, decimal = ".", zero_inicial = FALSE) {
  stopifnot(nrow(linha) == 1)
  z <- zero_inicial
  sep <- if (decimal == ",") "; " else ", "
  paste0("χ²(", formatC(linha$gl, format = "d"), ") = ", .pr_num(linha$qui2, 2, TRUE, decimal),
         sep, .pr_p(linha$p, decimal, z),
         sep, "CFI = ", .pr_num(linha$cfi, 3, z, decimal),
         sep, "TLI = ", .pr_num(linha$tli, 3, z, decimal),
         sep, "RMSEA = ", .pr_num(linha$rmsea, 3, z, decimal),
         sep, "IC 90% [", .pr_num(linha$rmsea_li, 3, z, decimal), sep,
         .pr_num(linha$rmsea_ls, 3, z, decimal), "]",
         sep, "SRMR = ", .pr_num(linha$srmr, 3, z, decimal))
}

# Versões dos pacotes, para a seção de análise de dados e para o LOG.
versoes_pacotes <- function(pacotes = NULL) {
  # Sem lista: os pacotes da psicometria mais os dos módulos de estatística carregados.
  if (is.null(pacotes)) {
    pacotes <- c("lavaan", "semTools", "mirt", "bootnet", "qgraph", "EGAnet",
                 "NetworkComparisonTest", "Boruta", "psych")
    if (exists(".ee_pacotes_estatistica")) pacotes <- unique(c(pacotes, get(".ee_pacotes_estatistica")))
  }
  ok <- pacotes[vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
  c(R = paste(R.version$major, R.version$minor, sep = "."),
    vapply(ok, function(p) as.character(utils::packageVersion(p)), ""))
}

# dMACS (v1.6): tamanho do efeito da não invariância por item (Nye & Drasgow, 2011; Wiley,
# cap. 27; Petersen, 2024, 15.10.3.4). Para indicadores contínuos o integral tem forma
# fechada: com Δ(η) = (ν_F - ν_R) + (λ_F - λ_R)'η e η ~ N(κ_F, Ψ_F) no grupo focal,
#   dMACS = sqrt[(Δν + Δλ'κ_F)² + Δλ'Ψ_F Δλ] / DP combinado do item,
# DP combinado = sqrt(((n_R - 1)s²_R + (n_F - 1)s²_F)/(n_R + n_F - 2)). A versão com sinal,
# Δν + Δλ'κ_F, é a diferença esperada no escore do item (focal menos referência) para
# pessoas com o mesmo traço; somada nos itens do fator é o ΔMean da escala. Rodar no
# modelo em que os parâmetros que diferem estão livres (configural ou invariância parcial
# saída de invariancia()$fits), com a métrica latente comparável entre os grupos.
# Régua citada por Nye et al. (2019): ,20 pequeno, ,40 médio, ,70 grande.
dmacs <- function(fit, grupo_ref = NULL, decimal = ".", idioma = "pt") {
  .pr_exige("lavaan")
  ng <- lavaan::lavInspect(fit, "ngroups")
  if (ng < 2) stop("dmacs() pede um modelo de grupos múltiplos (lavaan com group =)", call. = FALSE)
  if (length(lavaan::lavNames(fit, "ov.ord")))
    stop("itens ordinais: o dMACS daqui é para indicadores contínuos. Para itens ordinais, tri_dif() ou pinsearch::pin_effsize() no modelo do lavaan.", call. = FALSE)
  rot <- lavaan::lavInspect(fit, "group.label")
  ref <- if (is.null(grupo_ref)) 1L else if (is.character(grupo_ref)) match(grupo_ref, rot) else as.integer(grupo_ref)
  if (is.na(ref)) stop("grupo_ref não encontrado entre ", paste(rot, collapse = ", "), call. = FALSE)
  est <- lavaan::lavInspect(fit, "est"); dat <- lavaan::lavInspect(fit, "data")
  if (!is.list(dat) || is.matrix(dat)) dat <- list(dat)
  itens <- rownames(est[[1]]$lambda); fatores <- colnames(est[[1]]$lambda)
  dp_ref <- apply(dat[[ref]], 2, stats::sd, na.rm = TRUE); n_ref <- colSums(!is.na(dat[[ref]]))
  linhas <- list(); escalas <- list(); metrica_fixa <- FALSE
  for (g in setdiff(seq_len(ng), ref)) {
    eR <- est[[ref]]; eF <- est[[g]]
    nuR <- if (is.null(eR$nu)) rep(0, length(itens)) else as.numeric(eR$nu)
    nuF <- if (is.null(eF$nu)) rep(0, length(itens)) else as.numeric(eF$nu)
    kF <- if (is.null(eF$alpha)) rep(0, length(fatores)) else as.numeric(eF$alpha)
    PsiF <- eF$psi
    if (all(abs(kF) < 1e-8) && all(abs(diag(PsiF) - 1) < 1e-8)) metrica_fixa <- TRUE
    dp_f <- apply(dat[[g]], 2, stats::sd, na.rm = TRUE); n_f <- colSums(!is.na(dat[[g]]))
    dp_c <- sqrt(((n_ref - 1) * dp_ref^2 + (n_f - 1) * dp_f^2) / (n_ref + n_f - 2))
    dp_c <- dp_c[itens]
    for (i in seq_along(itens)) {
      dl <- eF$lambda[i, ] - eR$lambda[i, ]; dn <- nuF[i] - nuR[i]
      media <- dn + sum(dl * kF)
      quad <- media^2 + as.numeric(t(dl) %*% PsiF %*% dl)
      f_item <- fatores[which(eR$lambda[i, ] != 0 | eF$lambda[i, ] != 0)]
      linhas[[length(linhas) + 1]] <- data.frame(grupo_focal = rot[g], fator = paste(f_item, collapse = "+"), item = itens[i],
                                                 dmacs = sqrt(quad) / dp_c[[i]], delta_media = media, dp_combinado = dp_c[[i]])
    }
  }
  tab <- do.call(rbind, linhas); rownames(tab) <- NULL
  tab$classe <- cut(tab$dmacs, c(-Inf, .2, .4, .7, Inf), right = FALSE, labels = c("< ,20", "pequeno", "médio", "grande"))
  esc <- stats::aggregate(delta_media ~ grupo_focal + fator, data = tab, FUN = sum)
  names(esc)[3] <- "delta_media_escala"
  obs <- do.call(rbind, lapply(setdiff(seq_len(ng), ref), function(g) do.call(rbind, lapply(fatores, function(f) {
    its <- itens[est[[ref]]$lambda[, f] != 0 | est[[g]]$lambda[, f] != 0]
    data.frame(grupo_focal = rot[g], fator = f, dif_observada = mean(rowSums(dat[[g]][, its, drop = FALSE]), na.rm = TRUE) - mean(rowSums(dat[[ref]][, its, drop = FALSE]), na.rm = TRUE))
  }))))
  esc <- merge(esc, obs, all.x = TRUE)
  esc$proporcao_vies <- esc$delta_media_escala / esc$dif_observada
  grandes <- tab[tab$dmacs >= .4, , drop = FALSE]
  avisos <- character()
  if (metrica_fixa) {
    esc$proporcao_vies <- NA_real_
    avisos <- c(avisos, "Modelo com média 0 e variância 1 do traço em todos os grupos (configural com std.lv): toda diferença entre os grupos vira diferença de parâmetro de item, e o ΔMean reproduz a diferença observada por construção. Para separar viés de diferença real no traço, rodar no modelo de invariância parcial (itens âncora iguais, médias latentes livres).")
  }
  if (nrow(grandes)) avisos <- c(avisos, paste0("dMACS ≥ ,40 (médio ou maior pela régua de Nye et al., 2019) em: ", paste(unique(grandes$item), collapse = ", "),
                                               ". Conferir o conteúdo desses itens e se a diferença entre os grupos muda depois de liberá-los."))
  if (any(abs(esc$proporcao_vies) > .5 & is.finite(esc$proporcao_vies)))
    avisos <- c(avisos, "Em alguma escala o viés de medida responde por mais da metade da diferença observada entre os grupos (proporcao_vies): a comparação de médias brutas fica comprometida.")
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  mx <- tab[which.max(tab$dmacs), ]
  apa <- if (en) paste0("Item-level noninvariance effect sizes (dMACS; Nye & Drasgow, 2011), reference group ", rot[ref], ", ranged from ", .pr_num(min(tab$dmacs), 2, FALSE, decimal),
                        " to ", .pr_num(max(tab$dmacs), 2, FALSE, decimal), " (largest: ", mx$item, ").") else
    paste0("O tamanho do efeito da não invariância por item (dMACS; Nye & Drasgow, 2011), com ", rot[ref], " como referência, variou de ", .pr_num(min(tab$dmacs), 2, FALSE, decimal),
           " a ", .pr_num(max(tab$dmacs), 2, FALSE, decimal), " (maior: ", mx$item, ").")
  .ee_resultado("dMACS: tamanho do efeito da não invariância", n = sum(vapply(dat, nrow, 0)), tabela = tab,
                decisoes = c("Forma fechada do dMACS para indicadores contínuos (integral sobre a distribuição normal do traço no grupo focal, com as médias e (co)variâncias latentes estimadas).",
                             "DP combinado pelas variâncias dos dois grupos (o pacote dmacs 0.1.0, arquivado no CRAN, combina os DPs linearmente; a diferença é pequena).",
                             "delta_media: focal menos referência, na unidade do item; somado no fator é o ΔMean de Nye e Drasgow (2011), comparado à diferença observada do escore somado.",
                             "A métrica latente precisa ser comparável: no configural com std.lv em todos os grupos o traço tem média 0 e variância 1 nos dois, e o dMACS mede tudo o que difere nos parâmetros dos itens."),
                avisos = avisos, apa = apa, extra = list(tabelas_extra = list("Viés na escala (ΔMean)" = esc), escala = esc, grupo_ref = rot[ref]))
}

# Teste de equivalência do ajuste (v1.7; Counsell, Cribbie & Flora, 2020, como o
# equiv_chi() que Petersen, 2024, 15.10.3.5, usa): em vez de H0 "o modelo se ajusta",
# H0 é "o desajuste é pelo menos o de um RMSEA populacional ruim" (,08 por padrão).
# ε0 = gl × RMSEA0² / grupos; λ = (N - grupos) × ε0; p = P(χ²(gl, λ) <= χ² observado).
# p < α rejeita o desajuste ruim: o modelo (ou a piora entre modelos aninhados, com o Δχ²
# e o Δgl) é equivalente a um de ajuste aceitável. Serve para a sequência de invariância:
# p não significativo no passo diz que a invariância não foi demonstrada.
equivalencia_ajuste <- function(qui2, gl, n, grupos = 1, rmsea0 = .08, alfa = .05, decimal = ".", idioma = "pt") {
  if (length(qui2) != length(gl)) stop("qui2 e gl precisam ter o mesmo comprimento", call. = FALSE)
  eps0 <- gl * rmsea0^2 / grupos
  lambda <- (n - grupos) * eps0
  p <- stats::pchisq(qui2, gl, ncp = lambda)
  tab <- data.frame(qui2 = qui2, gl = gl, Fml = qui2 / (n - grupos), epsilon0 = eps0, lambda = lambda, p = p, equivalente = p < alfa)
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  apa <- vapply(seq_along(qui2), function(i) if (en) paste0("Equivalence test of fit (Counsell et al., 2020; RMSEA0 = ", .pr_num(rmsea0, 2, FALSE, decimal), "): χ²(", gl[i], ") = ",
    .pr_num(qui2[i], 2, TRUE, decimal), sep, .pr_p(p[i], decimal), if (p[i] < alfa) ", so the misfit is smaller than that of an RMSEA0 model." else ", so equivalence to acceptable fit was not shown.") else
    paste0("Teste de equivalência do ajuste (Counsell et al., 2020; RMSEA0 = ", .pr_num(rmsea0, 2, FALSE, decimal), "): χ²(", gl[i], ") = ", .pr_num(qui2[i], 2, TRUE, decimal), sep, .pr_p(p[i], decimal),
           if (p[i] < alfa) ", o desajuste é menor que o de um modelo com esse RMSEA." else ", a equivalência com um ajuste aceitável não foi demonstrada."), character(1))
  .ee_resultado("Teste de equivalência do ajuste", n = n, tabela = tab,
                decisoes = c("H0 de desajuste: RMSEA populacional >= RMSEA0; p = P(χ² não central <= observado), como o equiv_chi() de Counsell et al. (2020) usado por Petersen (2024).",
                             "Com o Δχ² e o Δgl de dois modelos aninhados, o teste diz se a piora da restrição é equivalente a desprezível. Com WLSMV ou MLR, usar o χ² escalonado e saber que a distribuição não central é aproximação."),
                apa = apa)
}

# Invariância ao longo de uma variável contínua ou ordinal (v1.7): testes baseados no escore
# (Merkle & Zeileis, 2013; Wang, Merkle & Zeileis, 2014), via strucchange::sctest num
# modelo lavaan de um grupo estimado por ML. Ordena os casos pela variável (idade, por
# exemplo) e testa se as contribuições ao escore dos parâmetros escolhidos flutuam mais que
# o acaso, sem cortar a variável em grupos. Funcional: "DM" (máximo duplo) e "CvM" para
# contínua; "WDMo" e "maxLMo" para ordinal; "LMuo" para categórica (o que Petersen usa).
invariancia_continua <- function(fit, ordenar_por, parametros = c("cargas", "interceptos", "todos"), funcional = NULL,
                                 decimal = ".", idioma = "pt") {
  .pr_exige("lavaan", "strucchange")
  parametros <- if (is.character(parametros)) match.arg(parametros) else parametros
  if (lavaan::lavInspect(fit, "ngroups") > 1) stop("invariancia_continua() é para modelo de um grupo", call. = FALSE)
  if (!lavaan::lavInspect(fit, "options")$estimator %in% c("ML")) stop("os testes baseados no escore pedem estimação por ML (estimator = \"ML\")", call. = FALSE)
  n <- lavaan::lavInspect(fit, "nobs")
  if (length(ordenar_por) != n) stop("ordenar_por precisa ter um valor por caso usado no modelo (", n, ")", call. = FALSE)
  if (anyNA(ordenar_por)) stop("ordenar_por tem ausentes: ajustar o modelo só nos casos com o moderador observado", call. = FALSE)
  pt <- lavaan::parTable(fit); pt <- pt[pt$free > 0, ]
  idx <- if (is.numeric(parametros)) parametros else switch(parametros,
    cargas = pt$free[pt$op == "=~"], interceptos = pt$free[pt$op == "~1" & pt$lhs %in% lavaan::lavNames(fit, "ov")], todos = pt$free)
  if (!length(idx)) stop("nenhum parâmetro livre do tipo pedido (interceptos pedem meanstructure = TRUE)", call. = FALSE)
  tipo <- if (is.ordered(ordenar_por)) "ordinal" else if (is.factor(ordenar_por) || is.character(ordenar_por)) "categorica" else "continua"
  if (tipo == "categorica") ordenar_por <- factor(ordenar_por)
  if (is.null(funcional)) funcional <- switch(tipo, continua = "DM", ordinal = "WDMo", categorica = "LMuo")
  ft <- funcional
  if (identical(funcional, "WDMo")) ft <- strucchange::ordwmax(ordenar_por)
  if (identical(funcional, "maxLMo")) ft <- strucchange::ordL2BB(ordenar_por, nproc = length(idx))
  if (identical(funcional, "LMuo")) ft <- strucchange::catL2BB(ordenar_por)
  st <- suppressWarnings(strucchange::sctest(fit, order.by = ordenar_por, parm = idx, functional = ft))
  nomes <- paste0(pt$lhs, pt$op, pt$rhs)[match(idx, pt$free)]
  tab <- data.frame(funcional = funcional, estatistica = unname(st$statistic), p = st$p.value, parametros = length(idx), tipo_moderador = tipo)
  en <- identical(idioma, "en"); sep <- .ee_sep(decimal)
  rot <- if (is.character(parametros)) c(cargas = "cargas", interceptos = "interceptos", todos = "todos os parâmetros")[[parametros]] else "parâmetros escolhidos"
  apa <- if (en) paste0("A score-based test of measurement invariance (Merkle & Zeileis, 2013; ", funcional, " statistic) along the ", tipo, " moderator for the ",
                        length(idx), " ", parametros, " parameters gave ", .ee_estat(funcional, st$statistic, NULL, decimal), sep, .pr_p(st$p.value, decimal), ".") else
    paste0("O teste de invariância baseado no escore (Merkle & Zeileis, 2013; estatística ", funcional, ") ao longo do moderador ", tipo, ", para ", length(idx), " ", rot, ", deu ",
           .ee_estat(funcional, st$statistic, NULL, decimal), sep, .pr_p(st$p.value, decimal), ".")
  .ee_resultado("Invariância ao longo de variável contínua ou ordinal (testes baseados no escore)", n = n, tabela = tab,
                decisoes = c(paste0("strucchange::sctest no modelo lavaan (ML), parâmetros testados: ", paste(nomes, collapse = ", "), "."),
                             "Funcional pelo tipo do moderador: DM para contínuo, WDMo para ordinal, LMuo para categórico (Merkle, Fan & Zeileis, 2014); outro pode ser pedido em funcional.",
                             "p < ,05: os parâmetros mudam ao longo do moderador; o teste não diz em que ponto, e o gráfico do processo de flutuação (plot(gefp(fit, order.by = ...))) mostra onde."),
                avisos = c(if (st$p.value < .05) "Parâmetros não invariantes ao longo do moderador: comparar médias ou correlações ao longo dele pede modelagem que acomode a mudança (MNLFA ou grupos com invariância parcial)." else
                  "p >= ,05 não prova invariância ao longo do moderador; o poder cai com n pequeno e com muitos parâmetros testados de uma vez."),
                apa = apa, extra = list(teste = st, parametros_testados = nomes))
}
