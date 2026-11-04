# est_comum.R: base comum das análises estatísticas da skill psicometria-r.
# Carregar sempre por scripts/carregar.R, que traz psicometria.R antes (usa .pr_num,
# .pr_p, .pr_exige e .pr_e_ordinal de lá).
#
# Contrato de todas as funções de análise (est_*.R):
#  - entrada: data.frame e nomes de variáveis como texto; nada é lido de arquivo aqui;
#  - casos incompletos nas variáveis da análise saem por deleção listwise, e o número
#    que saiu vai para o resultado (n_excluidos);
#  - saída: objeto de classe "pr_resultado" (ver .ee_resultado), com tabelas de
#    agregados, decisões tomadas e a frase APA; nenhuma linha de participante é
#    impressa (índices de casos atípicos ficam no objeto, fora do print);
#  - decisões de rota (Student ou Welch, correção de esfericidade, pós-hoc) seguem o
#    livro do usuário por padrão (rota = "classica") e ficam escritas em $decisoes;
#    a alternativa recomendada pela literatura sai sempre ao lado quando for barata.

.ee_versao <- "2.2.0"

# ---------------------------------------------------------------- objeto de resultado
.ee_resultado <- function(analise, n, n_excluidos = 0, tabela = NULL, efeito = NULL,
                          pressupostos = NULL, posthoc = NULL, decisoes = character(),
                          avisos = character(), apa = character(), extra = list()) {
  r <- c(list(analise = analise, n = n, n_excluidos = n_excluidos, tabela = tabela,
              efeito = efeito, pressupostos = pressupostos, posthoc = posthoc,
              decisoes = decisoes, avisos = avisos, apa = apa), extra)
  class(r) <- c("pr_resultado", "list")
  r
}

.ee_print_df <- function(d, titulo) {
  if (is.null(d) || !NROW(d)) return(invisible())
  cat("\n", titulo, "\n", sep = "")
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(x) signif(x, 4))
  print(as.data.frame(d), row.names = FALSE)
}

print.pr_resultado <- function(x, ...) {
  cat("== ", x$analise, " (n = ", x$n,
      if (isTRUE(x$n_excluidos > 0)) paste0("; ", x$n_excluidos, " excluídos por ausência") else "",
      ")\n", sep = "")
  if (is.data.frame(x$descritivas)) {
    dsc <- .ee_suprime_pequenos(as.data.frame(x$descritivas))
    .ee_print_df(dsc, "Descritivas:")
    if (!is.null(attr(dsc, "suprimidos")))
      cat("(grupo com n < ", .ee_n_min_print, " sem estatísticas no print, por privacidade; valores no objeto)\n", sep = "")
  }
  .ee_print_df(x$pressupostos, "Pressupostos:")
  .ee_print_df(x$tabela, "Resultado:")
  .ee_print_df(x$efeito, "Tamanho de efeito:")
  .ee_print_df(x$posthoc, if (is.character(x$titulo_posthoc)) x$titulo_posthoc else "Comparações:")
  # Tabelas extras com coluna n (médias por grupo, codificação de grupos) passam pela
  # mesma supressão de grupo pequeno das descritivas.
  if (is.list(x$tabelas_extra))
    for (nm in names(x$tabelas_extra)) {
      te <- x$tabelas_extra[[nm]]
      if (is.data.frame(te)) te <- .ee_suprime_pequenos(as.data.frame(te))
      .ee_print_df(te, paste0(nm, ":"))
      if (!is.null(attr(te, "suprimidos")))
        cat("(grupo com n < ", .ee_n_min_print, " sem estatísticas no print, por privacidade; valores no objeto)\n", sep = "")
    }
  if (length(x$decisoes)) cat("\nDecisões:\n", paste0("  - ", x$decisoes, collapse = "\n"), "\n", sep = "")
  if (length(x$avisos)) cat("\nAvisos:\n", paste0("  ! ", x$avisos, collapse = "\n"), "\n", sep = "")
  if (length(x$apa)) cat("\nAPA: ", paste(x$apa, collapse = "\n     "), "\n", sep = "")
  invisible(x)
}

relatar <- function(res) {
  stopifnot(inherits(res, "pr_resultado"))
  res$apa
}

# ---------------------------------------------------------------- não significativo
# p >= alfa não mostra ausência de efeito: os dados não destoam da hipótese nula, como
# não destoam de outros valores do intervalo, e o teste pode ter tido pouco poder
# (Greenland et al., 2016, equívocos 4, 6 e 8; Lakens, 2017). As funções de diferença e
# de associação chamam este aviso quando o teste principal não é significativo.
#  ic: trecho com a estimativa e o IC (fica de fora em teste unilateral, cujo IC é
#      aberto);
#  tost: tipo de teste_equivalencia() que cobre a estatística (NULL quando a skill não
#      tem teste de equivalência para ela);
#  rotulos: com vários termos ou pares, um rótulo por p (com o p e, se houver, o IC);
#      o aviso lista só os não significativos.
# Só em português, como os demais avisos.
.ee_aviso_nao_sig <- function(p, alfa = .05, ic = NULL, o_que = "diferença", tost = NULL, decimal = ".",
                              rotulos = NULL) {
  ns <- !is.na(p) & p >= alfa
  if (!any(ns)) return(character())
  ini <- if (is.null(rotulos)) {
    paste0(.pr_p(p[1], decimal), " só diz que os dados não destoam da ausência de ", o_que,
           ", e eles também não destoam de outros valores", if (!is.null(ic)) paste0(" (", ic, ")") else "")
  } else {
    paste0("em ", paste(rotulos[ns], collapse = "; "), ", os dados só não destoam da ausência de ", o_que,
           ", e também não destoam de outros valores")
  }
  paste0("Não significativo não quer dizer igual: ", ini, " (Greenland et al., 2016). ",
         if (!is.null(tost)) paste0("Para afirmar que o efeito é pequeno demais para importar, teste_equivalencia(tipo = \"", tost,
                                    "\") com limites fixados antes da análise pelo menor efeito de interesse (Lakens, 2017).") else
           paste0("A skill não tem teste de equivalência para esta estatística: relatar o IC do efeito e não escrever que não há ",
                  o_que, " (Lakens, 2017)."))
}

# ---------------------------------------------------------------- dados
.ee_dados <- function(dados, vars, numericas = character()) {
  if (!is.data.frame(dados)) stop("dados precisa ser um data.frame")
  falta <- setdiff(vars, names(dados))
  if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "))
  d <- as.data.frame(dados)[, vars, drop = FALSE]
  # Célula de texto vazia (ou só com espaços) é ausente, como no JASP; sem isso, o
  # read.csv cria uma categoria "" que entra na análise como grupo.
  for (v in vars) {
    x <- d[[v]]
    if (is.character(x)) { x[!is.na(x) & !nzchar(trimws(x))] <- NA; d[[v]] <- x }
    if (is.factor(x) && any(!nzchar(trimws(levels(x))))) {
      # Mantém a ordem original dos níveis (convenção de sinal, grupo de referência).
      lv <- levels(x); lv <- lv[nzchar(trimws(lv))]
      x <- as.character(x); x[!is.na(x) & !nzchar(trimws(x))] <- NA
      d[[v]] <- droplevels(factor(x, levels = lv))
    }
  }
  for (v in numericas) {
    if (!is.numeric(d[[v]])) stop("a variável ", v, " precisa ser numérica")
  }
  completos <- stats::complete.cases(d)
  list(d = d[completos, , drop = FALSE], n_excluidos = sum(!completos))
}

.ee_fator <- function(x, nome = "grupo") {
  f <- droplevels(as.factor(x))
  if (nlevels(f) < 2) stop("o fator ", nome, " precisa de pelo menos 2 níveis com dados")
  f
}

# ---------------------------------------------------------------- formatação APA
# Convenção: estatística de teste com 2 casas; p com 3 casas (p < .001); quantidades
# que não passam de 1 em módulo (p, r, R², η², ω², V, φ, AUC, poder) sem zero à
# esquerda; d, g, OR, B, IRR, f² e β padronizado com zero (o β pode passar de
# 1). Proporções de variância (η², ω², R², ΔR², f², ε²) com 3 casas, os demais
# efeitos com 2, exceto o módulo de associação, que usa 3 como os parágrafos-modelo do
# livro. decimal = "." por padrão (neutro entre idiomas); manuscrito em português usa
# decimal = ",", como o livro, e aí os itens da frase são separados por ponto
# e vírgula.
.ee_sep <- function(decimal) if (decimal == ",") "; " else ", "

.ee_gl <- function(gl, decimal = ".") {
  vapply(gl, function(g) {
    if (is.na(g)) return("NA")
    if (abs(g - round(g)) < 1e-8) return(format(round(g), big.mark = ""))
    .pr_num(g, 2, TRUE, decimal)
  }, character(1))
}

.ee_ic <- function(lo, hi, casas = 2, zero = TRUE, decimal = ".", nivel = .95, idioma = "pt") {
  rot <- if (idioma == "en") paste0(round(nivel * 100), "% CI") else paste0("IC ", round(nivel * 100), "%")
  paste0(rot, " [", .pr_num(lo, casas, zero, decimal), .ee_sep(decimal),
         .pr_num(hi, casas, zero, decimal), "]")
}

.ee_estat <- function(simbolo, valor, gl = NULL, decimal = ".", casas = 2, extra_gl = NULL) {
  g <- if (is.null(gl)) "" else paste0("(", paste(c(.ee_gl(gl, decimal), extra_gl),
                                                  collapse = .ee_sep(decimal)), ")")
  paste0(simbolo, g, " = ", .pr_num(valor, casas, TRUE, decimal))
}

.ee_efeito_txt <- function(simbolo, valor, lo = NA, hi = NA, limitado = TRUE, casas = 2,
                           decimal = ".", idioma = "pt", nivel = .95) {
  s <- paste0(simbolo, " = ", .pr_num(valor, casas, !limitado, decimal))
  if (!is.na(lo) && !is.na(hi)) s <- paste0(s, .ee_sep(decimal), .ee_ic(lo, hi, casas, !limitado, decimal, nivel, idioma))
  s
}

.ee_junta <- function(..., decimal = ".") paste(c(...), collapse = .ee_sep(decimal))

# ---------------------------------------------------------------- descritivas
# Assimetria e curtose no estimador do SPSS e do JASP (G1 e G2, curtose em excesso).
.ee_assimetria <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x); s2 <- sum((x - m)^2) / n
  g1 <- sum((x - m)^3) / n / s2^1.5
  g1 * sqrt(n * (n - 1)) / (n - 2)
}

.ee_curtose <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 4) return(NA_real_)
  m <- mean(x); s2 <- sum((x - m)^2) / n
  g2 <- sum((x - m)^4) / n / s2^2 - 3
  ((n + 1) * g2 + 6) * (n - 1) / ((n - 2) * (n - 3))
}

.ee_ep_assimetria <- function(n) sqrt(6 * n * (n - 1) / ((n - 2) * (n + 1) * (n + 3)))
.ee_ep_curtose <- function(n) 2 * .ee_ep_assimetria(n) * sqrt((n^2 - 1) / ((n - 3) * (n + 5)))

descritivas <- function(dados, vars, grupo = NULL) {
  todas <- c(vars, grupo)
  falta <- setdiff(todas, names(dados))
  if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "))
  um <- function(x) {
    x_ok <- x[!is.na(x)]; n <- length(x_ok)
    sw <- if (n >= 3 && n <= 5000 && stats::sd(x_ok) > 0) stats::shapiro.test(x_ok)$p.value else NA_real_
    ep <- if (n > 1) stats::sd(x_ok) / sqrt(n) else NA_real_
    tq <- if (n > 1) stats::qt(.975, n - 1) else NA_real_
    data.frame(n = n, ausentes = sum(is.na(x)), media = mean(x_ok), dp = stats::sd(x_ok), ep = ep,
               ic95_inf = mean(x_ok) - tq * ep, ic95_sup = mean(x_ok) + tq * ep,
               mediana = stats::median(x_ok), iiq = stats::IQR(x_ok),
               minimo = if (n) min(x_ok) else NA, maximo = if (n) max(x_ok) else NA,
               assimetria = .ee_assimetria(x_ok),
               ep_assimetria = if (n >= 3) .ee_ep_assimetria(n) else NA_real_,
               curtose = .ee_curtose(x_ok),
               ep_curtose = if (n >= 4) .ee_ep_curtose(n) else NA_real_,
               shapiro_p = sw)
  }
  out <- list()
  for (v in vars) {
    if (!is.numeric(dados[[v]])) stop("a variável ", v, " precisa ser numérica")
    if (is.null(grupo)) {
      out[[length(out) + 1]] <- cbind(variavel = v, um(dados[[v]]))
    } else {
      gch <- as.character(dados[[grupo]]); gch[!is.na(gch) & !nzchar(trimws(gch))] <- NA
      g <- if (is.factor(dados[[grupo]])) droplevels(factor(gch, levels = levels(dados[[grupo]]))) else droplevels(as.factor(gch))
      for (lv in levels(g)) {
        out[[length(out) + 1]] <- cbind(variavel = v, grupo = lv, um(dados[[v]][!is.na(g) & g == lv]))
      }
    }
  }
  res <- do.call(rbind, out)
  class(res) <- c("pr_descritivas", "data.frame")
  res
}

# Privacidade: no print, linha com n < 5 (grupo pequeno) mostra só n e ausentes; os
# valores continuam no objeto. Com n = 1, média, mínimo e máximo seriam o escore da
# pessoa.
.ee_n_min_print <- 5
.ee_suprime_pequenos <- function(d) {
  if (!is.data.frame(d) || !"n" %in% names(d)) return(d)
  peq <- !is.na(d$n) & d$n < .ee_n_min_print
  if (!any(peq)) return(d)
  manter <- c("variavel", "grupo", "n", "ausentes")
  for (cn in setdiff(names(d), manter)) if (is.numeric(d[[cn]])) d[[cn]][peq] <- NA
  attr(d, "suprimidos") <- sum(peq)
  d
}

print.pr_descritivas <- function(x, ...) {
  d <- .ee_suprime_pequenos(as.data.frame(unclass(x), stringsAsFactors = FALSE))
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(v) signif(v, 4))
  print.data.frame(d, row.names = FALSE)
  if (!is.null(attr(d, "suprimidos")))
    cat("(", attr(d, "suprimidos"), " linha(s) com n < ", .ee_n_min_print,
        " sem estatísticas no print, por privacidade; os valores ficam no objeto)\n", sep = "")
  invisible(x)
}

# Leitura de normalidade do livro (Vol. 1, cap. 5): avaliação contínua, nunca um
# liga-desliga. |assimetria| ou |curtose| abaixo de 1 = afastamento leve; acima de 2 =
# preocupa (heurística, não teste). O Shapiro-Wilk entra como um indício entre vários e
# perde sentido com n grande. Devolve texto por variável. Linha com n < 5 não ganha
# leitura (com tão poucos casos a forma não se avalia, e o texto exporia o grupo).
ler_normalidade <- function(desc, decimal = ".") {
  desc <- as.data.frame(unclass(desc), stringsAsFactors = FALSE)
  apply(desc, 1, function(l) {
    n <- as.numeric(l[["n"]])
    rotulo <- paste0(l[["variavel"]], if ("grupo" %in% names(l)) paste0(" [", l[["grupo"]], "]") else "")
    if (is.finite(n) && n < .ee_n_min_print) return(paste0(rotulo, ": n < ", .ee_n_min_print, ", sem leitura de normalidade"))
    a <- abs(as.numeric(l[["assimetria"]])); k <- abs(as.numeric(l[["curtose"]]))
    p <- suppressWarnings(as.numeric(l[["shapiro_p"]]))
    m <- max(a, k, na.rm = TRUE)
    faixa <- if (!is.finite(m)) "indeterminado" else if (m < 1) "afastamento leve" else if (m <= 2) "afastamento moderado" else "afastamento que preocupa"
    nota_n <- if (n >= 300) " (com n grande o Shapiro-Wilk acusa desvios irrelevantes: pesa o tamanho do afastamento)" else
      if (n < 30) " (com n pequeno o Shapiro-Wilk tem pouco poder: olhar o Q-Q plot)" else ""
    sp <- if (decimal == ",") "; " else ", "
    paste0(rotulo, ": ", faixa, "; assimetria ", .pr_num(as.numeric(l[["assimetria"]]), 2, TRUE, decimal),
           sp, "curtose ", .pr_num(as.numeric(l[["curtose"]]), 2, TRUE, decimal),
           if (!is.na(p)) paste0(sp, "Shapiro-Wilk ", .pr_p(p, decimal)) else "",
           nota_n)
  })
}

escores_z <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)

# Casos atípicos univariados. metodo "z": |z| acima do limite (padrão 3,29, p < .001
# bilateral); "mad": distância à mediana em MADs (Leys et al., 2013, limite padrão 2,5).
# O print mostra só contagens; os índices ficam em $indices, para uso local.
atipicos <- function(dados, vars, metodo = c("z", "mad"), limite = NULL) {
  metodo <- match.arg(metodo)
  if (is.null(limite)) limite <- if (metodo == "z") 3.29 else 2.5
  ind <- list(); tab <- list()
  for (v in vars) {
    x <- dados[[v]]
    s <- if (metodo == "z") escores_z(x) else (x - stats::median(x, na.rm = TRUE)) / stats::mad(x, na.rm = TRUE)
    i <- which(abs(s) > limite)
    ind[[v]] <- i
    tab[[v]] <- data.frame(variavel = v, metodo = metodo, limite = limite, n_atipicos = length(i),
                           prop = length(i) / sum(!is.na(x)))
  }
  .ee_resultado("Casos atípicos univariados", n = nrow(dados), tabela = do.call(rbind, tab),
                decisoes = "Caso atípico não sai do banco automaticamente: decidir por critério declarado antes da análise e relatar quantos saíram e por quê.",
                extra = list(indices = ind))
}

# Índice de mudança confiável e significância clínica (Jacobson & Truax, 1991):
# RCI = (pós - pré) / S_dif, S_dif = sqrt(2 * EP²), EP = DP_pre * sqrt(1 - confiabilidade).
# |RCI| > 1,96 = mudança confiável a 5%.
# Critérios de corte clínico de normalidade para recuperação:
#   Critério A: 2 DP da média clínica na direção funcional;
#   Critério B: 2 DP da média funcional/normativa;
#   Critério C: Ponto de intersecção ponderado pelas variâncias (M_clin*DP_norm + M_norm*DP_clin)/(DP_clin + DP_norm).
mudanca_confiavel <- function(pre, pos, dp_ref, confiabilidade, melhora = c("diminuir", "aumentar"),
                              media_clinica = NULL, dp_clinica = NULL,
                              media_norma = NULL, dp_norma = NULL,
                              criterio_recuperacao = c("nenhum", "A", "B", "C")) {
  melhora <- match.arg(melhora)
  criterio_recuperacao <- match.arg(criterio_recuperacao)
  ep <- dp_ref * sqrt(1 - confiabilidade); sdif <- sqrt(2 * ep^2)
  rci <- (pos - pre) / sdif
  sinal <- if (melhora == "diminuir") -1 else 1
  
  corte_clinico <- NULL
  decisao_corte <- NULL
  if (criterio_recuperacao != "nenhum") {
    if (criterio_recuperacao == "A") {
      if (is.null(media_clinica) || is.null(dp_clinica)) {
        stop("Critério A exige 'media_clinica' e 'dp_clinica' (2 DP da média clínica na direção funcional).")
      }
      corte_clinico <- if (melhora == "diminuir") media_clinica - 2 * dp_clinica else media_clinica + 2 * dp_clinica
      decisao_corte <- paste0("Corte Clínico A = ", round(corte_clinico, 2), " (M_clinica ", if (melhora == "diminuir") "- 2*DP" else "+ 2*DP", ").")
    } else if (criterio_recuperacao == "B") {
      if (is.null(media_norma) || is.null(dp_norma)) {
        stop("Critério B exige 'media_norma' e 'dp_norma' (2 DP da média normativa/funcional).")
      }
      corte_clinico <- if (melhora == "diminuir") media_norma + 2 * dp_norma else media_norma - 2 * dp_norma
      decisao_corte <- paste0("Corte Clínico B = ", round(corte_clinico, 2), " (M_norma ", if (melhora == "diminuir") "+ 2*DP" else "- 2*DP", ").")
    } else if (criterio_recuperacao == "C") {
      if (is.null(media_clinica) || is.null(dp_clinica) || is.null(media_norma) || is.null(dp_norma)) {
        stop("Critério C exige 'media_clinica', 'dp_clinica', 'media_norma' e 'dp_norma' (ponto ponderado entre as distribuições).")
      }
      corte_clinico <- (dp_clinica * media_norma + dp_norma * media_clinica) / (dp_clinica + dp_norma)
      decisao_corte <- paste0("Corte Clínico C = ", round(corte_clinico, 2), " (ponto de intersecção ponderado entre clínica e funcional).")
    }
  }
  
  if (is.null(corte_clinico)) {
    cls <- ifelse(is.na(rci), NA, ifelse(sinal * rci > 1.96, "melhora confiável",
                  ifelse(sinal * rci < -1.96, "piora confiável", "sem mudança confiável")))
  } else {
    passou_corte <- if (melhora == "diminuir") (pos < corte_clinico) else (pos > corte_clinico)
    cls <- ifelse(is.na(rci), NA,
           ifelse(sinal * rci < -1.96, "deteriorado",
           ifelse(sinal * rci <= 1.96, "inalterado",
           ifelse(passou_corte, "clinicamente recuperado", "melhora confiável sem recuperação"))))
  }
  
  decisoes <- c(
    paste0("S_dif = ", round(sdif, 3), "; mudança mínima confiável = ", round(1.96 * sdif, 2), " pontos (|RCI| > 1.96, p < .05)."),
    if (!is.null(decisao_corte)) decisao_corte else character(),
    paste0("O RCI é convenção de ensaio clínico em psicologia (Jacobson & Truax, 1991) e supõe o mesmo erro em toda a escala. ",
           "Alternativa psicométrica: erro padrão condicional ao nível do traço pela TRI (tri_ajustar(), curva de informação), ",
           "que muda ao longo da escala; ver references/fontes-e-decisoes.md.")
  )
  
  .ee_resultado("Índice de mudança confiável e significância clínica (Jacobson & Truax, 1991)",
                n = sum(!is.na(rci)),
                tabela = as.data.frame(table(classificacao = cls)),
                decisoes = decisoes,
                extra = list(rci = rci, classificacao = cls, corte_clinico = corte_clinico, criterio = criterio_recuperacao))
}

# ---------------------------------------------------------------- pacotes
.ee_pacotes_estatistica <- c("car", "afex", "emmeans", "effectsize", "lmerTest", "lme4",
                             "pwrss", "pwr", "pROC", "geepack", "metafor", "heplots",
                             "onewaytests", "rstatix", "ppcor", "lmtest", "sandwich", "boot",
                             "MASS", "psych", "GPArotation", "lavaan", "EGAnet", "bootnet",
                             "qgraph", "networktools", "mgm", "IsingFit",
                             "mirt", "semTools", "poLCA", "mclust", "survival",
                             "rankFD", "nparcomp", "nparLD", "MANOVA.RM", "netmeta", "meta", "bnlearn", "thurstonianIRT", "norm",
                             "nnet", "ordinal", "EFA.MRFA", "mice", "graphicalVAR", "mlVAR", "mokken", "mirtCAT", "strucchange")
# Só os evals usam (bancos de exemplo e conferências com pacotes de referência). Os
# evals da v1.4 pulam, com aviso, a conferência ao vivo quando o pacote falta.
.ee_pacotes_evals <- c("mlmRev", "metadat", "PMCMRplus", "TOSTER", "EFAtools", "irr", "BifactorIndicesCalculator",
                       "difR", "eRm", "tidyLPA", "meta", "mada", "naniar", "epitools",
                       "VGAM", "irrCAC", "pinsearch", "ShinyItemAnalysis", "rms", "DescTools")

checar_pacotes <- function(pacotes = .ee_pacotes_estatistica) {
  inst <- vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)
  data.frame(pacote = pacotes, instalado = inst,
             versao = vapply(pacotes, function(p) if (requireNamespace(p, quietly = TRUE))
               as.character(utils::packageVersion(p)) else NA_character_, character(1)),
             row.names = NULL)
}

# ---------------------------------------------------------------- v1.5: dados ausentes
# ausentes(): diagnóstico agregado dos dados ausentes, a aula "Dados faltantes: guia das
# melhores práticas". Sai a porcentagem por variável, o número de casos completos, os
# padrões de ausência (só os com n >= 5, por privacidade) e o teste MCAR de Little (1988),
# com média e covariância por EM (norm). As aulas imputam por EM (uma vez) quando o
# Little não rejeita; a skill registra que a imputação única por EM subestima os EP e
# recomenda FIML (ML/MLR no lavaan) ou imputação múltipla para a inferência. Little não
# significativo não prova MCAR (MAR e MCAR não se distinguem pelos dados).
ausentes <- function(dados, vars, decimal = ".", idioma = "pt") {
  falta <- setdiff(vars, names(dados)); if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "))
  d <- as.data.frame(dados[vars]); n <- nrow(d)
  M <- is.na(d)
  porv <- data.frame(variavel = vars, n_ausentes = colSums(M), pct_ausentes = 100 * colMeans(M), row.names = NULL)
  pad <- apply(M, 1, function(r) paste(as.integer(r), collapse = ""))
  tp <- sort(table(pad), decreasing = TRUE)
  padroes <- data.frame(padrao = names(tp), n = as.integer(tp), stringsAsFactors = FALSE)
  padroes$variaveis_ausentes <- vapply(padroes$padrao, function(p) { b <- strsplit(p, "")[[1]] == "1"; if (any(b)) paste(vars[b], collapse = ", ") else "(completo)" }, "")
  n_pad <- nrow(padroes); padroes_print <- padroes[padroes$n >= 5, ]
  lt <- NULL
  num <- vapply(d, is.numeric, TRUE)
  if (all(num) && any(M) && requireNamespace("norm", quietly = TRUE)) {
    x <- as.matrix(d); p <- ncol(x)
    s <- norm::prelim.norm(x)
    th <- norm::em.norm(s, showits = FALSE)
    pr <- norm::getparam.norm(s, th)
    mu <- pr$mu; S <- pr$sigma
    d2 <- 0; gl <- 0
    for (pp in unique(pad)) {
      i <- which(pad == pp); o <- which(strsplit(pp, "")[[1]] == "0")
      if (!length(o)) next
      yb <- colMeans(x[i, o, drop = FALSE]) - mu[o]
      d2 <- d2 + length(i) * drop(t(yb) %*% solve(S[o, o, drop = FALSE]) %*% yb)
      gl <- gl + length(o)
    }
    gl <- gl - p
    lt <- data.frame(teste = "Little (1988)", qui2 = d2, gl = gl, p = stats::pchisq(d2, gl, lower.tail = FALSE), padroes = n_pad)
  }
  completos <- sum(stats::complete.cases(d))
  pct_tot <- 100 * mean(M)
  en <- identical(idioma, "en")
  apa <- if (en) paste0(.pr_num(pct_tot, 1, TRUE, decimal), "% of the values were missing (", completos, " of ", n, " complete cases)",
                        if (!is.null(lt)) paste0("; Little's MCAR test: χ²(", lt$gl, ") = ", .pr_num(lt$qui2, 2, TRUE, decimal), ", ", .pr_p(lt$p, decimal)) else "", ".") else
    paste0(.pr_num(pct_tot, 1, TRUE, decimal), "% dos valores estavam ausentes (", completos, " de ", n, " casos completos)",
           if (!is.null(lt)) paste0("; teste MCAR de Little: χ²(", lt$gl, ") = ", .pr_num(lt$qui2, 2, TRUE, decimal), ", ", .pr_p(lt$p, decimal)) else "", ".")
  .ee_resultado("Dados ausentes", n = n, tabela = porv,
                decisoes = c("Little com média e covariância por EM (pacote norm), em variáveis numéricas.",
                             "Little não significativo não prova MCAR; significativo afasta MCAR, sem dizer se é MAR ou MNAR.",
                             "Rota para a análise: FIML (ML ou MLR) ou imputação múltipla. Imputação única por EM, ensinada nas aulas no SPSS, serve para estimar médias e covariâncias, não para testes: subestima os EP.",
                             "Padrões com n < 5 ficam fora do print (privacidade); o número total de padrões está no teste."),
                avisos = c(if (pct_tot > 20) "Mais de 20% de valores ausentes: a escolha do método pesa no resultado; fazer análise de sensibilidade." else character(),
                           if (!all(num)) "Variáveis não numéricas: o teste de Little não foi calculado." else character()),
                apa = apa, extra = list(tabelas_extra = c(if (!is.null(lt)) list(`Teste MCAR de Little` = lt), list(`Padrões de ausência (n >= 5)` = padroes_print)),
                                        little = lt, padroes = padroes))
}

# ---------------------------------------------------------------- v1.7: imputação múltipla
# imputacao_multipla(): imputação múltipla por equações encadeadas (mice; van Buuren &
# Groothuis-Oudshoorn, 2011), a rota que a aula de dados faltantes apresenta como avanço
# sobre a imputação única por EM e que Enders e Baraldi (2018, Wiley, cap. 6) recomendam
# junto com a máxima verossimilhança. Ajusta o modelo (linear ou GLM) em cada banco
# imputado e combina pelas regras de Rubin: estimativa média, variância total
# T = U + (1 + 1/m) B, gl de Barnard e Rubin (1999), fração de informação ausente (FMI).
# Sai ao lado o mesmo modelo em casos completos, para mostrar o quanto a deleção mudava.
# Os bancos imputados ficam em $imputacoes (objeto mids: tem linhas de pessoa; nunca
# imprimir nem levar à conversa).
.ee_rubin <- function(Q, U, n, k, nivel = .95) {
  m <- nrow(Q); qbar <- colMeans(Q); ubar <- colMeans(U); b <- apply(Q, 2, stats::var)
  t <- ubar + (1 + 1 / m) * b
  riv <- (1 + 1 / m) * b / ubar; lambda <- (1 + 1 / m) * b / t
  dfcom <- n - k
  dfold <- (m - 1) / lambda^2
  dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
  gl <- dfold * dfobs / (dfold + dfobs)
  gl[!is.finite(gl)] <- dfcom
  fmi <- (riv + 2 / (gl + 3)) / (riv + 1)
  ep <- sqrt(t); tt <- qbar / ep; q <- stats::qt(1 - (1 - nivel) / 2, gl)
  data.frame(termo = names(qbar), B = qbar, ep = ep, t = tt, gl = gl, p = 2 * stats::pt(-abs(tt), gl),
             ic_inf = qbar - q * ep, ic_sup = qbar + q * ep, riv = riv, lambda = lambda, fmi = fmi, row.names = NULL)
}

imputacao_multipla <- function(dados, vars, formula = NULL, familia = c("normal", "binomial", "poisson"), m = 20,
                               metodo = NULL, maxit = 10, semente = 2026, nivel = .95, decimal = ".", idioma = "pt") {
  familia <- match.arg(familia)
  .pr_exige("mice")
  falta <- setdiff(vars, names(dados)); if (length(falta)) stop("variáveis ausentes no banco: ", paste(falta, collapse = ", "), call. = FALSE)
  d <- as.data.frame(dados)[vars]
  for (v in vars) if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
  if (!is.null(formula)) {
    formula <- stats::as.formula(formula)
    fora <- setdiff(all.vars(formula), vars)
    if (length(fora)) stop("a fórmula usa variável fora de vars (entre no modelo de imputação): ", paste(fora, collapse = ", "), call. = FALSE)
  }
  n <- nrow(d); inc <- mean(!stats::complete.cases(d)); pct_var <- 100 * colMeans(is.na(d))
  if (inc == 0) stop("não há valores ausentes nessas variáveis: imputação desnecessária", call. = FALSE)
  avisos <- character(); decisoes <- character()
  if (m < 100 * inc) avisos <- c(avisos, paste0("m = ", m, " é menor que a % de casos incompletos (", round(100 * inc, 1),
                                                "%): White et al. (2011) sugerem m pelo menos igual a essa porcentagem para reduzir o erro de Monte Carlo."))
  args <- list(data = d, m = m, maxit = maxit, seed = semente, printFlag = FALSE)
  if (!is.null(metodo)) args$method <- metodo
  imp <- suppressWarnings(do.call(mice::mice, args))
  if (length(imp$loggedEvents) && nrow(imp$loggedEvents)) avisos <- c(avisos, paste0("O mice registrou ", nrow(imp$loggedEvents),
    " evento(s) (variável constante ou colinear tirada do modelo de imputação): conferir $imputacoes$loggedEvents."))
  metodos <- imp$method[nzchar(imp$method)]
  aus <- data.frame(variavel = vars, pct_ausentes = unname(pct_var), metodo = unname(ifelse(nzchar(imp$method[vars]), imp$method[vars], "(completa)")))
  decisoes <- c(decisoes, paste0("Imputação múltipla por equações encadeadas (mice), m = ", m, ", ", maxit, " iterações, semente ", semente,
                                 "; métodos por variável na tabela extra (pmm = predictive mean matching, padrão para numéricas)."),
                "Todas as variáveis de vars entram no modelo de imputação; auxiliares correlacionadas com a ausência ou com as variáveis tornam o MAR mais plausível (Enders & Baraldi, 2018).")
  res <- NULL; cc <- NULL; tab <- aus; apa <- character(); efeito <- NULL
  if (!is.null(formula)) {
    ajusta <- function(dd) if (familia == "normal") stats::lm(formula, data = dd) else stats::glm(formula, data = dd, family = familia)
    fits <- lapply(seq_len(m), function(i) ajusta(mice::complete(imp, i)))
    Q <- t(vapply(fits, stats::coef, numeric(length(stats::coef(fits[[1]])))))
    U <- t(vapply(fits, function(f) diag(stats::vcov(f)), numeric(ncol(Q))))
    k <- ncol(Q)
    res <- .ee_rubin(Q, U, n, k, nivel)
    if (familia != "normal") { res$exp_B <- exp(res$B); res$exp_ic_inf <- exp(res$ic_inf); res$exp_ic_sup <- exp(res$ic_sup) }
    f_cc <- ajusta(stats::na.omit(d[all.vars(formula)]))
    s_cc <- summary(f_cc)$coefficients
    cc <- data.frame(termo = rownames(s_cc), B_completos = s_cc[, 1], ep_completos = s_cc[, 2], p_completos = s_cc[, 4],
                     B_imputacao = res$B[match(rownames(s_cc), res$termo)], n_completos = stats::nobs(f_cc), row.names = NULL)
    r2 <- NULL
    if (familia == "normal") {
      rr <- vapply(fits, function(f) summary(f)$r.squared, 0)
      r2 <- tanh(mean(atanh(sqrt(rr))))^2
    }
    tab <- res
    sem <- res[res$termo != "(Intercept)", , drop = FALSE]
    efeito <- data.frame(indice = paste0(if (familia == "normal") "B: " else if (familia == "binomial") "OR: " else "IRR: ", sem$termo),
                         valor = if (familia == "normal") sem$B else sem$exp_B, ic_inf = if (familia == "normal") sem$ic_inf else sem$exp_ic_inf,
                         ic_sup = if (familia == "normal") sem$ic_sup else sem$exp_ic_sup, principal = TRUE)
    if (!is.null(r2)) efeito <- rbind(efeito, data.frame(indice = "R² combinado (Harel, z de Fisher)", valor = r2, ic_inf = NA, ic_sup = NA, principal = FALSE))
    if (any(res$fmi > .5)) avisos <- c(avisos, paste0("FMI acima de ,50 em ", paste(res$termo[res$fmi > .5], collapse = ", "),
                                                      ": metade ou mais da informação desse coeficiente vem do modelo de imputação; aumentar m e checar a sensibilidade ao MAR."))
    dif <- abs(cc$B_completos - cc$B_imputacao) / res$ep[match(cc$termo, res$termo)]
    if (any(dif > .5, na.rm = TRUE)) avisos <- c(avisos, paste0("Casos completos e imputação diferem em mais de meio EP em ", paste(cc$termo[dif > .5], collapse = ", "),
                                                                ": a deleção listwise enviesaria esse coeficiente se o mecanismo for MAR."))
    sep <- .ee_sep(decimal); en <- identical(idioma, "en")
    txt <- vapply(seq_len(nrow(sem)), function(i) paste0(sem$termo[i], ": ", if (familia == "normal")
      .ee_efeito_txt("B", sem$B[i], sem$ic_inf[i], sem$ic_sup[i], FALSE, 2, decimal, idioma, nivel) else
      .ee_efeito_txt(if (familia == "binomial") "OR" else "IRR", sem$exp_B[i], sem$exp_ic_inf[i], sem$exp_ic_sup[i], FALSE, 2, decimal, idioma, nivel),
      sep, .pr_p(sem$p[i], decimal)), character(1))
    apa <- if (en) paste0("Missing data (", .pr_num(100 * inc, 1, TRUE, decimal), "% of cases incomplete) were handled by multiple imputation by chained equations (mice; m = ", m,
                          "), with estimates pooled by Rubin's rules: ", paste(txt, collapse = "; "), ".") else
      paste0("Os dados ausentes (", .pr_num(100 * inc, 1, TRUE, decimal), "% dos casos incompletos) foram tratados por imputação múltipla por equações encadeadas (mice; m = ", m,
             "), com as estimativas combinadas pelas regras de Rubin: ", paste(txt, collapse = "; "), ".")
    decisoes <- c(decisoes, "Combinação pelas regras de Rubin com gl de Barnard e Rubin (1999), iguais ao mice::pool; FMI = fração de informação ausente. R² combinado pela média do z de Fisher de sqrt(R²) (Harel, 2009, como o mice::pool.r.squared).")
  } else {
    apa <- paste0("Imputação múltipla por equações encadeadas (mice; m = ", m, ") de ", length(vars), " variáveis com ", .pr_num(100 * inc, 1, TRUE, decimal),
                  "% dos casos incompletos. Analisar com with(res$imputacoes, ...) e combinar com mice::pool().")
  }
  .ee_resultado(paste0("Imputação múltipla", if (!is.null(formula)) paste0(" (", deparse(formula), ")") else ""), n = n, n_excluidos = 0,
                tabela = tab, efeito = efeito, decisoes = decisoes, avisos = avisos, apa = apa,
                extra = list(tabelas_extra = c(list(`Ausentes e método por variável` = aus), if (!is.null(cc)) list(`Casos completos vs imputação` = cc)),
                             imputacoes = imp, combinado = res, casos_completos = cc, pct_casos_incompletos = 100 * inc))
}
