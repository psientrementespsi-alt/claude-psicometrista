# ==============================================================================
# auditoria_higiene.R  (higiene-dados v1.3.0)
# Regras de validacao declarativas, codebook agregado e trilha de auditoria.
# Privacidade: nenhuma funcao devolve linha de participante; so contagens e
# resumos. Categorias com menos de 'n_min_celula' casos saem suprimidas.
# Inspiracao: validate (van der Loo & de Jonge, 2021, JSS 97(10)), codebook
# (Arslan, 2019, AMPPS), lumberjack/tidylog (registro de cada passo).
# ==============================================================================

#' Confere regras declaradas como expressoes do R sobre o banco
#'
#' @param df data.frame
#' @param regras vetor nomeado de expressoes em texto, ex. c(idade = "idade >= 18 & idade <= 99")
#' @return data.frame: regra, expressao, avaliados, violacoes, nao_avaliaveis, pct_violacao. Usa
#'   validate::confront quando o pacote existe; senao, avalia em base R. Nenhuma linha do banco e devolvida.
validar_regras <- function(df, regras) {
  if (is.null(names(regras)) || any(!nzchar(names(regras)))) names(regras) <- paste0("R", seq_along(regras))
  if (requireNamespace("validate", quietly = TRUE)) {
    v <- tryCatch({
      vv <- do.call(validate::validator, lapply(regras, str2lang))
      sm <- validate::summary(validate::confront(df, vv))
      data.frame(regra = names(regras), expressao = unname(regras), avaliados = sm$items - sm$nNA,
                 violacoes = sm$fails, nao_avaliaveis = sm$nNA, motor = "validate", stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (!is.null(v)) { v$pct_violacao <- round(100 * v$violacoes / pmax(1, v$avaliados), 1); return(v) }
  }
  out <- lapply(names(regras), function(nm) {
    r <- tryCatch(eval(str2lang(regras[[nm]]), envir = df, enclos = baseenv()), error = function(e) e)
    if (inherits(r, "error")) return(data.frame(regra = nm, expressao = regras[[nm]], avaliados = NA_integer_, violacoes = NA_integer_,
      nao_avaliaveis = NA_integer_, motor = paste("erro:", conditionMessage(r)), stringsAsFactors = FALSE))
    if (length(r) == 1L) r <- rep(r, nrow(df))
    data.frame(regra = nm, expressao = regras[[nm]], avaliados = sum(!is.na(r)), violacoes = sum(!r, na.rm = TRUE),
               nao_avaliaveis = sum(is.na(r)), motor = "base R", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out); out$pct_violacao <- round(100 * out$violacoes / pmax(1, out$avaliados), 1); out
}

#' Codebook agregado do banco (sem dado individual)
#'
#' @param df data.frame
#' @param rotulos vetor nomeado opcional com o rotulo de cada variavel
#' @param n_min_celula categorias com menos casos que isto aparecem como "<n" (padrao 5)
#' @param max_niveis numero maximo de categorias listadas por variavel
#' @param arquivo_md caminho opcional para gravar o codebook em Markdown
gerar_codebook <- function(df, rotulos = NULL, n_min_celula = 5, max_niveis = 8, arquivo_md = NULL) {
  linhas <- lapply(names(df), function(v) {
    x <- df[[v]]; n <- length(x); na <- sum(is.na(x)); nu <- length(unique(x[!is.na(x)]))
    tipo <- if (inherits(x, c("Date", "POSIXt"))) "data" else if (is.numeric(x)) "numerica" else if (is.factor(x)) "fator" else if (is.logical(x)) "logica" else "texto"
    livre <- tipo == "texto" && nu >= max(2 * n_min_celula, 0.5 * (n - na))
    resumo <- ""
    if (tipo == "numerica" && n - na > 0) {
      resumo <- sprintf("M = %.2f; DP = %.2f; Md = %.2f; min = %s; max = %s", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE),
                        stats::median(x, na.rm = TRUE), format(min(x, na.rm = TRUE)), format(max(x, na.rm = TRUE)))
      if (n - na < n_min_celula) resumo <- sprintf("resumo suprimido (menos de %d valores)", n_min_celula)
    } else if (tipo == "data" && n - na > 0) {
      resumo <- sprintf("de %s a %s", format(min(x, na.rm = TRUE), "%Y-%m"), format(max(x, na.rm = TRUE), "%Y-%m"))
    } else if (livre) {
      resumo <- "texto livre: conteudo omitido"
    } else if (n - na > 0) {
      tb <- sort(table(x, useNA = "no"), decreasing = TRUE)
      lab <- ifelse(tb < n_min_celula, sprintf("%s: <%d", names(tb), n_min_celula), sprintf("%s: %d", names(tb), as.integer(tb)))
      if (length(lab) > max_niveis) lab <- c(lab[seq_len(max_niveis)], sprintf("(+%d categorias)", length(lab) - max_niveis))
      resumo <- paste(lab, collapse = "; ")
    }
    data.frame(variavel = v, rotulo = if (!is.null(rotulos) && v %in% names(rotulos)) rotulos[[v]] else "",
               tipo = if (livre) "texto livre" else tipo, n_validos = n - na, n_ausentes = na,
               pct_ausente = round(100 * na / max(1, n), 1), n_valores_distintos = if (livre) NA_integer_ else nu,
               resumo = resumo, stringsAsFactors = FALSE)
  })
  cb <- do.call(rbind, linhas)
  if (!is.null(arquivo_md)) {
    md <- c(sprintf("# Codebook (gerado em %s; %d casos, %d variaveis)", format(Sys.Date()), nrow(df), ncol(df)), "",
            "| Variavel | Rotulo | Tipo | Validos | Ausentes (%) | Distintos | Resumo |", "|---|---|---|---:|---:|---:|---|",
            sprintf("| %s | %s | %s | %d | %d (%.1f) | %s | %s |", cb$variavel, cb$rotulo, cb$tipo, cb$n_validos, cb$n_ausentes,
                    cb$pct_ausente, ifelse(is.na(cb$n_valores_distintos), "", cb$n_valores_distintos), gsub("\\|", "/", cb$resumo)),
            "", sprintf("Categorias com menos de %d casos aparecem como <%d.", n_min_celula, n_min_celula))
    writeLines(md, arquivo_md, useBytes = TRUE)
  }
  cb
}

#' Trilha de auditoria: acrescenta um passo
#'
#' @param trilha data.frame anterior (ou NULL)
#' @param etapa nome curto da etapa
#' @param n_antes,n_depois casos antes e depois
#' @param decisao texto da regra aplicada (corte, criterio)
#' @param referencia fonte do criterio
registrar_passo <- function(trilha = NULL, etapa, n_antes, n_depois, decisao = "", referencia = "") {
  linha <- data.frame(ordem = if (is.null(trilha)) 1L else nrow(trilha) + 1L, etapa = etapa, n_antes = as.integer(n_antes),
                      n_depois = as.integer(n_depois), n_removido = as.integer(n_antes - n_depois), decisao = decisao,
                      referencia = referencia, momento = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE)
  if (is.null(trilha)) linha else rbind(trilha, linha)
}

#' Tabela de fluxo (Tabela 1 do JARS-Quant) em Markdown a partir da trilha
tabela_fluxo_md <- function(trilha, titulo = "Fluxo da amostra da coleta a amostra analitica") {
  c("**Tabela 1**", "", paste0("*", titulo, "*"), "", "| Etapa | Criterio | *n* removido | *n* restante |", "|---|---|---:|---:|",
    sprintf("| %s | %s | %d | %d |", trilha$etapa, gsub("\\|", "/", trilha$decisao), trilha$n_removido, trilha$n_depois), "",
    "*Nota.* Contagens geradas pela trilha de auditoria da higiene-dados; nenhum caso individual e exibido.")
}

#' Fluxograma de participantes em sintaxe Mermaid (JARS-Quant / CONSORT)
#'
#' Gera diagrama de fluxo direcionado em sintaxe Mermaid pronto para visualizacao e relato.
#'
#' @param trilha data.frame da trilha de auditoria do pipeline (ou lista retornada por executar_pipeline_higiene)
#' @param titulo Titulo do diagrama (comentario no bloco Mermaid)
#' @return Vetor de caracteres contendo a sintaxe Mermaid
fluxograma_mermaid <- function(trilha, titulo = "Fluxo de Participantes (JARS-Quant / CONSORT)") {
  if (is.list(trilha) && !is.data.frame(trilha) && !is.null(trilha$trilha)) {
    trilha <- trilha$trilha
  }
  if (!is.data.frame(trilha) || !all(c("etapa", "n_antes", "n_depois", "n_removido") %in% names(trilha))) {
    stop("trilha deve ser um data.frame de trilha de auditoria ou resultado de executar_pipeline_higiene()")
  }
  passos <- trilha[trilha$n_removido > 0, , drop = FALSE]
  n0 <- trilha$n_antes[1]
  nf <- utils::tail(trilha$n_depois, 1)
  
  linhas <- c(
    "```mermaid",
    "graph TD",
    sprintf("    %% %s", titulo),
    sprintf('    N0["Registros recebidos na coleta<br>(N = %d)"]', n0)
  )
  
  if (nrow(passos) == 0) {
    linhas <- c(linhas, sprintf('    N0 --> NF["Amostra analítica final<br>(N = %d)"]', nf))
  } else {
    no_atual <- "N0"
    for (i in seq_len(nrow(passos))) {
      p <- passos[i, ]
      no_depois <- paste0("N", i)
      no_excl <- paste0("E", i)
      crit <- gsub('["\n]', ' ', p$decisao)
      if (nzchar(p$referencia)) crit <- paste0(crit, " (", p$referencia, ")")
      
      linhas <- c(linhas,
        sprintf('    %s -.->|n = %d| %s["Excluídos: %s<br>Critério: %s"]',
                no_atual, p$n_removido, no_excl, p$etapa, crit)
      )
      if (i == nrow(passos)) {
        linhas <- c(linhas,
          sprintf('    %s --> NF["Amostra analítica final<br>(N = %d)"]', no_atual, nf)
        )
      } else {
        linhas <- c(linhas,
          sprintf('    %s --> %s["Amostra retida<br>(N = %d)"]', no_atual, no_depois, as.integer(p$n_depois))
        )
        no_atual <- no_depois
      }
    }
  }
  linhas <- c(linhas,
    "    classDef default fill:#f9f9f9,stroke:#333,stroke-width:1px;",
    "    classDef final fill:#d4edda,stroke:#28a745,stroke-width:2px;",
    "    classDef excl fill:#f8d7da,stroke:#dc3545,stroke-width:1px,stroke-dasharray: 5 5;",
    "    class NF final;"
  )
  for (i in seq_len(nrow(passos))) {
    linhas <- c(linhas, sprintf("    class E%d excl;", i))
  }
  linhas <- c(linhas, "```")
  class(linhas) <- c("pr_mermaid", "character")
  linhas
}

#' Paragrafo de Metodo (portugues) a partir da trilha e das decisoes da triagem
texto_metodo_higiene <- function(trilha, decisoes_cier = NULL, mcar = NULL, n_atipicos = NULL, estimador_mcd = NULL) {
  n0 <- trilha$n_antes[1]; nf <- utils::tail(trilha$n_depois, 1)
  passos <- trilha[trilha$n_removido > 0, ]
  frase_passos <- if (nrow(passos)) paste(sprintf("%s (%d)", tolower(passos$etapa), passos$n_removido), collapse = "; ") else "nenhuma exclusao"
  t <- sprintf("Dos %d registros recebidos, %d compuseram a amostra analitica. As exclusoes, na ordem em que foram aplicadas, foram: %s.", n0, nf, frase_passos)
  if (!is.null(decisoes_cier)) t <- paste0(t, " A triagem de respostas descuidadas seguiu familias de indices, e um caso so foi excluido quando assinalado em duas ou mais familias (Curran, 2016; Ward & Meade, 2023): ",
    paste(decisoes_cier, collapse = "; "), ".")
  p3 <- function(p) if (p < .001) "< .001" else paste("=", sub("^0", "", sprintf("%.3f", p)))
  if (!is.null(mcar) && !is.na(mcar$p_valor) && isTRUE(mcar$df > 0))
    t <- paste(t, sprintf("O teste de Little (1988) %s, χ²(%d) = %.2f, p %s.",
      if (mcar$p_valor >= .05) "nao rejeitou a hipotese de dados ausentes completamente ao acaso" else "rejeitou a hipotese de dados ausentes completamente ao acaso",
      as.integer(mcar$df), mcar$qui_quadrado, p3(mcar$p_valor)))
  if (!is.null(mcar) && isTRUE(mcar$df == 0)) t <- paste(t, "Nao restaram dados ausentes nos itens depois das exclusoes.")
  if (!is.null(n_atipicos)) {
    a <- suppressWarnings(as.numeric(sub(".*alpha = ([0-9.]+).*", "\\1", if (is.null(estimador_mcd)) "" else estimador_mcd)))
    est <- if (is.na(a)) "FastMCD" else paste0("FastMCD com alpha = ", sub("^0", "", sprintf("%.2f", a)))
    t <- paste(t, sprintf("A distancia robusta de Mahalanobis (%s; Leys et al., 2019), com corte em χ² a p < .001, assinalou %d %s, %s no banco para a analise de sensibilidade.",
      est, n_atipicos, if (n_atipicos == 1) "caso" else "casos", if (n_atipicos == 1) "retido" else "retidos"))
  }
  t
}
