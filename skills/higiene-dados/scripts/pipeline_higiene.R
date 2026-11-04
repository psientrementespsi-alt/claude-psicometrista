# ==============================================================================
# pipeline_higiene.R
# Script Mestre de Higienizacao, Triagem e Purificacao de Dados
# Skill: higiene-dados (v1.4.0)
# ==============================================================================

# Carrega os modulos auxiliares da skill com resolucao agnostica de ambiente
obter_caminho_script <- function() {
  # 1. Verifica se esta sendo chamado via source()
  for (i in rev(seq_along(sys.frames()))) {
    fr <- sys.frame(i)
    if (!is.null(fr$ofile)) {
      dir_candidato <- normalizePath(dirname(fr$ofile))
      if (file.exists(file.path(dir_candidato, "screening_cier.R"))) {
        return(dir_candidato)
      }
    }
  }
  # 2. Verifica se foi chamado diretamente via Rscript --file=...
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dir_candidato <- normalizePath(dirname(sub("^--file=", "", file_arg[1])))
    if (file.exists(file.path(dir_candidato, "screening_cier.R"))) {
      return(dir_candidato)
    }
  }
  # 3. Fallbacks de busca nos diretorios padrao do Gemini e Claude
  for (d in c(".", 
              "~/.gemini/config/skills/higiene-dados/scripts", 
              "~/.claude/skills/higiene-dados/scripts")) {
    if (file.exists(file.path(d, "screening_cier.R"))) return(normalizePath(d))
  }
  "."
}
caminho_script <- obter_caminho_script()
source(file.path(caminho_script, "screening_cier.R"))
source(file.path(caminho_script, "imputacao_mice.R"))
source(file.path(caminho_script, "auditoria_higiene.R"))
source(file.path(caminho_script, "ajuste_pessoa.R"))

#' Higieniza nomes de variaveis para snake_case e remove caracteres especiais
#' 
#' @param nomes Vetor de caracteres com os nomes brutos das colunas.
#' @return Vetor em formato snake_case limpo.
sanitizar_nomes_snake <- function(nomes) {
  # Substitui acentos comuns de portugues preservando letras
  nomes <- chartr("áéíóúâêîôûãõàèìòùäëïöüçÁÉÍÓÚÂÊÎÔÛÃÕÀÈÌÒÙÄËÏÖÜÇ",
                  "aeiouaeiouaoaeiouaeioucaeiouaeiouaoaeiouaeiouc", nomes)
  # Converte espacos, pontos e hifens em underscore
  nomes <- gsub("[^A-Za-z0-9_]", "_", nomes)
  # Remove underscores repetidos
  nomes <- gsub("_+", "_", nomes)
  # Remove underscores nas pontas e poem em minusculo
  nomes <- tolower(gsub("^_|_$", "", nomes))
  nomes
}

#' Executa o pipeline de higiene de dados (v2, higiene-dados 1.3.0)
#'
#' Etapas, nesta ordem, cada uma registrada na trilha de auditoria: registros vazios; recodificacao e valores
#' fora da escala (viram NA; contados em valores, nao em casos); regras declaradas (so relatadas); registros sem
#' resposta aos itens; duplicatas (so entre registros com resposta); omissao acima do limiar; teste de Little;
#' triagem C/IER por familias (repeticao, tempo, psychsyn, even-odd, Mahalanobis, atencao; exclusao com 2 ou mais);
#' atipicos multivariados por FastMCD (MCD75; assinalados, nao excluidos); MICE opcional.
#'
#' @param fatores lista nomeada com os itens de cada subescala (nomes como no banco bruto ou ja em snake_case).
#'   Usada no even-odd, no Mahalanobis por fator e no FastMCD sobre escores de subescala.
#' @param mcd_em "subescalas" (soma dos itens de cada fator; exige 'fatores') ou "itens"
#' @param regras vetor nomeado de expressoes para validar_regras(), avaliadas apos a recodificacao
#' @param regras_excluir nomes das regras cujos violadores saem do banco (criterio de inclusao, ex. consentimento);
#'   as demais so sao relatadas. NA na avaliacao de uma regra de exclusao nao exclui.
#' @param salvar_em pasta opcional: grava estado_<etapa>.rds a cada etapa (permite rodar em blocos e retomar)
#' @param etapas_ate nome da ultima etapa a rodar ("cier", "mcd" ...), para dividir a execucao em blocos
#' @param retomar_de caminho de um estado_<etapa>.rds gravado antes; o pipeline continua dali
#' @return lista: dados_purificados, log_descartes, tabela_fluxo_jars, trilha, tabela1_md, texto_metodo,
#'   regras, codebook, metricas_cier, teste_mcar, mcd, objeto_mice
executar_pipeline_higiene <- function(df_bruto,
                                      cols_itens,
                                      col_tempo = NULL,
                                      col_id = NULL,
                                      range_valido_itens = c(1, 5),
                                      limiar_omissao_caso = 0.20,
                                      executar_mice_se_ausente = FALSE,
                                      ignorar_todos_zeros = TRUE,
                                      mapa_niveis_likert = NULL,
                                      semente = 2026,
                                      fatores = NULL,
                                      alpha_mcd = 0.75,
                                      mcd_em = c("subescalas", "itens"),
                                      regras = NULL,
                                      regras_excluir = NULL,
                                      usar_mahad_cier = TRUE,
                                      p_mahad = 0.001,
                                      itens_atencao = NULL,
                                      esperado_atencao = NULL,
                                      incorreto_atencao = NULL,
                                      max_falhas_atencao = 0,
                                      rodar_mcar = TRUE,
                                      gerar_codebook_saida = TRUE,
                                      salvar_em = NULL,
                                      etapas_ate = NULL,
                                      retomar_de = NULL,
                                      verboso = TRUE) {
  mcd_em <- match.arg(mcd_em)
  diz <- function(...) if (isTRUE(verboso)) cat(sprintf(...), "\n", sep = "")
  ETAPAS <- c("sanitizar", "recodificar", "regras", "sem_resposta", "duplicatas", "omissao", "mcar", "cier", "mcd", "mice")
  salvar <- function(etapa, st) if (!is.null(salvar_em)) { dir.create(salvar_em, showWarnings = FALSE, recursive = TRUE)
    saveRDS(st, file.path(salvar_em, sprintf("estado_%s.rds", etapa))) }
  parar <- function(etapa) !is.null(etapas_ate) && identical(etapa, etapas_ate)

  PARAMS <- c("range_valido_itens", "limiar_omissao_caso", "executar_mice_se_ausente", "ignorar_todos_zeros", "mapa_niveis_likert",
              "semente", "alpha_mcd", "mcd_em", "regras", "regras_excluir", "usar_mahad_cier", "p_mahad", "itens_atencao",
              "esperado_atencao", "incorreto_atencao", "max_falhas_atencao", "rodar_mcar", "gerar_codebook_saida")
  if (!is.null(retomar_de)) {
    st <- readRDS(retomar_de); diz("Retomando depois da etapa '%s'", st$ultima)
    # os parametros da primeira execucao valem na retomada (evita trocar escala, cortes ou semente no meio)
    list2env(st$params, envir = environment())
  } else {
    st <- list(ultima = NULL, trilha = NULL, log = list(), cier = NULL, mcar = NULL, mcd = NULL, regras = NULL, mice = NULL)
    st$n_inicial <- nrow(df_bruto)
    st$cols <- sanitizar_nomes_snake(cols_itens)
    st$col_tempo <- if (!is.null(col_tempo)) sanitizar_nomes_snake(col_tempo) else NULL
    st$col_id <- if (!is.null(col_id)) sanitizar_nomes_snake(col_id) else NULL
    st$fatores <- if (!is.null(fatores)) lapply(fatores, sanitizar_nomes_snake) else NULL
    st$df <- df_bruto
    st$params <- mget(PARAMS, envir = environment())
  }
  feito <- function(e) !is.null(st$ultima) && match(e, ETAPAS) <= match(st$ultima, ETAPAS)
  fim <- function(e) { st$ultima <<- e; salvar(e, st); parar(e) }
  cols <- st$cols

  if (!feito("sanitizar")) {
    df <- st$df; colnames(df) <- sanitizar_nomes_snake(colnames(df))
    if (any(duplicated(colnames(df)))) stop("Nomes de coluna repetidos depois do snake_case: ", paste(unique(colnames(df)[duplicated(colnames(df))]), collapse = ", "))
    falt <- setdiff(cols, colnames(df)); if (length(falt)) stop("Itens ausentes no banco: ", paste(falt, collapse = ", "))
    n0 <- nrow(df)
    vaz <- rowSums(as.matrix(as.data.frame(lapply(df, function(v) !is.na(v) & trimws(as.character(v)) != "")))) == 0
    df <- df[!vaz, , drop = FALSE]
    st$trilha <- registrar_passo(st$trilha, "Registros vazios", n0, nrow(df), "todas as colunas vazias")
    st$df <- df; diz("[sanitizar] %d registros; %d vazios removidos", n0, sum(vaz))
    if (fim("sanitizar")) return(invisible(st))
  }
  if (!feito("recodificar")) {
    df <- st$df
    df[, cols] <- recodificar_likert_labels(df[, cols, drop = FALSE], mapa_niveis = mapa_niveis_likert)
    fora_total <- 0L
    for (cc in cols) { df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
      fora <- !is.na(df[[cc]]) & (df[[cc]] < range_valido_itens[1] | df[[cc]] > range_valido_itens[2] | df[[cc]] != round(df[[cc]]))
      fora_total <- fora_total + sum(fora); df[[cc]][fora] <- NA_real_ }
    if (!is.null(st$col_tempo)) df[[st$col_tempo]] <- suppressWarnings(as.numeric(df[[st$col_tempo]]))
    st$log$valores_fora_faixa <- fora_total
    st$trilha <- registrar_passo(st$trilha, "Recodificacao", nrow(df), nrow(df),
      sprintf("rotulos para %d-%d; %d valores fora da escala ou nao inteiros viraram NA", range_valido_itens[1], range_valido_itens[2], fora_total))
    st$df <- df; diz("[recodificar] %d valores fora da escala -> NA", fora_total)
    if (fim("recodificar")) return(invisible(st))
  }
  if (!feito("regras")) {
    if (!is.null(regras)) { st$regras <- validar_regras(st$df, regras); diz("[regras] %d regras; %d violacoes no total", nrow(st$regras), sum(st$regras$violacoes, na.rm = TRUE)) }
    if (!is.null(regras_excluir)) {
      if (is.null(regras) || !all(regras_excluir %in% names(regras))) stop("regras_excluir precisa nomear regras presentes em 'regras'.")
      df <- st$df; n0 <- nrow(df); viol <- rep(FALSE, n0)
      for (nm in regras_excluir) { v <- eval(str2lang(regras[[nm]]), envir = df, enclos = baseenv()); viol <- viol | (!is.na(v) & !v) }
      df <- df[!viol, , drop = FALSE]
      st$trilha <- registrar_passo(st$trilha, "Criterios de inclusao", n0, nrow(df), paste(sprintf("%s: %s", regras_excluir, regras[regras_excluir]), collapse = "; "))
      st$df <- df; diz("[regras] %d excluidos por criterio de inclusao", sum(viol))
    }
    if (fim("regras")) return(invisible(st))
  }
  if (!feito("sem_resposta")) {
    df <- st$df; n0 <- nrow(df); sem <- rowSums(!is.na(df[, cols, drop = FALSE])) == 0
    df <- df[!sem, , drop = FALSE]
    st$trilha <- registrar_passo(st$trilha, "Sem resposta aos itens", n0, nrow(df), "nenhum item respondido (abandono antes da escala)")
    st$df <- df; diz("[sem_resposta] %d removidos", sum(sem))
    if (fim("sem_resposta")) return(invisible(st))
  }
  if (!feito("duplicatas")) {
    df <- st$df; n0 <- nrow(df)
    subst <- setdiff(colnames(df), c(st$col_id, st$col_tempo))
    dup <- duplicated(df[, subst, drop = FALSE])
    if (isTRUE(ignorar_todos_zeros)) {   # protocolo todo no minimo nunca conta como duplicata (inventarios clinicos)
      M <- as.matrix(df[, cols, drop = FALSE])
      todos_min <- apply(M, 1, function(r) { r <- r[!is.na(r)]; length(r) > 0 && all(r == range_valido_itens[1]) })
      dup <- dup & !todos_min
    }
    df <- df[!dup, , drop = FALSE]
    st$log$duplicatas <- sum(dup)
    crit <- if (length(subst) > length(cols)) "registro identico em todas as colunas substantivas (exceto ID e tempo)" else "padrao de itens identico (todos no minimo protegidos)"
    if (sum(dup) > 0 && length(subst) - length(cols) < 3) {
      crit <- paste0(crit, "; ATENCAO: menos de 3 colunas alem dos itens, entao padroes iguais por acaso (ou retilineos) podem ter sido contados como duplicata: conferir com ID ou carimbo de tempo")
      warning("Duplicatas detectadas com menos de 3 colunas alem dos itens: conferir se sao registros repetidos de fato.", call. = FALSE)
    }
    st$trilha <- registrar_passo(st$trilha, "Duplicatas", n0, nrow(df), crit)
    st$df <- df; diz("[duplicatas] %d removidas", sum(dup))
    if (fim("duplicatas")) return(invisible(st))
  }
  if (!feito("omissao")) {
    df <- st$df; n0 <- nrow(df); po <- rowSums(is.na(df[, cols, drop = FALSE])) / length(cols)
    exc <- po > limiar_omissao_caso; df <- df[!exc, , drop = FALSE]
    st$log$omissao_excessiva <- sum(exc)
    st$trilha <- registrar_passo(st$trilha, "Omissao excessiva", n0, nrow(df), sprintf("mais de %.0f%% dos itens em branco", 100 * limiar_omissao_caso), "Newman (2014)")
    st$df <- df; diz("[omissao] %d removidos", sum(exc))
    if (fim("omissao")) return(invisible(st))
  }
  if (!feito("mcar")) {
    if (isTRUE(rodar_mcar)) { st$mcar <- teste_little_mcar(st$df[, cols, drop = FALSE]); diz("[mcar] %s", st$mcar$conclusao) }
    if (fim("mcar")) return(invisible(st))
  }
  if (!feito("cier")) {
    df <- st$df; n0 <- nrow(df); fat <- st$fatores
    tempo <- if (!is.null(st$col_tempo) && st$col_tempo %in% colnames(df)) df[[st$col_tempo]] else NULL
    att <- if (!is.null(itens_atencao)) df[, sanitizar_nomes_snake(itens_atencao), drop = FALSE] else NULL
    tr <- triagem_cier_completa(df[, cols, drop = FALSE], tempo_total_seg = tempo, n_itens_total = length(cols),
      corte_seg_por_item = 2, max_longstring_proporcao = 0.5, ignorar_todos_zeros = ignorar_todos_zeros, val_zero = range_valido_itens[1],
      fatores_evenodd = if (!is.null(fat) && length(fat) >= 2) fat else NULL,
      usar_mahad = usar_mahad_cier, p_mahad = p_mahad, fatores_mahad = fat,
      itens_atencao = att, esperado_atencao = esperado_atencao, incorreto_atencao = incorreto_atencao, max_falhas_atencao = max_falhas_atencao)
    exc <- tr$excluir_por_cier
    st$cier <- list(contagens = colSums(tr[, grep("^flag_|^familia_", names(tr)), drop = FALSE]), decisoes = attr(tr, "decisoes"),
                    avisos = attr(tr, "avisos"), n_excluidos = sum(exc))
    df <- df[!exc, , drop = FALSE]; st$log$respostas_descuidadas_cier <- sum(exc)
    st$trilha <- registrar_passo(st$trilha, "Respostas descuidadas", n0, nrow(df), "2 ou mais familias de indices C/IER", "Curran (2016); Ward & Meade (2023)")
    st$df <- df; diz("[cier] %d excluidos", sum(exc))
    if (fim("cier")) return(invisible(st))
  }
  if (!feito("mcd")) {
    df <- st$df; fat <- st$fatores
    X <- if (mcd_em == "subescalas" && !is.null(fat)) sapply(fat, function(ii) rowSums(df[, ii, drop = FALSE])) else df[, cols, drop = FALSE]
    md <- calc_mahad_robust(X, p_val = 0.001, semente = semente, alpha_mcd = alpha_mcd)
    df$flag_outlier_multivariado <- md$is_outlier
    st$mcd <- list(n_assinalados = sum(md$is_outlier), estimador = md$estimador, corte = md$corte_chi2,
                   base = if (mcd_em == "subescalas" && !is.null(fat)) "escores de subescala" else "itens")
    st$log$outliers_assinalados <- sum(md$is_outlier)
    st$trilha <- registrar_passo(st$trilha, "Atipicos multivariados (assinalados)", nrow(df), nrow(df),
      sprintf("%s sobre %s, p < .001; %d assinalados e retidos", md$estimador, st$mcd$base, sum(md$is_outlier)), "Leys et al. (2019)")
    st$df <- df; diz("[mcd] %d assinalados (%s)", sum(md$is_outlier), md$estimador)
    if (fim("mcd")) return(invisible(st))
  }
  if (!feito("mice")) {
    nas <- sum(is.na(st$df[, cols, drop = FALSE]))
    if (isTRUE(executar_mice_se_ausente) && nas > 0) st$mice <- executar_mice_psicometria(st$df[, cols, drop = FALSE], m = 20, seed = semente)
    if (fim("mice")) return(invisible(st))
  }
  tr <- st$trilha; nf <- nrow(st$df)
  fluxo <- data.frame(etapa = tr$etapa, n_excluido = tr$n_removido, n_restante = tr$n_depois,
                      pct_amostra_inicial = round(100 * tr$n_depois / st$n_inicial, 1), stringsAsFactors = FALSE)
  cb <- if (isTRUE(gerar_codebook_saida)) gerar_codebook(st$df) else NULL
  texto <- texto_metodo_higiene(tr, st$cier$decisoes, st$mcar, st$mcd$n_assinalados, st$mcd$estimador)
  fluxo_mmd <- fluxograma_mermaid(tr)
  diz("PURIFICACAO CONCLUIDA: N inicial = %d -> N analitico = %d", st$n_inicial, nf)
  list(dados_purificados = st$df, log_descartes = st$log, tabela_fluxo_jars = fluxo, trilha = tr,
       tabela1_md = tabela_fluxo_md(tr), fluxograma_mermaid = fluxo_mmd, texto_metodo = texto, regras = st$regras, codebook = cb,
       metricas_cier = st$cier, teste_mcar = st$mcar, mcd = st$mcd, objeto_mice = st$mice)
}
