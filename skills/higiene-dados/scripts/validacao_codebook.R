# ==============================================================================
# validacao_codebook.R  (higiene-dados, desativado na v1.5.0)
# Este modulo da v1.3.0 definia outra validar_regras() (que exigia lista) e outra
# gerar_codebook() com argumentos diferentes. Quem fazia source() de todos os
# scripts em ordem alfabetica sobrescrevia as versoes de auditoria_higiene.R e
# quebrava executar_pipeline_higiene(). As duas funcoes vivem so em
# auditoria_higiene.R; este arquivo apenas carrega aquele.
# ==============================================================================
local({
  arq <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cands <- c(if (length(arq)) dirname(normalizePath(arq)), tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL), ".", "scripts")
  d <- cands[file.exists(file.path(cands, "auditoria_higiene.R"))][1]
  if (is.na(d)) stop("validacao_codebook.R foi incorporado a auditoria_higiene.R; carregue esse arquivo.", call. = FALSE)
  sys.source(file.path(d, "auditoria_higiene.R"), envir = globalenv())
  message("validacao_codebook.R: validar_regras() e gerar_codebook() vem de auditoria_higiene.R (v1.5.0).")
})
