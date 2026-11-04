# carregar.R: ponto único de entrada da skill psicometria-r no R.
#   source("~/.local/share/psicometria-r/carregar.R")
# Carrega psicometria.R (AFC, invariância, TRI, redes) e os módulos de estatística.
# Cada módulo só exige os pacotes dele na hora da chamada (.pr_exige), então carregar
# tudo não falha por falta de pacote; checar_pacotes() lista o que falta.

local({
  # Acha o arquivo em qualquer frame de source(), inclusive dentro de suppressMessages().
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  ofiles <- Filter(function(f) grepl("carregar\\.R$", f), ofiles)
  aqui <- if (length(ofiles)) dirname(normalizePath(ofiles[[length(ofiles)]])) else NULL
  if (is.null(aqui) || !nzchar(aqui)) aqui <- path.expand("~/.local/share/psicometria-r")
  arquivos <- c("psicometria.R", "est_comum.R", "est_diferencas.R", "est_associacao.R",
                "est_processos.R", "est_modelos.R", "est_poder.R", "est_construcao.R",
                "est_redes.R", "est_tri.R", "est_vieses.R", "est_latentes.R", "est_epidemio.R",
                "est_naoparametricos.R", "est_mee.R", "est_roteamento.R",
                "est_bayes.R", "est_dif_nlr.R", "est_padroes.R",
                "est_rasch_cml.R", "est_medida.R")
  for (a in arquivos) {
    caminho <- file.path(aqui, a)
    if (file.exists(caminho)) sys.source(caminho, envir = globalenv()) else
      message("psicometria-r: arquivo ausente, pulado: ", a)
  }
  assign(".pr_dir", aqui, envir = globalenv())
})
invisible(TRUE)
