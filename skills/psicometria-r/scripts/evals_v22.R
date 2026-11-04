# evals_v22.R, psicometria-r 2.2.0: Rasch por CML (rota do eRm do curso), CSEM e equacionamento.
# O primeiro bloco reproduz a Tabela 1 do "Reporte Rasch dicotômico com eRm.docx" do curso (BLOT, N = 150).
args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
suppressMessages(source(file.path(pasta, "carregar.R")))
resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome)); resultados[[nome]] <<- ok; invisible(ok)
}
CURSO <- Sys.getenv("PR_CURSO_RAIZ", "~/Desktop/Psicometria-R")
BLOT <- file.path(CURSO, "03_Teoria_de_Resposta_ao_Item_TRI", "Introdução à Teoria de Resposta ao Item (TRI)", "Materiais_e_Scripts", "Modulo", "BLOT.sav")
if (file.exists(BLOT) && requireNamespace("haven", quietly = TRUE)) {
  b <- as.data.frame(haven::zap_labels(haven::read_sav(BLOT))); r <- rasch_cml(b, names(b)); it <- r$posthoc
  # valores da Tabela 1 do relatório do curso (eRm), itens 1, 6, 12, 13, 21, 30 e 35
  ref <- data.frame(i = c(1, 6, 12, 13, 21, 30, 35), dif = c(-0.767, -2.433, -1.766, 0.995, 2.294, 1.066, -0.293),
                    ep = c(0.256, 0.459, 0.354, 0.188, 0.194, 0.187, 0.227), out = c(0.673, 0.779, 0.233, 1.315, 1.659, 1.144, 0.723),
                    inf = c(0.966, 1.007, 0.664, 1.162, 1.248, 1.185, 0.915), out_t = c(-0.858, 0.031, -1.598, 2.002, 3.523, 1.011, -0.957),
                    inf_t = c(-0.156, 0.145, -1.316, 1.971, 2.529, 2.248, -0.630))
  teste("curso, Rasch BLOT: dificuldades e EP iguais ao eRm do relatório (tol. .0015)", max(abs(it$dificuldade[ref$i] - ref$dif)) < .0015 && max(abs(it$ep[ref$i] - ref$ep)) < .0015)
  teste("curso, Rasch BLOT: infit e outfit MNSQ iguais ao eRm (tol. .0015)", max(abs(it$outfit_msq[ref$i] - ref$out)) < .0015 && max(abs(it$infit_msq[ref$i] - ref$inf)) < .0015)
  teste("curso, Rasch BLOT: infit e outfit t iguais ao eRm (tol. .0015)", max(abs(it$outfit_t[ref$i] - ref$out_t)) < .0015 && max(abs(it$infit_t[ref$i] - ref$inf_t)) < .0015)
  teste("curso, Rasch BLOT: soma das dificuldades = 0 e 3 casos com escore extremo fora do ajuste", abs(sum(it$dificuldade)) < 1e-8 && r$tabela$valor[2] == 147)
} else cat("PULADO  banco BLOT do curso ausente\n")
set.seed(1); n <- 1500; th <- stats::rnorm(n); d <- seq(-1.5, 1.5, length = 10)
Y <- as.data.frame((matrix(stats::runif(n * 10), n) < stats::plogis(outer(th, d, "-"))) * 1)
rs <- rasch_cml(Y, names(Y))
teste("Rasch CML simulado: recupera as dificuldades (erro < .1) e Andersen não rejeita (p > .05)", max(abs(rs$posthoc$dificuldade - (d - mean(d)))) < .1 && rs$tabela$valor[5] > .05)
Y2 <- Y; Y2$V10 <- (stats::runif(n) < stats::plogis(2.5 * (th - d[10]))) * 1; r2 <- rasch_cml(Y2, names(Y2))
teste("Rasch CML: item com discriminação 2,5 faz o Andersen rejeitar e tem outfit < .7 (sobreajuste)", r2$tabela$valor[5] < .001 && r2$posthoc$outfit_msq[10] < .7)
# CSEM
Yb <- as.data.frame((matrix(stats::runif(400 * 20), 400) < .6) * 1); cb <- csem(Yb, names(Yb))
teste("CSEM binomial: sqrt(x(n-x)/(n-1)) e máximo no meio da escala", isTRUE(all.equal(cb$tabela$csem, sqrt((0:20) * (20 - 0:20) / 19))) && which.max(cb$tabela$csem) == 11)
ct <- csem(Yb, names(Yb), "thorndike"); dif <- rowSums(Yb[, seq(1, 20, 2)]) - rowSums(Yb[, seq(2, 20, 2)])
teste("CSEM de Thorndike: SEM geral = DP da diferença entre metades (eq. 72)", abs(ct$sem_geral - stats::sd(dif)) < 1e-12)
# Equacionamento
set.seed(2); x <- stats::rbinom(2000, 40, .55); y <- stats::rbinom(2000, 40, .60)
el <- equacionar(x, y, "linear", 0:40, 0:40)
teste("equacionamento linear: X convertida tem média e DP de Y", abs(el$posthoc$media[2] - mean(y)) < 1e-10 && abs(el$posthoc$dp[2] - stats::sd(y)) < 1e-10)
em <- equacionar(x, y, "media", 0:40, 0:40)
teste("equacionamento pela média: deslocamento constante = média(Y) - média(X)", isTRUE(all.equal(diff(em$tabela$equivalente_y), rep(1, 40))) && abs(em$tabela$equivalente_y[1] - (mean(y) - mean(x))) < 1e-10)
ee <- equacionar(x, y, "equipercentil", 0:40, 0:40); ei <- equacionar(x, x, "equipercentil", 0:40, 0:40)
obs <- sort(unique(x)); obs <- obs[obs > min(x) & obs < max(x)]
teste("equipercentil: X contra si mesma devolve a identidade nos escores observados e a média convertida fica a < .05 de média(Y)",
  max(abs(ei$tabela$equivalente_y[obs + 1] - obs)) < 1e-8 && abs(ee$posthoc$media[2] - mean(y)) < .05 && !is.unsorted(ee$tabela$equivalente_y))
cat(sprintf("\nTotal evals_v22 2.2.0: %d | Passaram: %d | Falharam: %d\n", length(resultados), sum(unlist(resultados)), sum(!unlist(resultados))))
