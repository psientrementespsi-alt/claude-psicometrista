# evals_vieses.R, testes do módulo de vieses de resposta (est_vieses.R, v1.4)
#
# Rodar da pasta da skill:  Rscript scripts/evals_vieses.R
# Dados simulados com semente fixa: 3 traços com 4 itens cada (2 positivos e 2 negativos,
# pares semânticos), aquiescência com variância .25 na resposta latente (.20 depois da
# padronização da resposta latente pelo WLSMV), 5 categorias. Nenhum dado de participante.

args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
source(file.path(pasta, "carregar.R"))
options(mc.cores = 1L)
resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome))
  resultados[[nome]] <<- ok
  invisible(ok)
}
t_inicio <- Sys.time()
quieto <- function(expr) suppressWarnings(suppressMessages(expr))
erro_de <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))

cat("== Aquiescência (curso, módulo de controle de aquiescência; Maydeu-Olivares & Coffman, 2006)\n")
set.seed(21); n <- 1500
th <- MASS::mvrnorm(n, rep(0, 3), matrix(c(1, .3, .2, .3, 1, .25, .2, .25, 1), 3))
acq <- stats::rnorm(n, 0, sqrt(.25))
sinal <- rep(c(1, -1, 1, -1), 3); fat <- rep(1:3, each = 4)
ys <- sapply(1:12, function(j) sinal[j] * .7 * th[, fat[j]] + acq + stats::rnorm(n, 0, sqrt(1 - .49)))
d <- as.data.frame(apply(ys, 2, function(x) as.integer(cut(x, c(-Inf, -1.2, -.4, .4, 1.2, Inf))))); names(d) <- paste0("y", 1:12)
fatores <- list(Abertura = paste0("y", 1:4), Consc = paste0("y", 5:8), Extro = paste0("y", 9:12))
polo <- stats::setNames(ifelse(sinal == 1, 1, 0), names(d))
r <- quieto(controle_aquiescencia(d, fatores, polo, minimo = 1, maximo = 5))
teste("índice = média das respostas sem inverter menos o ponto médio; correlação com a aquiescência simulada > .85", {
  abs(r$indice[1] - (mean(unlist(d[1, ])) - 3)) < 1e-12 && stats::cor(r$indice, acq) > .85 })
teste("interceptos randômicos recuperam a variância de aquiescência (.20 na métrica padronizada; tol. .03) e ajustam melhor que o modelo sem controle", {
  abs(r$variancia_aquiescencia$variancia - .20) < .03 && r$tabela$cfi[1] < .90 && r$tabela$cfi[2] > .99 })
teste("modelo de interceptos randômicos = sintaxe do curso no lavaan (mesmo qui-quadrado escalonado; tol. 1e-6)", {
  m <- paste0("Abertura =~ y1 + y2 + y3 + y4\nConsc =~ y5 + y6 + y7 + y8\nExtro =~ y9 + y10 + y11 + y12\n",
              "Acq =~ ", paste0("1*y", 1:12, collapse = " + "), "\nAbertura ~~ 0*Acq\nConsc ~~ 0*Acq\nExtro ~~ 0*Acq")
  f <- quieto(lavaan::cfa(m, d, ordered = TRUE, estimator = "WLSMV"))
  abs(lavaan::fitMeasures(f, "chisq.scaled") - r$tabela$qui2[2]) < 1e-6 })
teste("MIMIC do curso (itens no índice, fatores ortogonais ao índice): caminhos positivos em todos os itens, dos dois polos", {
  mt <- r$tabelas_extra$`MIMIC: itens no índice`; nrow(mt) == 12 && all(mt$b > 0) })
teste("recentralização sem arredondar: média de cada pessoa = 0 (tol. 1e-12); com arredondar, valores inteiros em [-2, 2]", {
  r2 <- quieto(controle_aquiescencia(d, fatores, polo, minimo = 1, maximo = 5, metodos = character(), arredondar = FALSE))
  max(abs(rowMeans(r2$dados_recentralizados))) < 1e-12 && all(unlist(r$dados_recentralizados) %in% -2:2) })
teste("análise paralela nos dados recentralizados sugere os 3 fatores simulados", {
  nf <- quieto(afe_numero_fatores(r$dados_recentralizados, names(d), n_iter = 10, hull = FALSE, map = FALSE)); nf$n_sugerido == 3 })
teste("escala desbalanceada gera aviso; fator sem item negativo para com mensagem", {
  p2 <- polo; p2["y2"] <- 1
  a <- any(grepl("não balanceada", quieto(controle_aquiescencia(d, fatores, p2, minimo = 1, maximo = 5, metodos = character()))$avisos))
  p3 <- polo; p3[c("y2", "y4")] <- 1
  a && grepl("sem item negativo", erro_de(controle_aquiescencia(d, fatores, p3, minimo = 1, maximo = 5, metodos = character()))) })
teste("privacidade: a impressão não traz linha de participante (nada de 1500 linhas, nenhum valor individual do índice)", {
  out <- capture.output(print(r))
  ind_txt <- vapply(r$indice[1:20], function(v) format(signif(v, 4)), "")
  length(out) < 120 && !any(vapply(ind_txt[abs(r$indice[1:20]) > 1e-8 & nchar(ind_txt) >= 6], function(v) any(grepl(v, out, fixed = TRUE)), TRUE)) })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
print(versoes_pacotes(c("lavaan", "psych")))
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
