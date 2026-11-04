# evals_construcao.R, testes do módulo de construção e validação de instrumentos
# (est_construcao.R)
#
# Rodar da pasta da skill:  Rscript scripts/evals_construcao.R
# Saída: uma linha PASSOU/FALHOU por teste e o resumo. Código de saída 1 se algum
# teste falhar. Nenhum dado de participante do usuário: notas de juízes digitadas a
# partir dos exemplos do curso, estruturas fatoriais simuladas com lavaan::simulateData,
# psych::sim.structure e MASS::mvrnorm com sementes fixas (dados ordinais com
# ordered_center = FALSE, como em evals_simulados.R) e os itens A1-A5 do psych::bfi.
# Com PR_DADOS_JASP apontando para a biblioteca de dados do JASP (a pasta que contém
# "6. Factor"; a busca é recursiva), rodam também os testes com o banco G Factor; sem a
# variável, eles são pulados com aviso. A análise paralela roda com 20 reamostras por
# velocidade; o padrão da função é 100. Tolerância declarada em cada teste.

args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
source(file.path(pasta, "carregar.R"))
# psych::fa.parallel e psych::polychoric, chamados direto aqui para comparar, abrem processos
# por fork; com BLAS OpenMP isso travou uma rodada. Um núcleo só evita o travamento.
options(mc.cores = 1L)
suppressPackageStartupMessages({ library(lavaan); library(psych) })

resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome))
  resultados[[nome]] <<- ok
  invisible(ok)
}
t_inicio <- Sys.time()
quieto <- function(expr) suppressWarnings(suppressMessages(expr))

cat("== CVC de Hernández-Nieto (curso, página 'Coeficiente de Validade de Conteúdo (CVC)')\n")
# formato longo: cada juiz aparece uma vez por critério; item 1 da planilha da aula
j3 <- data.frame(criterio = rep(c("clareza", "pertinencia", "relevancia"), each = 3),
                 i1 = c(4, 4, 4, 5, 4, 5, 5, 5, 4), i2 = c(5, 5, 5, 5, 5, 5, 5, 5, 5), i3 = c(3, 4, 2, 2, 3, 3, 4, 4, 4))
cvc <- validade_conteudo_cvc(j3, c("i1", "i2", "i3"), criterio = "criterio")
lin <- function(cr, it) cvc$tabela[cvc$tabela$criterio == cr & cvc$tabela$item == it, ]
teste("notas (5, 5, 4): média 4,67, CVCi .933, CVCc .896 (tol. .0005)",
      abs(lin("relevancia", "i1")$cvc_i - .9333) < .0005 && abs(lin("relevancia", "i1")$cvc_c - .8963) < .0005)
teste("item 1 da aula: clareza (4, 4, 4) = .763, inadequada; pertinência (5, 4, 5) e relevância .896, adequadas", {
  abs(lin("clareza", "i1")$cvc_c - .763) < .0005 && !lin("clareza", "i1")$adequado &&
    lin("pertinencia", "i1")$adequado && lin("relevancia", "i1")$adequado })
teste("Pe por número de juízes: 3 -> .0370, 4 -> .0039 (a aula transcreve .039 por engano), 2 -> .25", {
  pe <- function(J) validade_conteudo_cvc(data.frame(i = rep(5, J)), "i")$tabela$pe
  abs(pe(3) - .0370) < 5e-5 && abs(pe(4) - .0039) < 5e-5 && abs(pe(2) - .25) < 1e-12 })
teste("2 juízes com nota máxima: CVCc = .75 (inadequado) e aviso de mínimo de 3 juízes", {
  r <- validade_conteudo_cvc(data.frame(i = c(5, 5)), "i")
  abs(r$tabela$cvc_c - .75) < 1e-12 && !r$tabela$adequado && length(r$avisos) == 1 })
teste("CVC total = média dos CVCi menos a média dos Pe (= média dos CVCc), por critério", {
  t <- cvc$tabela
  all(abs(cvc$efeito$valor - tapply(t$cvc_c, t$criterio, mean)[c("clareza", "pertinencia", "relevancia")]) < 1e-12) })
teste("piloto com 30 participantes: Pe desprezível, CVC = média / 5 (curso, página do CVC com população-alvo)", {
  set.seed(1); p30 <- data.frame(q1 = sample(4:5, 30, TRUE), q2 = sample(3:5, 30, TRUE))
  r <- validade_conteudo_cvc(p30, c("q1", "q2"))
  all(abs(r$tabela$cvc_c - colMeans(p30) / 5) < 1e-12) })
teste("item adequado só se passa nos três critérios (lista de itens adequados)", identical(cvc$itens_adequados, "i2"))

cat("== IVC (curso, página 'Índice de Validade de Conteúdo (IVC)'; Polit & Beck, 2006)\n")
teste("(3, 4, 2) -> I-IVC .67, inaceitável", {
  r <- validade_conteudo_ivc(data.frame(i = c(3, 4, 2)), "i"); abs(r$tabela$i_ivc - .667) < .0005 && !r$tabela$adequado_curso })
iv <- data.frame(i1 = c(3, 4, 2, 3, 4), i2 = c(4, 2, 3, 1, 2), i3 = c(4, 4, 4, 3, 4))
ivc <- validade_conteudo_ivc(iv, c("i1", "i2", "i3"))
teste("(3, 4, 2, 3, 4) -> .80, aceitável no curso e não em Polit e Beck (5 juízes pedem 1,00)",
      ivc$tabela$i_ivc[1] == .8 && ivc$tabela$adequado_curso[1] && !ivc$tabela$adequado_polit_beck[1])
teste("(4, 2, 3, 1, 2) -> .40, inaceitável", ivc$tabela$i_ivc[2] == .4 && !ivc$tabela$adequado_curso[2])
teste("S-IVC/Ave = média dos I-IVC (.733) e S-IVC/UA = 1/3 (tol. .0005)",
      abs(ivc$efeito$valor[1] - mean(c(.8, .4, 1))) < 1e-12 && abs(ivc$efeito$valor[2] - 1 / 3) < 1e-12)
teste("kappa modificado com 4 de 5: Pc = .15625, k* = .763 (conta à mão; tol. .0005)",
      abs(ivc$tabela$kappa_mod[1] - (.8 - 5 / 32) / (1 - 5 / 32)) < 1e-12 && abs(ivc$tabela$kappa_mod[1] - .763) < .0005)
teste("6 juízes, 5 favoráveis: .83 passa nos dois critérios (Polit e Beck: >= .78 com 6 ou mais)", {
  r <- validade_conteudo_ivc(data.frame(i = c(4, 4, 3, 4, 3, 2)), "i"); r$tabela$adequado_curso && r$tabela$adequado_polit_beck })

cat("== RVC de Lawshe (curso, página 'Razão de Validade de Conteúdo (RVC)'; Ayre & Scally, 2014)\n")
cr <- rvc_critico(5:10)
teste("valores críticos da aula para 5 a 10 juízes: Ne 5, 6, 7, 7, 8, 9", all(cr$ne_critico == c(5, 6, 7, 7, 8, 9)))
teste("RVC crítico da aula: 1,00; 1,00; 1,00; .75; .78; .80 (tol. .005)",
      all(abs(cr$rvc_critico - c(1, 1, 1, .75, .78, .80)) < .005))
teste("método binomial exato de 5 a 40 juízes: P(X >= Ne) <= .05 < P(X >= Ne - 1)", {
  t <- rvc_critico(5:40)
  all(pbinom(t$ne_critico - 1, t$n_juizes, .5, lower.tail = FALSE) <= .05) &&
    all(pbinom(t$ne_critico - 2, t$n_juizes, .5, lower.tail = FALSE) > .05) })
rv <- data.frame(a = c(rep(3, 9), 1), b = c(rep(3, 8), 2, 2), c = c(3, 3, rep(1, 8)))
rvc <- validade_conteudo_rvc(rv, c("a", "b", "c"))
teste("exemplos da aula com 10 juízes: Ne 9 -> .80 atende; Ne 8 -> .60 não; Ne 2 -> -.60",
      all(abs(rvc$tabela$rvc - c(.8, .6, -.6)) < 1e-12) && identical(rvc$tabela$atende, c(TRUE, FALSE, FALSE)))
teste("4 juízes: sem valor crítico (NA) e aviso de mínimo de 5", {
  r <- validade_conteudo_rvc(data.frame(i = c(3, 3, 3, 3)), "i"); is.na(r$tabela$rvc_critico) && any(grepl("mínimo é 5", r$avisos)) })

cat("== V de Aiken (Aiken, 1985; IC de Penfield & Giacobbi, 2004)\n")
av <- validade_conteudo_aiken(data.frame(x = c(5, 5, 4, 4, 5), y = c(3, 4, 2, 5, 3)), c("x", "y"))
teste("(5, 5, 4, 4, 5) em 1 a 5 -> V = 18 / 20 = .90", abs(av$tabela$v[1] - .9) < 1e-12)
teste("IC escore igual ao intervalo de Wilson do prop.test (equivalência; tol. 1e-8)", {
  w <- prop.test(18, 20, correct = FALSE)$conf.int
  abs(av$tabela$ic_inf[1] - w[1]) < 1e-8 && abs(av$tabela$ic_sup[1] - w[2]) < 1e-8 })

cat("== AFE: estrutura ordinal conhecida (3 fatores, 4 itens cada, carga .70, r = .30, 4 categorias)\n")
it <- paste0("y", 1:12)
lim <- function(vars) paste(sprintf("%s | -1*t1 + 0*t2 + 1*t3", vars), collapse = "\n")
pop <- paste("f1 =~ .7*y1 + .7*y2 + .7*y3 + .7*y4", "f2 =~ .7*y5 + .7*y6 + .7*y7 + .7*y8",
             "f3 =~ .7*y9 + .7*y10 + .7*y11 + .7*y12", "f1 ~~ .3*f2 + .3*f3", "f2 ~~ .3*f3", lim(it), sep = "\n")
d3 <- simulateData(pop, sample.nobs = 600, seed = 2026, ordered_center = FALSE)
adq <- afe_adequacao(d3, it)
R3 <- adq$correlacao
teste("itens de 4 categorias tratados como ordinais: correlação policórica", adq$tipo_correlacao == "policórica" && adq$ordinal)
teste("KMO geral e por item iguais ao psych::KMO na mesma matriz (tol. 1e-10)", {
  k <- psych::KMO(R3); abs(adq$efeito$valor - k$MSA) < 1e-10 && all(abs(adq$tabela$msa - k$MSAi) < 1e-10) })
teste("Bartlett igual ao psych::cortest.bartlett (χ², gl e p; tol. 1e-6)", {
  b <- psych::cortest.bartlett(R3, n = 600)
  abs(adq$pressupostos$estatistica[2] - b$chisq) < 1e-6 && adq$pressupostos$gl[2] == b$df && abs(adq$pressupostos$p[2] - b$p.value) < 1e-10 })
teste("Bartlett pela fórmula: -(n - 1 - (2p + 5)/6)·ln|R| com gl = p(p - 1)/2 (tol. 1e-8)",
      abs(adq$pressupostos$estatistica[2] - (-(600 - 1 - (2 * 12 + 5) / 6) * log(det(R3)))) < 1e-8 && adq$pressupostos$gl[2] == 66)
nf3 <- afe_numero_fatores(d3, it, n_iter = 20, ega = TRUE)
teste("análise paralela por permutação (policórica) recupera 3 fatores nas duas bases", nf3$n_sugerido == 3 && nf3$n_outra_base == 3)
teste("EGA (EGAnet) recupera 3 dimensões", nf3$n_ega == 3)
teste("semente fixa: análise paralela reproduzível e sem mexer no gerador do usuário", {
  set.seed(99); antes <- runif(1); set.seed(99)
  a <- afe_numero_fatores(d3[1:200, ], it, n_iter = 5); depois <- runif(1)
  b <- afe_numero_fatores(d3[1:200, ], it, n_iter = 5)
  identical(a$tabela, b$tabela) && antes == depois })
ef3 <- quieto(afe(d3, it, n_fatores = 3))
L <- ef3$cargas
alvo <- rep(1:3, each = 4)
mapa <- apply(abs(L), 1, which.max)
teste("minres + oblimin: cada bloco de 4 itens carrega no mesmo fator, e fatores diferentes entre blocos",
      all(tapply(mapa, alvo, function(x) length(unique(x))) == 1) && length(unique(tapply(mapa, alvo, `[`, 1))) == 3)
teste("cargas principais recuperadas (.70 ± .10) e cargas secundárias abaixo de .15",
      { sec <- abs(L); sec[cbind(1:12, mapa)] <- NA
        all(abs(abs(L[cbind(1:12, mapa)]) - .7) < .10) && max(sec, na.rm = TRUE) < .15 })
teste("correlações entre fatores recuperadas (.30 ± .10)",
      all(abs(ef3$correlacoes_fatores[upper.tri(ef3$correlacoes_fatores)] - .3) < .10))
teste("nenhum alerta de carga cruzada ou baixa na estrutura limpa", all(ef3$tabela$alerta == ""))
teste("variância por fator soma a comunalidade total (tol. 1e-10)",
      abs(sum(ef3$efeito$valor[1:3]) * 12 - sum(ef3$comunalidades)) < 1e-10)
ew <- quieto(afe(d3, it, n_fatores = 3, extracao = "wlsmv"))
teste("extração WLSMV (lavaan::efa) concorda com minres nas cargas (diferença máxima < .05) e traz ajuste",
      max(abs(ew$cargas - L)) < .05 && any(grepl("^CFI", ew$efeito$indice)) && ew$efeito$valor[grepl("^RMSEA", ew$efeito$indice)] < .06)
teste("número de fatores pela análise paralela quando n_fatores = NULL", {
  r <- quieto(afe(d3, it, n_iter = 20)); ncol(r$cargas) == 3 && any(grepl("análise paralela", r$apa)) })

cat("== AFE: estrutura contínua (psych::sim.structure) e equivalência com psych::fa.parallel\n")
set.seed(2026)
fx <- matrix(0, 12, 3); fx[1:4, 1] <- .7; fx[5:8, 2] <- .7; fx[9:12, 3] <- .7
Phi <- matrix(.3, 3, 3); diag(Phi) <- 1
dc <- as.data.frame(psych::sim.structure(fx, Phi, n = 500)$observed); names(dc) <- paste0("v", 1:12)
nfc <- afe_numero_fatores(dc, names(dc), n_iter = 20)
teste("itens contínuos: correlação de Pearson e 3 fatores", nfc$tipo_correlacao == "Pearson" && nfc$n_sugerido == 3)
teste("mesma decisão do psych::fa.parallel (componentes, Pearson)", {
  set.seed(2026); invisible(utils::capture.output(pp <- quieto(psych::fa.parallel(dc, fa = "pc", n.iter = 20, plot = FALSE))))
  pp$ncomp == nfc$n_sugerido })
efc <- quieto(afe(dc, names(dc), n_fatores = 3))
teste("cargas recuperadas em dado contínuo (.70 ± .10)", all(abs(apply(abs(efc$cargas), 1, max) - .7) < .10))

cat("== AFE: alertas por item (carga cruzada, carga baixa, MSA baixo) e escala 0-100\n")
popx <- paste("f1 =~ .7*a1 + .7*a2 + .7*a3 + .7*a4 + .7*a5 + .5*x1 + .15*w1",
              "f2 =~ .7*b1 + .7*b2 + .7*b3 + .7*b4 + .7*b5 + .5*x1", "f1 ~~ .3*f2", "z ~~ 1*z", sep = "\n")
itx <- c(paste0("a", 1:5), paste0("b", 1:5), "x1", "w1", "z")
dx <- quieto(simulateData(popx, sample.nobs = 500, seed = 7))
dx <- as.data.frame(lapply(dx[itx], function(v) round(50 + 15 * v, 1)))
efx <- quieto(afe(dx, itx, n_iter = 20))
al <- setNames(efx$tabela$alerta, itx)
teste("escala 0-100 tratada como contínua (Pearson), como no PBAT", efx$tipo_correlacao == "Pearson")
teste("2 fatores pela análise paralela", ncol(efx$cargas) == 2)
teste("item com .50 nos dois fatores marcado como carga cruzada", al[["x1"]] == "carga cruzada")
teste("item com carga .15 e item de ruído marcados como carga baixa", al[["w1"]] == "carga baixa" && al[["z"]] == "carga baixa")
teste("item de ruído com o menor MSA, abaixo de .50, e aviso", {
  m <- efx$msa; m$item[which.min(m$msa)] == "z" && min(m$msa) < .50 && any(grepl("MSA abaixo de .50", efx$avisos)) })

cat("== Privacidade, ausentes e relato\n")
d3na <- d3; d3na$y1[c(3, 10, 50)] <- NA
ena <- quieto(afe(d3na, it, n_fatores = 3))
teste("deleção listwise contada em n_excluidos (3 casos)", ena$n == 597 && ena$n_excluidos == 3)
teste("nenhuma linha de participante no resultado: tabelas só com linhas por item, fator ou índice", {
  dfs <- Filter(function(x) is.data.frame(x) || is.matrix(x), ef3[setdiff(names(ef3), "modelo")])
  all(vapply(dfs, NROW, 0) <= 12) && length(utils::capture.output(print(ef3))) < 80 })
teste("print de todas as funções não falha e não passa de 80 linhas", {
  all(vapply(list(cvc, ivc, rvc, av, adq, nf3, ef3), function(r) { n <- length(utils::capture.output(print(r))); n > 5 && n < 80 }, TRUE)) })
teste("frase APA da AFE com KMO, Bartlett, extração, rotação, variância e correlações entre fatores", {
  a <- paste(ef3$apa, collapse = " ")
  grepl("KMO = .", a, fixed = TRUE) && grepl("χ²(66)", a, fixed = TRUE) && grepl("rotação oblimin", a) &&
    grepl("explicou", a) && grepl("correlações entre fatores variaram", a) })
teste("frase APA em inglês com vírgula decimal", {
  a <- quieto(afe(d3, it, n_fatores = 3, decimal = ",", idioma = "en"))$apa
  any(grepl("exploratory factor analysis", a)) && any(grepl("KMO = ,", a, fixed = TRUE)) })
teste("um fator: sem rotação e sem correlação entre fatores", {
  pop1 <- paste("f =~ .6*q1 + .6*q2 + .6*q3 + .6*q4 + .6*q5 + .6*q6",
                paste(sprintf("q%d | -0.5*t1 + 0.5*t2", 1:6), collapse = "\n"), sep = "\n")
  d1 <- simulateData(pop1, sample.nobs = 400, seed = 11, ordered_center = FALSE)
  n1 <- afe_numero_fatores(d1, paste0("q", 1:6), n_iter = 20)
  r1 <- quieto(afe(d1, paste0("q", 1:6), n_fatores = 1))
  n1$n_sugerido == 1 && ncol(r1$cargas) == 1 && !any(grepl("^r\\(", r1$efeito$indice)) && all(abs(r1$cargas - .6) < .12) })

cat("== Correções da revisão v1.1.0: dependência linear exata (C7)\n")
# Itens A1-A5 do psych::bfi (dados do pacote, não do usuário) com o escore total junto:
# posto 5 de 6 (qr da base R) e a função tem de parar nomeando a coluna, em vez de
# "system is computationally singular" (antes) ou de suavizar a matriz.
bfiA <- stats::na.omit(psych::bfi[paste0("A", 1:5)])
erro_de <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
teste("C7 escore total junto dos itens: afe_adequacao e afe param com mensagem que nomeia 'total' (posto 5 de 6 pelo qr)", {
  x <- bfiA; x$total <- rowSums(x)
  m1 <- erro_de(afe_adequacao(x, names(x), ordinal = FALSE)); m2 <- erro_de(quieto(afe(x, names(x), n_fatores = 1, ordinal = FALSE)))
  qr(scale(as.matrix(x)))$rank == 5 && grepl("linearmente dependentes", m1) && grepl("total é combinação exata", m1) &&
    grepl("escore total", m1) && grepl("linearmente dependentes", m2) && !grepl("computationally singular", m1) })
teste("C7 item duplicado e item invertido (7 - A2, ordinal): param na checagem dos dados brutos, antes da policórica", {
  x <- bfiA; x$A1_copia <- x$A1; y <- bfiA; y$A2_inv <- 7 - y$A2
  m1 <- erro_de(afe_adequacao(x, names(x))); m2 <- erro_de(afe_numero_fatores(y, names(y), n_iter = 5))
  grepl("A1_copia", m1) && grepl("linearmente dependentes", m1) && grepl("A2_inv", m2) && grepl("versão invertida", m2) })
teste("C7 n menor ou igual ao número de itens: mensagem própria, sem atribuir a dependência entre itens", {
  m <- erro_de(afe_adequacao(bfiA[1:5, ], paste0("A", 1:5), ordinal = FALSE))
  grepl("5 participantes com dados completos para 5 itens", m) && !grepl("linearmente dependentes", m) })
# Policórica não positiva definida sem dependência exata (n = 40, 8 itens de 4 categorias
# com categoria extrema rara; semente 4): a suavização continua valendo, com aviso, e o KMO
# é o do psych::KMO sobre a policórica do próprio psych (que suaviza do mesmo jeito).
teste("C7 policórica não positiva definida sem dependência exata: suavizada com aviso; KMO igual ao do psych (tol. 1e-10)", {
  set.seed(4); f <- rnorm(40); X <- sapply(1:8, function(j) .8 * f + .6 * rnorm(40))
  D <- as.data.frame(apply(X, 2, function(x) cut(x, c(-Inf, -1.5, 0, 1.6, Inf), labels = FALSE)))
  bruta <- quieto(psych::polychoric(D, smooth = FALSE)$rho)
  a <- quieto(afe_adequacao(D, names(D)))
  min(eigen(bruta)$values) < 0 && min(eigen(stats::cor(D))$values) > .05 &&
    any(grepl("suavizada com psych::cor.smooth", a$avisos)) &&
    abs(a$kmo$geral - psych::KMO(quieto(psych::polychoric(D)$rho))$MSA) < 1e-10 })

cat("== Correções da revisão v1.1.0: frase de fatorabilidade e caso Heywood (C12)\n")
# Seis normais independentes, n = 60 (semente 7; reprodução do revisor): KMO .45 e
# Bartlett χ²(15) = 14.94, p = .455, conferidos no psych::KMO e no psych::cortest.bartlett.
teste("C12 KMO .45 e Bartlett p = .455: a frase diz 'não foi indicada' com os números (psych confere KMO e χ²)", {
  set.seed(7); x <- as.data.frame(matrix(rnorm(60 * 6), 60, 6)); x[] <- lapply(x, function(v) round(v * 2 + 10, 1))
  a <- afe_adequacao(x, names(x)); k <- psych::KMO(stats::cor(x)); b <- psych::cortest.bartlett(stats::cor(x), n = 60)
  abs(a$kmo$geral - k$MSA) < 1e-10 && abs(b$chisq - 14.94) < .005 && abs(b$p.value - .455) < .0005 && !a$fatoravel &&
    grepl("não foi indicada: KMO = .45 (inaceitável", a$apa, fixed = TRUE) && grepl("χ²(15) = 14.94, p = .455", a$apa, fixed = TRUE) &&
    grepl("não rejeitou a matriz identidade", a$apa) && !grepl("foi indicada pelo", a$apa) })
# Três itens com r12 = r13 = .80 e r23 = .50 (dados com essa matriz exata, n = 300): KMO .53
# (psych::KMO), faixa miserável, com Bartlett significativo. Na solução de um fator com três
# indicadores, λ1² = r12·r13 / r23 = 1,28 (conta à mão): a carga de h1 passa de 1 (Heywood).
S3 <- matrix(c(1, .8, .8, .8, 1, .5, .8, .5, 1), 3)
set.seed(1); h3 <- as.data.frame(MASS::mvrnorm(300, rep(0, 3), S3, empirical = TRUE)); names(h3) <- paste0("h", 1:3)
teste("C12 KMO .53 (miserável) com Bartlett significativo: frase de não indicação e aviso da faixa miserável", {
  a <- afe_adequacao(h3, names(h3)); k <- psych::KMO(S3)
  abs(a$kmo$geral - k$MSA) < 1e-10 && round(k$MSA, 2) == .53 && !a$fatoravel &&
    grepl("não foi indicada: KMO = .53 (miserável", a$apa, fixed = TRUE) && grepl("que rejeitou a matriz identidade", a$apa) &&
    any(grepl("faixa miserável", a$avisos)) &&
    grepl("not supported: KMO = .53 (miserable", afe_adequacao(h3, names(h3), idioma = "en")$apa, fixed = TRUE) })
teste("C12 KMO >= .60 e Bartlett significativo (estrutura simulada, KMO .83): frase 'foi indicada'",
      adq$fatoravel && grepl("foi indicada pelo KMO = .83", adq$apa, fixed = TRUE))
teste("C12 caso Heywood (λ1 de h1 acima de 1 pela conta à mão): alerta, aviso e frase de solução imprópria em pt e en", {
  e <- quieto(afe(h3, names(h3), n_fatores = 1)); en <- quieto(afe(h3, names(h3), n_fatores = 1, idioma = "en"))
  fr <- e$apa[length(e$apa)]
  sqrt(.8 * .8 / .5) > 1 && e$solucao_impropria && grepl("^Heywood", e$tabela$alerta[1]) && abs(e$cargas[1, 1]) > 1 &&
    any(grepl("Caso Heywood", e$avisos)) && grepl("resultou em solução imprópria, com caso Heywood", fr) && grepl("em h1", fr) &&
    !grepl("sem cargas cruzadas", fr) && grepl("yielded an improper solution", en$apa[length(en$apa)]) })
teste("C12 solução sem Heywood continua com a frase descritiva", !ef3$solucao_impropria && grepl("sem cargas cruzadas", ef3$apa[length(ef3$apa)]))

cat("== Correções da revisão v1.1.0 com o banco G Factor do JASP (C7 e C12)\n")
DJ <- Sys.getenv("PR_DADOS_JASP")
jasp <- function(arquivo) {
  if (!nzchar(DJ)) return(NULL)
  f <- list.files(DJ, recursive = TRUE, full.names = TRUE)
  f <- f[basename(f) == arquivo]
  if (!length(f)) return(NULL)
  utils::read.csv(f[1], check.names = FALSE)
}
gf <- jasp("G Factor.csv")
if (is.null(gf)) cat("PULADO  testes com G Factor (PR_DADOS_JASP vazia ou banco não encontrado)\n") else {
  # Age = Years + Months/12 e cada Residuals.X é o resíduo de X sobre Age: R² = 1 nas
  # regressões da base R (conferência independente da dependência).
  teste("C7 G Factor com todas as colunas numéricas: para e nomeia Age e as colunas Residuals.* (R² = 1 no lm)", {
    it <- names(gf)[vapply(gf, is.numeric, TRUE)]
    m <- erro_de(afe_adequacao(gf, it))
    r2 <- c(summary(stats::lm(Age ~ Years + Months, gf))$r.squared,
            vapply(grep("^Residuals", names(gf), value = TRUE), function(v)
              quieto(summary(stats::lm(gf[[v]] ~ gf[[sub("Residuals.", "", v, fixed = TRUE)]] + gf$Age))$r.squared), 0))
    all(abs(r2 - 1) < 1e-8) && grepl("linearmente dependentes", m) && grepl("Age", m) &&
      all(vapply(grep("^Residuals", names(gf), value = TRUE), grepl, TRUE, x = m, fixed = TRUE)) && !grepl("mediano", m) })
  # 7 medidas, 2 fatores, minres + oblimin: carga de Weight 1.001 no psych::fa (aviso de
  # ultra-Heywood do próprio psych).
  teste("C12 G Factor, 2 fatores: Heywood em Weight (carga 1.001 no psych::fa) e frase de solução imprópria", {
    it7 <- c("Pitch", "Light", "Weight", "Classics", "French", "English", "Mathematics")
    e <- quieto(afe(gf, it7, n_fatores = 2, ordinal = FALSE))
    f <- quieto(psych::fa(stats::cor(stats::na.omit(gf[it7])), 2, n.obs = nrow(stats::na.omit(gf[it7])), fm = "minres", rotate = "oblimin"))
    fr <- e$apa[length(e$apa)]
    max(abs(unclass(f$loadings)["Weight", ])) > 1 && e$solucao_impropria && grepl("^Heywood", e$tabela$alerta[e$tabela$item == "Weight"]) &&
      grepl("solução imprópria, com caso Heywood", fr) && grepl("em Weight", fr) && !grepl("sem cargas cruzadas", fr) })
}

cat("== v1.4: Hull e MAP no número de fatores (Lorenzo-Seva, Timmerman & Kiers, 2011; Velicer, 1976)\n")
sim_fat <- function(n, nf, ipf, semente, carga = c(.45, .75), phi = .3, ordinal = FALSE) {
  set.seed(semente); p <- nf * ipf; L <- matrix(0, p, nf)
  for (k in seq_len(nf)) L[((k - 1) * ipf + 1):(k * ipf), k] <- stats::runif(ipf, carga[1], carga[2])
  Phi <- matrix(phi, nf, nf); diag(Phi) <- 1; S <- L %*% Phi %*% t(L); diag(S) <- 1
  X <- MASS::mvrnorm(n, rep(0, p), S)
  if (ordinal) X <- apply(X, 2, function(x) as.integer(cut(x, c(-Inf, -1, -.3, .4, 1.2, Inf))))
  d <- as.data.frame(X); names(d) <- paste0("i", seq_len(p)); d
}
h3 <- quieto(afe_numero_fatores(sim_fat(500, 3, 6, 31), paste0("i", 1:18), n_iter = 20))
h4 <- quieto(afe_numero_fatores(sim_fat(600, 4, 6, 32), paste0("i", 1:24), n_iter = 20))
ho <- quieto(afe_numero_fatores(sim_fat(500, 3, 6, 33, ordinal = TRUE), paste0("i", 1:18), n_iter = 20))
teste("Hull-CFI recupera a estrutura simulada: 3 fatores (contínuo), 4 fatores (contínuo) e 3 fatores (ordinal, policórica)",
      h3$n_hull == 3 && h4$n_hull == 4 && ho$n_hull == 3)
teste("tabela do Hull: 0 a J fatores, J = paralela + 1, gl do modelo nulo p(p - 1)/2 e CFI = 0 na solução nula", {
  t <- h3$tabelas_extra$hull
  identical(t$fatores, 0:(h3$n_sugerido + 1)) && t$gl[1] == 18 * 17 / 2 && t$cfi[1] == 0 &&
    all(t$gl[-1] == ((18 - t$fatores[-1])^2 - (18 + t$fatores[-1])) / 2) })
teste("CFI de cada solução = conta à mão com stats::factanal (F = objetivo da ML, qui2 = (N - 1)F; tol. 1e-8)", {
  d <- sim_fat(500, 3, 6, 31); R <- stats::cor(d); n <- 500; p <- 18
  cf <- vapply(1:4, function(q) {
    F <- unname(stats::factanal(covmat = R, factors = q, n.obs = n)$criteria["objective"])
    gl <- ((p - q)^2 - (p + q)) / 2; q2 <- F * (n - 1); q0 <- -(n - 1) * log(det(R))
    1 - max(q2 - gl, 0) / max(max(q2 - gl, 0), q0 - p * (p - 1) / 2) }, 0)
  all(abs(cf - h3$tabelas_extra$hull$cfi[2:5]) < 1e-8) })
if (requireNamespace("EFAtools", quietly = TRUE)) {
  teste("Hull igual ao EFAtools::HULL (ML, CFI e CAF, mesmo J): CFI e st com tol. 1e-6, mesma escolha", {
    ok <- TRUE
    for (cs in list(list(h3, sim_fat(500, 3, 6, 31), FALSE), list(h4, sim_fat(600, 4, 6, 32), FALSE), list(ho, sim_fat(500, 3, 6, 33, ordinal = TRUE), TRUE))) {
      r <- cs[[1]]; d <- cs[[2]]; J <- max(r$tabelas_extra$hull$fatores)
      R <- if (cs[[3]]) quieto(psych::polychoric(d, smooth = FALSE))$rho else stats::cor(d)
      e <- quieto(EFAtools::HULL(R, N = nrow(d), method = "ML", gof = c("CFI", "CAF"), n_fac_theor = J - 1))
      g <- function(nm) e$results[[which(vapply(e$results, `[[`, "", "name") == nm)]]
      ok <- ok && length(g("CFI")$y) == J + 1 && max(abs(g("CFI")$y - r$tabelas_extra$hull$cfi)) < 1e-6 &&
        max(abs(g("CFI")$st - r$tabelas_extra$hull$st_cfi), na.rm = TRUE) < 1e-6 &&
        e$n_fac_CFI == r$n_hull && e$n_fac_CAF == r$n_hull_caf
    }
    ok })
} else cat("PULADO  comparação com EFAtools::HULL (pacote não instalado; o teste da conta à mão acima cobre o CFI)\n")
teste("MAP igual ao psych::vss (média dos quadrados das parciais, tol. 1e-12) e recupera 3 componentes", {
  R <- stats::cor(sim_fat(500, 3, 6, 31))
  v <- psych::vss(R, n = 17, n.obs = 500, plot = FALSE)$map
  max(abs(v - h3$tabelas_extra$map$map[-1])) < 1e-12 && h3$n_map == 3 })
teste("MAP com m = 0: média dos quadrados das correlações fora da diagonal (Velicer, 1976)", {
  R <- stats::cor(sim_fat(500, 3, 6, 31)); R0 <- R; diag(R0) <- 0
  abs(h3$tabelas_extra$map$map[1] - sum(R0^2) / (18 * 17)) < 1e-14 })
teste("Hull não roda com menos de 6 itens (aviso) e a paralela segue", {
  r <- quieto(afe_numero_fatores(sim_fat(300, 1, 5, 34), paste0("i", 1:5), n_iter = 10))
  is.na(r$n_hull) && any(grepl("Hull não calculado", r$avisos)) && r$n_sugerido == 1 })
teste("n_teorico maior que a paralela sobe o teto do Hull (J = n_teorico + 1) e valor inválido para", {
  r <- quieto(afe_numero_fatores(sim_fat(500, 3, 6, 31), paste0("i", 1:18), n_iter = 10, n_teorico = 5))
  max(r$tabelas_extra$hull$fatores) == 6 && grepl("inteiro positivo", erro_de(afe_numero_fatores(sim_fat(200, 3, 6, 31), paste0("i", 1:18), n_iter = 5, n_teorico = 2.5))) })
teste("frase APA cita paralela, Hull e MAP; discordância entre paralela e Hull vira aviso", {
  r <- quieto(afe_numero_fatores(sim_fat(250, 1, 10, 14), paste0("i", 1:10), n_iter = 20))
  grepl("Hull", h3$apa) && grepl("MAP", h3$apa) &&
    (r$n_hull == r$n_sugerido || any(grepl("Hull-CFI sugere", r$avisos))) })

cat("== v1.4: teste-reteste e concordância entre avaliadores\n")
# Shrout e Fleiss (1979), tabela 2: 6 alvos, 4 juízes; valores publicados na tabela 4.
sf <- data.frame(j1 = c(9, 6, 8, 7, 10, 6), j2 = c(2, 1, 4, 1, 5, 2), j3 = c(5, 3, 6, 2, 6, 4), j4 = c(8, 2, 8, 6, 9, 7))
ci <- quieto(concordancia_juizes(sf, names(sf), "intervalar"))
teste("Shrout & Fleiss (1979): ICC(1,1) .17, ICC(2,1) .29, ICC(3,1) .71, ICC(1,4) .44, ICC(2,4) .62, ICC(3,4) .91 (tol. .005)", {
  v <- ci$tabela$icc; names(v) <- ci$tabela$tipo
  all(abs(v[c("ICC(1)", "ICC(A,1)", "ICC(C,1)", "ICC(1,k)", "ICC(A,k)", "ICC(C,k)")] - c(.17, .29, .71, .44, .62, .91)) < .005) })
set.seed(41); vv <- stats::rnorm(60, 50, 10)
rt <- data.frame(t1 = vv + stats::rnorm(60, 0, 4), t2 = vv + 5 + stats::rnorm(60, 0, 4))
rr <- quieto(confiabilidade_teste_reteste(rt, "t1", "t2"))
teste("teste-reteste: ICC(A,1) = conta à mão pelos quadrados médios de McGraw e Wong (1996) (tol. 1e-10)", {
  x <- as.matrix(rt); n <- nrow(x); k <- 2; mg <- mean(x)
  msr <- k * sum((rowMeans(x) - mg)^2) / (n - 1); msc <- n * sum((colMeans(x) - mg)^2) / (k - 1)
  mse <- (sum((x - mg)^2) - (n - 1) * msr - (k - 1) * msc) / ((n - 1) * (k - 1))
  abs(rr$efeito$valor[1] - (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)) < 1e-10 })
teste("teste-reteste: EPM = DP(todas) x raiz(1 - ICC), MMD95 = 1,96 x raiz(2) x EPM, Bland-Altman = média ± 1,96 DP das diferenças", {
  e <- rr$tabelas_extra[[1]]; epm <- stats::sd(c(rt$t1, rt$t2)) * sqrt(1 - rr$efeito$valor[1]); df <- rt$t2 - rt$t1
  abs(e$epm - epm) < 1e-10 && abs(e$mmd95 - 1.96 * sqrt(2) * epm) < 1e-10 &&
    abs(e$ba_inf - (mean(df) - 1.96 * stats::sd(df))) < 1e-10 && abs(e$ba_sup - (mean(df) + 1.96 * stats::sd(df))) < 1e-10 })
teste("teste-reteste com mudança média (+5 pontos): aviso de mudança entre ocasiões; frase com ICC(A,1) e régua de Koo e Li", {
  any(grepl("A média mudou", rr$avisos)) && grepl("ICC\\(A,1\\) = ", rr$apa) && grepl("Koo e Li", rr$apa) })
# kappa de Cohen à mão: tabela 2 x 2 [[20, 5], [10, 15]]: po = .70, pe = .50, kappa = .40.
kk <- data.frame(a = rep(c(1, 1, 2, 2), c(20, 5, 10, 15)), b = rep(c(1, 2, 1, 2), c(20, 5, 10, 15)))
teste("kappa de Cohen = .40 na tabela [[20, 5], [10, 15]] (po .70, pe .50; tol. 1e-12)",
      abs(quieto(concordancia_juizes(kk, c("a", "b"), "nominal"))$efeito$valor[1] - .40) < 1e-12)
teste("kappa ponderado quadrático = conta à mão (1 - soma w.o / soma w.e) numa tabela 3 x 3", {
  t3 <- matrix(c(15, 4, 1, 3, 12, 5, 0, 4, 16), 3, byrow = TRUE)
  idx <- which(t3 > 0, arr.ind = TRUE); d3 <- data.frame(a = rep(idx[, 1], t3[idx]), b = rep(idx[, 2], t3[idx]))
  o <- t3 / sum(t3); e <- outer(rowSums(o), colSums(o)); w <- outer(1:3, 1:3, function(i, j) (i - j)^2)
  abs(quieto(concordancia_juizes(d3, c("a", "b"), "ordinal"))$efeito$valor[1] - (1 - sum(w * o) / sum(w * e))) < 1e-10 })
# Fleiss: exemplo resolvido do verbete "Fleiss' kappa" (10 casos, 14 avaliadores, 5 categorias): kappa = .210.
cont <- matrix(c(0, 0, 0, 0, 14, 0, 2, 6, 4, 2, 0, 0, 3, 5, 6, 0, 3, 9, 2, 0, 2, 2, 8, 1, 1,
                 7, 7, 0, 0, 0, 3, 2, 6, 3, 0, 2, 5, 3, 2, 2, 6, 5, 2, 1, 0, 0, 2, 2, 3, 7), 10, byrow = TRUE)
fl <- as.data.frame(t(apply(cont, 1, function(l) rep(1:5, l)))); names(fl) <- paste0("r", 1:14)
kf <- quieto(concordancia_juizes(fl, names(fl), "nominal", reamostras = 200))
teste("kappa de Fleiss = .210 no exemplo resolvido (10 casos, 14 avaliadores; tol. .0005) e IC bootstrap contém a estimativa",
      abs(kf$efeito$valor[1] - .210) < .0005 && kf$efeito$ic_inf[1] < kf$efeito$valor[1] && kf$efeito$ic_sup[1] > kf$efeito$valor[1])
if (requireNamespace("irr", quietly = TRUE)) {
  teste("irr confere: kappa2, kappa2 quadrático, kappam.fleiss e icc (A,1) e (A,k) (tol. 1e-10)", {
    a <- abs(quieto(concordancia_juizes(kk, c("a", "b"), "nominal"))$efeito$valor[1] - irr::kappa2(kk)$value) < 1e-10
    b <- abs(kf$efeito$valor[1] - irr::kappam.fleiss(fl)$value) < 1e-10
    c <- abs(ci$efeito$valor[1] - irr::icc(sf, "twoway", "agreement", "single")$value) < 1e-10 &&
      abs(ci$efeito$valor[2] - irr::icc(sf, "twoway", "agreement", "average")$value) < 1e-10
    a && b && c })
} else cat("PULADO  conferência com o pacote irr (não instalado; as contas à mão acima cobrem os mesmos índices)\n")
teste("Fleiss com tipo ordinal avisa que trata categorias como nominais", {
  any(grepl("nominais", quieto(concordancia_juizes(fl[, 1:5], names(fl)[1:5], "ordinal", reamostras = 50))$avisos)) })

cat("== v1.4: análise clássica de itens e confiabilidade em EMA\n")
set.seed(2); ita_th <- stats::rnorm(500)
ita <- as.data.frame(sapply(1:8, function(j) as.integer(cut(.8 * ita_th + stats::rnorm(500, 0, .7) + c(0, .2, -.3, .5, 0, -.6, .3, .1)[j],
                                                              c(-Inf, -1, -.2, .6, 1.4, Inf))))); names(ita) <- paste0("i", 1:8)
ai <- analise_classica_itens(ita, names(ita))
teste("itens: RIT, RIR e alfa sem o item iguais ao psych::alpha (raw.r, r.drop, alpha.drop; tol. 1e-12)", {
  a <- quieto(psych::alpha(ita, warnings = FALSE))
  max(abs(a$item.stats$raw.r - ai$tabela$rit)) < 1e-12 && max(abs(a$item.stats$r.drop - ai$tabela$rir)) < 1e-12 &&
    max(abs(a$alpha.drop$raw_alpha - ai$tabela$alfa_sem_item)) < 1e-12 })
teste("itens: dificuldade corrigida = (média - mín)/(máx - mín) e ULI = diferença entre o terço superior e o inferior do total", {
  tot <- rowSums(ita); g <- cut(tot, unique(stats::quantile(tot, c(0, 1/3, 2/3, 1))), include.lowest = TRUE, labels = FALSE)
  u1 <- (mean(ita$i1[g == 3]) - mean(ita$i1[g == 1])) / (max(ita$i1) - min(ita$i1))
  abs(ai$tabela$dificuldade[1] - (mean(ita$i1) - 1) / 4) < 1e-12 && abs(ai$tabela$uli[1] - u1) < 1e-12 })
teste("itens: item invertido sem recodificar gera aviso de RIR negativo", {
  i2 <- ita; i2$i3 <- 6 - i2$i3; any(grepl("RIR negativo", analise_classica_itens(i2, names(i2))$avisos)) })
set.seed(1); ema <- expand.grid(id = 1:40, time = 1:10)
ema_p <- stats::rnorm(40)[ema$id]; ema_t <- stats::rnorm(10, 0, .3)[ema$time]; ema_pt <- stats::rnorm(nrow(ema), 0, .7)
for (i in 1:4) ema[[paste0("i", i)]] <- ema_p + ema_t + ema_pt + stats::rnorm(nrow(ema), 0, .5) + stats::rnorm(1, 0, .2)
ce <- quieto(confiabilidade_ema(ema, "id", "time", paste0("i", 1:4)))
teste("EMA: coeficientes iguais ao psych::multilevel.reliability (tol. 1e-12)", {
  r <- quieto(psych::multilevel.reliability(ema, grp = "id", Time = "time", items = paste0("i", 1:4), lme = FALSE, plot = FALSE, alpha = FALSE))
  max(abs(ce$tabela$valor - c(r$RkF, r$R1R, r$RkR, r$Rc))) < 1e-12 })
teste("EMA: Rc = σ²(pessoa × tempo) / [σ²(pessoa × tempo) + σ²(resíduo)/itens] pelos componentes; perto do populacional .887 (tol. .05)", {
  cp <- ce$tabelas_extra$`Componentes de variância`; v <- stats::setNames(cp$variancia, cp$componente)
  abs(ce$tabela$valor[4] - v[["ID x time"]] / (v[["ID x time"]] + v[["Residual"]] / 4)) < 1e-6 && abs(ce$tabela$valor[4] - .887) < .05 })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
if (!nzchar(DJ)) cat("Aviso: PR_DADOS_JASP vazia; os testes com o banco G Factor do JASP não rodaram.\n")
print(versoes_pacotes(c("psych", "GPArotation", "lavaan", "EGAnet", "EFAtools", "irr")))
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
