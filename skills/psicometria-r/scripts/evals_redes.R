# evals_redes.R, testes do módulo de redes psicométricas (est_redes.R)
#
# Rodar da pasta da skill:
#   PR_DADOS_JASP="/caminho/para/Data Library" Rscript scripts/evals_redes.R
# PR_DADOS_JASP aponta para a biblioteca de dados do JASP (a pasta que contém
# "12. Network Analysis"; qualquer pasta acima dela também serve, a busca é recursiva).
# Sem a variável, os testes de âncora do livro são pulados com aviso.
# Três tipos de teste, nesta ordem: (1) âncoras JASP 0.95.1 do livro do usuário (Vol. 7,
# banco BFI Network, com a linha md7:linha) e da tabela do JASP transcrita nas notas do
# curso (spec_T6, P10), que estimou a mesma rede; (2) verdade conhecida em dados
# simulados com semente fixa; (3) equivalência com bootnet, networktools e mgm chamados
# diretamente. Tolerância declarada em cada teste: para âncora, metade da última casa
# impressa pelo JASP (.0005).
# Reamostras: para caber em ~2 minutos, os bootstraps daqui usam de 20 a 120 reamostras
# (declarado em cada teste). Para relato, a função usa 1000 (padrão do livro).
# Saída: uma linha PASSOU/FALHOU por teste e o resumo; código de saída 1 se algo falhar.

args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
source(file.path(pasta, "carregar.R"))
if (!exists("rede_psicometrica")) sys.source(file.path(pasta, "est_redes.R"), envir = globalenv())

resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome))
  resultados[[nome]] <<- ok
  invisible(ok)
}
perto <- function(a, b, tol) isTRUE(all(abs(a - b) <= tol + 1e-12))
t_inicio <- Sys.time()

DJ <- Sys.getenv("PR_DADOS_JASP")
jasp <- function(arquivo) {
  if (!nzchar(DJ)) return(NULL)
  f <- list.files(DJ, recursive = TRUE, full.names = TRUE)
  f <- f[basename(f) == arquivo]
  if (!length(f)) return(NULL)
  utils::read.csv(f[1], check.names = FALSE)
}
pular <- function(banco) cat("PULADO  âncoras de", banco, "(PR_DADOS_JASP vazia ou banco não encontrado)\n")
todos <- list()  # resultados guardados para os testes de formatação

# ============================================================ 1. Âncoras JASP (livro, Vol. 7)
bfi <- jasp("BFI Network.csv")
if (is.null(bfi)) pular("BFI Network") else {
  cat("== Âncoras: BFI Network (Vol. 7, caps. 23 e 24; JASP 0.95.1, padrões do módulo Network)\n")
  it <- c(paste0("A", 1:5), paste0("C", 1:5), paste0("E", 1:5), paste0("N", 1:5), paste0("O", 1:5))
  rb <- rede_psicometrica(bfi, it)
  todos$bfi <- rb
  val <- function(ind) rb$tabela$valor[rb$tabela$indicador == ind]
  W <- rb$pesos
  teste("25 nós e 300 arestas possíveis, p(p-1)/2", val("nós") == 25 && val("arestas possíveis") == 300)
  teste("181 arestas não nulas de 300, 53 negativas, esparsidade .397 (JASP; contagem da Weights matrix em spec_T6 P10; o livro diz \"cerca de metade\")",
        val("arestas não nulas") == 181 && val("arestas negativas") == 53 && perto(val("esparsidade"), .397, .0005))
  teste("N = 2800 com ausentes par a par e tamanho de amostra \"Maximum\"", val("N usado na estimação") == 2800)
  teste("N1–N2 = .570, a aresta mais forte",
        perto(W["N1", "N2"], .570, .0005) && rb$tabelas_extra[[1]]$aresta[1] == "N1–N2")
  teste("N3–N4 = .302, a segunda maior do bloco", perto(W["N3", "N4"], .302, .0005))
  teste("C4–C5 = .301, A2–A3 = .274, A1–A2 = -.273, E2–E4 = -.274 (Weights matrix do JASP, spec_T6 P10; tol .0005)",
        perto(c(W["C4", "C5"], W["A2", "A3"], W["A1", "A2"], W["E2", "E4"]), c(.301, .274, -.273, -.274), .0005))
  fz <- rb$centralidade_z$forca; names(fz) <- it
  teste("força em z: C4 1.649, E2 1.389, E4 1.302, E5 1.289, C3 -1.706, A4 -1.502",
        perto(fz[c("C4", "E2", "E4", "E5", "C3", "A4")], c(1.649, 1.389, 1.302, 1.289, -1.706, -1.502), .0005) &&
          names(which.max(fz)) == "C4" && names(which.min(fz)) == "C3")
  # Tabela "Centrality measures per variable" do JASP (Betweenness, Closeness, Strength,
  # Expected influence), transcrita em spec_T6 P10; a coluna Strength coincide com.
  jt <- rbind(
    A1 = c(-1.096, -1.141, -1.058, -1.651), A2 = c(0.596, 0.143, 0.689, 0.152), A3 = c(1.020, 0.523, -0.359, 0.824),
    A4 = c(0.385, -0.063, -1.502, -0.135), A5 = c(-0.038, 0.739, -0.748, 0.788), C1 = c(-1.307, -1.178, -1.323, -0.240),
    C2 = c(1.548, 0.684, 1.222, 0.890), C3 = c(-1.413, -1.150, -1.706, -0.804), C4 = c(0.385, 0.086, 1.649, -0.599),
    C5 = c(0.596, 0.662, 0.149, -0.910), E1 = c(-1.413, -0.474, -0.203, -1.776), E2 = c(0.279, 1.108, 1.389, -1.216),
    E3 = c(0.491, 0.212, 0.190, 1.117), E4 = c(1.337, 1.449, 1.302, -0.971), E5 = c(-0.250, 1.215, 1.289, -0.378),
    N1 = c(-0.355, -1.767, 0.650, 1.723), N2 = c(-1.413, -2.340, 0.146, 1.281), N3 = c(1.125, -0.724, -0.710, 1.478),
    N4 = c(1.866, 0.546, 0.269, 1.015), N5 = c(0.279, 0.467, 0.284, 0.881), O1 = c(-1.413, -0.574, -0.986, -0.068),
    O2 = c(-0.990, -0.504, -1.353, -0.384), O3 = c(0.173, 0.323, 1.091, -0.543), O4 = c(0.279, 1.682, -0.329, 0.514),
    O5 = c(-0.673, 0.078, -0.042, -0.986))
  cz <- rb$centralidade_z
  teste("tabela de centralidade do JASP inteira: 25 itens x 4 índices em z (intermediação, proximidade, força, influência esperada; spec_T6 P10; tol .0005)",
        perto(as.matrix(cz[, c("intermediacao", "proximidade", "forca", "influencia_esperada")]), unname(jt[it, ]), .0005))
  teste("influência esperada em z: N1 1.723 no topo, N3 1.478, N2 1.281 (JASP, spec_T6 P10; tol .0005)",
        perto(cz$influencia_esperada[match(c("N1", "N3", "N2"), it)], c(1.723, 1.478, 1.281), .0005) &&
          it[which.max(cz$influencia_esperada)] == "N1")
  teste("frase APA com os elementos do livro: EBICglasso, γ = 0.50, policóricas, 181 de 300, N1 e N2 com peso .57, C4 no topo da força e frase de alcance",
        grepl("EBICglasso; γ = 0.50", rb$apa, fixed = TRUE) && grepl("correlações policóricas", rb$apa, fixed = TRUE) &&
          grepl("Das 300 arestas possíveis, 181 foram estimadas como não nulas", rb$apa, fixed = TRUE) &&
          grepl("ligou N1 e N2 (peso = .57)", rb$apa, fixed = TRUE) && grepl("de C4 (1.65)", rb$apa, fixed = TRUE) &&
          grepl("não autorizam inferências causais ou de intervenção", rb$apa, fixed = TRUE))
  rv <- rede_psicometrica(bfi, it, decimal = ",", redundancia = FALSE)
  teste("decimal = \",\": γ = 0,50 (com zero), peso = ,57 e densidade = ,60 (sem zero), z = 1,65",
        grepl("γ = 0,50", rv$apa, fixed = TRUE) && grepl("(peso = ,57)", rv$apa, fixed = TRUE) &&
          grepl("densidade = ,60; ", rv$apa, fixed = TRUE) && grepl("C4 (1,65)", rv$apa, fixed = TRUE))
  teste("redundância antes da rede: UVA sinaliza N1–N2 (wTO ≥ .25), a aresta que o livro diz poder ser redação parecida",
        !is.null(rb$redundancia) && any(rb$redundancia$node_i == "N1" & rb$redundancia$node_j == "N2") &&
          any(grepl("UVA (EGAnet) sinalizou", rb$avisos, fixed = TRUE)))
  teste("checklist do livro com N, nós e arestas possíveis e os três pendentes da estabilidade",
        nrow(rb$checklist) == 7 && sum(grepl("pendente: estabilidade_rede", rb$checklist$situacao)) == 3 &&
          grepl("N = 2800; 25 nós; 300 arestas possíveis", rb$checklist$situacao[6], fixed = TRUE))

  # Equivalência e padrões
  bn <- suppressWarnings(bootnet::estimateNetwork(bfi[it], default = "EBICglasso", corMethod = "cor_auto", tuning = 0.5,
                                                  missing = "pairwise", sampleSize = "maximum", verbose = FALSE))
  teste("equivalência: grafo igual ao bootnet::estimateNetwork com os padrões do JASP (pesquisa_github A.5; tol 1e-12)",
        perto(unname(W), unname(bn$graph), 1e-12))
  bd <- suppressWarnings(bootnet::estimateNetwork(bfi[it], default = "EBICglasso", verbose = FALSE))
  teste("os padrões do bootnet (Pearson, pairwise_average) dão outra rede: por isso a função declara os do JASP",
        sum(bd$graph[upper.tri(bd$graph)] != 0) != 181 || max(abs(bd$graph - bn$graph)) > .01)

  # Bootstrap curto no BFI (20 reamostras, só arestas, por velocidade)
  eb <- estabilidade_rede(rb, n_boot = 20, executar = "arestas")
  todos$bfi_boot <- eb
  cp <- eb$comparacao_arestas
  teste("BFI, bootstrap de arestas (20 reamostras): N1–N2 se distingue de N3–N4 e nenhum IC de outra aresta toca o dela",
        cp$aresta_1 == "N1–N2" && cp$aresta_2 == "N3–N4" && isTRUE(cp$distinguiveis) && cp$n_ic_sobrepostos == 0)
  teste("BFI, IC bootstrap de N1–N2 contém o peso estimado e fica acima de .50", {
    a <- eb$arestas_ic[1, ]; a$ic_inf <= a$peso && a$peso <= a$ic_sup && a$ic_inf > .50 })

  # Ponte com os cinco fatores como comunidades
  pb <- ponte(rb, comunidades = substr(it, 1, 1))
  todos$bfi_ponte <- pb
  nb <- networktools::bridge(W, communities = substr(it, 1, 1))
  teste("ponte no BFI igual ao networktools::bridge (força, influência esperada 1 passo, intermediação; tol 1e-12)",
        perto(pb$tabela$forca_ponte, unname(nb[["Bridge Strength"]]), 1e-12) &&
          perto(pb$tabela$influencia_esperada_ponte, unname(nb[["Bridge Expected Influence (1-step)"]]), 1e-12) &&
          perto(pb$tabela$intermediacao_ponte, unname(nb[["Bridge Betweenness"]]), 1e-12))
}

# ============================================================ 2. Verdade conhecida (simulação)
cat("== Verdade conhecida: GGM com rede plantada\n")
nos <- paste0("x", 1:10)
# Duas cadeias de 5 nós (correlação parcial .35 entre vizinhos) ligadas por x5–x6 (.25).
Pv <- matrix(0, 10, 10, dimnames = list(nos, nos))
for (i in c(1:4, 6:9)) Pv[i, i + 1] <- Pv[i + 1, i] <- .35
Pv[5, 6] <- Pv[6, 5] <- .25
Kv <- diag(10) - Pv
Sig <- stats::cov2cor(solve(Kv))
gera <- function(n, semente) {
  set.seed(semente)
  x <- as.data.frame(MASS::mvrnorm(n, rep(0, 10), Sig)); names(x) <- nos; x
}
x <- gera(3000, 1)
r <- rede_psicometrica(x, nos, redundancia = FALSE)
todos$sim <- r
Wr <- r$pesos; vt <- Pv[upper.tri(Pv)]; ve <- Wr[upper.tri(Wr)]
teste("recupera as 9 arestas plantadas (sensibilidade 1) e r ≥ .95 entre pesos estimados e verdadeiros (n = 3000)",
      all(ve[vt != 0] != 0) && stats::cor(ve, vt) >= .95)
teste("toda aresta plantada supera em |peso| qualquer aresta ausente na verdade", min(abs(ve[vt != 0])) > max(abs(ve[vt == 0])))
tb <- r$tabela
teste("contagens coerentes: 45 possíveis, não nulas = soma de W ≠ 0 na metade superior, densidade = não nulas/45",
      tb$valor[tb$indicador == "arestas possíveis"] == 45 &&
        tb$valor[tb$indicador == "arestas não nulas"] == sum(ve != 0) &&
        perto(tb$valor[tb$indicador == "densidade"], sum(ve != 0) / 45, 1e-12))
cb <- r$centralidade_bruta
ct <- qgraph::centrality(Wr, weighted = TRUE, signed = TRUE, all.shortest.paths = FALSE)
teste("centralidade bruta: força = Σ|w|, influência esperada = Σw, intermediação = qgraph/2 (como o JASP; tol 1e-12)",
      perto(cb$forca, unname(rowSums(abs(Wr))), 1e-12) && perto(cb$influencia_esperada, unname(rowSums(Wr)), 1e-12) &&
        perto(cb$intermediacao, unname(ct$Betweenness) / 2, 1e-12) && perto(cb$proximidade, unname(ct$Closeness), 1e-12))
teste("z do \"Normalized\" do JASP = scale() de cada índice (tol 1e-12)",
      perto(r$centralidade_z$forca, as.numeric(scale(cb$forca)), 1e-12) &&
        perto(r$centralidade_z$intermediacao, as.numeric(scale(cb$intermediacao)), 1e-12))
fzs <- r$centralidade_z$forca; names(fzs) <- nos
teste("força: nós internos das cadeias (x2-x4, x7-x9) acima das pontas x1 e x10, como na rede verdadeira",
      min(fzs[c("x2", "x3", "x4", "x7", "x8", "x9")]) > max(fzs[c("x1", "x10")]))

cat("== Ausentes e tamanho de amostra\n")
xa <- x[1:1000, ]
set.seed(7); for (v in nos) xa[[v]][sample(1000, 40)] <- NA
xa[5, ] <- NA  # um caso sem nenhuma resposta
ra <- rede_psicometrica(xa, nos, redundancia = FALSE)
rl <- rede_psicometrica(xa, nos, ausentes = "listwise", redundancia = FALSE)
teste("par a par + \"maximum\": N = casos com ao menos uma resposta (999), igual ao bootnet:::sampleSize_pairwise; caso vazio sai e é contado",
      ra$n == 999 && ra$n_excluidos == 1 && ra$n == bootnet:::sampleSize_pairwise(xa[-5, ], "maximum"))
teste("listwise: N = casos completos e n_excluidos = incompletos",
      rl$n == sum(stats::complete.cases(xa)) && rl$n_excluidos == sum(!stats::complete.cases(xa)))
bpa <- suppressWarnings(bootnet::estimateNetwork(xa[-5, ], default = "EBICglasso", corMethod = "cor_auto", tuning = 0.5,
                                                 missing = "pairwise", sampleSize = "maximum", verbose = FALSE))
teste("com ausentes, a rede par a par é a do bootnet com os mesmos argumentos (tol 1e-12)", perto(unname(ra$pesos), unname(bpa$graph), 1e-12))

cat("== Precisão e estabilidade (bootstrap)\n")
e_g <- estabilidade_rede(r, n_boot = 100)  # 100 reamostras em cada bootstrap, por velocidade
todos$sim_boot <- e_g
xs <- gera(150, 2)
rs <- rede_psicometrica(xs, nos, redundancia = FALSE)
e_p <- estabilidade_rede(rs, n_boot = 100, executar = "casos")
teste("CS(força) acima de .50 com n = 3000 e menor com n = 150 na mesma rede (100 reamostras de casos)",
      e_g$cs[["strength"]] > .50 && e_p$cs[["strength"]] < e_g$cs[["strength"]])
teste("CS pela definição (r > .70 em mais de 95% das reamostras, maior proporção descartada), refeito à mão a partir do bootTable", {
  bc <- e_g$boot_casos
  orig <- as.data.frame(bc$sampleTable); orig <- orig[orig$type == "strength", c("node1", "value")]
  bt <- as.data.frame(bc$bootTable); bt <- bt[bt$type == "strength", ]
  bt$prop <- 1 - bt$nPerson / bc$sample$nPerson
  rr <- tapply(seq_len(nrow(bt)), bt$name, function(ii) {
    s <- bt[ii, ]; c(prop = s$prop[1], r = stats::cor(s$value, orig$value[match(s$node1, orig$node1)])) })
  rr <- do.call(rbind, rr)
  P <- tapply(rr[, "r"] > .7, rr[, "prop"], mean)
  cs_mao <- max(c(0, as.numeric(names(P))[P > .95]))
  perto(cs_mao, e_g$cs[["strength"]], 1e-9) })
teste("cortes do CS: < .25 não interpretar; .25 a .50 cautela; > .50 preferível; teto marcado; sem patamar \"excelente\"", {
  l <- .ee_redes_leitura_cs(c(.20, .25, .44, .52, .75), .75)
  grepl("abaixo do mínimo", l[1]) && grepl("cautela", l[2]) && grepl("cautela", l[3]) && grepl("preferível de .50$", l[4]) &&
    grepl("teto", l[5]) && !any(grepl("excelente", l)) })
teste("decisões declaram n_boot, semente, caseMin/caseMax/caseN e a definição do CS com r ≥ .70 e 95%",
      any(grepl("100 reamostras", e_g$decisoes)) && any(grepl("semente 2026", e_g$decisoes)) &&
        any(grepl("caseMin = 0.05, caseMax = 0.75, caseN = 10", e_g$decisoes, fixed = TRUE)) &&
        any(grepl("95% de probabilidade, correlação de pelo menos .70", e_g$decisoes, fixed = TRUE)))
teste("ressalva do livro: IC de aresta regularizada é descritivo, não teste, na decisão e na frase APA",
      any(grepl("não teste", e_g$decisoes)) && grepl("descritivos da variabilidade", e_g$apa[1], fixed = TRUE))
teste("aviso quando n_boot < 1000 (padrão do livro)", any(grepl("para relato, 1000 ou mais", e_g$avisos)))

# Aresta dominante (x1–x2 = .50, demais .20) contra arestas iguais (.35): a diferença
# entre a mais forte e a segunda deve excluir o zero só no primeiro caso.
Pd <- matrix(0, 6, 6, dimnames = list(nos[1:6], nos[1:6]))
for (i in 1:5) Pd[i, i + 1] <- Pd[i + 1, i] <- .20
Pd[1, 2] <- Pd[2, 1] <- .50
set.seed(3); xd <- as.data.frame(MASS::mvrnorm(1000, rep(0, 6), stats::cov2cor(solve(diag(6) - Pd)))); names(xd) <- nos[1:6]
ed <- estabilidade_rede(rede_psicometrica(xd, nos[1:6], redundancia = FALSE), n_boot = 80, executar = "arestas")
Pe <- matrix(0, 6, 6, dimnames = list(nos[1:6], nos[1:6]))
for (i in 1:5) Pe[i, i + 1] <- Pe[i + 1, i] <- .35
set.seed(4); xe <- as.data.frame(MASS::mvrnorm(1000, rep(0, 6), stats::cov2cor(solve(diag(6) - Pe)))); names(xe) <- nos[1:6]
ee <- estabilidade_rede(rede_psicometrica(xe, nos[1:6], redundancia = FALSE), n_boot = 80, executar = "arestas")
teste("aresta dominante (.50 contra .20) se distingue da segunda; com arestas iguais (.35) a diferença inclui o zero (80 reamostras, n = 1000)",
      ed$comparacao_arestas$aresta_1 == "x1–x2" && isTRUE(ed$comparacao_arestas$distinguiveis) &&
        !isTRUE(ee$comparacao_arestas$distinguiveis) && ee$comparacao_arestas$n_ic_sobrepostos >= 1)
teste("frase APA da comparação: \"distinguiu-se\" num caso e \"não se distinguiu\" no outro",
      grepl("distinguiu-se da segunda maior", ed$apa[1], fixed = TRUE) && grepl("não se distinguiu da segunda maior", ee$apa[1], fixed = TRUE))
eca <- estabilidade_rede(r, n_boot = 40, indices = c("strength", "closeness", "betweenness"), executar = "casos", caso_n = 5)
teste("proximidade e intermediação entram quando pedidas, com a ressalva de Bringmann et al. (2019) (40 reamostras, 5 níveis)",
      all(c("closeness", "betweenness") %in% eca$tabela$indice_bootnet) && any(grepl("Bringmann", eca$decisoes)) &&
        length(eca$niveis_descarte) == 5)

cat("== Ising, mgm, previsibilidade, ponte e redundância\n")
gI <- matrix(0, 6, 6); for (i in 1:5) gI[i, i + 1] <- gI[i + 1, i] <- 1
set.seed(1); xb <- as.data.frame(IsingSampler::IsingSampler(3000, gI, rep(-0.5, 6))); names(xb) <- paste0("b", 1:6)
ri <- rede_psicometrica(xb, names(xb), estimador = "IsingFit", redundancia = FALSE)
todos$ising <- ri
Wi <- ri$pesos
teste("Ising: recupera a cadeia plantada de 6 nós binários (5 arestas, todas positivas, nenhuma falsa; n = 3000)",
      all((Wi[upper.tri(Wi)] != 0) == (gI[upper.tri(gI)] != 0)) && all(Wi[gI != 0] > 0))
xb2 <- xb; xb2[] <- lapply(xb, function(v) v + 1L)
teste("Ising com itens codificados 1/2: recodifica para 0/1 e dá a mesma rede",
      perto(rede_psicometrica(xb2, names(xb2), estimador = "IsingFit", redundancia = FALSE)$pesos, Wi, 1e-12))
teste("Ising recusa itens com mais de 2 categorias", inherits(tryCatch(rede_psicometrica(x, nos, estimador = "IsingFit"),
                                                                     error = function(e) e), "error"))
bI <- bootnet::estimateNetwork(xb, default = "IsingFit", tuning = 0.5, rule = "AND", split = "median", verbose = FALSE)
teste("equivalência: Ising igual ao bootnet::estimateNetwork(default = \"IsingFit\", tuning = 0.5) (tol 1e-12)",
      perto(unname(Wi), unname(bI$graph), 1e-12))
xm <- x[, 1:6]; xm$x6 <- as.integer(xm$x6 > 0)
rm6 <- rede_psicometrica(xm, names(xm), estimador = "mgm", redundancia = FALSE)
todos$mgm <- rm6
Wm <- rm6$pesos; cad <- cbind(1:5, 2:6)
teste("mgm com 5 gaussianos e 1 binário: tipos inferidos (g, c) e as 5 arestas da cadeia superam todas as outras",
      identical(rm6$tipos, c(rep("g", 5), "c")) && min(Wm[cad]) > max(abs(Wm[upper.tri(Wm)][-c(1, 3, 6, 10, 15)])))
x5k <- gera(5000, 5)
r5k <- rede_psicometrica(x5k, nos, redundancia = FALSE)
pv <- previsibilidade(r5k)
todos$prev <- pv
r2_teo <- 1 - 1 / diag(solve(Sig))
teste("previsibilidade: R² do mgm a .03 do R² verdadeiro de cada nó, 1 - 1/K_ii (n = 5000)", perto(pv$tabela$r2, unname(r2_teo), .03))
teste("previsibilidade: R² implicado pela rede EBICglasso a .05 do verdadeiro (encolhido pelo LASSO, nunca acima + .01)",
      perto(pv$tabela$r2_ggm_implicito, unname(r2_teo), .05) && all(pv$tabela$r2_ggm_implicito <= r2_teo + .01))
fit <- mgm::mgm(as.matrix(x5k), type = rep("g", 10), level = rep(1, 10), lambdaSel = "EBIC", lambdaGam = 0.5,
                ruleReg = "AND", k = 2, pbar = FALSE, signInfo = FALSE)
pr <- stats::predict(fit, as.matrix(x5k), errorCon = "R2")
teste("equivalência: previsibilidade igual ao mgm + predict chamados direto (γ = 0.5, regra AND; tol 1e-12)",
      perto(pv$tabela$r2, as.numeric(pr$errors$R2), 1e-12))
pvi <- previsibilidade(ri)
teste("previsibilidade em Ising: acurácia (CC) e acurácia normalizada (nCC) por nó, entre 0 e 1",
      all(!is.na(pvi$tabela$cc)) && all(pvi$tabela$cc > 0 & pvi$tabela$cc <= 1) && all(pvi$tabela$ncc >= 0 & pvi$tabela$ncc <= 1))
pn <- ponte(Pv, comunidades = list(a = nos[1:5], b = nos[6:10]))
teste("ponte na rede verdadeira: só x5 e x6 têm força de ponte (= .25, a aresta entre as cadeias), os demais zero",
      perto(pn$tabela$forca_ponte, c(0, 0, 0, 0, .25, .25, 0, 0, 0, 0), 1e-12))
teste("ponte recusa partição com uma comunidade só", inherits(tryCatch(ponte(Pv, rep("a", 10)), error = function(e) e), "error"))
xr <- x[1:1000, ]; set.seed(9); xr$x11 <- xr$x1 + stats::rnorm(1000, 0, .15)
rr <- rede_psicometrica(xr, c(nos, "x11"))
teste("redundância: UVA sinaliza um item quase duplicado (x1 e x11) antes da rede",
      !is.null(rr$redundancia) && any((rr$redundancia$node_i == "x1" & rr$redundancia$node_j == "x11") |
                                        (rr$redundancia$node_i == "x11" & rr$redundancia$node_j == "x1")))

# ============================================================ 3. Equivalência com bootnet
cat("== Equivalência com bootnet\n")
s_edge <- as.data.frame(summary(e_g$boot_arestas, statistics = "edge"))
teste("IC das arestas = quantis tipo 6 do bootTable, refeitos à mão (os mesmos do summary do bootnet; tol 1e-12)", {
  bt <- as.data.frame(e_g$boot_arestas$bootTable); bt <- bt[bt$type == "edge", ]
  q <- t(vapply(e_g$arestas_ic$id, function(id) stats::quantile(bt$value[bt$id == id], c(.025, .975), type = 6, names = FALSE), numeric(2)))
  perto(q[, 1], e_g$arestas_ic$ic_inf, 1e-12) && perto(q[, 2], e_g$arestas_ic$ic_sup, 1e-12) &&
    perto(e_g$arestas_ic$ic_inf, s_edge$q2.5[match(e_g$arestas_ic$id, s_edge$id)], 1e-12) })
teste("CS igual ao bootnet::corStability do mesmo objeto",
      perto(e_g$cs, bootnet::corStability(e_g$boot_casos, verbose = FALSE)[names(e_g$cs)], 1e-12))
e30 <- estabilidade_rede(r, n_boot = 30, executar = "arestas")
set.seed(2026)
b30 <- suppressMessages(bootnet::bootnet(r$modelo, nBoots = 30, type = "nonparametric", statistics = c("edge", "strength", "expectedInfluence"), verbose = FALSE))
s30 <- as.data.frame(summary(b30, statistics = "edge"))
teste("semente: bootnet::bootnet com set.seed(2026) e os mesmos argumentos dá os mesmos IC (30 reamostras; tol 1e-12)",
      perto(e30$arestas_ic$ic_inf, s30$q2.5[match(e30$arestas_ic$id, s30$id)], 1e-12) &&
        perto(e30$arestas_ic$ic_sup, s30$q97.5[match(e30$arestas_ic$id, s30$id)], 1e-12))
teste("diferença de |peso| entre as duas maiores = bootnet::differenceTest(measure = \"edge\") quando as duas são positivas", {
  cp <- e_g$comparacao_arestas
  a1 <- e_g$arestas_ic[e_g$arestas_ic$aresta == cp$aresta_1, ]; a2 <- e_g$arestas_ic[e_g$arestas_ic$aresta == cp$aresta_2, ]
  dt <- bootnet::differenceTest(e_g$boot_arestas, x = a2$no1, x2 = a2$no2, y = a1$no1, y2 = a1$no2, measure = "edge", verbose = FALSE)
  a1$peso > 0 && a2$peso > 0 && perto(c(dt$lower, dt$upper), c(cp$ic_inf, cp$ic_sup), 1e-12) })

# ============================================================ 4. Formatação, idioma e privacidade
cat("== Formatação, idioma e privacidade\n")
textos <- unlist(lapply(todos, function(o) c(o$apa, o$decisoes, o$avisos)))
teste("nenhum texto gerado usa travessão longo", !any(grepl("\u2014", textos)))
teste("nenhum texto gerado usa as palavras vetadas", !any(grepl("genuinamente|de fato|notavelmente|é importante notar|intricad|meticulos", textos, ignore.case = TRUE)))
teste("nenhum texto gerado cria patamar \"excelente\" para o CS nem chama centralidade de alvo de intervenção comprovado",
      !any(grepl("excelente|excellent", textos, ignore.case = TRUE)) &&
        !any(grepl("(?<!não )é alvo de intervenção", textos, perl = TRUE)) &&
        any(grepl("não é alvo de intervenção comprovado", textos, fixed = TRUE)))
ren <- rede_psicometrica(x, nos, idioma = "en", redundancia = FALSE)
een <- estabilidade_rede(ren, n_boot = 20, executar = "casos", caso_n = 4, idioma = "en")
teste("idioma = \"en\": frases em inglês com a frase de alcance e o CS",
      grepl("^The network of the 10 items", ren$apa) && grepl("do not support causal or intervention inferences", ren$apa, fixed = TRUE) &&
        grepl("CS(strength) = ", een$apa, fixed = TRUE))
teste("print funciona e nenhuma tabela impressa tem uma linha por participante", {
  out <- utils::capture.output(for (o in todos) print(o))
  tabs <- unlist(lapply(todos, function(o) c(list(o$tabela, o$efeito, o$posthoc), o$tabelas_extra)), recursive = FALSE)
  length(out) > 50 && all(vapply(tabs, function(t) is.null(t) || NROW(t) < 50, logical(1))) })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
if (!nzchar(DJ)) cat("Aviso: PR_DADOS_JASP vazia; as âncoras do livro não foram testadas.\n")
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
