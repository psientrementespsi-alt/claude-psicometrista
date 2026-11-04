# evals_simulados.R, testes da skill psicometria-r em dados simulados
#
# Rodar no R do computador do usuário, da pasta da skill:
#   Rscript scripts/evals_simulados.R
# Saída: uma linha PASSOU/FALHOU por teste e o resumo. Código de saída 1 se
# algum teste falhar. Nenhum dado real é usado; tudo vem de lavaan::simulateData
# e mirt::simdata com sementes fixas. Os dados ordinais são simulados com
# ordered_center = FALSE: o padrão do lavaan 0.7 centra cada resposta latente na
# média da própria amostra antes de cortar nos limiares, o que apagaria a
# variação de localização entre grupos e deixaria os testes de invariância fáceis.

args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
source(file.path(pasta, "psicometria.R"))
source(file.path(pasta, "est_comum.R"))  # checar_escore_somado (v1.3) usa est_comum.R
suppressPackageStartupMessages({
  library(lavaan); library(semTools); library(mirt); library(bootnet)
  library(EGAnet); library(NetworkComparisonTest)
})

resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome))
  resultados[[nome]] <<- ok
  invisible(ok)
}
t_inicio <- Sys.time()
# Nos cenários em que a hipótese nula é verdadeira (grupos invariantes), cada
# teste de diferença tem 5% de chance de rejeitar a .05 por acaso; com três
# passos, perto de 14%. Os evals usam alfa = .01 nesses cenários para não
# confundir erro tipo I esperado com defeito, e imprimem os p de cada passo.
ALFA_NULO <- .01
sim <- function(modelo, n, seed) simulateData(modelo, sample.nobs = n, seed = seed, ordered_center = FALSE)

# Dados: 3 fatores, 4 itens cada, carga padronizada .70 na resposta latente,
# correlação .30 entre fatores, 4 categorias (limiares -1, 0, 1).
itens <- paste0("y", 1:12)
limiares <- function(vars, desloc = 0) paste(sprintf("%s | %.2f*t1 + %.2f*t2 + %.2f*t3", vars,
                                                     -1 + desloc, 0 + desloc, 1 + desloc), collapse = "\n")
pop <- paste("f1 =~ .7*y1 + .7*y2 + .7*y3 + .7*y4",
             "f2 =~ .7*y5 + .7*y6 + .7*y7 + .7*y8",
             "f3 =~ .7*y9 + .7*y10 + .7*y11 + .7*y12",
             "f1 ~~ .3*f2 + .3*f3", "f2 ~~ .3*f3", limiares(itens), sep = "\n")
modelo <- "f1 =~ y1 + y2 + y3 + y4\n f2 =~ y5 + y6 + y7 + y8\n f3 =~ y9 + y10 + y11 + y12"
dados <- sim(pop, 1000, 2026)
N <- nrow(dados)

cat("== Auditoria e escolha do estimador\n")
aud <- auditar_itens(dados, itens)
teste("auditoria: uma linha por item, só agregados", nrow(aud) == length(itens) && all(aud$categorias == 4) &&
        all(vapply(aud, function(col) length(col) == length(itens), logical(1))))
teste("4 categorias -> WLSMV ordinal", { e <- sugerir_estimador(aud); e$estimador == "WLSMV" && e$ordinal })
cont <- as.data.frame(lapply(simulateData("f =~ .7*a + .7*b + .7*c", sample.nobs = 300, seed = 1),
                             function(x) round(50 + 15 * x)))
teste("escala 0-100 -> MLR contínuo", { e <- sugerir_estimador(auditar_itens(cont, c("a", "b", "c"))); e$estimador == "MLR" && !e$ordinal })
cinco <- dados; cinco[itens] <- lapply(dados[itens], function(x) x + (stats::runif(N) > .8))
teste("5 categorias simétricas, n = 1000 -> WLSMV principal com MLR de sensibilidade",
      { e <- sugerir_estimador(auditar_itens(cinco, itens)); e$estimador == "WLSMV" && identical(e$sensibilidade, "MLR") &&
          grepl("amostra pequena", e$motivo) })
teste("5 categorias simétricas, n = 150 -> MLR principal (Rhemtulla: ML com amostra pequena e limiares simétricos)",
      { e <- sugerir_estimador(auditar_itens(cinco[1:150, ], itens)); e$estimador == "MLR" && identical(e$sensibilidade, "WLSMV") })
assim <- paste("f =~ .7*z1 + .7*z2 + .7*z3 + .7*z4",
               paste(sprintf("z%d | 0.5*t1 + 1.0*t2 + 1.5*t3 + 2.0*t4", 1:4), collapse = "\n"), sep = "\n")
d_assim <- sim(assim, 1000, 3)
teste("5 categorias assimétricas -> WLSMV, motivo cita limiares assimétricos",
      { a <- auditar_itens(d_assim, paste0("z", 1:4)); e <- sugerir_estimador(a)
        all(a$categorias == 5) && e$estimador == "WLSMV" && grepl("assimétricos", e$motivo) })

cat("== AFC ordinal\n")
fit <- ajustar_afc(modelo, dados, itens)
std <- standardizedSolution(fit)
teste("cargas padronizadas recuperadas (.70 +- .08)", all(abs(std$est.std[std$op == "=~"] - .7) < .08))
teste("correlações entre fatores recuperadas (.30 +- .08)", all(abs(std$est.std[std$op == "~~" & std$lhs != std$rhs] - .3) < .08))
taj <- tabela_ajuste(fit, "3 fatores")
teste("tabela de ajuste traz versão escalonada e robusta", setequal(taj$versao, c("escalonada", "robusta")))
rob <- taj[taj$versao == "robusta", ]
teste("modelo correto: CFI robusto >= .95 e RMSEA robusto <= .06", rob$cfi >= .95 && rob$rmsea <= .06)
fit1 <- ajustar_afc("g =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 + y11 + y12", dados, itens)
rob1 <- tabela_ajuste(fit1, "1 fator")
rob1 <- rob1[rob1$versao == "robusta", ]
teste("modelo errado (1 fator) ajusta pior", rob1$cfi < rob$cfi - .05 && rob1$rmsea > rob$rmsea)

cat("== Confiabilidade e validade\n")
conf <- confiabilidade(fit, dados)
teste("ômega e alfa por fator, entre 0 e 1", nrow(conf) == 3 && all(conf$omega > 0 & conf$omega < 1 & conf$alfa > 0 & conf$alfa < 1))
om_lat <- unlist(compRelSEM(fit, ord.scale = FALSE))
teste("ômega da resposta latente perto do teórico (.794 +- .04)", all(abs(om_lat - 0.7935) < .04))
teste("ômega do escore observado abaixo do latente (categorização perde informação)", all(conf$omega < om_lat))
vd <- validade_discriminante(fit, dados)
teste("AVE perto do teórico (.49 +- .07)", all(abs(vd$ave - .49) < .07))
teste("HTMT < .85 e Fornell-Larcker atendido com fatores a .30", all(vd$htmt[upper.tri(vd$htmt)] < .85) && vd$fornell_larcker_ok)
pop90 <- paste("f1 =~ .7*y1 + .7*y2 + .7*y3 + .7*y4", "f2 =~ .7*y5 + .7*y6 + .7*y7 + .7*y8",
               "f1 ~~ .9*f2", limiares(paste0("y", 1:8)), sep = "\n")
d90 <- sim(pop90, 1000, 90)
vd90 <- validade_discriminante(ajustar_afc("f1 =~ y1 + y2 + y3 + y4\n f2 =~ y5 + y6 + y7 + y8", d90, paste0("y", 1:8)), d90)
teste("fatores a .90: HTMT > .85 e Fornell-Larcker reprovado", vd90$htmt[1, 2] > .85 && identical(vd90$fornell_larcker_ok, FALSE))
vd1 <- suppressMessages(validade_discriminante(fit1, dados))
teste("um fator só: Fornell-Larcker sai NA, não TRUE", is.na(vd1$fornell_larcker_ok))
teste("relato APA sem zero à esquerda", grepl("CFI = \\.9[0-9]{2}, TLI = \\.9", formatar_ajuste_apa(rob)) &&
        grepl("χ²\\([0-9]+\\)", formatar_ajuste_apa(rob)))
teste("relato com vírgula decimal e ponto e vírgula entre itens",
      { f <- formatar_ajuste_apa(rob, decimal = ","); grepl("CFI = ,9", f) && grepl("; TLI = ", f) && grepl("\\[,[0-9]{3}; ,[0-9]{3}\\]", f) })
teste("formato de p nos extremos e sem -0", .pr_p(.99996) == "p > .999" && .pr_p(.0004) == "p < .001" &&
        .pr_num(-0.0001, 3, FALSE) == ".000")
cat("  exemplo de relato (dados simulados):", formatar_ajuste_apa(rob), "\n")

cat("== Segunda ordem no formato do MPFI-24 (12 fatores de 2 itens, 6 categorias)\n")
# Tudo especificado na escala padronizada: carga de 1a ordem .75, de 2a ordem
# .70, correlação -.40 entre os fatores de 2a ordem, 5 limiares simétricos.
it24 <- paste0("m", 1:24)
fat12 <- paste0("s", 1:12)
pop2 <- c(
  sprintf("%s =~ .75*%s + .75*%s", fat12, it24[seq(1, 23, 2)], it24[seq(2, 24, 2)]),
  sprintf("flex =~ %s", paste0(".7*", fat12[1:6], collapse = " + ")),
  sprintf("inflex =~ %s", paste0(".7*", fat12[7:12], collapse = " + ")),
  "flex ~~ 1*flex", "inflex ~~ 1*inflex", "flex ~~ -.4*inflex",
  sprintf("%s ~~ .51*%s", fat12, fat12),
  sprintf("%s ~~ .4375*%s", it24, it24),
  sprintf("%s | -1.5*t1 + -0.75*t2 + 0*t3 + 0.75*t4 + 1.5*t5", it24))
mod2 <- paste(c(sprintf("%s =~ %s + %s", fat12, it24[seq(1, 23, 2)], it24[seq(2, 24, 2)]),
                sprintf("flex =~ %s", paste(fat12[1:6], collapse = " + ")),
                sprintf("inflex =~ %s", paste(fat12[7:12], collapse = " + "))), collapse = "\n")
d24 <- sim(paste(pop2, collapse = "\n"), 1000, 31)
e24 <- sugerir_estimador(auditar_itens(d24, it24))
teste("6 categorias simétricas -> WLSMV com MLR de sensibilidade", e24$estimador == "WLSMV" && identical(e24$sensibilidade, "MLR"))
fit24 <- ajustar_afc(mod2, d24, it24)
t24 <- tabela_ajuste(fit24, "2a ordem")
teste("2a ordem: ajuste robusto bom (CFI >= .95, RMSEA <= .06)",
      t24$cfi[t24$versao == "robusta"] >= .95 && t24$rmsea[t24$versao == "robusta"] <= .06)
s24 <- standardizedSolution(fit24)
teste("2a ordem: cargas de 2a ordem recuperadas (.70 +- .10)",
      all(abs(s24$est.std[s24$op == "=~" & s24$lhs %in% c("flex", "inflex")] - .7) < .10))
teste("2a ordem: correlação entre flex e inflex recuperada (-.40 +- .10)",
      abs(s24$est.std[s24$op == "~~" & s24$lhs == "flex" & s24$rhs == "inflex"] + .4) < .10)
# Valores teóricos na resposta latente para o composto de 12 itens:
# ômega hierárquico (só o fator geral) = 39,69 / 51,825 = .766;
# ômega total (todos os fatores comuns) = 46,575 / 51,825 = .899.
comp <- list(flex = it24[1:12], inflex = it24[13:24])
rel_h <- confiabilidade_compostos(fit24, comp, ord.scale = FALSE)
teste("2a ordem: ômega do composto no modelo de 2a ordem = hierárquico (.766 +- .04)",
      all(abs(rel_h$omega - .766) < .04))
mod_corr <- paste(sprintf("%s =~ %s + %s", fat12, it24[seq(1, 23, 2)], it24[seq(2, 24, 2)]), collapse = "\n")
fit_corr <- ajustar_afc(mod_corr, d24, it24)
rel_t <- confiabilidade_compostos(fit_corr, comp, ord.scale = FALSE)
teste("fatores correlacionados: ômega do composto = total (.899 +- .03)", all(abs(rel_t$omega - .899) < .03))
rel_obs <- confiabilidade_compostos(fit24, comp)
teste("ômega do composto observado (ordinal) abaixo do latente", all(rel_obs$omega < rel_h$omega))
sub_2a <- suppressMessages(confiabilidade(fit24, d24))
sub_corr <- confiabilidade(fit_corr, d24)
teste("subescalas: confiabilidade() roda em modelo de 2a ordem (12 subescalas, sem os fatores gerais)",
      nrow(sub_2a) == 12 && !any(c("flex", "inflex") %in% sub_2a$fator))
teste("subescalas: ômega no modelo de 2a ordem (parte do fator geral) < ômega no modelo correlacionado",
      all(sub_2a$omega < sub_corr$omega))
vd24 <- validade_discriminante(fit24, d24)
teste("2a ordem: validade_discriminante() roda nos 12 fatores de 1a ordem", length(vd24$ave) == 12 &&
        all(dim(vd24$htmt) == c(12, 12)) && isTRUE(vd24$segunda_ordem) && !is.na(vd24$fornell_larcker_ok))
teste("2a ordem: invariancia() recusa modelo com fator de ordem superior, com mensagem",
      { r <- tryCatch(invariancia(mod2, d24, it24, "s1"), error = function(e) conditionMessage(e))
        is.character(r) && grepl("ordem superior", r) })

cat("== Invariância ordinal (Wu e Estabrook, 2016)\n")
pop1 <- paste("f =~ .7*y1 + .7*y2 + .7*y3 + .7*y4 + .7*y5 + .7*y6", limiares(paste0("y", 1:6)), sep = "\n")
# y1 e y2 no grupo B: (a) todos os limiares deslocados em +.8 (deslocamento
# uniforme, equivale a mudar o intercepto da resposta latente); (b) só o
# primeiro limiar deslocado em +.8 (muda o espaçamento entre limiares).
uniforme <- paste("f =~ .7*y1 + .7*y2 + .7*y3 + .7*y4 + .7*y5 + .7*y6",
                  limiares(c("y1", "y2"), .8), limiares(paste0("y", 3:6)), sep = "\n")
espaco <- paste("f =~ .7*y1 + .7*y2 + .7*y3 + .7*y4 + .7*y5 + .7*y6",
                "y1 | -0.2*t1 + 0*t2 + 1*t3", "y2 | -0.2*t1 + 0*t2 + 1*t3",
                limiares(paste0("y", 3:6)), sep = "\n")
mod1 <- "f =~ y1 + y2 + y3 + y4 + y5 + y6"
it6 <- paste0("y", 1:6)
dois_grupos <- function(pop_b, s) rbind(cbind(sim(pop1, 600, s), g = "A"),
                                        cbind(sim(pop_b, 600, s + 1), g = "B"))
g_a <- dois_grupos(pop1, 11)
mostra_inv <- function(t) print(t[, c("modelo", "cfi", "rmsea", "delta_cfi", "delta_rmsea", "p_dif_qui2")], digits = 3)
inv_ok <- invariancia(mod1, g_a, it6, "g")$tabela; mostra_inv(inv_ok)
teste("grupos invariantes: |ΔCFI| <= .010 em todos os passos", all(abs(inv_ok$delta_cfi[-1]) <= .010))
teste("grupos invariantes: nenhum passo rejeitado (p > .01)", all(inv_ok$p_dif_qui2[-1] > ALFA_NULO))
inv_u <- invariancia(mod1, dois_grupos(uniforme, 13), it6, "g")$tabela; mostra_inv(inv_u)
teste("deslocamento uniforme dos limiares: passa nos limiares e cai nos interceptos (p < .001, ΔCFI < -.010)",
      inv_u$p_dif_qui2[2] > .05 && inv_u$p_dif_qui2[4] < .001 && inv_u$delta_cfi[4] < -.010)
inv_e <- invariancia(mod1, dois_grupos(espaco, 15), it6, "g")$tabela; mostra_inv(inv_e)
teste("mudança de espaçamento entre limiares: cai no passo dos limiares (p < .001)", inv_e$p_dif_qui2[2] < .001)
media <- paste(pop1, "f ~ 0.5*1", sep = "\n")  # média latente .5 maior no grupo B, itens invariantes
inv_m <- invariancia(mod1, dois_grupos(media, 17), it6, "g")$tabela; mostra_inv(inv_m)
teste("diferença real de média latente com itens invariantes: nenhum passo rejeitado (p > .01, |ΔCFI| <= .010)",
      all(inv_m$p_dif_qui2[-1] > ALFA_NULO) && all(abs(inv_m$delta_cfi[-1]) <= .010))
aud_g <- auditar_itens(g_a, it6, grupo = "g")
teste("auditoria por grupo traz a menor categoria dentro dos grupos", all(c("menor_categoria_no_grupo", "categoria_vazia_em_grupo") %in% names(aud_g)))
g_niv <- g_a; g_niv$g <- factor(g_niv$g, levels = c("A", "B", "C"))  # nível C sem casos
teste("auditoria por grupo ignora nível de grupo sem casos", !any(auditar_itens(g_niv, it6, grupo = "g")$categoria_vazia_em_grupo))
teste("validade_discriminante() recusa modelo multigrupo com mensagem",
      { r <- tryCatch(validade_discriminante(cfa(mod1, g_a, group = "g", ordered = it6), g_a), error = function(e) conditionMessage(e))
        is.character(r) && grepl("um grupo", r) })
tres <- paste("f =~ .7*y1 + .7*y2 + .7*y3 + .7*y4", paste(sprintf("y%d | -0.5*t1 + 0.5*t2", 1:4), collapse = "\n"), sep = "\n")
g3 <- rbind(cbind(sim(tres, 600, 41), g = "A"), cbind(sim(tres, 600, 42), g = "B"))
inv3 <- suppressMessages(invariancia("f =~ y1 + y2 + y3 + y4", g3, paste0("y", 1:4), "g"))$tabela
mostra_inv(inv3)
teste("3 categorias: o passo dos limiares não ganha grau de liberdade e fica sem p, sem teste inventado",
      inv3$delta_gl[2] == 0 && is.na(inv3$p_dif_qui2[2]) && all(inv3$p_dif_qui2[3:4] > ALFA_NULO))

cat("== Caminho contínuo (escala 0-100, MLR)\n")
pop_c <- "f =~ .7*a1 + .7*a2 + .7*a3 + .7*a4 + .7*a5 + .7*a6"
it_c <- paste0("a", 1:6)
para_100 <- function(d) { d[it_c] <- lapply(d[it_c], function(x) pmin(100, pmax(0, round(50 + 15 * x)))); d }
dc <- rbind(cbind(para_100(simulateData(pop_c, sample.nobs = 500, seed = 21)), g = "A"),
            cbind(para_100(simulateData(pop_c, sample.nobs = 500, seed = 22)), g = "B"))
e_c <- sugerir_estimador(auditar_itens(dc, it_c))
fit_c <- ajustar_afc("f =~ a1 + a2 + a3 + a4 + a5 + a6", dc, it_c, e_c$estimador, e_c$ordinal)
tc <- tabela_ajuste(fit_c, "0-100")
teste("MLR: estimador sugerido e ajuste robusto bom", e_c$estimador == "MLR" &&
        tc$cfi[tc$versao == "robusta"] >= .95)
inv_c <- invariancia("f =~ a1 + a2 + a3 + a4 + a5 + a6", dc, it_c, "g", ordinal = FALSE)$tabela; mostra_inv(inv_c)
teste("MLR: invariância contínua roda os 4 passos sem rejeição", identical(inv_c$modelo, c("configural", "metrica", "escalar", "estrita")) &&
        all(abs(inv_c$delta_cfi[-1]) <= .010) && all(inv_c$p_dif_qui2[-1] > ALFA_NULO))

cat("== TRI: modelo de resposta graduada (Samejima, 1969) e DIF\n")
set.seed(7)
a_true <- matrix(runif(8, 1, 2.5))
d_true <- t(apply(matrix(rnorm(8 * 3), 8), 1, function(r) sort(r, decreasing = TRUE)))
resp <- simdata(a_true, d_true, N = 1000, itemtype = "graded")
grm <- mirt(resp, 1, itemtype = "graded", verbose = FALSE)
a_est <- coef(grm, simplify = TRUE)$items[, "a1"]
teste("discriminações recuperadas (r > .85)", cor(a_est, a_true) > .85)
m2 <- M2(grm, type = "C2")
teste("M2 (C2) calculável, com RMSEA", is.finite(m2$RMSEA))
ifit <- itemfit(grm, fit_stats = "S_X2")
teste("S-X2 por item com p-valor", all(is.finite(ifit$p.S_X2)))
d_b <- d_true; d_b[1, ] <- d_b[1, ] - 1
grupo <- rep(c("A", "B"), each = 800)
resp_g <- rbind(simdata(a_true, d_true, N = 800, itemtype = "graded"),
                simdata(a_true, d_b, N = 800, itemtype = "graded"))
mg <- multipleGroup(resp_g, 1, group = grupo, itemtype = "graded", verbose = FALSE,
                    invariance = c(colnames(resp_g)[-1], "free_means", "free_var"))
dif <- DIF(mg, which.par = c("a1", "d1", "d2", "d3"), items2test = 1, verbose = FALSE)
teste("DIF plantado no item 1 detectado (p < .01)", dif$p[1] < .01)
mg2 <- multipleGroup(resp_g, 1, group = grupo, itemtype = "graded", verbose = FALSE,
                     invariance = c(colnames(resp_g)[4:8], "free_means", "free_var"))
dif2 <- DIF(mg2, which.par = c("a1", "d1", "d2", "d3"), items2test = 2:3, verbose = FALSE)
teste("itens sem DIF (2 e 3) não são acusados (p > .01)", all(dif2$p > ALFA_NULO))

perto <- function(a, b, tol) isTRUE(all(abs(a - b) <= tol + 1e-12))
cat("== Escore somado: paralelo contra congenérico (checar_escore_somado, v1.3)\n")
# Âncora: exemplo de McNeish e Wolf (2020) com o banco público Holzinger e Swineford
# (1939) do lavaan, itens x4 a x9 arredondados para o inteiro mais próximo, como no
# artigo. O arredondamento é meio para cima: com o round() do R, que leva o meio ao
# par, a média do item 2 sai 4,35 em vez de 4,46 e os qui-quadrados não batem. Valores
# do texto publicado: congenérico de 1 fator χ²(9) = 115,37, CFI = 0,81, RMSEA = 0,20
# [0,17; 0,23]; cargas padronizadas (Tabela 2) .82 .85 .79 .17 .18 .26; paralelo de 1
# fator CFI = 0,45, RMSEA = 0,23 [0,21; 0,25], SRMR = 0,24; 2 fatores congenérico
# χ²(8) = 14,74; 2 fatores paralelo sem covariância χ²(17) = 55,54, diferença χ²(9) =
# 40,80. O SRMR do artigo é o do Mplus (srmr_mplus no lavaan). Tolerância: metade da
# última casa impressa, mais 0,01 nos qui-quadrados (o artigo trunca 40,808 em 40,80).
# O χ²(19) = 5361,86 do paralelo de 1 fator no texto não é compatível com o RMSEA e o
# CFI que o próprio artigo dá para esse modelo; o lavaan dá 325,90, que reproduz os dois.
hs <- lavaan::HolzingerSwineford1939
hs6 <- as.data.frame(lapply(hs[paste0("x", 4:9)], function(x) floor(x + 0.5)))
mw1 <- checar_escore_somado(hs6, paste0("x", 4:9), estimador = "ML")
aj1 <- mw1$ajuste
teste("S-01 âncora McNeish e Wolf (2020): congenérico de 1 fator χ²(9) = 115,37, CFI = .81, RMSEA = .20 [.17, .23]",
      local({ l <- aj1[aj1$modelo == "congenérico", ]
        perto(l$qui2, 115.37, .005 + .01) && l$gl == 9 && perto(l$cfi, .81, .005) && perto(l$rmsea, .20, .005) &&
          perto(c(l$rmsea_li, l$rmsea_ls), c(.17, .23), .005) }))
teste("S-02 âncora: cargas padronizadas do congenérico = Tabela 2 do artigo (.82 .85 .79 .17 .18 .26)",
      perto(mw1$cargas$carga_padronizada, c(.82, .85, .79, .17, .18, .26), .005))
teste("S-03 âncora: paralelo de 1 fator CFI = .45, RMSEA = .23 [.21, .25], SRMR do Mplus = .24; restrições da soma rejeitadas",
      local({ l <- aj1[aj1$modelo == "paralelo", ]
        f <- mw1$fits$escala$paralelo
        perto(l$cfi, .45, .005) && perto(l$rmsea, .23, .005) && perto(c(l$rmsea_li, l$rmsea_ls), c(.21, .25), .005) &&
          perto(unname(lavaan::fitMeasures(f, "srmr_mplus")), .24, .005) && isTRUE(mw1$soma_rejeitada[["escala"]]) &&
          mw1$tabela$metodo[3] == "standard" }))
teste("S-04 âncora: modelo do artigo para subescalas (2 fatores, cargas 1, erros iguais dentro da subescala, sem covariância) dá χ²(17) = 55,54 e diferença χ²(9) = 40,80",
      local({ c2 <- lavaan::cfa("v =~ x4 + x5 + x6\n s =~ x7 + x8 + x9", hs6, estimator = "ML", std.lv = TRUE, meanstructure = TRUE)
        p2 <- lavaan::cfa(paste("v =~ 1*x4 + 1*x5 + 1*x6", "s =~ 1*x7 + 1*x8 + 1*x9", "v ~~ 0*s",
                                "x4 ~~ a*x4\n x5 ~~ a*x5\n x6 ~~ a*x6\n x7 ~~ b*x7\n x8 ~~ b*x8\n x9 ~~ b*x9", sep = "\n"),
                          hs6, estimator = "ML", meanstructure = TRUE)
        lr <- lavaan::lavTestLRT(c2, p2)
        perto(unname(lavaan::fitMeasures(c2, "chisq")), 14.74, .015) && perto(unname(lavaan::fitMeasures(p2, "chisq")), 55.54, .015) &&
          perto(lr[2, "Chisq diff"], 40.80, .015) && lr[2, "Df diff"] == 9 }))
mw2 <- checar_escore_somado(hs6, list(Verbal = c("x4", "x5", "x6"), Speeded = c("x7", "x8", "x9")), estimador = "ML")
teste("S-05 subescalas do artigo testadas uma a uma: 4 gl por subescala (paralelo de 3 itens contra congenérico saturado), decisão por subescala e nota da covariância",
      local({ d <- mw2$tabela[mw2$tabela$decisiva, ]
        all(d$gl_dif == 4) && all(names(mw2$soma_rejeitada) == c("Verbal", "Speeded")) && any(grepl("sem covariância", mw2$decisoes)) &&
          length(mw2$apa) == 2 && all(grepl("^Em (Verbal|Speeded), o modelo paralelo", mw2$apa)) }))
teste("S-06 subescala Verbal: cargas padronizadas quase iguais (amplitude < .06), mas o paralelo é rejeitado; a função mostra as duas coisas",
      local({ a <- mw2$amplitude[mw2$amplitude$escala == "Verbal", ]
        a$amplitude < .06 && isTRUE(mw2$soma_rejeitada[["Verbal"]]) }))

cat("== Escore somado: verdade conhecida\n")
pop_par <- paste("f =~ .7*x1 + .7*x2 + .7*x3 + .7*x4 + .7*x5", paste0("x", 1:5, " ~~ .51*x", 1:5, collapse = "\n"), sep = "\n")
pop_tau <- paste("f =~ .7*x1 + .7*x2 + .7*x3 + .7*x4 + .7*x5", "x1 ~~ .2*x1\n x2 ~~ .4*x2\n x3 ~~ .8*x3\n x4 ~~ 1.2*x4\n x5 ~~ 1.6*x5", sep = "\n")
pop_cong <- "f =~ .9*x1 + .8*x2 + .6*x3 + .4*x4 + .3*x5"
it5 <- paste0("x", 1:5)
s_par <- checar_escore_somado(simulateData(pop_par, sample.nobs = 600, seed = 31), it5)
s_tau <- checar_escore_somado(simulateData(pop_tau, sample.nobs = 600, seed = 32), it5)
s_cong <- checar_escore_somado(simulateData(pop_cong, sample.nobs = 600, seed = 33), it5)
pd <- function(r, i) r$tabela$p[i]
teste("S-07 paralelo verdadeiro (n = 600, MLR): nenhuma das três comparações rejeita (p > .01) e o aviso de poder sai",
      all(s_par$tabela$p > ALFA_NULO) && identical(s_par$soma_rejeitada[["escala"]], FALSE) &&
        any(grepl("não mostra que elas valem", s_par$avisos)) && all(s_par$tabela$metodo == "satorra.bentler.2001"))
teste("S-08 tau-equivalente verdadeiro (cargas iguais, erros diferentes): cargas iguais não rejeitadas, erros iguais rejeitados, paralelo rejeitado",
      pd(s_tau, 1) > ALFA_NULO && pd(s_tau, 2) < .001 && pd(s_tau, 3) < .001 && isTRUE(s_tau$soma_rejeitada[["escala"]]))
teste("S-09 congenérico verdadeiro (cargas .9 a .3): paralelo rejeitado (p < .001), recomendação de SEM ou escore fatorial com Croon",
      pd(s_cong, 3) < .001 && any(grepl("correção de Croon", s_cong$decisoes)) && any(grepl("Bartlett", s_cong$decisoes)))
lim4 <- function(v) paste(sprintf("%s | -1*t1 + 0*t2 + 1*t3", v), collapse = "\n")
it5y <- paste0("y", 1:5)
o_eq <- simulateData(paste("f =~ .7*y1 + .7*y2 + .7*y3 + .7*y4 + .7*y5", lim4(it5y), sep = "\n"), sample.nobs = 800, seed = 34, ordered_center = FALSE)
o_un <- simulateData(paste("f =~ .9*y1 + .8*y2 + .6*y3 + .4*y4 + .3*y5", lim4(it5y), sep = "\n"), sample.nobs = 800, seed = 35, ordered_center = FALSE)
s_oeq <- checar_escore_somado(o_eq, it5y); s_oun <- checar_escore_somado(o_un, it5y)
teste("S-10 ordinal (4 categorias, WLSMV theta): cargas iguais verdadeiras não rejeitadas, desiguais rejeitadas; método satorra.2000; uma comparação só",
      s_oeq$ordinal && nrow(s_oeq$tabela) == 1 && s_oeq$tabela$p > ALFA_NULO && s_oun$tabela$p < .001 &&
        s_oeq$tabela$metodo == "satorra.2000" && any(grepl("extensão ao caso ordinal é da skill", s_oeq$decisoes)))
teste("S-11 dois itens: contínuo testa só as variâncias residuais (1 gl); ordinal fica sem teste, com aviso e frase que não afirma nada",
      local({ c2 <- checar_escore_somado(simulateData(pop_cong, sample.nobs = 300, seed = 36), list(A = c("x1", "x2")))
        o2 <- checar_escore_somado(o_un, list(B = c("y1", "y2")))
        nrow(c2$tabela) == 1 && c2$tabela$gl_dif == 1 && any(grepl("igualdade das cargas não tem teste", c2$avisos)) &&
          is.null(o2$tabela) && is.na(o2$soma_rejeitada[["B"]]) && any(grepl("nenhuma restrição da soma tem teste", o2$avisos)) &&
          grepl("não têm teste", o2$apa) }))
teste("S-12 método lido do cabeçalho do lavTestLRT com aspas retas ou curvas (locale UTF-8 usa dQuote)",
      .pr_metodo_lrt("Scaled Chi-Squared Difference Test (method = \"satorra.bentler.2001\")") == "satorra.bentler.2001" &&
        .pr_metodo_lrt("Scaled and Shifted Chi-Squared Difference Test (method = \u201csatorra.2000\u201d)") == "satorra.2000" &&
        .pr_metodo_lrt("\nChi-Squared Difference Test\n") == "standard")
teste("S-13 ausentes: linha sem nenhum item sai (n_excluidos), FIML mantém as parciais; item constante e escala de 1 item param com mensagem",
      local({ dm <- simulateData(pop_cong, sample.nobs = 300, seed = 37); dm$x1[1:30] <- NA; dm[31:40, it5] <- NA
        r <- checar_escore_somado(dm, it5)
        dk <- dm; dk$x5 <- 1
        e1 <- tryCatch(checar_escore_somado(dk, it5), error = function(e) conditionMessage(e))
        e2 <- tryCatch(checar_escore_somado(dm, "x1"), error = function(e) conditionMessage(e))
        r$n == 290 && r$n_excluidos == 10 && grepl("sem variação", e1) && grepl("pelo menos 2 itens", e2) }))
teste("S-14 frase em inglês com decimal = ',' e frase em português com a referência do método",
      local({ e <- checar_escore_somado(hs6, list(Verbal = c("x4", "x5", "x6")), ordinal = FALSE, idioma = "en", decimal = ",")$apa
        grepl("^For Verbal, the parallel model", e) && grepl("Satorra & Bentler, 2001", e, fixed = TRUE) && grepl("90% CI [", e, fixed = TRUE) &&
          grepl("p = ,0", e, fixed = TRUE) && grepl("Satorra & Bentler, 2001", s_par$apa, fixed = TRUE) }))
teste("S-14b itens com até 7 valores distintos viram ordinais por padrão (WLSMV, Satorra, 2000); ordinal = FALSE ou estimador da família ML volta ao contínuo",
      local({ o <- checar_escore_somado(hs6, list(Verbal = c("x4", "x5", "x6")))
        isTRUE(o$ordinal) && identical(o$estimador, "WLSMV") && grepl("Satorra, 2000", o$apa, fixed = TRUE) &&
          !isTRUE(mw1$ordinal) && identical(mw1$estimador, "ML") }))
teste("S-15 textos sem travessão nem palavras vetadas; tabelas só com agregados (nenhuma com 30 linhas ou mais)",
      local({ rs <- list(mw1, mw2, s_par, s_tau, s_cong, s_oeq, s_oun)
        txt <- unlist(lapply(rs, function(x) c(x$apa, x$decisoes, x$avisos)))
        tabs <- unlist(lapply(rs, function(x) list(x$tabela, x$ajuste, x$cargas, x$amplitude)), recursive = FALSE)
        !any(grepl(enc2utf8("\u2014|\u2013|genuinamente|de fato|notavelmente|importante notar|intricado|meticuloso|cabe ressaltar|vale destacar|essencialmente|fundamentalmente|abrangente"),
                   enc2utf8(txt), useBytes = TRUE)) && all(vapply(tabs, function(d) is.null(d) || NROW(d) < 30, TRUE)) }))
teste("S-16 tabela_ajuste com ML e FIML: versão robusta traz o χ² padrão (antes saía NA)",
      local({ f <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5", simulateData(pop_cong, sample.nobs = 300, seed = 38), estimator = "ML", missing = "fiml")
        t <- suppressWarnings(tabela_ajuste(f)); l <- t[t$versao == "robusta", ]
        nrow(l) == 1 && perto(l$qui2, unname(lavaan::fitMeasures(f, "chisq")), 1e-8) && !is.na(l$p) }))

cat("== Redes (bootnet) e EGA\n")
rede <- estimateNetwork(dados[itens], default = "EBICglasso", corMethod = "cor_auto", verbose = FALSE)
W <- abs(rede$graph); fat <- rep(1:3, each = 4)
dentro <- W[outer(fat, fat, "==") & upper.tri(W)]; entre <- W[outer(fat, fat, "!=") & upper.tri(W)]
teste("arestas dentro do fator mais fortes que entre fatores", mean(dentro) > 3 * mean(entre))
boot <- bootnet(rede, nBoots = 100, type = "case", nCores = 4, statistics = "strength", verbose = FALSE)
cs <- corStability(boot, verbose = FALSE)
teste("coeficiente CS de força calculado (0 a .75)", is.finite(cs[["strength"]]) && cs[["strength"]] >= 0 && cs[["strength"]] <= .75)
bega <- bootEGA(dados[itens], iter = 100, type = "resampling", seed = 2026, plot.itemStability = FALSE, verbose = FALSE)
teste("bootEGA recupera 3 dimensões (mediana)", bega$summary.table$median.dim == 3)
est <- itemStability(bega, IS.plot = FALSE)
teste("estabilidade de item calculada para os 12 itens", length(est$item.stability$empirical.dimensions) == 12)

cat("== Comparação de redes (NCT)\n")
nct <- NCT(g_a[g_a$g == "A", it6], g_a[g_a$g == "B", it6], it = 100, progressbar = FALSE)  # Pearson, o padrão do NCT
teste("NCT entre grupos iguais: estrutura não difere (p > .01)", nct$nwinv.pval > ALFA_NULO)

cat("== v1.4: índices do modelo bifator (Rodriguez, Reise & Haviland, 2016)\n")
# Matriz populacional como sample.cov: o ajuste reproduz as cargas populacionais e os
# índices podem ser conferidos pela fórmula, sem depender do gerador aleatório.
bf_lg <- c(.7, .65, .6, .7, .55, .6, .5, .65, .6, .55, .7, .6)
bf_ls <- c(.4, .35, .45, .3, .5, .4, .45, .35, .3, .25, .35, .4)
bf_gr <- rep(1:3, each = 4)
bf_L <- cbind(bf_lg, sapply(1:3, function(k) ifelse(bf_gr == k, bf_ls, 0)))
bf_S <- bf_L %*% t(bf_L); diag(bf_S) <- 1; dimnames(bf_S) <- list(paste0("x", 1:12), paste0("x", 1:12))
bf_mod <- paste0("G =~ ", paste0("x", 1:12, collapse = " + "), "\n",
                 paste0("S", 1:3, " =~ ", vapply(1:3, function(k) paste0("x", which(bf_gr == k), collapse = " + "), ""), collapse = "\n"))
bf_fit <- lavaan::cfa(bf_mod, sample.cov = bf_S, sample.nobs = 1000, orthogonal = TRUE, std.lv = TRUE)
bf_uni <- lavaan::cfa(paste0("G =~ ", paste0("x", 1:12, collapse = " + ")), sample.cov = bf_S, sample.nobs = 1000, std.lv = TRUE)
bf <- indices_bifator(bf_fit, unidimensional = bf_uni)
teste("bifator: ωH, ωT, ECV e PUC iguais à fórmula com as cargas populacionais (tol. 1e-5)", {
  tot <- sum(bf_S); cs <- c(sum(bf_lg), tapply(bf_ls, bf_gr, sum))
  wh <- sum(bf_lg)^2 / tot; wt <- sum(cs^2) / tot; ecv <- sum(bf_lg^2) / (sum(bf_lg^2) + sum(bf_ls^2)); puc <- 1 - 3 * 6 / 66
  v <- setNames(bf$efeito$valor, bf$efeito$indice)
  abs(v[["ômega hierárquico (ωH)"]] - wh) < 1e-5 && abs(v[["ômega total (ωT)"]] - wt) < 1e-5 &&
    abs(v[["ECV"]] - ecv) < 1e-5 && abs(v[["PUC"]] - puc) < 1e-12 })
teste("bifator: ωHS da subescala 2 = (Σλs)² / soma da submatriz; H do geral = 1 / (1 + 1/Σ[λ²/(1 - λ²)]) (tol. 1e-5)", {
  i2 <- which(bf_gr == 2); hs2 <- sum(bf_ls[i2])^2 / sum(bf_S[i2, i2]); hg <- 1 / (1 + 1 / sum(bf_lg^2 / (1 - bf_lg^2)))
  abs(bf$tabela$omega_h[bf$tabela$fator == "S2"] - hs2) < 1e-5 && abs(bf$tabela$h[bf$tabela$fator == "G"] - hg) < 1e-5 })
teste("bifator: FD = raiz de diag(Λ' Σ^-1 Λ) com fatores ortogonais (tol. 1e-5)", {
  fd <- sqrt(diag(t(bf_L) %*% solve(bf_S) %*% bf_L))
  all(abs(bf$tabela$fd - fd) < 1e-5) })
if (requireNamespace("BifactorIndicesCalculator", quietly = TRUE)) {
  teste("bifator igual ao BifactorIndicesCalculator::bifactorIndices (ECV, PUC, ω, ωH, H, FD, ARPB, I-ECV; tol. 1e-6)", {
    b <- BifactorIndicesCalculator::bifactorIndices(bf_fit, UniLambda = bf_uni, standardized = TRUE)
    fl <- b$FactorLevelIndices; ml <- b$ModelLevelIndices
    all(abs(fl[bf$tabela$fator, "Omega"] - bf$tabela$omega) < 1e-6) && all(abs(fl[bf$tabela$fator, "OmegaH"] - bf$tabela$omega_h) < 1e-6) &&
      all(abs(fl[bf$tabela$fator, "ECV_SS"] - bf$tabela$ecv) < 1e-6) && all(abs(fl[bf$tabela$fator, "ECV_SG"] - bf$tabela$ecv_total) < 1e-6) &&
      all(abs(fl[bf$tabela$fator, "H"] - bf$tabela$h) < 1e-6) && all(abs(fl[bf$tabela$fator, "FD"] - bf$tabela$fd) < 1e-6) &&
      abs(ml[["PUC"]] - bf$efeito$valor[4]) < 1e-12 && abs(ml[["ARPB"]] - bf$arpb) < 1e-6 &&
      all(abs(b$ItemLevelIndices[, "IECV"] - bf$itens$i_ecv) < 1e-6) })
} else cat("PULADO  conferência com BifactorIndicesCalculator (não instalado; as fórmulas acima cobrem os índices)\n")
teste("bifator: modelo com fatores correlacionados para com mensagem de ortogonalidade", {
  f2 <- lavaan::cfa(bf_mod, sample.cov = bf_S, sample.nobs = 1000, std.lv = TRUE)
  m <- tryCatch({ indices_bifator(f2); "" }, error = function(e) conditionMessage(e))
  grepl("ortogonais", m) })
teste("bifator: ECV .73 e PUC .73 levam à frase de unidimensionalidade essencial; H dos específicos < .70 vira aviso", {
  isTRUE(bf$essencialmente_unidimensional) && grepl("essencialmente unidimensional", bf$apa) && any(grepl("H < ", bf$avisos)) })

cat("\n== Resumo\n")
ok <- unlist(resultados)
cat(sprintf("%d de %d testes passaram em %.0f s\n", sum(ok), length(ok),
            as.numeric(difftime(Sys.time(), t_inicio, units = "secs"))))
print(versoes_pacotes())
if (!all(ok)) { cat("Falharam:", paste(names(ok)[!ok], collapse = "; "), "\n"); quit(status = 1) }
