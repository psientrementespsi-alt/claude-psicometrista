# evals_higiene.R, higiene-dados 1.5.0 (1.5: ajuste de pessoa) (pipeline v2, auditoria, Mahalanobis e atencao como familias C/IER): recodificar_likert_labels(), teste_little_mcar() e executar_pipeline_higiene()
args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else if (file.exists("screening_cier.R")) "." else "scripts"
source(file.path(pasta, "screening_cier.R")); source(file.path(pasta, "imputacao_mice.R"))
source(file.path(pasta, "pipeline_higiene.R"))
resultados <- list()
teste <- function(nome, expr) { ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome)); resultados[[nome]] <<- ok; invisible(ok) }
rot <- c("Não aconteceu comigo nessa semana", "Aconteceu comigo algumas vezes na semana",
         "Aconteceu comigo em boa parte da semana", "Aconteceu comigo na maior parte do tempo da semana")
mapa <- setNames(0:3, rot)
set.seed(1); d <- data.frame(`item 1` = sample(rot, 50, TRUE), i2 = sample(rot, 50, TRUE), check.names = FALSE); d$i2[3] <- NA
teste("sem mapa e sem numero: para com erro (nao ordena alfabeticamente)",
      inherits(tryCatch(recodificar_likert_labels(d), error = function(e) e), "error"))
r <- recodificar_likert_labels(d, mapa_niveis = mapa)
teste("com mapa: 'Nao aconteceu comigo' = 0 e 'maior parte do tempo' = 3; NA preservado; nomes preservados",
      all(r[[1]][d[[1]] == rot[1]] == 0) && all(r[[1]][d[[1]] == rot[4]] == 3) && is.na(r$i2[3]) && identical(names(r), names(d)))
teste("com mapa: espacos extras nao quebram a correspondencia", recodificar_likert_labels(paste0(" ", rot[2], "  "), mapa) == 1)
teste("rotulo fora do mapa vira NA com aviso",
      { w <- NULL; v <- withCallingHandlers(recodificar_likert_labels(c(rot[1], "Talvez"), mapa), warning = function(x) { w <<- x; invokeRestart("muffleWarning") })
        !is.null(w) && identical(v, c(0, NA)) })
teste("rotulos com numero sem mapa: extrai o numero", identical(recodificar_likert_labels(c("1 - Nunca", "3 - Sempre")), c(1, 3)))
set.seed(42); n <- 300; S <- matrix(.5, 4, 4); diag(S) <- 1; Z <- MASS::mvrnorm(n, rep(0, 4), S)
X <- as.data.frame(Z); X[cbind(c(5, 17), c(2, 3))] <- NA
lt <- suppressWarnings(teste_little_mcar(X))
teste("padroes de um so caso entram: gl = soma das variaveis observadas por padrao - p (4 + 3 + 3 - 4 = 6)", identical(as.numeric(lt$df), 6) && lt$padroes == 3)
Y <- as.data.frame(Z); Y$V2[Y$V1 > 0.8] <- NA
teste("MAR forte (V2 ausente quando V1 alto) e rejeitado", teste_little_mcar(Y)$p_valor < .001)
set.seed(7); W <- as.data.frame(MASS::mvrnorm(n, rep(0, 4), S)); W[matrix(runif(n * 4) < .08, n)] <- NA
lw <- teste_little_mcar(W)
teste("MCAR simulado nao e rejeitado a 1%", lw$p_valor > .01)
if (requireNamespace("naniar", quietly = TRUE)) teste("igual ao naniar::mcar_test", abs(lw$qui_quadrado - naniar::mcar_test(W)$statistic) < .01) else cat("PULADO  igual ao naniar::mcar_test (naniar ausente)\n")

# Testes da versao 1.1.2: executar_pipeline_higiene() com protecao a protocolos de base
df_pipe_teste <- data.frame(
  id = 1:6,
  idade = c(20, 22, 25, 20, 30, 40),
  i1 = c(0, 0, 0, 1, 2, 3),
  i2 = c(0, 0, 0, 2, 1, 2),
  i3 = c(0, 0, 0, 3, 2, 1)
)
# Linhas 1 e 2 tem respostas 0 nos itens, mas idades diferentes (20 e 22). Linha 3 tem idade 25 e zeros nos itens.
# Se executar com ignorar_todos_zeros = TRUE e range_valido_itens = c(0, 3), nenhuma das 3 deve ser descartada como duplicata!
res_pipe <- suppressMessages(executar_pipeline_higiene(
  df_bruto = df_pipe_teste,
  cols_itens = c("i1", "i2", "i3"),
  col_id = "id",
  range_valido_itens = c(0, 3),
  ignorar_todos_zeros = TRUE
))

teste("pipeline protege protocolos de todos zeros de exclusao espuria por duplicata",
      res_pipe$log_descartes$duplicatas == 0 && nrow(res_pipe$dados_purificados) == 6)

# Teste com duplicata real em idade e itens
df_dup_real <- rbind(df_pipe_teste, df_pipe_teste[4, ])
res_dup <- suppressMessages(executar_pipeline_higiene(
  df_bruto = df_dup_real,
  cols_itens = c("i1", "i2", "i3"),
  col_id = "id",
  range_valido_itens = c(0, 3),
  ignorar_todos_zeros = TRUE
))
teste("pipeline detecta duplicatas reais em variaveis substantivas",
      res_dup$log_descartes$duplicatas == 1 && nrow(res_dup$dados_purificados) == 6)


cat("== v1.2.0: even-odd, triagem por familias e MCD75\n")
set.seed(11); n <- 200; th <- rnorm(n)
Y <- sapply(1:20, function(j) cut(th + rnorm(n, sd = .9), c(-Inf, -1, 0, 1, Inf), labels = FALSE) - 1)
colnames(Y) <- paste0("y", 1:20)
CR <- Sys.getenv("PR_CARELESS_R", "~/Desktop/psicomtria rrr/repositorios-higiene-dados/04_triagem_comportamental_e_careless/careless/R")
if (dir.exists(CR)) {
  ec <- new.env(); sys.source(file.path(CR, "evenodd.R"), envir = ec)
  ref_eo <- suppressWarnings(ec$evenodd(as.data.frame(Y), c(5, 5, 5, 5)))
  teste("even-odd igual ao careless::evenodd (subescalas contiguas, tol. 1e-10)", isTRUE(all.equal(unname(calc_evenodd(Y, c(5, 5, 5, 5))), unname(ref_eo), tolerance = 1e-10)))
} else cat("PULADO  even-odd contra o careless (repositorio ausente)\n")
ord <- sample(20); fat <- split(colnames(Y)[ord], rep(1:4, each = 5))
teste("even-odd aceita subescalas nao contiguas por nome (mesmo valor que reordenar as colunas)",
      isTRUE(all.equal(unname(calc_evenodd(Y, fat)), unname(calc_evenodd(Y[, unlist(fat)], c(5, 5, 5, 5))))))
Z <- Y; Z[1:4, ] <- 2; Z[5:6, ] <- 0
tri <- suppressMessages(triagem_cier_completa(Z, ignorar_todos_zeros = TRUE))
teste("repeticao e uma familia so: quem repete 2 sem outro indice nao e excluido", all(tri$familia_repeticao[1:4]) && !any(tri$excluir_por_cier[1:4]))
teste("protocolos com zero em todos os itens nao entram na familia repeticao", !any(tri$familia_repeticao[5:6]))
tri2 <- suppressMessages(triagem_cier_completa(Z, tempo_total_seg = c(rep(20, 4), rep(200, n - 4)), ignorar_todos_zeros = TRUE))
teste("repeticao e tempo < 2 s por item juntos excluem (2 familias)", all(tri2$excluir_por_cier[1:4]) && sum(tri2$excluir_por_cier) == 4)
teste("corte do longstring = metade da escala (Curran, 2016): 10 de 20 itens", grepl(">= 10 de 20", attr(tri, "decisoes")[1]))
avs <- character(); invisible(withCallingHandlers(triagem_cier_completa(Z, usar_psychsyn = TRUE, limite_informativo = 0), message = function(m) { avs <<- c(avs, conditionMessage(m)); invokeRestart("muffleMessage") }))
teste("indice que assinala acima do limite_informativo fica fora do criterio, com aviso", any(grepl("nao discrimina", avs)))
if (requireNamespace("robustbase", quietly = TRUE)) {
  set.seed(3); W <- MASS::mvrnorm(150, c(0, 0, 0), diag(3)); W[1:5, ] <- W[1:5, ] + 6
  mc <- calc_mahad_robust(W)
  teste("FastMCD: MCD75 por padrao (Leys et al., 2019) e os 5 atipicos plantados assinalados", grepl("alpha = 0.75", mc$estimador) && all(mc$is_outlier[1:5]))
} else cat("PULADO  FastMCD (robustbase ausente)\n")

# ---- 1.3.0 ----
set.seed(21); n <- 240; f1 <- rnorm(n); f2 <- rnorm(n)
I <- as.data.frame(cbind(sapply(1:6, function(j) pmin(4, pmax(0, round(2 + f1 + rnorm(n, 0, .7))))),
                         sapply(1:6, function(j) pmin(4, pmax(0, round(2 + f2 + rnorm(n, 0, .7)))))))
names(I) <- paste0("q", 1:12); I[1:5, ] <- matrix(0:4, 5, 12)
I[6:10, ] <- rbind(rep(c(0, 4), 6), rep(c(4, 0), 6), rep(c(0, 3), 6), rep(c(3, 0), 6), rep(c(1, 4), 6))
tmp <- c(rep(8, 10), runif(n - 10, 80, 400))
mm <- calc_mahad_cier(I)
teste("Mahalanobis C/IER: alternancia 0/4 e 4/0 (casos 6 e 7) sai com p < .001 contra qui-quadrado(12); retilineos nao", all(mm$flag_mahad[6:7]) && !any(mm$flag_mahad[1:5]))
mf <- calc_mahad_cier(I, fatores = list(1:6, 7:12))
teste("Mahalanobis por fator (Meade & Craig, 2012) com Bonferroni devolve p entre 0 e 1", all(mf$mahad_p >= 0 & mf$mahad_p <= 1, na.rm = TRUE))
att <- data.frame(a1 = c(rep(1, 5), rep(3, n - 5)), a2 = c(rep(5, 3), rep(1, n - 3)))
fa <- calc_falhas_atencao(att, esperado = c(3, NA), incorreto = list(NULL, 5))
teste("itens de atencao: instrucao direta e infrequencia somam falhas", all(fa[1:3] == 2) && all(fa[4:5] == 1) && all(fa[6:n] == 0))
tc <- triagem_cier_completa(I, tempo_total_seg = tmp, n_itens_total = 12, usar_psychsyn = FALSE, usar_mahad = TRUE)
teste("triagem: Mahalanobis + tempo exclui os alternados; repeticao + tempo exclui os retilineos", all(tc$excluir_por_cier[1:7]) && sum(tc$excluir_por_cier[-(1:10)]) <= 2)
dv <- data.frame(idade = c(17, 30, NA, 45, 120), sexo = c("F", "M", "F", NA, "F"))
vr <- validar_regras(dv, c(adulto = "idade >= 18", plausivel = "idade <= 99", sexo_valido = "sexo %in% c('F','M')"))
teste("validar_regras: contagens de violacao e de nao avaliaveis, sem linhas do banco", vr$violacoes[vr$regra == "adulto"] == 1 && vr$nao_avaliaveis[vr$regra == "adulto"] == 1 &&
  vr$violacoes[vr$regra == "plausivel"] == 1 && vr$violacoes[vr$regra == "sexo_valido"] == 1 && !any(c("idade", "sexo") %in% names(vr)))
cb <- gerar_codebook(data.frame(g = c(rep("a", 10), "raro"), txt = paste("resposta", 1:11), x = c(1:10, NA)))
teste("codebook: categoria com menos de 5 casos suprimida e texto livre omitido", grepl("raro: <5", cb$resumo[cb$variavel == "g"]) &&
  cb$tipo[cb$variavel == "txt"] == "texto livre" && !grepl("resposta 3", paste(cb$resumo, collapse = " ")) && cb$n_ausentes[cb$variavel == "x"] == 1)
dfp <- data.frame(I, tempo = tmp, consentiu = c(FALSE, rep(TRUE, n - 1)))
dfp[n + 1, ] <- NA; dfp[n + 2, ] <- dfp[100, ]; dfp$q3[50] <- 9
rp <- executar_pipeline_higiene(dfp, paste0("q", 1:12), col_tempo = "tempo", range_valido_itens = c(0, 4),
  fatores = list(A = paste0("q", 1:6), B = paste0("q", 7:12)), regras = c(consentimento = "consentiu == TRUE"),
  regras_excluir = "consentimento", verboso = FALSE)
tr <- rp$trilha
teste("pipeline v2: trilha encadeada (n_depois de cada etapa = n_antes da seguinte) e N final = linhas do banco",
  all(utils::head(tr$n_depois, -1) == tr$n_antes[-1]) && utils::tail(tr$n_depois, 1) == nrow(rp$dados_purificados))
teste("pipeline v2: vazio, criterio de inclusao, duplicata e valor fora da escala registrados",
  tr$n_removido[tr$etapa == "Registros vazios"] == 1 && tr$n_removido[tr$etapa == "Criterios de inclusao"] == 1 &&
  tr$n_removido[tr$etapa == "Duplicatas"] == 1 && rp$log_descartes$valores_fora_faixa == 1)
teste("pipeline v2: Tabela 1 em Markdown e texto de Metodo saem com os mesmos numeros da trilha",
  any(endsWith(rp$tabela1_md, sprintf("| %d |", utils::tail(tr$n_depois, 1)))) && grepl(sprintf("Dos %d registros", tr$n_antes[1]), rp$texto_metodo))
teste("fluxograma_mermaid: sintaxe valida gerada a partir da trilha e presente no resultado do pipeline",
  is.character(rp$fluxograma_mermaid) && any(grepl("^graph TD", rp$fluxograma_mermaid)) && any(grepl("Amostra analítica final", rp$fluxograma_mermaid)))
if (requireNamespace("mice", quietly = TRUE)) {
  df_incompleto <- data.frame(a = c(1, 2, NA, 4, 2, 6, 1, 8, 3, 5), b = c(5, 3, 4, NA, 1, 7, 2, 9, 8, 4), c = c(2, 4, 3, 1, 5, 2, 4, 3, 1, 5))
  imp_teste <- suppressWarnings(mice::mice(df_incompleto, m = 3, maxit = 2, printFlag = FALSE))
  exp_mi <- exportar_lavaan_mi(imp_teste)
  teste("exportar_lavaan_mi: gera lista com m data.frames completos prontos para semTools/lavaan.mi",
    is.list(exp_mi) && length(exp_mi) == 3 && all(vapply(exp_mi, function(x) !anyNA(x) && nrow(x) == 10, logical(1))))
} else cat("PULADO  exportar_lavaan_mi (pacote mice nao instalado)\n")
pasta_est <- file.path(tempdir(), "estado_eval"); unlink(pasta_est, recursive = TRUE)
st <- executar_pipeline_higiene(dfp, paste0("q", 1:12), col_tempo = "tempo", range_valido_itens = c(0, 4), fatores = list(A = paste0("q", 1:6), B = paste0("q", 7:12)),
  regras = c(consentimento = "consentiu == TRUE"), regras_excluir = "consentimento", salvar_em = pasta_est, etapas_ate = "omissao", verboso = FALSE)
rr <- executar_pipeline_higiene(NULL, NULL, retomar_de = file.path(pasta_est, "estado_omissao.rds"), verboso = FALSE)
teste("pipeline v2 em blocos: parar em 'omissao' e retomar da o mesmo resultado que rodar direto",
  st$ultima == "omissao" && identical(rr$trilha$n_depois, rp$trilha$n_depois) && identical(rr$mcd$n_assinalados, rp$mcd$n_assinalados))
F_DASS <- Sys.getenv("PR_DASS_CSV", "")
if (file.exists(F_DASS) && requireNamespace("robustbase", quietly = TRUE)) {
  d <- read.csv(F_DASS, fileEncoding = "UTF-8", check.names = FALSE, colClasses = "character", na.strings = c("", "NA"))
  it <- sprintf("dass_%02d", 1:21); raw <- d[, 72:92]; names(raw) <- it
  sb <- list(depressao = c(3,5,10,13,16,17,21), ansiedade = c(2,4,7,9,15,19,20), estresse = c(1,6,8,11,12,14,18))
  outras <- d[, -(72:92)]; outras <- outras[, !grepl("^Survey Timestamp|^Complete\\?", names(outras))][, -1]
  rd <- executar_pipeline_higiene(data.frame(outras, raw, consentiu = grepl("^Concordo", d[[grep("^Após ler o termo", names(d))[1]]]), check.names = FALSE), it,
    range_valido_itens = c(0, 3), mapa_niveis_likert = mapa, fatores = lapply(sb, function(v) sprintf("dass_%02d", v)),
    regras = c(consentimento = "consentiu == TRUE"), regras_excluir = "consentimento", verboso = FALSE, gerar_codebook_saida = FALSE)
  teste("DASS-21 do projeto: pipeline v2 reproduz N analitico = 116 e 5 atipicos MCD75 nas subescalas (fase 1 de 29/10/2026)",
    nrow(rd$dados_purificados) == 116 && rd$mcd$n_assinalados == 5)
} else cat("PULADO  reproducao do banco DASS-21 (arquivo ausente)\n")
# ---- v1.5.0: ajuste de pessoa (Niessen, Meijer & Tendeiro, 2016) ----
source(file.path(pasta, "ajuste_pessoa.R"))
gA <- calc_guttman_poly(data.frame(N23 = c(1, 4), N74 = c(2, 3), N81 = c(1, 1)), categorias = 1:4,
  p_passos = c(.49, .58, .05, .85, .77, .12, .53, .58, .12))
teste("G^p reproduz o Apendice A de Niessen et al. (2016): respondentes 12 e 161 com 0 e 6 erros, passos na ordem da Tabela A2",
  identical(as.numeric(gA), c(0, 6)) && identical(attr(gA, "ordem_passos"), c("N74-1", "N74-2", "N23-2", "N81-2", "N81-1", "N23-1", "N74-3", "N81-3", "N23-3")))
gD <- calc_guttman_poly(rbind(c(1, 1, 1, 0, 0), c(1, 0, 0, 1, 1)), 1:5, categorias = 0:1, p_passos = c(.9, .7, .5, .3, .1))
teste("G^p dicotomico: [1,1,1,0,0] tem 0 erros e [1,0,0,1,1] tem 4 (exemplo da secao 3.1)", identical(as.numeric(gD), c(0, 4)))
set.seed(11); nn <- 400; th <- rnorm(nn); aa <- runif(8, 1, 2.2); bb <- t(apply(matrix(rnorm(24), 8), 1, sort))
sim <- sapply(1:8, function(j) { pp <- sapply(1:3, function(k) plogis(aa[j] * (th - bb[j, k]))); rowSums(pp > runif(nn)) })
sim[1:40, ] <- matrix(sample(0:3, 320, TRUE), 40); colnames(sim) <- paste0("q", 1:8); Dp <- as.data.frame(sim)
if (requireNamespace("mirt", quietly = TRUE)) {
  fitp <- mirt::mirt(Dp, 1, "graded", verbose = FALSE); thp <- mirt::fscores(fitp, verbose = FALSE)[, 1]
  teste("lz^p igual ao Zh do mirt::personfit() com o mesmo theta EAP (dif. < 1e-8)",
    max(abs(calc_lz_poly(Dp) - mirt::personfit(fitp, Theta = matrix(thp))$Zh)) < 1e-8)
} else cat("PULADO  lz^p contra mirt::personfit (mirt ausente)\n")
ap <- suppressMessages(calc_ajuste_pessoa(Dp, list(esc = names(Dp)), calcular_lz = FALSE))
teste("ajuste de pessoa: com 10% de respondentes aleatorios, a maioria dos assinalados pelo G^p_m e aleatoria",
  mean(which(ap$flag_ajuste_pessoa) <= 40) > .5 && attr(ap, "resumo")$n_assinalados <= ceiling(.06 * nn))
tp <- suppressMessages(triagem_cier_completa(Dp, tempo_total_seg = c(rep(5, 40), rep(120, nn - 40)), n_itens_total = 8,
  usar_psychsyn = FALSE, flag_ajuste_pessoa = ap))
teste("triagem: ajuste de pessoa conta como familia; so exclui junto com outra (tempo)",
  all(tp$excluir_por_cier == (tp$flag_ajuste_pessoa + tp$flag_tempo_rapido + tp$familia_repeticao >= 2)) && any(tp$flag_ajuste_pessoa & !tp$excluir_por_cier))
ev <- new.env(); suppressMessages(sys.source(file.path(pasta, "validacao_codebook.R"), envir = ev))
teste("validacao_codebook.R nao sobrescreve mais validar_regras() com outra assinatura", identical(names(formals(validar_regras)), c("df", "regras")))
cat(sprintf("\nTotal higiene-dados 1.5.0: %d | Passaram: %d | Falharam: %d\n", length(resultados), sum(unlist(resultados)), sum(!unlist(resultados))))

