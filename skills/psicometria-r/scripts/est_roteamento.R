# ==============================================================================
# est_roteamento.R (psicometria-r v2.0.0)
# Roteamento Algoritmico e Recomendacao de Testes Estatisticos
# Baseado na arvore de decisao do curso:
# 'A Escolha do Teste Estatistico Apropriado (Final).xlsx'
# ==============================================================================

#' Recomenda o Teste Estatistico Apropriado a partir do Delineamento e Pressupostos
#'
#' Implementa a logica da matriz de decisao do curso Psicometria-R para guiar o pesquisador
#' na selecao do teste parametrico ideal, alternativas nao-parametricas, funcao correspondente
#' na skill e justificativa metodologica detalhada.
#'
#' @param desfecho_tipo "continua" (intervalar/razao) ou "categorica" (binaria, nominal, ordinal)
#' @param n_desfechos "um" (univariado) ou "dois_ou_mais" (multivariado)
#' @param preditor_tipo "categorico", "continuo", "ambos" ou "nenhum"
#' @param n_preditores "um" ou "dois_ou_mais"
#' @param n_categorias_preditor "duas" ou "mais_de_duas" (relevante se preditor for categorico)
#' @param delineamento "entre" (grupos independentes), "pareado" (medidas repetidas) ou "misto"
#' @param normalidade TRUE se a distribuicao do desfecho/residuos for normal, FALSE caso violada
#' @param homocedasticidade TRUE se as variancias entre grupos forem homogeneas, FALSE caso heterocedasticas
#' @param objetivo Opcional: "diferenca", "associacao", "predicao", "estrutura_latente", "sobrevida", "epidemiologia"
#' @return Objeto da classe 'pr_recomendacao' contendo o teste indicado, funcao da skill e justificativa
recomendar_teste <- function(desfecho_tipo = c("continua", "categorica"),
                             n_desfechos = c("um", "dois_ou_mais"),
                             preditor_tipo = c("categorico", "continuo", "ambos", "nenhum"),
                             n_preditores = c("um", "dois_ou_mais"),
                             n_categorias_preditor = c("duas", "mais_de_duas"),
                             delineamento = c("entre", "pareado", "misto"),
                             normalidade = TRUE,
                             homocedasticidade = TRUE,
                             objetivo = NULL) {
  
  desfecho_tipo <- match.arg(desfecho_tipo)
  n_desfechos <- match.arg(n_desfechos)
  preditor_tipo <- match.arg(preditor_tipo)
  n_preditores <- match.arg(n_preditores)
  n_categorias_preditor <- match.arg(n_categorias_preditor)
  delineamento <- match.arg(delineamento)
  
  teste <- ""
  funcao <- ""
  alternativa <- ""
  justificativa <- ""
  pacote_referencia <- "stats / psicometria-r"
  
  # --- BLOCO 1: DESFECHOS MULTIPLOS (MANOVA / SEM) ---
  if (n_desfechos == "dois_ou_mais") {
    if (desfecho_tipo == "continua") {
      if (preditor_tipo == "categorico") {
        if (normalidade) {
          teste <- if (n_preditores == "um") "MANOVA (Análise de Variância Multivariada)" else "MANOVA Fatorial"
          funcao <- "manova_grupos(dados, vds, grupo)"
          alternativa <- "manova_np(dados, vds, grupo)"
          justificativa <- "Multiplas variaveis continuas dependentes avaliadas simultaneamente para controlar o erro Tipo I familiar."
        } else {
          teste <- "MANOVA Não-Paramétrica com Bootstrap de Wald (MANOVA.RM)"
          funcao <- "manova_np(dados, vds, grupo, iter = 1000)"
          alternativa <- "Testes univariados individuais com correcao de Bonferroni/Holm"
          justificativa <- "Violacao de normalidade multivariada ou homocedasticidade das matrizes de covariancia (Box M)."
        }
      } else if (preditor_tipo == "ambos") {
        teste <- "MANCOVA (Análise Multivariada de Covariância)"
        funcao <- "manova_grupos(dados, vds, grupo, covariaveis)"
        alternativa <- "Modelagem por Equações Estruturais (MEE)"
        justificativa <- "Multiplas VDs continuas com controle estatistico de covariaveis continuas."
      } else {
        teste <- "Modelagem por Equações Estruturais (MEE) ou Redes Psicométricas"
        funcao <- "mee(modelo, dados) ou rede_psicometrica(dados)"
        alternativa <- "Analise de Regressao Multivariada"
        justificativa <- "Multiplas variaveis latentes ou observadas inter-relacionadas simultaneamente."
      }
    } else {
      teste <- "Modelos de Classes Latentes (LCA) ou Redes Ising"
      funcao <- "classes_latentes(dados, itens) ou rede_psicometrica(dados, modelo = 'IsingFit')"
      alternativa <- "Regressoes Logisticas Multivariadas"
      justificativa <- "Multiplos desfechos categoricos observados."
    }
    
  # --- BLOCO 2: DESFECHO UNICO CATEGORICO (QUI-QUADRADO / LOGISTICA) ---
  } else if (desfecho_tipo == "categorica") {
    if (preditor_tipo == "categorico") {
      if (delineamento == "entre") {
        teste <- "Teste Qui-Quadrado de Independência (ou Exato de Fisher)"
        funcao <- "qui_quadrado(dados, var1, var2)"
        alternativa <- "qui_quadrado(dados, var1, var2, exato = TRUE)"
        justificativa <- "Associacao entre duas variaveis categoricas medidas em grupos independentes. Use Fisher se frequencias esperadas < 5."
      } else {
        if (n_categorias_preditor == "duas") {
          teste <- "Teste de McNemar para Dados Binários Pareados"
          funcao <- "mcnemar(dados, t1, t2)"
          alternativa <- "Equações de Estimativas Generalizadas (GEE) logísticas"
          justificativa <- "Comparacao da proporcao de respostas binarias antes e depois (medidas repetidas em 2 tempos)."
        } else {
          teste <- "Teste Q de Cochran para Medidas Repetidas Politômicas/3+ Tempos"
          funcao <- "q_cochran(dados, colunas_tempo)"
          alternativa <- "GEE binomial"
          justificativa <- "Avaliacao de mudanca em respostas dicotomicas ao longo de 3 ou mais momentos."
        }
      }
    } else if (preditor_tipo == "continuo") {
      if (n_preditores == "um") {
        teste <- "Regressão Logística Binária Simples / Curva ROC"
        funcao <- "regressao_logistica(dados, desfecho, preditor) ou curva_roc(dados, desfecho, preditor)"
        alternativa <- "Correlacao Ponto-Bisserial (.ee_assoc_binaria)"
        justificativa <- "Predicao da probabilidade de ocorrencia de um evento binario por uma variavel continua."
      } else {
        teste <- "Regressão Logística Múltipla (ou Multinomial se > 2 categorias)"
        funcao <- "regressao_logistica(dados, desfecho, preditores) ou modelo_multinomial(dados, desfecho, preditores)"
        alternativa <- "Modelos Lineares Generalizados (GLM)"
        justificativa <- "Predicao multivariada de desfecho categorico por conjunto de preditores continuos."
      }
    } else {
      teste <- "Regressão Logística com Preditores Mistos (GLM Binomial)"
      funcao <- "modelo_glm(dados, desfecho, preditores, familia = 'binomial')"
      alternativa <- "Arvores de Decisao / Random Forest"
      justificativa <- "Desfecho categorico com combinacao de preditores continuos e categoricos."
    }
    
  # --- BLOCO 3: DESFECHO UNICO CONTINUO (T / ANOVA / REGRESSAO / CORRELACAO) ---
  } else {
    # 3.1 Preditor Continuo
    if (preditor_tipo == "continuo") {
      if (n_preditores == "um") {
        if (normalidade) {
          teste <- "Correlação Linear de Pearson ou Regressão Linear Simples"
          funcao <- "correlacoes(dados, c(var1, var2), metodo = 'pearson') ou regressao_linear(dados, y, x)"
          alternativa <- "correlacoes(dados, c(var1, var2), metodo = 'spearman')"
          justificativa <- "Associacao bivariada linear entre duas variaveis continuas normalmente distribuidas."
        } else {
          teste <- "Correlação Não-Paramétrica de Spearman ou Tau de Kendall"
          funcao <- "correlacoes(dados, c(var1, var2), metodo = 'spearman') ou correlacoes(dados, c(var1, var2), metodo = 'kendall')"
          alternativa <- "Regressao robusta com bootstrapping"
          justificativa <- "Relacao monotona entre variaveis continuas quando o pressuposto de normalidade bivariada e violado."
        }
      } else {
        teste <- "Regressão Linear Múltipla (MQO)"
        funcao <- "regressao_linear(dados, y, preditores)"
        alternativa <- "Regressao Robusta HC3 / GLM Gaussiano com bootstrap"
        justificativa <- "Modelagem preditiva linear de uma variavel continua com multiplos preditores continuos."
      }
      
    # 3.2 Preditores Mistos (Continuo + Categorico)
    } else if (preditor_tipo == "ambos") {
      if (delineamento == "entre") {
        teste <- "ANCOVA (Análise de Covariância) ou Regressão Múltipla com Moderação"
        funcao <- "ancova(dados, y, grupo, covariavel) ou moderacao(dados, y, x, w)"
        alternativa <- "GLM com interacao"
        justificativa <- "Comparacao de grupos com controle estatistico de covariavel continua concomitante."
      } else {
        teste <- "Modelos Lineares Hierárquicos / Multinível (LMM) ou GEE"
        funcao <- "multinivel_linear(dados, formula, cluster)"
        alternativa <- "ANOVA Mista Split-Plot"
        justificativa <- "Delineamento longitudinal com medidas repetidas e preditores de nivel 1 e 2."
      }
      
    # 3.3 Preditor Categorico
    } else if (preditor_tipo == "categorico") {
      # 3.3.1 Um preditor categorico
      if (n_preditores == "um") {
        # Duas categorias
        if (n_categorias_preditor == "duas") {
          if (delineamento == "entre") {
            if (normalidade) {
              teste <- if (homocedasticidade) "Teste t de Student para Amostras Independentes" else "Teste t de Welch (Variâncias Desiguais)"
              funcao <- "teste_t_independente(dados, y, grupo)"
              alternativa <- "teste_mann_whitney(dados, y, grupo)"
              justificativa <- "Comparacao de medias de 2 grupos independentes. A skill aplica correcao de Welch automaticamente se heterocedastico."
            } else {
              teste <- "Teste U de Mann-Whitney (Wilcoxon Rank-Sum)"
              funcao <- "teste_mann_whitney(dados, y, grupo)"
              alternativa <- "teste_t_independente(dados, y, grupo, bootstrap = TRUE) ou modelo_glm(dados, y, grupo)"
              justificativa <- "Comparacao de postos entre 2 condicoes independentes com assimetria ou violacao de normalidade."
            }
          } else { # pareado
            if (normalidade) {
              teste <- "Teste t de Student para Amostras Pareadas"
              funcao <- "teste_t_pareado(dados, t1, t2)"
              alternativa <- "teste_wilcoxon_pareado(dados, t1, t2)"
              justificativa <- "Comparacao de medias intra-sujeitos entre 2 momentos ou condicoes pareadas normais."
            } else {
              teste <- "Teste dos Postos Sinalizados de Wilcoxon Pareado"
              funcao <- "teste_wilcoxon_pareado(dados, t1, t2)"
              alternativa <- "modelo_gee(dados, y ~ tempo, id)"
              justificativa <- "Comparacao nao-parametrica intra-sujeitos quando as diferencas nao seguem distribuicao normal."
            }
          }
        # Mais de duas categorias
        } else {
          if (delineamento == "entre") {
            if (normalidade) {
              teste <- if (homocedasticidade) "ANOVA de Um Fator (One-Way ANOVA)" else "ANOVA de Welch com Post-hoc de Games-Howell"
              funcao <- "anova_um_fator(dados, y, grupo)"
              alternativa <- "teste_kruskal_wallis(dados, y, grupo)"
              justificativa <- "Comparacao de medias entre 3 ou mais grupos independentes. Em caso de heterocedasticidade, usa Welch."
            } else {
              teste <- "Teste de Kruskal-Wallis (ANOVA por Postos)"
              funcao <- "teste_kruskal_wallis(dados, y, grupo)"
              alternativa <- "modelo_glm(dados, y, grupo) ou np_fatorial(dados, y ~ grupo)"
              justificativa <- "Comparacao nao-parametrica de 3 ou mais grupos independentes sob violacao de normalidade."
            }
          } else { # pareado
            if (normalidade) {
              teste <- "ANOVA de Medidas Repetidas de Um Fator"
              funcao <- "anova_medidas_repetidas(dados, colunas_tempo)"
              alternativa <- "teste_friedman(dados, colunas_tempo) ou modelo_gee(dados, y ~ tempo, id)"
              justificativa <- "Comparacao de medias do mesmo grupo ao longo de 3 ou mais momentos (com correcao de Greenhouse-Geisser para esfericidade)."
            } else {
              teste <- "Teste de Friedman para Medidas Repetidas"
              funcao <- "teste_friedman(dados, colunas_tempo)"
              alternativa <- "np_longitudinal(dados, formula, id)"
              justificativa <- "Comparacao nao-parametrica por postos para medidas repetidas em 3 ou mais momentos."
            }
          }
        }
      # 3.3.2 Dois ou mais preditores categoricos (Fatorial)
      } else {
        if (delineamento == "entre") {
          if (normalidade) {
            teste <- "ANOVA Fatorial Entre-Sujeitos (Two-Way / Three-Way ANOVA)"
            funcao <- "anova_fatorial(dados, y, c(fator1, fator2))"
            alternativa <- "np_fatorial(dados, y ~ fator1 * fator2)"
            justificativa <- "Avaliacao de efeitos principais e termos de interacao entre 2 ou mais fatores entre-sujeitos."
          } else {
            teste <- "ANOVA Fatorial Não-Paramétrica por Postos (rankFD)"
            funcao <- "np_fatorial(dados, y ~ fator1 * fator2)"
            alternativa <- "Regressao Linear Multipla com estimador sanduiche robusto HC3"
            justificativa <- "Delineamento fatorial entre-sujeitos sob violacao de normalidade multivariada ou homocedasticidade."
          }
        } else if (delineamento == "pareado") {
          if (normalidade) {
            teste <- "ANOVA Fatorial de Medidas Repetidas (Within-Subjects)"
            funcao <- "anova_medidas_repetidas(dados, itens_ou_tempos)"
            alternativa <- "np_longitudinal(dados, formula, id)"
            justificativa <- "Delineamento fatorial em que todos os fatores sao medidos intra-sujeitos."
          } else {
            teste <- "Delineamento Longitudinal Não-Paramétrico de Brunner-Munzel (nparLD)"
            funcao <- "np_longitudinal(dados, formula, id)"
            alternativa <- "Modelos Lineares Generalizados Misto (GLMM)"
            justificativa <- "Delineamento fatorial de medidas repetidas com postos relativos (RTE) sem suposicao de esfericidade."
          }
        } else { # misto
          teste <- "ANOVA Fatorial Mista (Split-Plot) ou Modelo Linear Misto (LMM)"
          funcao <- "anova_medidas_repetidas(dados, tempos, entre = grupo) ou multinivel_linear(dados, y ~ tempo * grupo, cluster)"
          alternativa <- "np_longitudinal(dados, formula, id, modelo = 'F1-LD-F1')"
          justificativa <- "Delineamento com pelo menos um fator entre-sujeitos e pelo menos um fator intra-sujeitos (medidas repetidas)."
        }
      }
    }
  }
  
  res <- list(
    teste_recomendado = teste,
    funcao_skill = funcao,
    alternativa_nao_parametrica = alternativa,
    justificativa = justificativa,
    parametros_informados = list(
      desfecho_tipo = desfecho_tipo,
      n_desfechos = n_desfechos,
      preditor_tipo = preditor_tipo,
      n_preditores = n_preditores,
      n_categorias_preditor = n_categorias_preditor,
      delineamento = delineamento,
      normalidade = normalidade,
      homocedasticidade = homocedasticidade,
      objetivo = objetivo
    )
  )
  class(res) <- c("pr_recomendacao", "list")
  res
}

#' Imprime o Objeto de Recomendacao Estatistica
print.pr_recomendacao <- function(x, ...) {
  cat("==============================================================================\n")
  cat(" RECOMENDAÇÃO DE ANÁLISE ESTATÍSTICA (Psicometria-R v2.0)\n")
  cat("==============================================================================\n")
  cat("  -> Teste Indicado:   ", x$teste_recomendado, "\n", sep = "")
  cat("  -> Função na Skill:  ", x$funcao_skill, "\n", sep = "")
  if (nzchar(x$alternativa_nao_parametrica)) {
    cat("  -> Alternativa:      ", x$alternativa_nao_parametrica, "\n", sep = "")
  }
  cat("  -> Justificativa:    ", x$justificativa, "\n", sep = "")
  cat("==============================================================================\n")
  invisible(x)
}
