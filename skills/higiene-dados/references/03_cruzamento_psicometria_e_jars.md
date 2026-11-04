# Cruzamento Linha a Linha: Corpus de Higiene com a Skill `psicometria-r`

## Sumário

- 7 (linha 12)
- 84 (linha 89)

**Finalidade:** Mapear a intersecção exata entre os 30 textos do corpus de higiene e os arquivos estruturais da skill `psicometria-r` (`references/*.md` e scripts `est_*.R`), fornecendo a Claude o blueprint completo para a futura implementação, **sem alterar nenhum arquivo da skill nesta sessão**.

---

## 1. Mapeamento por Módulo e Arquivo de `references/` da Skill

### 1.1 `references/fontes-e-decisoes.md` (Regras de Decisão, Pressupostos e Integridade)
- **Status Atual na v1.3:** A skill estabelece regras estritas sobre dados agregados, relato de ICs e pressupostos como gradiente de indícios (não teste rígido de significância).
- **Lacuna Existente:** Falta um protocolo formal e codificado para a etapa pré-analítica de triagem de respondentes desatentos (*data screening*) e uma regra declarada para a escolha entre FIML e Imputação Múltipla.
- **Obras do Corpus que Fundamentam a Expansão:**
  - *Aguinis et al. (2013)* (`10.1177/1094428112470848`): A árvore de decisão em 14 etapas para diagnosticar e relatar *outliers* com dupla modelagem (com e sem o caso).
  - *DeSimone, Harms, & DeSimone (2015)* (`10.1002/job.1962`): A regra prática de triagem multinível.
  - *Ward & Meade (2023)* (`10.1146/annurev-psych-040422-045007`): A regra de ouro de que **nenhum participante é excluído por um único critério atencional isolado**; exigem-se no mínimo duas violações independentes.
  - *van der Loo & de Jonge (2018)* (`10.1002/9781118897126`): A separação entre regras de validação (DSL do `validate`) e a manipulação dos dados.
- **Regras Sugeridas para Claude Incorporar:**
  1. *Regra Pré-Analítica de Triagem:* Todo banco de dados comportamental deve passar por triagem de desatenção antes da AFE ou AFC. O log de participantes removidos deve relatar o motivo exato (speeding, longstring ou consistência) e o percentual da amostra total.
  2. *Regra de Relato de Outliers:* Nenhum outlier com valor biologicamente plausível é expurgado sem análise de sensibilidade que compare as cargas fatoriais e correlações com e sem ele.

---

### 1.2 `references/construcao-instrumentos.md` (Construção, Itens e AFE)
- **Status Atual na v1.3:** Cobre análise clássica de itens, dificuldade, discriminação e retenção fatorial na AFE.
- **Lacuna Existente:** Não possui diretrizes para limpeza de respostas invariantes (*straightlining*) em matrizes de Likert nem para detecção de tempos excessivamente rápidos por página.
- **Obras do Corpus que Fundamentam a Expansão:**
  - *Curran (2016)* (`10.1016/j.jesp.2015.07.006`): Índices de sinônimos e antônimos psicométricos pós-hoc usando pares de itens com $|r| > .60$.
  - *Meade & Craig (2012)* (`10.1037/a0028085`): Métrica de invariabilidade (*longstring*) e consistência par-ímpar.
  - *Huang et al. (2015)* (`10.1007/s10869-014-9357-6`): Uso de escalas de infrequência no blueprint do teste.
  - *Boehmke (2016)* (`10.1007/978-3-319-45599-0`): Tratamento de itens com variância zero ou quasi-zero antes da fatorabilidade (KMO e Bartlett).
- **Blueprint para o Script `est_construcao.R`:**
  - Adição da função `triagem_desatentos()` acionável antes de `afe_adequacao()`, calculando `careless::longstring()` e `careless::irv()`.

---

### 1.3 `references/afc-ordinal.md` (AFC Categórica, Estimadores e Dados Ausentes)
- **Status Atual na v1.3:** Discute estimadores WLSMV vs MLR e menciona dados ausentes, mas sem aprofundar os pressupostos de estimação.
- **Lacuna Existente:** Falta a distinção formal entre FIML (ótimo para MLR com indicadores contínuos) e Imputação Múltipla com PMM/polr via `mice` (necessária para indicadores categóricos/ordinais quando se quer evitar o pressuposto MCAR do pairwise do WLSMV).
- **Obras do Corpus que Fundamentam a Expansão:**
  - *Schafer & Graham (2002)* (`10.1037/1082-989X.7.2.147`) e *Graham (2009)* (`10.1146/annurev.psych.58.110405.085530`): Teoria dos mecanismos MCAR, MAR e MNAR e uso de variáveis auxiliares.
  - *Lang & Little (2018)* (`10.1007/s11121-016-0644-5`): Árvore de escolha formal entre FIML e MI.
  - *van Buuren & Groothuis-Oudshoorn (2011)* (`10.18637/jss.v045.i03`): Implementação prática de FCS via pacote `mice`.
  - *Tierney & Cook (2023)* (`10.18637/jss.v105.i07`): Diagnóstico exploratório de ausência via pacote `naniar` e teste de MCAR de Little (`mcar_test()`).
  - *Wiley Handbook (2018, Cap. 6 e Cap. 9)*: Alerta técnico de que a estimação pairwise com matriz policórica (WLSMV) assume dados omissos completamente aleatórios (MCAR) nos pares bivariados; se o mecanismo for MAR, a Imputação Múltipla via `mice` antes da AFC ordinal é a abordagem defensável.
- **Blueprint para `psicometria.R` (`ajustar_afc`):**
  - Implementar verificação automática com `naniar::mcar_test()`. Se o teste falhar e houver indicadores categóricos ordinais com missingness significativo ($> 5\%$), sugerir rota com `mice` antes do `lavaan`.

---

### 1.4 `references/associacao.md` (Correlações, Outliers e Curva ROC)
- **Status Atual na v1.3:** Possui a função `curva_roc` e correlações lineares, mencionando outliers de passagem.
- **Lacuna Existente:** Não possui implementação nativa do MAD para diagnóstico robusto univariado nem de Mahalanobis-MCD para correlações bivariadas e multivariadas.
- **Obras do Corpus que Fundamentam a Expansão:**
  - *Leys et al. (2013)* (`10.1016/j.jesp.2013.03.013`): A fórmula clássica de $MAD = 1.4826 \cdot \text{mediana}(|X_i - \text{mediana}(X)|)$ e o ponto de corte ótimo de 2,5 MAD.
  - *Leys et al. (2019)* (`10.5334/irsp.289`): Mahalanobis robusto com MCD.
  - *Hubert, Debruyne, & Rousseeuw (2018)* (`10.1002/wics.1421`): Algoritmo FastMCD via pacote `robustbase::covMcd`.
  - *Petersen (2024, Cap. 3 e Cap. 10)*: Impacto de outliers de alavancagem em coeficientes de calibração clínica e regressão logística diagnóstica.
- **Blueprint para `est_associacao.R`:**
  - Criação da função auxiliar `detectar_outliers_mad()` e inclusão do argumento `outlier_mad = TRUE` em análises correlacionais e descritivas da skill.

---

### 1.5 `references/invariancia.md` (Invariância de Medida e Equidade)
- **Status Atual na v1.3:** Trata da sequência de invariância configural, métrica, escalar e estrita.
- **Lacuna Existente:** Não adverte sobre o impacto da contaminação diferencial: respondentes desatentos ou dados ausentes não tratados em apenas um dos subgrupos (ex: grupo online vs grupo presencial) criam não-invariância métrica e escalar espúria.
- **Obras do Corpus que Fundamentam a Expansão:**
  - *Chmielewski & Kucker (2020)* (`10.1177/1948550619875149`): Taxas diferenciais de desatenção entre amostras universitárias e plataformas online (MTurk/Prolific) quebram a equivalência de medida.
  - *Cook & Pitoniak (2026, Cap. 6)*: Equidade e padronização das condições de resposta como pré-requisito da invariância.
  - *Petersen (2024, Cap. 15)*: Diferenciação entre viés de medida e viés de resposta por desatenção.

---

### 1.6 `references/relato-apa.md` (Tabelas e Frases Padrão APA 7)
- **Status Atual na v1.3:** Fornece frases para Método e Resultados de testes inferenciais e psicométricos.
- **Lacuna Existente:** Falta a frase padrão recomendada pelas diretrizes metodológicas da APA 7 (JARS-Quant) para a seção de Participantes e Triagem Pré-Analítica.
- **Obras do Corpus que Fundamentam a Expansão:**
  - *Lang & Little (2018)* e *Aguinis et al. (2013)*: Frases padronizadas de relato de dados omissos e outliers.
  - *Ward & Meade (2023)*: Padrão de relatório de exclusão de respostas desatentas.
- **Template APA 7 para Claude Adicionar:**
  > *"Antes das análises primárias, os dados foram submetidos a triagem atencional sistemática (Ward & Meade, 2023). Participantes foram sinalizados por velocidade de resposta excessivamente rápida (< 1,5 segundos/item), invariabilidade excessiva de resposta (longstring) e inconsistência intra-individual em sinônimos psicométricos. Seguindo o critério conservador de requerer no mínimo duas violações independentes, N = X participantes (Y%) foram excluídos. O diagnóstico de dados ausentes indicou mecanismo [MCAR/MAR], tratado via [FIML no lavaan / Imputação Múltipla com m = 20 no pacote mice]. Casos atípicos univariados foram avaliados pelo Desvio Absoluto Mediano (MAD > 2,5; Leys et al., 2013), com análises de sensibilidade demonstrando estabilidade das estimativas substantivas."*

---

## 2. Blueprint de Funções Prontas para a Implementação do Claude

Quando o Claude for integrar este material à skill `psicometria-r`, as seguintes funções representam a estrutura ideal:

```r
# ==============================================================================
# BLUEPRINT PARA CLAUDE: Módulo est_higiene.R (ou expansão de est_comum.R)
# ==============================================================================

# 1. Triagem Comportamental Integrada
triagem_desatentos <- function(dados, colunas_itens, col_tempo = NULL, 
                                corte_longstring = 10, corte_tempo_item = 1.5) {
  # Implementação baseada em careless (Curran 2016; Ward & Meade 2023)
  # Retorna flags e data frame limpo
}

# 2. Diagnóstico Robusto de Outliers
diagnosticar_outliers_mad <- function(x, corte = 2.5) {
  # Implementação baseada em Leys et al. (2013)
  # Retorna vetor lógico e limites empíricos
}

# 3. Diagnóstico e Teste de Ausentes
diagnosticar_ausentes <- function(dados, colunas) {
  # Implementação baseada em naniar (Tierney & Cook 2023; Schafer & Graham 2002)
  # Executa mcar_test() e sumariza padrões
}

# 4. Auditoria e Validação Lógica
validar_consistencia <- function(dados, regras_validator) {
  # Implementação baseada em validate (van der Loo & de Jonge 2018, 2021)
  # Confronta regras e devolve tabela de violações
}
```
