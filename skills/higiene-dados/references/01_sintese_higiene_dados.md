# Manual de Síntese Sistemática: Higiene, Triagem e Tratamento de Dados em R

## Sumário

- 8 (linha 18)
- 44 (linha 54)
- 105 (linha 115)
- 175 (linha 185)
- 240 (linha 250)
- 285 (linha 295)
- 303 (linha 313)

**Corpus Base:** 30 obras fundamentais (3 livros integrais de referência em R, 26 artigos metodológicos peer-reviewed e 1 resenha técnica) cruzadas com os 6 tratados de psicometria da pasta `livros-completos-md/`.  
**Finalidade:** Fornecer a base teórica, as regras formais, os pontos de corte e os scripts reproduzíveis em R para a etapa pré-analítica, deixando o material 100% preparado para posterior adaptação pelo Claude na skill `psicometria-r`.

---

## 1. Visão Geral da Arquitetura Pré-Analítica

A literatura contemporânea em psicometria e ciência de dados (Wickham, 2014; van der Loo & de Jonge, 2018; Ward & Meade, 2023; Aguinis et al., 2013; Schafer & Graham, 2002) converge em um princípio intransigente: **nenhuma análise inferencial ou psicométrica (AFE, AFC, TRI, invariância, redes) deve ser rodada sobre dados brutos não higienizados**.

A negligência na triagem pré-analítica gera quatro distorções severas:
1. **Inflação espúria de correlações e fatores artificiais:** Respondentes desatentos (*careless/insufficient effort responding - C/IER*) produzem fatores espúrios de método ou de aquiescência, distorcendo índices de ajuste global ($\chi^2$, CFI, RMSEA) e gerando dimensionalidade espúria (Curran, 2016; Ward & Meade, 2023).
2. **Distorção de estimadores por contaminação:** A presença de 5% de *outliers* não detectados reduz a eficiência da média aritmética e do desvio-padrão a zero, enquanto correlações de Pearson colapsam ou se tornam falsamente significativas (Leys et al., 2013, 2019; Rousseeuw & Hubert, 2011).
3. **Erros de Tipo I e Tipo II por tratamentos ingênuos de ausência:** A exclusão por lista (*listwise deletion*) induz viés sistemático nas cargas fatoriais caso o mecanismo não seja MCAR (*Missing Completely at Random*), além de descartar poder estatístico precioso (Schafer & Graham, 2002; Lang & Little, 2018).
4. **Falta de reprodutibilidade e corrupção de variáveis:** Modificações pontuais em planilhas sem rastreamento automático (*audit trail*) impossibilitam a auditoria científica e violam as diretrizes de integridade da APA 7 e FAIR (Broman & Woo, 2018; van der Loo, 2021; Peng, 2011).

O fluxo rigoroso de higiene divide-se em 5 fases sequenciais:

```
[Dados Brutos]
      │
      ▼
[Fase 1: Estruturação Tidy e Validação Técnica] (Wickham 2014; Broman & Woo 2018; validate)
      │
      ▼
[Fase 2: Triagem Comportamental / Careless Responding] (Ward & Meade 2023; careless; IRV, Longstring)
      │
      ▼
[Fase 3: Diagnóstico e Manejo de Outliers] (Leys et al. 2013, 2019; MAD; Mahalanobis-MCD)
      │
      ▼
[Fase 4: Diagnóstico e Tratamento de Dados Ausentes] (Schafer & Graham 2002; naniar; FIML / mice)
      │
      ▼
[Fase 5: Documentação Automatizada e Rastreamento] (Arslan 2019; codebook; lumberjack)
      │
      ▼
[Banco Limpo Pronto para Modelagem Psicométrica]
```

---

## 2. Fase 1: Arquitetura de Dados Tabulares e Validação Formal

### 2.1 Princípios de Tidy Data e Planilhas Estruturadas
Conforme estabelecido por **Wickham (2014)** (`10.18637/jss.v059.i10`) e **Broman & Woo (2018)** (`10.1080/00031305.2017.1375989`):
1. Cada variável deve residir em sua própria coluna.
2. Cada observação (participante) deve ocupar exatamente uma linha.
3. Cada unidade observacional deve formar uma tabela dedicada.
4. Nomes de variáveis devem usar caracteres alfanuméricos minúsculos sem acentos, com separadores em *snake_case* (ex: `item_01`, `tempo_total_seg`).
5. Células mescladas, linhas em branco, cores como significado de dado e cabeçalhos em múltiplas linhas são expressamente proibidos em dados brutos.

### 2.2 Validação Declarativa com o Pacote `validate`
**Van der Loo & de Jonge (2018, 2021)** (`10.1002/9781118897126`, `10.18637/jss.v097.i10`) estruturaram a validação formal em R via linguagem declarativa específica de domínio (DSL). As regras de consistência lógica e limites numéricos devem ser separadas do código de execução:

```r
library(validate)

# Definição formal das regras de integridade psicométrica
regras_integridade <- validator(
  # Definição do grupo de itens para validação em lote
  G := var_group(item_01, item_02, item_03, item_04, item_05,
                 item_06, item_07, item_08, item_09, item_10,
                 item_11, item_12, item_13, item_14, item_15,
                 item_16, item_17, item_18, item_19, item_20),
  # Limites de escala Likert (ex: 1 a 5)
  itens_intervalo = in_range(G, min = 1, max = 5),
  
  # Tipagem estrita
  id_unico = is_unique(id_participante),
  idade_valida = idade >= 18 & idade <= 100,
  
  # Consistência lógica de tempos (ex: tempo mínimo plausível de leitura)
  tempo_minimo = tempo_total_seg >= 120,
  
  # Ausência de pontuações impossíveis
  nao_negativo = tempo_total_seg > 0
)

# Confronto dos dados contra as regras
confronto <- confront(dados_brutos, regras_integridade)
summary(confronto)
barplot(confronto, main = "Violações de Regras de Integridade")
```

### 2.3 Rastreamento Automático com o Pacote `lumberjack`
**Van der Loo (2021)** (`10.18637/jss.v098.i01`) demonstrou que a limpeza de dados em R deve gerar registros automáticos de mutação sem intervenção manual:

```r
library(lumberjack)

# Registro de mutações durante o pipeline de limpeza
dados_limpos <- dados_brutos %L>%
  start_log(log = cellwise$new(key = "id_participante")) %L>%
  # Passo 1: Impor NAs em códigos de erro (ex: 99, 999)
  mutate(across(starts_with("item_"), ~ ifelse(. %in% c(99, 999), NA, .))) %L>%
  # Passo 2: Conversão de tipos
  mutate(sexo = as.factor(sexo)) %L>%
  dump_log(file = "auditoria_limpeza.csv", stop = TRUE)
```

---

## 3. Fase 2: Triagem Comportamental (*Careless Responding / Screening*)

A triagem comportamental visa diagnosticar participantes que responderam com esforço insuficiente (*Insufficient Effort Responding - IER*), gerando dados corrompidos por desatenção, fazendas de cliques ou bots (Ward & Meade, 2023; Curran, 2016; Meade & Craig, 2012; Chmielewski & Kucker, 2020).

### 3.1 Métricas de Triagem Pós-Hoc
Segundo a revisão de **Ward & Meade (2023)** (`10.1146/annurev-psych-040422-045007`) e os estudos de **Curran (2016)** e **Meade & Craig (2012)**:

1. **Invariabilidade de Resposta (*Longstring*):**
   - *Conceito:* Contagem máxima de respostas idênticas dadas consecutivamente ao longo da escala (ex: marcar `4` em 15 itens consecutivos).
   - *Cálculo:* `careless::longstring(dados_itens)`.
   - *Ponto de corte:* Avaliado contra a distribuição empírica ou simulação de respostas aleatórias (regra prática: mais de 6 a 10 respostas consecutivas idênticas em escalas com itens invertidos).
2. **Variância Intra-Individual de Resposta (*IRV - Intra-individual Response Variability*):**
   - *Conceito:* O desvio-padrão das respostas de um indivíduo através de uma bateria de itens.
   - *Fórmula:* $IRV_i = \sqrt{\frac{1}{K-1} \sum_{j=1}^{K} (X_{ij} - \bar{X}_i)^2}$, onde $K$ é o número de itens.
   - *Interpretação:* $IRV \approx 0$ indica respondentes invariantes (*straightliners*); valores excessivamente altos associados a tempos rápidos indicam cliques aleatórios alternados.
3. **Sinônimos e Antônimos Psicométricos:**
   - *Conceito:* Correlação intra-individual calculada sobre pares de itens com correlação de amostra extremamente alta ($r > .60$, sinônimos) ou negativa ($r < -.60$, antônimos).
   - *Interpretação:* O respondente atencioso deve pontuar de forma consistente nesses pares. Correlações intra-individuais negativas em sinônimos psicométricos revelam forte desatenção (Curran, 2016; careless::psychsyn).
4. **Consistência Par-Ímpar (*Even-Odd Consistency*):**
   - *Conceito:* Correlação dividida das metades de uma mesma subescala ajustada pela fórmula de Spearman-Brown para cada indivíduo.
5. **Velocidade de Resposta (*Speeders*):**
   - *Ponto de corte de tempo:* DeSimone, Harms, & DeSimone (2015) e Ward & Meade (2023) estabelecem que respostas com menos de **1,5 a 2,0 segundos por item** (ou menos de 1 segundo para itens muito curtos) são fisicamente incompatíveis com a leitura e processamento cognitivo do enunciado.

### 3.2 Script de Screening no R via Pacote `careless`

```r
library(careless)
library(dplyr)

# 1. Longstring
ls_scores <- longstring(dados_itens)

# 2. Variância Intra-Individual (IRV)
irv_scores <- irv(dados_itens)

# 3. Sinônimos Psicométricos (itens com r > .60; Meade & Craig, 2012)
syn_scores <- psychsyn(dados_itens, critval = 0.60)

# 4. Antônimos Psicométricos (itens com r < -0.60; Curran, 2016)
ant_scores <- psychant(dados_itens, critval = -0.60)

# 5. Distância de Mahalanobis Multivariada (D²) baseada no perfil de itens (sem plot interativo)
mah_scores <- mahad(dados_itens, plot = FALSE, flag = FALSE)

# Matriz Consolidada de Triagem
triagem <- tibble(
  id = dados$id,
  tempo_item = dados$tempo_total_seg / ncol(dados_itens),
  longstring = ls_scores,
  irv = irv_scores,
  syn = syn_scores,
  ant = ant_scores,
  d2 = mah_scores
)

# Regra de Exclusão Conservadora (Ward & Meade, 2023: nunca excluir por apenas 1 critério isolado):
# Excluir apenas participantes que falharem em pelo menos DOIS critérios independentes:
triagem <- triagem %>%
  mutate(
    flag_speed = tempo_item < 1.5,
    flag_long = longstring > 10,  # ajustar conforme tamanho da escala
    flag_irv = irv < quantile(irv, 0.05, na.rm = TRUE),
    flag_syn = syn < 0.0,
    total_flags = flag_speed + flag_long + flag_irv + flag_syn,
    excluir_careless = total_flags >= 2
  )
```

---

## 4. Fase 3: Diagnóstico e Manejo de Casos Atípicos (*Outliers*)

### 4.1 O Erro da Regra da Média $\pm$ 2 ou 3 Desvios-Padrão
**Leys, Ley, Klein, Bernard, & Licata (2013)** (`10.1016/j.jesp.2013.03.013`) e **Leys et al. (2019)** (`10.5334/irsp.289`) demonstraram matematicamente por que a regra clássica de excluir valores além de 2 ou 3 desvios-padrão da média é ineficaz:
1. A média aritmética e o desvio-padrão não são resistentes a anomalias: o próprio caso extremo infla a média e multiplica o desvio-padrão (*mascaramento*).
2. O ponto de corte afasta-se proporcionalmente ao tamanho do caso extremo, falhando em detectar outliers óbvios em amostras pequenas ou moderadas.

### 4.2 O Desvio Absoluto Mediano (MAD)
A alternativa univariada formal e robusta é o **Median Absolute Deviation (MAD)**:

$$\text{MAD} = b \cdot \text{mediana}(|X_i - \text{mediana}(X)|)$$

Onde $b = 1.4826$ é a constante de correção assintótica para equivalência ao desvio-padrão sob distribuição normal.

- **Critérios de Corte Recomendados na Literatura:**
  - **Conservador (Miller, 1991; Leys et al., 2013):** $[\text{Mediana} - 3 \cdot \text{MAD}, \text{Mediana} + 3 \cdot \text{MAD}]$.
  - **Moderado / Padrão Recomendado (Leys et al., 2013):** $[\text{Mediana} - 2.5 \cdot \text{MAD}, \text{Mediana} + 2.5 \cdot \text{MAD}]$.
  - **Nota Metodológica sobre o Corte 2.24:** O limiar de 2.24 unidades provém da regra paramétrica de desvios-padrão de Aguinis et al. (2013, p. 280; citando Martin & Roberts, 2010) para delimitar 2.5% em cada cauda sob normalidade ($p < .025$ bicaudal), e não deve ser confundido com os multiplicadores clássicos de MAD de Miller (1991) e Leys et al. (2013), que são 2.5 e 3.0.

### 4.3 Casos Atípicos Multivariados: Mahalanobis-MCD
Para variáveis combinadas ou escores de subescalas, a distância clássica de Mahalanobis sofre do mesmo problema de mascaramento porque a matriz de covariância clássica é sensível a pontos de alavancagem (Rousseeuw & Hubert, 2011; Hubert et al., 2018).

A solução é a **Distância de Mahalanobis com Determinante de Covariância Mínima (FastMCD)** implementada via `robustbase::covMcd`:

$$RD(X_i) = \sqrt{(X_i - \hat{\mu}_{\text{MCD}})^T \hat{\Sigma}_{\text{MCD}}^{-1} (X_i - \hat{\mu}_{\text{MCD}})}$$

- *Ponto de corte:* Valores de $RD^2 > \chi^2_{p, 0.999}$ ou corte ajustado de Green & Martin.

### 4.4 Script R para Detecção de Outliers Univariados e Multivariados

```r
library(robustbase)

# 1. Função Univariada via MAD (Leys et al., 2013)
detectar_outliers_mad <- function(x, corte = 2.5) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, constant = 1.4826, na.rm = TRUE)
  lim_inf <- med - corte * mad_val
  lim_sup <- med + corte * mad_val
  return(x < lim_inf | x > lim_sup)
}

# 2. Função Multivariada via FastMCD (Rousseeuw & Hubert, 2011; Hubert et al., 2018)
detectar_outliers_mcd <- function(dados_matriz, alpha = 0.75, conf = 0.999) {
  mcd_res <- covMcd(dados_matriz, alpha = alpha)
  d2_robusto <- mahalanobis(dados_matriz, center = mcd_res$center, cov = mcd_res$cov)
  corte_chi2 <- qchisq(conf, df = ncol(dados_matriz))
  return(list(
    outliers = d2_robusto > corte_chi2,
    distancias = sqrt(d2_robusto),
    corte = sqrt(corte_chi2)
  ))
}
```

### 4.5 Árvore de Decisão de Aguinis et al. (2013)
Conforme **Aguinis, Gottfredson, & Joo (2013)** (`10.1177/1094428112470848`):
1. **Erro de Registro / Codificação:** Se for erro identificável (ex: digitação de 55 em escala de 1 a 5), corrige-se pelo protocolo original ou converte-se para `NA`.
2. **Outlier Interessante / Legítimo da População:** Se o caso reflete um indivíduo genuíno com traços extremos:
   - **NÃO DELETAR cegamente.**
   - Rodar a análise primária com o caso incluído.
   - Rodar análise de sensibilidade com método robusto (ou sem o caso) e relatar **ambos** no manuscrito, declarando o impacto nas estimativas de parâmetros e nas conclusões substantivas.

---

## 5. Fase 4: Diagnóstico e Tratamento de Dados Ausentes (*Missing Data*)

### 5.1 Taxonomia dos Mecanismos (Schafer & Graham, 2002; Graham, 2009)
1. **MCAR (Missing Completely at Random):** A ausência é estocasticamente independente de variáveis observadas e não observadas ($P(M|Y_{\text{obs}}, Y_{\text{mis}}) = P(M)$).
   - *Diagnóstico:* Teste de MCAR de Little (`naniar::mcar_test()`). Se $p > .05$, não há evidência contra MCAR.
2. **MAR (Missing at Random):** A ausência depende de variáveis observadas, mas não do valor omisso em si ($P(M|Y_{\text{obs}}, Y_{\text{mis}}) = P(M|Y_{\text{obs}})$).
   - *Tratamento exigido:* FIML (*Full Information Maximum Likelihood*) ou Imputação Múltipla (MI). Ambos produzem estimativas não enviesadas sob MAR.
3. **MNAR (Missing Not at Random):** A probabilidade de ausência depende do próprio valor omisso (ex: participantes deprimidos deixam de responder perguntas sobre humor disfórico). Exige modelagem de seleção de Heckman ou modelos de classes latentes com ausência informada.

### 5.2 Diagnóstico Visual com o Pacote `naniar`
**Tierney & Cook (2023)** (`10.18637/jss.v105.i07`):

```r
library(naniar)

# Mapeamento do percentual de ausência por item e por participante
vis_miss(dados_brutos)
gg_miss_upset(dados_brutos)

# Teste formal de Little
teste_little <- mcar_test(dados_brutos)
print(teste_little)
```

### 5.3 Escolha Metodológica: FIML vs Imputação Múltipla (`mice`)
Conforme a síntese de **Lang & Little (2018)** (`10.1007/s11121-016-0644-5`) e **van Buuren (2011, 2018)**:

| Critério | FIML (ex: `lavaan::cfa(..., missing = "fiml")`) | Imputação Múltipla (ex: `mice`) |
|---|---|---|
| **Ambiente de Modelagem** | AFC, SEM, Modelos de Trajetória, Regressão Latente. | Qualquer modelo (não lineares, regressões múltiplas, AFE). |
| **Itens Ordinais / Categóricos** | Limitado sob WLSMV (requer pairwise com pressuposto MCAR; MLR com integração numérica fica lento com $>5$ fatores). | **Ideal via FCS:** permite modelos de imputação politômicos ordenados (`polr` ou `pmm`) respeitando categorias discretas. |
| **Variáveis Auxiliares** | Devem ser integradas diretamente na sintaxe do modelo SEM. | Podem ser incluídas na matriz de preditores do `mice` sem alterar o modelo substantivo final. |
| **Padrão de Pooling** | Parâmetros e erros-padrão estimados diretamente em etapa única. | $m$ bancos imputados (regra prática: $m \ge \text{percentual de ausência}$, ex: $m = 20$), combinados pelas Regras de Rubin. |

### 5.4 Imputação Não Paramétrica via `missForest`
Para bancos de dados complexos com interações e relações não lineares, **Stekhoven & Bühlmann (2012)** (`10.1093/bioinformatics/btr597`) forneceram o algoritmo `missForest`, que não assume normalidade nem linearidade e performa imputação iterativa baseada em florestas aleatórias:

```r
library(missForest)
imputacao_rf <- missForest(dados_mistos, maxiter = 10, ntree = 100)
dados_completos <- imputacao_rf$ximp
```

---

## 6. Fase 5: Documentação Automatizada, Dicionários e Reprodutibilidade

Conforme **Arslan (2019)** (`10.1177/2515245919838783`):
O código de preparação de dados deve culminar na criação automática de um dicionário formal de dados (*codebook*) contendo metadados das escalas, confiabilidades alfa/ômega preliminares, distribuições visuais e percentual de valores omissos:

```r
library(codebook)

# Atribuição de rótulos e metadados formais
var_label(dados$item_01) <- "Sinto-me energizado pela manhã."
val_labels(dados$item_01) <- c("Discordo totalmente" = 1, "Discordo" = 2, "Neutro" = 3, "Concordo" = 4, "Concordo totalmente" = 5)

# Geração automática do livro de códigos interativo
codebook(dados)
```

---

## 7. Roteiro Prático de Execução para a Skill `psicometria-r`

Este referencial fornece a Claude a base exata para a criação das seguintes funções modulares na skill:

1. `triagem_comportamental(dados, colunas_itens, tempo_col = NULL)`:
   - Extrai `longstring`, `irv`, `psychsyn`, `psychant` e calcula tempo médio por item.
   - Aplica a regra conservadora de exclusão (mínimo de 2 violações).
2. `validar_regras(dados, arquivo_regras_yaml)`:
   - Confronta os dados com a gramática declarativa do pacote `validate` e gera sumário de violações.
3. `diagnosticar_outliers(dados, metodo = c("mad", "mcd", "ambos"), corte = 2.5)`:
   - Retorna flags de casos extremos pelo MAD (univariado) e Mahalanobis-MCD (multivariado), gerando gráfico diagnóstico.
4. `diagnosticar_ausentes(dados, colunas)`:
   - Roda o teste de Little (`naniar::mcar_test`), mapeia percentual por variável e gera recomendação formal: FIML (se AFC contínua) ou MICE (se itens ordinais/AFE).
5. `gerar_dicionario(dados, output_html = "codebook.html")`:
   - Documenta automaticamente o banco limpo via pacote `codebook`.
