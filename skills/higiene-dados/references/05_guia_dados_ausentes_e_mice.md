# Guia de Diagnóstico e Tratamento de Dados Ausentes

Este guia estabelece os procedimentos teóricos e computacionais para avaliar mecanismos de perda e aplicar métodos modernos de imputação múltipla em pesquisas em psicologia.

---

## 1. Mecanismos de Omissão (Rubin, 1976; Schafer & Graham, 2002)

1. **MCAR (Missing Completely at Random):**
   - A probabilidade de um dado ser ausente é estritamente independente de qualquer variável observada ou não observada no estudo.
   - *Implicação:* Análise de casos completos (*listwise deletion*) é não viesada, embora cause perda severa de poder estatístico.
   - *Teste:* Teste de MCAR de Little (1988). Sob $H_0$, os dados são MCAR. Se $p > .05$, o mecanismo MCAR é plausível.

2. **MAR (Missing at Random):**
   - A probabilidade de omissão depende de variáveis observadas no banco de dados (ex.: idade, escolaridade, traço prévio), mas não dos valores hipotéticos da própria variável ausente após condicionamento nas variáveis observadas.
   - *Implicação:* *Listwise deletion* gera estimativas severamente enviesadas.
   - *Tratamento adequado:* Full Information Maximum Likelihood (FIML) ou Imputação Múltipla por Equações Encadeadas (MICE).

3. **MNAR (Missing Not at Random):**
   - A probabilidade de omissão depende diretamente do valor da própria variável omitida (ex.: participantes deprimidos graves que deixam de responder itens de ideação suicida).
   - *Tratamento:* Exige modelos de seleção ou *pattern-mixture models*. MICE e FIML podem ser complementados com análise de sensibilidade com deslocamento de parâmetros (*tipping point analysis*).

---

## 2. Diagnóstico Prático de Ausentes

1. **Inspeção de Padrões:**
   - Visualização da matriz de omissão (ex.: via `naniar::vis_miss()` ou tabelas de coocorrência).
   - Identificação de padrões de atrito (*attrition*) em desenhos longitudinais.
2. **Cálculo da Proporção de Perda:**
   - Taxa por variável: itens com $> 30\%$ a $40\%$ de ausência devem ser investigados criticamente antes de entrar no modelo.
   - Taxa por participante: indivíduos com mais de $20\%$ de omissão na bateria central geralmente são excluídos da amostra analítica com justificativa documentada.
3. **Variáveis Auxiliares (Graham, 2009):**
   - Variáveis que não fazem parte do modelo substantivo de interesse, mas que estão correlacionadas com a probabilidade de omissão ou com as variáveis com dados ausentes ($r > .30$).
   - A inclusão de variáveis auxiliares amplia a sustentabilidade do pressuposto MAR, reduzindo o viés de estimação.

---

## 3. Imputação Múltipla com MICE (van Buuren, 2018)

### 3.1 Fundamentos do MICE
O MICE (*Multivariate Imputation by Chained Equations*, também conhecido como FCS - *Fully Conditional Specification*) imputa dados ausentes variável por variável, especificando um modelo condicional univariado para cada variável em função de todas as outras.

### 3.2 Escolha dos Métodos de Imputação
- **Variáveis Contínuas / Escores:** Regressão linear bayesiana (`norm`), correspondência preditiva média (`pmm`). O método `pmm` preserva a distribuição original dos dados e evita imputações fora dos limites reais.
- **Variáveis Binárias:** Regressão logística (`logreg`).
- **Variáveis Ordinais (Itens Likert):** Regressão logística ordinal com chances proporcionais (`polr`).
- **Variáveis Nominais Poliatômicas:** Regressão logística multinomial (`polyreg`).

### 3.3 Número de Imputações ($m$)
A recomendação histórica de $m = 3$ a $5$ conjuntos de dados é inadequada. A prática contemporânea recomenda:
$$m \ge 20 \quad \text{a} \quad 50$$
ou $m \approx 100 \times \text{proporção máxima de omissão}$ (ex.: $25\%$ de omissão $\implies m \ge 25$).

### 3.4 Regras de Rubin para Agrupamento (*Pooling*)
As estimativas pontuais dos $m$ modelos imputados são agregadas pela média simples:
$$\bar{Q} = \frac{1}{m} \sum_{l=1}^m \hat{Q}^{(l)}$$
A variância total combina a variância intra-imputação ($\bar{U}$) e a variância entre-imputações ($B$):
$$T = \bar{U} + \left(1 + \frac{1}{m}\right) B$$

---

## 4. Articulação com a Psicometria no R

- **Para AFC com Itens Ordinais (WLSMV / DWLS):** A estimação sobre matriz policórica no `lavaan` com exclusão pairwise assume formalmente **MCAR**. Se os dados violarem MCAR, a rota recomendada é imputar os itens ordinais via MICE (`polr` ou `pmm`) em $m$ bancos, ajustar a AFC em cada banco via WLSMV e agrupar com `semTools::runMI()`.
- **Para AFC Contínua ou Modelos de Regressão:** Usar FIML com estimador robusto (`estimator = "MLR"`, `missing = "fiml"`), que é eficiente e consistente sob **MAR**.
