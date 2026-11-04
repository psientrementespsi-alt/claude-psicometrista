# MANUAL SISTEMÁTICO DE RELATO ESTATÍSTICO SEGUNDO APA 7 E JARS-QUANT

## Sumário

- 2 (linha 10)
- 16 (linha 24)
- 34 (linha 42)
- 50 (linha 58)
- 149 (linha 157)
## NORMAS DE ESCRITA, TABELAS DE PUBLICAÇÃO E CONVENÇÕES METODOLÓGICAS

**Autor/Organizador:** Antigravity AI Orchestrator  
**Destino:** `~/Desktop/psicomtria rrr/referencial-relato-estatistico/`  
**Data:** 29 de Outubro de 2026  
**Corpus de Fundamentação:**
- *Publication Manual of the American Psychological Association* (APA 7th ed., 2020)
- *Journal Article Reporting Standards for Quantitative Research in Psychology* (JARS-Quant; Appelbaum et al., 2018; DOI: `10.1037/amp0000191`)
- *Statistical Methods in Psychology Journals: Guidelines and Explanations* (Wilkinson & Task Force, 1999; DOI: `10.1037/0003-066X.54.8.594`)
- *Reporting Effect Sizes in Original Psychological Research* (Pek & Flora, 2018; DOI: `10.1037/met0000126`)
- *The Prevalence of Statistical Reporting Errors in Psychology* (Nuijten et al., 2016; DOI: `10.3758/s13428-015-0664-2`)

---

## 1. PRINCÍPIOS FUNDAMENTAIS E REGRAS DE OURO DA APA 7

### 1.1 Formatação Tipográfica de Símbolos Estatísticos
1. **Uso Obrigatório de Itálico:** Todas as letras latinas usadas como símbolos estatísticos devem ser grafadas em itálico ($N$, $n$, $M$, $SD$, $Mdn$, $t$, $F$, $r$, $p$, $z$, $b$, $d$).
2. **Sem Itálico em Letras Gregas e Subscritos:** Letras gregas ($\alpha$, $\beta$, $\omega$, $\chi^2$, $\eta^2$, $\Delta$) **nunca** levam itálico. Subscritos e sobrescritos textuais também não levam itálico, a menos que correspondam a variáveis estatísticas (ex: $\eta_p^2$, $\omega_h$, $R^2$, mas $S_{\text{diff}}$).
3. **A Regra do Zero à Esquerda (Leading Zero):**
   - **Omitir o zero inicial** para estatísticas que **não podem ultrapassar 1.0**: proporções, probabilidades, correlações e coeficientes de confiabilidade (escrever $.05$, $.85$, $.001$, e **nunca** $0.05$ ou $0.85$). Isso se aplica estritamente a: $p$, $r$, $\rho$, $\tau$, $\beta$ padronizado, $\alpha$, $\omega$, $\text{CFI}$, $\text{TLI}$.
   - **Manter o zero inicial** para estatísticas que **podem teoricamente exceder 1.0**: médias, desvios-padrão, valores de teste, distâncias e tamanhos de efeito não limitados (escrever $0.45$, $1.20$, $0.08$). Isso se aplica a: $M$, $SD$, $t$, $F$, $\chi^2$, $z$, $d$ de Cohen, $\text{RMSEA}$, $\text{SRMR}$, $B$ não padronizado, $OR$ (Odds Ratio).
4. **Casas Decimais Padronizadas:**
   - **Duas casas decimais:** Para quase todas as estatísticas descritivas e inferenciais ($M = 24.50$, $SD = 3.12$, $t(98) = 2.45$, $F(2, 85) = 4.10$, $d = 0.52$).
   - **Três casas decimais:** Para valores de $p$ ($p = .024$, $p = .008$) e índices de ajuste que exigem alta precisão em psicometria ($\text{RMSEA} = 0.048$, $\text{CFI} = .962$).
   - **Valores Extremos de $p$:** Nunca escrever $p = .000$. O menor valor reportável é $p < .001$. Quando o valor for marginal ou limítrofe, reportar o número exato com três casas (ex: $p = .049$ ou $p = .051$), e nunca relatar apenas "$p < .05$" ou "$p > .05$".
5. **Espaçamento e Pontuação:**
   - Inserir espaço antes e depois de operadores matemáticos ($=$, $<$, $>$, $+$, $-$). Exemplo correto: $t(48) = 2.31, p = .025$; incorreto: $t(48)=2.31, p=.025$.
   - O número de graus de liberdade deve ser colocado entre parênteses imediatamente após a letra do teste, sem espaço: $t(98)$, $F(2, 45)$, $\chi^2(12)$.

---

## 2. REQUISITOS DO JARS-QUANT (APA TASK FORCE)

Segundo o relatório oficial do APA Publications and Communications Board Task Force (Appelbaum et al., 2018), todo manuscrito empírico quantitativo deve relatar:

1. **Poder Amostral e Determinação do $N$:** Relatar na subseção de Participantes ou Delineamento se o tamanho amostral foi planejado *a priori* via cálculo de poder estatístico (ex: G*Power, simulação de Monte Carlo ou precisão de RMSEA), especificando o tamanho de efeito mínimo clinicamente relevante, o erro Tipo I ($\alpha = .05$) e o poder desejado ($1 - \beta \ge .80$).
2. **Fluxo de Participantes e Dados Omissos:**
   - Descrever o número total de participantes recrutados, o número de desistentes e o número exato de casos excluídos com a justificativa técnica explícita para cada exclusão (critérios de integridade, tempo de resposta excessivamente rápido, invariância de respostas).
   - Indicar a porcentagem de dados ausentes por variável e globalmente.
   - Especificar o teste de mecanismo de ausência realizado (ex: teste de MCAR de Little) e o método formal de tratamento adotado (Imputação Múltipla com $m$ replicações ou FIML com MLR). A exclusão por lista (*listwise*) só é admitida se comprovadamente MCAR e com perda inferior a 5%.
3. **Tamanhos de Efeito e Intervalos de Confiança (Mandatórios):**
   - Wilkinson & Task Force (1999) e Pek & Flora (2018) estabelecem que **o relato de $p$-valor isolado é cientificamente inaceitável**.
   - Toda estimativa pontual (médias, diferenças, coeficientes de regressão, cargas fatoriais) deve ser acompanhada de uma medida de tamanho de efeito apropriada e seu respectivo **Intervalo de Confiança de 95%** ($95\% \text{ CI } [LL, UL]$).
   - Se o intervalo de confiança for assimétrico (como em efeitos indiretos mediados via bootstrap), explicitar que se trata de bootstrap percentil ou BCa com o número de reamostragens ($B \ge 5.000$).

---

## 3. CONVENÇÕES DE RELATO POR TIPO DE ANÁLISE

### 3.1 Estatística Descritiva e Caracterização Demográfica
- **Variáveis Contínuas:** Relatar Média ($M$) e Desvio-Padrão ($SD$). Se a distribuição for marcadamente assimétrica, relatar Mediana ($Mdn$) e Intervalo Interquartílico ($IQR$).
- **Variáveis Categóricas:** Relatar contagem absoluta ($n$) e porcentagem relativa ($\%$).
- **Proteção de Privacidade:** Em tabelas cruzadas, qualquer célula com contagem $n < 5$ deve ter sua contagem suprimida ou agrupada para evitar reidentificação de participantes.
- **Modelo de Frase em Português:**
  > "A amostra foi composta por 320 participantes adultos ($M = 28.45$ anos, $SD = 7.12$, variando de 18 a 65 anos), sendo a maioria do sexo feminino ($n = 210$, $65.6\%$) e solteira ($n = 185$, $57.8\%$). O tempo mediano de resposta foi de 14.50 minutos ($IQR = 6.20$)."
- **Modelo de Frase em Inglês:**
  > "The sample comprised 320 adult participants ($M = 28.45$ years, $SD = 7.12$, range: 18–65), predominantly female ($n = 210$, $65.6\%$) and unmarried ($n = 185$, $57.8\%$). The median completion time was 14.50 minutes ($IQR = 6.20$)."

---

### 3.2 Testes de Comparação de Médias ($t$ de Student e Welch)
- **Especificação:** Indicar se foi adotada a correção de Welch para variâncias desiguais ($t$ de Welch com graus de liberdade corrigidos/fracionários).
- **Tamanho de Efeito:** Relatar $d$ de Cohen ou $g$ de Hedges (para amostras pequenas $N < 20$ ou desiguais) e o IC de 95% do tamanho de efeito.
- **Modelo de Frase (Welch):**
  > "O grupo intervenção apresentou escores significativamente maiores de bem-estar ($M = 42.10$, $SD = 5.30$) em comparação ao grupo controle ($M = 35.40$, $SD = 7.80$), $t(82.41) = 4.82$, $p < .001$, $d = 1.01$, $95\% \text{ CI } [0.57, 1.45]$."
- **Equivalente Não Paramétrico (Mann-Whitney):**
  > "Uma comparação via teste $U$ de Mann-Whitney indicou que o grupo clínico reportou maiores níveis de ansiedade ($Mdn = 18.00$, $IQR = 6.00$) do que o grupo não clínico ($Mdn = 10.00$, $IQR = 4.00$), $U = 412.00$, $p = .003$, correlação bisserial de postos $r_{\text{rb}} = .42$, $95\% \text{ CI } [.16, .63]$."

---

### 3.3 Análise de Variância (ANOVA)
- **Especificação:** Indicar o teste de homogeneidade de variâncias (Levene). Em caso de violação, relatar a ANOVA de Welch ($F_{\text{Welch}}$).
- **Tamanho de Efeito:** Relatar $\eta_p^2$ (eta ao quadrado parcial) ou $\omega^2$ (ômega ao quadrado, menos enviesado).
- **Comparações Post-Hoc:** Indicar o método de controle do erro familiar (Tukey HSD para grupos equilibrados, Games-Howell para variâncias heterogêneas, Bonferroni ou Holm).
- **Modelo de Frase:**
  > "Uma ANOVA univariada de uma via revelou efeito principal significativo da condição experimental sobre o desempenho, $F(2, 147) = 8.45$, $p < .001$, $\eta_p^2 = .10$, $95\% \text{ CI } [.03, .18]$. Comparações post-hoc com correção de Tukey indicaram que a condição A ($M = 78.20$, $SD = 9.10$) superou significativamente a condição B ($M = 70.10$, $SD = 8.50$, $p = .002$, $d = 0.92$) e o controle ($M = 68.40$, $SD = 10.20$, $p < .001$, $d = 1.01$), sem diferença significativa entre B e o controle ($p = .620$)."

---

### 3.4 Correlações Lineares e Não Paramétricas
- **Especificação:** Relatar o coeficiente ($r$, $\rho$ ou $\tau$), o número de graus de liberdade ($N - 2$), o valor de $p$ e o IC de 95%.
- **Modelo de Frase:**
  > "Observou-se correlação positiva moderada e estatisticamente significativa entre suporte social e autoestima, $r(198) = .42$, $p < .001$, $95\% \text{ CI } [.30, .53]$. A correlação parcial controlando para os efeitos de idade permaneceu robusta, $r_{\text{p}}(197) = .39$, $p < .001$."

---

### 3.5 Regressão Linear Múltipla e Regressão Logística
- **Regressão Linear:** Relatar o ajuste global do modelo ($R^2$, $R^2_{\text{adj}}$, $F(k, N - k - 1)$, $p$) e, para cada preditor, o coeficiente não padronizado com seu erro-padrão ($B$, $SE$), o coeficiente padronizado ($\beta$), a estatística $t$ e o IC de 95% de $B$.
- **Regressão Logística:** Relatar a razão de chances ($OR$ ou $\text{Exp}(B)$), seu IC de 95%, o teste de Wald ($z$), o teste de Hosmer-Lemeshow e os pseudo-$R^2$ (Cox & Snell, Nagelkerke).
- **Modelo de Frase (Linear):**
  > "O modelo de regressão linear explicou $34.5\%$ da variância total da satisfação com a vida, $R^2 = .35$, $R^2_{\text{adj}} = .33$, $F(3, 196) = 34.82$, $p < .001$. Neuroticismo associou-se negativamente à satisfação ($B = -0.38$, $SE = 0.06$, $\beta = -.41$, $t = -6.12$, $p < .001$, $95\% \text{ CI } [-0.50, -0.26]$), enquanto extroversão exibiu associação positiva ($B = 0.24$, $SE = 0.05$, $\beta = .30$, $t = 4.55$, $p < .001$, $95\% \text{ CI } [0.14, 0.34]$)."
- **Modelo de Frase (Logística):**
  > "A probabilidade de adesão ao tratamento foi significativamente prevista pela aliança terapêutica, $B = 0.65$, $SE = 0.18$, Wald $z = 3.61$, $p < .001$, $OR = 1.92$, $95\% \text{ CI } [1.35, 2.73]$, indicando que cada incremento de um ponto na escala dobra praticamente as chances de conclusão da terapia."

---

### 3.6 Mediação e Moderação (Modelos de Processo / PROCESS / SEM)
- **Mediação:** Relatar os efeitos direto ($c'$), indireto ($a \times b$) e total ($c$). **Nunca usar o teste de Sobel**, pois ele assume normalidade do produto de dois coeficientes (violação matemática grave). O relato canônico exige **Intervalo de Confiança com Bootstrap Percentil ou BCa** (mínimo de 5.000 replicações).
- **Moderação:** Relatar a interação produto, a mudança na variância explicada ($\Delta R^2$, $F$-change) e a análise de declives simples (*simple slopes*) nos níveis $-1 SD$, Média e $+1 SD$ do moderador (ou via técnica de Johnson-Neyman).
- **Modelo de Frase (Mediação):**
  > "A análise de mediação indicou efeito total significativo do estresse precoce sobre a depressão adulta, $c = 0.45$, $SE = 0.08$, $t = 5.62$, $p < .001$, $95\% \text{ CI } [0.29, 0.61]$. O efeito indireto transmitido pela desregulação emocional foi estatisticamente significativo, $ab = 0.18$, $SE_{\text{boot}} = 0.04$, $95\% \text{ CI}_{\text{boot}} [0.10, 0.27]$ ($5.000$ reamostragens percentis). Controlando-se para o mediador, o efeito direto permaneceu significativo, porém atenuado, $c' = 0.27$, $SE = 0.08$, $t = 3.38$, $p = .001$, $95\% \text{ CI } [0.11, 0.43]$, caracterizando uma mediação parcial."

---

### 3.7 Modelagem Psicométrica (AFE, AFC e MEE)
- **AFE (Análise Fatorial Exploratória):**
  - Informar o teste KMO e a esfericidade de Bartlett.
  - Justificar o método de retenção fatorial (Análise Paralela com matriz policórica ou MRFA; não usar a regra de Kaiser dos autovalores $> 1$ isoladamente).
  - Indicar o estimador (WLSMV / ULS para ordinais, MLR para contínuos) e a rotação oblíqua (oblimin, geomin, promax).
  - Relatar a tabela com cargas fatoriais (omitindo cargas $< .30$ ou $< .40$), comunalidades ($h^2$) e correlações entre fatores.
- **AFC (Análise Fatorial Confirmatória):**
  - Especificar o software e pacote (`lavaan` no R), o estimador (`WLSMV` com matriz policórica para itens ordinais; `MLR` com `missing = "ML"` sob MAR).
  - Relatar a bateria global de índices de ajuste com seus pontos de corte recomendados:
    - $\chi^2$ escalado (Satorra-Bentler / Yuan-Bentler) com graus de liberdade e $p$-valor.
    - $\text{CFI} \ge .95$ (aceitável $\ge .90$)
    - $\text{TLI} \ge .95$ (aceitável $\ge .90$)
    - $\text{RMSEA} \le .06$ com Intervalo de Confiança de 90% ($\text{90\% CI } [LL, UL]$) e probabilidade de ajuste estrito ($p_{\text{close}} > .05$)
    - $\text{SRMR} \le .08$ (para contínuos) ou $\text{WRMR} \le 1.00$
- **Confiabilidade Composta:**
  - Relatar o Ômega Total de McDonald ($\omega_t$) e o Ômega Categórico (Green & Yang, 2009). Se houver modelo bifatorial, relatar o Ômega Hierárquico do fator geral ($\omega_h$), os ômegas hierárquicos das subescalas ($\omega_{hs}$) e a Variância Comum Explicada ($\text{ECV}$).
- **Modelo de Frase (AFC):**
  > "A estrutura tridimensional proposta apresentou excelente ajuste aos dados empíricos sob estimação robusta WLSMV, $\chi^2(167) = 312.45$, $p < .001$, $\text{CFI} = .978$, $\text{TLI} = .975$, $\text{RMSEA} = 0.046$, $\text{90\% CI } [0.038, 0.054]$, $p_{\text{close}} = .782$, $\text{SRMR} = 0.051$. Todas as cargas fatoriais padronizadas foram estatisticamente significativas ($p < .001$), variando de $\lambda = .58$ a $.89$. A confiabilidade composta foi elevada para todos os fatores: Fator 1 ($\omega = .88$), Fator 2 ($\omega = .85$) e Fator 3 ($\omega = .82$)."

---

### 3.8 Invariância Fatorial de Medida (Entre Grupos ou Longitudinal)
- **Sequência Hierárquica:** Configural $\to$ Métrica (cargas iguais) $\to$ Escalar (limiares/interceptos iguais) $\to$ Estrita (resíduos iguais).
- **Critérios de $\Delta$ Índices (Cheung & Rensvold, 2002; Chen, 2007):**
  - $\Delta\text{CFI} \le .010$ acompanhado de $\Delta\text{RMSEA} \le .015$ ou $\Delta\text{SRMR} \le .030$ (para métrica) e $\le .010$ (para escalar).
  - Informar que a significância de $\Delta\chi^2$ não deve ser usada como critério único de rejeição devido à sua sensibilidade extrema ao tamanho amostral.
- **Modelo de Frase:**
  > "A invariância de medida foi confirmada entre homens ($n = 150$) e mulheres ($n = 170$). O modelo métrico não apresentou degradação prática em relação ao configural ($\Delta\text{CFI} = .003$, $\Delta\text{RMSEA} = .001$, $\Delta\text{SRMR} = .004$). O modelo escalar subsequente sustentou a premissa de equivalência de limiares ($\Delta\text{CFI} = .006$, $\Delta\text{RMSEA} = .002$, $\Delta\text{SRMR} = .003$), viabilizando a comparação substantiva de médias latentes entre os grupos."

---

### 3.9 Redes Psicométricas Transversais (GGM / bootnet)
- **Especificação:** Indicar o estimador (EBICglasso com hiperparâmetro de ajuste $\gamma = 0.50$ para favorecer especificidade), o tipo de matriz de entrada (correlações policóricas para itens ordinais) e o layout visual (Fruchterman-Reingold).
- **Estabilidade e Acurácia (Isvoranu et al., 2022; Epskamp et al., 2018):**
  - Relatar a estabilidade dos índices de centralidade via bootstrap não paramétrico e o coeficiente de estabilidade de correlação ($\text{CS-coefficient}$).
  - Indicar se $\text{CS}(\text{cor} = 0.7) \ge .25$ (mínimo aceitável) ou $\ge .50$ (recomendado).
- **Comparação entre Grupos:** Relatar os resultados do Network Comparison Test (NCT), indicando a invariância da estrutura da rede ($M$, $p$) e da força global ($S$, $p$).
- **Modelo de Frase:**
  > "A rede de sintomas estimou 24 arestas não nulas de um total de 45 possíveis (esparsidade de $46.7\%$). O sintoma 'anedonia' apresentou a maior força de nó ($r_s = 1.42$), seguido por 'tristeza persistente' ($r_s = 1.25$). A análise de estabilidade por bootstrap indicou excelente replicabilidade da ordem de centralidade ($\text{CS-coefficient} = .65$). O Network Comparison Test não indicou diferenças significativas na estrutura de conectividade entre os sexos ($M = 0.18$, $p = .340$), tampouco na força global da rede ($S = 0.42$, $p = .510$)."

---

## 4. ESTRUTURAÇÃO DE TABELAS NO PADRÃO APA 7

As tabelas no padrão APA 7 devem obedecer estritamente a estas regras de design:
1. **Ausência de Linhas Verticais:** É terminantemente proibido o uso de grades ou linhas verticais separando colunas.
2. **Apenas Três Linhas Horizontais Principais:**
   - Uma linha no topo da tabela.
   - Uma linha separando os cabeçalhos das colunas do corpo de dados.
   - Uma linha na base da tabela, logo acima das notas de rodapé.
   - Linhas horizontais intermediárias só são admitidas para separar agrupamentos de cabeçalhos de múltiplos níveis (*spanner heads*).
3. **Numeração e Título:**
   - **Tabela X** (em negrito, alinhada à esquerda).
   - *Título da Tabela em Itálico* (em linha separada, caixa alta e baixa, sem ponto final).
4. **Notas de Rodapé:**
   - *Nota.* Texto explicativo geral, abreviações e detalhes do modelo.
   - *$p < .05$. **$p < .01$. ***$p < .001$.
