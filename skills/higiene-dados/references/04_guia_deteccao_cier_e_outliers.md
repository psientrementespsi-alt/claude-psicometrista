# Guia Prático de Detecção de C/IER e Outliers Robustos

Este guia consolida as diretrizes metodológicas, formulações matemáticas e rotinas computacionais para identificar respostas descuidadas ou insuficientes (Careless / Insufficient Effort Responding - C/IER) e observações atípicas multivariadas em dados psicométricos.

---

## 1. Triagem de Respostas Descuidadas (C/IER)

Conforme Curran (2016), Meade & Craig (2012) e Ward & Meade (2023), respostas descuidadas geram correlações espúrias, distorcem a dimensionalidade fatorial (gerando fatores artificiais de método) e atenuam a confiabilidade.

### 1.1 Métodos Diretos e Pós-hoc

1. **Longstring Index:**
   - *Definição:* O número máximo de respostas consecutivas com o mesmo valor de escala dado por um respondente em um bloco de itens.
   - *Cálculo:* Avalia-se o comprimento da maior sequência ininterrupta.
   - *Ponto de corte:* Valores que excedem a extensão máxima esperada sob resposta atenta (frequentemente > 6 a 10 itens idênticos consecutivos, dependendo da proporção de itens invertidos). A presença de itens invertidos é essencial para que o longstring capture desatenção e não construto real.

2. **Variância Intraindividual de Resposta (IRV):**
   - *Fórmula:*
     $$IRV_i = \sqrt{\frac{1}{K-1} \sum_{j=1}^K (X_{ij} - \bar{X}_i)^2}$$
     onde $K$ é o número de itens, $X_{ij}$ é a resposta do indivíduo $i$ ao item $j$, e $\bar{X}_i$ é a média das respostas do indivíduo $i$.
   - *Interpretação:* $IRV = 0$ indica *straightlining* perfeito (resposta invariante em toda a escala). Valores extremamente baixos ou anômalos no percentil inferior indicam negligência.

3. **Pares de Sinônimos e Antônimos Psicométricos:**
   - *Identificação dos pares:* Selecionam-se itens na amostra com correlações empíricas elevadas ($|r| > .60$).
   - *Cálculo:* Calcula-se a correlação intraindividual entre as respostas do participante nos itens do par.
   - *Interpretação:* Respondentes atentos devem apresentar correlação positiva alta nos sinônimos e negativa nos antônimos. Correlações próximas de zero ou invertidas sinalizam leitura sem apreensão de sentido.

4. **Piso Cognitivo de Tempo de Resposta (Page / Item Time):**
   - *Regra empírica:* 2 segundos por item, medido por página (Huang et al., 2012, citado por Curran, 2016, como corte deliberadamente conservador; Ward & Meade, 2023, citando Bowling et al., 2016).
   - *Cálculo:* Tempo total do questionário dividido pelo número de itens, ou tempos por página registrados pela plataforma de pesquisa.
   - *Ponto de corte:* tempo inferior a 2 s por item forma a família tempo (`corte_seg_por_item = 2` em `triagem_cier_completa()`).
5. **Ajuste de pessoa (Niessen, Meijer & Tendeiro, 2016):**
   - *G^p:* passos de item ordenados pela popularidade condicional; um erro é cada passo menos popular dado quando um mais popular não foi dado. Soma nas subescalas (G^p_m).
   - *lz^p:* log-verossimilhança padronizada do padrão sob o GRM; valores baixos indicam padrão improvável. Soma nas subescalas (lz^p_m).
   - *Ponto de corte:* empírico (quantil da amostra), por isso a família só exclui junto com outra; relatar resultados com e sem os excluídos, como pede o artigo.

### 1.2 Regra de Decisão Multimétodo
A literatura contemporânea desaconselha o descarte com base em um único índice isolado. Recomenda-se a regra de convergência: o participante só é descartado ou marcado para exclusão se violar pelo menos **2 critérios independentes** de C/IER (ex.: violar tempo mínimo e apresentar longstring anômalo).

---

## 2. Detecção Robusta de Outliers (Leys et al., 2013, 2019)

### 2.1 Rejeição da Regra Média $\pm$ 2 ou 3 SD
A média e o desvio-padrão são estimadores com ponto de ruptura de $0\%$: um único valor extremo arbitrário pode deslocar a média para o infinito e inflar o desvio-padrão. Consequentemente:
- O desvio-padrão inflado aumenta o intervalo de aceitação, deixando de detectar atípicos verdadeiros (*masking effect*).
- Não atípicos podem ser erroneamente classificados como atípicos (*swamping effect*).

### 2.2 Desvio Absoluto Mediano (MAD)
O MAD possui ponto de ruptura de $50\%$.
- *Fórmula:*
  $$MAD = b \cdot \text{mediana}(|X_i - \text{mediana}(X)|)$$
  Para distribuição normal, a constante de consistência é $b = 1.4826$.
- *Pontos de corte recomendados (Leys et al., 2013, p. 765):*
  - Moderado: $|X_i - \text{mediana}(X)| / MAD > 2.5$
  - Conservador: $|X_i - \text{mediana}(X)| / MAD > 3.0$

### 2.3 Detecção Multivariada Robusta: FastMCD
A distância de Mahalanobis clássica $D^2(X) = (X - \bar{X})' S^{-1} (X - \bar{X})$ utiliza o vetor de médias $\bar{X}$ e a matriz de covariância amostral $S$, que são vulneráveis a atípicos.
- *Solução:* Estimar o centroide robusto $\mu_{\text{MCD}}$ e a matriz de dispersão $\Sigma_{\text{MCD}}$ através do Minimum Covariance Determinant (FastMCD; Rousseeuw & Driessen, 1999).
- *Distância Robusta:*
  $$RD_i^2 = (X_i - \hat{\mu}_{\text{MCD}})' \hat{\Sigma}_{\text{MCD}}^{-1} (X_i - \hat{\mu}_{\text{MCD}})$$
- *Ponto de corte:* Quantil $\chi^2$ com $p$ graus de liberdade no nível $\alpha = .001$ ou corte ajustado com correção de consistência (Hubert et al., 2018).

---

## 3. Análise de Sensibilidade e Relato (Aguinis et al., 2013)

Conforme a taxonomia de Aguinis et al. (2013):
1. Distinguir se o atípico é decorrente de erro técnico de digitação/registro (deve ser corrigido ou convertido em NA), erro de amostragem/população inadequada (deve ser excluído) ou variação legítima do fenômeno psicológico.
2. Não deletar sumariamente os casos legítimos extremos.
3. Criar uma variável booleana `flag_outlier` no banco.
4. Conduzir a modelagem psicométrica completa na amostra total ($N_{\text{total}}$) e na amostra retida ($N_{\text{retido}}$).
5. Se as conclusões substantivas (estrutura fatorial, cargas, índices de ajuste) não se alterarem, reportar a estabilidade dos resultados. Se se alterarem, discutir a influência dos casos influentes de forma transparente.
