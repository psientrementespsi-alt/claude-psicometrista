# Template APA 7: Seção de Resultados

Este template fornece o esqueleto canônico para redação da seção de Resultados em psicologia e ciências do comportamento, cobrindo análises preliminares, testes de hipóteses principais, tamanhos de efeito e análise de sensibilidade.

---

# Resultados

## Análises Preliminares e Triagem Amostral
Apresentar os dados de fluxo amostral decorrentes da higienização, estatísticas descritivas básicas ($M, SD$, assimetria e curtose) e matriz de correlação preliminar bivariada.

> "A triagem pré-analítica de dados resultou na exclusão de $n = [EXCLUSOES]$ participantes por violações de consistência e respostas descuidadas, consolidando uma amostra analítica final de $N = [N\_FINAL]$ participantes (conforme detalhado na Tabela 1). A inspeção dos valores de assimetria ($|Sk| < 1.5$) e curtose ($|Ku| < 2.0$) indicou adequação à distribuição unimodal. A Tabela 2 apresenta as médias, desvios-padrão e correlações de Pearson entre todas as variáveis do estudo. Conforme hipotetizado, a variável [VAR\_X] associou-se positivamente a [VAR\_Y], $r([DF]) = .[VALOR]$, $p < .001$."

## Estrutura Fatorial e Evidências Psicométricas
Apresentar os resultados da modelagem fatorial (AFE, AFC ou Bifator), parâmetros locais e confiabilidade.

> "A Análise Fatorial Confirmatória da escala demonstrou índices excelentes de ajuste global aos dados sob estimação robusta WLSMV: $\chi^2([DF]) = [VALOR]$, $p = [P\_VAL]$, $\text{CFI} = .[VALOR]$, $\text{TLI} = .[VALOR]$, $\text{RMSEA} = 0.[VALOR]$ (IC 90% $[0.[MIN], 0.[MAX]]$) e $\text{SRMR} = 0.[VALOR]$. Todas as cargas fatoriais padronizadas foram elevadas e estatisticamente significativas ($p < .001$), variando entre $\lambda = .[MIN]$ e $\lambda = .[MAX]$ (Tabela 3). Os coeficientes de consistência interna confirmaram a precisão da medida, com $\alpha = .[VALOR]$ e $\omega = .[VALOR]$ para o Fator Geral."

## Testes de Hipóteses Principais

### Exemplo A: Comparação entre Grupos (Testes de Diferença)
> "Para testar o efeito da condição experimental, foi conduzido um teste *t* de Student para amostras independentes. Os participantes submetidos à intervenção relataram níveis significativamente superiores de bem-estar subjetivo ($M = 4.25$, $SD = 0.58$) em comparação aos participantes do grupo controle ($M = 3.60$, $SD = 0.72$), $t([DF]) = [T\_STAT]$, $p < .001$, $d = 0.99$, IC 95% $[0.65, 1.33]$, correspondendo a um tamanho de efeito de grande magnitude."

### Exemplo B: Modelo de Regressão e Processos (Mediação)
> "Para avaliar se o efeito de [PREDITOR] sobre [DESFECHO] é mediado por [MEDIADOR], foi conduzida uma análise de mediação linear com 5.000 reamostragens por bootstrap. O preditor exerceu impacto direto significativo sobre o mediador, $a = [VALOR]$, $SE = [SE]$, $t([DF]) = [T]$, $p < .001$, e o mediador associou-se significativamente ao desfeito controlando pelo preditor, $b = [VALOR]$, $SE = [SE]$, $t([DF]) = [T]$, $p < .001$. O efeito indireto do preditor sobre o desfecho via mediador foi estatisticamente significativo, $ab = [VALOR]$, $SE = [SE]$, IC 95% bootstrap $[[MIN], [MAX]]$. Como o intervalo de confiança não intercepta o valor zero, a hipótese de mediação foi corroborada."

## Análise de Sensibilidade para Casos Atípicos
Seguindo a recomendação de Aguinis et al. (2013), reportar a re-estimação do modelo após controle dos atípicos multivariados.

> "Para verificar a robustez dos achados diante de observações extremas, foi realizada uma análise de sensibilidade re-estimando o modelo final com a exclusão dos $n = [N\_OUT]$ casos identificados como atípicos multivariados pelo estimador FastMCD. Os parâmetros estimados mantiveram estabilidade substantiva: o efeito indireto manteve-se significativo ($ab = [VALOR]$, IC 95% $[[MIN], [MAX]]$) e os índices de ajuste fatorial não sofreram alterações relevantes ($\Delta\text{CFI} < .005$). Esses resultados evidenciam que os efeitos observados não foram impulsionados por observações atípicas."
