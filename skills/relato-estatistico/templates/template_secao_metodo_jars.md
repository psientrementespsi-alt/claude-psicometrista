# Template JARS-Quant: Seção de Método

Este template fornece a estrutura completa exigida pelas diretrizes JARS-Quant (Appelbaum et al., 2018) para a seção de Método em estudos quantitativos e psicométricos.

---

# Método

## Participantes e Delineamento Amostral
Descrever a população-alvo, os critérios de inclusão e exclusão determinados a priori, o procedimento de recrutamento (conveniência, amostragem probabilística, painel online) e o poder estatístico a priori calculado (ex.: tamanho de efeito mínimo detectável para $1 - \beta = .80$ e $\alpha = .05$). Relatar dados sociodemográficos agregados ($N$ total, percentuais de gênero/etnia, idade com média e desvio-padrão).

> "A amostra inicial foi composta por $N = [INSERIR]$ participantes recrutados via [PROCEDIMENTO]. O cálculo a priori de poder estatístico indicou que uma amostra mínima de $N = [INSERIR]$ participantes seria necessária para detectar um tamanho de efeito de magnitude [MAGNITUDE: ex., d = 0.35 ou r = .20] com poder de 80% e $\alpha = .05$. Os participantes apresentaram média de idade de $M = [INSERIR]$ anos ($SD = [INSERIR]$, amplitude de [MIN] a [MAX] anos), sendo [PERCENTUAL]% identificados com o gênero feminino e [PERCENTUAL]% com o gênero masculino."

## Instrumentos e Medidas
Para cada escala ou tarefa comportamental utilizada:
1. Nome da escala, autores e referência original de validação.
2. Descrição conceitual do construto mensurado.
3. Quantidade de itens e formato de resposta (ex.: escala Likert de 5 pontos variando de 1 = 'Discordo totalmente' a 5 = 'Concordo totalmente').
4. Sentido da pontuação e identificação de itens com pontuação invertida.
5. Evidências prévias de validade e coeficientes de confiabilidade observados no presente estudo ($\alpha$ de Cronbach e $\omega$ de McDonald com IC 95%).

> "A **Escala de [NOME DA MEDIDA]** ([CITAÇÃO AUTORES, ANO]), adaptada para o contexto brasileiro por [CITAÇÃO, ANO], avalia o construto de [CONSTRUTO]. O instrumento é composto por [K] itens organizados em [M] fatores: [NOME DO FATOR 1] ([K1] itens, ex.: '[EXEMPLO DE ITEM]') e [NOME DO FATOR 2] ([K2] itens). As respostas são registradas em escala Likert de [ESCALA: ex., 5 pontos de 1 ('Discordo totalmente') a 5 ('Concordo totalmente')]. No presente estudo, a consistência interna da escala foi de $\alpha = [VALOR]$ e $\omega = [VALOR]$ para o Fator 1 e $\alpha = [VALOR]$ e $\omega = [VALOR]$ para o Fator 2."

## Procedimentos de Higienização de Dados e Preparação Pré-Analítica
Relatar as etapas de higienização conforme o protocolo da skill `higiene-dados`:
1. Remoção de duplicatas e dados estruturalmente inválidos.
2. Triagem de respostas descuidadas e insuficientes (C/IER; Curran, 2016; Ward & Meade, 2023), especificando a regra de convergência (tempo de resposta inferior a 1,5s por item concomitante a longstring extremo ou IRV nulo).
3. Avaliação de padrões de perda de dados e teste de Little para MCAR.
4. Identificação e tratamento de atípicos multivariados via Minimum Covariance Determinant (FastMCD; Leys et al., 2018) com previsão de análise de sensibilidade.

## Plano de Análise Estatística
Descrever a sequência planejada de análises inferenciais e psicométricas:
1. Softwares e pacotes utilizados com versões e citações formais (ex.: R versão 4.x, pacote `lavaan` para AFC, `bootnet` para redes).
2. Estimadores empregados (ex.: WLSMV sobre matriz de correlações policóricas para itens ordinais; FIML com MLR para dados contínuos sob MAR).
3. Critérios de decisão e pontos de corte para índices globais de ajuste ($\text{CFI} \ge .95$, $\text{TLI} \ge .95$, $\text{RMSEA} \le .06$ com IC 90%, $\text{SRMR} \le .08$).
4. Critérios para avaliação de invariância de medida ($\Delta\text{CFI} \le .010$ e $\Delta\text{RMSEA} \le .015$).
5. Procedimento de bootstrap (5.000 reamostragens) para intervalos de confiança de efeitos indiretos em mediação.
