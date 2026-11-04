# Estudos epidemiológicos: razão de prevalência, contagens e sobrevida

Funções em `scripts/est_epidemio.R` (e `modelo_glm()` de `est_modelos.R` para
contagens); testes em `scripts/evals_latentes_epidemio.R` (`tendencia_joinpoint()` e
`taxa_padronizada()`, v1.5, em `scripts/evals_complementos.R`). Base: o módulo de
análises para estudos epidemiológicos do curso (regressão de Poisson, binomial
negativa, Poisson robusta, tábua de sobrevida, log-rank e Cox).

## Qual função

| Pergunta | Função |
|---|---|
| Desfecho binário comum em estudo transversal ou coorte: razão de prevalência (RP) | `razao_prevalencia(dados, desfecho, preditores)` |
| Contagem (número de eventos): taxa, razão de taxas (IRR) | `modelo_glm(..., familia = "poisson")`; com superdispersão, `familia = "binomial_negativa"` (`references/modelos.md`) |
| Tempo até o evento com censura | `sobrevida(dados, tempo, evento, grupo, covariaveis)` |

## Razão de prevalência (Poisson robusta)

Com desfecho comum, a razão de chances da logística se afasta da RP e
superestima a associação quando lida como risco (Barros & Hirakata, 2003). A
Poisson com ligação log e erro padrão sanduíche estima a RP direto (Zou, 2004).
O EP robusto é o HC0, o mesmo da GEE independente com um caso por cluster; os
evals conferem com o geepack e com a conta à mão. Prevalência prevista acima de
1 para alguém vira aviso (limitação conhecida da Poisson para desfecho
binário; a RP segue válida). Evento: a mesma regra da `regressao_logistica()`
(0/1 = 1, fator = segundo nível, texto exige `evento`).

## Sobrevida

`sobrevida()` devolve, numa chamada:

- Kaplan-Meier: n, eventos, mediana com IC (transformação log, padrão do
  survival) por grupo; mediana NA quando a curva não cai abaixo de 50%.
- Tábua de sobrevida nos `tempos` pedidos (em risco, eventos, sobrevida e IC).
- Log-rank (Mantel-Cox) entre os grupos.
- Cox com `grupo` e `covariaveis`: HR com IC de Wald, empates por Efron, e
  riscos proporcionais pelos resíduos de Schoenfeld (`cox.zph`); violação vira
  aviso (estratificar ou interação com o tempo). Menos de 10 eventos por
  parâmetro vira aviso.

`evento` 0/1 ou lógico; com outro código (o `lung` do survival usa 2 = óbito),
`evento_codigo = 2`. Os evals reproduzem o log-rank χ²(1) = 10,33 e o HR de 0,588
[0,424; 0,816] do banco `lung`.

## Tendência temporal por pontos de inflexão (v1.5)

`tendencia_joinpoint(dados, tempo, casos, populacao)` ou `taxa =`: o modelo do programa
Joinpoint (Kim et al., 2000) que o curso de séries temporais usa, em R. Log-linear
contínuo, pontos procurados nos anos observados, pelo menos 2 observações da ponta e 2
entre pontos, máximo de pontos pela regra do Joinpoint (até 11 anos, 1; 12 a 16, 2...),
número de pontos pelo BIC do Joinpoint. APC por segmento, AAPC do período (média das
inclinações ponderada pelo comprimento dos segmentos), IC paramétrico (t com n - 2k - 2
gl) ou por quantis empíricos (`ic = "empirico"`, reamostra os resíduos e reajusta a
localização dos pontos). Tendência crescente ou decrescente quando o IC não inclui zero,
estacionária nos demais casos, como no relato das aulas.

Âncora: com a incidência de tuberculose 2010-2019 da aula (masculino, feminino e
ambos), os pontos (2015, 2016, 2015), as APC e as AAPC batem com a saída do Joinpoint
4.9.1 do relato da aula na casa decimal relatada (masculino: -1,0 e 2,7, AAPC 0,7). Os IC
do relato da aula são mais largos que os dois métodos da skill (masculino 2010-2015:
[-2,3; 0,3] no Joinpoint; [-1,7; -0,2] paramétrico; [-2,0; -0,3] empírico, 4.999 reamostras), e a sessão
salva (.jps/.jpo) não deixa ler o método de IC nem a opção de variância usados. Para IC
idêntico ao do programa, rodar o Joinpoint; com a skill, declarar o método de IC no relato.
Permutação (Kim et al., 2000) e BIC ponderado do Joinpoint ficam fora; autocorrelação
dos erros também.

## Taxa padronizada (v1.5)

`taxa_padronizada(dados, casos, populacao, padrao, grupo)`: padronização direta (média
das taxas por estrato de idade ponderada pela população padrão), com IC gama de Fay e
Feuer (1997); igual ao `epitools::ageadjust.direct` nos evals. É a "análise ajustada por
idade" das aulas de séries temporais. Declarar a população padrão (Censo, OMS).

## Fora das funções

Modelos que a aula só cita, sem código: regressão zero-inflada e zero-truncada (contagem
com excesso de zeros que tem significado próprio; `pscl::zeroinfl` e `pscl::hurdle`,
comparados com a binomial negativa do `modelo_glm()` por AIC e pelo teste de Vuong, que é
controverso para modelos zero-inflados) e riscos competitivos, quando a censura é
informativa ou há mais de um desfecho que impede o outro (`cmprsk`/`tidycmprsk`, incidência
cumulativa de Fine e Gray). Rodar à mão, no computador do usuário, e conferir com o
material da disciplina antes de relatar.

Acesso a dados abertos (DATASUS, Tabwin, Tabnet, IBGE, INEP, IPEA, INMET), que o curso
também ensina, não é análise estatística da skill. Para ler arquivos .dbc do DATASUS no R
há o pacote `read.dbc` (o `extrair_dbc.R` das aulas); a limpeza vira bloco registrado do
script (`references/script-reproduzivel.md`). Microdados individuados do DATASUS são dado
de pessoa: valem as regras de privacidade da skill (rodar no computador, só agregados na
conversa).
