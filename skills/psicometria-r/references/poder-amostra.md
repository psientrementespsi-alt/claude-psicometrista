# Poder e tamanho amostral (reproduzindo o G*Power 3.1)

## Sumário

- 11 (linha 24)
- 52 (linha 65)
- 94 (linha 107)
- 119 (linha 132)
- 153 (linha 166)
- 164 (linha 177)
- 174 (linha 187)
- 205 (linha 218)
- 249 (linha 262)
- 294 (linha 307)

Funções em `scripts/est_poder.R`; testes em `scripts/evals_poder.R` (108 testes na
v1.4, cerca de 13 s). Todas devolvem `pr_resultado` (subclasse `pr_poder`, com impressão em lista).

Índice: 1. Qual função usar · 2. O que toda função faz · 3. Tamanho de efeito e
conversões · 4. Medidas repetidas e a opção "as in Cohen (1988)" · 5. Perda prevista ·
6. Post hoc e sensibilidade · 7. Como relatar · 8. Armadilhas · 9. Fontes e divergências ·
10. AFC e MEE (v1.4)

## 1. Qual função usar

O cálculo é sempre para um teste e um efeito específicos (livro): primeiro
fixar a análise final, depois o efeito que ela vai testar.

| Análise planejada | Função | G*Power (Test family / Statistical test) |
|---|---|---|
| Correlação de Pearson | `poder_correlacao(r)` | Exact / Correlation: Bivariate normal model |
| t independente | `poder_teste_t(d)` | t tests / Means: Difference between two independent means (two groups) |
| t pareado (pré e pós) | `poder_teste_t(dz, "pareado")` | t tests / Means: Difference between two dependent means (matched pairs) |
| t de uma amostra | `poder_teste_t(d, "uma_amostra")` | t tests / Means: Difference from constant (one sample case) |
| ANOVA de um fator | `poder_anova(f, niveis = k)` | F tests / ANOVA: Fixed effects, omnibus, one-way |
| ANOVA fatorial (um termo) | `poder_anova(f, niveis = c(a, b), efeito = "A:B")` | F tests / ANOVA: Fixed effects, special, main effects and interactions |
| ANCOVA | `poder_anova(..., covariaveis = c)` | F tests / ANCOVA: Fixed effects, main effects and interactions |
| Medidas repetidas: intra, entre, interação | `poder_anova_mr(f, medidas, grupos, efeito)` | F tests / ANOVA: Repeated measures, within factors / between factors / within-between interaction |
| Regressão: R² do modelo | `poder_regressao("global", p, r2 = )` | F tests / Linear multiple regression: Fixed model, R² deviation from zero |
| Regressão hierárquica: ΔR² do bloco | `poder_regressao("incremento", p, r2_reduzido = , delta_r2 = , testados = q)` | F tests / Linear multiple regression: Fixed model, R² increase |
| Regressão: um coeficiente | `poder_regressao("coeficiente", p, r2_parcial = )` | t tests / Linear multiple regression: Fixed model, single regression coefficient |
| Regressão logística (um preditor) | `poder_logistica(or, p0)` | z tests / Logistic regression |
| Qui-quadrado | `poder_quiquadrado(w, gl)` | χ² tests / Goodness-of-fit tests: Contingency tables |
| AFC: ajuste do modelo pelo RMSEA (v1.4) | `poder_afc_rmsea(rmsea0, rmsea1, gl)` | fora do G*Power (seção 10) |
| AFC: largura do IC do RMSEA (v1.4) | `n_aipe_rmsea(rmsea, gl, largura)` | fora do G*Power (seção 10) |
| AFC/MEE: um parâmetro por modelos aninhados (v1.4) | `poder_afc_aninhados(...)` | fora do G*Power (seção 10) |
| AFC/MEE: largura do IC de um parâmetro (v1.4) | `n_aipe_parametro(modelo, Sigma, parametro, largura)` | fora do G*Power (seção 10) |

Fora do escopo: mediação, moderação, modelos mistos e redes, que pedem simulação
(`simr::powerSim` para modelos mistos; para rede, o livro indica simular cenários com o
netSimulator do bootnet). A AFC entrou na v1.4 pelas quatro abordagens do
capítulo 5 do Wiley Handbook (seção 10); simulação de Monte Carlo para MEE complexa
continua fora. GLM e GEE: as aulas "Cálculo amostral - GLM/GEE" usam um programa web
(o G*Power não cobre) que pede as médias por condição, o DP, a correlação entre medidas
repetidas e a estatística de teste multivariada, que é o fluxo do GLIMMPSE (Kreidler et
al., 2013); a skill não tem função para isso. Alternativa no R: simular o desenho com os
valores esperados e ajustar `modelo_gee()` a cada réplica, contando as rejeições, com a
ressalva de rota fora do escopo testado.

Amostragem (aulas "Métodos de amostragem", "Amostragem probabilística" e "não
probabilística"): o cálculo de N supõe amostra aleatória da população; com amostra por
conveniência ou bola de neve o N garante poder para o efeito na amostra, não
representatividade. Dizer o método de amostragem no Método e não generalizar além dele.

## 2. O que toda função faz

- Sem `n`: resolve o N mínimo para o `poder` pedido (padrão .80, α = .05 bicaudal).
  Com `n`: devolve o poder alcançado (post hoc, com aviso; seção 6).
- Imprime, em três blocos: "G*Power, entrada" (os campos na ordem da tela, para refazer
  a conta), "G*Power, saída" (λ ou δ, valor crítico, gl, N, poder real) e "Skill"
  (N mínimo inteiro, N com grupos ou células iguais, N com perda). A saída usa as casas
  do G*Power: λ, δ, valores críticos e poder com 6 casas (seção 11.3 do manual: poder
  0.99999, que com 4 casas virava "1"); gl e N inteiros. No t de uma amostra o campo é
  "Effect size d", como no G*Power (manual, seção 20); dz é só do pareado.
- Motor: distribuições não centrais do R com os parâmetros do manual do G*Power
  (λ = f²·N para ANOVA, ANCOVA e regressão; δ = d·√(n1n2/(n1+n2)) ou dz·√N no t;
  λ = w²·N no χ²). Correlação exata e logística usam o pwrss 1.3.1
  (`power.exact.onecor`, `power.z.logistic`). Os evals conferem a equivalência com
  `pwr::pwr.t.test`, `pwr.t2n.test`, `pwr.anova.test`, `pwr.f2.test`,
  `pwr.chisq.test` e com `pwrss::power.f.ancova`, `power.f.mixed.anova` e
  `power.f.regression`.
- Arredondamento, como o G*Power: ANOVA de um fator e medidas repetidas em múltiplo do
  número de grupos (f = .25: 159, 180 e 200 para 3, 4 e 5 grupos; mínimos 158, 179 e
  196); fatorial, ANCOVA, regressão, correlação e logística no N mínimo (244 na
  interação do 3 × 5 do curso com ηp² = .06, que não é múltiplo de 15; a skill mostra ao
  lado o N com células iguais, 255); t independente com n por grupo arredondado para
  cima. Com `razao` (N2/N1) diferente de 1, a regra do G*Power: resolve n1* contínuo com
  n2 = razão × n1 e devolve teto(n1*) e teto(razão × n1*), cada grupo por conta própria
  (`poder_teste_t(.5, razao = 1.571428)` dá 53 + 82 = 135; buscar o menor n1 inteiro
  daria 53 + 84). A regra está conferida no manual em dois exemplos que não são do t de
  duas médias (seção 9). Com efeito muito grande o menor desenho (pelo menos 2 casos em
  cada grupo) já passa do poder pedido: a função devolve esse desenho, com aviso
  (`poder_teste_t(3, razao = 3)`: 2 + 6 = 8).
- Casas do efeito na frase: o menor número de casas (2 a 4) com que o efeito
  arredondado devolve o mesmo N. Em `poder_anova_mr(eta2p = .10, medidas = 4, rho = NA,
  poder = .90)`, f = .3333 sai como "0.3333" e o N é 129, porque 0.333 daria 130 e 0.33
  daria 132 (a armadilha da aula de medidas repetidas). Na regressão o f² começa em 3
  casas. A conferência usa o mesmo arredondamento da frase (meio para cima):
  `poder_teste_t(.625)` escreve "d = 0.625" e 84, porque 0.63 daria 82.
- α e poder na frase sem zero à esquerda e com as casas do valor usado na conta: α .05,
  .01, .001, .005; poder .80, .90, .975, .999; um valor sem forma curta, como o
  Bonferroni .05/3, sai com três algarismos significativos (.0167).
- `justificativa`: a origem do efeito entra na frase entre parênteses; sem ela, a frase
  sai com "[origem do tamanho de efeito: preencher]". O livro exige o efeito e a origem
  antes do N; o curso pede a justificativa antes do valor.

## 3. Tamanho de efeito e conversões

A peça difícil é o efeito, não a conta. Ordem do livro: literatura da
área, meta-análise, estudos comparáveis, piloto; convenção de magnitude só como ponto de
partida declarado; superestimar o efeito é o erro mais caro. O curso acrescenta o menor
efeito de interesse (SESOI), em unidades brutas convertidas pelo DP esperado.

| Conversão | Função | Fórmula |
|---|---|---|
| η² ou ηp² para f | `efeito_f_de_eta2()` | f = √(η² / (1 − η²)) |
| f para η² | `efeito_eta2_de_f()` | η² = f² / (1 + f²) |
| d de médias e DP | `efeito_d_de_medias(m1, m2, dp1, dp2, n1, n2)` | DP = √((DP1² + DP2²)/2), o "Determine" do G*Power; com n, DP combinado por gl |
| f de médias (um fator) | `efeito_f_de_medias(medias, dp, n)` | f = DP das médias (ponderado por n) / DP dentro |
| r para d e d para r | `efeito_d_de_r()`, `efeito_r_de_d(d, n1, n2)` | d = 2r/√(1 − r²); r = d/√(d² + a), a = (n1+n2)²/(n1n2) |
| d de t | `efeito_d_de_t(t, n1, n2)` | d = t·√(1/n1 + 1/n2) |
| f² de R² | `efeito_f2_de_r2()` | f² = R² / (1 − R²); vale para R² parcial |
| f² de ΔR² | `efeito_f2_de_delta_r2(r2_reduzido, delta_r2)` | f² = ΔR² / (1 − (R²0 + ΔR²)) |
| OR de probabilidades | `efeito_or_de_probabilidades(p1, p0)` | (p1/(1 − p1)) / (p0/(1 − p0)) |
| dz de d entre momentos | `efeito_dz_de_d(d, r)` ou `r_pareado =` | dz = d / √(2(1 − r)) |

Magnitudes de Cohen (1988; 1992) como ponto de partida: d .20/.50/.80; r .10/.30/.50;
f .10/.25/.40; f² .02/.15/.35; w .10/.30/.50; η² .01/.06/.14 é a derivação de f. O
curso cita R² ≈ .13 (médio) e ≈ .26 (grande). Para OR não há corte universal: depende
da incidência (seção 8).

## 4. Medidas repetidas e a opção "as in Cohen (1988)"

`poder_anova_mr(f, medidas, grupos = 1, efeito = "intra" | "entre" | "interacao", rho = 0.5, epsilon = 1)`.

- Modo padrão (`rho` numérico) = o padrão do G*Power ("as in GPower 3.0"), descrito em
  Faul et al. (2007) e implementado no `pwrss::power.f.mixed.anova`: λ = f²·N·m·ε/(1 − ρ)
  para intra e interação; λ = f²·N·m/(1 + (m − 1)ρ) para entre; gl intra (m − 1)ε e
  (N − k)(m − 1)ε; interação (k − 1)(m − 1)ε e (N − k)(m − 1)ε; entre k − 1 e N − k.
  Reproduz o exemplo do artigo (f = .2572, N = 90, 3 grupos, 3 medidas, ρ = .3: λ =
  25.52, gl 2 e 174, F crítico 3.048, poder .997; interação .653; entre .488). Aqui o f
  é o do efeito sem a correlação; ρ entra pelo campo "Corr among rep measures" (padrão
  da tela 0,5).
- `rho = NA`: o f já incorpora a correlação (rota `rho.within = NA` do pwrss): λ =
  f²·N·ε. É a rota que reproduz os dois exemplos do curso (132 no intra com f = .33, 4
  medidas e poder .90; 309 na interação 3 × 2 com f = .204 e poder .90), calculados no
  G*Power com Options = "as in Cohen (1988)".
- **A opção "as in Cohen (1988)" não está implementada como tal.** O manual do G*Power
  (versão de 2023) não tem capítulo de medidas repetidas e Faul et al. (2007) não
  descrevem essa opção; a fórmula dela não foi encontrada em fonte primária. A
  coincidência com λ = f²·N·ε vale nos dois exemplos do curso, que é pouco para
  afirmar equivalência geral. A função não oferece argumento com esse nome, avisa
  quando `rho = NA` e pede conferir o N na tela do G*Power antes de escrever "as in
  Cohen (1988)" no método. A opção "as in SPSS" também não está implementada.
- Um ηp² publicado de ANOVA de medidas repetidas (SPSS, JASP) já embute a correlação:
  com ele, `rho = NA` é a leitura coerente; no modo padrão o f precisa ser o do efeito
  sem a correlação.
- N em múltiplos do número de grupos (o 309 do curso; mínimo inteiro 308). Esfericidade:
  com 2 medidas não se aplica (ε = 1, nada a relatar); com 3 ou mais, ε = 1 é a escolha
  liberal (intra, 4 medidas, `rho = NA`, f = .33 e poder .90: 132 com ε = 1, 162 com
  .75, 217 com .50).
- "Number of groups" é o número de grupos entre sujeitos (desenho só intra: 1). A aula
  usa 4 no exemplo intra; ali não muda o N (132 nos dois), mas pode mudar em outros
  desenhos.

## 5. Perda prevista

`perda = 0.20` acrescenta "N para coleta" = teto(N / (1 − perda)), em múltiplo do número
de grupos quando o desenho pede grupos iguais (one-way com f = .23 e 3 grupos: 186, que
com 20% vira 234). Esse
ajuste **não vem das fontes do usuário**: o curso não ensina ajuste por perda, e o livro
diz apenas que é boa prática inflar o N alvo em 10% a 20% como decisão do pesquisador
. Feito como N × (1 + perda), o inflado fica
abaixo do planejado depois da perda (62 × 1,2 = 75; 75 × 0,8 = 60). A função mostra a
conta do livro no aviso e usa a divisão. `ajustar_perda(n, perda)` faz só essa conta.

## 6. Post hoc e sensibilidade

Com `n` informado a função devolve o poder alcançado e sempre avisa: poder calculado com
o efeito observado na própria amostra é função direta do p (no limiar de 5%, perto de
50%) e não acrescenta informação (livro e; o curso, na página "Tipos de
análise de poder", diz o mesmo, embora a página "Cálculo de Tamanho Amostral" proponha o
contrário: vale a primeira). O uso legítimo é com efeito definido antes (SESOI,
meta-análise). Análise de sensibilidade (menor efeito detectável com N fixo) não tem
função própria: rodar a função com `n` para uma grade de efeitos, ou ver Lakens (2022).

## 7. Como relatar

```r
source("~/.local/share/psicometria-r/carregar.R")
poder_correlacao(.35, justificativa = "magnitude de achado anterior comparável")
poder_anova(eta2p = .06, niveis = c(3, 5), efeito = "A:B", justificativa = "...")
poder_anova_mr(eta2p = .04, medidas = 2, grupos = 3, efeito = "interacao", rho = NA, poder = .90)
poder_regressao("incremento", preditores = 5, r2_reduzido = .20, delta_r2 = .05, testados = 3)
poder_logistica(or = 1.52, p0 = .30, perda = .15)
```

Frase gerada por `poder_anova(f = .23, niveis = 3, perda = .20, justificativa = "DP 6,15
do estudo normativo e médias 12, 12 e 15 como diferença mínima de interesse")` (one-way
do curso, com perda de 20%): "O tamanho amostral foi determinado
a priori com os parâmetros do G*Power 3.1 (Faul et al., 2007) para uma ANOVA de um fator
com 3 grupos, com f = 0.23 (DP 6,15 do estudo normativo e médias 12, 12 e 15 como
diferença mínima de interesse), α = .05 e poder de .80, resultando em um mínimo de 186
participantes (62 por grupo). Prevendo perda de 20%, a meta de coleta é de 234
participantes." O formato segue a frase-modelo do livro com uma
diferença: α e poder saem sem zero à esquerda (com as casas do valor usado), pela regra
do próprio livro para quantidades que não passam de 1; a frase-modelo escreve
"poder de 0,80". A frase cita o G*Power: rodar no
G*Power as configurações impressas e conferir o N antes de usar; se o relato for do R,
citar o pacote (`citation("pwrss")`). Duas exceções em que a frase não cita o G*Power,
porque o N do resultado não é o da tela: `poder_correlacao(..., metodo = "fisher")`
("pela aproximação z de Fisher"; o bloco "G*Power, saída" continua com o N exato, 61
para ρ = .35, e o 62 de Fisher fica no bloco da skill) e a logística com R² other X > 0
e preditor binário (armadilha 11). O relato completo pede: análise e efeito-alvo (na fatorial, qual
termo), efeito com a origem, α e caudas, poder, parâmetros do desenho (grupos, medidas,
ρ, ε, preditores testados e totais, p0 e distribuição do preditor), N e meta de coleta.

## 8. Armadilhas

1. One-way do curso: a aula relata "≈ 189"; com f = .23 e 3 grupos o G*Power dá 186 (62
   por grupo), e o texto-modelo da própria aula já cita 186. Teste: 186.
2. Fatorial 3 × 5: a aula escreve f ≈ .225, mas ηp² = .06 dá f = .2526, o único valor
   compatível com o N = 244 da tela (com .225 seriam 305). Erro de transcrição.
3. Preditor único: o "R² parcial = .067" da aula é o f² (de R²p = .0625). Digitar .067
   como R² parcial dá 112, não 120. R² parcial = β² só vale com preditores não
   correlacionados.
4. t pareado: o G*Power usa dz. Um d, g ou Δ calculado com o DP dos escores só é igual a
   dz quando r = .50 entre as medidas: `poder_teste_t(.30, "pareado", r_pareado = r,
   poder = .90)` dá N = 166, 119 e 72 com r = .30, .50 e .70 (com o poder padrão .80,
   125, 90 e 55). Informar `r_pareado` ou o tipo de d.
5. Correlação: `pwr::pwr.r.test` dá 85 onde o G*Power exato dá 84 (ρ = .30). O livro faz
   a conta de Fisher no texto (62) e mostra a tela exata (61); a função usa o exato e
   mostra Fisher ao lado (seção 9).
6. ΔR² e coeficiente: `pwr::pwr.f2.test` usa λ = f²(u + v + 1), que só é f²·N no modelo
   global; no incremento diverge (f² = .15, 2 de 5 preditores: G*Power 68, pwr 71).
7. Logística: o N depende muito de p0 (`poder_logistica(1.52, p0)`, preditor normal,
   Demidenko com correção, poder .80: 898 com p0 = .05; 227 com .30; 196 com .50). p0 deve refletir a prevalência esperada do desfecho. Os cortes de OR de Chen et
   al. (2010) são para exposição binária; com preditor contínuo o OR é por 1 DP.
8. ANCOVA: o N não muda com as covariáveis só se o ηp² assumido já for o ajustado; se
   veio de estudo sem a covariável, f ajustado = f / √(1 − r²).
9. ηp² de outro desenho (outros fatores, covariáveis, medidas repetidas) não é
   intercambiável com o do desenho planejado.
10. Regra 50 + 8k (citada pelo curso) é só checagem: com f² = .15 e poder .80,
    `poder_regressao("global", k, f2 = .15)` dá 85 para k = 4 (regra: 82) e 118 para
    k = 10 (regra: 130).
11. Logística com R² other X > 0 e a opção padrão (Demidenko com correção de variância):
    o manual do G*Power documenta N/(1 − R²) com a correção mantida (eq. 14-15), mas o
    único exemplo com R² > 0 (seção 29.4: OR = 1,5, p0 = .5, poder .95, R² = .1, preditor
    normal) imprime 395, que é o N sem a correção dividido por 1 − R² (354,74 / 0,9); a
    fórmula documentada daria 375. A skill segue o exemplo, mostra o valor da fórmula
    documentada no bloco da skill e avisa que a regra vem de um só exemplo. O caminho
    sem correção com R² other X reproduz os dois exemplos de análise de compromisso de
    Faul et al. (2009), com R² = .40 (z crítico 3,454879 e 2,146971, reproduzidos nos
    evals pelo poder post hoc da função com 1 − poder = α); o que fica apoiado numa só
    tela é a combinação "com correção" mais R² > 0 com preditor normal (395); com
    correção, preditor binário e R² > 0 não há tela. Com o curso
    (OR = 1,52, p0 = .30, R² = .20) isso dá 293, e não 284. Nessa última combinação
    (com correção, binário, R² > 0) o N sai no bloco da skill, com aviso de que não foi conferido
    contra o G*Power, e a frase descreve o procedimento (Demidenko, 2007, sem correção;
    fator de Hsieh et al., 1998) sem citar o programa. Conferir na tela antes de relatar.

## 9. Fontes e divergências

Livro (Biblioteca Análise de Dados, Vol. 1, cap. 7): sistema de quatro peças e origem do
efeito; poder relativo a teste e efeito; G*Power para correlação
; números-âncora r = .35: Fisher 62 e 82, exato 61 e 81; r = .10: Fisher 783,
exato 782, poder .61 com N = 500; conta de Fisher;
perda de 10% a 20%; erros; frase-modelo;
interação como efeito mais exigente. Todos reproduzidos nos evals.

Divergência livro × G*Power na correlação: o texto chama Fisher de padrão das
ferramentas, mas a tela do G*Power usa a distribuição exata (61, não 62). A
skill usa o exato por padrão (é o que o G*Power devolve e o que o curso usa) e imprime o
N de Fisher em `$decisoes`; `metodo = "fisher"` reproduz a conta do texto.

Curso (Psicometria Online Academy, páginas "Cálculo de Tamanho Amostral", "Tipos de
análise de poder", "Tamanho de Efeito para ANOVAs no G*Power" e as páginas de cada
procedimento no G*Power, de "Correlação: Cálculo Amostral" a "Regressão Logística -
Cálculo Amostral"): todos os exemplos foram recalculados; os reproduzidos estão nos
evals, com as exceções das armadilhas 1 a 3. O ΔR² da aula não tem N impresso; 168 é
cálculo da skill, conferido contra o pwrss. O 284 da logística com R² other X = .20
(sensibilidade do material) era conta do pwrss com a correção mantida, não tela do
G*Power (armadilha 11).

Manual do G*Power 3.1 (versão de 1/6/2023), conferido nos evals: logística, seção 29.4
(337, 355 e 317 com preditor normal; 1437, 2158, 2368 e 1498 com preditor binário; 1437
também pelo Hsieh com preditor binário; poder alcançado do Hsieh 0.950486 com 317 e
0.950068 com 1437, pela fórmula de Hsieh et al. resolvida para o poder, porque o pwrss
devolve o poder pedido; 395 com R² other X = .1, armadilha 11); ANOVA fatorial, seção 11.3
(λ 53.935690, F crítico 3.123907, poder 0.99999); t de uma amostra, seção 20.3 (30 e
1492); regra de alocação com razão diferente de 1, seção 15
(duas variâncias, razão 2: 148 + 295) e seção 17 (diferença entre inclinações, teste t,
razão 1,571428: 163 + 256, δ = 2.811598). Nenhum exemplo do manual é do t de duas médias
com razão diferente de 1; ali a regra é a mesma, mas sem tela conferida.

Literatura: Faul, Erdfelder, Lang e Buchner (2007), DOI 10.3758/BF03193146 (G*Power 3 e
o exemplo de medidas repetidas); Faul, Erdfelder, Buchner e Lang (2009), DOI
10.3758/BRM.41.4.1149 (logística com R² other X pelo caminho de Demidenko sem correção:
os dois exemplos de análise de compromisso, reproduzidos nos evals); Cohen (1988), DOI da reedição 10.4324/9780203771587, e
Cohen (1992), DOI 10.1037/0033-2909.112.1.155 (convenções de magnitude); Lakens (2022),
DOI 10.1525/collabra.33267 (justificativa de N e sensibilidade). Citados pelo curso e
usados nas opções, com DOI conferido em `referencias-verificadas.md`: Demidenko (2007) e
Hsieh et al. (1998) (métodos da logística no G*Power e no pwrss) e Chen, Cohen e Chen
(2010) (OR).
Pacotes: pwrss 1.3.1 e pwr 1.3-0 (`citation()` no R da sessão).

## 10. AFC e MEE (v1.4)

Kelley e Lai (2018, cap. 5 do Wiley Handbook, na pasta Psicometria-R) separam quatro
perguntas; cada uma tem uma função, e os exemplos do capítulo são reproduzidos nos
evals. O tamanho amostral só precisa responder à pergunta que o estudo faz.

1. **Poder para o RMSEA do modelo** (MacCallum, Browne & Sugawara, 1996):
   `poder_afc_rmsea(rmsea0, rmsea1, gl)`. λ = (N − 1)·gl·ε²; com ε0 < ε1 é o teste de
   ajuste aproximado (H0 ε ≤ ε0; o clássico é ,05 contra ,08), com ε0 > ε1 o de ajuste
   não aproximado (H0 ε ≥ ε0, rejeitar mostra que o ajuste é melhor que o corte), e
   ε0 = 0 o de ajuste exato (pouco informativo: o ajuste perfeito já é implausível).
   `grupos` divide λ pelo número de grupos, como o semTools.
2. **AIPE do RMSEA** (Kelley & Lai, 2011): `n_aipe_rmsea(rmsea, gl, largura, nivel =
   .90)`, o N para a largura esperada do IC não passar da pedida. Largura esperada
   quer dizer que metade dos estudos com esse N terá IC mais largo.
3. **Poder para um parâmetro por modelos aninhados**: `poder_afc_aninhados()` com os
   RMSEA dos dois modelos (MacCallum, Browne & Cai, 2006), com os F de máxima
   verossimilhança, ou com Σ populacional e as duas sintaxes lavaan (Satorra & Saris,
   1985); λ = (N − 1)·ΔF com Δgl.
4. **AIPE de um parâmetro**: `n_aipe_parametro(modelo, Sigma, parametro, largura)`,
   pela escala do erro padrão assintótico da máxima verossimilhança (a garantia,
   *assurance*, do MBESS não está implementada).

O gl do item 1 é o do modelo que vai ser testado (p(p + 1)/2 menos os parâmetros
livres). Σ populacional sai de correlações e desvios-padrão plausíveis ou da matriz
implicada por cargas propostas; o capítulo mostra como montar.

Números do capítulo e da skill: 1.268 e 1.116 nos dois testes do RMSEA (o capítulo dá
1.269 e 1.113 pelo MBESS, que para a busca quando o poder fica a ,001 do alvo; o
semTools dá os exatos, iguais aos da skill); 802, 1.059 e 1.544 no AIPE do RMSEA; 336,
268 e 365 nos aninhados por RMSEA (o capítulo dá 364 no último, pelo mesmo motivo);
1.007 no AIPE de λ1. As funções devolvem `pr_resultado` com uma tabela de campos (não a
tela do G*Power, que não tem esses testes) e frase APA.
