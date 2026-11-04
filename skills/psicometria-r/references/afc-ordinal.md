# AFC com itens ordinais e quase contínuos

## Sumário

- 12 (linha 23)
- 54 (linha 65)
- 95 (linha 106)
- 131 (linha 142)
- 147 (linha 158)
- 201 (linha 212)
- 286 (linha 297)
- 338 (linha 349)

Índice: Escolha do estimador · Código · Índices de ajuste · Cargas e solução · Confiabilidade · Escore somado · Modelo bifator · Validade convergente e discriminante

Funções: `auditar_itens()`, `sugerir_estimador()`, `ajustar_afc()`,
`tabela_ajuste()`, `confiabilidade()`, `confiabilidade_compostos()`,
`validade_discriminante()`, `formatar_ajuste_apa()`, `checar_escore_somado()`, todas em
`scripts/psicometria.R` (a última usa também `est_comum.R`: carregar por `carregar.R`).
Modelo estrutural (regressões entre latentes, efeito indireto, multigrupo, interação
latente, curva de crescimento): `references/mee.md` (v1.5).

## Escolha do estimador

A regra de `sugerir_estimador()` vem de Rhemtulla, Brosseau-Liard e Savalei
(2012), simulação com 2 a 7 categorias, quatro tamanhos de amostra, duas
distribuições latentes e cinco padrões de limiares. O resumo do artigo diz:
cargas e erros-padrão saíram mais precisos com cat-LS (a família do WLSMV),
sobretudo com menos de 5 categorias; correlações entre fatores e ajuste saíram
igualmente bem com ML; cat-LS foi mais sensível a amostra pequena e a violação
da normalidade das variáveis latentes contínuas; ML foi mais sensível a
limiares assimétricos e ficou enviesado com cargas altas. Recomendação dos
autores: cat-LS com menos de 5 categorias, ML com 5 ou mais quando a amostra é
pequena e os limiares são aproximadamente simétricos, e qualquer um dos dois
com 6 ou 7 categorias em muitas condições.

Li (2016), outra simulação, aponta na mesma direção para as cargas: WLSMV menos
enviesado e mais preciso que MLR em quase todas as condições. Com amostra
pequena ou distribuição latente moderadamente não normal, WLSMV produziu
superestimação moderada das correlações entre fatores, e com N = 200 o MLR deu
erros-padrão melhores quando as latentes eram não normais.

Como a skill aplica isso:

| Situação | Estimador principal | Sensibilidade | De onde vem |
|---|---|---|---|
| Menos de 5 categorias | WLSMV, `ordered` | nenhuma | recomendação dos autores |
| 5 a 7 categorias, assimetria acima de 1 | WLSMV | MLR, com ressalva | ML enviesa com limiar assimétrico |
| 5 a 7 categorias, assimetria até 1, n até 200 | MLR | WLSMV | recomendação dos autores (ML com 5 ou mais categorias, amostra pequena, limiares simétricos) |
| 5 categorias, assimetria até 1, n acima de 200 | WLSMV | MLR | escolha da skill: cat-LS estima cargas com mais precisão |
| 6 ou 7 categorias, assimetria até 1, n acima de 200 | WLSMV | MLR | os dois aceitáveis segundo os autores; WLSMV por escolha da skill |
| 8 ou mais valores distintos (0-100, VAS) | MLR com FIML | nenhuma | escala quase contínua |

Duas operacionalizações são da skill, não do artigo: a assimetria do item como
aproximação da assimetria dos limiares, e "amostra pequena" como n até 200 (o
artigo não dá um número; Li, 2016, achou os problemas do WLSMV com N = 200). O
texto do manuscrito justifica o estimador pelo número de categorias, pela
distribuição e pelo n, citando Rhemtulla et al. (2012), e diz qual foi a
análise de sensibilidade.

Com n até 200 e WLSMV (menos de 5 categorias, ou limiares assimétricos), a
função acrescenta aviso: é o cenário em que as duas simulações acharam os
problemas do cat-LS.

## Código

```r
source("~/.local/share/psicometria-r/carregar.R")
itens <- paste0("item", 1:24)
aud <- auditar_itens(dados, itens)        # agregados por item
est <- sugerir_estimador(aud)             # lista: estimador, ordinal, motivo, sensibilidade
fit <- ajustar_afc(modelo, dados, itens, est$estimador, est$ordinal)
tabela_ajuste(fit, "modelo teórico")      # versões escalonada e robusta, rotuladas
```

`ajustar_afc()` decide sozinha entre ordinal e contínuo se `ordinal` não for
informado (até 7 valores distintos conta como ordinal, a mesma regra do
`cor_auto`), e o estimador padrão segue: WLSMV para ordinal, MLR para
contínuo. Com itens ordinais usa `parameterization = "theta"` (a mesma da
invariância, para os números baterem entre as etapas).

Ausentes. Com WLSMV o lavaan só oferece deleção listwise ou pairwise, e a
função usa pairwise, a mesma da `invariancia()`, para o n não mudar entre as
etapas. Deleção pairwise pressupõe ausência completamente ao acaso (MCAR); com
ausência ligada a outras variáveis do estudo (atrito seletivo, abandono por
gravidade de sintoma), FIML e imputação múltipla são as rotas defensáveis
(Enders & Baraldi, 2018, cap. 6 do Wiley Handbook, na pasta Psicometria-R). Para
itens ordinais com WLSMV, isso quer dizer imputação múltipla antes da AFC
(`mice`, com os resultados combinados por `semTools::runMI()` ou
`lavaan.mi`), dita no relato; com MLR, FIML (`missing = "fiml"`). O teste de MCAR
de Little não prova MCAR: não rejeitar não autoriza a pairwise (o curso inverte a
leitura da H0 numa aula; ver `fontes-e-decisoes.md`). Relatar a proporção de
ausentes por item (vem da auditoria) e a escolha feita.

Para modelo de segunda ordem (como o MPFI-24), a sintaxe lavaan é a de sempre:

```r
modelo <- "
  aceitacao =~ i1 + i2
  ...
  flexibilidade =~ aceitacao + presente + ...
  inflexibilidade =~ esquiva + ...
"
```

## Índices de ajuste

Relatar χ² escalonado com gl e p, CFI, TLI, RMSEA com IC 90% e SRMR (a lista
exigida pela fabrica-de-artigo). Com WLSMV e MLR o lavaan 0.7.3 dá duas
versões de CFI, TLI e RMSEA: a escalonada, tradicional, e a robusta, mais
nova. A documentação do `fitMeasures()` recomenda a robusta sempre que
existir. Com dado ordinal a robusta só é calculada se a matriz policórica for
positiva definida; caso contrário sai NA, o que segundo a mesma documentação
não é raro em amostra pequena. `tabela_ajuste()` avisa quando isso acontece.

Savalei (2021) mostra por que isso importa: as computações correntes dos
índices com dado categórico indicam ajuste bem melhor do que os mesmos dados
indicariam sem a categorização, o que leva a aceitar modelo ruim com mais
frequência. As computações propostas no artigo, desenvolvidas para DWLS,
ficaram mais próximas dos valores com dado contínuo em todas as condições
estudadas, exceto com itens binários e N = 200, onde todos os índices
categóricos foram mal.

Sobre os cortes. Hu e Bentler (1999) derivaram para o método ML os valores de
referência que a fabrica-de-artigo cita (CFI e TLI perto de .95, RMSEA perto
de .06, SRMR perto de .08). Xia e Yang (2019) compararam, para a mesma matriz
policórica e o mesmo modelo, os índices de ML com os de DWLS e ULS: DWLS e ULS
deram RMSEA menor e CFI e TLI maiores que ML em todas as condições, escalonados
ou não, de modo que aplicar os cortes convencionais com DWLS tende a não
descobrir desajuste. Na prática: com WLSMV, CFI de .96 não autoriza "bom
ajuste" sem ressalva. A saída é relatar os valores, preferir a versão robusta,
olhar o SRMR e os resíduos locais (`lavResiduals(fit)`, só correlações
residuais, que são agregados) e comparar modelos concorrentes.

McNeish e Wolf (2023) propõem cortes dinâmicos, derivados por simulação para o
modelo e os dados em análise, com um aplicativo Shiny que automatiza o
processo, e uma abordagem específica para modelo de um fator. O pacote R
`dynamic` não estava no CRAN em 10/10/2026 (conferido na instalação); usar o
aplicativo e conferir no artigo se o modelo e o tipo de dado em questão são
cobertos antes de citar.

## Cargas e solução

`lavaan::standardizedSolution(fit)` dá cargas padronizadas com erro-padrão e
IC. Conferir antes de interpretar:

- solução imprópria (variância negativa, correlação acima de 1): o lavaan
  avisa; modelo com caso Heywood não vai para o relato sem tratamento;
- carga abaixo de .40: o modo PSICOMETRICO da fabrica-de-artigo pede
  justificativa (a tabela de referência `statistical_reporting_standards.md`
  da mesma skill fala em .50, ideal .70; dizer qual critério foi usado);
- carga padronizada negativa: item invertido sem recodificar; `confiabilidade()`
  avisa, e alfa e ômega de pesos unitários saem errados até recodificar;
- fator com 2 itens: identificado só pelas correlações com outros fatores ou
  pelo fator de segunda ordem, e a confiabilidade dele fica limitada (caso do
  MPFI-24, com 12 fatores de 2 itens).

## Confiabilidade

`confiabilidade(fit, dados)` devolve, por fator, ômega e alfa.

O ômega vem de `semTools::compRelSEM()`. Com itens ordinais, o padrão
(`ord.scale = TRUE`) é a confiabilidade do escore somado observado, o método
de Green e Yang (2009), que admite relação não linear entre fator e item. É
esse que vai para o relato, porque é o escore que o usuário do instrumento vai
somar. `compRelSEM(fit, ord.scale = FALSE)` dá a confiabilidade das respostas
latentes contínuas, que é maior (os evals mostram isso) e não descreve o
escore observado; não relatar como se fosse.

O alfa é calculado nos escores observados e vai junto porque a fabrica-de-artigo
exige alfa e ômega lado a lado (Dunn et al., 2014; McNeish, 2018). McNeish
(2018) resume por que o alfa costuma subestimar a confiabilidade quando os
pressupostos dele não valem e apresenta alternativas (ômega total, ômega total
de Revelle, maior limite inferior, coeficiente H).

Intervalo de confiança (v1.6): `confiabilidade(fit, dados, ic = TRUE, reamostras =
1000)` acrescenta IC por bootstrap percentílico dos casos para ômega e alfa, reajustando
o mesmo modelo (estimador, ordinal e ausentes iguais) em cada reamostra, a rota que
Kelley e Pornprasertmanit (2016) recomendam para o ômega; Martinková e Hladká (2023,
4.4.4) fazem o mesmo para a confiabilidade por componentes de variância. Reamostra que
não converge sai e é contada (mensagem; acima de 5% o IC fica otimista). Com WLSMV leva
minutos. Nos evals, o IC do ômega bate com o IC delta do parâmetro definido no lavaan e
o do alfa com o de Feldt (`psych::alpha.ci`), as duas com tolerância de ,015. Relatar
como "ω = ,82, IC 95% [,80; ,84]" e dizer que o IC é por bootstrap.

Modelo de segunda ordem. No semTools 0.5.10 os argumentos `higher=` e
`omit.factors=` do `compRelSEM()` estão obsoletos (o pacote avisa); fatores de
ordem superior são detectados sozinhos e tratados como a única fonte de escore
verdadeiro, com os resíduos dos fatores de primeira ordem contando como erro.
Duas consequências, ambas conferidas nos evals contra valores teóricos:

- `confiabilidade()` num modelo de segunda ordem devolve, para cada
  subescala, só a parte da variância devida ao fator geral (nos evals, perto
  de .30 para subescalas de 2 itens), e não a confiabilidade da subescala.
  Para a subescala, ajustar o modelo de fatores correlacionados e rodar
  `confiabilidade()` nele. A função avisa quando detecta fator de ordem
  superior.
- `confiabilidade_compostos(fit, list(flex = itens_flex, inflex = itens_inflex))`
  dá a confiabilidade do escore global. Rodado no modelo de segunda ordem, o
  valor é um ômega hierárquico (só o fator geral como escore verdadeiro;
  teórico .766 no eval, estimado dentro de .04); rodado no modelo de fatores
  correlacionados, é o ômega total do composto (teórico .899, estimado dentro
  de .03). Relatar os dois, com o nome certo.
- Pela documentação e pelo código do semTools, `compRelSEM(fit, true =
  list(s1 = "s1"))` daria o ômega da subescala dentro do próprio modelo de
  segunda ordem. A skill não testou esse caminho; até testar, usar o modelo de
  fatores correlacionados.

`confiabilidade()` é para modelo de um grupo (com modelo multigrupo ela para e
pede para rodar por grupo).

## Escore somado

`checar_escore_somado(dados, itens, ordinal = NULL, estimador = NULL, alfa = .05,
decimal = ".", idioma = "pt")`

Somar (ou tirar a média dos) itens é pontuar com um modelo. McNeish e Wolf (2020)
mostram que a soma corresponde ao escore fatorial do modelo paralelo: cargas não
padronizadas iguais (não precisam valer 1), variâncias residuais iguais, interceptos dos
itens livres. A AFC usual é congenérica, com cargas e resíduos livres. Validar com o
congenérico e pontuar com a soma é o que os autores chamam de "bait-and-switch": a
evidência vem de um modelo e o escore, de outro que ninguém testou. A função testa esse
outro modelo.

- `itens`: vetor (uma escala) ou lista nomeada (uma entrada por subescala). Cada
  subescala vai num modelo de um fator próprio. No exemplo do artigo as subescalas
  entram juntas, com erros iguais dentro de cada uma e sem covariância entre os fatores,
  porque com covariância o escore fatorial deixa de ser transformação linear da soma de
  cada subescala. Testar uma a uma dá a decisão por subescala; a correlação entre
  subescalas fica para a AFC de validade. Para escore total, todos os itens num vetor.
- Itens contínuos (MLR por padrão; `estimador = "ML"` reproduz o artigo), 3 itens ou
  mais: congenérico, tau-equivalente (cargas iguais) e paralelo (cargas e resíduos
  iguais), com três comparações. A decisiva é paralelo contra congenérico, a do artigo.
  As intermediárias mostram qual restrição falha; essa sequência é da skill, o artigo
  compara paralelo e congenérico direto.
- 2 itens: o congenérico não é identificado, e sai só paralelo contra tau-equivalente
  (resíduos iguais, 1 gl).
- Itens ordinais (até 7 valores distintos, a mesma regra de `.pr_e_ordinal()`, ou
  `ordinal = TRUE`): WLSMV com parametrização theta, em que as variâncias residuais das
  respostas latentes ficam fixas em 1. Cargas iguais contra livres já testa o paralelo
  dessas respostas latentes; a soma das categorias observadas ainda as trata como
  igualmente espaçadas, o que o teste não alcança. Esse caso é extensão da skill: o
  artigo trata os itens como contínuos. Com 2 itens ordinais não há teste, e a frase diz
  isso sem afirmar nada sobre a soma. Itens com poucas categorias tratados como
  contínuos: `ordinal = FALSE` ou um estimador da família ML.
- Diferença de qui-quadrado pelo `lavTestLRT()`, com o método lido do cabeçalho que o
  próprio lavaan imprime: ML, padrão; MLR, escalonada (Satorra & Bentler, 2001); WLSMV,
  escalonada e deslocada (Satorra, 2000). Diferença escalonada negativa é refeita com
  Satorra e Bentler (2010), com registro. O cabeçalho usa aspas curvas em locale UTF-8
  (`dQuote`) e retas em C; a leitura aceita as duas (bug achado na v1.3, coberto nos
  evals).
- Saída (`pr_resultado`): `$tabela` (uma linha por comparação: Δχ², Δgl, p, método,
  resultado, coluna `decisiva`), `$ajuste` (χ², gl, p, CFI, TLI, RMSEA com IC 90% e SRMR
  de cada modelo, na mesma versão em todos os modelos da escala), `$cargas` (não
  padronizadas, padronizadas e resíduos do modelo menos restrito), `$amplitude` (maior
  menos menor carga padronizada), `$soma_rejeitada` (por escala), `$fits`,
  `$metodo_lrt`, `$ordinal`, `$estimador`, `$n_por_escala`, `$decisoes`, `$avisos`,
  `$apa` (PT ou EN, `decimal = ","`). Linha sem nenhum item da escala sai
  (`$n_excluidos`); o FIML usa as linhas parciais. Item constante ou escala de 1 item
  param com mensagem.

Leitura do resultado:

- Rejeitar o paralelo não obriga a abandonar a soma em qualquer caso. O próprio artigo
  diz que um congenérico com cargas padronizadas quase iguais pode ser aproximado pelo
  paralelo, e por isso a função mostra o teste e a amplitude das cargas lado a lado. A
  decisão vai para o relato com os dois números. Quando a soma for rejeitada e a
  diferença importar, a saída indicada é o modelo de mensuração dentro de uma SEM
  (abordagem simultânea) ou o escore fatorial com a correção de Croon; para o escore
  fatorial, Bartlett quando ele for desfecho e regressão quando for preditor, como
  McNeish e Wolf (2020) resumem de Skrondal e Laake (2001). No lavaan 0.7.3,
  `lavPredict(fit, method = "Bartlett")` ou `"regression"` (itens contínuos) e `"EBM"`
  ou `"ML"` (ordinais); a função `fsr()` que o artigo cita para a correção de Croon não
  existe mais nessa versão (conferido).
- Não rejeitar não mostra que as restrições valem: com amostra pequena o teste tem pouco
  poder (Greenland et al., 2016). Conferir o ajuste do paralelo e a amplitude.
- Correlação alta entre soma e escore fatorial não substitui o teste. No exemplo do
  artigo, a soma da subescala Speeded e o escore de Bartlett correlacionaram .985
  (R² = .970), e a conclusão sobre a diferença entre as escolas mudou de direção e de
  significância conforme o escore usado.

Âncora (evals, `evals_simulados.R`): o exemplo do artigo, com o banco
HolzingerSwineford1939 do lavaan, itens x4 a x9 arredondados para o inteiro mais
próximo, se reproduz no R. O arredondamento tem de ser meio para cima
(`floor(x + 0.5)`); o `round()` do R leva o meio ao par, e a média de x5 sai 4,35, longe
dos 4,46 da Tabela 1 do artigo, o que muda os qui-quadrados. Com ML: congenérico de um fator χ²(9) = 115,37,
CFI = .81, RMSEA = .20 [.17, .23]; cargas padronizadas .82, .85, .79, .17, .18, .26
(Tabela 2 do artigo); paralelo de um fator CFI = .45, RMSEA = .23 [.21, .25], SRMR = .24
(o SRMR do artigo é o do Mplus, `srmr_mplus` no lavaan); dois fatores congenérico
χ²(8) = 14,74; dois fatores paralelo sem covariância χ²(17) = 55,54, diferença χ²(9) =
40,80. O χ²(19) = 5361,86 que o texto dá para o paralelo de um fator não é compatível
com o CFI e o RMSEA que o próprio artigo dá para esse modelo; o lavaan dá 325,90, que
reproduz os dois. Na função, por subescala: Verbal tem cargas padronizadas quase iguais
(amplitude .05) e mesmo assim o paralelo é rejeitado, Δχ²(4) = 17,18, p = .002; Speeded,
Δχ²(4) = 10,69, p = .030.

## Modelo bifator (v1.4)

Quando os fatores de primeira ordem se correlacionam muito, a pergunta é se o
escore total (fator geral) basta ou se as subescalas medem algo além dele. O
bifator põe um fator geral em todos os itens e fatores específicos ortogonais;
no ajuste ele quase sempre ganha (absorve mais variância), então a decisão sai
dos índices (Rodriguez et al., 2016a; Reise et al., 2018, cap. 22 do Wiley
Handbook, na pasta). O curso ensina o bifator na AFC (ômega hierárquico pelas
fórmulas de Rios e Wells, numa planilha) e lê o ECV > ,85 do FACTOR como critério
de unidimensionalidade.

```r
m_bif <- "
  G  =~ i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9
  S1 =~ i1 + i2 + i3
  S2 =~ i4 + i5 + i6
  S3 =~ i7 + i8 + i9"
fit_bif <- ajustar_afc(m_bif, dados, itens, orthogonal = TRUE, std.lv = TRUE)
fit_uni <- ajustar_afc("G =~ i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9", dados, itens, std.lv = TRUE)
indices_bifator(fit_bif, unidimensional = fit_uni)
```

O que sai: ωT e ωH do escore total; por subescala, ωS e ωHS; ECV (geral), ECV de
cada específico dentro da subescala e em relação ao total; I-ECV por item; PUC; H
de Hancock e Mueller (2001) e determinação fatorial por fator; ARPB, quando o
modelo de um fator entra em `unidimensional`. Variância do escore somado pela
matriz implicada pelo modelo. Iguais ao `BifactorIndicesCalculator` nos evals.

Réguas, declaradas na saída com a fonte: ωH ≥ ,80 para ler o escore total como
medida do fator geral (Rodriguez et al., 2016a); ECV > ,70 com PUC > ,70 para
tratar a variância comum como essencialmente unidimensional (Rodriguez et al.,
2016a); com PUC < ,80, ECV > ,60 e ωH > ,70 (Reise et al., 2013); H ≥ ,70 para
fator bem definido (Hancock & Mueller, 2001); determinação ≥ ,90 para usar escore
fatorial (Gorsuch, 1983); ARPB abaixo de 10 a 15% (Muthén et al., 1987). As
réguas de ECV e PUC chegaram à skill por fonte secundária (ver a âncora em
`referencias-verificadas.md`): ler Rodriguez et al. antes de citar número de
corte num manuscrito. ωHS não tem régua fixa: compará-lo com o ωS mostra quanto
da confiabilidade da subescala vem emprestada do fator geral.

Itens ordinais (WLSMV): as cargas são das respostas latentes, e o ω aqui é o da
soma das respostas latentes. O ω da soma observada das categorias (Green & Yang,
2009) sai do `BifactorIndicesCalculator::bifactorIndices(fit)`.

Bifator ou segunda ordem: o de segunda ordem é aninhado no bifator com
restrições de proporcionalidade; os dois respondem a perguntas parecidas, mas só
o bifator separa a variância específica de cada subescala. Para o MPFI-24, ver
`references/pbat-mpfi.md`.

Outros modelos que a pasta discute sem função própria: MTMM (traço-método,
cap. do Wiley Handbook sobre validade) e interação latente; usar lavaan direto,
com a documentação aberta, e dizer no relato que a rota saiu do escopo testado.

## Validade convergente e discriminante

`validade_discriminante(fit, dados)` (modelo de um grupo) devolve AVE, raiz da AVE, correlações
entre fatores, HTMT e o critério de Fornell-Larcker (raiz da AVE de cada fator
maior que as correlações dele), sempre entre fatores de primeira ordem
alinhados por nome. Com um fator só, o critério sai NA (não há o que
comparar). Em modelo de segunda ordem, o HTMT é calculado sobre a parte de
primeira ordem do modelo. Os evals conferem os três casos e um par de fatores
correlacionados a .90, em que HTMT e Fornell-Larcker reprovam.

- AVE de .50 para cima é a convenção que a fabrica-de-artigo usa.
- HTMT abaixo de .85 (Henseler et al., 2015). O HTMT foi proposto e testado no
  contexto de SEM baseada em variância (PLS); na simulação dos autores,
  Fornell-Larcker e cargas cruzadas falharam em detectar falta de validade
  discriminante em situações comuns, e o HTMT foi melhor. O uso em AFC é
  prática corrente, e o `semTools::htmt()` aceita itens ordinais.
- Com fatores de segunda ordem e fatores de 2 itens, AVE e HTMT dos fatores de
  primeira ordem dizem pouco; relatar com essa ressalva.
