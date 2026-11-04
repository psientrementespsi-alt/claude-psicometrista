# Redes psicométricas, EGA e comparação de redes

## Sumário

- 12 (linha 29)
- 26 (linha 43)
- 86 (linha 103)
- 129 (linha 146)
- 145 (linha 162)
- 153 (linha 170)
- 168 (linha 185)
- 199 (linha 216)
- 238 (linha 255)
- 252 (linha 269)
- 289 (linha 306)
- 309 (linha 326)
- 342 (linha 359)
- 352 (linha 369)

Índice: 1. Quando usar | 2. Estimar: `rede_psicometrica()` | 3. Precisão e estabilidade:
`estabilidade_rede()` | 4. Ponte e previsibilidade | 5. Tamanho de amostra | 6. Como relatar |
7. NCT | 8. EGA e bootEGA | 9. Boruta | 10. Rede não é modelo fatorial | 11. Armadilhas |
12. Fontes e divergências

Pacotes: bootnet, qgraph, networktools (ponte), mgm (previsibilidade, dados mistos), IsingFit
(binários), EGAnet (UVA, EGA, bootEGA), NetworkComparisonTest, Boruta. Funções em
`scripts/est_redes.R`; testes em `scripts/evals_redes.R`.

## 1. Quando usar

A pergunta pede rede quando o foco está em elementos nomeados (sintomas, itens, processos), nas
relações diretas entre eles, em comorbidade ou em pontes, e quando se aceita
o pressuposto de que os elementos interagem diretamente; quem pensa em variável latente
reformula a pergunta como correlação entre fatores. Rede transversal descreve
associação, não causa.

| Dado dos nós | Modelo | `estimador` |
|---|---|---|
| ordinal (Likert), contínuo, 0-100 | GGM: correlações parciais regularizadas | `"EBICglasso"` (padrão) |
| binário (0/1) | Ising | `"IsingFit"` |
| misto (contínuo, categórico, contagem) | modelo gráfico misto | `"mgm"` (com `tipos`) |

## 2. Estimar a rede: `rede_psicometrica()`

```r
res <- rede_psicometrica(dados, itens, estimador = "EBICglasso", tuning = 0.5, cor = "cor_auto",
                         ausentes = "pairwise", tamanho_amostra = "maximum", tipos = NULL,
                         redundancia = TRUE, n_fortes = 10, decimal = ".", idioma = "pt")
res; relatar(res)                       # tabelas, decisões, avisos e frase APA
res$pesos; res$centralidade_bruta; res$centralidade_z
plot(res$modelo, layout = "spring")     # figura (bootnet/qgraph)
```

Padrões = módulo Network do JASP 0.95.1, que o livro usa (; `jaspNetwork`, A.5 de
`pesquisa_github.md`). Os do bootnet são outros, então tudo vai explícito:

| Argumento | Skill (JASP) | bootnet 1.9.1 | Fonte |
|---|---|---|---|
| `tuning` (γ do EBIC) | 0.5 | 0.5 (0.25 no IsingFit e no mgm) |; Epskamp & Fried (2018) |
| `cor` | `"cor_auto"` | `"cor"` (Pearson) | policórica para ordinal |
| `ausentes` / `tamanho_amostra` | `"pairwise"` / `"maximum"` | `"pairwise"` / `"pairwise_average"` | |
| centralidade | escores z ("Normalized") | bruta | |

O que a função decide sozinha:
- `cor_auto` trata item com até 7 categorias como ordinal (policórica) e acima disso como
  contínuo (Pearson). Em 0-100 (PBAT) a rede fica em Pearson, a mesma regra descrita na
  adaptação romena do PBAT (Popa et al., 2025).
- `"maximum"`: o N do EBIC é o número de casos com ao menos uma resposta. Com ausentes,
  `"maximum"` e `"pairwise_average"` podem dar redes diferentes; declarar. Caso sem resposta
  nenhuma sai e conta em `n_excluidos`. IsingFit e mgm só aceitam listwise (declarado); item
  binário 1/2 vira 0/1; no mgm sem `tipos`, item com 2 valores vira `"c"` e o resto `"g"`.
- Centralidade igual à do JASP: `qgraph::centrality(W, weighted = TRUE, signed = TRUE,
  all.shortest.paths = FALSE)`, intermediação dividida por 2 e z = (x - média) / DP amostral.
- `redundancia = TRUE` roda `EGAnet::UVA()` antes e avisa pares com wTO ≥ .25. Aresta muito
  forte entre itens vizinhos pode ser redação parecida, distinção que a rede sozinha não
  resolve: decidir antes se mantém, combina ou remove o item, e declarar.
- O aviso de rede densa do qgraph sai traduzido. `threshold = TRUE` (em `...`) ganha
  especificidade à custa de sensibilidade; Isvoranu e Epskamp (2023), com 564 mil bases
  simuladas, mostram que essa troca sempre existe e que o método se escolhe pela pergunta.
  Decidir antes de olhar o resultado; não trocar de método até a rede ficar bonita.

Saída: `tabela` (nós, arestas possíveis p(p-1)/2, não nulas em número exato, positivas,
negativas, densidade, esparsidade, N, incompletos, γ), arestas mais fortes, centralidade em z,
`$checklist` e a frase APA. No banco BFI Network do livro (2800 casos, 25 itens), que a
função reproduz até a terceira casa do JASP:

> A rede dos 25 itens foi estimada como modelo gráfico gaussiano (arestas = correlações
> parciais), com regularização LASSO e seleção do modelo pelo EBIC (EBICglasso; γ = 0.50), a
> partir de correlações policóricas (cor_auto), com exclusão par a par de ausentes e tamanho de
> amostra máximo (N = 2800). Das 300 arestas possíveis, 181 foram estimadas como não nulas
> (densidade = .60, 128 positivas e 53 negativas). A aresta de maior peso ligou N1 e N2 (peso
> = .57), seguida de N3–N4 (peso = .30). Na força (escores z), os maiores valores foram de C4
> (1.65), E2 (1.39), E4 (1.30), e o menor foi de C3 (-1.71). Os resultados descrevem
> associações condicionais em dados transversais e não autorizam inferências causais ou de
> intervenção.

Índices: força é o principal e o mais estável; influência esperada quando aresta
negativa representa oposição; proximidade e intermediação não medem importância do
nó em redes psicológicas (Bringmann et al., 2019; Hallquist et al., 2021) e só entram com CS
adequado. A ordem de centralidade só se interpreta depois da seção 3. Posição na figura não é
centralidade: o layout é algorítmico.

## 3. Precisão e estabilidade: `estabilidade_rede()`

Obrigatório antes de interpretar; são duas perguntas, e é comum relatar uma e esquecer a outra
.

```r
est <- estabilidade_rede(res, n_boot = 1000, indices = c("strength", "expectedInfluence"),
                         semente = 2026, n_boot_casos = 1000, caso_min = 0.05, caso_max = 0.75,
                         caso_n = 10, executar = c("arestas", "casos"), comunidades = NULL)
est$arestas_ic; est$comparacao_arestas; est$cs
plot(est$boot_arestas, order = "sample", labels = FALSE)   # figura do suplementar
```

Precisão das arestas: bootstrap não paramétrico (participantes reamostrados com reposição),
1000 reamostras, semente fixada antes de cada bootstrap, IC 95% por percentis
(quantil tipo 6, o do `summary` do bootnet).
- IC de rede regularizada é descritivo da variabilidade, não teste: o LASSO encolhe pesos para
  zero, e o IC incluir o zero não prova que a aresta é nula.
- IC sobrepostos: não se pode dizer que uma aresta é mais forte. A função compara a
  mais forte com a segunda pela diferença de |peso| nas reamostras (IC 95%; igual ao
  `differenceTest` quando as duas são positivas) e conta as outras com IC sobreposto ao dela.
  Outros pares: `bootnet::differenceTest(est$boot_arestas, x = "A", x2 = "B", y = "C", y2 =
  "D", measure = "edge")`; força: `differenceTest(est$boot_arestas, "A", "B", "strength")`.
  Sem correção para comparações múltiplas.
- Ordinal pede o não paramétrico: o paramétrico gera dados da rede estimada, supondo que ela
  é a verdadeira, e amostra dado contínuo normal.

Estabilidade da centralidade: bootstrap de descarte de casos com `caseMin = 0.05`, `caseMax =
0.75`, `caseN = 10` declarados (padrões do bootnet; o livro só ilustra "10%, 20%").
CS-coefficient: maior proporção de casos que se pode descartar mantendo, com 95% de
probabilidade, correlação de pelo menos .70 com a centralidade original (; Epskamp,
Borsboom & Fried, 2018). Cortes: não abaixo de .25 e, de preferência, acima de .50, que os
autores chamam de um tanto arbitrários. Não há patamar de .70: .70 é a correlação da definição
. CS de .30 fica acima do mínimo e abaixo do desejável; CS baixo num índice
não invalida a rede. O CS só assume valores da grade testada, com teto em
`caso_max` (a função marca). Proximidade e intermediação entram em `indices` quando o texto
as citar; índices de ponte (`"bridgeStrength"` etc.) exigem `comunidades`.

Reamostras: os evals usam de 20 a 100 por velocidade; para relato, 1000 ou mais (a função
avisa abaixo disso). Com `n_nucleos > 1`, conferir a reprodutibilidade com a semente. A frase
gerada traz o IC da aresta mais forte, a diferença com a segunda, quantas têm IC sobreposto ao
dela e o CS de cada índice com a leitura do corte.

## 4. Ponte e previsibilidade

`ponte(res, comunidades)`: índices de Jones, Ma e McNally (2021) pelo `networktools::bridge`
(força de ponte, influência esperada de ponte de 1 e 2 passos, intermediação e proximidade de
ponte), brutos e z da força de ponte. Comunidades definidas antes do resultado
(teoria, subescalas ou EGA declarado), como vetor na ordem dos nós ou `list(fator = itens)`.
Sem corte para "nó ponte": conferir no artigo antes de adotar um. Estabilidade:
`estabilidade_rede(res, indices = "bridgeStrength", comunidades = ...)`.

`previsibilidade(res, tuning = NULL)`: proporção da variância do nó explicada pelos vizinhos,
"uma espécie de R² do nó" (Haslbeck & Waldorp, 2018), pelo mgm (EBIC com o γ da
rede, regra AND, casos completos). Nó gaussiano: R²; categórico: acurácia (CC) e acurácia
normalizada além da categoria modal (nCC). Em rede EBICglasso sai também `r2_ggm_implicito` =
1 - 1/diagonal da precisão regularizada, coerente com a figura e encolhido pelo LASSO. Nó
central com previsibilidade baixa sugere determinantes fora do modelo.

## 5. Tamanho de amostra

Arestas possíveis = p(p-1)/2: 10 nós, 45; 20, 190; 25, 300; 30, 435. Não existe
número mágico nem razão fixa por parâmetro. Estratégias: comparar com
estudos prévios (mesmo instrumento, nós e estimador), preferir versões curtas e simular antes
da coleta com `bootnet::netSimulator` ou com o pacote powerly (Constantin, Schuurman
& Vermunt, 2026). A função relata N, nós e arestas possíveis juntos e não julga o N por corte.

## 6. Como relatar

Elementos do parágrafo-modelo e do checklist do livro, que as funções
montam em `$checklist`: (1) estrutura antes dos números (GGM, EBICglasso, γ, correlação,
ausentes, N); (2) número exato de arestas não nulas sobre as possíveis (181 de 300, não
"cerca de metade"); (3) aresta mais forte com peso e IC; (4) arestas menores com IC
sobrepostos ditas indistinguíveis; (5) CS de cada índice citado, interpretando só o estável;
(6) N em relação aos nós; (7) frase de alcance: associações condicionais em dados
transversais, sem inferência causal ou de intervenção.

Formato: peso de aresta GGM e CS sem zero à esquerda (.57; .25); γ e z com zero
(0.50; 1.65); com `decimal = ","`, vírgula decimal e ponto e vírgula entre itens. Pesos de
Ising e mgm não são correlações e saem com zero. Gráficos de bootstrap e versão colorida da
rede vão para o suplementar.

## 7. Comparar redes entre grupos (NCT)

```r
library(NetworkComparisonTest)
set.seed(2026)
nct <- NCT(dados[g == "A", itens], dados[g == "B", itens], it = 1000,
           test.edges = TRUE, edges = "all", p.adjust.methods = "holm")
nct$nwinv.pval    # estrutura (diferença máxima de aresta)
nct$glstrinv.pval # força global
```

O NCT (van Borkulo et al., 2023) testa por permutação, entre duas amostras independentes, a
estrutura da rede, a força de cada aresta e a força global. Comparar centralidade entre grupos
exige teste formal, não dois gráficos lado a lado; o split do JASP só compara a olho.

Cuidado com o estimador. Com dados brutos, o NCT 2.2.4 estima as redes com correlação de
Pearson por padrão (o código do pacote, `NCT_estimator_GGM`, usa `cor()`), e o próprio pacote
comenta que o uso de `cor_auto` com o NCT não foi validado. Com itens ordinais, a rede testada
pelo NCT é então outra que a rede policórica relatada pelo bootnet: declarar Pearson no texto
(ou passar objetos do `estimateNetwork` ao NCT, conferindo na ajuda) e tratar ausentes antes,
porque `cor()` com NA faz a análise falhar. No PBAT, em 0-100, as duas redes já são de
Pearson. Na simulação dos autores, o erro tipo I ficou perto do nominal e o poder foi
suficiente quando o tamanho da amostra e a diferença entre as redes eram substanciais.

Leitura correta de p > .05 no NCT: ausência de evidência de diferença, não demonstração de que
as redes são iguais. Com amostra modesta ou diferença pequena o teste não tem poder para achar
a diferença. Relatar assim, com o tamanho de cada grupo. (A adaptação romena do PBAT lê o NCT
não significativo como evidência de que as redes clínica e não clínica são invariantes e
justifica o poder pelo N total de 637; o grupo clínico, porém, tinha 152 casos e a rede dele
teve CS de .20, o que pede a ressalva; Popa et al., 2025.)

## 8. EGA e bootEGA (dimensionalidade)

```r
library(EGAnet)
ega <- EGA(dados[itens])
set.seed(2026)
boot <- bootEGA(dados[itens], iter = 500, type = "resampling", seed = 2026)
boot$summary.table                 # mediana e frequência do número de dimensões
dimensionStability(boot)           # consistência estrutural por dimensão
itemStability(boot)                # estabilidade de cada item
```

EGA estima o número de dimensões com rede regularizada e algoritmo de comunidades. Golino e
Epskamp (2017) acharam desempenho comparável à análise paralela em várias condições e superior
na estrutura de quatro fatores com correlação .70 entre eles; Golino et al. (2020), numa
simulação mais ampla, acharam EGA e EGAtmfg tão precisos quanto a análise paralela, o melhor
dos métodos tradicionais, e com as melhores propriedades em amostra grande.

O bootEGA (Christensen & Golino, 2021) gera a distribuição amostral do EGA: mediana e
frequência do número de dimensões, consistência estrutural (com que frequência cada dimensão se
repete exatamente) e estabilidade de item, que separa item mal alocado, item multidimensional e
item redundante. O artigo deriva por simulação diretrizes para estabilidade aceitável de item:
na discussão, estabilidade de item de .70 ou mais; na seção de resultados, os pontos ótimos
variaram de .651 (GLASSO com reamostragem) a .721 (TMFG paramétrico), e valores entre .65 e
.75 marcam onde os itens começam a perder estabilidade. Com GLASSO e reamostragem, .70 é
conservador. (Âncora lida na página do artigo na MDPI por um revisor independente em
10/10/2026; conferir a frase no PDF antes de citação literal.) O corte de "mais de 50% das
réplicas" das notas do curso não é adotado (seção 12).

`type`: o padrão do EGAnet é `"parametric"` (gera dados normais multivariados a partir da
matriz de correlação); `"resampling"` reamostra os casos. Christensen e Golino (2021, seção 4)
recomendam a reamostragem quando a distribuição subjacente é desconhecida ou não normal, que é
o caso de item ordinal; o mesmo motivo leva Epskamp et al. (2018) a recomendar o bootstrap não
paramétrico no bootnet.

Redundância de itens: `UVA()` do EGAnet identifica pares de itens localmente dependentes antes
do EGA e antes da rede de itens (`rede_psicometrica()` já roda). Itens de conteúdo quase
idêntico (comuns em versões curtas e traduções) puxam dimensões espúrias e arestas infladas.

## 9. Boruta (relevância de item para desfecho)

```r
library(Boruta)
set.seed(2026)
br <- Boruta(phq9_total ~ ., data = dados[c("phq9_total", itens)], maxRuns = 500)
attStats(br)
```

Boruta (Kursa & Rudnicki, 2010) seleciona as variáveis relevantes para prever um desfecho
comparando a importância de cada uma, num random forest, com cópias embaralhadas delas. Não é
análise de estrutura: diz quais itens se relacionam com o desfecho no grupo, não como os itens
se organizam. É a análise central da literatura do PBAT (ver `references/pbat-mpfi.md`).

## 10. Rede não é modelo fatorial (e vice-versa)

Rede regularizada entre itens e AFC respondem perguntas diferentes, e a escolha deve vir do
modelo teórico (variável latente comum ou interação entre componentes), não do resultado que
agrada. Os dois paradigmas são lentes complementares, e tratar a rede como refutação da
psicometria clássica é erro conceitual. No framework do usuário (network approach), a
rede é o modelo de primeira linha para sintomas e processos quando a pergunta é sobre relações
diretas entre eles; para validade estrutural de instrumento com escore somado, a AFC continua
sendo a evidência que periódico e banca esperam. Quando as duas entram no mesmo estudo, dizer
o papel de cada uma. Com itens de um mesmo fator latente, a força pode só repetir a carga
fatorial (Hallquist et al., 2021).

**Por que o ajuste não decide.** Os dois tipos de modelo podem gerar os mesmos dados. Com
itens binários, Kruis e Maris (2016) mostram que o modelo de Ising (rede), um modelo de TRI
multidimensional (causa comum, com tantas dimensões latentes quanto os autovalores não nulos
da matriz de interações; o Rasch é o caso mais simples) e um modelo de seleção por colisor
(efeito comum) são matematicamente equivalentes: o mesmo resultado admite as três leituras.
van Bork et al. (2021) generalizam: toda rede tem um modelo de variáveis latentes
equivalente, e vice-versa, mas esses pares raramente são a comparação que interessa (uma rede
com arestas zeradas equivale a um modelo de p − 1 fatores com restrições de
proporcionalidade nas cargas, sem leitura substantiva), e os modelos que saem de hipóteses
de interesse em geral não são equivalentes. Para a skill, isso quer dizer que AFC com bom
ajuste não refuta a rede, rede estável não refuta o fator, e comparar o ajuste dos dois não
escolhe o mecanismo.

**Quando os dados separam.** Fator único e rede esparsa (com arestas zero) não têm membro
equivalente em comum, e o fator único implica três coisas que a rede esparsa não implica
(van Bork et al., 2021): na população, nenhuma correlação parcial é zero (salvo carga
padronizada 0 ou ±1), nenhuma tem sinal oposto ao da correlação bivariada e nenhuma é mais
forte que ela. O teste PCL (partial correlation likelihood) dos autores compara a proporção
de parciais que trocam de sinal ou crescem em módulo com as distribuições simuladas sob o
fator único e sob a rede esparsa de melhor ajuste; o código está no Apêndice A do artigo.
A skill não implementa o PCL. Contar essas parciais na amostra não é teste, porque o acaso
produz algumas mesmo com fator único verdadeiro. Quando a pergunta for "fator único ou rede
esparsa", citar o teste, rodá-lo à parte e dizer que a decisão vale para esses dois modelos,
não para fator contra rede em geral.

## 10b. Rede bayesiana (modelo gráfico probabilístico dirigido, v1.5)

O módulo "Probabilistic Graphical Models" do curso (scripts com `bnlearn` e `pcalg`)
aprende um grafo acíclico dirigido dos dados. `rede_bayesiana(dados, vars, algoritmo)`
faz isso com o bnlearn (Scutari, 2010): hc (padrão), tabu, pc.stable ou mmhc, e, em vez
de uma rede única (o que as aulas fazem), a rede média de `reamostras` bootstrap
(`boot.strength` + `averaged.network`; limiar de Scutari & Nagarajan, 2013). Cada aresta
sai com força (proporção de redes em que aparece) e direção (proporção com aquele
sentido). Variáveis com até 7 valores viram rede discreta (como o código das aulas, que
converte inteiros em fator); contínuas, gaussiana; misturar os dois tipos pede
discretizar. `comparar_algoritmos = TRUE` faz a validação cruzada das aulas (`bn.cv`) para
os quatro algoritmos.
- Direção de aresta aprendida de dados transversais não é efeito causal: só a classe de
  equivalência é identificada, e arestas com direção perto de 0,5 têm sentido
  indeterminado. A aula "PGMs versus MEE" compara com a MEE: na MEE a estrutura vem da
  teoria e é testada; aqui ela é estimada dos dados e é hipótese.
- Algoritmo PC e efeitos de intervenção do `pcalg` (`pc`, `ida`) ficam fora das funções.
- Evals: no banco `gaussian.test` do bnlearn (o da aula), o esqueleto da rede média é o da
  rede geradora.

## 10c. Redes temporais de EMA e diário (v1.7)

`rede_temporal(dados, vars, id = NULL, tempo = NULL, dia = NULL, bip = NULL,
destendenciar = TRUE, gamma = .5, aleatorios = "orthogonal")`, pelo modelo GVAR de lag 1
dos capítulos 9 e 10 de Isvoranu et al. (2022). Dados longos, uma linha por medida, já
ordenáveis por pessoa, dia e bipe.

- **Uma pessoa** (sem `id`): graphicalVAR (Abegaz & Wit, 2013), regularizado com seleção
  por EBIC. Rede temporal em correlações parciais dirigidas (PDC, linha = t - 1, coluna =
  t, com as autorregressões) e contemporânea em correlações parciais (PCC), o que fica
  correlacionado no mesmo bipe depois de tirar o temporal. `dia` corta os pares que cruzam
  a noite e `bip` respeita medidas faltantes, como no tutorial 10.1 do livro (não usar
  `dia` com uma medida por dia: o pacote apaga tudo). Aviso abaixo de 100 pontos no tempo
  (Epskamp, Waldorp et al., 2018, recuperaram redes esparsas de 8 nós com 100).
- **Várias pessoas** (`id`): mlVAR em dois passos (lme4). Tabela temporal com efeito
  fixo, EP, p e DP do efeito aleatório (quanto a aresta varia entre pessoas); redes
  contemporânea e entre pessoas nas tabelas extras. Aviso com menos de 20 medidas por
  pessoa (a centralização na média da pessoa pede pelo menos cerca de 20, cap. 10) e com
  menos de 30 pessoas. `aleatorios = "orthogonal"` é o padrão (rápido, recomendado com
  muitas variáveis); `"correlated"` estima as covariâncias dos aleatórios.
- **Destendenciar** (padrão): tira a tendência linear no tempo de cada variável, por
  pessoa. Epskamp, van Borkulo et al. (2018) mostraram que destendenciar tudo ou só o
  significativo dá o mesmo e que não destendenciar baixa a especificidade; quando a
  mudança no tempo é a pergunta (efeito de intervenção), `destendenciar = FALSE` e
  modelo de mudança, não rede estacionária.

Nos evals, a rede de uma pessoa é igual à chamada direta do graphicalVAR e recupera a PDC
simulada (r > ,90); a multinível é igual ao mlVAR direto e recupera os efeitos fixos
simulados (r > ,95). Leitura: aresta temporal é precedência no tempo dentro da pessoa,
não causa; aresta ausente numa rede regularizada não mostra ausência de efeito. A rede de
uma pessoa é resultado sobre um participante: com dado real, na conversa só a estrutura,
sem identificar a pessoa, e os dados e as redes por pessoa do mlVAR ficam no objeto.

## 11. Armadilhas

- Aresta transversal como seta causal; centralidade alta como alvo de intervenção comprovado
  (hipótese; erro de paradigma que persiste com tudo estável).
- Layout lido como centralidade; espessuras comparadas a olho; intermediação sem CS
; CS baixo como invalidação; IC regularizado como teste;
  centralidade comparada entre grupos sem teste.
- Rede sem γ, correlação, ausentes e N declarados: padrões do bootnet e do JASP dão redes
  diferentes no mesmo banco (testado em `evals_redes.R`).

## 12. Fontes e divergências

Livro do usuário, Vol. 7 (md7): paradigma e alcance (69, 79, 101), estimação e padrões do JASP
(206, 216, 312), centralidade (233-258), ponte (250), previsibilidade (314), confiabilidade
(361-401), relato (431-435, 520), formato (560-562). Notas do curso (análise de redes,
`spec_T6.md`): estimador por tipo de dado (P4, P9), centralidade (P5), precisão (P6), tamanho
de amostra (P7), estabilidade (P8), aulas no JASP (P9, P10); EGA em `notion_g4.md`, seção 11.
Literatura: Epskamp, Borsboom e Fried (2018); Epskamp e Fried (2018); Isvoranu e Epskamp
(2023); Bringmann et al. (2019); Hallquist et al. (2021); Jones et al. (2021); Haslbeck e
Waldorp (2018); Constantin et al. (2026); van Borkulo et al. (2023); Golino e Epskamp (2017);
Golino et al. (2020); Christensen e Golino (2021); Kursa e Rudnicki (2010); Popa et al. (2025);
Kruis e Maris (2016); van Bork et al. (2021).

Reproduzido nos evals (BFI Network, JASP 0.95.1): 181 arestas não nulas de 300 (53 negativas,
esparsidade .397), N1–N2 = .570, N3–N4 = .302, força em z de C4 (1.649) a C3 (-1.706) e a
tabela de centralidade inteira do JASP (25 itens x 4 índices), todos a .0005.

Divergências decididas:
- "Cerca de metade das 300 arestas": o JASP dá 181 de 300 (60%); a função dá o número
  exato, como o livro pede. O ,53 do pré-cálculo (; corpo em, com a
  nota de transposição de) e o IC [,49; ,56] dele não se combinam
  com o ,57 do JASP.
- Do curso, não adotados: CS ≥ .70 "excelente" (P8 a P10; a fonte só tem .25 e .50);
  nó central como alvo de intervenção (P1, P8); N mínimo de 300 a 500 e ganhos a cada
  50 casos (P7, sem fonte); estabilidade de item no bootEGA acima de 50% (`notion_g4`
  seção 11; vale ≥ .70, Christensen & Golino, 2021); justificativa do não paramétrico por viés
  do LASSO e "19.000 replicações" (P6, sem fonte; valem e); C3 como central
  (P10; C3 tem a menor força).
- γ = 0.5 também no IsingFit e no mgm, como no JASP (o bootnet usa 0.25 nos dois).
- O JASP faz o case-dropping com força, proximidade e intermediação; a skill usa força e
  influência esperada por padrão e as outras sob pedido.
- Bootstrap em Ising: o curso cita outros métodos para binários sem nomeá-los (P6); a função
  roda o não paramétrico e avisa.
