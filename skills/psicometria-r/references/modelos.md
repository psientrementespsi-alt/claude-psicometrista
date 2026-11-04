# Modelos avançados: GLM, ordinal e multinomial, GEE, multinível e meta-análise

Funções em `scripts/est_modelos.R`; testes em `scripts/evals_modelos.R` (136 sem
`PR_DADOS_JASP` na v1.5, com meta-análise em rede e extração de dados; cerca de 25 s). O livro do usuário não cobre estes modelos (o GEE aparece como
"horizonte, fora do escopo" e). As regras vêm das notas do curso,
conferidas contra a literatura; os erros do curso ficam na seção 9 e não são
reproduzidos.

Índice: 1. Qual análise | 2. Funções | 3. GLM | 4. GEE | 5. Multinível linear |
6. Multinível logístico | 7. Meta-análise (7c em rede; 7d extração de dados) | 8. Âncoras |
9. Fontes e divergências

## 1. Qual análise

| Desenho e pergunta | Função |
|---|---|
| Uma medida por pessoa; Y de contagem, 0/1 ou positiva assimétrica; n pequeno | `modelo_glm()` |
| Duas ou mais medidas por pessoa, efeito médio populacional | `modelo_gee()` |
| Pessoas dentro de grupos, efeitos de nível 1 e 2, variação entre grupos | `multinivel_linear()` |
| Idem com desfecho 0/1 | `multinivel_logistico()` |
| Combinar estudos | `calcular_efeitos_meta()` e `meta_analise()` |
| Efeito muda com característica dos estudos | `meta_subgrupos()`, `meta_regressao()` |
| Combinar prevalências (v1.4) | `meta_proporcao()` |
| Combinar sensibilidade e especificidade de um teste (v1.4) | `meta_acuracia()` |
| Razão de prevalência, sobrevida (v1.4) | `razao_prevalencia()`, `sobrevida()` (`references/epidemiologia.md`) |
| Subgrupos latentes de pessoas (v1.4) | `classes_latentes()`, `perfis_latentes()` (`references/classes-perfis.md`) |
| Desfecho com 3 ou mais categorias ordenadas (v1.6) | `modelo_ordinal()` |
| Desfecho com 3 ou mais categorias sem ordem (v1.6) | `modelo_multinomial()` |
| Efeito da média do grupo além da posição individual (v1.6) | `efeito_contextual()` |

Medidas repetidas: GEE dá efeito populacional e exige MCAR; o modelo misto dá efeito
condicional ao sujeito e vale sob MAR (o curso não tem regra; decide a pergunta).

## 2. Funções

```r
escolher_familia(dados, vd)
modelo_glm(dados, vd, preditores, familia = "normal", ligacao = NULL, referencia = NULL,
           fatores = NULL, escala = c("ml", "pearson"), exposicao = NULL, comparar = NULL)
comparar_familias(dados, vd, preditores, familias = c("normal", "gama"), ligacoes = NULL, fatores = NULL)
modelo_gee(dados, vd, preditores, id, tempo, familia = "normal", ligacao = NULL, correlacao = "unstructured",
           ep = c("robusto", "mancl_derouen"), referencia = NULL, fatores = NULL, comparar = NULL)
comparar_correlacoes_gee(dados, vd, preditores, id, tempo, familia, ligacao, correlacoes, fatores = NULL)
dados <- centralizar_multinivel(dados, vars, cluster, tipo = c("geral", "grupo"))
multinivel_linear(dados, vd, cluster, nivel1, nivel2, inclinacao, interacao, referencias,
                  fatores = NULL, modelo_final = NULL, estimador_final = c("REML", "ML"))
multinivel_logistico(dados, vd, cluster, nivel1, nivel2, inclinacao, interacao, referencias,
                     fatores = NULL, evento = NULL, modelo_final = NULL, nAGQ = 1)
calcular_efeitos_meta(dados, medida = c("SMD", "MD", "OR", "RR", "COR"), estudo,
                      media1, dp1, n1, media2, dp2, n2, eventos1, eventos2, r, n)
meta_analise(dados, yi = "yi", vi = "vi", estudo, medida, metodo = "REML",
             hartung_knapp = FALSE, sensibilidade = c("PM", "DL"), vies = TRUE, k_min_vies = 10)
meta_subgrupos(dados, subgrupo, yi, vi, estudo, medida, metodo = "REML", hartung_knapp = FALSE)
meta_regressao(dados, moderadores, yi, vi, medida, metodo = "REML", hartung_knapp = FALSE,
               fatores = NULL, permutacao = FALSE, semente = 2026, reamostras = 1000)
modelo_ordinal(dados, vd, preditores, ligacao = c("logit", "probit"), referencia = NULL, fatores = NULL)   # v1.6
modelo_multinomial(dados, vd, preditores, referencia_vd = NULL, referencia = NULL, fatores = NULL)       # v1.6
efeito_contextual(resultado, var)                                                                       # v1.6
```

Todas aceitam `nivel` (o rótulo do IC na frase acompanha: `nivel = .90` escreve "IC
90%"), `decimal` e `idioma`; `preditores` aceita `"grupo * tempo"`.
Devolvem `pr_resultado` (`modelo_glm` e `modelo_gee` com a subclasse `pr_modelo`, que só
muda o print), exceto `centralizar_multinivel` (banco com colunas novas,
invisível, com print só de agregados) e `calcular_efeitos_meta` (estudos com `yi`, `vi`,
`sei`, `n_total`).

Grupo gravado como número (1/2/3, comum em bancos do SPSS e do JASP) entra como
covariável contínua, com 1 gl, e pode inverter a decisão. No Viagra, com o Wald χ² da
função (escala ML) nas duas codificações, dose numérica dá χ²(1) = 3,89, p = .049, e dose
como fator dá χ²(2) = 5,37, p = .068; com o F clássico nas duas, F(1, 28) = 3,63,
p = .067, e F(2, 27) = 2,416, p = .108, nenhuma significativa (a inversão aparece no
Wald, que fica liberal com n pequeno). O livro avisa o mesmo no JASP. Toda
função de modelo aceita `fatores = "grupo"` e avisa quando um preditor numérico fora de
`fatores` tem de 3 a 7 valores distintos, todos inteiros (a mesma regra da associação e
dos processos; 0/1 não avisa, porque é a codificação que o livro manda para a
dicotômica).

## 3. GLM

Distribuição e ligação (curso, "Funções de Ligação"): contínua simétrica,
normal-identidade (B = diferença média); positiva assimétrica, gama-identidade
(diferença na escala original, opção do curso) ou gama-log (exp(B) = razão de médias);
0/1, binomial-logit (OR); contagem, Poisson-log (IRR); contagem superdispersa,
binomial negativa-log (`MASS::glm.nb`). `escolher_familia()` sugere pelo tipo de Y. Sem
`ligacao`, a gama usa identidade com só fatores (mesmas médias ajustadas da log) e log
com preditor contínuo (a identidade pode gerar média negativa). Zeros impedem a gama: a
função recusa e explica que somar constante muda o modelo.

- Escala da normal e da gama por máxima verossimilhança (`escala = "ml"`), como o SPSS
  GENLIN do curso: reproduz os EPs das médias e o Wald do curso. Com n pequeno o Wald
  fica liberal; `escala = "pearson"` usa n - p (padrão do R). Na normal, o F clássico
  (= ANOVA Tipo III) sai na mesma tabela como sensibilidade ("Respondendo Revisores").
- Termos: Wald χ² Tipo III com contrastes soma-zero, um gl por coluna do termo.
- Coeficientes contra a referência (primeiro nível; o SPSS usa o último:
  `referencia = list(grupo = "controle")`), IC de Wald, exp(B) com IC em log e logit.
- Médias marginais (`emmeans`, escala da resposta) e pares com Bonferroni, IC dos pares
  também ajustado; com log os pares são razões, com logit são OR. `$posthoc` traz a
  coluna `n` (casos do nível ou da célula; no GEE, sujeitos): nível com n < 5 sai do
  print sem média e IC (com n = 1 seriam o escore da pessoa), o valor fica no objeto e
  `$avisos` lembra que o par e o coeficiente desse nível, somados à média do outro
  nível, devolvem o escore. A frase só traz os pares, nunca a média de um nível.
- Poisson: dispersão de Pearson (χ²/gl, p); com p < .05 avisa e recomenda binomial
  negativa, com o AIC das duas ao lado. `exposicao` entra como `offset(log(...))`.
- AIC = -2 logLik + 2k com k contando a escala (na gama, forma por ML como no SPSS).
  `comparar_familias()`: AIC, ΔAIC e peso de Akaike, só para o mesmo Y, casos e
  preditores; avisa se mistura distribuição contínua com discreta.

Efeito principal: B (identidade) ou exp(B) (log/logit) de cada coeficiente, com IC.

```r
modelo_glm(ratos, "congelamento", "grupo", referencia = list(grupo = "AP5 CONT"))$apa[1]
# Generalized Linear Model (distribuição normal, ligação identidade); AIC = 286.27:
# grupo, Wald χ²(2) = 32.28, p < .001.
```

Fora: proporção agregada. Normal-identidade é a própria ANOVA. Ordinal e multinomial
têm função própria desde a v1.6 (seção 3b).

## 3b. Desfecho ordinal e nominal (v1.6)

Martinková e Hladká (2023, 6.7) apresentam os dois modelos como modelos de item; aqui
servem a qualquer desfecho com 3 ou mais categorias. `escolher_familia()` aponta para eles
quando o desfecho é categórico com mais de 2 níveis.

`modelo_ordinal()` ajusta o logito cumulativo de chances proporcionais (`MASS::polr`; probit
com `ligacao = "probit"`). Parametrização do PLUM do SPSS: logit P(Y ≤ j) = limiar_j − Xβ,
de modo que B > 0 desloca para categorias mais altas e OR = exp(B) é a razão das chances
de estar acima de qualquer corte. Sai: coeficientes com IC de Wald, limiares, Wald χ²
Tipo III por termo (contrastes soma-zero), teste global por razão de verossimilhança,
pseudo-R² de Cox e Snell, Nagelkerke e McFadden, e o teste das linhas paralelas como no
PLUM: razão de verossimilhança contra o modelo com um coeficiente por corte
(`ordinal::clm` com `nominal =`), global e por variável (igual ao `ordinal::nominal_test`).
Desfecho numérico com até 10 valores vira fator ordenado pela ordem dos valores; fator sem
ordem declarada usa a ordem dos níveis, com aviso. Linhas paralelas rejeitadas com n
grande pedem olhar os OR por corte (logísticas binárias em cada corte) antes de trocar de
modelo; as saídas são a multinomial ou o modelo de chances proporcionais parciais
(`ordinal::clm(..., nominal = ~ x)`). Âncora: `MASS::housing` (Venables e Ripley), com
InflHigh = 1,2888, TypeTerrace = −1,0910 e limiar Low|Medium = −0,4961, conferidos nos evals.

`modelo_multinomial()` ajusta o logito de categoria-base (`nnet::multinom`): um conjunto
de coeficientes por categoria contra a referência (`referencia_vd`, padrão o primeiro
nível), RRR = exp(B) com IC de Wald, testes dos termos por razão de verossimilhança
(modelo sem o termo, como os "Likelihood Ratio Tests" do NOMREG), pseudo-R² e acerto de
classificação contra a categoria mais comum. Célula vazia no cruzamento com um fator e EP
acima de 10 avisam separação. Os coeficientes batem com o `VGAM::vglm` nos evals.

```r
modelo_ordinal(dados, "gravidade", c("grupo", "idade"), fatores = "grupo")$apa
modelo_multinomial(dados, "diagnostico", c("idade", "sexo"), referencia_vd = "controle")$apa
```

## 4. GEE

Formato longo, uma linha por sujeito e tempo; a função ordena por `id` e `tempo` e cada
linha ausente sai sozinha. `tempo` dentro de `preditores` entra como fator; passe-o como
fator na ordem cronológica (texto é ordenado alfabeticamente, com aviso).

- Correlação de trabalho: padrão do curso é a não estruturada; se não convergir ou der
  matriz não positiva definida, tenta AR(1) e depois intercambiável, com aviso. Com 2
  tempos as três são a mesma.
- AR(1) no geepack usa a ordem dos tempos, não a distância: com intervalos desiguais a
  função avisa (resolve a contradição do curso entre "intervalos irregulares" e "AR(1)
  supõe espaçamento regular").
- EP sanduíche por padrão: consistente mesmo com a correlação errada (Liang e Zeger,
  1986). Com menos de 40 sujeitos ele infla o erro tipo I (Mancl e DeRouen, 2001): a
  função calcula a correção de Mancl e DeRouen (ausente no geepack) e a põe ao lado;
  `ep = "mancl_derouen"` a torna principal. Nos evals, com 12 sujeitos e α = .05, o
  robusto rejeitou 9,5% e o corrigido 5,0%.
- Ausência: GEE padrão só vale sob MCAR; sob MAR, GEE ponderado (Robins et al., 1995) ou
  imputação múltipla (Schafer e Graham, 2002). Avisa sempre que sai linha.
- `comparar_correlacoes_gee()`: QIC (Pan, 2001) na mesma família e ligação, CIC (Hin e
  Wang, 2009) ao lado; QICu para preditores. QIC entre famílias não é usado.
- Binomial negativa não existe no geepack (só com θ fixo): a função recusa.

Problemas do geepack 1.3.13 contornados: "unstructured" com tempos faltando no meio
derruba a sessão do R (segfault) e, sem `waves`, pareia correlações pela posição da
linha; a função usa `corstr = "userdefined"` com `zcor = genZcor(..., corstrv = 4)`. O
`geepack::QIC` não divide a quasi-verossimilhança pela escala (na normal vira soma de
quadrados em unidades de Y²) e escreve a da gama como -y/(μ - log μ); a função usa
-y/μ - log μ e a escala do modelo de independência, e bate com o geepack quando a
escala é 1.

Interação sempre relatada; pós-hoc de interação não significativa não é confirmatório.

## 5. Multinível linear

Sequência do curso (0-1-2-3-4, prática no JASP, reporte): M0 nulo; M1 nível 1; M2 nível
2; M3 inclinação aleatória de `inclinacao`; M4 interação `inclinacao × interacao`,
mantendo a inclinação aleatória. A página "Ajuste de Modelos Multiníveis" inverte M2 e
M3; segue-se a ordem das outras três.

- Entrada: `nivel2` precisa ser constante no cluster (erro); `nivel1` sem variação no
  cluster e contínuo não centralizado geram aviso.
- Comparações em ML (partes fixas mudam): parâmetros, logLik, desvio = -2 logLik, AIC,
  BIC, Δ-2LL, Δgl, p do LRT, R² marginal e condicional por modelo. M3 testa variância na
  fronteira: mistura 50:50 de χ²(gl - 1) e χ²(gl) ao lado (Self e Liang, 1987; Stram e
  Lee, 1994).
- Modelo final: último passo com melhora significativa (no M3, pela mistura) ou
  `modelo_final`; reestimado em REML (Pinheiro e Bates, 2000) com gl de Satterthwaite
  (Luke, 2017). O curso não fala de ML/REML; `estimador_final = "ML"` reproduz o curso.
- ICC do nulo em REML (principal) e ML, sem limiar de entrada (nota se < .05). Efeito
  de desenho = 1 + (m - 1) ICC, m = tamanho médio; "deff < 2 dispensa multinível" não
  vale com preditor de nível 2 ou clusters < 10 (Lai e Kwok, 2015).
- R² de Nakagawa e Schielzeth (2013), inclinações pela extensão de Nakagawa et al.
  (2017); confere com `performance::r2_nakagawa`.
- Inclinação aleatória: τ11, ρ01, faixa de 95% (γ10 ± 1,96 DP) e de 68% (± 1 DP, do
  curso).
- Diagnóstico: singularidade, avisos do otimizador, n de clusters (EPs de nível 2
  viesados com 50 grupos ou menos: Maas e Hox, 2005), clusters < 10 casos, assimetria e
  curtose dos resíduos.

- Inclinação aleatória só de preditor numérico; fator para na entrada (usar dummy 0/1).

`centralizar_multinivel()`: padrão `tipo = "geral"`, o do curso: `X_cgm` (nível 2
centrado na média entre clusters). Para nível 1 cria ao lado, sem custo, a alternativa
de Enders e Tofighi (2007): `X_cwc` (efeito intragrupo puro) e `X_media_cgm` (média do
cluster centrada, efeito entre grupos), recomendada para efeito de nível 1 puro e para
interação entre níveis. `tipo = "grupo"` cria só estas e recusa variável de nível 2.
Devolve o banco de forma invisível; o print mostra só n, clusters e resumo das colunas.

Efeito contextual (v1.6; a aula de regressão multinível estima o "efeito contextual" com a
média da escola como preditor de nível 2): com `X_cwc` e `X_media_cgm` no mesmo modelo,
`efeito_contextual(resultado, "X")` devolve b_media − b_cwc com EP pela covariância dos
efeitos fixos. É o quanto a média do grupo acrescenta além da posição da pessoa dentro
dele; com X centrado na média geral e a média do cluster no modelo, o coeficiente da média
já é esse número (os evals conferem a igualdade).

```r
multinivel_linear(Exam, "normexam", "school", nivel1 = c("standLRT", "sex"), nivel2 = c("schavg", "schgend"),
  inclinacao = "standLRT", interacao = "schgend", referencias = list(sex = "F", schgend = "mixed"))$apa[1]
# Modelo nulo: ICC = .168 (REML; .166 em ML), ou seja, 16.8% da variância de normexam
# está entre os clusters (school); efeito de desenho = 11.34.
```

Estrutura aleatória máxima (Barr et al., 2013) ou parcimoniosa (Matuschek et al.,
2017) é debate aberto: declarar a escolha.

## 6. Multinível logístico

Mesma sequência com `lme4::glmer` (binomial, logit), tudo por ML (Laplace; `nAGQ` > 1
só com intercepto aleatório, e então o LRT do M3 usa o modelo anterior reajustado por
Laplace, para não misturar aproximações; GLMM não tem REML). OR = exp(b) com IC de Wald; ICC
latente = τ00/(τ00 + π²/3) (Goldstein et al., 2002); R² com variância residual π²/3.
O texto diz que os OR são condicionais ao cluster (maiores em módulo que os
populacionais da logística simples ou do GEE quando τ00 é grande) e que τ00 muda de
escala ao entrar preditores (não comparar τ00 entre modelos como na linear).

## 7. Meta-análise

`calcular_efeitos_meta()` (metafor::escalc): SMD = g de Hedges (correção exata); MD; OR
e RR em log (0,5 nas tabelas com zero; estudos sem evento nos dois grupos saem); COR em
z de Fisher, com volta para r.

`meta_analise()`:

- Efeitos aleatórios por padrão, escolhidos pela hipótese sobre os efeitos verdadeiros,
  não pelo I² (curso, "Fixo vs Aleatórios"); τ² por REML (Veroniki et al., 2016; Langan
  et al., 2019). Ao lado: Hartung-Knapp (IntHout et al., 2014; principal com
  `hartung_knapp = TRUE`), Paule-Mandel, DerSimonian-Laird e efeito comum.
- Heterogeneidade: Q (p < .10 como no curso, pouco poder); I² principal = (Q - gl)/Q
  truncado em 0, a fórmula do curso e do pacote meta; I² do metafor (por τ²) com IC ao
  lado; τ² com IC Q-profile; faixas sobrepostas do Cochrane Handbook só como leitura (as faixas vêm do
  Handbook, não de Higgins et al., 2003, que propõe o I²).
- Intervalo de predição = μ ± t(k - 2) · √(τ² + EP²) (Riley et al., 2011; IntHout et
  al., 2016), com k >= 3; o padrão do metafor (normal) é mais estreito (BCG: .15 a 1,55
  com z; .13 a 1,78 com t(11)). Com k < 5, aviso de instabilidade. `metodo = "EE"` vira
  "efeito comum", sem τ² nem intervalo de predição. Leave-one-out com faixa e estudos que mudam a significância;
  `extra$estudos` traz IC e peso (%) de cada estudo para o forest plot.
- Assimetria do funil só com k >= 10 (Sterne et al., 2011): Egger clássico (Egger et
  al., 1997) para MD e COR; para OR e RR o principal é o de Peters et al. (2006),
  ln(efeito) em 1/N com pesos S·F/N (= `meta::metabias(method.bias = "Peters")`), com
  Harbord pelo escore (Harbord et al., 2006) para OR; ambos pedem as contagens de
  `calcular_efeitos_meta()`. Para SMD, regressão em 1/√n. Assimetria indica efeitos de estudos pequenos,
  não prova viés de publicação. Trim-and-fill só como sensibilidade (Duval e Tweedie,
  2000; falha com heterogeneidade: Terrin et al., 2003; Peters et al., 2007).

`meta_subgrupos()`: efeitos aleatórios com τ² próprio por subgrupo (padrão do meta e do
curso) e Q entre subgrupos; τ² comum (meta-regressão) ao lado. A modificação de efeito
se lê no teste entre subgrupos; avisa subgrupos com menos de 5 estudos e, com
Hartung-Knapp, que o Q entre fica mais liberal (EP de HK sem ajuste ad hoc pode encolher).

`meta_regressao()`: QM, QE, τ² e I² residuais e R² como no metafor (origem dos números
do curso); exp(B) e % de mudança para OR e RR; Hartung-Knapp e permutação (semente
2026) opcionais; avisa menos de 10 estudos por covariável. B com 4 casas quando |B| <
0,1. O R² é descritivo: cai ao acrescentar moderador com k pequeno (75,62% para 53,84%
no curso); comparar modelos por LRT em ML ou AIC.

```r
meta_analise(calcular_efeitos_meta(bcg, "RR", estudo = "estudo", eventos1 = "tpos", n1 = "n1",
             eventos2 = "cpos", n2 = "n2"), decimal = ",")$apa[1]
# Meta-análise de efeitos aleatórios (REML; k = 13): RR = 0,49; IC 95% [0,34; 0,70];
# z = -3,97; p < ,001; intervalo de predição 95% [0,13; 1,78].
```

### 7b. Meta-análise de prevalência e de acurácia (v1.4)

`meta_proporcao(dados, eventos, n, estudo, transformacao = "PLO")`, pelo módulo de
meta-análise de prevalência do curso: IC de cada estudo pelo método exato de
Clopper-Pearson, efeitos aleatórios (REML) sobre a proporção transformada e volta
para proporção, com intervalo de predição, Q, I² e τ². `transformacao`: `"PLO"`
(logito, padrão), `"PFT"` (duplo arco-seno de Freeman-Tukey, volta pela média
harmônica dos n) ou `"PR"` (bruta). O logito é o padrão porque a volta do
Freeman-Tukey pode distorcer a estimativa combinada (Schwarzer et al., 2019); com
muitos estudos de prevalência extrema, o GLMM do `meta::metaprop(method = "GLMM")` é a
alternativa. Igual ao `meta::metaprop` (inverso da variância, REML) nos evals, nas
três transformações.

`meta_acuracia(dados, vp, fp, fn, vn)`: modelo bivariado de Reitsma et al. (2005)
sobre os logitos de sensibilidade e especificidade (REML, `metafor::rma.mv`, efeitos
aleatórios correlacionados), com sensibilidade e especificidade combinadas e IC, RV+,
RV− e razão de chances diagnóstica (IC pelo método delta), τ² e correlação, e o teste
de assimetria de Deeks et al. (2005) (p < ,10). Correção de 0,5 em todas as células
quando algum estudo tem célula zero, como o `mada::reitsma`, com o qual os evals
conferem o banco AuditC (sensibilidade ,891 [,808; ,941]; especificidade ,780). O curso
tem só a aula de abertura de meta-análise de acurácia; o modelo segue o padrão da
literatura citada.

## 7c. Meta-análise em rede (v1.5)

`meta_rede(dados, estudo, tratamento, n, media, dp | eventos, referencia)`, a rota do
curso "Metanálise de Rede" (script `Franchini2012.R`, netmeta). Entrada por braço, no
formato longo (uma linha por braço de estudo); o `pairwise()` (que saiu do netmeta 3
para o pacote meta) monta os contrastes, inclusive de estudos com 3 braços.
- Modelo comum e aleatório lado a lado; o relato usa o aleatório por padrão (as aulas
  usam só o comum). Heterogeneidade: τ² (DerSimonian-Laird, padrão do netmeta) e I².
- Inconsistência: decomposição do Q (dentro dos desenhos = heterogeneidade; entre
  desenhos = inconsistência) e separação de nós (`netsplit`: direto contra indireto nas
  comparações que têm os dois). Q entre desenhos com p < ,10 vira aviso: transitividade
  em dúvida, ler o netsplit antes do ranking.
- Ranking por P-score (Rücker & Schwarzer, 2015), descritivo; `menor_melhor` diz o
  sentido. Não substitui os IC.
- Gráficos da rede, netheat, league table e funil comparativo das aulas: rodar à parte
  sobre `$modelo` (`netgraph`, `netheat`, `netleague`, `funnel`).
- Evals: com o banco Franchini 2012 da aula (no netmeta), efeitos iguais ao netmeta no
  formato largo; rede de dois tratamentos igual à meta-análise par a par do metafor
  (comum = FE; aleatório = DL); P-scores com média 0,5.

## 7d. Extração de dados (v1.5)

`extrair_media_dp(tipo, ...)`, a aula de extração de dados (Cochrane Handbook, cap. 6),
estudo a estudo:
- `"mediana_iqr"`, `"mediana_amplitude"`, `"mediana_iqr_amplitude"`: média de Luo et al.
  (2018) (padrão) ou de Wan et al. (2014), DP de Wan et al. (2014), pelas rotinas do pacote
  meta (tabela de Wan para n até 201). O exemplo da aula (n = 15; 9,6; 14,7; 30,7) dá média
  18,66 e DP 17,03.
- `"ep"` e `"ic"`: DP pelo EP ou pelo IC da média (t com n < 60, z acima, como o Handbook).
- `"combinar"`: dois grupos num só (fórmula do Handbook); igual ao banco unido nos evals.
- `"mudanca"`: DP da mudança pré-pós com a correlação r imputada; `"correlacao"`: r
  implícita num estudo que relata os três DP. Declarar a origem de r e testar outros
  valores na sensibilidade.
Converter mediana em média supõe distribuição perto da normal; com assimetria forte (o
motivo de o estudo ter relatado mediana), rodar a meta-análise com e sem os estudos
convertidos.

## 8. Âncoras

| Fonte | Números | Situação |
|---|---|---|
| Livro, Viagra, com `fatores = "dose"` | F(2, 27) = 2,416, p = ,108; ANCOVA dose F(2, 26) = 4,142, covariável F(1, 26) = 4,959 | reproduzidos (pulados sem `PR_DADOS_JASP`) |
| "GLM - Descrição dos Resultados" | Wald χ²(2) = 32,27; AIC 286,27; β 112,65 e 133,25; EP 18,04 e 16,14 | reproduzidos com dados montados das médias, n e σ²_ML (dos EPs), AP5 CONT = 74,25 |
| Multinível, JASP e reporte (`mlmRev::Exam`) | ICC .166 (ML) e .168 (REML), τ00 .172, σ² .848; M1 .559, -.171, Δ 1680,6; M2 .348, .194, .151, Δ 14,7; M3 Var .015, DP .122, r .48, faixa .427 a .672; M4 -.016, -.021, Δ 0,2, p .89; JASP .563 e .343 | reproduzidos (.5595 e .5498 truncados para .559 e .549 no curso) |
| Idem | Δ-2LL M3 vs M2 = 40,199 | não na última casa: R dá 40,181 (três otimizadores, contr.sum); tolerância .03 |
| Meta, DM de mudança (Cho, Fritz, Morone) | REML 1,7713 [1,0013; 2,5414], τ² .2407, Q 4,17, I² 52,0%; comum 1,68 [1,16; 2,20]; 8 estimadores de τ²; pesos | reproduzidos dos N, médias e DPs do forest plot |
| Meta, BCG (`metadat::dat.bcg`) | RR .49 [.34; .70], τ² .3132, I² 92,1%; subgrupos, Q entre 1,86 (p .3943); latitude B -.0291 (EP .0072), τ² .0764, I² 68,39%, R² 75,62%, QE 30,7331, QM 16,3571; alocação; latitude + alocação | reproduzidos |

Sem dados disponíveis (sem âncora): GEE do Experimento 2; bancos Richardson2018,
Balduzzi2019 e GomesNeto2017 do curso (nomes de banco, não referências).

## 9. Fontes e divergências

Livro: GEE fora do escopo; multinível e meta-análise ausentes. Vale
o curso onde não contraria a literatura.

Curso, GLM e GEE (erros não reproduzidos): AP5 CONT 75,25 é 74,25 (os β excedem a
diferença em 1,00; com 75,25 o Wald seria 31,76) e "GLM - Prática" (Salina ≈ 86,
β ≈ 102) é incoerente com o documento; "Wald (2, 23)" e "Wald (2, 24)" viram Wald
χ²(gl); QIC entre famílias (31,53 vs 16308,22) é artefato de escala, o QIC escolhe
correlação dentro da família (Pan, 2001); "o link não afeta o ajuste" é falso fora de
modelos saturados só com fatores; somar constante para a gama muda o modelo; binomial
negativa sem teste de superdispersão; sem EP sanduíche nem correção para 26 a 27 ratos
(função traz Mancl e DeRouen); sem discussão de MCAR (função avisa); contradição AR(1)
vs intervalos irregulares (seção 4).

Curso, multinível: desvio = 2 × logLik no checklist de "Ajuste" (correto: -2 logLik);
ordem M2/M3 invertida em "Ajuste"; ML vs REML não discutido; LRT da variância sem
correção de fronteira; "Var(u1j) p < 0,01" no reporte (variância não tem p no lme4; o p
vem do LRT); média geral como padrão (mantida: não é erro, responde outra pergunta;
o CWC de Enders e Tofighi, 2007, sai ao lado); ICC .168 vs
.166 é REML vs ML, não erro; a faixa de 68% (.427 a .672) confere com os valores exatos
(.5498 ± .1224), a suspeita de erro na extração não procede; "não é obrigatório manter a
inclinação aleatória no M4" (JASP) contraria o reporte, que a mantém; logística sem
nota de OR condicional e de reescala de τ00.

Curso, meta-análise: I² do curso = (Q - gl)/Q do pacote meta, o do metafor usa τ²
(92,1% vs 92,22% no BCG; a função mostra os dois); Egger em desfecho binário tem erro
tipo I inflado (função usa Peters e Harbord); Hartung-Knapp, intervalo de predição
e IC de τ² ausentes nas aulas (calculados); R² caindo na meta-regressão é artefato com k
pequeno; subgrupos com efeito comum na prática (função usa efeitos aleatórios).

Literatura (DOI conferido): Liang e Zeger (1986); Pan (2001); Hin e Wang (2009); Mancl
e DeRouen (2001); Robins et al. (1995); Schafer e Graham (2002); Maas e Hox (2005);
Enders e Tofighi (2007); Luke (2017); Self e Liang (1987); Stram e Lee (1994); Pinheiro
e Bates (2000); Lai e Kwok (2015); Goldstein et al. (2002); Nakagawa e Schielzeth
(2013); Nakagawa et al. (2017); Barr et al. (2013); Matuschek et al. (2017); Veroniki
et al. (2016); Langan et al. (2019); IntHout et al. (2014, 2016); Riley et al. (2011);
Higgins et al. (2003); Egger et al. (1997); Sterne et al. (2011); Peters et al. (2006,
2007); Harbord et al. (2006); Duval e Tweedie (2000); Terrin et al. (2003).
