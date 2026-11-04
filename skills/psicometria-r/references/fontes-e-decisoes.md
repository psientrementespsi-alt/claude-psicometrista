# Fontes, decisões transversais e erros das fontes que a skill não reproduz

Índice: 1. Hierarquia das fontes · 2. Rotas: clássica e alternativa · 3. Descritivas,
normalidade e atípicos · 4. Convenções de relato · 5. Erros das fontes não reproduzidos ·
6. Onde a literatura discorda do livro · 7. Afirmações e o que elas exigem

## 1. Hierarquia das fontes

1. **Livro do usuário**: Biblioteca Análise de Dados (Análise de Dados na Ciência com
   IA), 7 volumes, exemplos no JASP 0.95.1 com bancos da biblioteca do JASP. É a fonte
   dos padrões das funções. As referências `mdX:linha` nos arquivos de módulo apontam
   para os volumes (md1 a md7).
2. **Notas do curso Psicometria Online Academy** (páginas do Notion compartilhadas com
   o usuário): valem onde o livro é omisso (não paramétricos, McNemar e Cochran, ROC,
   cálculo amostral no G*Power, GLM, GEE, multinível, meta-análise, construção de
   instrumentos). Material de terceiros: a skill parafraseia, não copia.
3. **Pasta Psicometria-R do usuário** (v1.4): seis livros-base completos (Petersen,
   2024; Martinková & Hladká, 2023; Isvoranu et al., 2022; Irwing et al., 2018, o Wiley
   Handbook; Cook & Pitoniak, 2026; Rust et al., 2021), 58 recortes por módulo da
   skill, o mapa dos módulos, a auditoria CoVe de fidelidade e o curso Psicometria
   Online Academy completo (transcrições, sintaxes R e Rmd, inputs de Mplus). Os
   livros decidem a psicometria onde o livro de análise de dados é omisso (número de
   fatores, bifator, poder para AFC, viés preditivo, análise de itens, TRI, EMA); as
   transcrições e o código do curso decidem os procedimentos que o curso ensina e os
   livros não (aquiescência, desejabilidade, escolha forçada, classes latentes,
   epidemiologia). Material de terceiros: a skill parafraseia, não copia.
4. **Literatura com DOI conferido** (`referencias-verificadas.md`): decide quando livro
   ou curso contradizem a evidência metodológica, e entra como alternativa quando o
   livro escolheu outra rota.

Conflito entre livro e curso: vale o livro. Livro omisso e curso errado diante da
literatura: vale a literatura, com a nota no módulo. Erro de conta, de transcrição ou
de atribuição em qualquer fonte não é reproduzido (seção 5).

## 2. Rotas: clássica e alternativa

O livro expõe a rota clássica e registra a alternativa como defensável (Caixa de teto
13.3). A skill faz o mesmo:

| Situação | Padrão (livro) | Alternativa, sempre calculada ao lado |
|---|---|---|
| t independente | Levene decide Student ou Welch | `rota = "welch"`: Welch sempre (Delacre, Lakens & Leys, 2017) |
| ANOVA de um fator | Levene decide ANOVA padrão ou Welch com Games-Howell | Welch e Brown-Forsythe ao lado; `rota = "welch"` |
| Medidas repetidas | Mauchly; com violação, Greenhouse-Geisser se ε ≤ ,75 e Huynh-Feldt acima | `esfericidade = "epsilon"`: decide só pelo ε, sem o Mauchly como porta |
| Pós-hoc | só com ômnibus significativo; Bonferroni (3 ou 4 comparações), Tukey (todos os pares), Games-Howell (variâncias desiguais) | Holm nos não paramétricos (Dunn, Conover), com a correção do curso ao lado |
| ANCOVA com pressupostos violados | bootstrap com 1.000 reamostras | o bootstrap não corrige heterogeneidade das inclinações (teste fator × covariável sempre) |
| Mediação | bootstrap percentílico, 10.000 reamostras, semente 2026 | 5.000 (curso e PROCESS) como argumento |
| Centralização no multinível (livro omisso) | média geral, como o curso | `tipo = "grupo"` (CWC), que separa efeito dentro e entre clusters (Enders & Tofighi, 2007) |
| Kruskal-Wallis e Friedman (livro omisso só nomeia os testes) | Dunn e Conover com Holm (literatura: Holm, 1979) | Bonferroni, que o curso usa, como argumento |

Nenhuma rota troca de teste em silêncio: a escolha e o motivo vão para `$decisoes`.

## 3. Descritivas, normalidade e atípicos

- Assimetria e curtose no estimador do JASP e do SPSS (G1 e G2, curtose em excesso), com
  erro-padrão; `descritivas()` também dá IC 95% da média, mediana, IIQ e Shapiro-Wilk.
- Normalidade é avaliação contínua, nunca liga-desliga: |assimetria|
  ou |curtose| abaixo de 1 é afastamento leve, acima de 2 preocupa (heurística do livro).
  O Shapiro-Wilk ganha poder com o n e acusa desvio irrelevante em amostra grande; com n
  pequeno, quase nunca acusa. Gráfico Q-Q pesa mais que o p. `ler_normalidade()` aplica
  essa leitura.
- O que precisa ser normal depende da análise: resíduos na regressão e na ANOVA, a
  diferença no t pareado, a variável dentro de cada grupo no t independente.
- Na dúvida entre Pearson e Spearman, rodar os dois e comparar; se divergirem, olhar a
  dispersão. Nada de trocar de teste em silêncio.
- Atípicos: `atipicos()` por |z| > 3,29 ou por MAD (Leys et al., 2013, limite 2,5).
  Caso atípico não sai do banco automaticamente; o critério é declarado antes.
- Mudança clínica individual em pré e pós: `mudanca_confiavel()` (Jacobson & Truax,
  1991), com DP e confiabilidade de referência declarados. O RCI é convenção de ensaio
  clínico e supõe o mesmo erro de medida em toda a escala; a alternativa psicométrica
  é o erro padrão condicional ao nível do traço pela TRI (`tri_ajustar()`, curva de
  informação), que muda ao longo da escala (auditoria CoVe da pasta, módulo 01).
- Dados ausentes (v1.5): `ausentes(dados, vars)` resume a porcentagem por variável, os
  casos completos, os padrões (só os com n ≥ 5 no print) e o teste MCAR de Little (1988),
  igual ao `naniar::mcar_test` nos evals. Little não significativo não prova MCAR;
  significativo afasta MCAR sem dizer se é MAR ou MNAR. Rota para a análise: FIML (ML ou
  MLR no lavaan e nos modelos mistos) ou imputação múltipla; a imputação única por EM
  das aulas, no SPSS, estima médias e covariâncias, mas subestima os EP em testes. Com
  mais de 20% de ausentes, análise de sensibilidade.
- Imputação múltipla (v1.7): `imputacao_multipla(dados, vars, formula, familia, m)`
  imputa por equações encadeadas (mice; PMM nas numéricas), ajusta o modelo linear ou o
  GLM em cada banco e combina pelas regras de Rubin com os gl de Barnard e Rubin (1999),
  iguais ao `mice::pool` nos evals, com a FMI por coeficiente, o R² combinado e o mesmo
  modelo em casos completos ao lado. Variáveis auxiliares entram em `vars` (tornam o MAR
  mais plausível; Enders & Baraldi, 2018); m pelo menos igual à % de casos incompletos
  (White et al., 2011), com aviso. Os bancos imputados ficam em `$imputacoes` para
  outras análises (`with()` e `mice::pool()`); têm linhas de pessoa e não se imprimem.
  A aula de dados faltantes faz a imputação múltipla no SPSS e diz que ela pouco difere
  da EM; a diferença está nos EP, que só a múltipla carrega.

## 4. Convenções de relato

Livro, Apêndice A e Vol. 2:

- símbolos estatísticos em itálico no manuscrito (*M*, *DP*, *t*, *F*, *p*, *r*);
- vírgula decimal em português; descritivas com duas casas;
- p exato com duas ou três casas, *p* < ,001 abaixo disso, nunca *p* = ,000;
- sem zero à esquerda no relato formal de quantidade que não passa de 1 em módulo
: p, r, ρ, τ, R², η², ω², V, φ, AUC, poder. O zero fica em d, g, OR, B,
  f² e β padronizado, que podem passar de 1. Os parágrafos de relato do
  Vol. 3 trazem "0,477" e "0,077" com zero, e a frase-modelo de cálculo amostral traz
  "poder de 0,80", contra a regra do próprio livro; a skill segue a regra;
- tamanho de efeito com IC; N efetivo e mecanismo de exclusão de casos declarados.

Casas: 2 para estatística de teste e para efeitos do tipo d, g e OR; 3 para p e para
proporções de variância (η², ω², R², ΔR², f², ε²). O módulo de associação usa 3 casas
em B, β, OR, r, ρ, τ, V, φ, w e AUC, como os parágrafos-modelo do livro (;
), porque com OR abaixo de ,10 duas casas perdem informação (0,072
viraria 0,07); processos também usa 3, como o livro.

Separador decimal: o padrão das funções é `decimal = "."`, neutro entre idiomas, para
que a mesma saída sirva a manuscrito em inglês. Manuscrito em português pede
`decimal = ","`, que segue o livro e troca a vírgula entre os itens da frase
por ponto e vírgula. A IA que monta o relato em português passa `decimal = ","`.

## 5. Erros das fontes não reproduzidos

Livro do usuário:

| Local | Problema | O que a skill faz |
|---|---|---|
| | 0,200 chamado de OR dos homens; é a chance (odds) dos homens | OR = 0,100 (chance 90% menor) |
| | "cerca de metade das 300 arestas" | o JASP dá 181 de 300 (60%); a função conta as arestas não nulas |
| | N1-N2 = ,53 (337 e corpo) com IC [,49; ,56] (542, 560), do pré-cálculo | o JASP dá ,57; não usar ,53 e o IC como gabarito do ,57 |
| | Henze-Zirkler = 14,459 | o número não se reproduz; a conclusão sim: nos resíduos da MANOVA, que são o objeto do pressuposto, HZ = 0,92, p = ,023 (nas VDs brutas, 0,69, p = ,108); a skill testa os resíduos |
| | g de Hedges "preferível quando os grupos têm tamanhos bem diferentes" | a correção J de Hedges depende só dos graus de liberdade (amostra pequena), não do desbalanceamento; a skill dá o g ao lado do d sem aviso por razão de tamanhos |
| Vol. 4 x Vol. 5 | η²p como padrão em desenho com mais de um fator (Vol. 4) contra ω²p (Vol. 5) | ω²p em fatorial e ANCOVA entre sujeitos; η²p em desenho com fator intra |
| (caixa de curtose) | fórmula da caixa não é o estimador do JASP (0,51 contra 0,53); o próprio livro registra a diferença em | estimador do JASP |
| (t pareado) | d de 0,47 sem fórmula declarada | é o d com o DP combinado das duas medidas, √((DP1² + DP2²)/2), o principal; o d_av de Lakens (média simples dos DP), o d_z e o d_rm saem ao lado |

Curso Psicometria Online (resumo; detalhes nos módulos):

- transcrição: "|z| > 0,96" (é 1,96); t de amostra única −0,69 (é −1,69); diferenças de
  média da ANOVA de um fator que não batem com as médias; AP5 CONT 75,25 (é 74,25); f
  "0,225" (é 0,2526); N de 189 na ANOVA de um fator com f = ,23 (o mínimo é 186);
- rótulo: ω²p = 0,037 na ANCOVA, que é η²p; "R² parcial 0,067", que é f²; g de Hedges da
  MANOVA fatorial (0,12 e 0,37; recalculados 0,22 e 0,49); faixas de d "0,20-0,39 e
  0,40-0,79" atribuídas a Cohen, que não são dele; faixas de r_rb atribuídas a Grissom e
  Kim, não achadas em fonte publicada;
- método: g de Hedges como correção de heterocedasticidade (não é); proporção mediada
  como c'/c (é ab/c, com ressalvas); Sobel descrito como liberal (é conservador); desvio
  = 2 logLik (é −2 logLik); QIC para comparar famílias de distribuição (diferença de
  escala, não de ajuste); Wald com dois graus de liberdade; Egger para OR (Peters ou
  Harbord); CS ≥ ,70 como "excelente" e nó central como alvo de intervenção (o livro
  chama isso de erro); N fixo de 300 a 500 para redes; estabilidade de item do bootEGA
  acima de 50% (a literatura usa ≥ ,70); "Number of groups = 4" em desenho só intra
  (é 1); OR de Chen para preditor binário aplicados a preditor contínuo;
- método (v1.4, conferido nas transcrições): teste de MCAR de Little lido com a H0
  invertida ("a hipótese nula é de que os dados não são completamente aleatórios"; a
  H0 do teste é MCAR, e não rejeitar não prova MCAR); NCT com p > ,05 lido como
  "os grupos são invariantes" (não significativo não é igual; ver o aviso de não
  significância das funções e o TOST); reteste pela Pearson, que não enxerga mudança
  sistemática entre as ocasiões (a skill usa o ICC(A,1) como principal e põe a
  Pearson ao lado).
- método e código (v1.5, conferido nos scripts das aulas): imputação única por EM do
  SPSS como solução para o ausente depois do Little não significativo (serve para
  estimar médias e covariâncias, mas subestima os EP em testes; a skill recomenda FIML ou
  imputação múltipla, `ausentes()`); AFE e AFC em metades aleatórias do mesmo banco como
  "confirmação" (as próprias aulas extras desaconselham; `construcao-instrumentos.md`);
  `nparLD.R` com argumentos que o nparLD 2.3 removeu (`description`, `plot.CI`) e
  `mctp` com `asy.method = "log.odds"`, cujo estimador sai em escala transformada
  (`naoparametricos-complexos.md`); `Latent growth curve com SEM.R` com a sintaxe rodada
  diferente da descrita nos comentários (`mee.md`, seção 7); `Franchini2012.R` só com o
  modelo de efeito comum, sem olhar a heterogeneidade antes (`modelos.md`, 7c);
  `pairwise()` chamado do netmeta, de onde saiu na versão 3 (está no pacote meta); rede
  bayesiana única, sem reamostragem, lida pela direção das setas
  (`redes-ega.md`, 10b); intervalos de confiança do Joinpoint no relato da aula sem o
  método declarado (`epidemiologia.md`).
- pacotes de referência (v1.6, conferidos ao montar os evals): o `dmacs` 0.1.0, que
  Petersen (2024) usa e que saiu do CRAN, combina os DPs dos grupos linearmente (não as
  variâncias), multiplica o integral por √(variância do traço focal), o que só se anula
  quando essa variância vale 1, e no caso ordinal calcula o escore esperado sem a escala
  do resíduo da resposta latente; a skill usa a forma fechada para indicadores contínuos
  (igual ao `pinsearch::dmacs`) e deixa itens ordinais para `tri_dif()`
  (`invariancia.md`). O `irr::kripp.alpha` dá valores levemente diferentes do irrCAC e da
  fórmula de Krippendorff com 3 ou mais juízes; a referência dos evals é o irrCAC
  (`construcao-instrumentos.md`, 3c). O "Intercept" do `rms::validate` é o intercepto da
  recalibração junto com a inclinação, não a calibração geral (offset) que a skill relata
  (`associacao.md`, seção 6).
- aula de limitações do kappa (v1.6): o curso mostra o paradoxo (90% de concordância e
  κ = 0) e encerra dizendo que "há outras medidas" sem nomear; a skill acrescenta o AC1 de
  Gwet e o alfa de Krippendorff e avisa o paradoxo (`concordancia_juizes()`).

## 6. Onde a literatura discorda do livro

Regras do livro que a literatura não sustenta, mantidas como padrão por serem a fonte do
usuário, com a alternativa calculada: escolha entre Student e Welch pelo Levene
(Delacre et al., 2017); Mauchly como porta da correção e o corte de ,75 do ε, regra de bolso que Barcikowski e
Robey (1984) e Girden (1992) formularam a partir dos resultados de Huynh e Feldt (1976),
como descrevem Blanca et al. (2023). Nos resíduos ajustados do
qui-quadrado o livro já manda endurecer o limiar de 1,96 na lógica de Bonferroni quando
as células são muitas, como Beasley e Schumacker (1995); a skill aplica.
Pós-hoc da ANOVA de um fator: o livro aceita Bonferroni com 3 ou 4 comparações e prefere
Tukey quando o número cresce; o automático da skill usa Tukey (ou Games-Howell
com variâncias desiguais) em qualquer número de grupos e registra a alternativa em
`$decisoes`.
Regras do curso que a skill não adota por contrariar a literatura: Durbin-Watson por
faixa fixa de 1,5 a 2,5 (a skill dá o teste com p),
centralização para reduzir multicolinearidade, Spearman automático com dado não normal
(Bishara & Hittner, 2012). Detalhes e referências em cada arquivo de módulo, seção
"Fontes e divergências".

## 7. Afirmações e o que elas exigem

Frases comuns em manuscrito que o resultado usual não sustenta sozinho. A coluna do meio diz
o que a afirmação pede; a da direita, onde a skill confere.

| Afirmação | O que ela exige | Onde a skill confere |
|---|---|---|
| "Os grupos não diferem", "não há associação" | p ≥ α não basta: os dados também não destoam de outros valores do IC (Greenland et al., 2016). Para afirmar ausência de efeito que importe, equivalência com limites fixados antes (Lakens, 2017) | aviso "Não significativo não quer dizer igual" em `$avisos` das funções de diferença e de associação; `teste_equivalencia()` |
| "A diferença é pequena demais para importar" | limites do menor efeito de interesse, justificados antes dos dados, e IC 90% dentro deles | `teste_equivalencia()`, com os quatro resultados de Lakens (2017) em `$resultado` |
| "O escore total mede o construto" | a soma supõe o modelo paralelo; a AFC congenérica não é evidência para ela (McNeish & Wolf, 2020) | `checar_escore_somado()` (`references/afc-ordinal.md`, Escore somado) |
| "Cada subescala pode ser somada" | paralelo em cada subescala; no modelo do escore, fatores sem covariância | `checar_escore_somado(dados, list(...))`, uma decisão por subescala |
| "Não rejeitamos o paralelo, logo a soma vale" | não rejeitar não mostra que as restrições valem; olhar poder, ajuste do paralelo e amplitude das cargas | `$avisos` e `$amplitude` de `checar_escore_somado()` |
| "O construto é uma rede, não um fator" (ou o contrário) | razão teórica sobre o mecanismo; modelos equivalentes têm o mesmo ajuste (Kruis & Maris, 2016; van Bork et al., 2021). Fator único contra rede esparsa tem teste próprio (PCL) | `references/redes-ega.md`, seção 10; o PCL não está implementado |
