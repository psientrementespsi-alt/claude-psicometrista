# Construção e validação de instrumentos: validade de conteúdo e AFE

## Sumário

- 12 (linha 24)
- 51 (linha 63)
- 98 (linha 110)
- 212 (linha 224)
- 225 (linha 237)
- 320 (linha 332)
- 356 (linha 368)
- 393 (linha 405)
- 417 (linha 429)

Funções em `scripts/est_construcao.R`; testes em `scripts/evals_construcao.R` (80
testes na v1.4, 82 com `PR_DADOS_JASP` apontando para os dados do JASP; cerca de 60 s).
Complementa a AFC de `scripts/psicometria.R` (auditoria, estimador, ajuste,
confiabilidade, validade discriminante), que não muda.

Índice: 1. Onde cada função entra · 2. Validade de conteúdo · 3. Análise fatorial
exploratória · 3b. Análise clássica de itens · 3c. Fidedignidade (reteste, juízes,
EMA) · 4. Depois da AFE · 5. Como relatar · 6. Armadilhas · 7. Fontes e divergências

## 1. Onde cada função entra

Etapas da construção no curso (página "Construção de Medidas de Autorrelato"): revisão
da literatura, definição constitutiva e operacional, itens, comitê de juízes, população-
alvo, e só então as evidências de estrutura interna, precisão e relações externas
(página "Confiabilidade Teste-Reteste": conteúdo, estrutura interna, precisão, relações
externas).

| Etapa | Pergunta | Função |
|---|---|---|
| Juízes: clareza, pertinência, relevância (1 a 5) | O item comunica e representa o construto? | `validade_conteudo_cvc()` |
| Juízes: relevância (1 a 4) | Quantos juízes acham o item relevante? | `validade_conteudo_ivc()` |
| Juízes: essencialidade (1 a 3) | O item é essencial, além de relevante? | `validade_conteudo_rvc()`, `rvc_critico()` |
| Juízes: qualquer escala ordenada | Alternativa com intervalo de confiança | `validade_conteudo_aiken()` |
| População-alvo (piloto) | Clareza, adequação e compreensão | `validade_conteudo_cvc()` com N participantes |
| Estrutura interna exploratória | Quantos fatores e quais itens em cada um? | `auditar_itens()`, `afe_adequacao()`, `afe_numero_fatores()`, `afe()` |
| Estrutura confirmatória (outra amostra) | O modelo se sustenta? | `ajustar_afc()` e o resto de `references/afc-ordinal.md` |
| Itens pela TCT | Dificuldade e discriminação de cada item | `analise_classica_itens()`; múltipla escolha: `analise_distratores()` (v1.6) |
| Itens pela TRI | Discriminação, informação, DIF | `tri_ajustar()`, `tri_dif()` (`references/tri.md`) |
| Precisão | Confiabilidade por fator | `confiabilidade()` no modelo da AFC; `indices_bifator()` se houver fator geral |
| Precisão no tempo e entre avaliadores | Estabilidade, concordância, mudança intraindividual | `confiabilidade_teste_reteste()`, `concordancia_juizes()`, `confiabilidade_ema()` |
| Precisão com várias fontes de erro (v1.6) | Quanto pesa cada faceta (itens, avaliadores, ocasiões) e quantas condições bastam | `generalizabilidade()` |
| Planejar o tamanho do teste (v1.6) | Quantos itens para chegar a uma confiabilidade | `spearman_brown()` |
| Normas (v1.7) | Como converter o escore bruto em percentil, z, T, padrão e estanino, por estrato | `normas()` |
| Itens pela TRI não paramétrica (v1.7) | Os itens formam uma escala de Mokken? | `escala_mokken()` (`references/tri.md`) |
| Validade de critério com medida imprecisa ou amostra selecionada (v1.6) | Quanto seria a correlação sem erro ou sem seleção | `corrigir_correlacao()` |

Antes dos itens, a matriz de especificação (blueprint): domínios do construto,
proporção de itens por domínio e formato de resposta, decididos pela definição
operacional e conferidos pelos juízes (Cook & Pitoniak, 2026, cap. 7, na pasta
Psicometria-R). A RVC de Lawshe também é apresentada por Martinková e Hladká (2023,
cap. 2).

Na adaptação, o professor não usa índice de concordância (avaliação qualitativa dos
juízes, porque o original já tem evidências; página "Avaliação da Síntese - Experts").
A planilha "Ficha - Validade de Conteúdo [Adaptação]" do próprio curso sugere avaliação
estruturada; quem quiser quantificar as equivalências (semântica, idiomática,
experiencial, conceitual) pode usar o CVC ou o IVC por critério de equivalência.

## 2. Validade de conteúdo

Formato dos dados: um juiz por linha, um item por coluna, notas numéricas. Juiz com nota
faltando sai por deleção listwise (contado em `n_excluidos`). Para o CVC com três
critérios, formato longo: uma coluna diz o critério (`criterio = "criterio"`), e cada
juiz aparece uma vez por critério.

| Índice | Escala | Cálculo | Corte padrão | De onde |
|---|---|---|---|---|
| CVC (Hernández-Nieto) | 1 a 5 | CVCi = média / 5; Pe = (1/J)^J; CVCc = CVCi − Pe; total = média dos CVCi − média dos Pe | CVCc ≥ .80 por item e critério | curso, página do CVC |
| I-IVC | 1 a 4 | juízes com 3 ou 4 / J | ≥ .80 (curso); ao lado, Polit e Beck (2006): 1,00 com 3 a 5 juízes, ≥ .78 com 6 ou mais | curso, página do IVC |
| S-IVC/Ave e S-IVC/UA | | média dos I-IVC; proporção de itens com I-IVC = 1 | Ave ≥ .90 (planilha do curso; Polit & Beck, 2006); UA sem corte | idem |
| Kappa modificado (k*) | | (I-IVC − Pc)/(1 − Pc), Pc = C(J, A)·.5^J | > .74 excelente, .60 a .74 bom, .40 a .59 regular | Polit, Beck e Owen (2007) |
| RVC (Lawshe) | 1 a 3 | (Ne − J/2)/(J/2), Ne = notas "essencial" | valor crítico binomial exato por J | Ayre e Scally (2014); curso |
| V de Aiken | mín. a máx. | Σ(nota − mín.)/(J·(c − 1)) com IC escore | sem corte nas fontes (argumento `corte`) | Aiken (1985); Penfield e Giacobbi (2004) |

O que as funções decidem:

- CVC: item adequado só se CVCc ≥ corte em todos os critérios (`$itens_adequados`).
  Com 2 juízes Pe = .25 e nem notas máximas passam (.75): a função avisa; o curso pede
  pelo menos 3. Pe com 3 juízes = .0370; com 4 = .0039 (a transcrição da aula traz
  .039, erro). Diagnóstico do curso: clareza baixa aponta redação; pertinência baixa,
  utilidade ou redundância; relevância baixa, representação do construto.
- IVC: as duas decisões saem lado a lado (`adequado_curso`, `adequado_polit_beck`). Elas
  discordam com 3 a 5 juízes: 4 de 5 (.80) passa no curso e não em Polit e Beck. O
  curso decide item a item e não usa o total (a média esconde item ruim).
- RVC: `rvc_critico(J)` aplica o método de Ayre e Scally: menor Ne com P(X ≥ Ne) ≤ .05
  numa binomial com p = .5 (unilateral). Reproduz a tabela da aula (5 a 10 juízes: Ne 5,
  6, 7, 7, 8, 9; RVC 1,00, 1,00, 1,00, .75, .78, .80). Com menos de 5 juízes não existe
  valor crítico (a função devolve NA e avisa). Acima de 10 juízes os valores saem do
  mesmo método e não foram conferidos contra a tabela impressa do artigo: conferir
  antes de citar. Efeito resumo: IVC de Lawshe, média do RVC dos itens mantidos.
- Aiken: o IC escore é o intervalo de Wilson da proporção V com J·(c − 1) ensaios (os
  evals conferem contra `prop.test`). Comparar o limite inferior com um corte declarado.

Relevância (IVC) e essencialidade (RVC) são perguntas diferentes: item correto mas
redundante passa no IVC e cai no RVC. Nenhum índice substitui os comentários dos juízes,
que dizem o que mudar, e nenhum garante carga fatorial, unidimensionalidade ou
compreensão pela população-alvo.

Piloto com população-alvo: mesma fórmula do CVC com N participantes (critérios clareza,
adequação e compreensão). O termo (1/N)^N foi derivado para juízes; aplicado a 30
participantes ele some, e o CVC vira média / 5 (extensão do curso, sem respaldo
formal). A pergunta "o item precisa ser modificado?" (sim/não) fica fora do CVC: é uma
proporção por item, que o curso lê com cautela abaixo de .90 (heurística do
professor).

## 3. Análise fatorial exploratória

Fluxo: `auditar_itens()` (categorias raras, assimetria) → `afe_adequacao()` →
`afe_numero_fatores()` → `afe()` → AFC em outra amostra. `afe()` sozinha roda as três
etapas e junta as frases.

Tipo de correlação. Itens com até 7 valores distintos são ordinais e entram com
correlação policórica (`psych::polychoric`); acima disso, Pearson. É a mesma regra de
`.pr_e_ordinal()` usada por `ajustar_afc()` e pelo `cor_auto`, então AFE e AFC tratam os
itens do mesmo jeito. Escala 0-100 (PBAT) vai para Pearson.

Antes de qualquer correlação, as três funções da AFE conferem os dados brutos e param
com mensagem em português quando a matriz seria singular por construção: item constante;
n menor ou igual ao número de itens; dependência linear exata entre colunas (posto do
`qr` dos dados padronizados menor que o número de itens, ou menor autovalor da Pearson
abaixo de 1e-8), com os nomes das colunas que são combinação exata das demais. Causas
habituais: escore total ou de subescala junto dos itens, item duplicado, item junto da
sua versão invertida, colunas derivadas (resíduos, médias). Suavizar uma matriz dessas
daria KMO e cargas sem sentido. Só a policórica, estimada par a par, pode sair não
positiva definida sem dependência exata nos dados (comum em amostra pequena ou com
categoria rara); nesse caso ela é suavizada com `psych::cor.smooth` (o mesmo que o
`psych::polychoric` faria por padrão), com aviso, e a suavização entra no relato.

Adequação (`afe_adequacao`). KMO geral e MSA por item calculados pela matriz
anti-imagem (iguais ao `psych::KMO` nos evals), com os rótulos de Kaiser (1974): .90
maravilhoso, .80 meritório, .70 mediano, .60 medíocre, .50 miserável, abaixo de .50
inaceitável. Item com MSA < .50 gera aviso (retirar um por vez e recalcular). Bartlett:
χ² = −(n − 1 − (2p + 5)/6)·ln|R|, gl = p(p − 1)/2 (igual ao `psych::cortest.bartlett`);
com n grande quase sempre rejeita, e sobre a policórica é aproximado: pesa pouco. A
frase APA só diz que a fatorabilidade "foi indicada" com KMO ≥ .60 e Bartlett com
p < .05 (convenção da skill: a faixa miserável de Kaiser, .50 a .60, não basta); nos
demais casos diz que não foi indicada, com o KMO, o rótulo de Kaiser e o Bartlett, e a
função avisa (KMO abaixo de .50, KMO na faixa .50 a .60, Bartlett sem rejeitar a matriz
identidade). `$fatoravel` guarda a decisão.

Número de fatores (`afe_numero_fatores`). Análise paralela (Horn, 1965) com dados
permutados: cada coluna é embaralhada de forma independente (Buja & Eyuboglu, 1992), o
que mantém a distribuição de cada item e desfaz as correlações, e a correlação aleatória
é do mesmo tipo da observada (policórica com policórica). Retém fatores enquanto o
autovalor observado passa do percentil 95 dos aleatórios (Glorfeld, 1995). Duas bases
saem sempre: autovalores da matriz de correlação ("componentes", o critério original de
Horn e o principal por padrão) e da matriz reduzida com SMC na diagonal ("fatores").
Quando discordam, a função avisa: decidir pela teoria e pela leitura das soluções
vizinhas e relatar as duas. Padrões: 100 reamostras (os evals usam 20), semente 2026
(restaura o gerador do usuário ao sair). Com policórica, cada reamostra recalcula a
matriz: 100 reamostras com 24 itens levam cerca de um minuto. `ega = TRUE` acrescenta a
EGA (EGAnet, padrões do pacote), que o curso indica como complemento à análise paralela
(página da EGA); estabilidade com bootEGA em `references/redes-ega.md`.

Hull e MAP (v1.4). `hull = TRUE` (padrão) acrescenta o método de Hull (Lorenzo-Seva et
al., 2011): ajusta soluções de 0 a J fatores por máxima verossimilhança (`stats::factanal`
sobre a mesma matriz), registra o ajuste e os graus de liberdade, descarta as soluções
que não melhoram o ajuste e as que ficam abaixo da reta entre as vizinhas (contorno
convexo superior; Ceulemans & Kiers, 2006) e escolhe a de maior razão st, o cotovelo do
scree formalizado. J = número da paralela (ou `n_teorico`, se maior) + 1, limitado pelo
limite de Ledermann e com mínimo de 3: a paralela fixa o teto, contra superfatoração.
Índice principal: CFI (com χ² da máxima verossimilhança; Hull-CFI); secundário: CAF, a
parte comum explicada (1 − KMO da matriz residual; Hull-CAF), o índice indicado quando
não há χ². Os evals conferem CFI, st e a escolha com o `EFAtools::HULL` (tol. 1e-6) e o
CFI com a conta à mão. Com policórica o χ² é aproximado, mas o Hull usa só a posição
relativa das soluções; o capítulo 11 do Wiley Handbook mostra o Hull-CFI do modelo
ordinal concordando com o Hull-CFI por ULS sobre policóricas. `map = TRUE` (padrão)
acrescenta o MAP de Velicer (1976) e o revisado com a quarta potência (Velicer et al.,
2000), iguais ao `psych::vss` nos evals.

O que o capítulo 11 do Wiley Handbook (Timmerman et al., 2018), na pasta Psicometria-R,
recomenda: considerar Hull e paralela, cada um na variante coerente com o modelo e o
estimador, e decidir pela interpretabilidade e pela plausibilidade teórica das
soluções indicadas; o Hull foi melhor que a paralela com amostras grandes e muitos
itens por fator, e a paralela melhor com N = 100 e 5 itens por fator (a função lembra
isso quando a amostra ou a razão itens/fator caem nessa faixa); o MAP foi desaconselhado
para fatores comuns, porque piora com unicidades altas (Lorenzo-Seva et al., 2011), e
por isso entra como critério auxiliar. Discordância entre paralela e Hull vira aviso.
A PA-MRFA (Timmerman & Lorenzo-Seva, 2011), preferida pelo capítulo para itens
ordinais, entra com `mrfa = TRUE` (v1.6; pacote EFA.MRFA, opcional): compara a % de
variância comum de cada fator, pela análise fatorial de posto mínimo, com a do acaso,
sobre policóricas quando os itens são ordinais. Discordância entre PA-MRFA e paralela
vira aviso, e o capítulo pesa mais a PA-MRFA com itens ordinais e cargas moderadas.

Extração e rotação (`afe`).

- `extracao = "minres"` (padrão): `psych::fa` sobre a matriz policórica ou de Pearson.
  Devolve o RMSR (resíduos fora da diagonal); os índices de ajuste do psych para
  minres sobre policórica são aproximados e ficam de fora.
- `extracao = "wlsmv"`: `lavaan::efa`, com WLSMV para itens ordinais e MLR para
  contínuos (os mesmos estimadores da AFC). Traz CFI, TLI, RMSEA (IC 90%) e SRMR, na
  versão robusta quando existe, pela `tabela_ajuste()` de `psicometria.R`. Com WLSMV os
  cortes de Hu e Bentler (1999) tendem a não acusar desajuste (Xia & Yang, 2019; ver
  `references/afc-ordinal.md`). Nos evals as cargas das duas extrações diferem menos de
  .05.
- Rotação oblíqua por padrão (`oblimin`, via GPArotation): fatores psicológicos raramente
  são independentes, e quando são a oblíqua devolve o mesmo que a ortogonal. Rotação
  ortogonal pedida pelo usuário gera aviso. Com 1 fator não há rotação.
- Fatores ordenados pela variância explicada e com sinal tal que a soma das cargas seja
  positiva. Variância por fator = soma de λ·(λΦ) da matriz padrão, que soma a
  comunalidade total mesmo com fatores correlacionados (é o "Sum of sq (obliq)
  loadings" do lavaan).
- Tabela por item: cargas do padrão, h², u², fator principal e alerta. Alertas
  (convenção da skill; o curso não dá corte): "carga cruzada" com |λ| ≥ .30 em dois ou
  mais fatores; "carga baixa" com maior |λ| < .40 (a mesma leitura de .40 da AFC em
  `references/afc-ordinal.md`); "Heywood" com h² ou |λ| ≥ 1, ou quando a extração
  avisa caso Heywood ou variância negativa. Os cortes são argumentos (`cruzada_min`,
  `carga_min`).
- Caso Heywood torna a solução imprópria: a função avisa com os itens envolvidos, marca
  `$solucao_impropria` e a frase APA diz que a solução é imprópria, com os números entre
  parênteses só para rastreio, em vez de descrever as cargas como estimativas.
- Exploratória não é sem teoria (página "Estrutura Interna"): aceitar uma solução só se
  os agrupamentos fizerem sentido para o construto; não aceitar qualquer solução do
  software.

Privacidade: tudo o que sai é por item, por fator ou matriz de correlações. O objeto do
modelo (`$modelo`, lavaan ou psych) fica no resultado e fora do print; o do lavaan
guarda os dados, então não enviar o objeto inteiro para fora do computador.

## 3b. Análise clássica de itens (v1.4)

`analise_classica_itens(dados, itens)`, pela TCT de Martinková e Hladká (2023, cap. 5,
na pasta Psicometria-R): dificuldade (proporção de acerto ou endosso em dicotômicos;
(média − mínimo)/(máximo − mínimo) em ordinais), RIT (correlação item-total), RIR
(item-resto, sem o próprio item no total), ULI (diferença de dificuldade entre o terço
superior e o inferior do escore total) e alfa sem o item. Régua de Ebel (1954), pelo
livro: ULI por terços ≥ ,20, salvo item muito fácil ou muito difícil; ULI baixo em item
de dificuldade intermediária vira aviso, e RIR negativo denuncia item invertido sem
recodificar. RIT, RIR e alfa sem o item iguais ao `psych::alpha` nos evals. Rodar por
subescala quando a estrutura é multidimensional; alfa que sobe sem um item não decide
exclusão sozinho.

## 3c. Fidedignidade além da consistência interna (v1.4)

**Teste-reteste**: `confiabilidade_teste_reteste(dados, "t1", "t2")`. O índice
principal é o ICC(A,1), duas vias, concordância absoluta, medida única (McGraw & Wong,
1996), que Koo e Li (2016) indicam para reteste; o ICC(C,1) e a Pearson saem ao lado. O
curso usa a Pearson, que não enxerga mudança sistemática entre as ocasiões: a diferença
entre ICC(C,1) e ICC(A,1) e o teste da mudança média mostram essa mudança. Saem também o
erro padrão de medida (EPM = DP de todas as medidas × raiz(1 − ICC)), a mudança mínima
detectável a 95% (1,96 × raiz(2) × EPM; Weir, 2005) e os limites de concordância de
Bland e Altman (1986). Régua de Koo e Li aplicada ao IC, como eles recomendam: abaixo de
,50 fraca, ,50 a ,75 moderada, ,75 a ,90 boa, acima de ,90 excelente; n < 30 vira aviso.
O relato diz o intervalo entre as ocasiões e se o construto é traço ou estado.

**Concordância entre juízes ou avaliadores**: `concordancia_juizes(dados, avaliadores,
tipo)`. Nominal: kappa de Cohen (1960) com 2 avaliadores, Fleiss (1971) com 3 ou mais
(IC por bootstrap dos casos); ordinal: kappa ponderado quadrático (Cohen, 1968) com 2;
intervalar: ICC(A,1) e ICC(A,k). Régua de Landis e Koch (1977) declarada como régua. O
curso usa o kappa na triagem de revisões sistemáticas e o kappa de Fleiss e o ICC na
análise de juízes de desejabilidade social. Conferido nos evals com o irr, com os
valores publicados de Shrout e Fleiss (1979) e com o exemplo resolvido do kappa de
Fleiss (κ = ,210).

Desde a v1.6 a função devolve também, fora do índice principal, o AC1 de Gwet (2008)
(nominal; AC2 com pesos quadráticos quando `tipo = "ordinal"`) e o alfa de Krippendorff
(Hayes & Krippendorff, 2007) na métrica do tipo (nominal, ordinal pelos postos acumulados,
intervalar), os dois com IC por bootstrap dos casos. O motivo é a aula "Limitações do
coeficiente Kappa" do curso: com quase todos os artigos incluídos (ou excluídos), a
concordância esperada pelo acaso fica alta e o kappa despenca mesmo com concordância de
90% (exemplo da aula: 9 sim-sim e 1 sim-não dá κ = 0), e com todos na mesma categoria o
kappa nem é calculável. Esse é o paradoxo de Feinstein e Cicchetti (1990); a função avisa
quando a concordância simples passa de 80% e o κ fica abaixo de ,40 (ou é incalculável),
e manda relatar o AC1 e a concordância simples ao lado. Para que o AC1 enxergue a
categoria que ninguém usou, passar os juízes como fator com os níveis possíveis (ou
`categorias = c("sim", "não")`). AC1, AC2 e Krippendorff nominal e intervalar iguais ao
irrCAC nos evals; o Krippendorff ordinal conferido pela identidade com o intervalar
sobre os postos médios (o `irr::kripp.alpha` dá valores levemente diferentes do irrCAC
com 3 ou mais juízes e não é usado como referência).

**Teoria da generalizabilidade** (v1.6): `generalizabilidade(dados, pessoa, escore,
facetas)`, com uma faceta (pessoa × item, ou pessoa × avaliador) ou duas (pessoa × item
× ocasião), delineamento cruzado e uma observação por célula; dados longos, ou largos com
`itens =` para uma faceta. Estudo G: componentes por REML (`lme4`), como Martinková e
Hladká (2023, 4.4 e 4.5) e Petersen (2024, cap. 5); o resíduo absorve a interação de
maior ordem. Estudo D: Eρ² = σ²p/(σ²p + σ²rel), para ordenar pessoas, e Φ =
σ²p/(σ²p + σ²abs), para decidir por ponto de corte (soma os efeitos principais das
facetas ao erro); `estudo_d = list(item = c(5, 10, 20), ocasiao = 1:3)` projeta outros
tamanhos, a generalização da profecia de Spearman-Brown. Âncoras nos evals: o exemplo
pessoa × item × ocasião de Petersen (Eρ² = ,8731; Φ = ,6374, iguais ao `gtheory`, pacote
arquivado no CRAN que a skill não exige) e o HCI de Martinková (σ²p = ,0237, σ²i = ,0270,
σ²e = ,1885; Φ = ,6874 com 20 itens). Delineamentos aninhados (itens dentro de pessoa,
avaliadores diferentes por pessoa) ficam fora; para EMA, `confiabilidade_ema()`.

**Profecia de Spearman-Brown** (v1.6): `spearman_brown(confiabilidade, k_atual, k_novo)`
projeta a confiabilidade com outro número de itens, e `alvo =` devolve quantos itens
paralelos levam à confiabilidade desejada (arredondado para cima). Exemplo do livro de
Martinková e Hladká (2023, 4.2.5.1): 20 itens com ,75 levados a 30 dão ,8182; para ,85
são 38 itens. Supõe itens novos tão bons quanto os atuais.

**Correção de correlação** (v1.6): `corrigir_correlacao(r, confiabilidade_x,
confiabilidade_y, dp_irrestrito, dp_restrito, n)`. Atenuação: r/sqrt(rxx ryy), o que
seria a correlação entre escores verdadeiros (Martinková & Hladká, 2023, 4.2.3);
restrição de amplitude por seleção direta em X: caso 2 de Thorndike (1949), como no
`psych::rangeCorrection`; exemplo do livro (2.4.3): candidatos admitidos com X > 58, r
restrito ,3633, corrigido ,6548. Com as duas correções, a restrição vem primeiro. Valor
acima de 1 avisa (confiabilidade subestimada ou erros correlacionados). Correção serve
para interpretar validade de critério; o teste de hipótese continua sobre o r observado.

**Normas** (v1.7): `normas(dados, escore, grupos = NULL, percentil = "medio")` monta,
para cada escore bruto, o percentil, z (média e DP do estrato), T = 50 + 10z, escalar =
10 + 3z, padrão = 100 + 15z e estanino pelas faixas cumulativas 4, 11, 23, 40, 60, 77, 89
e 96% (Petersen, 2024, cap. 1, tabela 1.1). Percentil `"medio"` (padrão: abaixo + metade
dos iguais, sobre n, que não dá 0 nem 100 no extremo) ou `"petersen"` (o
`dplyr::percent_rank` do livro). `grupos = c("faixa_etaria", "escolaridade")` dá norma
por estrato, com aviso abaixo de 30 casos; `normalizar = TRUE` acrescenta z e T
normalizados (qnorm do percentil) para escore assimétrico; `valores =` monta a tabela
para escores que a amostra não teve. A tabela é agregada (sem frequência por escore no
print). Converter o escore de um paciente é consulta à tabela, feita no computador. A
norma descreve a população que a amostra representa: o relato diz como a amostra
normativa foi composta.

**Análise de distratores** (v1.6): `analise_distratores(dados, itens, gabarito)`, para
itens de múltipla escolha com as respostas em letras (Martinková & Hladká, 2023, 5.5):
proporção que escolhe cada opção em cada grupo do escore (terços, como o
ShinyItemAnalysis) e discriminação da opção (grupo superior menos inferior). Gabarito que
não cresce com o escore, distrator que cresce e distrator quase nunca escolhido viram
aviso.

**Medidas intensivas (EMA, diário)**: `confiabilidade_ema(dados, id, tempo, itens)`,
pela teoria da generalizabilidade de Shrout e Lane (2012), via
`psych::multilevel.reliability`: RkF (confiabilidade entre pessoas da média de todos
os tempos), R1R e RkR (tempos sorteados) e Rc, a confiabilidade da mudança dentro da
pessoa (Cranford et al., 2006), que é a que importa para análise intraindividual
(redes temporais, mudança no PBAT). Dados em formato largo (uma linha por pessoa e
ocasião); com dados desbalanceados os componentes saem do lme4.

## 4. Depois da AFE

- AFC em outra amostra com o modelo sugerido: `ajustar_afc()`, `tabela_ajuste()`. AFE e
  AFC na mesma amostra não confirmam nada. Metade aleatória do mesmo banco também não
  (aulas "Devo dividir a minha amostra?" e "Problemas em dividir a amostra"): as metades
  carregam os mesmos vieses de recrutamento, a confirmação tende a repetir a exploratória e
  cada análise perde metade do poder. Sem amostra independente, relatar a AFE (ou a AFC
  de um modelo teórico) e deixar a replicação como limitação; se dividir mesmo assim,
  dizer que é validação cruzada interna, não confirmação.
- Precisão por fator, depois da estrutura (ordem do curso: estrutura, fatores, cargas,
  ajuste, confiabilidade): `confiabilidade()` dá ômega e alfa por fator. O curso usa
  confiabilidade ≥ .70 como referência não mecânica (página "Consistência Interna") e
  prefere a confiabilidade composta; com itens ordinais a skill relata o ômega do escore
  observado (Green & Yang, 2009), como em `references/afc-ordinal.md`. Não calcular um
  coeficiente único para escala multidimensional.
- Validade convergente e discriminante: `validade_discriminante()` (AVE, HTMT,
  Fornell-Larcker). O curso define discriminante como correlação nula com construto
  não relacionado; a skill oferece os critérios de MEE e documenta a diferença.
- Teste-reteste: o curso usa a correlação de Pearson entre T1 e T2; a skill relata o
  ICC(A,1) com EPM (`confiabilidade_teste_reteste()`, seção 3c), porque Pearson ignora
  mudança sistemática de nível entre as aplicações.
- Duas metades (aula "Fidedignidade Split-Half"): `confiabilidade_metades()` (v1.5) dá,
  para ímpares contra pares ou primeira contra segunda metade, a correlação entre as
  metades, Spearman-Brown e Guttman λ4 (principal, porque não supõe metades com a mesma
  variância); com `metodo = "todas"`, média, mínimo e máximo de λ4 em todas as divisões ao
  meio (com número par de itens a média é o alfa, Cronbach, 1951; o máximo superestima com
  muitos itens). O ômega segue como estimativa principal de consistência interna.
- Evidências baseadas nos processos de resposta (aula "Evidências de validade baseadas no
  padrão de resposta"): entrevista cognitiva, grupo focal, rastreamento ocular e tempo de
  resposta. São evidências qualitativas ou de delineamento, sem função na skill; o tempo
  de resposta entra como critério de limpeza declarado antes da análise (respostas rápidas
  demais), com o número de casos excluídos no relato.
- PBAT: os autores tratam o instrumento como banco de itens, não como escala; AFE no
  PBAT entra como descrição da organização por valência, em Pearson (0-100), nunca como
  validação de escore somado (`references/pbat-mpfi.md`).

## 5. Como relatar

```r
source("~/.local/share/psicometria-r/carregar.R")
validade_conteudo_cvc(juizes, itens, criterio = "criterio")   # formato longo
validade_conteudo_ivc(juizes_relevancia, itens)
validade_conteudo_rvc(juizes_essencialidade, itens)
res <- afe(dados, itens)               # adequação + análise paralela + extração
res <- afe(dados, itens, n_fatores = 3, extracao = "wlsmv")
relatar(res)
```

Frases geradas (dados simulados dos evals, 12 itens de 4 categorias, n = 600):
"A fatorabilidade da matriz de correlações policóricas foi indicada pelo KMO = .83 (MSA
dos itens de .81 a .85) e pelo teste de Bartlett, χ²(66) = 2180.45, p < .001." e "A
análise fatorial exploratória (extração por resíduos mínimos sobre correlações
policóricas, rotação oblimin) com 3 fatores explicou 49.0% da variância; as cargas
principais variaram de .64 a .73, sem cargas cruzadas; as correlações entre fatores
variaram de .28 a .37." Com `extracao = "wlsmv"` a frase termina com CFI, TLI, RMSEA e
SRMR. Quando a matriz não é fatorável (seis variáveis normais independentes, n = 60, dos
evals): "A fatorabilidade da matriz de correlações de Pearson não foi indicada: KMO = .45
(inaceitável pelos rótulos de Kaiser, 1974; MSA dos itens de .36 a .57) e teste de
Bartlett, χ²(15) = 14.94, p = .455, que não rejeitou a matriz identidade." Com caso
Heywood (três itens com r12 = r13 = .80 e r23 = .50, um fator, dos evals): "A análise
fatorial exploratória (extração por resíduos mínimos sobre correlações de Pearson) com 1
fator resultou em solução imprópria, com caso Heywood (carga ou comunalidade igual ou
acima de 1) em h1; os valores dessa solução (71.8% da variância; cargas principais de
.74 a 1.03) não devem ser interpretados como estimativas." CVC: "A validade de conteúdo foi avaliada por 3 juízes numa escala de 1 a 5. O
coeficiente de validade de conteúdo total (CVC; Hernández-Nieto, 2002) foi clareza =
.76, pertinencia = .79, relevancia = .87; 1 de 3 itens alcançaram CVCc ≥ .80 em todos
os critérios." (os nomes dos critérios vêm da coluna dos dados).

O método do manuscrito declara: número e perfil dos juízes, escala e critérios, índice e
corte com a fonte; para a AFE, tipo de correlação e por quê, KMO e Bartlett, critério
de retenção (análise paralela, base, reamostras, percentil, semente), extração, rotação,
cortes de carga e o que foi feito com itens problemáticos.

## 6. Armadilhas

1. Planilha de CVC da aula: o total de clareza usa L4:L17 (exclui o item 1), pertinência
   e relevância usam todas as linhas. A função usa todos os itens.
2. Planilha de IVC: o texto diz "delineada para 3 juízes" e a fórmula divide por 5. A
   função divide pelo número real de juízes.
3. Cortes diferentes por fonte: CVC .80 no curso, .70 aceito por alguns autores (o
   procedimento é o mesmo usado por Cassepp-Borges et al., 2010); I-IVC .80 no curso,
   1,00 ou .78 em Polit e Beck (2006). Declarar o corte e a fonte.
4. Ayre e Scally (2014) recalcularam pelo método exato os valores críticos da tabela
   original de Lawshe (1975), e é essa versão que a aula usa. Usar `rvc_critico()` e
   citar Ayre e Scally.
5. S-IVC/UA cai quando o número de juízes cresce (um juiz desfavorável basta).
6. KMO e Bartlett não dizem quantos fatores nem se a estrutura é a esperada.
7. Bases da análise paralela podem discordar; relatar as duas e decidir pela teoria.
8. Retirar item por carga baixa ou cruzada é decisão do autor, um de cada vez, com a
   AFE refeita a cada retirada e relatada.
9. Escore total ou de subescala no mesmo banco dos itens (ou item duplicado, item e sua
   versão invertida, colunas de resíduos): a matriz fica singular. As funções da AFE param
   e nomeiam a coluna; passar só os itens em `itens`.
10. Caso Heywood (carga ou comunalidade ≥ 1): solução imprópria, que não se relata como
    estrutura. Rever o número de fatores (em geral, fatores demais), itens quase
    redundantes ou o tamanho da amostra.

## 7. Fontes e divergências

Livro do usuário: os sete volumes não tratam de validade de conteúdo nem de AFE (busca
nas extrações; o Vol. 7 cobre redes). Os padrões deste módulo vêm do curso, e onde o
curso é omisso, da literatura ou de convenção da skill declarada como tal.

Curso (Psicometria Online Academy): páginas "Construção de Medidas de Autorrelato",
"Índice de Validade de Conteúdo (IVC)", "Razão de Validade de Conteúdo (RVC)",
"Coeficiente de Validade de Conteúdo (CVC)", "Coeficiente de Validade de Conteúdo com
População-Alvo", "Avaliação da Síntese - Experts", "Estrutura Interna", "Consistência
Interna", "Confiabilidade Teste-Reteste" e a página da EGA. Números reproduzidos nos
evals: CVC (5, 5, 4) = .896 e (4, 4, 4) = .763; Pe de 2, 3 e 4 juízes; IVC .67, .80 e
.40 dos exemplos; tabela de RVC crítico de 5 a 10 juízes e os RVC .80, .60 e −.60 com 10
juízes. Erros da fonte não reproduzidos: Pe de 4 juízes ".039" (correto .0039); faixas
das planilhas (armadilhas 1 e 2). O curso não dá corte para cargas, KMO ou retenção de
fatores; os cortes usados aqui são da literatura (KMO) ou da skill (cargas).

Divergências decididas: IVC com corte do curso por padrão e Polit e Beck ao lado; CVC
com .80 do curso; RVC pelo método exato (que reproduz a aula) em vez de tabela fixa;
análise paralela com base em componentes por padrão e a base de fatores ao lado;
discriminante pelos critérios de MEE além da definição do curso (já em `psicometria.R`).

Literatura (DOIs conferidos no Crossref, referência completa em `referencias-verificadas.md`):
Polit e Beck (2006); Polit, Beck e Owen (2007); Lynn (1986); Lawshe (1975); Ayre e
Scally (2014); Aiken (1985); Penfield e Giacobbi (2004); Kaiser (1974); Bartlett (1950);
Horn (1965); Buja e Eyuboglu (1992); Glorfeld (1995). Sem DOI: Hernández-Nieto (2002,
livro) e Cassepp-Borges, Balbinotti e Teodoro (2010, capítulo), citados pelo curso. Já
na lista da skill: Golino e Epskamp (2017), Green e Yang (2009), Hu e Bentler (1999),
Xia e Yang (2019), International Test Commission (2018). Pacotes: psych, GPArotation,
lavaan, EGAnet (`citation()` e `versoes_pacotes()` no R da sessão).
