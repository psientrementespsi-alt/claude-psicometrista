# Invariância de medida

## Sumário

- 16 (linha 27)
- 70 (linha 81)
- 108 (linha 119)
- 124 (linha 135)
- 152 (linha 163)
- 174 (linha 185)
- 183 (linha 194)
- 191 (linha 202)

Índice: Sequência · Critérios · Quando um passo falha · Invariância longitudinal · Comparação entre grupos sem invariância escalar · Viés de medida, DIF e viés preditivo

Função: `invariancia(modelo, dados, itens, grupo, ordinal = NULL, estimador =
NULL)`, em `scripts/psicometria.R`. Sem `ordinal`, item com até 7 valores
distintos conta como ordinal (WLSMV) e o resto como contínuo (MLR). Para a
invariância usar o mesmo estimador da AFC, passar o que `sugerir_estimador()`
devolveu: `invariancia(modelo, dados, itens, grupo, ordinal = est$ordinal,
estimador = est$estimador)`; com MLR, a sequência é a contínua. A função não
repassa `group.partial`: invariância parcial se monta com `measEq.syntax()`
direto (seção abaixo). Devolve `tabela` (um modelo por linha, com CFI, TLI,
RMSEA, SRMR, Δ de cada um e p do teste de diferença de χ² escalonado contra o
passo anterior) e `fits` (os objetos lavaan).

## Sequência

Itens ordinais, com a identificação de Wu e Estabrook (2016):
configural, limiares, cargas (com limiares), interceptos (com limiares e
cargas). Itens contínuos ou quase contínuos (0-100): configural, métrica
(cargas), escalar (cargas e interceptos), estrita (mais resíduos). Na
correspondência com o vocabulário do modo PSICOMETRICO da fabrica-de-artigo,
"cargas" é a métrica e "interceptos" é a escalar.

Por que limiares primeiro no caso ordinal. Wu e Estabrook (2016) mostram que
as condições usuais de identificação do modelo configural fixam a escala das
respostas latentes de jeitos diferentes, e que, ao impor invariância, essas
condições podem virar restrições e gerar modelos não equivalentes, com
conclusões diferentes; eles derivam condições de identificação que não
dependem da parametrização escolhida para o modelo de base. No
`semTools::measEq.syntax()`, que implementa isso com
`ID.cat = "Wu.Estabrook.2016"` (o padrão), a documentação diz que os autores
recomendaram igualar os limiares primeiro e que isso permite liberar
restrições de identificação que deixam de ser necessárias, a começar pelo
intercepto da resposta latente de cada item nos grupos além do primeiro.

Consequência prática, mostrada nos evals da skill: quando todos os limiares de
um item se deslocam juntos num grupo (deslocamento uniforme, o equivalente
ordinal de um intercepto diferente), o passo dos limiares passa e a falha
aparece no passo dos interceptos (no eval: p = .85 nos limiares, ΔCFI = -.113
e p < .001 nos interceptos). Quando muda o espaçamento entre os limiares, a
falha aparece já no passo dos limiares. Ler a tabela com isso em mente: passo
dos limiares aprovado não quer dizer que as categorias funcionam igual nos
grupos até os interceptos serem testados.

`invariancia()` usa `ID.fac = "std.lv"` (consistente com Wu e Estabrook,
segundo a mesma documentação), `parameterization = "theta"` e, com itens
ordinais, `missing = "pairwise"`, o mesmo da `ajustar_afc()`.

Itens com 2 ou 3 categorias. Pela documentação do `measEq.syntax()`, com item
binário não há teste independente dos limiares, e com 3 categorias a igualdade
dos limiares é pressuposta. Nos evals, com 3 categorias, o passo dos limiares
não ganha grau de liberdade; a função deixa o p desse passo em branco (NA) e
avisa, em vez de produzir um teste vazio.

Modelo com fator de ordem superior (o MPFI-24, por exemplo). A função recusa:
o `measEq.syntax()` troca a identificação para `ID.fac = "ul"` quando detecta
fator de ordem superior, e os passos de cargas e interceptos só alcançam os
indicadores observados (cargas de segunda ordem pedem "regressions" e
interceptos de fatores latentes, "means"). O caminho: testar a medida no
modelo de fatores correlacionados de primeira ordem com `invariancia()` e,
para o nível de segunda ordem, montar à mão os passos com `group.equal`
incluindo "regressions" e depois "means", conferindo os argumentos na ajuda do
`measEq.syntax()`. Esse segundo caminho não está coberto pelos evals.

Diferença real de média latente não é falta de invariância: nos evals, com
média latente .5 maior num grupo e itens invariantes, nenhum passo rejeita
(p > .01, |ΔCFI| ≤ .010).

## Critérios

Usar três informações juntas e relatar as três:

1. **ΔCFI.** Cheung e Rensvold (2002) examinaram 20 índices em simulação de dois
   grupos e recomendaram ΔCFI, ΔGamma hat e ΔMcNCI para avaliar invariância,
   por serem independentes da complexidade do modelo e do tamanho da amostra.
   O critério de queda de .01 no CFI é o mais usado (Putnick & Bornstein, 2016,
   seção sobre critérios de ajuste).
2. **ΔRMSEA e ΔSRMR.** Chen (2007), em duas simulações, achou o SRMR mais
   sensível a falta de invariância nas cargas do que nos interceptos ou
   resíduos, CFI e RMSEA igualmente sensíveis aos três níveis, e mudanças nos
   índices maiores com grupos de tamanho igual. Os valores propostos (queda de
   .01 no CFI junto com aumento de .015 no RMSEA, e de .030 no SRMR para
   cargas ou .015 para interceptos e resíduos) estão aqui pela síntese de
   Putnick e Bornstein (2016), que os apresenta para amostra com poder
   adequado, grupos de tamanho igual e não invariância mista (umas cargas
   maiores e outras menores num grupo). Grupos clínico e não clínico costumam
   ter tamanhos bem diferentes; nesse caso, dizer que a condição não vale.
   Conferir no original antes de citar os números como de Chen (2007)
   [âncora secundária].
3. **Teste de diferença de χ² escalonado** (`lavTestLRT`, coluna `p_dif_qui2`).
   Sensível a diferenças pequenas com n grande, motivo pelo qual parte dos
   pesquisadores passou a usar os Δ de índices (Putnick & Bornstein, 2016), mas
   continua sendo a única peça com distribuição conhecida. Com três passos a
   .05, uma rejeição ao acaso acontece perto de 14% das vezes quando os grupos
   são invariantes; olhar o tamanho da diferença nos parâmetros antes de
   concluir.

Ressalva para WLSMV. Os critérios de Δ foram derivados com dado contínuo e ML.
Savalei (2021) mostra que os índices categóricos correntes indicam ajuste
melhor que o dos mesmos dados sem categorização; por extensão (inferência
desta skill, não achado de um estudo de invariância), os Δ derivados com dado
contínuo não têm calibração garantida para WLSMV. Tratar como heurística,
usar a versão robusta dos índices em todos os passos (`invariancia()` faz isso
quando ela existe em todos) e dar peso ao teste de diferença e ao tamanho das
diferenças nos parâmetros.

## Quando um passo falha

Não seguir para o próximo passo como se nada tivesse acontecido: as restrições
se acumulam sobre o desajuste (no eval de espaçamento, a falha nos limiares
contaminou os passos seguintes). Putnick e Bornstein (2016) descrevem três
saídas: localizar a fonte liberando restrições uma a uma (para trás) ou
acrescentando (para frente) até chegar a um modelo de invariância parcial;
retirar os itens não invariantes e testar de novo; ou concluir que o
construto não é invariante e não comparar os grupos.

Invariância parcial com `measEq.syntax()`: o argumento `group.partial` recebe
os parâmetros liberados na sintaxe do lavaan, por exemplo
`group.partial = c("item3 | t1", "item3 | t2", "item3 | t3")` para limiares ou
`"f1 =~ item3"` para carga. Toda liberação vai relatada com a justificativa
(conteúdo do item, tradução, contexto cultural), não só o índice.

## Tamanho do efeito da não invariância: `dmacs()` (v1.6)

O teste diz se há diferença; o dMACS (Nye & Drasgow, 2011) diz quanto ela pesa na
unidade do item: a raiz da diferença quadrática média entre os escores esperados do item
nos dois grupos, para pessoas com o mesmo traço, ao longo da distribuição do traço no
grupo focal, dividida pelo DP combinado do item. Petersen (2024, 15.10.3.4) roda o
`dmacs::lavaan_dmacs` no configural; o Wiley Handbook (cap. 27) traz a fórmula e o ΔMean
de escala.

```r
inv <- invariancia(modelo, dados, itens, "grupo")
dmacs(inv$fits$configural, grupo_ref = "referencia")   # ou o modelo parcial
```

Saem, por item: `dmacs`, `delta_media` (focal menos referência, com sinal, na unidade do
item) e a classe pela régua citada por Nye et al. (2019): ,20 pequeno, ,40 médio, ,70
grande (lida em fonte secundária; conferir o artigo antes de usar como critério de
corte). Por fator: a soma dos `delta_media` (ΔMean), comparada com a diferença observada
do escore somado. Para indicadores contínuos a função usa a forma fechada do integral
(sem grade numérica), que bate com o `pinsearch::dmacs` (exato) e com o pacote `dmacs`,
arquivado no CRAN, a menos do DP combinado (o pacote combina os DPs linearmente; a
skill, as variâncias). Leitura: no configural com std.lv em todos os grupos o traço tem
média 0 e variância 1 nos dois, toda diferença vira diferença de item e o ΔMean reproduz
a diferença observada por construção (a função avisa e não calcula a proporção). Para
separar viés de diferença real no traço, rodar no modelo de invariância parcial (âncoras
iguais, médias latentes livres). Itens ordinais ficam fora: `tri_dif()` ou
`pinsearch::pin_effsize()`.

## Equivalência do ajuste e invariância ao longo de variável contínua (v1.7)

**Teste de equivalência** (`equivalencia_ajuste(qui2, gl, n, grupos, rmsea0 = .08)`;
Counsell et al., 2020, como o `equiv_chi()` que Petersen, 2024, 15.10.3.5, usa). Inverte
a pergunta do χ²: H0 é "o desajuste é pelo menos o de um RMSEA populacional de ,08", e p
< ,05 mostra ajuste equivalente a aceitável. Aplicado ao χ² de cada passo e ao Δχ² com o
Δgl entre passos, diz se a invariância foi demonstrada, não só se não foi rejeitada.
Reproduz os três exemplos do livro (configural p = ,147; métrico p = ,076; diferença p =
,1375, nenhum demonstrado). Com WLSMV ou MLR, usar o χ² escalonado e saber que a
distribuição não central é aproximação.

**Variável contínua ou ordinal** (`invariancia_continua(fit, ordenar_por, parametros)`):
testes baseados no escore (Merkle & Zeileis, 2013; Merkle, Fan & Zeileis, 2014; Wang et
al., 2014, que Petersen cita) num modelo lavaan de um grupo por ML. Ordena os casos pelo
moderador (idade, anos de escolaridade) e testa se as contribuições ao escore das cargas
(ou interceptos, ou de todos) flutuam além do acaso, sem cortar o moderador em faixas.
Funcional DM para contínuo, WDMo para ordinal, LMuo para categórico. Nos evals, igual ao
`strucchange::sctest` direto; carga que cresce com a idade dá p < ,001 e a invariante,
p > ,05. Moderador com ausente para a função (ajustar o modelo nos casos com o moderador).
p ≥ ,05 não prova invariância; p < ,05 não diz onde muda (gráfico do processo de
flutuação, `plot(strucchange::gefp(fit, order.by = ...))`, no computador).

## Invariância longitudinal

`invariancia_longitudinal(dados, ondas, construto)` (v1.5, `references/mee.md`, seção
6) faz a sequência entre ondas com `longFacNames` e `long.equal` do `measEq.syntax()`,
resíduos do mesmo item correlacionados entre ondas, os mesmos passos e critérios daqui.
Nos evals, dados invariantes passam em todos os passos e um intercepto deslocado na onda
2 derruba o passo escalar. É o passo prévio da curva de crescimento com variável latente
por onda (aula "Latent growth curve com SEM").

## Comparação entre grupos sem invariância escalar

Comparar médias latentes entre grupos pressupõe ter passado pelos passos
configural, métrico e escalar (Putnick & Bornstein, 2016, seção sobre
diferenças de médias latentes); no caso ordinal, limiares e interceptos
entram nessa exigência. Sem isso, diferença de média não vai para o
manuscrito como diferença no construto, e o mesmo vale para escore somado.

## Viés de medida, DIF e viés preditivo

Três perguntas diferentes sobre equidade (Petersen, 2024, cap. 15; Cook & Pitoniak,
2026, cap. 6, na pasta Psicometria-R):

- **Invariância de medida** (esta página): os parâmetros do modelo de medida (cargas,
  limiares, interceptos) são iguais entre grupos? Falha aqui impede comparar médias
  latentes.
- **DIF** (`tri_dif()`, `references/tri.md`): olhar item a item a mesma pergunta, pela
  TRI ou pela regressão logística. Não somar DIF e invariância como evidências
  independentes.
- **Viés preditivo** (`vies_preditivo()`, `references/associacao.md`, seção 7b): o
  escore prevê um critério externo do mesmo jeito nos grupos (Cleary, 1968)? Um teste
  pode ser invariante e ainda prever de forma diferente, e o contrário.

MIMIC (o grupo como covariável que prediz o fator e, item a item, o próprio item) é a
alternativa quando os grupos são pequenos para a análise multigrupo ou a variável é
contínua (idade). A skill não tem função própria; o MIMIC de aquiescência de
`controle_aquiescencia()` mostra a sintaxe (`references/vieses-resposta.md`).
