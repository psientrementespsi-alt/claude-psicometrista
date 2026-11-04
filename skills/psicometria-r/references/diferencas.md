# Diferenças: comparação de grupos e condições (est_diferencas.R)

## Sumário

- 13 (linha 31)
- 35 (linha 53)
- 65 (linha 83)
- 112 (linha 130)
- 135 (linha 153)
- 162 (linha 180)
- 184 (linha 202)
- 197 (linha 215)
- 216 (linha 234)
- 232 (linha 250)
- 242 (linha 260)
- 249 (linha 267)
- 260 (linha 278)
- 314 (linha 332)
- 322 (linha 340)

Módulo de testes de diferença da skill. Padrões seguem o livro do usuário (Biblioteca
Análise de Dados, Vol. 4 e 5, JASP 0.95.1); a literatura entra como argumento ou ao lado.
Testes: `scripts/evals_diferencas.R` (197 testes; 82 com os bancos do livro, todas as âncoras
reproduzidas; sem `PR_DADOS_JASP` rodam 115, dos quais 6 conferem contra o TOSTER).

Índice: 1. Mapa de decisão · 2. Saída comum · 3. Testes t · 4. Mann-Whitney e Wilcoxon ·
5. ANOVA de um fator · 6. Medidas repetidas, fatorial intra e mista · 7. ANOVA fatorial ·
8. ANCOVA · 9. MANOVA · 10. Kruskal-Wallis e Friedman · 11. Rótulos de magnitude ·
12. Equivalência (TOST) e não significativo · 13. Armadilhas · 14. Fontes e divergências

## 1. Mapa de decisão

| Desenho | Função |
|---|---|
| 2 grupos de pessoas distintas, VD contínua | `teste_t_independente()` |
| Mesmas pessoas, 2 momentos | `teste_t_pareado()` |
| 1 grupo contra valor de referência | `teste_t_uma_amostra()` (Wilcoxon de uma amostra ao lado) |
| 2 grupos, normalidade gravemente violada com n pequeno | `teste_mann_whitney()` |
| 2 momentos, idem | `teste_wilcoxon_pareado()` |
| 3+ grupos, 1 fator | `anova_um_fator()` (nunca bateria de t:) |
| 3+ condições nas mesmas pessoas | `anova_medidas_repetidas(medidas)` |
| 2+ fatores intra | `anova_medidas_repetidas(medidas, intra = list(...))` |
| Entre × intra (grupo × tempo) | `anova_medidas_repetidas(medidas, entre = "grupo")` |
| 2+ fatores entre grupos | `anova_fatorial()` |
| Grupos ajustando por covariável contínua | `ancova()` |
| 2+ VDs com coerência teórica | `manova_grupos()` |
| 3+ grupos / 3+ condições sem normalidade | `teste_kruskal_wallis()` / `teste_friedman()` |
| Mostrar que a diferença (ou a correlação) é pequena demais para importar | `teste_equivalencia()` (seção 12) |

Não paramétrico só com violação grave e n pequeno; antes, Welch,
correção de gl ou bootstrap (curso, "Testes Não-Paramétricos").

## 2. Saída comum

Todas devolvem `pr_resultado`: `$tabela`, `$efeito` (coluna `principal` = o índice que o
livro manda relatar; `termo` quando há vários efeitos), `$pressupostos` (checagem,
estatística, gl, p, leitura por magnitude), `$posthoc`, `$decisoes`, `$avisos`, `$apa`.
Extras ficam no primeiro nível do objeto: `$descritivas` (n, M, DP por grupo, momento ou
célula), `$medias` e `$inclinacoes` (ANCOVA), `$univariadas` e `$box_m` (MANOVA),
`$efeitos_simples_F` (fatorial), `$modelo`. Todas aceitam `decimal = ","` (itens passam a
ser separados por "; ") e `idioma = "en"`; `alfa = .05` e `nivel = .95`. Casos incompletos
nas variáveis da análise saem por listwise e são contados em `$n_excluidos`.
Convenção das comparações par a par: rótulo "A - B" = A menos B, na ordem dos níveis.
Privacidade: o print da base omite as estatísticas de linhas de `$descritivas` com n < 5
(grupos, momentos e células de desenhos mistos), e a frase segue a mesma regra: grupo ou
momento com n < 5 sai da frase sem M, DP ou Md ("n < 5, descritivas omitidas"), com registro
em `$avisos` (t independente, pareado, uma amostra, Mann-Whitney, Wilcoxon pareado; na ANCOVA,
a média ajustada). No t independente saem M e DP dos dois grupos: com M e DP do outro grupo,
a diferença e o t devolvem por álgebra M e DP do grupo pequeno (com n = 2, os dois escores).
O que a regra não fecha: as estatísticas do teste no print (diferença, gl de Welch, Δ de
Glass) com as descritivas do outro grupo ainda permitem essa conta; no t de uma amostra, a
diferença M − mu devolve M; na ANCOVA, as diferenças ajustadas devolvem a média ajustada.
Com grupo tão pequeno, juntar categorias ou não relatar o teste. t, ANOVA, Mann-Whitney e
Kruskal-Wallis exigem ao menos 2 casos por grupo. VD sem variação para com mensagem que diz
o motivo, no lugar do erro cru do R: VD toda igual (testes t, ANOVA de um fator, fatorial,
ANCOVA, Mann-Whitney e Kruskal-Wallis), igual dentro de todos os grupos (t independente,
ANOVA de um fator, fatorial e Mann-Whitney; o Kruskal-Wallis roda, é separação perfeita) e
diferenças todas iguais (pareado). Medidas repetidas, MANOVA e Friedman não têm essa
checagem. Os IC da frase usam o rótulo de `nivel`; com `decimal = ","`, os itens de
`$decisoes` e da frase são separados por "; " (inclusive Levene, Mauchly e inclinações da
ANCOVA).

## 3. Testes t

`teste_t_independente(dados, vd, grupo, rota = "classica", centro_levene = "media", indice = "d", alternativa = "bilateral", controle = NULL)`

- Student e Welch sempre calculados. `rota = "classica"` (livro): Levene p < alfa
  leva à linha de Welch; senão Student. `rota = "welch"`: Welch sempre (Delacre et al.,
  2017; Zimmerman, 2004; Caixa 13.3). A outra linha fica com `escolhido = FALSE`.
- Levene: o livro chama de Levene e imprime F = 2,362 no Directed Reading, que é o Levene
  centrado na **média**. O painel de teste t do JASP 0.95.1 vem por padrão com a variante
  centrada na mediana, rotulada "Brown-Forsythe" (2,342 no mesmo banco). Padrão
  `centro_levene = "media"` (o número do livro); a outra versão sai em `$decisoes`.
- Normalidade: Shapiro-Wilk por grupo como um indício com leitura de assimetria
  e curtose (Vol. 1); Shapiro dos resíduos, que é o que o painel do JASP mostra.
- Efeito: d com DP combinado (Student) ou com a média simples das variâncias (Welch), como
  o JASP; IC pelo t não central. g de Hedges (J gama exata) e Δ de Glass (`controle`) ao
  lado. `indice = "g"` torna o g principal; com gl < 20 um aviso lembra que o d superestima
  o efeito em 4% ou mais (a correção depende só dos gl; ver divergência 13).
- Sinal: diferença = 1º nível − 2º (a ordem do JASP); a frase bilateral usa |t| e diz qual
  grupo teve média maior. `alternativa = "maior"` testa 1º nível > 2º. Com teste unilateral a
  frase declara H1 ("unilateral (H1: média de Control maior que a de Treat)"), dá t, diferença
  e d com o sinal de 1º − 2º (os de `$tabela` e `$efeito`), IC unilaterais ("IC 95%
  unilateral [-17.34, ∞]") e marca "no sentido oposto ao de H1" quando os dados vão contra a
  hipótese (no Directed Reading com `"maior"`: t(42) = −2.27, p = .986). O pareado faz o mesmo
  com depois − antes (H1: aumento ou diminuição) e o de uma amostra com M − mu.
- Levene degenerado: com 2 casos num grupo (os dois desvios ao centro são iguais), variância
  zero num grupo ou ajuste perfeito (F que é ruído numérico, como F(1, 2) = 8,1e31 com 2
  casos por grupo), o Levene não decide a rota: vale a linha de Welch, o motivo vai para
  `$decisoes`, `$avisos` e a leitura do pressuposto, e F e p saem como NA quando são ruído.
  O mesmo vale na ANOVA de um fator (Welch e Games-Howell).
- Frase: `Um teste t de Student para amostras independentes indicou média maior em Treat
  (M = 51.48, DP = 11.01) do que em Control (M = 41.52, DP = 17.15), t(42) = 2.27, p = .029;
  diferença média 9.95, IC 95% [1.09, 18.82], d = 0.68, IC 95% [0.07, 1.29] (efeito moderado).`

`teste_t_pareado(dados, antes, depois)`: diferença = depois − antes (positivo = aumento,
); o JASP mostra o t com sinal oposto. Normalidade da diferença. Quatro d:

| Índice | Fórmula | Papel |
|---|---|---|
| d na métrica original (principal) | M_dif / √((DP₁² + DP₂²)/2) | o do livro: 0,47 no Weight Gain (Caixa 13.2) |
| d_av de Lakens | M_dif / ((DP₁ + DP₂)/2) | quase igual; Lakens (2013) |
| d_rm | d_z · √(2(1 − r)) | o "Correct for correlation" do JASP: 0,45, não 0,47 |
| d_z | M_dif / DP_dif | padrão do JASP (2,71); infla com r alto; serve para poder |

IC do principal pelo método "av" do effectsize (leva o pareamento em conta).
`teste_t_uma_amostra(dados, vd, mu)`: d = (M − mu)/DP; o livro trata só conceitualmente
; o Wilcoxon de uma amostra sai na tabela.

## 4. Mann-Whitney e Wilcoxon

- `teste_mann_whitney(dados, vd, grupo, sinal_rrb = "direcional", regua_rrb = "r", exato = NULL)`:
  U = U do 1º grupo (o W do R e do JASP), z com correção de empates, Hodges-Lehmann com IC,
  descritivas com mediana e média dos postos (não é teste de mediana).
- `teste_wilcoxon_pareado(dados, antes, depois, ...)`: zeros descartados antes de ranquear
  (método de Wilcoxon), V = soma dos postos positivos de depois − antes, z sem correção de
  empates (como o JASP), r_rb pareada = 2V/maxw − 1 (Kerby, 2014): pondera pelos postos,
  não é proporção de pessoas que mudaram. O Wilcoxon de uma amostra (em
  `teste_t_uma_amostra`) segue a mesma regra, com os valores iguais a mu descartados.
- Unilateral: a frase declara H1 ("valores de A tendem a ser maiores que os de B"; no
  Wilcoxon, "aumento" ou "diminuição"), mantém z e r_rb com sinal e marca "no sentido oposto
  ao de H1" quando é o caso.
- Sinal do r_rb: `"direcional"` (padrão): positivo = 1º grupo maior (MW) ou aumento
  (Wilcoxon), o sinal do effectsize e do JASP ≥ 0.96. O JASP 0.95.1 usa 1 − 2U/(n₁n₂) no
  independente (sinal invertido, corrigido no próprio JASP em 2026) e antes − depois no
  pareado; `sinal_rrb = "jasp_0951"` reproduz a tela antiga. Magnitude idêntica.
- p exato: `exato = NULL` reproduz o R 4.5.1 do JASP 0.95.1 (Mann-Whitney exato só com n < 50
  por grupo sem empates; Wilcoxon só com n < 50 diferenças não nulas, sem empates e sem
  zeros; senão normal com correção de continuidade e de empates). O R 4.6 mudou o
  `wilcox.test` (exato condicional com empates; no pareado, zeros nos postos e V muda);
  `exato = TRUE` chama esse `wilcox.test(exact = TRUE)` do R instalado; `FALSE` força a normal.

## 5. ANOVA de um fator

`anova_um_fator(dados, vd, grupo, rota = "classica", posthoc = "auto", exigir_omnibus = TRUE, centro_levene = "media")`

- Tabela com três linhas: ANOVA clássica (SQ Tipo III, contrastes de soma zero, SQ do
  fator e do resíduo), Welch (`oneway.test`) e Brown-Forsythe F* (fórmula do
  `onewaytests::bf.test`, implementada sem dependência). Rota clássica: Levene
  (centrado na média, o do painel de ANOVA do JASP) sem evidência → clássica + Tukey;
  com evidência → Welch + Games-Howell. `rota = "welch"`: Delacre et al. (2019).
- Efeito: ω² principal (Caixa 14.1), η² (= η²p com um fator) ao lado; IC pelo F
  não central bilateral (o padrão do effectsize é unilateral). ω² vem da tabela clássica
  também na rota Welch, como no JASP. Na rota Welch o ω² clássico usa o QM do erro combinado,
  que supõe variâncias iguais, e pode contradizer o F de Welch (F de Welch significativo com
  ω² = ,000). Por isso, nessa rota: terceira linha em `$efeito` com o ω² aproximado pelo F de
  Welch, (F − 1)gl1/(F·gl1 + gl2 + 1) com os gl de Welch (`effectsize::F_to_omega2`; é o
  que `effectsize::omega_squared()` devolve para `oneway.test(var.equal = FALSE)`, marcado
  pelo próprio pacote como aproximação); a frase declara a origem de cada um ("ω² da ANOVA
  clássica, como no JASP: ...; pelo F de Welch, ω² aproximado = .157"); nota em `$decisoes`;
  aviso sempre, e aviso de contradição quando o Welch é significativo e o IC do ω² clássico
  inclui zero.
- Pós-hoc só com ômnibus significativo: `"auto"` = Tukey ou Games-Howell pela
  rota; ou `"bonferroni"`, `"holm"`, `"nenhum"`. Colunas: dif, EP, t, gl, p e IC ajustados,
  `d_qme` (dif/√QME, a coluna "Cohen's d" do JASP, a da frase) e `g_par` (g de Hedges só com
  os dois grupos do par, o índice que o corpo cita). Games-Howell manual
  (fórmula do JASP) = `rstatix::games_howell_test`. Ômnibus significativo sem par isolado
  não é contradição; a frase registra o menor p ajustado.

## 6. Medidas repetidas, fatorial intra e mista

`anova_medidas_repetidas(dados, medidas, intra = NULL, entre = NULL, esfericidade = "classica", limiar_epsilon = .75, posthoc = "bonferroni", exigir_omnibus = TRUE)`

- Formato largo, uma coluna por célula intra. `intra = NULL` cria o fator "condicao";
  `intra = list(bloco = c("A","B"), condicao = c("C","I","N"))` segue a ordem das células do
  JASP (último fator varia mais rápido). `entre` = fator(es) entre sujeitos. Motor:
  `afex::aov_car`, SQ Tipo III.
- Esfericidade (fatores com 3+ níveis): Mauchly W, χ², gl, p; ε de GG, HF (truncado em 1,
  como o JASP) e limite inferior. Rota do livro: o Mauchly decide se corrige; GG
  quando o ε de GG arredondado a duas casas é ≤ ,75, HF acima (fronteira inclusiva: o livro
  adota GG com ε = ,752). `esfericidade = "epsilon"` decide só pelo ε;
  `limiar_epsilon = .60` segue Blanca et al. (2023). Tabela com gl e p sem e com correção.
- Efeito: η²p principal com qualquer fator intra; η²G (afex) e ω² ao lado (com um
  só fator intra, a fórmula generalizada do JASP: 0,020 no Stroop; nos demais, ω²p do
  effectsize).
- Pós-hoc Bonferroni (Caixa 14.3) no modelo multivariado do afex: erro por par, gl = n − 1
  (o t(18) de). Com um fator: `d_par` (DP combinado do par: o ,35 do corpo),
  `d_jasp` (DP médio das condições: o −0,366 da tela) e `d_z` (≈ 1,06, inflado). Interação
  de 2 fatores significativa: efeitos simples com Bonferroni (ou Holm) sobre o conjunto.
- Mista: Levene por medida entre os grupos (média, como o JASP) e M de Box (critério p < ,001).

## 7. ANOVA fatorial entre grupos

`anova_fatorial(dados, vd, fatores, tipo_sq = 3, posthoc = "tukey", centro_levene = "media")`

- SQ Tipo III com soma zero (pré-cálculo do livro e padrão do JASP); `tipo_sq = 2`
  sem interação relevante (Langsrud, 2003). Células vazias param a função.
- ω²p principal e η²p com IC; R²; Levene entre as células (o do JASP); Shapiro.
- A frase abre pela interação. Interação de 2 fatores significativa: efeitos
  simples com t no gl do erro, IC não ajustado, d = dif/√QME, p sem ajuste (o que o livro
  imprime) e p ajustado sobre o conjunto (Caixa 16.A): Bonferroni, ou Holm com
  `posthoc = "holm"`; `"nenhum"` não os calcula. `$efeitos_simples_F` traz os F conjuntos.
  Sem interação: Tukey nos efeitos principais com 3+ níveis. Interação de 3ª ordem: aviso.

## 8. ANCOVA

`ancova(dados, vd, fator, covariaveis, posthoc = "bonferroni", bootstrap = FALSE, n_boot = 1000, semente = 2026, tipo_ic_boot = "bca", grupos_aleatorizados = NA)`

- Um fator, 1+ covariáveis, SQ Tipo III. ω²p principal e η²p; não parciais em
  `$efeito_nao_parcial` só para conferência.
- Checagens: homogeneidade das inclinações (fator × covariável; Caixa 15.A), com a
  inclinação por grupo em `$inclinacoes` (parcial, com as demais covariáveis, quando há 2+:
  é a que o teste avalia) e r bivariado; diferença dos grupos na covariável; Levene da VD
  bruta (o do corpo) e dos resíduos (o do painel do JASP); Shapiro dos resíduos.
- Médias brutas e ajustadas (`$medias`). Comparações das ajustadas, só para o fator: p e IC
  ajustados (os da tela), IC não ajustado (o do corpo), d = dif/√QME.
- `bootstrap = TRUE`: 1000 reamostragens de casos estratificadas por grupo (Caixa 15.C,
), IC BCa (percentil se falhar); não conserta inclinações heterogêneas; o gerador
  aleatório do usuário é restaurado. Sem bootstrap e com normalidade ou homogeneidade
  violadas, um aviso sugere essa rota.
- Avisos: inclinações violadas (modelar a interação); Levene só nos resíduos; sem
  `grupos_aleatorizados = TRUE`, Miller e Chapman (2001) e "ajustando por", não "controlando".

## 9. MANOVA

`manova_grupos(dados, vds, fatores, tipo_sq = 3, alfa_box = .001)`

- Quatro estatísticas por termo; Pillai é a principal, fixada antes (Olson, 1976). F pelas aproximações do `summary.manova` (as da tela do JASP; o Hotelling-Lawley
  do corpo, 3,19 com gl 50,19, é a de McKeon), com p = nº de VDs mesmo quando o efeito tem
  menos gl; conferido com 2 grupos × 3 VDs e 3 grupos × 4 VDs.
- Efeito multivariado η²p = V/s, s = mín(nº de VDs, gl do efeito): o "≈ ,17".
- M de Box (`heplots::boxM`, = JASP) com um só critério, `alfa_box = .001` (usual pela
  sensibilidade do teste; o livro não declara critério), na leitura e na frase.
  Normalidade multivariada: Henze-Zirkler dos resíduos (o objeto do pressuposto; = `MVN::hz`;
  0,920, p = ,023 no Viagra) e, para conferir a tela, o Shapiro-Wilk multivariado do JASP
  sobre as VDs brutas (W = ,955, p = ,200). Box violado: manter Pillai e registrar.
- Univariadas com η²p em `$univariadas`; a frase só as relata com multivariado
  significativo. `tipo_sq = 1` reproduz o JASP (sequencial); com um fator coincidem.

## 10. Kruskal-Wallis e Friedman

- `teste_kruskal_wallis(dados, vd, grupo, correcao = "holm", n_boot = 1000, semente = 2026)`:
  H, gl, p (2+ casos por grupo); ε²_R = H/(N − 1) principal com IC bootstrap percentil;
  η²_H (o "rank η²" do JASP e do curso) ao lado. Dunn (o do JASP; = `rstatix::dunn_test`)
  só com ômnibus significativo, p de Bonferroni e de Holm; r_rb do Mann-Whitney do par.
- `teste_friedman(dados, medidas, nomes, correcao = "holm")`: χ², gl, p, W de Kendall =
  χ²/(n(k − 1)) principal; Conover (o do JASP; = `PMCMRplus::frdAllPairsConoverTest`) com
  Bonferroni e Holm; r_rb do Wilcoxon pareado ("A - B": positivo = A maior).

## 10b. Delineamentos complexos sem normalidade (v1.5)

Dois ou mais fatores, medidas repetidas cruzadas com grupos ou várias VDs, com dado
ordinal ou assimétrico: `np_fatorial()`, `np_longitudinal()` e `manova_np()`, em
`references/naoparametricos-complexos.md` (rankFD, nparLD e MANOVA.RM, a rota do curso
"Testes não paramétricos para delineamentos complexos").

## 11. Rótulos de magnitude

| Índice | Cortes (pequeno, médio, grande) | Fonte |
|---|---|---|
| d, g | ,20 / ,50 / ,80 como limiares inferiores; "próximo do" a menos de ,05 do seguinte | Cohen (1988); livro |
| η², η²p, ω², ω²p, ε² | ,01 / ,06 / ,14 | derivados de f de Cohen; aplicar a ω² é extrapolação usual |
| r_rb, W | ,10 / ,30 / ,50 (`regua_rrb = "grissom_kim"`: ,125/,304/,465) | convenção de r (JASP, effectsize) |

Os cortes reproduzem as leituras do livro: 0,47 "pequeno próximo do moderado",
ω² = ,025 pequeno, ω²p = ,17 grande, ,008 desprezível.

## 12. Equivalência (TOST) e não significativo

**Não significativo não quer dizer igual.** p ≥ α só diz que os dados não destoam da
hipótese nula; também não destoam dos outros valores do IC, e o teste pode ter tido pouco
poder (Greenland et al., 2016, equívocos 4, 6 e 8: aceitar H0, dizer que não houve efeito,
ler o efeito como pequeno). Toda função deste módulo, quando o teste principal não é
significativo, põe em `$avisos` uma nota que começa por "Não significativo não quer dizer
igual", com a estimativa e o IC no teste bilateral (no unilateral o IC é aberto e fica de
fora). Com vários termos ou pares, a nota lista só os não significativos. A nota aponta
`teste_equivalencia()` quando a skill tem equivalência para a estatística (t independente,
pareado e de uma amostra); nos demais (Mann-Whitney, Wilcoxon, ANOVAs, ANCOVA, MANOVA,
Kruskal-Wallis, Friedman) manda relatar o IC do efeito e não escrever que não há diferença.
O texto não muda a frase APA: a frase continua dizendo o que o teste disse.

**`teste_equivalencia(dados, tipo, limites, ...)`**: dois testes unilaterais (TOST, como em
Lakens, 2017), H0 inferior Δ ≤ ΔL e H0 superior Δ ≥ ΔU. Há equivalência quando os dois
rejeitam, o mesmo que o IC de 1 − 2α (90% com α = .05) ficar dentro dos limites; relata-se
o maior dos dois p (Lakens, 2017, nota 2 e exemplo). O teste usual (Δ = 0) sai ao lado, e os
dois juntos dão um dos quatro resultados de Lakens (2017, Figura 1), em `$resultado`:
equivalente e não diferente de zero; diferente de zero e equivalente; diferente de zero e
não equivalente; indeterminado.

| `tipo` | Estimativa | Conta (Lakens, 2017) |
|---|---|---|
| `"independente"` (`vd`, `grupo`) | diferença 1º nível − 2º nível, como `teste_t_independente()` | Welch com gl de Satterthwaite (eq. 3 e 4); `welch = FALSE` dá Student |
| `"pareado"` (`antes`, `depois`) | depois − antes, como `teste_t_pareado()` | eq. 6 |
| `"uma_amostra"` (`vd`, `mu`) | M − mu | eq. 7 |
| `"correlacao"` (`x`, `y`) | r de Pearson | z de Fisher com EP = 1/√(N − 3), contra a normal (eq. 8) |

Saída: `$tabela` (os dois unilaterais e o usual: estatística, gl, p), `$efeito` (estimativa
com IC 90%, o principal, e IC 95%), `$limites` (na unidade da estimativa), `$equivalente`,
`$diferente_de_zero`, `$resultado`, `$p_tost`, mais `$decisoes`, `$avisos` e `$apa` (PT ou
EN, `decimal = ","`).

- **Limites**: `c(inferior, superior)` ou um número (±). Vêm do menor efeito de interesse,
  fixado antes de ver os dados; a função não sugere valor e registra isso em `$decisoes`.
- **Uma amostra**: os limites são relativos a mu, como na eq. 7 (`limites = c(-2, 2)` com
  `mu = 20`). O `TOSTER::t_TOST` 0.8.6 pede os limites na escala da média (`eqb = c(18, 22)`
  com `mu = 20`) e dá o IC da média, não da diferença; t e p são os mesmos (conferido nos
  evals).
- **`escala = "d"`**: limites em d de Cohen (d_z no pareado) são convertidos pelo DP da
  amostra (média das variâncias com Welch, DP combinado com Student, DP das diferenças no
  pareado). O limite passa a depender dos dados, e a documentação do `TOSTER::t_TOST`
  avisa que limites padronizados dão resultado enviesado ("Raw is strongly recommended as
  SMD bounds will produce biased results"); a função registra o aviso. Padrão: `"bruta"`.
- **Correlação**: limites em r, entre −1 e 1; o usual ao lado é o t do `cor.test`, o mesmo
  de `correlacoes()` (`references/associacao.md`).
- **IC mais largo que os limites**: a amostra não tem como mostrar equivalência, por maior
  que seja a concordância com zero; sai aviso em `$avisos`, e no resultado indeterminado a
  frase não fala em igualdade.
- Conferência: t, gl, p e IC batem com `TOSTER::t_TOST` (independente com Welch e Student,
  pareado, uma amostra, limites em d), e p, z e IC batem com `TOSTER::z_cor_test`
  (correlação), tolerância 1e-8.

## 13. Armadilhas

- `carregar.R` com `source()` direto; fator numérico (dose 1/2/3) vira fator (rotular antes).
- p > ,05 num teste de pressuposto não "confirma" o pressuposto.
- Pós-hoc com ômnibus não significativo é pescar resultado; `exigir_omnibus =
  FALSE` só para comparações planejadas e declaradas.
- Medidas repetidas são listwise; com ausentes, modelos mistos ou GEE.

## 14. Fontes e divergências

Livro (Vol. 4-5): Student/Welch; JASP; d, g, Δ
; d pareado; ANOVA; esfericidade (,
); não paramétricos; ω² vs η²p; pós-hoc
; ANCOVA;
fatorial; MANOVA; recomputação.

Divergências decididas (detalhes nas seções 3 a 10):

1. Levene do t: livro = média (2,362); painel t do JASP = mediana (2,342). Padrão: média.
2. d pareado: livro 0,47 = DP combinado; "Correct for correlation" do JASP = d_rm (0,45).
3. Pós-hoc de um fator: tela d = −0,449 (√QME); corpo g = ,46/,39 = g do par. Saem os dois.
4. Pós-hoc da RM: tela −0,366 (DP médio das condições); corpo ,35 = `d_par`, o da frase.
5. Métrica com mais de um fator: a Caixa 14.1 anuncia η²p no Vol. 5, que adota
   ω²p na ANCOVA e na fatorial. Vale o Vol. 5 entre grupos; com fator
   intra, η²p. O curso prefere ω²p na mista; prevalece o livro.
6. Dunn e Conover: curso = Bonferroni até 3 grupos; Holm (1979) controla o mesmo erro com
   poder igual ou maior. Padrão: Holm; as duas colunas saem.
7. W de Kendall: o curso recomenda não interpretá-lo; a extração do curso registra que é o
   efeito usual do Friedman (concordância na ordenação). Mantido.
8. Cortes de r_rb atribuídos pelo curso a Grissom & Kim (2012) não foram achados em fonte
   publicada: padrão = régua de r.
9. Hotelling-Lawley: reproduz-se a tela (3,101; gl 50), não o corpo (3,19; gl 50,19).
10. Henze-Zirkler = 14,459 não se reproduz como número: pela fórmula do MVN, as
    VDs brutas dão 0,69 (p = ,108) e os resíduos 0,92 (p = ,023). Nos resíduos, o objeto do
    pressuposto, o HZ rejeita, como conclui o livro; é o que o módulo relata.
11. `wilcox.test` do R 4.6 ≠ R 4.5.1 do JASP 0.95.1: regra antiga fixada (seção 4).
12. Brown-Forsythe F*: ausente no livro; sai na tabela por ser barato (Delacre et al., 2019).
13. Erro do livro não reproduzido: e dizem que o g de Hedges é preferível
    com grupos de tamanhos bem diferentes. A correção J depende só dos gl (n₁ + n₂ − 2), não
    do desbalanceamento. O módulo não emite aviso por desbalanceamento nem troca o índice do
    pós-hoc por isso; o g fica como opção para amostras pequenas (aviso com gl < 20).

Curso (Psicometria Online Academy): páginas de teste t (independente e pareado), ANOVA
one-way, ANCOVA, medidas repetidas, fatorial de medidas repetidas, MANOVA, Mann-Whitney,
Kruskal-Wallis, Wilcoxon e Friedman.

Literatura: Delacre, Lakens & Leys (2017; correção 2022); Delacre, Leys, Mora & Lakens
(2019); Zimmerman (2004); Lakens (2013); Lakens (2017); Greenland et al. (2016);
Goertzen & Cribbie (2010, pela eq. 8 de Lakens, 2017); Dunlap et al. (1996); Greenhouse &
Geisser (1959); Huynh & Feldt (1976); Blanca et al. (2023); Games & Howell (1976); Holm
(1979); Miller & Chapman (2001); Langsrud (2003); Olson (1976); Cohen (1988). DOIs em
`referencias-verificadas.md`.
