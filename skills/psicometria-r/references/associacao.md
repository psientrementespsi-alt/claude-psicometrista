# Associação: correlação, regressão e desfecho categórico (est_associacao.R)

Módulo de análises de associação da skill. Padrões seguem o livro do usuário (Biblioteca
Análise de Dados, Vol. 1-3, JASP 0.95.1); a literatura entra como argumento ou ao lado.
Testes: `scripts/evals_associacao.R` (134 testes; 54 usam os bancos do JASP e reproduzem o livro;
sem `PR_DADOS_JASP` rodam 80).

Índice: 1. Mapa de decisão · 2. Saída comum · 3. Correlações · 4. Regressão linear ·
5. Hierárquica · 6. Logística · 7. ROC · 8. Qui-quadrado · 9. McNemar e Cochran ·
10. Armadilhas · 11. Fontes e divergências

## 1. Mapa de decisão

| Pergunta | Função |
|---|---|
| Duas ou mais contínuas variam juntas? | `correlacoes()` (Pearson e Spearman lado a lado) |
| Contínua × dicotômica verdadeira | `correlacao_ponto_bisserial()` |
| X e Y controlando Z | `correlacao_parcial()` |
| "X se associa mais a Y do que a Z?" | `comparar_correlacoes()` (independentes ou dependentes) |
| Quanto o desfecho contínuo acompanha preditores | `regressao_linear()` |
| Contribuição "acima e além" de blocos | `regressao_hierarquica()` |
| Desfecho dicotômico | `regressao_logistica()` (nunca a linear) |
| Escore discrimina um desfecho? Ponto de corte | `curva_roc()` |
| Duas categóricas independentes | `qui_quadrado()`; uma contra distribuição esperada: `qui_quadrado_aderencia()` |
| Binário nos mesmos sujeitos, 2 momentos / 3+ | `mcnemar()` / `q_cochran()` |

Relação não monotônica (U invertido): nenhum coeficiente descreve; olhar a dispersão
. Não categorizar contínua para caber no qui-quadrado.

## 2. Saída comum

Todas devolvem `pr_resultado` (subclasse `pr_associacao`): `$tabela`, `$efeito` (coluna
`principal` = o que o livro manda relatar), `$pressupostos`, `$posthoc`, `$decisoes`,
`$avisos`, `$apa`, mais `$tabelas_extra` (modelo global, matriz de confusão, células,
chances) impressas entre o resultado e as decisões. Objetos de modelo em `$modelo`,
`$modelos`, `$roc`; índices de casos atípicos em `$extra_diagnostico` (nunca impressos).
Todas aceitam `decimal = ","` (separador de itens vira "; ") e `idioma = "en"`.
Texto vazio ("") vira ausente antes da deleção listwise; nível de fator sem casos fica.

Formato: `casas = 3` para B, β, OR, r, ρ, τ, V, φ, w e AUC, como os parágrafos-modelo
dos Vol. 2-3 (, r(101) = -,441B = -0,485, β = -0,321OR = 0,275V = 0,077); o Apêndice A do Vol. 1 usa duas,
inconsistência do livro. Testes com 2 casas, p com 3. β com zero (pode passar de 1,
); r, ρ, τ, R², V, φ, AUC e p sem. Valor que arredondaria para zero ganha casas
até 2 algarismos significativos (B = -0.0000089, não 0.000). Abaixo de 1e-10 em módulo o
valor é zero numérico (resto de ponto flutuante) e sai como zero: r = .000, não
r = .000000000000000014.

Evento (código 1) quando `evento` não é informado, regra única do módulo: 0/1 → 1;
lógico → TRUE; outro numérico → maior valor; fator → segundo nível; texto → erro pedindo
`evento` (a ordem alfabética poria "controle" como evento em "caso"/"controle"). Exceções
declaradas: ponto-bisserial (`codigo_1`) e OR do 2×2 usam a ordem alfabética do texto.

## 3. Correlações

`correlacoes(dados, vars, metodos = c("pearson", "spearman"), ausentes = "pairwise", nivel = .95, regua = "cohen", limiar_divergencia = .10, casas = 3)`

- Pearson e Spearman sempre lado a lado por padrão; `"kendall"`
  opcional (amostras pequenas ou muitos empates; aviso automático se algum par
  tem n < 30). Nenhuma troca silenciosa: a frase relata o Pearson e diz se cada outro
  coeficiente convergiu. Spearman diverge com |r − ρ| ≥ .10 ou significância diferente;
  Kendall só pela significância (τ é menor que r por construção). Divergência gera aviso:
  voltar à dispersão, não escolher o mais conveniente.
- IC: Pearson pelo `cor.test` (z de Fisher); Spearman z de Fisher com EP = 1/√(n − 3);
  Kendall assintótico de Hollander-Wolfe-Chicken. Os três iguais ao JASP 0.95.1.
- Ausentes par a par (padrão do JASP), n por par na tabela; `ausentes = "listwise"`.
- Pressupostos: normalidade de cada variável como grau (W, assimetria, curtose, n).
- Magnitude: `regua = "cohen"` (.10/.30/.50, livro), `"schober"` (.10/.40/.70/.90, livro)
  ou `"gignac"` (.10/.20/.30, literatura). Réguas de partida, não leis.
- Frase: `Exam e Anxiety associaram-se negativamente, com magnitude moderada pela régua de
  Cohen (1988), r(101) = -.441, p < .001, IC 95% [-.585, -.271], r² = 19.4% de variância
  compartilhada; Spearman convergiu, ρ(101) = -.405, ...` + frase de fronteira causal.
- Ponto-bisserial: Pearson com a dicotômica 0/1 (`codigo_1` escolhe quem é 1); tabela de
  codificação e médias por grupo; t equivalente = t de Student. Dicotomia artificial pede
  o bisserial, não este (curso, "Correlação Ponto Bi Serial"). A tabela de grupos fica em
  `$tabelas_extra` com a coluna `n`: o print apaga média e DP de grupo com n < 5 (com
  n = 1 seria o escore da pessoa), a frase não traz média nem DP por grupo e `$avisos`
  registra o grupo pequeno. Limite, dito também no aviso: o t equivalente mais média e
  DP do outro grupo ainda permitem recalcular a média do grupo pequeno (com n = 1, o
  escore); não publicar a análise assim.
- Parcial: estimativa e p do `ppcor` (gl = n − 2 − k, = JASP); IC de Fisher com
  EP = 1/√(n − 3 − k) (o JASP não dá IC); bivariada ao lado e queda em r² (p.p.). Não
  distingue confundidor de mediador.
- `comparar_correlacoes(r1, n1, r2, n2)`: independentes, z de Fisher + IC de Zou (2007)
  da diferença. `comparar_correlacoes(r1, n, r2, r12 = r_entre_as_nao_comuns)`:
  dependentes sobrepostas, z de Steiger (1980) com Meng et al. (1992) ao lado (= cocor).
  Âncora do curso: −.653 vs −.719, r23 = .518, N = 390 → |z| = 2.057, p = .040.
- Não significativo não quer dizer igual: `correlacoes()` (por par), ponto-bisserial,
  parcial, `comparar_correlacoes()` e `qui_quadrado()` põem em `$avisos` a nota de
  Greenland et al. (2016) com o coeficiente e o IC (regra comum em `references/diferencas.md`,
  seção 12). A nota aponta `teste_equivalencia(tipo = "correlacao")` quando o coeficiente
  relatado é o Pearson e `tipo = "independente"` no ponto-bisserial (o r_pb equivale ao t
  entre os grupos). Spearman, Kendall, parcial, diferença entre correlações e V ou φ não
  têm TOST na skill: a nota manda relatar o IC e não escrever que não há associação.
- `teste_equivalencia(dados, "correlacao", limites, x = , y = )`: TOST do r de Pearson pelo
  z de Fisher com EP = 1/√(N − 3), eq. 8 de Lakens (2017), que a toma de Goertzen e Cribbie
  (2010). Limites em r, do menor efeito de interesse, fixados antes. p, z do limite decisivo
  e IC 90% iguais aos do `TOSTER::z_cor_test`. Com rho = 0 e limites ±.10, n = 2000 mostra
  equivalência e n = 30 não: o IC 90% fica mais largo que os limites, e a função avisa.

## 4. Regressão linear (simples e múltipla)

`regressao_linear(dados, vd, preditores, referencias = NULL, fatores = NULL, nivel = .95, tipo_residuo = "studentizado", ep_robusto = FALSE, casas = 3)`

- Enter. Desfecho dicotômico é recusado com indicação da logística.
- Preditor de texto/fator/lógico vira k − 1 indicadoras contra a referência (primeiro
  nível ou `referencias = list(var = "nível")`, ou `c(var = "nível")`); nunca código
  numérico.
- Grupo gravado como número (escolaridade 1/2/3, formato comum vindo do SPSS ou de
  planilha): declarar `fatores = "escolaridade"`, que vira indicadoras com o menor código
  como referência (ou a de `referencias`). Sem isso, preditor numérico com 3 a 7 valores
  distintos, todos inteiros, entra como contínuo (um coeficiente, 1 gl) e `$avisos` diz
  para declarar `fatores`; 0/1 e escalas (mais de 7 valores ou com decimais) não avisam.
  Item Likert como preditor também dispara o aviso, que nesse caso pode ser ignorado.
  Mesma regra e mesmo argumento em `regressao_hierarquica()` e `regressao_logistica()`;
  `referencias` para variável numérica fora de `fatores` dá erro.
- Categoria de preditor categórico com menos de 2 casos: erro claro antes do ajuste (com 1
  caso a indicadora reproduz o valor da pessoa e o resíduo dela fica indefinido); agrupar
  categorias ou retirar o caso.
- Coeficientes: B, EP, IC (t), β = B·DP(x)/DP(y) (= JASP), t, p, tolerância, VIF. β de
  fator com 2 níveis = recodificar 0/1 como escala; com 3+ níveis fica em branco.
- Modelo (`$tabelas_extra`): R, R², R² ajustado, RMSE vs RMSE nulo (DP do desfecho), F,
  gl, p, f². Efeito principal: R² ajustado e f² = R²/(1 − R²) (.02/.15/.35).
  Os dois saem sem IC, e `$decisoes` diz "R² ajustado e f² sem IC: não calculado" (na
  hierárquica, "R² ajustado, ΔR² e f²"); quem monta o relato leva a nota junto.
- Diagnósticos (`$pressupostos`), todos sem excluir casos:

| Checagem | Como | Regra |
|---|---|---|
| Multicolinearidade | `car::vif`; fator 3+ níveis: GVIF^(1/gl) | alerta VIF > 5 ou 10; tolerância < .20 ou .10 |
| Independência | DW + p bilateral `lmtest::dwtest` (= JASP) | perto de 2; só com ordem significativa das linhas |
| Normalidade dos resíduos | SW nos resíduos brutos + assimetria/curtose | grau, não veto; resíduos, não desfecho |
| Homocedasticidade | Breusch-Pagan studentizado + gráfico | p < .05: aviso e sugestão `ep_robusto = TRUE` (HC3) |
| Resíduos padronizados | `rstandard` (JASP) ou `"zresid"` (SPSS) | nenhum |z| > 3; até 1% > 2,5; cerca de 5% > 2 (curso) |
| Influência | Cook | D > 1 (e 4/n como alerta) |
| Multivariados nos preditores | Mahalanobis nas colunas do modelo | χ²(k) crítico a p < .001 |

- Frase: `O modelo foi estatisticamente significativo, F(2, 100) = 13.18, p < .001,
  explicando 19.3% da variância de Exam (R² ajustado = .193, f² = 0.264).` e uma por
  preditor: `Anxiety: B = -0.485, IC 95% [-0.863, -0.107], β = -0.321, t(100) = -2.55, p = .012.`
  Modelo não significativo gera aviso para não aprofundar coeficientes (curso).

## 5. Regressão hierárquica

`regressao_hierarquica(dados, vd, blocos = list(contexto = "Gender", comportamento = "Revise", ansiedade = "Anxiety"), referencias = NULL, fatores = NULL, ...)`

- Blocos na ordem da hipótese, fixada antes de rodar; mesma amostra em todos os
  degraus. Tabela M0 (nulo) a Mk como o JASP: R², R² aj, RMSE, ΔR², F de mudança
  (`anova(anterior, atual)`, dois modelos de cada vez), gl, p, f² do bloco.
- f² do bloco = ΔR²/(1 − R² do degrau que inclui o bloco): reproduz .189 e .063 (;
). Coeficientes só do modelo final (`$tabelas_extra`); os intermediários ficam em
  `$coeficientes_degraus`, para consulta. Diagnósticos no modelo final.
- Frase por degrau: `No degrau 3, a inclusão de Anxiety (ansiedade) acrescentou
  explicação significativa acima dos blocos anteriores, ΔR² = .050, F(1, 99) = 6.28,
  p = .014, f² = 0.063.` e o modelo final. Bloco sem nome vira "bloco 1" ("block 1" com
  `idioma = "en"`).

## 6. Regressão logística binária

`regressao_logistica(dados, vd, preditores, evento = NULL, referencias = NULL, fatores = NULL, corte = .5, ic = "wald", box_tidwell = TRUE, validar = FALSE, reamostras = 500, semente = 2026, casas = 3)`

- Evento pela regra da seção 2 (texto exige `evento = "sim"`). A decisão sai escrita, com
  as referências dos fatores ("evento é o quê, comparado com quem?"). N efetivo
  e excluídos abrem a frase.
- Tabela: B, EP, z de Wald, Wald = z² (gl 1), p, OR, IC de Wald (`confint.default`,
  = JASP), mudança percentual da chance. `ic = "perfil"` acrescenta o IC por perfil.
- Modelo vs nulo: deviances, χ²(gl = parâmetros acrescentados), AIC, BIC. Pseudo-R²
  (fórmulas do JASP): Nagelkerke e Tjur principais; McFadden e Cox-Snell ao lado.
- Classificação no corte (classe 1 se p > corte, como o JASP): matriz de confusão,
  acurácia contra a taxa base, sensibilidade, especificidade, precisão, VPN, F1, Brier e
  Brier escalado, AUC com IC de DeLong.
- Pressupostos: VIF, EPV (alerta < 10), separação (avisos do glm ou |B| > 10),
  Box-Tidwell para contínuos (x·ln x; desloca x se ≤ 0), resíduos de Pearson |z| > 3 e
  Cook > 1. Com separação, a frase troca o OR afetado por "OR não estimável por máxima
  verossimilhança; usar Firth ou métodos exatos" e não imprime pseudo-R².
- Leitura do OR: < 1 → "(1 − OR)·100% menor"; 1 a 2 → "% maior"; ≥ 2 → "vezes a da
  referência" ou "multiplicada por OR a cada unidade". Nunca "vezes menos".
- Frase: `PClass: 2nd vs 1st: OR = 0.275, IC 95% [0.165, 0.457], p < .001 (chance 72.5%
  menor que em 1st).` e a classificação: `78.6% ... ante 58.6% da categoria mais frequente`.
- Validação interna (v1.6): `validar = TRUE, reamostras = 500` faz o bootstrap de
  otimismo de Harrell et al. (1996): reajusta o modelo em cada reamostra, mede AUC, Brier,
  inclinação de calibração e intercepto de calibração geral (offset) na reamostra e no
  banco original, e desconta a média da diferença. No próprio banco a inclinação vale 1 e
  o intercepto 0 por construção, por isso só o valor corrigido informa: inclinação < ,90
  avisa sobreajuste (o encolhimento esperado dos coeficientes em dados novos; Petersen,
  2024, cap. 9). Steyerberg et al. (2001) mostram o bootstrap como a validação interna
  mais eficiente; validação externa continua sendo a prova. Nos evals, perto do
  `rms::validate` (AUC tol. ,02; inclinação tol. ,08, com reamostras diferentes).

## 7. Curva ROC

`curva_roc(dados, desfecho, preditor, evento = NULL, maior_indica_evento = TRUE, cortes = NULL, casas = 3)`

- Evento pela mesma regra da logística (seção 2). `pROC` com direção fixada (o
  `direction = "auto"` pode inverter a AUC); AUC com IC de DeLong; corte de Youden (todos,
  se empatarem; ±Inf descartados); sensibilidade e especificidade nos `cortes`
  informados (positivo = escore ≥ corte). AUC ≤ .50: aviso para conferir o evento e
  `maior_indica_evento`, e a frase do ponto de corte não é gerada.
- Régua de Hosmer et al. (2013): .70 aceitável, .80 excelente, .90 excepcional. O curso
  usa rótulos próprios para as mesmas faixas ("razoável/boa", "muito boa", "excelente").
- Corte e índice a otimizar são decisão da pergunta clínica.
- v1.4, `prevalencia` e `custo_fn_fp`. VPP e VPN dependem da prevalência: os da
  amostra só valem onde a prevalência for a mesma. Com `prevalencia` (proporção de
  casos na população onde o teste vai ser usado), VPP = Se·π / [Se·π + (1 − Sp)(1 − π)]
  e VPN = Sp(1 − π) / [Sp(1 − π) + (1 − Se)π] saem ao lado dos da amostra; com
  prevalência baixa o VPP cai mesmo com sensibilidade e especificidade altas. RV+ =
  Se/(1 − Sp) e RV− = (1 − Se)/Sp, que não dependem da prevalência, saem para todo
  corte. Com `custo_fn_fp` (quanto custa um falso negativo em relação a um falso
  positivo), o corte que maximiza Se + r·Sp, r = (1 − π)/(custo·π) (pROC,
  `best.weights`), entra na tabela ao lado do Youden: é o corte por utilidade
  esperada que Petersen (2024, cap. 9, na pasta Psicometria-R) demonstra. O custo é
  decisão clínica declarada, não estimativa.

## 7b. Viés preditivo de Cleary (v1.4)

`vies_preditivo(dados, criterio, preditor, grupo)`: o escore prevê o critério do mesmo
jeito nos grupos? (Cleary, 1968; Petersen, 2024, cap. 15). Regressão moderada em
etapas (Lautenschlager & Mendoza, 1986): ômnibus (critério ~ escore contra critério ~
escore × grupo); se rejeitar, inclinações (validade diferencial); se as inclinações
não diferem, interceptos (o mesmo escore prevê critério diferente). O escore entra
centrado na média geral. Saem F, gl, p, ΔR² e f² de cada etapa, a reta de cada grupo
na escala original e o resíduo médio de cada grupo em relação à reta comum (positivo,
a reta comum subestima o grupo). Avisos: variância residual muito diferente entre
grupos (pressuposto da regressão moderada) e poder baixo para inclinações (Aguinis
et al., 2010). Viés preditivo não é viés de medida: ver `references/invariancia.md`.

## 8. Qui-quadrado

`qui_quadrado(dados, linhas, colunas, evento = NULL, referencia = NULL, exato = "auto", correcao_continuidade = FALSE, alfa_residuos = .05, semente = 2026, simulacoes = 10000, casas = 3)`
(aceita também uma tabela ou matriz de contagens em `dados`).

- Ordem de leitura do livro: esperadas → χ² e p → efeito → resíduos.
- Esperadas: no máximo 20% < 5 e nenhuma < 1. Violada com `exato = "auto"`:
  p lido = Fisher exato (troca visível em `$decisoes`); r×c grande demais para o exato
  vira Monte Carlo com semente local (não mexe no gerador da sessão).
- χ² de Pearson sem correção (JASP); 2×2 mostra Yates e Fisher ao lado e declara qual
  foi lido.
- 2×2: φ com o sinal do JASP; OR = chance do `evento` (coluna) na outra linha ÷ chance
  na `referencia` (linha), IC de Woolf, inverso ao lado, chances e probabilidades por
  grupo, OR condicional de Fisher. Com célula zero, OR e IC saem com a correção de
  Haldane (+0,5 em todas as células), e a frase e `$decisoes` dizem que o OR é corrigido e
  quanto seria sem a correção (infinito, zero ou indeterminado); o OR condicional de
  Fisher fica ao lado. r×c: V de Cramér sem correção (JASP) e w = V·√(min −
  1) = √(χ²/N) para comparar com .10/.30/.50. V, φ e w com IC bilateral
  (`effectsize`, qui-quadrado não central), que o JASP não dá.
- 2×2, risco (v1.6): razão de riscos com IC log-normal, diferença de riscos com IC de
  Newcombe (1998, híbrido de Wilson; igual ao `DescTools::BinomDiffCI(method = "score")`)
  e NNT = 1/|DR| com o IC invertido de Altman (1998): se o IC da DR cruza o zero, o do NNT
  passa pelo infinito e sai como "x (aumento) a ∞ a y (redução)". Em delineamento
  transversal o RR é razão de prevalências e o NNT não tem leitura de tratamento.
- Célula observada < 5 gera aviso: a tabela de células traz contagens, e cruzamento de
  variáveis demográficas com células < 5 não se publica sem agrupar categorias.
- Resíduos ajustados (Haberman) por célula: limiar 1,96 (livro) e Bonferroni para o
  número de células ao lado; a frase diz quais sobrevivem ao ajuste.
- Frase: `Houve associação estatisticamente significativa entre Physical Activity e Fruit
  Consumption, χ²(4, N = 1184) = 14.15, p = .007, V de Cramér = .077, IC 95% [...]
  (w = V·√2 = .109, efeito pequeno pela régua de Cohen, 1988).` + resíduos.
- Aderência: `chisq.test(x, p)`, w = √(χ²/N) com IC; esperadas violadas → p por Monte
  Carlo. As categorias vêm dos níveis do fator ou dos nomes de `proporcoes`: categoria
  sem observação entra com contagem 0 (tirá-la mudaria a hipótese). A frase decide a
  significância com alfa = 1 − `nivel`, como o resto do módulo.

## 9. McNemar e Q de Cochran

- Evento pela regra da seção 2 (fator: 2º nível da primeira medida; texto exige `evento`).
- `mcnemar(dados, antes, depois, evento = NULL, teste = "assintotico")`: tabela de pares,
  χ² sem correção (lido por padrão), com correção, exato binomial e mid-p sempre lado a
  lado; aviso se pares discordantes < 25. Efeitos: OR pareado b/c com IC e diferença de
  proporções com IC.
- `q_cochran(dados, vars, evento = NULL, ajuste = "bonferroni")`: Q ~ χ²(k − 1);
  proporção por medida; pós-hoc McNemar par a par com p ajustado (Bonferroni, como o
  curso; `"holm"` opcional) e, ao lado, o z par a par do SPSS que o curso relata; a frase
  segue `decimal` e `idioma`.
  Conferência: Σz² = k·Q/2, identidade que os números do curso cumprem (53,2).

## 10. Armadilhas

- `anova(m1, m2, m3)` com três modelos usa o resíduo do maior: F de mudança 19,88 em vez
  de 18,89 no Exam Anxiety. `confint(glm)` no R ≥ 4.4 é por perfil e não bate com o JASP.
- Chance de um grupo não é OR: 0,200 é a chance dos homens no Titanic; o OR
  homens/mulheres é 0,100, 90% menor (erro do livro em).
- r nunca em porcentagem (r² sim); Tjur é diferença em pontos percentuais, não relativa.
- Grupo gravado como 1/2/3 entra como contínuo se não for declarado em `fatores`: com
  efeito não monotônico (grupo 2 mais alto), o coeficiente linear sai nulo (B = -0,428,
  p = ,245) e o efeito real some (com indicadoras, F(2, 146) = 18,10, p < ,001). O aviso
  em `$avisos` existe para isso; não ignorar quando a variável é grupo.
- Correlações dependentes pedem o N efetivo comum aos três pares (o exemplo do curso só
  reproduz com N = 390, não com o N = 524 da matriz).
- Seleção automática de preditores (backward, forward, stepwise, por p ou por AIC; o
  curso mostra `olsrr` e `step()` e ressalva que tirar variável só pelo p nem sempre é
  recomendado): os coeficientes que sobram ficam inflados, os EP e os p ficam pequenos
  demais e o R² otimista (Harrell et al., 1996). Escolher os preditores pela teoria; se
  houver seleção, relatá-la como exploratória e medir o otimismo com `validar = TRUE`
  sabendo que o bootstrap só é honesto se repetir a seleção dentro de cada reamostra
  (a função valida o modelo final, não o procedimento de seleção).

## 11. Fontes e divergências

Livro: correlação; regressão; hierárquica
; logística; qui-quadrado. Curso (Psicometria Online
Academy), páginas de Correlação (e Fisher r-to-z, Ponto Bi Serial, Parcial), Regressão
Linear Simples, Identificação de Problemas (regressão múltipla), Atividade prática de
Regressão Múltipla (JASP), Regressão Logística Binária, Qui-Quadrado 2x2 e 4X2, McNemar
& Q de Cochran.

Âncoras: todos os números testados dos Vol. 1-3 batem com o JASP 0.95.1 (meia casa
impressa). O SW dos resíduos do livro (.965) é o dos resíduos brutos (`rstandard` daria
.963); o 9,99 do OR vem das contagens; a célula vazia de PClass do Titanic é ausente.

Divergências decididas:

| Tema | Livro / JASP | Curso | Literatura | Decisão |
|---|---|---|---|---|
| Pearson vs Spearman | rodar os dois e comparar | Spearman se não normal | Spearman automático não é a melhor opção (Bishara & Hittner, 2012) | os dois sempre, aviso de divergência |
| Faixas de r | Cohen; Schober (réguas) | Cohen, Schober | Gignac & Szodorai (2016); Funder & Ozer (2019) | `regua`, Cohen padrão |
| IC da parcial | JASP sem IC | sem IC | Fisher com n − 3 − k | IC de Fisher acrescentado |
| Correlações dependentes | "versões adaptadas" | psychometrica.de ("Fisher") | Steiger (1980); Meng et al. (1992) | Steiger lido, Meng ao lado |
| Resíduo padronizado | JASP `rstandard` | ZRESID (SPSS) | origem das regras 3/2,5/2 não resolvida | `tipo_residuo`, JASP padrão |
| Durbin-Watson | "perto de 2" | 1,5 a 2,5 | sem faixa fixa; teste formal (Durbin & Watson, 1950, 1951) | DW + p bilateral, sem faixa |
| VIF/tolerância | 5 ou 10; .20 ou .10 | > 10; < .10 ou .20 | 10 é convenção, não critério de invalidade (O'Brien, 2007) | alerta > 5, forte > 10 |
| VIF de fator 3+ níveis | JASP mostra GVIF^(1/(2gl)) | omisso | omisso | GVIF^(1/gl), escala do VIF |
| Mahalanobis | omisso | 11/15/25 | χ² crítico p < .001 (Tabachnick & Fidell, não conferido) | χ² crítico; tabela do curso citada |
| Cook | JASP > 1 | > 1 | Cook (1977); 4/n sem fonte primária | > 1; 4/n como alerta |
| Heterocedasticidade | consequência, sem alternativa | EP robusto | omisso | Breusch-Pagan + `ep_robusto` (HC3), não padrão |
| f² do bloco | "modelo completo" ambíguo | omisso | Cohen (1988) | R² do degrau (reproduz o livro) |
| IC do OR (logística) | Wald (JASP) | omisso | perfil é o padrão do R ≥ 4.4 | Wald; `ic = "perfil"` |
| Linearidade do logito | Box-Tidwell como horizonte | ausente | omisso | Box-Tidwell calculado |
| EPV, separação | separação citada | ausentes | regra usual de 10 EPV | checados |
| Faixas de AUC | omisso | rótulos próprios | Hosmer et al. (2013) .70/.80/.90 | Hosmer et al. |
| Yates | JASP sem correção | não citado | Yates conservador (Haviland, 1990; Campbell, 2007) | sem correção, Yates ao lado |
| Fisher r×c | JASP 0.95.1 só 2×2 | sugere Fisher | Cochran (1954) | exato r×c, Monte Carlo se pesado |
| Resíduos ajustados | 1,96; Bonferroni na caixa | 1,96 | sem correção não se sustenta com muitas células (Haberman, 1973; Beasley & Schumacker, 1995) | 1,96 e Bonferroni lado a lado |
| V em tabela maior | V·√(min − 1) | Cohen geral | w = V·√(k − 1) (Cohen, 1988) | V e w rotulados |
| McNemar | fora do escopo; não no JASP 0.95.1 | com e sem correção | mid-p ou sem correção (Fagerland et al., 2013) | sem correção lido; os 4 ao lado |
| Pós-hoc do Q | omisso | z do SPSS + Bonferroni | omisso | McNemar + Bonferroni; z do SPSS ao lado |
| Hosmer-Lemeshow | ausente no JASP | ausente | omisso | não calculado |
| Casas decimais | 2 no Apêndice A; 3 nos parágrafos dos Vol. 2-3 | 2 a 3 | APA: o necessário | 3 nas frases (`casas`) |
| Zero à esquerda do β | mantido | omisso | APA mantém se pode passar de 1 | mantido (como a base, est_comum.R) |
| IC de V e φ | JASP sem IC | sem IC | efeito com IC | IC bilateral do effectsize |

Referências (DOI conferido no Crossref): Cohen (1988), 10.4324/9780203771587; Gignac &
Szodorai (2016), 10.1016/j.paid.2016.06.069; Funder & Ozer (2019), 10.1177/2515245919847202;
Bishara & Hittner (2012), 10.1037/a0028087; O'Brien (2007), 10.1007/s11135-006-9018-6;
Durbin & Watson (1950), 10.2307/2332391, e (1951), 10.2307/2332325; Cook (1977),
10.1080/00401706.1977.10489493; Hosmer, Lemeshow & Sturdivant (2013), 10.1002/9781118548387;
Cochran (1954), 10.2307/3001616; Haviland (1990), 10.1002/sim.4780090403; Campbell (2007),
10.1002/sim.2832; Haberman (1973), 10.2307/2529686; Beasley & Schumacker (1995),
10.1080/00220973.1995.9943797. Referência completa das demais (Schober et al., 2018;
Steiger, 1980; Meng et al., 1992; Zou, 2007; DeLong et al., 1988; Youden, 1950; Hollander
et al., 2015; Fagerland et al., 2013; Tabachnick & Fidell, 2019; Greenland et al., 2016;
Lakens, 2017; Goertzen & Cribbie, 2010) em `referencias-verificadas.md`.
