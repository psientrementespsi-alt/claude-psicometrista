# Vieses de resposta: aquiescência, desejabilidade social e escolha forçada

Funções em `scripts/est_vieses.R`; testes em `scripts/evals_vieses.R` e, para
`controle_desejabilidade()` e `escolha_forcada()` (v1.5), em `scripts/evals_complementos.R`. Base: o
módulo de controle de vieses do curso Psicometria Online Academy (aquiescência,
desejabilidade social, escolha forçada), com o código das aulas na pasta
Psicometria-R.

Índice: 1. Aquiescência com `controle_aquiescencia()` · 2. Qual controle usar ·
3. Desejabilidade social · 4. Escolha forçada (TRI thurstoniana) · 5. Relato

## 1. Aquiescência com `controle_aquiescencia()`

Aquiescência é a tendência a concordar com o item independente do conteúdo. Ela
infla a correlação entre itens do mesmo polo, achata a dos polos opostos e pode
fazer itens positivos e negativos parecerem fatores diferentes. Em criança e
adolescente ela é maior (Soto et al., 2008, citado no curso).

```r
fatores <- list(Abertura = c("y1", "y2", "y3", "y4"), Consc = c("y5", "y6", "y7", "y8"))
polo <- c(y1 = 1, y2 = 0, y3 = 1, y4 = 0, y5 = 1, y6 = 0, y7 = 1, y8 = 0)   # 1 positivo, 0 negativo
r <- controle_aquiescencia(dados, fatores, polo, minimo = 1, maximo = 5)
r$tabela                     # ajuste: sem controle, interceptos randômicos, MIMIC
r$variancia_aquiescencia     # variância do fator de aquiescência e p
r$dados_recentralizados      # para a AFE seguinte; nunca imprimir
```

O que a função faz, na ordem do curso:

1. **Índice de aquiescência**: média das respostas brutas (sem inverter os
   negativos) de cada pessoa, centrada no ponto médio da escala (0 = sem
   aquiescência). Numa escala balanceada o traço se cancela e sobra a tendência a
   concordar. `pares` restringe o índice aos itens com par semântico.
2. **Recentralização** (ipsatização): resposta centrada menos o índice da pessoa.
   O curso arredonda e trunca no intervalo da escala para seguir com análise
   ordinal (`arredondar = TRUE`). Depois da recentralização todos têm a mesma
   média: o escore não compara pessoas pelo nível médio. As aulas usam a função
   `recode_for_acq` do LabAPE, carregada da internet com `source()`; a skill
   reimplementa a conta e não baixa código.
3. **Interceptos randômicos** (Maydeu-Olivares & Coffman, 2006): fator de método
   com carga fixa 1 em todos os itens, na métrica original (negativos sem
   inverter), ortogonal aos fatores de conteúdo. A variância dele é a variância
   de aquiescência; os evals conferem que a sintaxe gerada dá o mesmo
   qui-quadrado da sintaxe das aulas.
4. **MIMIC**: cada item regredido no índice, com os fatores de conteúdo
   ortogonais ao índice (sintaxe do curso).

O balanceamento é conferido antes: fator sem item negativo para a função (o
índice confundiria traço e aquiescência); número diferente de positivos e
negativos vira aviso (controle parcial). Primi et al. (2019) comparam os métodos
clássicos e modernos de medir e corrigir aquiescência.

## 2. Qual controle usar

- Para descrever e decidir se há aquiescência: variância do fator nos
  interceptos randômicos e a mudança de ajuste e de cargas entre os modelos.
- Para a estrutura fatorial: AFE ou AFC nos dados recentralizados, ou o modelo
  de interceptos randômicos direto.
- Para escores: o escore de conteúdo do modelo de interceptos randômicos, ou o
  escore somado de uma escala balanceada (que já cancela a aquiescência). Escore
  recentralizado não compara pessoas pelo nível médio.

## 3. Desejabilidade social

O curso trabalha três peças, todas com funções da skill ou com lavaan direto:

- **Juízes de desejabilidade**: juízes avaliam o quanto cada item é desejável;
  a concordância sai de `concordancia_juizes()` (ICC ou Fleiss, como nas aulas,
  que usam `irr::kappam.fleiss` e o ICC).
- **Escala neutralizada** (Bäckström et al., 2009): itens reescritos em
  formulação neutra; compara-se a estrutura e o fator geral entre a versão
  original e a neutralizada (AFC com `ajustar_afc()`, invariância com
  `invariancia()`).
  As aulas também testam a média das avaliações contra o ponto médio da escala de
  desejabilidade (`teste_t_uma_amostra()`, mu = 5) e a correlação entre polo do
  traço e desejabilidade (`correlacoes()`), rodada a rodada de juízes.
- **Quádruplas e fator de desejabilidade** (v1.5): `controle_desejabilidade(dados,
  itens, polo_ds)` compara, no mesmo banco, traço só; traço + interceptos randômicos
  (aquiescência); traço + fator de desejabilidade ortogonal com o sinal de cada carga
  preso pelo polo de desejabilidade do item (desejável > ,01; indesejável < -,01), como
  nas aulas "Estrutura interna das quádruplas". `polo_ds` vem dos juízes (+1, -1, ou 0
  para item neutro, que fica fora do fator). Com `criterio`, a validade do traço contra
  a medida externa (o ACME das aulas) sai no modelo sem e com controle. Nos evals
  simulados, o controle aproxima as cargas do traço e o β do critério dos valores da
  população, e todas as cargas de desejabilidade saem com o sinal do polo.
- **MIMIC com uma medida de desejabilidade**: como o MIMIC de aquiescência,
  trocando o índice pelo escore da escala de desejabilidade. Sem função própria:
  sintaxe lavaan análoga à de `controle_aquiescencia()` (itens `~ escore_ds`,
  fatores de conteúdo ortogonais ao escore), conferida com `ajustar_afc()`. A
  comparação multigrupo das versões quádrupla e neutralizada sai de `invariancia()`.

## 4. Escolha forçada (TRI thurstoniana)

Blocos de itens em que a pessoa ordena ou escolhe; os escores brutos são
ipsativos e não comparam pessoas. O modelo thurstoniano de Brown e
Maydeu-Olivares (2011) recupera escores normativos. As aulas usam o pacote
`thurstonianIRT` (Bürkner, 2019): `set_block()` para definir os blocos,
`make_TIRT_data()` para o formato longo de comparações binárias,
`fit_TIRT_lavaan(estimator = "ULSMV")` e `predict()` para os escores; a
precisão vem da confiabilidade empírica pela matriz de informação, com o código
da aula.

`escolha_forcada(dados, blocos, formato)` (v1.5) faz essa rota: `blocos` é uma lista com
`itens`, `tracos` e `sinais` por bloco; `formato = "ranks"` recebe os postos dentro de
cada bloco (a recodificação para comparações binárias que as aulas fazem com o pacote
RecodeFCit, do GitHub, fica com o `make_TIRT_data`), `"pairwise"` recebe as comparações
já binárias. Saem o ajuste (ULSMV), as correlações entre traços, as cargas por comparação
e a confiabilidade empírica por traço pela fórmula das aulas (informação de cada
comparação no θ estimado, EP² = diag[(I + Φ⁻¹)⁻¹], ρ = var(θ̂)/[var(θ̂) + média(EP²)]),
com os θ calculados e descartados. Cada traço só é identificado até o sinal: a função
compara as cargas com o sinal declarado dos itens e reflete o traço que convergiu
invertido (as correlações saem já no sentido declarado). Com itens todos de mesmo sinal
os escores ficam ipsativos e a recuperação dos traços cai (num teste com 3 traços e
sinais iguais, correlações estimadas de ,77, ,57 e ,03 para ,30, 0 e -,20 simuladas); use
itens de sinais mistos. Blocos de 3 ou mais itens têm redundâncias entre comparações, e
Brown e Maydeu-Olivares (2012) corrigem os graus de liberdade por elas; a skill não aplica
a correção, então χ² e índices que dependem dos gl ficam aproximados. Nos evals (5 traços,
30 itens, sinais mistos), as correlações entre traços e a confiabilidade empírica batem
com a simulação.

## 5. Relato

- Dizer qual controle foi usado, com o estimador, e mostrar o ajuste com e sem
  controle (tabela de `controle_aquiescencia()`).
- Relatar a variância de aquiescência com o p e o tamanho das cargas de conteúdo
  antes e depois.
- Recentralização: dizer o arredondamento e o truncamento.
- Frase pronta em `$apa`; o índice e os dados recentralizados não vão para a
  conversa nem para o manuscrito, só agregados.
