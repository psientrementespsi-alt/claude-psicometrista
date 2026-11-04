# Modelagem por equações estruturais além da AFC (est_mee.R, v1.5)

Fontes: cursos "Análise Fatorial e Modelagem por Equações Estruturais" (JASP), "MPLUS" e
"R: Análise Fatorial e MEE" da pasta Psicometria-R (scripts `MEE_MRP.R`,
`Moderação_categórica.R`, `Interação de Fatores Latentes.R`, `Sintaxe_MEE.Rmd`) e as
sintaxes `Latent growth curve.R` e `Latent growth curve com SEM.R` do acervo (pasta 07).
Literatura: Little et al. (2006), Marsh et al. (2004), Schoemann e Jorgensen (2021),
Cheung e Rensvold (2002), Wu e Estabrook (2016).

Sumário: 1. Sequência · 2. `mee()` · 3. `mee_multigrupo()` · 4. `interacao_latente()` ·
5. `curva_crescimento()` · 6. `invariancia_longitudinal()` · 7. Divergências · 8. Relato

## 1. Sequência

Como nas aulas: cada instrumento na AFC (`ajustar_afc()`), depois o modelo de mensuração
com todos os fatores correlacionados, depois o estrutural. Mudança de ajuste entre o de
mensuração e o estrutural mostra o custo das restrições estruturais. Estimador e versão dos
índices como na AFC (`afc-ordinal.md`): WLSMV com item ordinal, MLR com FIML nos demais.

## 2. `mee(modelo, dados, bootstrap = 0)`

- Sintaxe lavaan completa: `=~` para a medida, `~` para as regressões, rótulos (`a*X`) e
  `:=` para efeitos definidos (indireto `a*b`, total `a*b + c`).
- Tabela: cada regressão e efeito definido com b, EP, z, p (da estimativa não
  padronizada) e β com IC. Efeito principal: o primeiro efeito definido (o indireto), ou o
  primeiro caminho se não houver.
- IC: método delta por padrão; `bootstrap = 1000` (ou mais) dá o IC percentil das
  estimativas padronizadas, a rota das aulas (`bootstrapLavaan` com `standardizedSolution`).
  Efeito indireto pede bootstrap: o produto a × b é assimétrico. Semente fixa; reamostras
  sem convergência descartadas e contadas.
- Extras: ajuste, R² das endógenas, cargas padronizadas. Solução imprópria (carga > 1, R²
  fora de [0, 1]) vira aviso.
- `pct_mediada := ind/total` das aulas fica fora do padrão: a proporção mediada é instável
  com efeito total pequeno e muda de sinal com efeitos inconsistentes (`processos.md`).

## 3. `mee_multigrupo(modelo, dados, grupo)`

Moderação por variável categórica (aula "Moderação (moderadora categórica)").
- Passos: configural → métrica (cargas) → escalar (cargas + limiares no ordinal, ou +
  interceptos no contínuo) → regressões iguais. Teste de diferença (escalonado com WLSMV e
  MLR) e ΔCFI entre vizinhos.
- As aulas igualam variâncias e covariâncias latentes antes das regressões.
  `igualar_latentes = TRUE` reproduz; o padrão da skill compara as regressões a partir do
  escalar, porque a igualdade de caminhos não padronizados só pede a medida invariante, e
  variâncias latentes diferentes entre grupos não são moderação.
- Teste por caminho: o modelo com todas as regressões iguais contra o mesmo modelo com só
  aquele caminho livre (`group.partial`), como nas aulas. Tabela de β por grupo no modelo de
  referência.

## 4. `interacao_latente(dados, x, z, y)`

Interação entre fatores latentes por indicadores-produto (aula "Moderação entre variáveis
latentes"): `semTools::indProd` com centralização residual (Little et al., 2006),
produtos pareados pela ordem dos itens (Marsh et al., 2004), fatores com variância 1 e
estrutura de médias, estimação ML/MLR. Saem o coeficiente da interação (b e β), o ΔR² de Y
contra o modelo sem a interação e as inclinações simples de Y em X a -1, 0 e +1 DP de Z
(`probe2WayRC`). Exige o mesmo número de itens em X e Z; com números diferentes, montar
parcelas ou escolher os pares pelos maiores pesos. Item ordinal com poucas categorias pede
LMS (Mplus) e sai do escopo testado.

## 5. `curva_crescimento(dados, ondas)`

Curva de crescimento latente (`lavaan::growth`), banco largo (uma coluna por onda).
- `forma = "linear"`: cargas da inclinação fixas nos `tempos` (0, 1, 2... por padrão; usar
  os intervalos reais se as ondas não forem equidistantes).
- `forma = "livre"` (base latente): primeiras cargas fixas em 0 e 1, demais livres; a
  inclinação média é a mudança da onda 1 para a 2 e as cargas livres dizem a forma da
  curva. Pede 4 ou mais ondas para ter gl.
- `residuos_adjacentes = TRUE` correlaciona resíduos de ondas vizinhas (a segunda curva das
  aulas); com 4 ondas e base latente o modelo fica saturado (gl = 0), e a função avisa.
- `residuos_iguais = TRUE` iguala as variâncias residuais: o modelo fica equivalente ao
  modelo misto com intercepto e inclinação aleatórios (conferido contra `lmer` nos evals).
- `preditores`: covariáveis do intercepto e da inclinação.
- Saem médias, variâncias e covariância de intercepto e inclinação; efeito principal: a
  inclinação média.

## 6. `invariancia_longitudinal(dados, ondas)`

Passo prévio da curva com variável latente em cada onda (segunda ordem) e de qualquer
comparação de médias latentes entre ondas. `ondas` é uma lista com os itens de cada onda,
na mesma ordem. Passos como em `invariancia()` (ordinal: limiares → cargas → interceptos,
Wu & Estabrook, 2016; contínuo: métrica → escalar → estrita), com `longFacNames` do
`measEq.syntax` e resíduos do mesmo item correlacionados entre ondas. Critério: ΔCFI ≤ -0,01
lido junto com o teste de diferença.

A curva de segunda ordem (fatores por onda → intercepto e inclinação) fica fora das
funções: montar a partir da sintaxe do modelo escalar (`$modelos$interceptos` ou
`$modelos$escalar`, `lavaan::partable`), e rodar com `lavaan::sem`, com a ressalva de rota
fora do escopo testado.

## 7. Divergências com as aulas

- `Latent growth curve com SEM.R`: o roteiro manda trocar as linhas `0*1` dos interceptos
  por `NA*1` e todas as variâncias dos fatores para o rótulo `psi.1_1`, mas a sintaxe que
  roda (`mod_curva`) mantém `0*1` e os rótulos `psi.1_1`, `psi.2_2`, `psi.3_3`. O modelo
  rodado não é o descrito; quem reproduzir precisa decidir qual dos dois quer e dizer.
- `Latent growth curve.R` lê os dados de uma URL da UCLA (`gpa_wide.csv`); a skill não
  depende de download.
- `MEE_MRP.R` roda o bootstrap com `parallel = "multicore"` e `ncpus = 8`, que no Windows não
  paraleliza e no Mac depende dos núcleos; a skill roda sequencial, com semente.

## 8. Relato

Ajuste com a versão dos índices, β de cada caminho com IC (dizer se bootstrap percentil e
quantas reamostras), R² das endógenas. Efeito indireto: β com IC bootstrap, sem rótulo
completa/parcial. Multigrupo: a sequência de passos com ΔCFI e o teste de diferença, e o
teste de cada caminho. Interação latente: b, β, p, ΔR² e as inclinações simples. Curva:
médias e variâncias de intercepto e inclinação, com a forma da curva e o tratamento dos
resíduos.
