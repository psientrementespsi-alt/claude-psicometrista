# Classes latentes (LCA) e perfis latentes (LPA)

Funções em `scripts/est_latentes.R`; testes em `scripts/evals_latentes_epidemio.R`
(`classes_distal()`, v1.5, em `scripts/evals_complementos.R`).
Base: o módulo de análise de classes e perfis latentes do curso (Psicometria
Online Academy), que roda LCA no poLCA e LPA no tidyLPA.

## Qual função

| Indicadores | Função | Motor |
|---|---|---|
| categóricos (sim/não, Likert curto tratado como categoria) | `classes_latentes(dados, itens, k_max = 5, nrep = 10)` | poLCA |
| contínuos (escores de escalas) | `perfis_latentes(dados, vars, k_max = 5, modelo = "EEI")` | mclust |

`modelo` do LPA: `"EEI"` (variâncias iguais entre perfis, covariâncias zero; o
modelo 1 do tidyLPA, o das aulas), `"VVI"` (variâncias livres), `"EEE"` e
`"VVV"` (com covariâncias). Começar pelo EEI; modelos mais livres pedem N maior
e convergem pior. `padronizar = TRUE` põe os indicadores em escore z, útil
quando as escalas têm amplitudes diferentes.

## O que sai

Para cada solução de 1 a `k_max`: log-verossimilhança, número de parâmetros,
AIC, BIC, SABIC, entropia relativa, `prob_min`/`prob_max` (as do tidyLPA),
`avepp_min` (menor média da posterior entre os alocados em cada classe, a
tabela "Average Latent Class Probabilities" do Mplus), proporção da menor e da
maior classe, convergência e o LMR entre k - 1 e k classes (Lo et al., 2001, na
implementação do `tidyLPA::calc_lrt`). No LPA, `blrt = 500` acrescenta o BLRT do
mclust (demorado). As probabilidades de resposta por classe ficam em `$perfis`
(LCA) e as médias por perfil em `$medias` (LPA). A classe de cada pessoa fica
em `$classe` e nunca é impressa.

Os evals conferem: AIC e BIC iguais aos do poLCA; LL, BIC, entropia, prob_min e
n_min iguais aos do `tidyLPA::estimate_profiles(models = 1)`; recuperação das
estruturas simuladas.

## Como decidir o número de classes (critérios do curso)

1. Comparativos: menores AIC, BIC e SABIC. Eles tendem a cair sempre com mais
   classes; quando o BIC ainda cai na maior solução testada, a função avisa para
   ampliar `k_max`.
2. Testes de razão de verossimilhança: LMR (e BLRT, se pedido) entre k - 1 e k.
   O curso considera o BLRT o mais robusto, mas quase sempre significativo na
   prática.
3. Adequação, não comparação: entropia ≥ ,80; probabilidade de classificação
   alta; menor classe ≥ 5% da amostra.
4. Interpretabilidade das probabilidades ou médias por classe e a teoria. A
   skill lista candidatas e soluções adequadas em `$decisoes`; a escolha final é
   do pesquisador e vai justificada no relato.

Soluções locais: o poLCA parte de `nrep` pontos aleatórios e fica com o maior
log-verossimilhança; aviso de máximo não replicado ou sem convergência pede
`nrep` e `maxiter` maiores, como nas aulas. Semente fixada (`semente = 2026`).

## Usar as classes depois

Salvar a classe modal e usá-la como variável de grupo ignora o erro de
classificação. O curso apresenta as abordagens de três passos com correção de
viés (BCH e máxima verossimilhança, no tidySEM) e avisa que, com entropia muito
baixa, até elas podem continuar enviesadas.

`classes_distal(resultado, dados, distal, k)` (v1.5) faz os três passos com a correção
BCH (Bolck et al., 2004; Vermunt, 2010) sobre a saída de `classes_latentes()` ou
`perfis_latentes()`: matriz D de erro de classificação, D[t, s] = P(atribuído a s |
classe t), pesos BCH w* = W D⁻¹ (podem ser negativos) e médias (variável contínua) ou
proporções (categórica) por classe, com EP sanduíche por pessoa, Wald global e pares com
Holm. A coluna `media_modal` mostra a média da atribuição ingênua só para comparação. Nos
evals, com entropia ,56, a BCH chega mais perto da diferença simulada que a atribuição
modal; com classificação perfeita, as duas coincidem. Com entropia abaixo de ,60 a
função avisa. O EP não carrega a incerteza do passo 1 (como no tidySEM e no Latent GOLD).
Covariáveis que predizem a classe (o modelo `BCH(..., model = "Y ~ X")` das aulas) e o
método de máxima verossimilhança ficam fora: rodar no tidySEM com a ressalva de escopo.

## Relato

Tabela com todas as soluções testadas (não só a escolhida), o critério que
decidiu, entropia e tamanho da menor classe da solução final, e o gráfico de
perfis (probabilidades ou médias por classe), que não tem dado individual.
Frase pronta em `$apa`.
