# Preparação de dados, figuras e tradução de Mplus para lavaan

Base: a "Segunda Camada em R" do curso (linguagem e scripts, tidyverse, ggplot2)
e os inputs de Mplus das aulas. Nada aqui tem função própria; tudo entra como
bloco registrado do script (`references/script-reproduzivel.md`).

## Preparação dos dados

- O arquivo bruto não é editado à mão: recodificação, inversão de itens,
  exclusões e criação de escores são blocos registrados, com o critério escrito
  no código (`subset(dados, tempo_resposta >= 300)`), nunca lista de IDs.
- Inversão de item Likert: `dados$y2_inv <- (min + max) - dados$y2`, com min e
  max da escala, não os observados. Conferir depois com
  `analise_classica_itens()` (RIR negativo denuncia item não invertido).
- Fatores com rótulos explícitos (`factor(x, levels = c(1, 2), labels =
  c("controle", "intervenção"))`) antes de qualquer função de grupo; a
  referência das funções é o primeiro nível ou o argumento `referencias`.
- Ausentes: contar por variável antes (`colSums(is.na(dados))`, que é agregado)
  e decidir a rota (listwise, FIML, imputação) antes de ver o resultado; ver
  `references/afc-ordinal.md` (seção de ausentes).
- tidyverse ou R base dão o mesmo resultado; o que importa é o bloco estar no
  script. `dplyr::mutate()` e `across()` servem para recodificar muitos itens.

## Figuras

- Gráfico de agregados (médias com IC, curvas de informação, perfis latentes,
  floresta de meta-análise, curvas de sobrevida) pode voltar para a conversa.
- Gráfico com um ponto por pessoa (dispersão, boxplot com pontos, mapa de Wright
  com pessoas) fica no computador do usuário: grava o arquivo e não vem para a
  conversa.
- Para manuscrito em APA 7: `ggplot2::theme_classic()`, fonte sem serifa de 10 a
  12 pt, eixos com nome e unidade, legenda dentro da figura quando couber, sem
  título dentro da imagem (o título vai no texto da figura) e
  `ggsave("figura1.png", width = 16, height = 10, units = "cm", dpi = 300)`.
- Barras de erro: dizer se são IC 95% ou EP, no texto da figura.
- Cores: paleta segura para daltonismo (`scale_colour_viridis_d()`) ou tons de
  cinza; não depender só de cor para separar grupos. O curso de ggplot2 mostra as
  paletas de publicação dos pacotes `scico` e `ggsci`; vale o mesmo critério.
- Raincloud (aula "Raincloud plot"; Allen et al., 2019): densidade, boxplot e pontos
  brutos por grupo. Tem um ponto por pessoa, então é gráfico do computador do usuário;
  para a conversa, a versão sem os pontos (densidade e boxplot).
- Diagrama de AFC ou MEE (aulas "Visualização com semPlot"): `semPlot::semPaths(fit,
  "std", style = "lisrel", intercepts = FALSE, thresholds = FALSE, residuals = FALSE,
  structural = TRUE, layout = matriz)`, como no `MEE_MRP.R` (Epskamp, 2015). É figura de
  agregados (cargas e caminhos). Para o manuscrito, conferir pela diagram-auditor-apa7
  quando for diagrama de mediação.
- Gráfico de efeitos relativos dos não paramétricos complexos (RTE por tempo e grupo,
  com IC do mctp.rm), como nas aulas: pontos e linhas por grupo, eixo y de 0 a 1, linha
  de referência em 0,5.

## Mplus para lavaan

| Mplus | lavaan |
|---|---|
| `F1 BY y1-y4;` | `F1 =~ y1 + y2 + y3 + y4` |
| `F1 BY y1* y2-y4; F1@1;` | `F1 =~ NA*y1 + y2 + y3 + y4` e `F1 ~~ 1*F1` (ou `std.lv = TRUE`) |
| `y1@1` (carga fixa) | `1*y1` |
| `F1 WITH F2;` / `F1 WITH F2@0;` | `F1 ~~ F2` / `F1 ~~ 0*F2` |
| `y ON x;` | `y ~ x` |
| `[y1];` (intercepto ou média) | `y1 ~ 1` |
| `y1-y4 (a);` (igualdade) | mesmo rótulo: `a*y1 + a*y2 + a*y3 + a*y4` |
| `CATEGORICAL ARE y1-y12;` | `ordered = c("y1", ..., "y12")` |
| `ESTIMATOR = WLSMV;` | `estimator = "WLSMV"` (parametrização delta nos dois por padrão) |
| `MODEL INDIRECT: y IND x;` | rótulos nos caminhos e `ind := a*b` |
| `F1-F3 BY y1-y12 (*1);` (ESEM, geomin) | `efa("efa1")*F1 + efa("efa1")*F2 + efa("efa1")*F3 =~ y1 + ... + y12`, `rotation = "geomin"` |
| `OUTPUT: STDYX MODINDICES;` | `summary(fit, standardized = TRUE)`; `modindices(fit)` |
| `TYPE = MIXTURE;` (LCA/LPA) | `classes_latentes()` / `perfis_latentes()` |
| `GROUPING IS sexo (1 = m 2 = f);` | `group = "sexo"` (invariância: `invariancia()`) |

Diferenças que aparecem: qui-quadrado e índices do WLSMV costumam bater na
segunda casa, mas o CFI e o TLI podem diferir um pouco pela forma de calcular o
modelo de base; o Mplus usa FIML por padrão com dados contínuos e o lavaan só
com `missing = "fiml"`. Ao migrar um input das aulas, comparar o número de
parâmetros livres e os graus de liberdade antes dos índices.
