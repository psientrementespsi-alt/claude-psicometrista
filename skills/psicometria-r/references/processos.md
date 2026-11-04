# Processos: moderação, mediação e processo condicional

## Sumário

- 18 (linha 30)
- 36 (linha 48)
- 86 (linha 98)
- 125 (linha 137)
- 161 (linha 173)
- 183 (linha 195)
- 196 (linha 208)
- 228 (linha 240)
- 237 (linha 249)

Funções em `scripts/est_processos.R`; testes em `scripts/evals_processos.R`. Estimação por
MQO e bootstrap percentílico próprio, no formato dos modelos do PROCESS (numeração de
Hayes, 2022), sem usar a macro.

Índice
1. Qual análise (regra de decisão)
2. Funções e argumentos
3. Moderação (modelo 1)
4. Mediação simples, paralela e serial (modelos 4 e 6)
5. Mediação moderada (modelos 7, 14 e 15)
6. Pressupostos e o que a função checa
7. Como relatar
8. Armadilhas
9. Fontes e divergências

## 1. Qual análise (regra de decisão)

A técnica nasce da pergunta. Antes de rodar, desenhar o diagrama
conceitual: quem é X, M, W e Y e por que cada seta existe; se não dá para desenhar, a
hipótese ainda não está pronta.

| Pergunta | Marca no texto | Função | Modelo PROCESS |
|---|---|---|---|
| X afeta Y? | sem mecanismo nem condição | `regressao_linear` (outro módulo) | |
| Por qual caminho? | "porque", "por meio de" | `mediacao()` | 4 (um ou mais M em paralelo), 6 (serial) |
| Depende de quê? | "mais forte entre", "depende de" | `moderacao()` | 1 |
| O mecanismo depende de W? | combinação das duas | `mediacao_moderada()` | 7 (W em a), 14 (W em b), 15 (W em b e c') |

Paralela ou serial é decisão teórica: os mediadores se afetam?. Onde W entra
(a, b, c') também vem da teoria e muda a história contada. Número de modelo
não substitui teoria. Y dicotômica não é coberta: a função recusa e
orienta para regressão logística.

## 2. Funções e argumentos

```r
moderacao(dados, y, x, w, covariaveis = NULL, centralizar = FALSE,
          valores_w = c("media_dp", "percentis"), johnson_neyman = TRUE,
          ep = c("classico", "HC3"), nivel = .95, jn_prop_min = .05,
          casas = 3, decimal = ".", idioma = "pt")
mediacao(dados, x, m, y, covariaveis = NULL, tipo = c("paralela", "serial"),
         reamostras = 10000, semente = 2026, nivel = .95, padronizado = FALSE,
         contrastes = TRUE, ep = c("classico", "HC3"), proporcao_mediada = FALSE,
         delineamento = c("transversal", "longitudinal", "experimental"),
         casas = 3, decimal = ".", idioma = "pt")
mediacao_moderada(dados, x, m, y, w, modelo = 7, covariaveis = NULL,
                  valores_w = c("media_dp", "percentis"), centralizar = FALSE,
                  reamostras = 10000, semente = 2026, nivel = .95, ep = c("classico", "HC3"),
                  johnson_neyman = TRUE, jn_prop_min = .05,
                  delineamento = c("transversal", "longitudinal", "experimental"),
                  casas = 3, decimal = ".", idioma = "pt")
```

Comum às três: célula de texto vazia ("" ou só espaços, comum em CSV) vira ausente;
deleção listwise nas variáveis da análise (`n_excluidos`, com nota em `$decisoes` quando
houve vazios); X numérico ou fator de dois níveis (vira 0/1, primeiro nível = 0);
covariáveis numéricas ou fatores (indicadores) entram em todas as equações e ficam na
média nas predições; `valores_w` aceita também um vetor numérico de valores teóricos ou
clínicos; `casas` controla coeficientes e efeitos na frase (o livro usa 3 casas,
), e um valor não nulo que arredondaria para zero sai com 2
algarismos significativos (ex.: índice -0,00049; limite de IC -0,00030), com nota em
`$decisoes` sugerindo reescalar a variável; t e F saem com 2. Saída `pr_resultado`: `tabela` (coeficientes de todas as equações, com coluna
`caminho` na mediação), `efeito`, `posthoc` (efeitos condicionais), `pressupostos`,
`decisoes`, `avisos`, `apa`.

Grupo gravado como número (1/2/3, comum em banco vindo do SPSS ou de planilha): X, W ou
covariável numérica com 3 a 7 valores distintos, todos inteiros, entra como contínua e
gera aviso em `$avisos` (mesma regra de `est_associacao.R`; códigos numéricos diretos para
3 ou mais categorias são proibidos no livro). X ou W 0/1 continuam aceitos sem
aviso; item Likert dispara o aviso, que nesse caso pode ser ignorado. O que fazer depende
do papel: W da `moderacao()` e covariável viram fator antes (`factor()`), com indicadores;
X de mediação e de mediação moderada e W da mediação moderada com 3 ou mais categorias não
são cobertos (rodar com dois grupos por vez).

Número de reamostras na frase: separador de milhar só quando não se confunde com o
decimal ("10.000 reamostras" com `decimal = ","`; "10,000 resamples" em inglês com
`decimal = "."`); nos casos mistos, sem separador ("10000 resamples" ao lado de
"ab = 0,241").

Índices sem IC (ΔR², f² e R² da moderação; ΔR² das interações na mediação moderada)
ganham a nota "sem IC: não calculado" em `$decisoes` (regra 6 da skill), para levar ao
relato.

## 3. Moderação (modelo 1)

Modelo Y = b0 + b1·X + b2·W + b3·X·W (+ covariáveis); X, W e o produto entram sempre
juntos. W contínua, binária ou categórica (indicadores, primeira categoria de
referência; o teste da moderação vira o F conjunto dos produtos). W categórica gravada
como 1/2/3 precisa virar fator antes; numérica, ela seria sondada em média ± 1 DP de
códigos, sem leitura, e a função avisa (seção 2).

- Teste da moderação: b3 com t e p; ΔR² ao retirar o(s) produto(s) e F(gl1, gl2); f² =
  ΔR²/(1 - R²). Com EP clássico, F = t² do produto. `principal` = ΔR².
- Centralização (`centralizar = TRUE`): só X e W contínuas (binárias e categóricas
  ficam). Muda apenas b1 e b2; b3, EP(b3), R² e F não mudam. A função
  ajusta sempre a outra parametrização, guarda em `$coeficientes_alternativos` e escreve
  em `$decisoes` o b1 da outra versão e a conferência de invariância. Centralizar não
  conserta multicolinearidade essencial; o VIF alto do produto é artefato (;
  Echambadi & Hess, 2007; Dalal & Zickar, 2012), por isso a função não checa VIF.
- Sondagem (`$posthoc`): efeito condicional θ(w) = b1 + b3·w com EP, t, p e IC t, em
  média e ±1 DP de W (padrão do livro), nos percentis 16/50/84
  (`valores_w = "percentis"`, padrão do PROCESS e do curso; `quantile` tipo 7) ou em
  valores informados. Valor fora da faixa observada gera aviso (extrapolação).
- Johnson-Neyman (W contínua): fronteiras por solução analítica de
  (b1 + b3·w)² = t²crít·Var(θ(w)), isto é, a equação de segundo grau
  (b3² - t²V33)w² + 2(b1b3 - t²V13)w + (b1² - t²V11) = 0, guardando as raízes dentro da
  faixa observada. `$johnson_neyman` traz `pontos` (W, efeito, t, % de casos abaixo e
  acima), `regioes` (faixas, significância, sinal, n e % de casos) e `grade` (20 pontos
  entre mínimo e máximo de W mais as fronteiras, como o PROCESS). Região significativa com
  menos de `jn_prop_min` dos casos gera aviso (sugestão do curso).
- `ep = "HC3"`: EP robusto (igual a `sandwich::vcovHC(type = "HC3")`) em coeficientes,
  efeitos condicionais e J-N; o F da interação passa a ser o de Wald robusto, como no
  PROCESS, e o ΔR² fica descritivo.
- `$grade_grafico`: Y previsto para X (0/1, ou média ±1 DP, ou percentis) cruzado com os
  valores de W, covariáveis na média. Uma linha por nível de W no gráfico; linhas
  paralelas = sem moderação.
- Padrão descritivo (`$padrao`): "amplifica", "atenua" ou "inverte", lido nos valores
  sondados; é leitura, não teste.

Se b3 não é significativo, a frase APA para no teste da interação; os efeitos
condicionais ficam no objeto só como descrição.

## 4. Mediação simples, paralela e serial (modelos 4 e 6)

Equações por MQO: cada M ~ X (+ mediadores anteriores, se serial) (+ covariáveis);
Y ~ X + todos os M (+ covariáveis); e o modelo do efeito total Y ~ X (+ covariáveis).
Caminhos rotulados em `tabela$caminho`: a (a1, a2...), d21 (serial: M1 → M2), b, c', c.

- Efeito indireto: a·b na simples; a_j·b_j por mediador na paralela; na serial todos os
  caminhos ordenados (com 2 M: a1·b1, a2·b2, a1·d21·b2). Indireto total = soma.
  `principal` = efeitos específicos (foco recomendado).
- Inferência: IC percentílico por bootstrap de casos, 10.000 reamostras e semente 2026
  por padrão. Os índices são sorteados como no `boot::boot`, e o IC usa a
  mesma interpolação do `boot.ci(type = "perc")`: mesma semente, mesmo IC (testado). Com
  n × reamostras acima de 2e7 a função sorteia uma reamostra por vez (poupa memória) e
  avisa em `$avisos`: o IC segue reprodutível pela semente, mas não igual ao do
  `boot::boot`. Efeitos
  direto e total com IC pela distribuição t (como o PROCESS).
- Identidade c = c' + indireto total conferida e escrita em `$decisoes` e
  `$identidade`; em MQO com os mesmos casos ela fecha exatamente.
- Contrastes (`contrastes = TRUE`, com 2+ efeitos específicos): diferença bruta entre
  pares, com IC bootstrap. Dois efeitos significativos não implicam diferença
  significativa. Com sinais opostos (ambos com IC sem o zero) a diferença
  bruta mistura direção e magnitude; a função avisa.
- `padronizado = TRUE`: X contínuo, completamente padronizado ab·DP_X/DP_Y; X binário,
  parcialmente padronizado ab/DP_Y; DP recalculados em cada reamostra. O padrão segue o
  livro: coeficientes brutos.
- Sobel: só para mediação simples, em `$sobel`, como comparação histórica.
- `proporcao_mediada = TRUE` (só simples): ab/c, sem IC, com aviso; não calculada quando
  ab e c têm sinais opostos. Pede n grande (MacKinnon et al., 1995).
- Mediação inconsistente: aviso só quando indireto e direto têm sinais opostos e os
  dois IC excluem o zero; |c'| pode superar |c| e c pode ser não significativo com efeito
  indireto real. c' de sinal oposto mas compatível com zero (o caso comum de
  c' perto de zero) não gera aviso. O mesmo critério vale para indiretos específicos de
  sinais opostos.
- `delineamento`: frase de alcance causal em `$avisos` (transversal é o padrão;
  Maxwell & Cole, 2007).

## 5. Mediação moderada (modelos 7, 14 e 15)

Um mediador, W numérica (contínua ou binária) ou fator de dois níveis.

| Modelo | Equação de M | Equação de Y | Indireto condicional | Índice |
|---|---|---|---|---|
| 7 | X, W, X·W | X, M | (a1 + a3·W)·b | a3·b |
| 14 | X | X, M, W, M·W | a·(b1 + b3·W) | a·b3 |
| 15 | X | X, M, W, X·W, M·W | a·(b1 + b3·W); direto c'1 + c'3·W | a·b3 |

- Índice de mediação moderada e indiretos condicionais com IC bootstrap percentílico
  (valores de W fixos nas reamostras). O índice é a inclinação do indireto em W: a
  diferença entre dois indiretos condicionais é índice × distância entre os W. O livro
  exige indiretos condicionais com incerteza e um índice da variação ao longo de W
; ambos têm `principal = TRUE`.
- `$posthoc`: o caminho moderado (a no 7; b no 14 e 15) em cada valor de W, com t; no 15,
  também os efeitos diretos condicionais. `$testes_interacao`: b, t, p, ΔR² e F de cada
  produto. J-N no caminho moderado.
- Interação significativa num caminho é necessária, mas não suficiente: a
  conclusão vem do IC do índice. Centralizar (`centralizar = TRUE`) atinge só as
  contínuas que formam produtos e não muda índice nem indiretos condicionais.

## 6. Pressupostos e o que a função checa

| Checagem | Como | Leitura automática |
|---|---|---|
| Heterocedasticidade | Breusch-Pagan estudentizado (igual a `lmtest::bptest`) por equação | p < .05: sugerir `ep = "HC3"`; o IC bootstrap não depende disso |
| Normalidade dos resíduos | Shapiro-Wilk (3 ≤ n ≤ 5000) | com n ≥ 300 o teste acusa desvios irrelevantes; olhar Q-Q plot |
| Colinearidade perfeita | posto da matriz | erro com os termos envolvidos |
| Célula de texto vazia | "" ou espaços viram NA antes do listwise | nunca vira categoria nem referência |
| Y dicotômica, mediador binário | contagem de valores | erro com orientação |

Linearidade, independência e confundidores não medidos não são testáveis aqui; a frase de
delineamento lembra o limite causal.

## 7. Como relatar

```r
med <- mediacao(dados, x = "Condition", m = "MediaInfluence", y = "Intention", decimal = ",")
relatar(med)
```

Frase gerada (Presumed Media Influence, biblioteca do JASP):

> O efeito indireto de Condition sobre Intention por meio de MediaInfluence foi ab = 0,241
> (IC 95% [0,005; 0,519]; bootstrap percentílico, 10.000 reamostras, semente 2026), com
> a = 0,477 (p = ,045) e b = 0,506 (p < ,001). O efeito direto foi c' = 0,254; EP = 0,256;
> t(120) = 0,99; p = ,322. O efeito total foi c = 0,496; EP = 0,278; t(121) = 1,79; p = ,077.

```r
mod <- moderacao(dados, y = "Justification", x = "Framing", w = "Skepticism", decimal = ",")
```

> A interação entre Framing e Skepticism na predição de Justification foi significativa,
> b = 0,201; EP = 0,055; t(207) = 3,64; p < ,001; ΔR² = ,048; f² = 0,064. O efeito
> condicional de Framing sobre Justification foi b = -0,292 (p = ,068) em Skepticism =
> 1,345 (-1 DP); b = 0,117 (p = ,297) em Skepticism = 3,378 (média); b = 0,526 (p = ,001)
> em Skepticism = 5,411 (+1 DP). Pela técnica de Johnson-Neyman, o efeito de Framing foi
> significativo quando Skepticism ≤ 1,171 (negativo; 6,6% dos casos) ou Skepticism ≥
> 3,934 (positivo; 32,2% dos casos).

Na mediação moderada a frase traz o teste de cada produto, o índice com IC e os indiretos
condicionais. Relatar sempre significância, magnitude e precisão; no Método,
número de reamostras, semente, tipo de IC e tratamento de ausentes. Redação para
transversal: associativa, por exemplo "a associação foi explicada, em termos estatísticos,
por M".

## 8. Armadilhas

- Rotular "completa/parcial" ou exigir c significativo.
- Ler b1 e b2 como efeitos principais com interação; b1 é o efeito de X com W = 0
. Testar o produto sem X e W, ou usar subgrupos.
- Concluir mediação moderada pela interação num caminho, sem índice.
- Proporção mediada como c'/c (erro do curso); J-N com poucos casos na região; W fora da
  faixa observada; mudar a ordem serial sem teoria.

## 9. Fontes e divergências

Livro (Vol. 6, `biblioteca/md6.md`): caminhos e identidade;
coeficientes brutos; bootstrap percentílico, 10.000, semente 2026;
três informações do relato; sem completa/parcial; c não precisa ser
significativo; Sobel histórico; múltiplos mediadores e contraste
; modelo com produto; centralização; ΔR² e f²
; média ±1 DP; spotlight e J-N; processo condicional
; limite causal por delineamento.

Curso (Psicometria Online Academy; `fontes/spec_T5.md` e `notion_g3.md`, `notion_g4.md`):
páginas Moderação, Prática de Moderação, Descrição da Moderação, Mediação, Prática e
Reporte da Mediação, Mediação Paralela e Serial (teoria e R), Processo Condicional
(modelo 15) e Moderação no R (HC3).

Literatura (DOI conferido): Bauer & Curran (2005); Hayes & Matthes (2009); Spiller et al.
(2013); Echambadi & Hess (2007); Dalal & Zickar (2012); Preacher & Hayes (2008); Hayes &
Scharkow (2013); Fritz, Taylor & MacKinnon (2012); MacKinnon et al. (2002); MacKinnon,
Warsi & Dwyer (1995); Rucker et al. (2011); Maxwell & Cole (2007). Hayes (2022), livro
sem DOI (ISBN 9781462549030), para numeração dos modelos e índice de mediação
moderada.

| Tema | Livro | Curso / literatura | Decisão |
|---|---|---|---|
| Reamostras | 10.000, semente 2026 | PROCESS e curso: 5.000 | Padrão 10.000; argumento `reamostras` |
| Tipo de IC | percentílico | curso relata BCa em texto | Só percentílico (Fritz et al., 2012; Hayes & Scharkow, 2013); rótulo sempre explícito |
| Valores de W | média ±1 DP | curso e PROCESS: P16/50/84; Spiller et al. (2013): valores focais ou J-N | Padrão média ±1 DP; `"percentis"` e vetor numérico; aviso quando o valor cai fora da faixa |
| Centralização | opcional, só leitura | curso centraliza por padrão | Padrão escala original (como); a outra versão sai sempre em `$coeficientes_alternativos` |
| EP robusto | não trata | curso usa HC3 (PROCESS hc = 3) | `ep = "HC3"` opcional; F da interação vira Wald robusto |
| Efeitos condicionais | sempre na demonstração | curso: só com b3 significativo | Calculados sempre; frase só com b3 significativo |
| Proporção mediada | não trata | curso: c'/c numa página (erro), ab/c noutra | ab/c, opcional, só simples, com aviso (MacKinnon et al., 1995) |
| Sobel | histórico | curso: "infla erro tipo I" (erro; é conservador, MacKinnon et al., 2002) | Só em `$sobel`, fora da inferência |
| Contraste na serial | não trata | curso: sem contraste na serial | Permitido (o PROCESS permite), com a ressalva de leitura |
| Numeração dos modelos | não cita | curso usa PROCESS 1, 4, 6, 7, 14, 15 | Numeração de Hayes só para interoperabilidade |
| Percentis | não trata | PROCESS: algoritmo próprio, não conferido | `quantile` tipo 7; pode diferir na 3ª casa do PROCESS com W de muitos valores |

Reprodução do livro (evals): mediação (a, b, c, c', ab, p de b, c e c', valores sem
arredondar, Sobel z = 1,89 e p = ,059) e moderação (b1, b2, b3 e p, R², ΔR², f², b1
centralizado, efeitos condicionais e p) batem até a última casa impressa. Os p do livro
são de teste t por MQO: o módulo Process do JASP (lavaan, MV e z) daria p = ,072 para c,
contra ,077 no livro. O IC bootstrap do livro, [0,006; 0,521], não é reproduzido
exatamente: com semente 2026 e o esquema do `boot::boot` sai [0,005; 0,519], diferença
dentro do erro de Monte Carlo (o próprio livro relata [0,023; 0,541] noutra execução); o
gerador usado no livro não é conhecido. Do curso, com o mesmo banco: indiretos paralelos
0,203 e 0,189, total 0,392, contraste 0,014 com IC contendo zero, e efeito serial 0,049.
Johnson-Neyman do banco Climate Change Disasters: fronteiras em ceticismo 1,171 e 3,934.
