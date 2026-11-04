# TRI: ajuste, comparação de modelos e DIF

Funções em `scripts/est_tri.R` (mirt; Chalmers, 2012); testes em
`scripts/evals_tri.R`. Todas devolvem `pr_resultado`.

Índice: 1. Quando usar · 2. Qual modelo · 3. `tri_ajustar()` · 4. Réguas e o que o
curso ensina · 5. Confiabilidade e precisão · 6. `tri_comparar()` · 7. DIF com
`tri_dif()` · 8. Relação com a AFC e a invariância · 9. TRI multidimensional · 10. Fora das funções (CAT,
equacionamento, ancoragem)

## 1. Quando usar

A pergunta é sobre o item (discriminação, em que faixa do traço ele informa, se
funciona igual entre grupos) ou sobre a precisão ao longo do traço, não só sobre
a estrutura fatorial. Unidimensionalidade é pressuposto dos modelos daqui:
testar antes com AFC de um fator ou EGA (`references/redes-ega.md`), e não usar
o ajuste da TRI como prova de unidimensionalidade. Escala multidimensional:
rodar por subescala.

## 2. Qual modelo

| `modelo` | Itens | Modelo | itemtype do mirt |
|---|---|---|---|
| `"rasch"` | dicotômicos | Rasch (1PL com discriminação comum 1) | `Rasch` |
| `"2PL"` | dicotômicos | 2 parâmetros (Birnbaum) | `2PL` |
| `"3PL"` | dicotômicos | 3 parâmetros, com acerto ao acaso | `3PL` |
| `"pcm"` | ordinais | crédito parcial de Masters (1982), família Rasch | `Rasch` |
| `"rsm"` | ordinais | escala de avaliação de Andrich (1978), família Rasch | `rsm` |
| `"gpcm"` | ordinais | crédito parcial generalizado de Muraki (1992) | `gpcm` |
| `"graded"` | ordinais | resposta graduada de Samejima (1969), o padrão | `graded` |

O curso ensina Rasch dicotômico e politômico no WinSteps (mapa de Wright, infit
e outfit, limiares); a skill roda os mesmos modelos no mirt por máxima
verossimilhança marginal, então os números não batem com o WinSteps (estimação
conjunta) casa a casa: comparar a leitura, não a terceira casa. Na dúvida entre
família Rasch e modelos com discriminação livre, `tri_comparar()` (seção 6).

## 3. `tri_ajustar()`

```r
r <- tri_ajustar(dados, itens, modelo = "graded", decimal = ",")
r$parametros      # a, b (ou limiares) com EP; rsm em forma de interceptos
r$ajuste_itens    # S-X² com p corrigido (BH); infit e outfit na família Rasch
r$informacao      # pico de informação e faixas de θ com precisão >= ,80 e >= ,90
relatar(r)
```

O que sai:

- **Ajuste global**: M2 para dicotômicos; M2* (Cai & Hansen, 2013) para
  politômicos; C2 quando o M2* não tem graus de liberdade (comum com poucos itens
  de muitas categorias). A tabela diz qual foi usada, com RMSEA (IC 90%), SRMSR,
  CFI e TLI. Família M2: Maydeu-Olivares e Joe (2006).
- **Ajuste por item**: S-X² de Orlando e Thissen (2000), com p corrigido por
  Benjamini-Hochberg (argumento `corrigir`). Na família Rasch, infit e outfit
  (MNSQ) com a régua de Linacre (2002), 0,5 a 1,5 (`regua_mnsq = "linacre"`), ou
  a régua por tamanho de amostra que o curso mostra (`"amostra"`).
- **Dependência local**: Q3 relativo de Christensen et al. (2017), Q3 menos a
  média dos Q3 fora da diagonal, com limite ,20 (`q3_limite`). Pares acima
  saem em `$pares_q3`.
- **Dimensionalidade na família Rasch**: autovalor do primeiro contraste da
  análise de componentes dos resíduos padronizados; 2 ou mais vira aviso.
- **Limiares desordenados** (PCM e GPCM): aviso por item. Não é desajuste: diz
  que uma categoria nunca é a mais provável; pensar em juntar categorias com
  base no conteúdo.
- **Ajuste de pessoas**: só agregado (porcentagem com Zh < -1,96, Drasgow et
  al., 1985; porcentagem com infit > 1,5). θ individual nunca é impresso.
- **Avisos de amostra**: N abaixo do mínimo citado no curso (Şahin & Anıl,
  2017): 150 na família Rasch, 200 nos demais, 500 no 3PL.
- **Discriminação**: a > 4 vira aviso (possível caso Heywood da TRI); a < ,40,
  item fraco (régua da skill, declarada como tal).

`ausentes = "fiml"` usa toda resposta disponível (o mirt estima com ausentes);
o S-X² não sai com ausentes, e a função avisa.

## 4. Réguas e o que o curso ensina

- MNSQ (infit, outfit): 0,5 a 1,5 produtivo para medida (Linacre, 2002). ZSTD é
  o teste de significância do MNSQ e depende do N: com N grande sai significativo
  para desvios pequenos, então o MNSQ vem primeiro e o ZSTD só desempata.
- Q3 relativo > ,20: dependência local a investigar (conteúdo repetido, item
  derivado de outro).
- Primeiro contraste dos resíduos ≥ 2: indício de segunda dimensão.
- Confiabilidade: a skill não fixa corte para a marginal; o uso pretendido
  decide, e a faixa de θ com precisão ≥ ,80 e ≥ ,90 mostra onde o escore serve.

## 5. Confiabilidade e precisão

Confiabilidade marginal (Green et al., 1984): σ²/(σ² + E[EP²(θ)]), com a
esperança tomada na distribuição latente estimada N(μ, σ²). Sai também a
empírica (variância dos EAP dividida pela variância dos EAP mais a média dos
EP²).

**Divergência com o mirt.** `mirt::marginal_rxx()` integra em N(0, 1) mesmo
quando a variância latente é estimada (família Rasch): infla com σ² < 1 e
encolhe com σ² > 1. Nos evals, com σ² = 2,39: mirt ,46; skill ,79; empírica
,80. A skill usa a fórmula com a variância estimada.

Precisão local: `$informacao` traz a faixa de θ (na métrica da variância
latente) em que a precisão σ²/(σ² + EP²) passa de ,80 e de ,90. Escala clínica
que informa só no extremo alto mede mal a população geral, e o contrário também.
A curva (`plot(r$modelo_mirt, type = "info")`) não tem dado de participante e
pode voltar para a conversa.

## 6. `tri_comparar()`

```r
tri_comparar(dados, itens, modelos = c("pcm", "gpcm", "graded"))
tri_comparar(dados, itens, modelos = c("rasch", "2PL", "3PL"))
```

Log-verossimilhança, número de parâmetros, AIC, BIC e BIC ajustado; razão de
verossimilhança só entre aninhados (Rasch ⊂ 2PL ⊂ 3PL; RSM ⊂ PCM ⊂ GPCM). GRM e
GPCM não são aninhados: comparar por critério de informação. Modelo mais
simples com ajuste aceitável e itens interpretáveis costuma ganhar; a família
Rasch tem propriedades de medida (escore total suficiente) que justificam
mantê-la com diferença pequena de ajuste.

## 7. DIF com `tri_dif()`

```r
tri_dif(dados, itens, "sexo", modelo = "graded", metodo = "mirt")
tri_dif(dados, itens, "sexo", metodo = "logistica", purificar = TRUE)
tri_dif(dados, itens, "sexo", metodo = "mh")   # v1.6, itens 0/1
```

**Rota mirt** (padrão): modelos de grupos múltiplos com média e variância
latentes livres no grupo focal. Sem `ancoras`, duas etapas: (1) a partir do
modelo todo restrito, libera um item por vez (`scheme = "drop"`) e marca os
candidatos pelo p corrigido; (2) usa os itens não marcados como âncora e testa
os candidatos acrescentando restrições (`"add"`). Com `ancoras`, testa os demais
contra elas. Na família Rasch sai o contraste de localização (focal menos
referência, em logits), lido pela classificação em logits que o curso usa no
WinSteps (abaixo de ,43 desprezível, ,43 a ,64 moderado, acima de ,64 grande);
âncoras ficam como "âncora (sem DIF)". RSM é
testado como PCM. Âncora contaminada enviesa o teste dos outros itens: relatar
como as âncoras foram escolhidas.

**Rota logística** (Swaminathan & Rogers, 1990; Zumbo, 1999): regressão
logística (dicotômicos) ou ordinal (`MASS::polr`) do item no escore de
pareamento, no grupo e na interação; testes uniforme, não uniforme e conjunto
(2 gl), com p corrigido. Escore purificado em rodadas (o item testado sempre
entra no próprio escore, como no difR). Tamanho: ΔR² de Nagelkerke pelas
classes de Jodoin e Gierl (2001), A < ,035, B ,035 a ,070, C > ,070, derivadas
para itens dicotômicos; para ordinais, o ΔR² de McFadden ≥ ,02 de Choi et al.
(2011) sai ao lado. DIF relevante = significativo e ΔR² ≥ ,035. Os evals
conferem a estatística e o ΔR² com `difR::difLogistic` (tol. 1e-8).

**Rota Mantel-Haenszel** (v1.6; `metodo = "mh"`, só itens dicotômicos; Martinková &
Hladká, 2023, 9.3.2): estratos pelo escore total (com o item), purificado por padrão;
αMH (chance de acerto da referência, primeiro nível do grupo, sobre a do focal no mesmo
escore), χ² de MH com correção de continuidade e ΔMH = −2,35 ln αMH, escala delta do
ETS, pela régua do difR: |Δ| < 1 A, 1 a 1,5 B, ≥ 1,5 C. DIF relevante = significativo e
|Δ| ≥ 1. Só enxerga DIF uniforme; DIF não uniforme pede a rota logística. Os evals
conferem χ² e αMH com `difR::difMH`, com e sem purificação (tol. 1e-8).

Tamanho importa mais que o p: com n grande, diferença irrelevante dá p
pequeno. Metade ou mais dos itens com DIF sugere diferença de traço mal
modelada ou âncoras contaminadas, não teste enviesado.

## 8. Relação com a AFC e a invariância

DIF na TRI e invariância de limiares, cargas e interceptos na AFC ordinal
respondem à mesma pergunta. Se o estudo já faz a sequência de
`references/invariancia.md`, o DIF localiza item a item; não relatar as duas
análises como evidências independentes. Viés preditivo (relação do escore com
um critério externo) é outra pergunta: `vies_preditivo()`.

## 9. TRI multidimensional (v1.5)

`tri_multidimensional(dados, itens, dimensoes, modelo)`, a rota das aulas "GRM
multidimensional (exploratório/confirmatório)" e "PCM multidimensional" do curso
R: TRI.
- `dimensoes = 3`: exploratório, rotação oblimin (padrão do mirt e das aulas),
  cargas rotacionadas e correlações entre dimensões. Só GRM, GPCM ou 2PL: Rasch
  exploratório não tem discriminação livre para rotacionar.
- `dimensoes = list(DEP = c(...), ANS = c(...))`: confirmatório com dimensões
  correlacionadas (`COV`), inclusive PCM (Rasch) multidimensional, com infit e outfit.
- Estimação EM com até 2 dimensões; MHRM com 3 ou mais, como nas aulas (semente fixa).
- Saem: ajuste global (M2, M2* ou C2), cargas e h², correlações entre dimensões,
  parâmetros na parametrização inclinação-intercepto (a1..ak, d; o mirt não converte
  para a, b com mais de uma dimensão), ajuste dos itens com p corrigido, confiabilidade
  empírica por dimensão (θ EAP calculado e descartado) e comparação com o
  unidimensional (AIC, BIC, SABIC e razão de verossimilhança).
- Bifator da TRI (`bfactor`, aula do PCM) fica fora: o bifator pela AFC está em
  `indices_bifator()`.
- Os evals recuperam a estrutura simulada (cada item na dimensão certa, correlação ,40
  entre dimensões) e conferem a confiabilidade empírica com o quadrado da correlação
  entre θ verdadeiro e EAP.

## 9b. TRI não paramétrica de Mokken (v1.7)

`escala_mokken(dados, itens, limite = .30, aisp = TRUE, monotonicidade = TRUE, iio =
FALSE)`, pelo capítulo 15 do Wiley Handbook (Meijer & Tendeiro, 2018) e pelo pacote
mokken (van der Ark, 2012): H da escala com IC, Hi e Hij por item com EP, a seleção
automática de itens (AISP) com limite c, a monotonicidade manifesta (crit > 80 é
violação séria) e, com `iio = TRUE`, a ordenação invariante dos itens (HT ≥ ,30). Régua
de Sijtsma e Molenaar (2002), citada no capítulo: H < ,30 não escalável, ,30 fraca, ,40
média, ,50 forte. Ressalva do capítulo, que a função repete no aviso: item que se ajusta
ao Rasch com discriminação perto de 1 pode ter Hi < ,30, e a AISP descarta itens que um
modelo paramétrico aceitaria; Mokken complementa a TRI paramétrica e a AFC, não decide
sozinho. Escores precisam ser inteiros; o deslocamento é o mesmo para todos os itens
(recodificar item a item mudaria H). H, Hi, EP e a partição da AISP iguais ao mokken
nos evals; dados de Guttman perfeitos dão H = 1. Sequência da análise em Sijtsma e van
der Ark (2017).

## 9c. Simulação pós-hoc de CAT (v1.7)

`simular_cat(dados, itens, modelo = "2PL", ep_max = .30, max_itens = NULL)`, como em
Martinková e Hladká (2023, 10.2.6) e Petersen (2024, cap. 20): calibra o modelo (ou usa
`fit =`), refaz o teste de cada pessoa como se fosse adaptativo com as respostas que ela
deu (mirtCAT; primeiro item de máxima informação, θ por MAP, seleção por máxima
informação, parada com EP ≤ ,30, que corresponde a confiabilidade perto de ,91) e
devolve só agregados: itens usados (média, mediana, extremos), % de itens economizados,
% que atingiu o EP, correlação, RMSE e viés entre o θ do CAT e o do teste completo, e a
exposição de cada item. Nenhum θ por pessoa fica no objeto. Pós-hoc superestima a
concordância com o teste inteiro (as respostas são as mesmas) e não mede efeitos de
ordem; banco com pouca informação nos extremos deixa gente sem atingir o EP (aviso).
Nos evals, o número médio de itens é igual ao do mirtCAT direto. Serve para planejar
forma curta ou banco adaptativo, não para aplicar CAT a pacientes.

## 10. Fora das funções

O módulo 05 da pasta Psicometria-R (Cook & Pitoniak, 2026; Wiley Handbook, cap.
27) cobre ainda equacionamento, só em texto de livro, sem aula com código. A
skill não tem função validada para ele. Caminho, sem conferência nos evals:
equacionamento com o pacote `equate`
(TCT) ou `plink`/`SNSequate` (TRI). Usar só com a documentação do pacote aberta e
dizer no relato que a rota saiu do escopo testado da skill. Fixar parâmetros de
itens (ancoragem de calibração), como no script "Exemplo para fixar parâmetros no
mirt" das aulas, também fica fora: `mirt(..., pars = "values")`, editar `value` e
`est`, e rodar de novo.
