# Pacote de resultados e relato APA 7

## Sumário

- 8 (linha 16)
- 38 (linha 46)
- 91 (linha 99)
- 105 (linha 113)
- 112 (linha 120)

Índice: Conteúdo do pacote · Frases de resultado · Diferenças em relação ao modo PSICOMETRICO da fabrica-de-artigo · Referências de software · O que não entra no pacote

O pacote é o que esta skill entrega para a fabrica-de-artigo (modo
PSICOMETRICO). Todo número sai do R da sessão; nada é digitado à mão.

## Conteúdo do pacote

1. **Análise de dados (para o Método).** Software e versões
   (`versoes_pacotes()`), estimador com a justificativa devolvida por
   `sugerir_estimador()`, tratamento de ausentes, colapso de categorias (se
   houve), versão dos índices de ajuste (robusta ou escalonada, e por quê),
   critérios de invariância, número de reamostragens e sementes.
2. **Tabela de ajuste por modelo** (`tabela_ajuste()` de cada modelo
   concorrente, só a versão declarada): modelo, χ², gl, p, CFI, TLI, RMSEA
   [IC 90%], SRMR.
3. **Tabela de cargas** (`lavaan::standardizedSolution()`): carga padronizada,
   erro-padrão, IC 95%, por fator.
4. **Correlações entre fatores, AVE e HTMT** (`validade_discriminante()`).
5. **Confiabilidade** (`confiabilidade()` por fator; `confiabilidade_compostos()` para o composto de segunda ordem, se houver): ômega do escore observado e alfa
   por fator, mais composto de segunda ordem se houver.
6. **Invariância** (`invariancia()$tabela`): um passo por linha, com Δ e p.
7. **TRI, rede, EGA**, quando entrarem: parâmetros dos itens, ajuste (M2* ou
   C2, S-X²), DIF; rede pelo `$checklist` de `rede_psicometrica()` (estimador e γ,
   ausentes e tamanho amostral da estimação, arestas não nulas sobre as possíveis,
   arestas mais fortes, centralidade em z) mais o CS por índice e os IC de arestas
   de `estabilidade_rede()`, com a frase de alcance causal; bootEGA com mediana de
   dimensões e estabilidade.
8. **Análises estatísticas** (diferenças, associação, processos, modelos, poder):
   de cada `pr_resultado`, a frase `$apa` (ou `relatar(res)`), a tabela `$tabela`,
   o índice principal de `$efeito` com IC, os pressupostos lidos e as decisões de
   rota de `$decisoes`. Cálculo amostral com a tela do G*Power para replicar
   (`poder_*()`).
9. **Decisões do autor pendentes**, marcadas `[DECISÃO DO AUTOR]` (o marcador
   que a fabrica-de-artigo já usa): colapso, modelo final, itens retirados.

## Frases de resultado

`formatar_ajuste_apa(linha)` monta a frase de ajuste a partir de uma linha de
`tabela_ajuste()`:

```r
taj <- tabela_ajuste(fit, "três fatores")
v <- if ("robusta" %in% taj$versao) "robusta" else "escalonada"  # robusta pode sair NA
formatar_ajuste_apa(taj[taj$versao == v, ])
```

Saída nos evals da skill (dados simulados, modelo de três fatores correto):
`χ²(51) = 66.87, p = .067, CFI = .993, TLI = .991, RMSEA = .021, IC 90%
[.000, .035], SRMR = .026`. Serve só para mostrar o formato; os números de um
relato vêm da execução com os dados do estudo.

Regras de formato aplicadas pela função, coerentes com
`fabrica-de-artigo/references/statistical_reporting_standards.md`:

- sem zero à esquerda em p, CFI, TLI, RMSEA e SRMR (os que não passam de 1, e
  o TLI pela mesma convenção de relato do CFI, embora ele possa passar de 1);
  χ² com zero;
- p com três casas, `p < .001` abaixo disso e `p > .999` acima;
- `decimal = ","` troca o separador decimal para manuscrito em português e
  passa a separar os itens da frase (e os limites do IC) por ponto e vírgula,
  para a vírgula não ter dois papéis; o uso de zero à esquerda com vírgula
  varia entre periódicos brasileiros (`zero_inicial = TRUE` liga), então
  conferir a norma do periódico-alvo.

Frase do estimador, a adaptar com o que `sugerir_estimador()` devolveu:

> Como os itens têm quatro categorias de resposta, a AFC usou o estimador
> WLSMV sobre correlações policóricas (Rhemtulla et al., 2012). Os índices de
> ajuste relatados são as versões robustas calculadas pelo lavaan (versão X).

Frase de ressalva dos cortes com WLSMV:

> Os pontos de corte convencionais (Hu & Bentler, 1999) foram derivados para
> estimação por máxima verossimilhança; com DWLS, RMSEA tende a ser menor e
> CFI e TLI maiores para o mesmo grau de desajuste (Xia & Yang, 2019), e por
> isso os valores são interpretados como indicativos.

Funções da v1.4 com frase pronta em `$apa` (`relatar(res)`): `afe_numero_fatores()`
(paralela, Hull e MAP numa frase), `analise_classica_itens()`,
`confiabilidade_teste_reteste()`, `concordancia_juizes()`, `confiabilidade_ema()`,
`indices_bifator()`, `tri_ajustar()`, `tri_comparar()`, `tri_dif()`,
`poder_afc_rmsea()`, `n_aipe_rmsea()`, `poder_afc_aninhados()`, `n_aipe_parametro()`,
`curva_roc()`, `vies_preditivo()`, `controle_aquiescencia()`, `classes_latentes()`,
`perfis_latentes()`, `razao_prevalencia()`, `sobrevida()`, `meta_proporcao()` e
`meta_acuracia()`. As frases de poder para AFC deixam a origem dos valores
populacionais como "[preencher]", como as do G*Power; as de bifator e TRI trazem a
fonte da fórmula. Réguas citadas nas frases são as declaradas em `$decisoes`.

## Diferenças em relação ao modo PSICOMETRICO da fabrica-de-artigo

O modo PSICOMETRICO foi alinhado com esta skill em 11/10/2026, mas três pontos
de vocabulário precisam passar certos para o texto:

- invariância ordinal: "cargas" corresponde à invariância métrica e
  "interceptos" à escalar; não há passo estrito na sequência ordinal daqui;
- ômega: dizer qual (ômega do escore observado de Green e Yang, 2009, o padrão
  do `compRelSEM` com itens ordinais; ômega hierárquico ou total do composto
  global, conforme o modelo), em vez de só "ômega de McDonald";
- PBAT: banco de itens de uso idionômico, sem alfa, ômega, AVE nem AFC de
  validação de escore (ver `references/pbat-mpfi.md`), o que o checklist de
  escala do modo PSICOMETRICO não cobre.

## Referências de software

Sair de `citation("lavaan")`, `citation("semTools")` etc. no R da sessão, que
dá a referência que os autores do pacote pedem, com a versão instalada. Para o
lavaan, a referência do artigo é Rosseel (2012); para o mirt, Chalmers (2012);
para o bootnet, Epskamp, Borsboom e Fried (2018). R em si: `citation()`.

## O que não entra no pacote

Linha de participante, ID, tabela cruzada com célula menor que 5, saída
completa de `summary()` com dados descritivos individuais. Figura de rede e
curva de informação podem entrar, porque não carregam dado individual.
