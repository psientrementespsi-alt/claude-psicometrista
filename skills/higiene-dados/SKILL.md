---
name: "higiene-dados"
description: "Sanitização de planilhas brutas, validação declarativa de regras de dados, triagem de respostas descuidadas e insuficientes (C/IER: longstring, IRV, psychsyn/psychant, tempo de resposta, ajuste de pessoa G^p e lz^p), detecção robusta de atípicos (MAD univariado, FastMCD multivariado), diagnóstico e tratamento de dados ausentes (teste de Little, MICE com variáveis auxiliares, missForest) e geração de codebook e fluxograma de exclusão (JARS-Quant). Prepara o banco purificado para a psicometria-r."
metadata:
  version: '1.5.0'
---

# SKILL: Higiene e Tratamento de Dados (v1.5.0)

Executa a sanitização técnica, a validação lógica, a triagem de respostas descuidadas, a detecção de atípicos e o tratamento de dados ausentes em bancos de pesquisa em psicologia e ciências do comportamento. Transforma uma planilha bruta de coleta em um banco analítico auditado e purificado, gerando o log completo de descartes exigido pelas diretrizes APA 7 e JARS-Quant antes de qualquer modelagem na `psicometria-r`.

<required_reading>
0. Se uma skill `regras-comuns` estiver instalada, leia-a antes de começar; se não estiver, siga as regras desta skill.
1. Antes de executar qualquer etapa, leia o arquivo correspondente em `references/`:
   - `references/01_sintese_higiene_dados.md` (fundamentos teóricos e metodológicos)
   - `references/04_guia_deteccao_cier_e_outliers.md` (fórmulas e pontos de corte para C/IER e atípicos)
   - `references/05_guia_dados_ausentes_e_mice.md` (diagnóstico MCAR/MAR e imputação MICE)
</required_reading>

## IRON RULES (desta skill, somadas às do hub)

1. **Dado de participante não sai do computador.** A higiene roda estritamente no R local da máquina do usuário. Nunca copiar planilha bruta para a nuvem, nunca colar linhas individuais de participantes na conversa e nunca exibir identificadores pessoais (IDs de e-mail, nomes, IPs).
2. **Imutabilidade da base bruta.** O arquivo original de coleta é intocável (`read-only`). O script de higiene gera um novo arquivo serializado (`dados_purificados.rds`) e um registro detalhado de auditoria (`log_descartes.rds`).
3. **Trilha de auditoria reprodutível (JARS-Quant).** Todo descarte de participante ou valor precisa ser justificado por um critério explícito a priori, registrado com contagem exata no fluxograma de inclusão e exclusão. Descarte ad hoc sem critério matemático documentado é terminantemente proibido.
4. **Análise de sensibilidade para atípicos (Regra de Aguinis).** Participantes classificados como atípicos extremos não são deletados silenciosamente. O banco purificado recebe uma flag indicadora (`flag_outlier`), permitindo que a `psicometria-r` execute e relate a modelagem com e sem os atípicos.
5. **Nenhum pacote ausente interrompe o pipeline.** Os scripts da skill contam com implementações puras e vetorizadas em R para as métricas centrais (MAD, MCD, IRV, longstring, psychsyn, even-odd, Mahalanobis, Little MCAR, regras, codebook), integrando-se nativamente a pacotes especializados (`janitor`, `validate`, `careless`, `naniar`, `mice`, `robustbase`) quando disponíveis.

---

## Economia de tokens

1. Scripts de `scripts/` são executados com `source()`, nunca lidos como contexto. Para a assinatura de uma função, use `args(nome)` no R ou `grep -n "nome <- function" scripts/*.R`.
2. Leia só a referência do módulo do pedido. Referência longa abre com sumário e números de linha: leia o sumário e depois só a seção necessária.
3. Não releia arquivo já lido nesta conversa; reaproveite o que está no contexto.
4. Saída do R enxuta: resumos e agregados, nunca objetos inteiros (o que também protege dado de participante).
5. Todo cálculo, conferência ou tabela determinística sai de função da skill; não refazer conta à mão em tokens.
6. Modo rápido, para exploração e iteração: só números e decisões, sem a verificação narrativa completa. Na entrega final, a verificação completa e o checklist da skill são obrigatórios.
7. Modelo: execução mecânica (rodar script, formatar tabela) cabe num modelo leve; interpretação e decisão metodológica pedem o modelo mais forte. Em lote via API, Batch API (50% de desconto) e cache dos materiais relidos.

## O Pipeline de Higiene em 5 Fases

```
   [Planilha Bruta: .csv / .xlsx / .sav]
                     │
                     ▼
  ┌─────────────────────────────────────┐
  │ FASE 1: Sanitização Estrutural     │ -> janitor, tipos, duplicatas
  └──────────────────┬──────────────────┘
                     │
                     ▼
  ┌─────────────────────────────────────┐
  │ FASE 2: Validação Declarativa       │ -> regras, critérios de inclusão, ranges
  └──────────────────┬──────────────────┘
                     │
                     ▼
  ┌─────────────────────────────────────┐
  │ FASE 3: Triagem C/IER (Careless)    │ -> famílias: repetição, tempo, psychsyn, even-odd, Mahalanobis, atenção, ajuste de pessoa
  └──────────────────┬──────────────────┘
                     │
                     ▼
  ┌─────────────────────────────────────┐
  │ FASE 4: Atípicos Robustos (Outliers)│ -> MAD univariado, FastMCD (MCD75)
  └──────────────────┬──────────────────┘
                     │
                     ▼
  ┌─────────────────────────────────────┐
  │ FASE 5: Diagnóstico e Imputação     │ -> Little MCAR, MICE, lavaan.mi, codebook
  └──────────────────┬──────────────────┘
                     │
                     ▼
   [dados_purificados.rds + log_descartes.rds] ──> Para a psicometria-r
```

### Pipeline v2 (`executar_pipeline_higiene()`, 1.3.0 e 1.4.0)
Uma chamada roda, nesta ordem, cada etapa registrada na trilha de auditoria (`registrar_passo()`): registros vazios; recodificação de rótulos e valores fora da escala ou não inteiros para NA (contados em valores); regras e critérios de inclusão; registros sem resposta aos itens; duplicatas (só entre quem respondeu; protocolo todo no mínimo nunca conta, e com menos de três colunas além dos itens o pipeline avisa que padrões iguais podem ser coincidência); omissão acima de `limiar_omissao_caso`; Little; C/IER por famílias; FastMCD75 sobre escores de subescala (`fatores`, `mcd_em = "subescalas"`) ou itens, só assinalando; MICE opcional. Devolve `dados_purificados`, `trilha`, `tabela1_md` (Tabela 1 do fluxo em Markdown), `fluxograma_mermaid` (diagrama de fluxo de participantes em Mermaid, no formato JARS-Quant/CONSORT, gerado da mesma trilha por `fluxograma_mermaid()`), `texto_metodo` (parágrafo de Método em português com os números da trilha), `regras`, `codebook` (`gerar_codebook()`: só agregados, categoria com menos de 5 casos como "<5", texto livre omitido), `metricas_cier`, `teste_mcar` e `mcd`. Nada devolve linha de participante.
Para o limite de 60 s do r-btw: `salvar_em` grava `estado_<etapa>.rds` a cada etapa, `etapas_ate` para numa etapa e `retomar_de` continua com os parâmetros da primeira execução (escala, cortes e semente não mudam no meio). O eval refaz o banco DASS-21 do projeto e reproduz N = 116 e 5 atípicos MCD75 da fase 1 de 29/10/2026.

### Fase 1: Sanitização Estrutural
- Padronização de nomes de variáveis em `snake_case` com `janitor::clean_names()`.
- Remoção de linhas de metadados de plataformas online (ex.: cabeçalhos duplicados do Qualtrics).
- Remoção de linhas e colunas 100% vazias (`janitor::remove_empty()`).
- Conversão explícita de tipos: variáveis categóricas para `factor`, itens psicométricos para `numeric` ou `ordered`, datas para `Date`.
- Identificação e isolamento de duplicatas de resposta (mesmo IP, mesmo ID ou respostas idênticas em toda a extensão do questionário).

### Fase 2: Validação Declarativa de Regras
- Definição formal de limites admissíveis para cada escala (ex.: itens Likert de 1 a 5).
- Verificação de consistência lógica (ex.: idade de universitário >= 16; tempo de término posterior ao de início).
- Confrontação de regras com `validar_regras(df, regras)` (em `scripts/auditoria_higiene.R`): regras escritas como expressões do R (`c(adulto = "idade >= 18")`), avaliadas pelo `validate` quando instalado e em R base quando não; devolve só contagens (avaliados, violações, não avaliáveis). No pipeline, `regras_excluir` nomeia as regras que são critério de inclusão (consentimento, idade) e tiram o caso; as demais só são relatadas.
- Valores fora de escopo são recodificados para `NA` ou sinalizados para inspeção.

### Fase 3: Triagem de Respostas Descuidadas (C/IER)
Conforme as recomendações de Curran (2016) e Ward & Meade (2023), adota-se múltiplos indicadores:
1. **Longstring:** Contagem da maior sequência de respostas consecutivas idênticas. Corte padrão: metade da escala (Curran, 2016, regra de base declarada pelo autor como ponto de partida). Itens invertidos quebram sequências legítimas; em inventário clínico unipolar, `ignorar_todos_zeros = TRUE` protege quem não tem sintoma.
2. **Variância Intraindividual de Resposta (IRV):** Identifica variabilidade nula (straightlining, $IRV = 0$). Como IRV = 0 implica longstring igual ao número de itens, longstring e IRV formam uma família só (repetição) e nunca somam duas infrações.
3. **Sinônimos e Antônimos Psicométricos:** Correlação intraindividual em pares de itens com correlação amostral $|r| > .60$; correlação intrapessoal negativa nos sinônimos é o corte conservador de Curran (2016). **Even-odd** (`calc_evenodd()`, igual ao `careless::evenodd` nos evals e aceitando subescalas não contíguas): Spearman-Brown entre médias pares e ímpares das subescalas, menor que zero como corte; com menos de 4 subescalas a correlação intrapessoal é instável e a função avisa. Índice que assinala mais de 20% dos casos calculáveis não discrimina na amostra e sai do critério com aviso (heurística da skill, `limite_informativo`).
4. **Piso Cognitivo de Tempo de Resposta:** Tempo por item abaixo de 2 segundos (Huang et al., 2012, citado por Curran, 2016, como corte conservador). Sem início e fim registrados, a família tempo não é calculada e isso é declarado.
5. **Mahalanobis como índice de C/IER** (`calc_mahad_cier()`, `usar_mahad = TRUE`): D² clássico nos itens contra χ² com k gl (Meade & Craig, 2012; DeSimone et al., 2015). DeSimone et al. sugerem os 5% superiores e lembram que a distribuição χ² supõe itens normais; a skill usa p < .001, mais conservador com itens ordinais. Com `fatores`, um D² por subescala e o menor p com Bonferroni, como Meade e Craig fizeram com 60 itens por fator. Pega resposta aleatória ou alternada; retilíneo perto da média não sai, e é a família repetição que o pega. É distinto do FastMCD da Fase 4, que procura atípico substantivo e só assinala.
6. **Itens de atenção** (`calc_falhas_atencao()`, `itens_atencao`): itens de instrução direta (resposta esperada) e de infrequência (respostas que contam como falha; Huang et al., 2015). Mais de `max_falhas_atencao` falhas forma uma família.
7. **Ajuste de pessoa** (`calc_ajuste_pessoa()` em `scripts/ajuste_pessoa.R`, passado à triagem em `flag_ajuste_pessoa`): erros de Guttman politômicos G^p, com passos de item ordenados pela popularidade condicional do Apêndice A, e lz^p sob o modelo de resposta graduada (mirt), somados nas subescalas (extensão multiteste; Niessen, Meijer & Tendeiro, 2016). Serve para questionários de tamanho médio, em que even-odd e psychsyn têm poucos pares. O artigo usa corte empírico (5%, 10% ou 20% mais extremos), então a família assinala cerca de 5% por construção e só exclui junto com outra. G^p é o padrão por não exigir modelo; lz^p exige o GRM ajustado e itens sem ausentes. Escala com menos de 5 itens gera aviso de baixa sensibilidade.
8. **Critério de Exclusão:** Excluir apenas participantes com infração em pelo menos duas famílias distintas (repetição, tempo, psychsyn, even-odd, Mahalanobis, atenção, ajuste de pessoa; Ward & Meade, 2023). `triagem_cier_completa()` devolve as famílias, as decisões em `attr(, "decisoes")` e os avisos em `attr(, "avisos")`.

### Fase 4: Detecção Robusta de Atípicos (Outliers)
Conforme demonstrado por Leys et al. (2013, 2019):
- **Proibição:** A regra da média $\pm$ 2 ou 3 desvios-padrão é terminantemente proibida, pois sofre do efeito de mascaramento (*masking*) e inflação da própria variância.
- **Outliers Univariados:** Utilização do Desvio Absoluto Mediano (MAD), com constante de escala $b = 1.4826$:
  $$MAD = 1.4826 \cdot \text{mediana}(|X_i - \text{mediana}(X)|)$$
  Corte conservador: $|X_i - \text{mediana}(X)| / MAD > 3.0$ (ou 2.5 sob suspeita prévia).
- **Outliers Multivariados:** Distância de Mahalanobis robusta com o Minimum Covariance Determinant (FastMCD) reponderado, padrão MCD75 (ponto de ruptura .25; Leys et al., 2019), com χ² a p = .001; MCD50 só com `alpha_mcd = 0.5`, quando se suspeita de mais de 25% de atípicos. Com matriz singular (itens com mais da metade das respostas no mesmo valor) a grade sobe para .85 e .90. Quando os atípicos vão definir a análise de sensibilidade de uma regressão, não incluir o desfecho no espaço do MCD: a exclusão passa a selecionar pelo desfecho e infla o R² por construção; a rota da `psicometria-r` (Cook e Mahalanobis nos preditores) responde a outra pergunta e sai ao lado.
- **Tratamento:** Não descartar imediatamente. Marcar na coluna `flag_outlier` e planejar análise de sensibilidade.

### Fase 5: Diagnóstico e Tratamento de Dados Ausentes
Conforme Schafer & Graham (2002) e van Buuren (2018):
- **Diagnóstico:** Inspeção de padrões de omissão e cálculo da taxa de ausência por item e por participante. Participantes com $> 20\%$ de omissão na escala principal são candidatos a exclusão do modelo substantivo.
- **Mecanismos:**
  - Teste de Little para Missing Completely at Random (MCAR).
  - Se $p > .05$, o teste não rejeita MCAR, o que não prova MCAR nem autoriza sozinho a deleção pairwise; com poucas omissões o poder é muito baixo (psicometria-r, `afc-ordinal.md` e `fontes-e-decisoes.md`).
  - Se $p \le .05$, MCAR fica afastado sem dizer se é MAR ou MNAR; trabalhar com MAR condicional às variáveis demográficas e auxiliares e discutir MNAR.
- **Tratamento:**
  - Proibição de imputação pela média ou dedutiva simples.
  - Para modelagem contínua ou regressão: Imputação Múltipla por Equações Encadeadas (MICE) com $m = 20$ a $50$ imputações, incorporando variáveis auxiliares preditoras de omissão.
  - Para levar a imputação ao lavaan: `exportar_lavaan_mi(imp)` converte o objeto `mids` numa lista de bancos imputados pronta para `semTools::runMI()` / `lavaan.mi`, com pooling pelas regras de Rubin (1987) e a fração de informação ausente (FMI) relatada.
  - Para AFC categórica na `psicometria-r`: se MAR for sustentável, encaminhar para FIML/MLR ou MICE; se WLSMV for usado, alertar que a deleção pairwise assume MCAR estrito.

---

## Estrutura de Arquivos da Skill

Os scripts ficam no Mac em `~/.gemini/config/skills/higiene-dados/` (cópias em `~/.claude/skills/` e `~/skills/`); a versão completa também está na página da skill no Notion.

- `SKILL.md`: Instruções operacionais centrais e catálogo de regras.
- `references/`:
  - `01_sintese_higiene_dados.md`: Síntese das 30 obras fundamentais de higiene e tratamento de dados.
  - `02_taxonomia_credit_e_fontes.md`: Matriz de contribuição científica e referências completas com DOI.
  - `03_cruzamento_psicometria_e_jars.md`: Alinhamento metodológico com a skill de psicometria e JARS-Quant.
  - `04_guia_deteccao_cier_e_outliers.md`: Fórmulas matemáticas, algoritmos e pontos de corte para C/IER e atípicos.
  - `05_guia_dados_ausentes_e_mice.md`: Protocolo de diagnóstico de ausentes e rotinas MICE.
  - `06_acervo_local_e_uso.md`: Índice dos papers e repositórios da pasta `psicomtria rrr` e qual decisão da skill cada um sustenta.
- `scripts/`:
  - `pipeline_higiene.R`: Script mestre para execução automatizada e geração do relatório de purificação.
  - `screening_cier.R`: Longstring, IRV, psychsyn/psychant, even-odd, Mahalanobis como C/IER, itens de atenção, triagem por famílias e Mahalanobis robusto (FastMCD).
  - `imputacao_mice.R`: Funções para diagnóstico de ausentes, imputação encadeada e `exportar_lavaan_mi()`.
  - `ajuste_pessoa.R`: `calc_guttman_poly()`, `calc_lz_poly()` e `calc_ajuste_pessoa()` (Niessen et al., 2016).
  - `validacao_codebook.R`: desativado na 1.5.0; só carrega `auditoria_higiene.R` (a versão antiga sobrescrevia `validar_regras()` com outra assinatura).
  - `auditoria_higiene.R`: `validar_regras()`, `gerar_codebook()`, `registrar_passo()`, `tabela_fluxo_md()`, `fluxograma_mermaid()`, `texto_metodo_higiene()`.
  - `evals_higiene.R`: 37 testes (na 1.5, G^p contra o Apêndice A e o exemplo dicotômico de Niessen et al., lz^p igual ao `mirt::personfit()`, família de ajuste na triagem; antes, fluxograma Mermaid e exportação para o lavaan.mi, além de recodificação, Little, duplicatas, even-odd contra o careless, triagem por famílias, MCD75, Mahalanobis e atenção, regras, codebook, pipeline v2 encadeado e em blocos, reprodução do banco DASS-21).
- `templates/`:
  - `template_tabela_descarte_apa7.md`: Template de tabela de fluxo de participantes para publicação.
