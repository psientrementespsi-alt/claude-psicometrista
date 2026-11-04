---
name: "relato-estatistico"
description: "Redação e formatação de resultados estatísticos e psicométricos em psicologia no padrão APA 7 e JARS-Quant. Converte saídas e modelos do R (lavaan, mirt, bootnet, lm, glm, PROCESS) em prosa acadêmica fluida e tabelas de publicação prontas para Word, PDF ou Markdown, com verificação automática de consistência interna via statcheck."
metadata:
  version: '1.5.0'
---

# SKILL: Relato Estatístico APA 7 e JARS-Quant (v1.5.0)

Traduz saídas numéricas, modelos estatísticos e análises psicométricas rodadas no R para texto acadêmico e tabelas rigorosamente aderentes às normas da 7ª edição do Manual de Publicação da American Psychological Association (APA 7) e às diretrizes do Journal Article Reporting Standards for Quantitative Research (JARS-Quant; Appelbaum et al., 2018). Garante a consistência interna das estatísticas inferenciais através do protocolo `statcheck` (Nuijten et al., 2016).

<required_reading>
0. Se uma skill `regras-comuns` estiver instalada, leia-a antes de começar; se não estiver, siga as regras desta skill.
1. Antes de redigir qualquer seção ou tabela, consulte o módulo correspondente em `references/`:
   - `references/01_manual_relato_apa7_jars.md` (normas completas de estilo e JARS-Quant)
   - `references/03_blueprints_fraseado_apa7.md` (estruturas sintáticas para cada família de modelos)
   - `references/04_diretrizes_tabelas_e_figuras.md` (normas estritas de diagramação de tabelas)
   - `references/05_auditoria_consistencia_statcheck.md` (protocolo de conferência de NHST)
   - `references/06_acervo_local_e_uso.md` (repositórios e papers locais; convenções que divergem entre as skills)
</required_reading>

## IRON RULES (desta skill, somadas às do hub)

1. **Zero inventado: todo número sai do R.** É terminantemente proibido preencher graus de liberdade, tamanhos de efeito ou intervalos de confiança de cabeça ou por aproximação. O que não foi computado na sessão é declarado como não calculado.
2. **Regras de Ouro do Zero à Esquerda (APA 7):**
   - **Sem zero antes da vírgula/ponto** para grandezas matematicamente limitadas a 1:
     - Valores de significância ($p < .001$, $p = .042$).
     - Coeficientes de correlação ($r = .45$, $\rho = .38$).
     - Cargas fatoriais e coeficientes padronizados ($\lambda = .78$, $\beta = .32$).
     - Índices de confiabilidade ($\alpha = .87$, $\omega = .89$).
     - Índices de ajuste com teto unitário ($\text{CFI} = .965$, $\text{TLI} = .952$, $R^2 = .24$).
   - **Com zero antes da vírgula/ponto** para grandezas que podem ultrapassar 1:
     - Médias e desvios-padrão ($M = 0.45$, $SD = 1.20$).
     - Estatísticas de teste ($t = 2.45$, $F = 4.12$, $\chi^2 = 15.30$, $z = 1.96$).
     - Tamanhos de efeito não limitados ($d = 0.52$, $g = 0.48$).
     - Índices de resíduo e erro quadrático ($\text{RMSEA} = 0.048$, $\text{SRMR} = 0.035$).
   - Esta é a convenção padrão (perfil `apa7_relato`). Para o periódico-alvo, o perfil do Modo 5 decide; divergências com a `psicometria-r` (RMSEA e SRMR sem zero em `relato-apa.md`; β com zero em `associacao.md`) são decisão do autor, resolvidas pela escolha do perfil e marcadas `[DECISÃO DO AUTOR]` no pacote (`references/06_acervo_local_e_uso.md`).
3. **Tipografia e Itálico Estatístico:**
   - Letras latinas que representam estatísticas ou símbolos matemáticos são em itálico: *$M$, $SD$, $t$, $F$, $p$, $r$, $d$, $z$, $N$, $n$*.
   - Letras gregas **nunca** são em itálico: $\alpha, \beta, \omega, \chi^2, \eta_p^2, \lambda$.
   - Siglas multiletras **nunca** são em itálico: CFI, TLI, RMSEA, SRMR, AIC, BIC, SE, IC.
4. **Proibição do Teste de Sobel em Mediação:**
   - O teste de normalidade do produto de Sobel é terminantemente proibido devido à assimetria inerente da distribuição do efeito indireto ($ab$). O efeito indireto deve ser obrigatoriamente relatado com Intervalo de Confiança obtido por bootstrap (5.000 reamostragens, percentil ou BCa).
5. **Diagramação de Tabelas APA 7:**
   - Exatamente 3 linhas horizontais visíveis: no topo da tabela, abaixo da linha de cabeçalhos e no rodapé inferior (antes da nota).
   - Zero linhas verticais (proibição absoluta de separadores de colunas verticais).
   - Todo número alinhado pelo ponto decimal.
6. **Auditoria Statcheck Obrigatória:**
   - Qualquer frase que reporte testes de hipótese nula (*NHST*: $t, F, \chi^2, r, z$) acompanhados de $p$ deve ser conferida pela regra do statcheck (seção Modo 4). Discrepâncias grosseiras (ex.: relatar $p < .05$ quando a estatística resulta em $p = .12$) barram a entrega do texto.
7. **Ômega com itens ordinais:** o que vai para o relato é o ômega do escore observado (Green & Yang, 2009) com IC por bootstrap; `omega_delta` da `psicometria-r` descreve a variável-resposta latente e só entra rotulado como tal.

---

## Economia de tokens

1. Scripts de `scripts/` são executados com `source()`, nunca lidos como contexto. Para a assinatura de uma função, use `args(nome)` no R ou `grep -n "nome <- function" scripts/*.R`.
2. Leia só a referência do módulo do pedido. Referência longa abre com sumário e números de linha: leia o sumário e depois só a seção necessária.
3. Não releia arquivo já lido nesta conversa; reaproveite o que está no contexto.
4. Saída do R enxuta: resumos e agregados, nunca objetos inteiros (o que também protege dado de participante).
5. Todo cálculo, conferência ou tabela determinística sai de função da skill; não refazer conta à mão em tokens.
6. Modo rápido, para exploração e iteração: só números e decisões, sem a verificação narrativa completa. Na entrega final, a verificação completa e o checklist da skill são obrigatórios.
7. Modelo: execução mecânica (rodar script, formatar tabela) cabe num modelo leve; interpretação e decisão metodológica pedem o modelo mais forte. Em lote via API, Batch API (50% de desconto) e cache dos materiais relidos.

## Os 5 Modos Operacionais da Skill

### Modo 1: Descritivas e Comparações de Grupos
- **Testes parametrizados e robustos:** Teste *t* de Student para amostras independentes e pareadas, Teste *t* de Welch (correção para variâncias heterogêneas), Teste U de Mann-Whitney / Wilcoxon.
- **ANOVA e Comparações Múltiplas:** ANOVA unidirecional (*One-Way*) e fatorial, correção de Welch e Brown-Forsythe, pós-hoc de Tukey HSD ou Games-Howell.
- **Tamanhos de Efeito Obrigatórios:** *d* de Cohen ou *g* de Hedges (com IC 95%), $\eta^2$ ou $\eta_p^2$ para ANOVA, $V$ de Cramér para tabelas de contingência.
- **Formato canônico:**
  > O grupo intervenção apresentou escores significativamente superiores ($M = 4.82$, $SD = 0.65$) em comparação ao grupo controle ($M = 3.91$, $SD = 0.81$), $t(98) = 6.18$, $p < .001$, $d = 1.24$, IC 95% [0.81, 1.66].

### Modo 2: Modelagem e Processos (Regressão, Mediação e Moderação)
- **Regressão Linear e Logística:** Relato de $R^2$, $R^2_{\text{adj}}$, $F(df_1, df_2)$, coeficientes não padronizados ($B$) com erro-padrão ($SE$) e coeficientes padronizados ($\beta$); efeito principal $R^2_{\text{adj}}$ e $f^2$, declarados sem IC quando não calculado. Em regressão logística: razão de chances (*Odds Ratio* - $OR$) com IC 95%.
- **Mediação (Processo de Hayes):** Decomposição formal em efeito direto ($c'$), efeito indireto ($ab$) e efeito total ($c$). Relato do IC de bootstrap:
  > O efeito indireto do estresse percebido sobre a exaustão emocional via ruminação foi estatisticamente significativo, $ab = 0.28$, $SE = 0.06$, IC 95% bootstrap [0.17, 0.41], confirmando a mediação indireta.

### Modo 3: Psicometria Avançada
- **Análise Fatorial Exploratória (AFE):** Método de extração (PAF/ML/ULS), critério de retenção (Análise Paralela / Hull), tipo de rotação (oblimin / geomin), variância explicada total e tabela de cargas fatoriais com comunalidades ($h^2$).
- **Análise Fatorial Confirmatória (AFC):** Especificação do estimador (WLSMV para itens ordinais, MLR para contínuos) com justificativa; versão dos índices declarada (robusta ou escalonada). Índices globais: $\chi^2(df)$, $p$, CFI, TLI, RMSEA (com IC 90%), SRMR, mais resíduos locais; cargas padronizadas com IC 95%; AVE e HTMT quando fatores correlacionam alto.
- **Modelos Bifator:** Coeficientes de confiabilidade hierárquica ($\omega_h$ para o fator geral, $\omega_{hs}$ para fatores específicos), *Explained Common Variance* (ECV) e *Percentage of Uncontaminated Correlations* (PUC).
- **Invariância de Medida:** Sequência de restrições (configural, métrica/fraca, escalar/forte, estrita) e critérios de invariância comparativa: $\Delta\text{CFI} \le .010$ e $\Delta\text{RMSEA} \le .015$ (Chen, 2007; Cheung & Rensvold, 2002).
- **Redes Psicométricas:** Estimador (`EBICglasso` com $\gamma = 0.25$ ou $0.50$), centralidade de força/influência esperada com estabilidade do índice $CS$ ($CS > 0.50$ via *bootnet*), e comparação de redes via `NetworkComparisonTest` (NCT).

### Modo 4: Auditoria de Consistência Statcheck
- Extração de todas as tuplas estatísticas do texto no formato: `estatistica(df) = valor, p = valor` (t, F, χ², r e z), em inglês ou português: vírgula decimal, ponto e vírgula entre itens, χ² em Unicode, sinal de menos tipográfico e itálico de markdown são normalizados.
- Regra de decisão portada do pacote statcheck (Nuijten et al., 2016): o $p$ é recalculado nos dois extremos do intervalo de arredondamento da estatística relatada e comparado nas casas decimais do $p$ relatado; erro grosseiro é a inconsistência que muda a decisão a α = .05 (com p = .05 contado como significativo, o padrão do pacote). `evals_relato.R` confere os veredictos contra o código do statcheck do repositório.
- Classificação:
  - **Consistente:** $p$ calculado coincide com o valor reportado dentro do arredondamento.
  - **Inconsistência Leve:** discrepância decimal sem mudança na decisão de significância estatística ($\alpha = .05$).
  - **Erro Grosseiro (Gross Inconsistency):** o texto afirma significância ($p < .05$), mas o recálculo estatístico aponta não significância ($p \ge .05$), ou vice-versa. Erros grosseiros exigem correção imediata antes da entrega ao usuário.

### Modo 5: Formatação por periódico, frases do R, efeitos e .docx (1.3.0 a 1.5.0)
- **Perfil por manuscrito.** O zero à esquerda e o separador decimal saem de um perfil escolhido para o periódico-alvo: `fmt_apa(x, tipo, perfil)`, `fmt_p_apa()`, `fmt_ic_apa()`, `listar_perfis_apa()` (em `scripts/formatador_apa.R`). Perfis: `apa7_relato` (padrão, regra 2), `apa7_estrito` (APA 6.36 literal, como o `gt1` do papaja: zero em β, TLI, RMSEA e SRMR), `psicometria_r` (sem zero em RMSEA e SRMR, zero em β), `apa7_ptbr` (0,31 com ponto e vírgula entre itens), `apa7_ptbr_sem_zero` (,31) e perfis de periódico (`rbtc`, `avaliacao_psicologica`, `prc`) que herdam de um perfil-base. Perfil de periódico só traz regra conferida na norma publicada; o que a norma não diz fica herdado e o perfil sai com `verificado = FALSE` e a URL da norma. `aviso_perfil_apa()` imprime a ressalva, que vai junto da entrega. `perfil_apa_personalizado()` monta o perfil do autor.
- **Frase sai do objeto.** `frases_apa(obj, perfil)` (em `scripts/frases_apa.R`) escreve a frase de `t.test`, `cor.test`, `chisq.test`, `lm` (modelo e cada coeficiente, com β para preditor contínuo), `lavaan` (índices robustos quando existem, senão escalonados, declarados no atributo `versao_indices`; modelo saturado é recusado) e, desde a 1.4.0, `glm` (logística e Poisson: B, EP, z, p e OR ou RP com IC 95%), `lmerMod`/`glmerMod` do lme4 (efeitos fixos com EP, t ou z, p e IC, variâncias τ₀₀ e σ² e ICC), `coxph` do survival (HR com IC 95%), `rma` do metafor (efeito combinado com IC, z, p e heterogeneidade Q, I² e τ²) e `roc` do pROC (AUC com IC e corte de Youden com sensibilidade e especificidade). A saída passa no `auditar_texto_statcheck()` em todos os perfis.
- **Efeito e IC.** `auditar_efeitos_texto()` marca cada frase com teste que não traz tamanho de efeito ou IC; `conferir_d_de_t()` (com `hedges = TRUE` para o g_s, eq. 4), `conferir_eta2p_de_F()` e, desde a 1.5.0, `conferir_omega2p_de_F()` (eq. 15, com o N total) reconstituem d, g, η²p e ω²p a partir da estatística relatada, no intervalo de arredondamento (Lakens, 2013), e acusam valor incompatível. O eval reproduz o exemplo dos filmes do artigo (t(18) = 2,52, d_s = 1,13, g_s = 1,08, ω²p = ,21).
- **Redes (1.5.0).** `frases_rede(rede, cs)` ou `frases_apa(rede)` com o objeto de `bootnet::estimateNetwork()`: nós, estimador e γ, N, arestas diferentes de zero sobre as possíveis, densidade, sinais e peso absoluto médio; com a saída de `bootnet::corStability()`, cada CS(cor = 0,7) classificado pelos cortes de Epskamp, Borsboom e Fried (2018): abaixo de 0,25 não interpretar a ordem de centralidade, entre 0,25 e 0,50 com cautela, acima de 0,50 estável. Os autores chamam esses cortes de algo arbitrários; a frase cita a fonte em vez de vender o corte como regra.
- **Word sem pacote extra.** `escrever_tabelas_docx(arquivo, tabelas)` (em `scripts/docx_apa.R`) grava uma ou mais tabelas em .docx só com R base: três linhas horizontais, nenhuma vertical, Times New Roman 12, rótulo em negrito, título em itálico, nota abaixo, `*texto*` em itálico e alinhamento por tab decimal. `tabela_para_docx()` converte a saída do gerador.

---

## Estrutura de Arquivos da Skill

Os scripts ficam no Mac em `~/.gemini/config/skills/relato-estatistico/` (cópias em `~/.claude/skills/` e `~/skills/`); a versão completa também está na página da skill no Notion.

- `SKILL.md`: Catálogo de regras operacionais e instruções de estilo.
- `references/`:
  - `01_manual_relato_apa7_jars.md`: Manual aprofundado de relato APA 7 e JARS-Quant.
  - `02_taxonomia_literatura_e_pacotes.md`: Mapeamento de literatura metodológica e pacotes R de relato.
  - `03_blueprints_fraseado_apa7.md`: Biblioteca de frases e esqueletos canônicos para cada teste.
  - `04_diretrizes_tabelas_e_figuras.md`: Regras de diagramação e tipografia de tabelas.
  - `05_auditoria_consistencia_statcheck.md`: Protocolo detalhado do auditor statcheck.
  - `06_acervo_local_e_uso.md`: Índice dos repositórios (statcheck, papaja, apaTables, rempsyc, modelsummary, gtsummary, report, equatiomatic) e papers, e as convenções que divergem entre as skills.
- `scripts/`:
  - `gerador_tabelas_apa.R`: Script R para compilar e formatar tabelas no padrão APA 7 (`estilo = "linhas"` para as três linhas horizontais; a tabela de correlações relata o n por par quando há ausentes).
  - `auditor_statcheck.R`: Script R para ler strings ou arquivos de texto e auditar NHST pela regra do statcheck.
  - `formatador_apa.R`: formatador único com perfis por periódico.
  - `frases_apa.R`: frases APA geradas de htest, lm, glm, lavaan, lme4, coxph, metafor, pROC e bootnet (`frases_rede()`).
  - `checagem_efeitos.R`: presença de efeito e IC; d, g, η²p e ω²p reconstituídos.
  - `docx_apa.R`: tabela APA em .docx com R base.
  - `evals_relato.R`: 30 testes (na 1.5, o exemplo dos filmes de Lakens, 2013, para g e ω²p e o exemplo de CS de Epskamp et al., 2018, para a frase de rede; antes, auditor com 24 veredictos conferidos contra o statcheck, perfis, frases de htest, lm, lavaan, glm, lme4, coxph, metafor e ROC, efeitos, .docx).
- `templates/`:
  - `template_secao_metodo_jars.md`: Template estruturado para a seção de Método.
  - `template_secao_resultados_apa7.md`: Template estruturado para a seção de Resultados.
