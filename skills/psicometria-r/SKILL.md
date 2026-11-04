---
name: "psicometria-r"
description: "Análise de dados no R do usuário, com a estatística da psicologia (descritivas e ausentes, testes t, ANOVAs, ANCOVA, MANOVA, não paramétricos fatoriais, correlação, regressões (inclusive ordinal e multinomial), ROC, qui-quadrado, moderação, mediação, GLM, GEE, multinível, meta-análise em rede, de prevalência e de acurácia, sobrevida, joinpoint, poder no G*Power e para AFC) e a psicometria (AFE, análise de itens, AFC ordinal, MEE e curva de crescimento, ômega e bifator, invariância e dMACS, TRI e MIRT, DIF, validade de conteúdo, reteste, concordância (kappa e AC1), teoria G, EMA, aquiescência, desejabilidade, escolha forçada, classes e perfis latentes, redes, EGA e redes bayesianas). Reproduz o JASP do livro e entrega tabelas e frases APA. Use para rodar, conferir, escolher ou planejar análise estatística ou psicométrica, calcular N, ler output do R, JASP ou Mplus, ou desenhar a validação do PBAT. Redigir o manuscrito é com fabrica-de-artigo; revisar, com verificacao-artigo; literatura, com pesquisa-academica."
metadata:
  version: '2.2.0'
---

# SKILL: Análise de dados e psicometria no R (v2.2)

Roda e interpreta, no R do computador do usuário, as análises de um estudo em
psicologia e entrega um pacote de resultados pronto para a fabrica-de-artigo:
estatística descritiva e inferencial (diferenças, inclusive não paramétricas em
delineamentos complexos, associação, processos, modelos para medidas repetidas e
dados aninhados, meta-análise, epidemiologia, poder, fator de Bayes) e psicometria (construção de
itens, AFE, AFC e MEE, bifator, confiabilidade, invariância, TRI e DIF, vieses de
resposta, classes latentes, redes, EGA e redes bayesianas, pontos de corte).

<required_reading>
0. Se uma skill `regras-comuns` estiver instalada, leia-a antes de começar; se não estiver, siga as regras desta skill.
1. Antes de rodar uma análise, leia o arquivo de `references/` do módulo dela
   (tabela "Módulos" abaixo). As decisões das funções estão explicadas lá.
2. Na primeira análise da sessão, leia `references/script-reproduzivel.md`.
</required_reading>

## IRON RULES (desta skill, somadas às do hub)

1. **Dado de participante não sai do computador do usuário.** A análise roda no
   R local. Nunca copiar planilha de participantes para o container da nuvem,
   nunca colar linhas na conversa, nunca rodar `head()`, `print(dados)`,
   `View()` ou qualquer coisa que devolva linha individual ou ID. O que volta
   para a conversa é agregado: n, frequência por categoria, médias, cargas,
   índices, correlações. Ao levar para a conversa ou para o manuscrito uma
   tabela cruzada de variáveis demográficas, célula menor que 5 é suprimida (as
   contagens do `qui_quadrado()` não são mascaradas, então a supressão é tarefa
   de quem relata); o print de `descritivas()` por grupo já apaga as estatísticas
   de grupo com n < 5. As funções da skill só imprimem agregados; as que criam
   variável por pessoa (índices de atípicos, `$indice` e `$dados_recentralizados`
   de `controle_aquiescencia()`, `$classe`, `$posterior` e `$linhas` das classes e perfis latentes,
   `$rci` de `mudanca_confiavel()`) guardam no objeto, e esses elementos nunca são
   impressos nem levados à conversa. Gráfico com um ponto por pessoa fica no
   computador. As ferramentas do btw que descrevem objetos da sessão
   (`btw_tool_env_describe_data_frame`, `btw_tool_env_describe_environment`)
   devolvem linhas e ficaram fora do servidor r-btw do Mac; se alguém religar o
   grupo "env", elas continuam proibidas em dado de participante.
2. **Todo número do relato sai do output do R desta sessão.** Nada de estimar,
   arredondar de cabeça ou completar índice ausente. O que o R não calculou é
   "não calculado", com o motivo.
3. **Nunca inventar argumento de função.** Usar as funções da skill; fora delas,
   conferir a ajuda (`btw_tool_docs_help_page` no servidor r-btw, ou `args()`).
   O lavaan 0.7.3 documenta argumentos em snake_case e aceita a grafia antiga
   com ponto; conferir nomes e padrões novos antes de reaproveitar código antigo.
4. **A fonte do usuário decide a rota.** Ordem (`references/fontes-e-decisoes.md`,
   seção 1): o livro do usuário (Biblioteca Análise de Dados, 7 volumes, JASP
   0.95.1); as notas e o curso Psicometria Online Academy; a pasta Psicometria-R
   (seis livros-base de psicometria, recortes por módulo, auditoria CoVe, curso
   completo com código); a literatura com DOI conferido, que entra como argumento
   ou ao lado. Toda decisão sai escrita em `$decisoes`. Erro de conta, de
   transcrição ou de método das fontes não é reproduzido (lista no mesmo
   arquivo, seção 5). Se o usuário mandar o artigo original do instrumento, o
   manual ou o protocolo, a estrutura, o formato de resposta e a pontuação vêm
   dali.
5. **Pressuposto não é liga-desliga.** Normalidade e homogeneidade são lidas por
   magnitude, gráfico e n, como no livro; p do Shapiro-Wilk ou do Levene é um
   indício entre vários. A alternativa (Welch, não paramétrico, bootstrap,
   correção de esfericidade) sai ao lado e a escolha fica justificada.
6. **Tamanho de efeito com IC; não significativo não é igual.** O índice
   principal é o que o livro manda relatar (coluna `principal` em `$efeito`) e
   sai com IC sempre que há método estabelecido; onde a função não dá IC, o relato
   diz que não foi calculado. Rótulos de magnitude são régua, não lei, e a régua
   usada é declarada com a fonte. p ≥ α não mostra ausência de efeito: para
   afirmar equivalência, `teste_equivalencia()` com limites justificados, ou o
   fator de Bayes (`bf_teste_t()`, `bf_correlacao()`) com a prior declarada.
7. **Estimador declarado e justificado; versão dos índices declarada.** Na AFC
   com WLSMV ou MLR, a tabela nunca mistura índices escalonados e robustos, e o
   texto diz qual foi usada. Corte de índice de ajuste é referência, não teste
   (Hu & Bentler, 1999, derivados para ML; Xia & Yang, 2019).
8. **Modificação de modelo só com justificativa teórica**, sempre relatada.
9. **Semente fixada e versões registradas** em toda análise com reamostragem ou
   partida aleatória (bootstrap, Monte Carlo, bootEGA, NCT, redes, análise
   paralela, LCA, TRI): argumento `semente` (padrão 2026, o do livro) e
   `versoes_pacotes()` no LOG.
10. **Referência só da lista verificada** (`references/referencias-verificadas.md`:
    entradas com DOI conferido no Crossref, inclusive as da v2.1, mais as citadas
    sem DOI, com a âncora de cada afirmação). Fora dela: pesquisa-academica e
    `verify_citations.py`. Referência de software sai de `citation("pacote")`.
11. **Todo número do relato tem código salvo que o refaz.** O script `.R` fica ao
    lado dos dados e passa por uma reexecução em R limpo antes da entrega
    (`references/script-reproduzivel.md`).
12. **Rota fora do escopo testado é dita.** Pacote ou modelo sem função e sem
    eval na skill (equacionamento com âncora (NEAT) ou por TRI, Rasch politômico por CML, redes com variação no tempo e GVAR de painel, MTMM,
    curva de crescimento de segunda ordem, interação latente com item ordinal, três
    passos com covariáveis, poder para GLM e GEE, ANOVA e regressão bayesianas) só roda com a documentação do pacote aberta, e o relato e o
    LOG dizem que a rota saiu do escopo testado.

## Economia de tokens

1. Scripts de `scripts/` são executados com `source()`, nunca lidos como contexto. Para a assinatura de uma função, use `args(nome)` no R ou `grep -n "nome <- function" scripts/*.R`.
2. Leia só a referência do módulo do pedido. Referência longa abre com sumário e números de linha: leia o sumário e depois só a seção necessária.
3. Não releia arquivo já lido nesta conversa; reaproveite o que está no contexto.
4. Saída do R enxuta: resumos e agregados, nunca objetos inteiros (o que também protege dado de participante).
5. Todo cálculo, conferência ou tabela determinística sai de função da skill; não refazer conta à mão em tokens.
6. Modo rápido, para exploração e iteração: só números e decisões, sem a verificação narrativa completa. Na entrega final, a verificação completa e o checklist da skill são obrigatórios.
7. Modelo: execução mecânica (rodar script, formatar tabela) cabe num modelo leve; interpretação e decisão metodológica pedem o modelo mais forte. Em lote via API, Batch API (50% de desconto) e cache dos materiais relidos.

## Roteamento

| Pedido | Skill |
|---|---|
| Rodar, escolher, conferir ou planejar análise estatística ou psicométrica; calcular N | esta |
| Higiene do banco bruto antes da análise (C/IER, atípicos, ausentes, fluxograma) | higiene-dados |
| Formatar o relato por periódico, tabela .docx, auditoria statcheck | relato-estatistico |
| Redigir Método e Resultados | fabrica-de-artigo (modo EMPIRICO ou PSICOMETRICO), recebendo o pacote de resultados daqui |
| Verificar manuscrito antes da submissão | verificacao-artigo |
| Buscar referências, conferir DOI | pesquisa-academica (`verify_citations.py`) |
| Diagrama de mediação | diagram-auditor-apa7 (auditar) |
| Pôster com os resultados | poster-estrategista |

## Onde a análise roda

1. **Servidor MCP r-btw no Claude Desktop do Mac** (preferido). Executar com
   `btw_tool_run_r`; conferir argumentos com `btw_tool_docs_help_page`. Cada
   chamada faz uma coisa só e devolve só agregados. Toda chamada que mexe nos
   dados passa por `registrar()` (`references/script-reproduzivel.md`). Cada chamada
   do r-btw tem limite de 60 s: trabalho longo vai em blocos que gravam `.rds`, ou
   num `Rscript` disparado com `system2(..., timeout = 50)` para não travar a sessão.
2. **Shell do computador** (sem r-btw): script `.R` ao lado dos dados, `Rscript`.
3. **Container da nuvem**: só com dados simulados ou públicos (evals, bancos da
   biblioteca do JASP, exemplos de pacote). Nunca com dado de participante.

Carregar tudo com `source("~/.local/share/psicometria-r/carregar.R")`. Conferir
que `.ee_versao` é "2.2.0"; se não, copiar `scripts/` de novo. `checar_pacotes()`
lista o que falta; pacote ausente só é instalado com pedido explícito do usuário.

## Script reproduzível (resumo)

A abertura cria `analise_<data>.R` ao lado dos dados, copia as funções da skill
para a pasta e define `registrar()`; a leitura do banco entra com o md5 do
arquivo; cada análise entra como bloco de texto cru rotulado por tipo
(descritivo, confirmatório, exploratório) e destino (relatado ou não); o
fechamento registra o bloco `numeros`, anexa a sessão e reexecuta o script em R
limpo (`status` 0 e `all.equal` TRUE liberam o Passo 3). Código completo, regras
e testes em `references/script-reproduzivel.md`.

## Módulos

| Pergunta | Funções principais | Ler antes |
|---|---|---|
| Qual teste usar para o meu delineamento? | `recomendar_teste()` | `referencias-verificadas.md` |
| Como são os dados? Ausentes, normalidade, atípicos, mudança clínica | `descritivas()`, `ausentes()`, `imputacao_multipla()`, `ler_normalidade()`, `atipicos()`, `mudanca_confiavel()` (RCI e critérios de recuperação A, B e C de Jacobson e Truax) | `fontes-e-decisoes.md`, `figuras-dados.md` |
| Grupos ou condições diferem? São equivalentes? | `teste_t_independente()`, `teste_t_pareado()`, `teste_t_uma_amostra()`, `teste_equivalencia()`, `anova_um_fator()`, `anova_medidas_repetidas()`, `anova_fatorial()`, `ancova()`, `manova_grupos()`, `teste_mann_whitney()`, `teste_wilcoxon_pareado()`, `teste_kruskal_wallis()`, `teste_friedman()` | `diferencas.md` |
| Quanto a evidência favorece H1 ou H0? (fator de Bayes) | `bf_teste_t()` (uma amostra, pareado, independente), `bf_correlacao()` | `diferencas.md`, `associacao.md`; priors padrão do JASP |
| Fatorial, longitudinal ou multivariado sem normalidade | `np_fatorial()`, `np_longitudinal()`, `manova_np()` | `naoparametricos-complexos.md` |
| Variáveis se associam? Quem prediz o desfecho? O escore prevê igual nos grupos? | `correlacoes()`, `correlacao_parcial()`, `correlacao_ponto_bisserial()`, `comparar_correlacoes()`, `regressao_linear()`, `regressao_hierarquica()`, `regressao_logistica()` (com `validar = TRUE`), `curva_roc()`, `vies_preditivo()`, `qui_quadrado()` (RR, DR e NNT no 2 × 2), `qui_quadrado_aderencia()`, `mcnemar()`, `q_cochran()` | `associacao.md` |
| Como (mediação), quando (moderação), como e quando | `moderacao()`, `mediacao()`, `mediacao_moderada()` | `processos.md` |
| Desfecho não normal, ordinal ou nominal, medidas repetidas, dados aninhados, síntese de estudos | `escolher_familia()`, `modelo_glm()`, `modelo_ordinal()`, `modelo_multinomial()`, `modelo_gee()`, `multinivel_linear()`, `multinivel_logistico()`, `efeito_contextual()`, `meta_analise()`, `meta_subgrupos()`, `meta_regressao()`, `meta_proporcao()`, `meta_acuracia()`, `meta_rede()`, `extrair_media_dp()` | `modelos.md` |
| Estudo epidemiológico: razão de prevalência, sobrevida, tendência, taxa padronizada | `razao_prevalencia()`, `sobrevida()`, `tendencia_joinpoint()`, `taxa_padronizada()` | `epidemiologia.md` |
| Subgrupos latentes de pessoas e variável externa | `classes_latentes()`, `perfis_latentes()`, `classes_distal()` | `classes-perfis.md` |
| Quantos participantes? | `poder_correlacao()`, `poder_teste_t()`, `poder_anova()`, `poder_anova_mr()`, `poder_regressao()`, `poder_logistica()`, `poder_quiquadrado()`, `poder_afc_rmsea()`, `n_aipe_rmsea()`, `poder_afc_aninhados()`, `n_aipe_parametro()`, conversões `efeito_*()` | `poder-amostra.md` |
| Itens e escala: conteúdo, itens, estrutura exploratória, fidedignidade, tamanho do teste, correção de r | `validade_conteudo_cvc()`, `validade_conteudo_ivc()`, `validade_conteudo_rvc()`, `validade_conteudo_aiken()`, `analise_classica_itens()`, `analise_distratores()`, `afe_adequacao()`, `afe_numero_fatores()`, `afe()`, `confiabilidade_metades()`, `confiabilidade_teste_reteste()`, `concordancia_juizes()`, `confiabilidade_ema()`, `generalizabilidade()`, `spearman_brown()`, `corrigir_correlacao()`, `normas()` | `construcao-instrumentos.md` |
| Ponto de corte por julgamento de especialistas | `corte_angoff()` (Angoff modificado), `corte_bookmark()` (caderno ordenado pela TRI com critério RP) | `construcao-instrumentos.md`; Cook e Pitoniak (2026, cap. 12) |
| Estrutura confirmatória, confiabilidade, bifator, escore somado, validade | `auditar_itens()`, `sugerir_estimador()`, `ajustar_afc()`, `tabela_ajuste()`, `confiabilidade()`, `confiabilidade_bootstrap_blocos()`, `confiabilidade_compostos()`, `indices_bifator()`, `checar_escore_somado()`, `validade_discriminante()` | `afc-ordinal.md` |
| Modelo estrutural: efeito indireto, moderação por grupo, interação latente, curva de crescimento | `mee()`, `mee_multigrupo()`, `interacao_latente()`, `curva_crescimento()` | `mee.md` |
| Invariância entre grupos, ondas ou ao longo de idade; equivalência do ajuste; tamanho da não invariância; viés de medida e preditivo | `invariancia()`, `invariancia_longitudinal()`, `invariancia_continua()`, `equivalencia_ajuste()`, `dmacs()`, `vies_preditivo()`, `tri_dif()` | `invariancia.md`, `mee.md` |
| Rasch dicotômico pela rota do eRm do curso: CML, infit e outfit, Andersen, Wald, PCA dos resíduos | `rasch_cml()` | `tri.md` |
| Erro-padrão condicional de medida; equacionamento de formas | `csem()` (binomial de Lord, Thorndike, TRI), `equacionar()` (média, linear, equipercentil) | `construcao-instrumentos.md`; Cook e Pitoniak (2026, caps. 5 e 11) |
| Itens pela TRI, DIF, TRI multidimensional, Mokken, CAT pós-hoc | `tri_ajustar()`, `tri_comparar()`, `tri_dif()`, `dif_nlr()` (regressão logística generalizada com acerto ao acaso e desatenção, rota do difNLR), `tri_multidimensional()`, `escala_mokken()`, `simular_cat()` | `tri.md` |
| Aquiescência, desejabilidade social, escolha forçada | `controle_aquiescencia()`, `controle_desejabilidade()`, `escolha_forcada()` | `vieses-resposta.md` |
| Redes, centralidade, EGA, NCT, rede bayesiana, rede temporal de EMA (N = 1 e multinível) | `rede_psicometrica()`, `estabilidade_rede()`, `ponte()`, `previsibilidade()`, `rede_bayesiana()`, `rede_temporal()` | `redes-ega.md` |
| PBAT ou MPFI-24 | conforme o desenho | `pbat-mpfi.md` |

As funções de análise devolvem um objeto `pr_resultado`: `$pressupostos`,
`$tabela`, `$efeito` (com a coluna `principal`), `$posthoc`, `$decisoes`,
`$avisos`, `$apa` (frase pronta; `relatar(res)`) e, quando há, `$tabelas_extra`,
com os argumentos comuns `decimal = ","` para manuscrito em português (itens da
frase separados por ponto e vírgula), `idioma = "en"` e `nivel` do IC. As funções
de poder do G*Power devolvem `pr_poder`, com a tela do G*Power para replicar; as
de AFC devolvem `pr_resultado` com uma tabela de campos. Exceções: `descritivas()`
devolve tabela; `ler_normalidade()`, texto; `efeito_*()`, número;
`centralizar_multinivel()`, o banco (invisível); `rvc_critico()`,
`calcular_efeitos_meta()`, `extrair_media_dp()`, `ajustar_perda()` e `checar_pacotes()`, tabelas;
`versoes_pacotes()`, texto. As de `psicometria.R` (AFC, confiabilidade,
invariância) devolvem tabelas, listas e o objeto lavaan, e o relato do ajuste sai
de `formatar_ajuste_apa()`; `confiabilidade(ic = TRUE, metodo = "bootstrap")` dá o IC percentílico do ômega do escore observado e do alfa, a rota do relato com itens ordinais (`afc-ordinal.md`). Com 1.000 reamostras o bootstrap passa do limite de 60 s do r-btw: `confiabilidade_bootstrap_blocos(fit, dados, pasta)` roda em blocos gravados na pasta (pelo menos um bloco por chamada, até `max_segundos`) e continua na chamada seguinte, com o mesmo plano de reamostras e o mesmo IC do bootstrap direto; pasta com plano de outro modelo, banco, semente ou B é recusada. O padrão `metodo = "delta"` dá `omega_delta` (variável-resposta latente, igual ao observado só com itens contínuos) com EP e IC Wald, mais o alfa com IC de Feldt;
`indices_bifator()`, `checar_escore_somado()` e `dmacs()` devolvem `pr_resultado`.
O fator de Bayes (`est_bayes.R`) sai em R base por integração numérica, com a prior padrão do JASP (Cauchy r = sqrt(2)/2 no t, beta esticada kappa = 1 na correlação), robustez em outras larguras e rótulo de força como régua; `dif_nlr()` (`est_dif_nlr.R`) estima por máxima verossimilhança e LRT, enquanto o difNLR usa por padrão mínimos quadrados não lineares; `corte_angoff()` e `corte_bookmark()` (`est_padroes.R`) devolvem o corte, o EP entre juízes e a tabela de itens ou o caderno ordenado.

## Fluxo

**Passo 0. Delineamento.** Antes de qualquer código: pergunta (se, como,
quando, como e quando; diferença, equivalência ou associação), desenho (entre ou
intra sujeitos, aninhamento, medidas repetidas, tempo até evento), tipo de cada
variável, n e, se o estudo ainda vai coletar, o cálculo amostral (`poder_*()`;
para AFC, a pergunta do capítulo 5 do Wiley que o estudo faz); `recomendar_teste()`
confere a escolha contra a árvore do curso. Em estudo
psicométrico, objeto (conteúdo, estrutura, bifator, adaptação, invariância, TRI,
rede, vieses, ponto de corte) e estrutura original da fonte do usuário; COSMIN em
`fabrica-de-artigo/references/modo-psicometrico.md`. PBAT (mestrado) ou MPFI-24:
ler `references/pbat-mpfi.md`.

**Passo 1. Dados.** `descritivas()` e `ler_normalidade()` por variável e grupo;
ausentes declarados (as funções fazem deleção listwise e informam quantos saíram;
correlações e redes usam pairwise por padrão, como o JASP; rota para ausentes na
AFC em `afc-ordinal.md`); atípicos por critério declarado antes da análise. Em
escala, `auditar_itens()` e `analise_classica_itens()`.

**Passo 2. Análise.** A função do módulo, com a rota padrão e a alternativa ao
lado quando houver. Ler `$pressupostos`, `$decisoes` e `$avisos` antes do
resultado.

**Fechamento do script.** Bloco `numeros`, sessão e reexecução em R limpo. Sem
`status` 0 e `all.equal` TRUE, não há Passo 3.

**Passo 3. Pacote de resultados.** Tabelas e frases de `$apa`, no formato de
`references/relato-apa.md`, com estimador ou rota, versão dos índices, sementes e
versões dos pacotes. Esse pacote é o insumo da fabrica-de-artigo.

**Passo 4. LOG** no template do hub, com `versoes_pacotes()`, sementes, decisões
de rota, colapso de categorias e modificações de modelo em `notas_revisao`, o
caminho do script, o `tools::md5sum()` dele e do banco e o resultado da
reexecução.

## Antes de entregar

- [ ] Nenhuma linha de participante saiu do computador do usuário; nenhum vetor por pessoa (`$indice`, `$classe`, `$posterior`, `$rci`, dados recentralizados) foi impresso
- [ ] Rota de cada análise declarada (clássica do livro ou alternativa) e alternativa calculada ao lado quando a função oferece
- [ ] Tamanho de efeito principal com IC em todo teste; régua de magnitude declarada com a fonte
- [ ] Resultado não significativo sem leitura de igualdade; equivalência só com TOST e limites justificados, ou fator de Bayes com a prior e a robustez relatadas
- [ ] Pós-hoc só com ômnibus significativo e com a correção declarada
- [ ] Na AFE: paralela e Hull relatados, discordância discutida pela interpretabilidade; MAP só como auxiliar
- [ ] Na AFC: estimador justificado, índices completos com versão declarada e ressalva de Xia e Yang (2019) com WLSMV; escore somado só com o modelo paralelo testado (`checar_escore_somado()`)
- [ ] Com fator geral plausível: bifator com ωH, ωHS, ECV e PUC, réguas declaradas como réguas
- [ ] Invariância com sequência, critério e teste de diferença; passo que falhou investigado; DIF e invariância não somados como evidências independentes
- [ ] TRI: modelo comparado, ajuste global com a estatística declarada (M2, M2* ou C2), faixa de precisão de θ relatada
- [ ] Reteste pelo ICC(A,1) com EPM; concordância com kappa ou ICC conforme a escala
- [ ] Ponto de corte: método, número de juízes, EP entre juízes e, no Bookmark, o RP declarado
- [ ] Classes ou perfis: todas as soluções na tabela, critério que decidiu e entropia da escolhida
- [ ] Rede com CS por índice e centralidade só interpretada se o CS permitir; frase de alcance causal; rede temporal com pontos por pessoa, destendência e leitura de precedência, não causa
- [ ] Ausentes: mecanismo discutido; imputação múltipla com m, método, FMI e casos completos ao lado
- [ ] Mediação sem rótulo completa/parcial e com cautela causal conforme o desenho; efeito indireto na MEE com IC bootstrap
- [ ] Não paramétrico em delineamento complexo: ATS como teste principal e efeito relativo com IC
- [ ] Meta-análise em rede: heterogeneidade e inconsistência (Q entre desenhos, netsplit) lidas antes do ranking
- [ ] Joinpoint: método do IC declarado (paramétrico ou empírico); escolha forçada com itens de sinais mistos; BCH com a entropia relatada
- [ ] Rota fora do escopo testado dita no relato e no LOG
- [ ] Semente e versões dos pacotes no LOG
- [ ] Script `analise_<data>.R` ao lado dos dados, com a pasta `psicometria-r/` copiada, o md5 do banco no bloco de leitura e todo bloco que gerou número do relato rotulado por tipo e destino
- [ ] Reexecução em R limpo com status 0 e `all.equal` TRUE, ou divergência explicada no LOG
- [ ] Nenhum ID de participante escrito no código do script
- [ ] Toda referência da lista verificada ou conferida com `verify_citations.py`

## Evals

Testes, cobertura de cada suíte e última execução: `references/evals.md`. Ler só quando o pedido for rodar, conferir ou ampliar os testes.

