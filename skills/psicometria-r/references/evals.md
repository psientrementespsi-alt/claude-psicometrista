# Evals da psicometria-r


`evals_v22.R` (2.2, 11 testes) reproduz a Tabela 1 do relatório de Rasch do curso (BLOT, N = 150: dificuldade, EP, infit e outfit MNSQ e t do eRm, todos a .0015) e confere CSEM e equacionamento pelas propriedades das fórmulas. `evals_v21.R` (2.1, 12 testes, mais 2 pulados sem BayesFactor e difNLR) confere o fator de Bayes do t contra Monte Carlo independente e o da correlação contra a aproximação de Jeffreys, o `dif_nlr()` 2PL contra o LRT do glm (igual a 1e-6) e a recuperação de DIF plantado, e a aritmética do Angoff e do Bookmark; com BayesFactor ou difNLR instalados, compara com `ttestBF()`, `correlationBF()` e `difNLR()`. `evals_curso_didatico.R` (2.0, 15 testes) confere `recomendar_teste()` contra a árvore da planilha do curso, os critérios A, B e C de `mudanca_confiavel()` e os 7 bancos `.sav` clássicos do módulo 08 (ANOVA, ANOVA de medidas repetidas, MANOVA, regressão logística, qui-quadrado, entre outros). `evals_curso.R` (1.9, 8 testes) confere a skill contra os bancos `.sav` e os scripts do curso Psicometria-R (pasta em `Sys.getenv("PR_CURSO")`): teste t independente igual ao `t.test` do script e ao t = 3,230 e d = 0,28 do texto; `comparar_correlacoes()` coerente com o `psych::r.test` do script; `ajustar_afc()` igual ao `lavaan::cfa()` do script no LOT-R e no MRP (a skill usa a parametrização theta, o script a delta; índices iguais até a sexta casa); bootstrap em blocos igual ao direto.

Um arquivo por módulo em `scripts/` (`evals_simulados.R` para a psicometria,
`evals_tri.R`, `evals_vieses.R`, `evals_latentes_epidemio.R`, `evals_naoparametricos.R`,
`evals_mee.R`, `evals_complementos.R` para as demais funções da v1.5, `evals_v16.R`
e `evals_v17.R` para as da v1.6 e da v1.7, `evals_v18.R` para a confiabilidade da v1.8, e
`evals_<módulo>.R` para os outros), com três tipos de teste: números do livro do
usuário reproduzidos nos bancos da biblioteca do JASP (lidos de
`Sys.getenv("PR_DADOS_JASP")`, pulados sem essa variável); verdade conhecida em
dados simulados; e equivalência com outra implementação (car, afex, effectsize,
rstatix, bootnet, geepack, metafor, psych, G*Power, EFAtools, irr,
BifactorIndicesCalculator, semTools, MBESS pelos números do Wiley, difR, eRm,
poLCA, tidyLPA, meta, mada, survival e, na v1.5, rankFD, nparLD, MANOVA.RM, lme4,
netmeta, naniar, epitools, bnlearn e a saída do Joinpoint 4.9.1 da aula; na v1.6,
VGAM, ordinal, irrCAC, gtheory, dmacs, pinsearch, difR, rms, DescTools, ShinyItemAnalysis
e EFA.MRFA, mais os exemplos de Martinková e Hladká e de Petersen; na v1.7, mice,
graphicalVAR, mlVAR, mokken, mirtCAT e strucchange). Conferência
com pacote ausente é pulada com aviso, e a conta à mão ou o valor publicado cobre o
mesmo ponto.

Última execução completa no container (29/10/2026, R 4.6.1, sem os bancos do
JASP): 867 testes, todos passando (diferenças 114, associação 87, processos 44,
modelos 136, poder 108, construção 80, redes 39, psicometria 76, TRI 21, vieses 8,
latentes e epidemiologia 13, não paramétricos 16, MEE 23, complementos 30, v1.6 46,
v1.7 26); no Mac, em 29/10/2026, `evals_v18.R` (6) e `evals_v21.R` (12) passando. Última execução completa no Mac (04/11/2026, R 4.6.1,
21 suítes, sem `PR_DADOS_JASP`): 1.167 testes passando, nenhum falhando, 2 pulados (gtheory e dmacs, arquivados no CRAN)
por pacote ausente (BayesFactor, difNLR, eRm, naniar e outros) ou por banco do JASP
ausente. Rodar de
novo no computador do usuário depois de instalar ou atualizar pacotes, por exemplo
`Rscript ~/.local/share/psicometria-r/evals_mee.R`.

Achados dos evals que viraram regra de relato: deslocamento uniforme de limiares
aparece no passo dos interceptos (`invariancia.md`); o ômega por subescala num
modelo de segunda ordem é relativo ao fator geral (`afc-ordinal.md`); com itens
de 3 categorias o passo dos limiares fica sem teste; o R 4.6 mudou o
`wilcox.test` com empates e zeros, e a skill fixa a regra do JASP 0.95.1
(`diferencas.md`); `mirt::marginal_rxx()` ignora a variância latente estimada na
família Rasch (`tri.md`); o MBESS aproxima o N mínimo no poder do RMSEA
(`poder-amostra.md`); o `mctp` com log-odds das aulas devolve o estimador em escala
transformada, e o nparLD 2.3 removeu os argumentos das aulas
(`naoparametricos-complexos.md`); os IC do Joinpoint do relato da aula não se
reproduzem, só as estimativas (`epidemiologia.md`); TRI thurstoniana com itens de mesmo
sinal recupera mal os traços (`vieses-resposta.md`); o kappa despenca com prevalência
desigual, como na aula de limitações do kappa, e o AC1 não (`construcao-instrumentos.md`);
no próprio banco a inclinação de calibração da logística vale 1 por construção, e só a
corrigida pelo bootstrap informa (`associacao.md`); o pacote dmacs usado por Petersen
combina DPs em vez de variâncias e escala mal o integral quando a variância latente
focal não é 1 (`invariancia.md`); no Mokken, recodificar item a item para começar em 0
muda H e a ordenação dos itens, e a skill desloca todos juntos (`tri.md`); no configural
com std.lv o ΔMean reproduz a diferença observada por construção (`invariancia.md`); o
texto do curso imprime DP = 3,71 onde o banco dá 3,72 e mistura p uni e bicaudal na
comparação de correlações dependentes; no DIF 2PL a regressão
logística generalizada coincide com o glm, e com acerto ao acaso a diferença de c
pode aparecer como diferença de b (identificação fraca), então o LRT conjunto é o que se lê (`tri.md`).
