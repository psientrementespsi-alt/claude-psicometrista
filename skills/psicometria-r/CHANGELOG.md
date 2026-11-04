# CHANGELOG, psicometria-r

## 2.2.0 (29/10/2026)

- Novo `est_rasch_cml.R`: `rasch_cml()`, Rasch dicotômico por máxima verossimilhança condicional em R base, a rota do eRm dos scripts do curso: dificuldades com soma zero e EP, theta por ML, infit e outfit (MNSQ e t), teste de Andersen (com a exclusão de itens extremos num subgrupo, como o LRtest), Wald por item, PCA dos resíduos e confiabilidade de separação. Reproduz a Tabela 1 do relatório de Rasch do curso (BLOT) a .0015.
- Novo `est_medida.R`: `csem()` (binomial de Lord, Thorndike e TRI; Cook & Pitoniak, 2026, cap. 5, eq. 71-78) e `equacionar()` (média, linear e equipercentil com continuização por posto percentil; Cook & Pitoniak, 2026, cap. 11, eq. 15-16; Kolen & Brennan, 2014).
- `regressao_linear()`: diagnóstico de influência ganhou DFFITS e DFBETAS (Belsley et al., 1980), que os scripts de regressão do curso usam.
- Novo `evals_v22.R` (11 testes).
- `.ee_versao` = "2.2.0".
- Conferência de 04/11/2026: as 21 suítes de evals rodaram no Mac (1.167 testes passando, nenhum falhando, 2 pulados: gtheory e dmacs, arquivados no CRAN; os demais pulados foram cobertos instalando os pacotes de comparação).
- SKILL.md unificado em 04/11/2026: o texto do Mac recebeu as linhas que só a versão da conta do Claude tinha (rotas para as outras skills, notas de uso), e a conta recebe as desta versão.

## 2.1.0 (29/10/2026)

- Novo `est_bayes.R`: `bf_teste_t()` (uma amostra, pareado e independente; JZS de Rouder et al., 2009, prior Cauchy r = sqrt(2)/2 do JASP, com robustez em r = 1 e sqrt(2)) e `bf_correlacao()` (prior beta esticada, kappa = 1 do JASP; Ly et al., 2016), em R base por integração numérica, sem o BayesFactor.
- Novo `est_dif_nlr.R`: `dif_nlr()`, DIF por regressão logística generalizada (modelos 2PL, 3PLc, 3PLd e 4PL; eq. 9.9 de Martinková & Hladká, 2023), por máxima verossimilhança e LRT, sem o difNLR. Com 2PL o LRT é o do glm.
- Novo `est_padroes.R`: `corte_angoff()` (Angoff modificado, com EP e IC entre juízes e a tabela de itens) e `corte_bookmark()` (caderno ordenado por RP, theta de corte por juiz e escore esperado pela curva característica do teste), a partir de Cook & Pitoniak (2026, cap. 12).
- Novo `evals_v21.R`: 12 testes passando; 2 comparações puladas porque BayesFactor e difNLR não estão instalados no Mac.
- `referencias-verificadas.md`: Rouder et al. (2009), Ly et al. (2016) e Drabinová & Martinková (2017) conferidas no Crossref; Cook & Pitoniak (2026) entre as citadas sem DOI.
- `.ee_versao` = "2.1.0".

## 2.0.0 (29/10/2026)

- Novo `scripts/est_roteamento.R`: função algorítmica declarativa `recomendar_teste()`, implementando 100% da árvore de decisão da planilha didática do curso (`A Escolha do Teste Estatistico Apropriado (Final).xlsx`) para guiar a seleção de testes univariados e multivariados, rotas paramétricas vs não-paramétricas, alternativas de bootstrap e funções exatas da skill.
- Expansão de `mudanca_confiavel()` em `est_comum.R`: suporte formal aos Critérios de Recuperação Clínica A, B e C de Jacobson & Truax (1991), classificando os participantes em 4 categorias (*deteriorado*, *inalterado*, *melhora confiável sem recuperação*, *clinicamente recuperado*).
- `.ee_versao` = "2.0.0".

## 1.9.0 (29/10/2026)

- Nova `confiabilidade_bootstrap_blocos()`: o bootstrap da confiabilidade em blocos gravados numa pasta, para caber no limite de 60 s do r-btw. O plano de reamostras usa a mesma sequência de `sample.int` que `confiabilidade(metodo = "bootstrap")`, então o IC final é o mesmo; pasta com plano de outro modelo, banco, semente ou B é recusada.
- `.ee_versao` = "1.9.0".

## 1.8.1 (29/10/2026)

- **`confiabilidade(ic = TRUE, metodo = "delta")` corrigida.** A 1.8.0 aplicava o gradiente de (soma das cargas)^2 / ((soma)^2 + soma dos residuos) sem a variancia do fator e sem os parametros fixos (primeira carga do item marcador e residuos dos itens ordinais na parametrizacao delta), e centrava o IC no omega de `compRelSEM()`, que e outra quantidade. Na DASS-21 (N = 116, WLSMV) o IC da depressao saiu [.881, 1.000]. Agora o omega entra como parametro definido num refit com `std.lv = TRUE` (theta = 1 - lambda^2 nos itens ordinais, theta rotulado nos continuos), e o delta method do lavaan da `omega_delta`, `omega_ep` e o IC Wald centrado nele. A coluna `omega` (compRelSEM, escore observado) continua sem mudanca. Modelos de segunda ordem, com carga cruzada ou covariancia residual ficam sem IC delta, com aviso.
- IC de Feldt do alfa passa a usar o n listwise de cada fator, e nao `nrow(dados)`.
- `evals_v18.R`: o teste do delta confere `omega_delta` contra o parametro definido escrito a mao, em itens continuos (onde `omega_delta` = `compRelSEM`) e ordinais.
- SKILL.md: a frase sobre `confiabilidade(ic = TRUE)` dizia bootstrap, mas desde a 1.8.0 o padrão é `metodo = "delta"`; agora diz o que cada método entrega e que o ômega do escore observado com itens ordinais pede `metodo = "bootstrap"`.

## 1.8.0 (29/10/2026)

Quarta expansao sistematica e saneamento de robustez de execucao:
- **Confiabilidade com Metodo Delta Instantaneo (`confiabilidade()`)**:
  - Introducao do argumento `metodo = c("delta", "bootstrap")`, sendo `"delta"` o padrao.
  - Calculo analitico em forma fechada do erro-padrao assintotico e intervalo de confianca de 95% de Wald para o ômega de McDonald (Raykov, 2001; Kelley & Cheng, 2012) a partir da matriz de covariancia assintotica `vcov(fit)` em milissegundos (< 0.1s), sem necessidade de bootstrap.
  - Calculo exato via distribuicao F para o intervalo de confianca de 95% do alfa de Cronbach (Feldt, 1965; van Zyl et al., 2000).
  - Correcao estrutural do modo `"bootstrap"`: sanitizacao automatica de `fit@call` contra variaveis de escopo fechado, eliminando falhas de convergencia por simbolos nao avaliados.
- **Saneamento de Chamada em `ajustar_afc()`**:
  - Reconstrucao da chamada `lavaan::cfa` com avaliacao explicita dos argumentos (`estimator`, `ordered`, `missing`, `parameterization`), garantindo compatibilidade total com `lavaan::update()` e `semTools::compRelSEM()`.
- Versao da skill atualizada para 1.8.0.

## 1.7.0 (29/10/2026)

Terceira conferência da pasta Psicometria-R: pacotes chamados nos códigos dos livros e das
aulas que a skill ainda não usava (mice, graphicalVAR, mlVAR, mokken, mirtCAT,
strucchange, o equiv_chi de Petersen), cruzados com os capítulos que os usam. O que tinha
código e exemplo virou função testada; o PBAT ganha a rede temporal que o
`pbat-mpfi.md` marcava como lacuna. `.ee_versao` passa a "1.7.0".

Funções novas:

- **Imputação múltipla** (est_comum.R): `imputacao_multipla()`, mice com regras de Rubin
  e gl de Barnard-Rubin, FMI, R² combinado e casos completos ao lado (igual ao
  mice::pool).
- **Rede temporal** (est_redes.R): `rede_temporal()`, GVAR de lag 1 com graphicalVAR (uma
  pessoa) ou mlVAR (multinível), destendência e avisos de pontos no tempo (Isvoranu et
  al., 2022, cap. 9 e 10).
- **Mokken** (est_tri.R): `escala_mokken()`, H, Hi e Hij com EP, AISP, monotonicidade e
  IIO (Wiley, cap. 15).
- **CAT pós-hoc** (est_tri.R): `simular_cat()`, mirtCAT com parada por EP, só agregados
  (Martinková e Hladká, 10.2.6; Petersen, cap. 20).
- **Normas** (est_construcao.R): `normas()`, percentil, z, T, escalar, padrão e estanino
  por estrato (Petersen, cap. 1).
- **Invariância** (psicometria.R): `equivalencia_ajuste()`, teste de equivalência do χ²
  de Counsell et al. (2020), que reproduz os três exemplos de Petersen (15.10.3.5), e
  `invariancia_continua()`, testes baseados no escore ao longo de idade ou outra
  variável contínua ou ordinal (strucchange).

Referências: 14 novas com DOI conferido (272 no total). Evals: `evals_v17.R` (26 testes).
CAT e redes temporais saem da lista de rotas fora do escopo testado.

## 1.6.0 (29/10/2026)

Segunda revisão linha a linha da pasta Psicometria-R: varredura dos 174 pacotes e das
funções chamadas nos scripts das aulas, dos termos de método nas 655 transcrições, dos
sumários dos seis livros (Martinková e Hladká, Wiley Handbook, Cook e Pitoniak, Petersen,
Isvoranu et al., Rust et al.), do mapa modular e da auditoria CoVe. O que o curso ou os
livros fazem com código e a skill não fazia virou função com eval contra o exemplo do
livro ou um pacote de referência. `.ee_versao` passa a "1.6.0".

Funções novas:

- **Desfecho ordinal e nominal** (est_modelos.R): `modelo_ordinal()` (MASS::polr,
  parametrização do PLUM, Wald tipo III, teste das linhas paralelas global e por variável
  pelo ordinal::clm) e `modelo_multinomial()` (nnet::multinom, RRR, razão de
  verossimilhança por termo como o NOMREG). `escolher_familia()` aponta para elas.
- **Efeito contextual** (est_modelos.R): `efeito_contextual()`, b entre grupos menos b
  intragrupo com EP, a pergunta da aula de regressão multinível.
- **Teoria da generalizabilidade** (est_construcao.R): `generalizabilidade()`, estudos G
  (REML) e D (Eρ² e Φ, com projeção de tamanhos), uma ou duas facetas cruzadas;
  reproduz Petersen (cap. 5) e Martinková e Hladká (4.5).
- **Profecia de Spearman-Brown e correções de correlação** (est_construcao.R):
  `spearman_brown()` e `corrigir_correlacao()` (atenuação; restrição de amplitude, caso 2
  de Thorndike), com os exemplos de Martinková e Hladká (4.2.5.1 e 2.4.3).
- **Análise de distratores** (est_construcao.R): `analise_distratores()`, igual ao
  ShinyItemAnalysis.
- **dMACS** (psicometria.R): `dmacs()`, forma fechada para indicadores contínuos, ΔMean
  de escala e aviso quando a métrica latente fixa torna a decomposição circular.

Funções ampliadas:

- `concordancia_juizes()`: AC1/AC2 de Gwet e alfa de Krippendorff com IC bootstrap, e
  aviso do paradoxo do kappa com os três casos da aula de limitações do kappa; o kappa
  incalculável não derruba mais a função.
- `confiabilidade(ic = TRUE)`: IC por bootstrap dos casos para ômega e alfa, reajustando
  o mesmo modelo.
- `tri_dif(metodo = "mh")`: Mantel-Haenszel com purificação e ΔMH do ETS, igual ao difR.
- `regressao_logistica(validar = TRUE)`: validação interna por bootstrap de otimismo (AUC,
  Brier, inclinação e intercepto de calibração).
- `qui_quadrado()` 2 × 2: RR, DR com IC de Newcombe e NNT com IC de Altman.
- `afe_numero_fatores(mrfa = TRUE)`: PA-MRFA de Timmerman e Lorenzo-Seva (EFA.MRFA), a
  preferida pelo capítulo 11 do Wiley Handbook.

Referências: 14 novas com DOI conferido (258 no total) e Thorndike (1949) sem DOI.
Notas novas: seleção automática de preditores (associacao.md), zero-inflada e riscos
competitivos (epidemiologia.md), defeitos dos pacotes de referência (fontes-e-decisoes.md).
Evals: `evals_v16.R` (46 testes).

## 1.5.0 (29/10/2026)

Revisão sistemática da pasta Psicometria-R (11 pastas, 1.934 arquivos, 655 aulas
transcritas e os scripts das aulas) contra a skill, aula por aula e script por script. As
lacunas com código nas aulas viraram funções testadas; o que não tem código nas aulas
(CAT, equacionamento) segue como rota fora do escopo testado. `.ee_versao` passa a
"1.5.0".

Funções novas:

- **Não paramétricos em delineamentos complexos** (est_naoparametricos.R, novo):
  `np_fatorial()` (rankFD: ATS, WTS, efeitos relativos; mctp), `np_longitudinal()`
  (nparLD: LD-F1, F1-LD-F1, F2-LD-F1; mctp.rm) e `manova_np()` (MANOVA.RM: MATS com
  bootstrap paramétrico, univariados com Holm).
- **MEE além da AFC** (est_mee.R, novo): `mee()` (modelo estrutural com efeitos
  definidos e IC bootstrap percentil), `mee_multigrupo()` (moderação categórica, com
  teste por caminho), `interacao_latente()` (indicadores-produto com centralização
  residual, inclinações simples), `curva_crescimento()` (linear ou de base latente; com
  resíduos iguais equivale ao modelo misto) e `invariancia_longitudinal()`.
- **TRI multidimensional** (est_tri.R): `tri_multidimensional()` (GRM, GPCM, 2PL e PCM;
  exploratório com rotação ou confirmatório; EM ou MHRM; confiabilidade empírica por
  dimensão; comparação com o unidimensional).
- **Meta-análise** (est_modelos.R): `meta_rede()` (netmeta: comum e aleatório,
  decomposição do Q, netsplit, P-score) e `extrair_media_dp()` (mediana para média e DP
  por Luo e Wan, EP e IC para DP, grupos combinados, DP da mudança).
- **Classes latentes**: `classes_distal()` (três passos com correção BCH).
- **Redes**: `rede_bayesiana()` (bnlearn, rede média de reamostras bootstrap).
- **Construção**: `confiabilidade_metades()` (Spearman-Brown, Guttman λ4, todas as divisões).
- **Epidemiologia**: `tendencia_joinpoint()` (APC e AAPC; pontos e estimativas iguais aos
  do Joinpoint 4.9.1 da aula) e `taxa_padronizada()` (igual ao epitools).
- **Vieses de resposta**: `controle_desejabilidade()` (quádruplas e fator de
  desejabilidade com sinal pelo polo) e `escolha_forcada()` (TRI thurstoniana com
  thurstonianIRT, confiabilidade empírica pela fórmula das aulas, reflexão de traço
  invertido).
- **Dados ausentes**: `ausentes()` (padrões e teste MCAR de Little, igual ao naniar).
- `classes_latentes()` e `perfis_latentes()` guardam as probabilidades posteriores e as
  linhas usadas (no objeto, fora do print), para o BCH.

Achados que viraram regra: o nparLD 2.3 removeu os argumentos das aulas; o `mctp` com
log-odds devolve o estimador em escala transformada (a skill usa Fisher); a ATS de dois
grupos do rankFD não é o Brunner-Munzel ao quadrado; o `pairwise()` saiu do netmeta 3;
os IC do Joinpoint no relato da aula não se reproduzem por nenhum dos dois métodos
(estimativas sim); a sintaxe da curva de segunda ordem das aulas difere do roteiro
comentado; a TRI thurstoniana com itens de mesmo sinal recupera mal os traços.

References novas: `naoparametricos-complexos.md`, `mee.md`. Atualizadas: tri, modelos,
classes-perfis, epidemiologia, vieses-resposta, redes-ega, construcao-instrumentos,
fontes-e-decisoes, figuras-dados, diferencas, poder-amostra, invariancia, afc-ordinal e
referencias-verificadas (24 DOIs novos, todos VERIFIED). Evals novos:
`evals_naoparametricos.R`, `evals_mee.R`, `evals_complementos.R`; blocos novos em
`evals_tri.R` e `evals_modelos.R`.

## 1.4.0 (29/10/2026)

Atualização pela pasta Psicometria-R do usuário (seis livros-base, 58 recortes por
módulo, mapa, auditoria CoVe e o curso Psicometria Online Academy completo, com
transcrições, sintaxes R e inputs de Mplus). Cada função nova foi conferida contra um
pacote de referência, um valor publicado ou a conta à mão. `.ee_versao` passa a "1.4.0".

Funções novas:

- **Número de fatores**: `afe_numero_fatores()` ganha Hull (CFI e CAF; Lorenzo-Seva
  et al., 2011), igual ao `EFAtools::HULL`, e MAP de Velicer (original e revisado),
  igual ao `psych::vss`, com a leitura do cap. 11 do Wiley Handbook; `n_teorico`.
- **Itens e fidedignidade** (est_construcao.R): `analise_classica_itens()`
  (dificuldade, RIT, RIR, ULI, alfa sem o item; Martinková & Hladká, 2023);
  `confiabilidade_teste_reteste()` (ICC(A,1) principal, EPM, mudança mínima detectável,
  Bland-Altman); `concordancia_juizes()` (Cohen, ponderado, Fleiss, ICC);
  `confiabilidade_ema()` (RkF, R1R, RkR, Rc; Shrout & Lane, 2012).
- **Bifator** (psicometria.R): `indices_bifator()` (ωH, ωHS, ECV, PUC, H, FD, I-ECV,
  ARPB), igual ao BifactorIndicesCalculator.
- **TRI** (est_tri.R, novo): `tri_ajustar()` (Rasch, 2PL, 3PL, PCM, RSM, GPCM, GRM;
  M2/M2*/C2, S-X², infit e outfit, Q3 relativo, contraste dos resíduos, precisão por
  faixa de θ), `tri_comparar()` e `tri_dif()` (mirt em duas etapas de âncoras;
  regressão logística purificada igual ao difR). Achado: `mirt::marginal_rxx()` ignora a
  variância latente estimada na família Rasch; a skill usa a fórmula com ela.
- **Poder para AFC** (est_poder.R): `poder_afc_rmsea()`, `n_aipe_rmsea()`,
  `poder_afc_aninhados()`, `n_aipe_parametro()`; reproduzem os números do cap. 5 do
  Wiley Handbook (com o N mínimo exato onde o MBESS aproxima).
- **ROC e viés preditivo** (est_associacao.R): `curva_roc()` com `prevalencia`,
  `custo_fn_fp`, VPP e VPN por Bayes, RV+ e RV−; `vies_preditivo()` (Cleary, step-down).
- **Vieses de resposta** (est_vieses.R, novo): `controle_aquiescencia()` (índice,
  recentralização, interceptos randômicos e MIMIC do curso).
- **Classes e perfis latentes** (est_latentes.R, novo): `classes_latentes()` (poLCA) e
  `perfis_latentes()` (mclust), iguais ao poLCA e ao tidyLPA.
- **Epidemiologia** (est_epidemio.R, novo): `razao_prevalencia()` (Poisson robusta) e
  `sobrevida()` (Kaplan-Meier, log-rank, Cox, Schoenfeld).
- **Meta-análise** (est_modelos.R): `meta_proporcao()` (igual ao `meta::metaprop`) e
  `meta_acuracia()` (igual ao `mada::reitsma`).
- `mudanca_confiavel()` ganha a ressalva do RCI como convenção e a alternativa pelo erro
  condicional da TRI.

Documentação: `tri.md` reescrito; novos `vieses-resposta.md`, `classes-perfis.md`,
`epidemiologia.md`, `figuras-dados.md` (preparação de dados, figuras APA, tradução de
Mplus para lavaan) e `script-reproduzivel.md` (seção movida do SKILL.md, que ficou mais
curto); atualizados construcao-instrumentos, afc-ordinal (bifator; ausentes pela
pasta), poder-amostra (seção 10), associacao, invariancia (viés de medida, DIF e viés
preditivo), fontes-e-decisoes (hierarquia com a pasta; erros do curso: MCAR de Little
com H0 invertida, NCT não significativo lido como igualdade, reteste pela Pearson),
pbat-mpfi (EMA, redes intraindividuais, bifator do MPFI-24) e modelos.
`referencias-verificadas.md`: 59 entradas novas com DOI (todas VERIFIED no
`verify_citations.py` 3.2.0, sem retratação) e 6 sem DOI.

Evals novos ou ampliados: `evals_tri.R`, `evals_vieses.R`,
`evals_latentes_epidemio.R`; ampliados construção, psicometria, poder, associação e
modelos. Conferências ao vivo com pacotes que o Mac pode não ter (EFAtools, irr,
BifactorIndicesCalculator, difR, eRm, tidyLPA, meta, mada, semTools) são puladas com
aviso quando o pacote falta; a conta à mão ou o valor publicado cobre o mesmo ponto.

## 1.3.0 (26/10/2026)

Três mudanças de conteúdo, cada uma conferida na fonte primária, e quatro bugs
achados no caminho. `.ee_versao` passa a "1.3.0".

- **Não significativo não quer dizer igual.** Todo teste de diferença e os de
  associação (correlações, ponto-bisserial, parcial, comparação de correlações,
  qui-quadrado) põem em `$avisos` uma nota quando o teste principal não é
  significativo, com a estimativa e o IC (Greenland et al., 2016, equívocos 4, 6
  e 8), e apontam o TOST quando ele existe para a estatística. Função nova
  `teste_equivalencia()` (est_diferencas.R): TOST independente (Welch ou Student),
  pareado, de uma amostra e de correlação, pelas eq. 3, 4, 6, 7 e 8 de Lakens
  (2017), com o teste usual ao lado e os quatro resultados combinados; limites
  obrigatórios, sem sugestão de valor; limites em d convertidos pelo DP da
  amostra, com o aviso de viés da documentação do TOSTER.
- **Escore somado.** Função nova `checar_escore_somado()` (psicometria.R): por
  escala ou subescala, congenérico, tau-equivalente e paralelo, com a comparação
  decisiva paralelo contra congenérico (McNeish & Wolf, 2020); diferença de χ²
  pelo `lavTestLRT()` com o método do estimador (Satorra & Bentler, 2001; Satorra,
  2000; refação por Satorra & Bentler, 2010, se negativa); caso ordinal com
  WLSMV theta (cargas iguais contra livres); 2 itens; ajuste, cargas e amplitude;
  frase PT e EN; recomendação de SEM ou escore fatorial quando rejeitado.
- **Rede não é fator.** `references/redes-ega.md`, seção 10: por que o ajuste
  não escolhe entre rede e fator (Kruis & Maris, 2016; van Bork et al., 2021) e
  quando os dados separam fator único de rede esparsa (teste PCL, citado e não
  implementado). `references/fontes-e-decisoes.md`, seção 7 nova: afirmação, o
  que ela exige e onde a skill confere.

CoVe direto na fonte (texto publicado ou documentação do pacote), com o que
mudou por causa dela:

- McNeish e Wolf (2020): o exemplo do artigo (HolzingerSwineford1939, x4 a x9
  arredondados) se reproduz no lavaan com ML, inclusive as cargas da Tabela 2 e
  a diferença χ²(9) = 40,80 do modelo de subescalas. O arredondamento tem de ser
  meio para cima; com o `round()` do R os números não batem. O χ²(19) = 5361,86
  do paralelo de um fator no texto é incompatível com o CFI e o RMSEA do próprio
  artigo; o lavaan dá 325,90, que reproduz os dois. O SRMR do artigo é o do
  Mplus. O artigo não usa o tau-equivalente como passo e trata os itens como
  contínuos: a sequência de três comparações e o caso ordinal ficaram declarados
  como da skill. O artigo testa as subescalas juntas e sem covariância; a função
  testa uma a uma e diz isso. A frase sobre R² = .970 ficou com o contexto do
  artigo (subescala Speeded, diferença entre escolas).
- Lakens (2017): eq. 3, 4, 6, 7 e 8, IC de 1 − 2α, relato do maior p e os quatro
  resultados conferidos no texto. Na eq. 7 os limites de uma amostra são
  relativos a mu; o `TOSTER::t_TOST` 0.8.6 pede os limites na escala da média.
  Os evals conferem os dois jeitos.
- Documentação do TOSTER 0.8.6: o aviso sobre limites padronizados ("SMD bounds
  will produce biased results") está no `t_TOST`.
- Greenland et al. (2016): numeração dos equívocos 4, 6 e 8 conferida no texto.
- Kruis e Maris (2016) e van Bork et al. (2021): equivalência e as três
  implicações do fator único conferidas no texto; o PCL ficou como não
  implementado.
- Lista verificada: 10 entradas novas com DOI (161 no total), todas com
  `is_retracted` falso no OpenAlex, e Croon (2002) em "Citadas sem DOI". Satorra
  e Bentler (2010) tem corrigenda de 2025 ligada no `updated-by`, registrada na
  nota de correções. Goertzen e Cribbie (2010) e Skrondal e Laake (2001), citados
  em texto de saída, entraram na lista. No lavaan 0.7.3 a função `fsr()` que o
  artigo cita para a correção de Croon não existe mais; o texto indica
  `lavPredict()`.

Bugs corrigidos:

- `tabela_ajuste()` com ML e FIML: a versão robusta saía com χ² NA (procurava o
  escalonado, que não existe com ML); agora traz o χ² padrão.
- Correlação numericamente zero (1e-17) saía com dezenas de casas na frase; abaixo
  de 1e-10 em módulo o valor sai como zero.
- O cabeçalho do `lavTestLRT()` usa aspas curvas em locale UTF-8 e retas em C; a
  leitura do método aceita as duas.
- `psych::polychoric` abre processos por fork (`mcmapply`); com BLAS OpenMP uma
  rodada de evals travou sem erro. A skill passa a pedir um núcleo nas chamadas da
  policórica (`.ee_policorica()`, mesmo resultado), e o `evals_construcao.R` fixa
  `mc.cores = 1` para as comparações diretas com o psych.

Evals (26/10/2026, R 4.5.3, com os bancos do JASP): diferenças 197/197 (17
novos), associação 134/134 (8 novos), psicometria 70/70 (17 novos, com as
âncoras do artigo), processos 62/62, modelos 122/122, poder 96/96, construção
59/59, redes 57/57. TOSTER entrou em `.ee_pacotes_evals`; sem ele, os testes de
equivalência com outra implementação são pulados.

## 1.2.0 (25/10/2026)

Script reproduzível para a rota do r-btw, em que o `btw_tool_run_r` roda o
código no ambiente global e a documentação não fala em registro do que foi
executado. Só o SKILL.md mudou; references e scripts seguem os da 1.1.0, e a
conferência do `.ee_versao` ficou escrita como "1.1.0" porque o card não guarda
a chave de versão do cabeçalho.

- IRON RULE 11: todo número do relato tem código salvo que o refaz.
- Seção "Script reproduzível": `registrar()` roda cada bloco e só o grava no
  `analise_<data>.R` depois que ele roda sem erro; a abertura copia as funções
  da skill para `psicometria-r/` ao lado dos dados; o bloco de leitura guarda o
  md5 do banco; o rótulo diz o tipo (descritivo, confirmatório, exploratório) e
  o destino (relatado ou não).
- Fechamento antes do Passo 3: bloco `numeros`, `sessionInfo()` como comentário
  no fim do script e reexecução em R limpo (`Rscript --no-init-file`); o pacote
  de resultados só sai com `status` 0 e `all.equal` TRUE ou com a divergência
  explicada no LOG.
- LOG com o caminho e o md5 do script e do banco e o resultado da reexecução;
  três itens novos no checklist.
- renv, Git, WORCS e o fluxo de Peikert e Brandmaier (2021) como opcionais para
  o mestrado e o OSF.

Revisão contra Sandve et al. (2013) e Obels et al. (2020) e duas passadas de
CoVe mudaram, antes da publicação: cópia das funções junto do script, md5 do
banco, arquivo bruto sem edição à mão, rótulo por tipo e destino, reanálise de
parecerista com novo fechamento, figura em bloco registrado, README e codebook
no depósito do OSF, Git como opcional e o limite do que a reexecução garante; o
Quarto deixou de ser atribuído a Peikert e Brandmaier (2021), a nota de teste
passou a descrever cada teste em separado, a macro do SPSS passou a ser de um
artigo da amostra de Obels, e a seção passou a dizer que aplica parte das dez
regras de Sandve. Os seis DOIs saíram VERIFIED no `verify_citations.py`.

Teste do template em 25/10/2026, no container (R 4.3.3) e no Mac pelo r-btw (R
4.6.1, macOS 26.5.2), só com dados simulados. No Mac, o diretório de trabalho
volta para `/` a cada chamada, e o bloco de leitura parava no md5; `registrar()`
passou a fixar a pasta com `setwd()` e `on.exit()`, e o `file.copy()` da
abertura ficou silencioso. Depois da correção, a reexecução e a segunda sessão
fecharam com `status` 0 e `all.equal` TRUE, e a trava do md5 parou o banco
alterado.

Revisão final no Notion (26/10/2026): a frase de Obels sobre descritiva,
confirmatória e exploratória ficou restrita aos scripts de SPSS, com o "muitas
vezes" do artigo, e ganhou a recomendação deles de dizer qual código é de
análise que não entra no artigo; o Fluxo passou a aceitar, como o fechamento e o
checklist, a divergência explicada no LOG.

## 1.1.0 (11/10/2026)

Expansão para a estatística mais usada em psicologia, construída a partir das
fontes do usuário (livro Biblioteca Análise de Dados, 7 volumes, JASP 0.95.1;
notas do curso Psicometria Online Academy, 94 páginas lidas) e conferida com a
literatura (151 entradas com DOI verificado no Crossref) e com o código do
JASP e do G*Power.

- Entrada única `scripts/carregar.R`; base comum `est_comum.R` (objeto
  `pr_resultado`, frase APA, descritivas com assimetria e curtose do JASP,
  leitura de normalidade do livro, atípicos por z e MAD, índice de mudança
  confiável).
- Módulos novos: `est_diferencas.R` (t, ANOVA de um fator, medidas repetidas,
  fatorial, mista, ANCOVA, MANOVA, Mann-Whitney, Wilcoxon, Kruskal-Wallis,
  Friedman), `est_associacao.R` (correlações, parcial, ponto-bisserial,
  comparação de correlações, regressão linear e hierárquica com diagnósticos,
  logística, ROC, qui-quadrado, McNemar, Q de Cochran), `est_processos.R`
  (moderação com Johnson-Neyman, mediação simples, paralela e serial, mediação
  moderada nos modelos 7, 14 e 15), `est_modelos.R` (GLM, GEE com QIC,
  multinível linear e logístico, meta-análise), `est_poder.R` (tamanho amostral
  com os números do G*Power), `est_construcao.R` (validade de conteúdo, AFE) e
  `est_redes.R` (rede com os padrões do JASP, estabilidade, ponte,
  previsibilidade).
- Referências novas por módulo e `references/fontes-e-decisoes.md` (hierarquia
  das fontes, rotas clássica e alternativa, erros das fontes não reproduzidos).
- `redes-ega.md` alinhado ao livro: γ, ausentes e tamanho amostral da
  estimação declarados, IC de aresta como descritivo, CS por índice, ponte,
  previsibilidade, frase de alcance causal e checklist de relato.
- `.pr_num` arredonda meio para cima, como o JASP (o `round()` do R dava 155,04
  para 155,045); `versoes_pacotes()` inclui os pacotes dos módulos.
- Evals: 180 + 126 + 62 + 122 + 96 + 59 + 57 testes novos (diferenças,
  associação, processos, modelos, poder, construção, redes), dos quais 178 usam os
  bancos da biblioteca do JASP (diferenças 82, associação 54, processos 18, modelos 4,
  construção 2, redes 18) e são pulados sem eles; quase todos reproduzem números do
  livro, e os demais testam comportamento (os 2 de construção usam o banco G Factor,
  que o livro não usa). Mais os 53 da psicometria. Todos passam no R 4.6.1
  (14/10/2026), com e sem os bancos.
- Documentação: `relato-apa.md` ganhou o checklist da rede e o item de análises
  estatísticas no pacote de resultados; `afc-ordinal.md`, `tri.md` e o cabeçalho
  de `psicometria.R` passam a carregar por `carregar.R`; IRON RULES 4 a 6 novas
  (fonte do usuário decide a rota, pressuposto não é liga-desliga, tamanho de
  efeito com IC); índice no topo dos arquivos de referência longos.

Revisão independente (CoVe com cinco revisores, 11/10/2026) mudou, antes da
publicação:

- diferenças: F, gl e p da MANOVA; Wilcoxon pareado com zeros; pós-hoc da
  fatorial; Henze-Zirkler nos resíduos; célula de texto vazia vira ausente; o
  aviso que repetia o erro do livro sobre o g de Hedges e desbalanceamento saiu;
  com 3 ou 4 comparações, `$decisoes` lembra que o livro aceita Bonferroni;
- associação: 3 casas como os parágrafos-modelo do livro; β com zero à
  esquerda; IC de V e φ; aderência com categoria de contagem zero;
  categoria de evento única na ROC;
- modelos: preditor numérico com poucos valores avisa e aceita `fatores =`;
  intervalo de predição com t(k − 2) (Riley et al., 2011); teste de Peters com os
  pesos do artigo; `centralizar_multinivel()` invisível e com padrão na média
  geral, como o curso; inclinação aleatória de fator recusada com mensagem; LRT
  com a mesma aproximação no multinível logístico com nAGQ > 1;
- poder: N da logística com R² dos outros preditores pela regra do exemplo do
  manual do G*Power (395, não 375), com aviso; alocação desigual arredondada como
  o G*Power (teto de n1* e de razão · n1*); α na frase com as casas do valor usado;
  com `metodo = "fisher"` a frase não cita o G*Power; α e poder sem zero à
  esquerda;
- construção: dependência linear exata entre itens (escore total, item
  duplicado) para com mensagem; frase da adequação segue KMO e Bartlett; caso
  Heywood dito na frase da AFE;
- comum: descritivas com erro-padrão da assimetria e da curtose e supressão, no
  print, de grupo com n < 5.

Verificação final (CoVe com quatro revisores, 14/10/2026) mudou, antes da
publicação:

- privacidade: a supressão de grupo com n < 5 vale também para tabelas extras com
  coluna n (médias marginais do GLM e do GEE, codificação de grupos do
  ponto-bisserial) e para as frases APA (t, Mann-Whitney, Wilcoxon, ANCOVA,
  ponto-bisserial), que omitem M, DP e Md e avisam; `ler_normalidade()` não lê
  grupo com n < 5; os avisos dizem que estatísticas de teste com um grupo minúsculo
  ainda permitem reconstruir valores, e recomendam agrupar categorias;
- grupo gravado como número (3 a 7 valores inteiros) avisa também na regressão
  linear, hierárquica e logística (que ganham `fatores =`), na moderação e na
  mediação; a regra é a mesma em todos os módulos (0/1 não avisa);
- diferenças: na rota Welch, o ω² aproximado pelo F de Welch sai ao lado do
  clássico, com aviso quando os dois contradizem; teste unilateral diz H1 na frase;
  VD constante e Levene degenerado param ou avisam com mensagem clara; separador
  e vírgula decimal coerentes em `$decisoes`;
- associação: magnitude no idioma pedido; OR com correção de Haldane dito na frase;
  aderência com o α de `nivel`; R² ajustado, ΔR² e f² sem IC registrados como não
  calculados; categoria com menos de 2 casos para com mensagem;
- modelos: `nivel` chega a todos os rótulos de IC (meta-análise, comparações do GLM
  e do GEE); no GEE, sujeito que determina sozinho um parâmetro (nível com um só
  sujeito) deixa a correção de Mancl e DeRouen indefinida, o que agora vira aviso e
  EP robusto, em vez de erro (no R do Mac a função parava);
- poder: conferência do N com o mesmo arredondamento da frase; Hsieh com preditor
  binário; desenho mínimo com efeito grande devolvido com aviso; poder na frase com
  as casas do valor usado; bloco do G*Power com 6 casas e rótulo "Effect size d" no
  t de uma amostra; dois exemplos de compromisso de Faul et al. (2009) nos evals;
- nível de fator com célula vazia (texto "") não muda mais a ordem dos níveis.

Achado que virou regra: o R 4.6 mudou o `wilcox.test` (p exato condicional com
empates; zeros mantidos nos postos do pareado); a skill fixa a regra do R 4.5.1,
a mesma do JASP 0.95.1.

## 1.0.0 (10/10/2026)

Criação. Funções em `scripts/psicometria.R` (auditoria de itens, inclusive
por grupo, sugestão de estimador, AFC, tabela de ajuste com versões
escalonada e robusta, confiabilidade por fator e de compostos, validade
discriminante, invariância ordinal e contínua, relato APA), evals em
`scripts/evals_simulados.R` (53 testes, 53 passando no R 4.6.1 do Mac),
references por análise e lista de 35 referências conferidas com
`verify_citations.py` v3.2.

Revisão independente antes da publicação (10-11/10/2026) mudou: a regra de
estimador com 5 a 7 categorias (MLR principal com amostra pequena e limiares
simétricos, como recomendam Rhemtulla et al., 2012); `validade_discriminante()`
só entre fatores de primeira ordem, alinhada por nome e com NA para um fator
só; `invariancia()` recusa modelo de ordem superior e não inventa p em passo
sem grau de liberdade; ausentes iguais na AFC e na invariância; formato APA
com vírgula decimal; evals com `ordered_center = FALSE`, diferença real de
média latente, fatores a .90, itens de 3 categorias e itens sem DIF; ferramentas
`env` do btw fora do servidor r-btw; PBAT como instrumento do mestrado.

Achados dos evals que viraram regra de relato:

- com a identificação de Wu e Estabrook (2016), deslocamento uniforme dos
  limiares de um item passa no passo dos limiares e cai no dos interceptos;
- no semTools 0.5.10, `compRelSEM(higher = )` está obsoleto, e num modelo de
  segunda ordem o ômega por subescala é relativo ao fator geral; o ômega do
  escore global sai com `W=` e muda de sentido conforme o modelo (hierárquico
  na segunda ordem, total nos fatores correlacionados).

Âncora do corte de estabilidade de item do bootEGA (≥ .70; Christensen &
Golino, 2021) resolvida pela página do artigo na MDPI; conferir no PDF antes
de citação literal.
