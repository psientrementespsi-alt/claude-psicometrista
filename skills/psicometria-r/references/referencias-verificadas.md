# Referências verificadas da skill

## Sumário

- 53 (linha 74)
- 105 (linha 126)
- 117 (linha 138)
- 140 (linha 161)
- 162 (linha 183)
- 195 (linha 216)
- 227 (linha 248)
- 246 (linha 267)
- 291 (linha 312)
- 309 (linha 330)
- 348 (linha 369)
- 356 (linha 377)
- 387 (linha 408)
- 410 (linha 431)
- 433 (linha 454)
- 443 (linha 464)
- 452 (linha 473)
- 473 (linha 494)

Índice: Psicometria e AFC · Invariância · TRI · Redes e EGA · Diferenças entre grupos · Associação e regressão · Processos · Modelos avançados · Poder e tamanho amostral · Construção de instrumentos · Relato · v1.5: módulos novos · Citadas sem DOI · Não verificadas

Conferidas em 11/10/2026 (Faul et al., 2009, em 14/10/2026; as entradas da v1.3,
em 26/10/2026, também com `is_retracted` do OpenAlex, falso em todas). Cada DOI foi
consultado com `curl` no endpoint `https://api.crossref.org/works/<DOI>` (JSON
bruto, sem intermediário); título, primeiro autor e ano bateram com a citação nas
161 entradas com DOI da v1.3 (uma delas é a correção de Delacre et al., 2022), 59 da v1.4, 24 da v1.5, 14 da v1.6 e 14 da v1.7 (29/10/2026, mesma rotina; 28/10/2026 para as da v1.4, `verify_citations.py` 3.2.0: todas VERIFIED, nenhuma retratação; correções registradas em McGraw & Wong, 1996, Koo & Li, 2016, e Rodriguez et al., 2016b), mais 5 DOIs
de erratas e correções citados só nesta nota, também consultados.
Autores, volume, número e páginas vêm do Crossref; quando o Crossref só dá a
primeira página (registros do JSTOR) ou lista só o primeiro autor (Peters et al.,
2006), o restante vem da pesquisa de literatura de 11/10/2026
(`fontes/pesquisa_literatura.md`). Páginas dos artigos do JSS vêm da meta
`citation_firstpage`/`citation_lastpage` do site do periódico (conferência de
10/10/2026, com `verify_citations.py` v3.2, mantida). Livros sem DOI: ISBN
conferido no Open Library quando possível.

**Campo `updated-by` lido em todos os DOIs.** Nenhuma retratação, retirada
(withdrawal), remoção ou manifestação de preocupação (expression of concern).
Correções e erratas registradas, que não invalidam o trabalho mas pedem leitura
antes de citar números:

- Samejima (1969): errata em *Psychometrika, 35*(1), 139 (1970),
  https://doi.org/10.1007/BF02290599 (ligada no `updated-by`).
- Grégoire et al. (2020): errata coletiva do *Journal of Contextual Behavioral
  Science, 21*, 223–226 (2021), https://doi.org/10.1016/j.jcbs.2021.03.005, que
  cobre dezenas de artigos dos volumes 15 a 18 (conteúdo não lido; o formato
  sugere correção editorial em lote, não do resultado).
- Funder e Ozer (2019): corrigendum em *Advances in Methods and Practices in
  Psychological Science, 3*(4), 509 (2020), https://doi.org/10.1177/2515245920979282.
- Satorra e Bentler (2010): corrigenda em *Psychometrika, 90*(1), 415 (2025),
  https://doi.org/10.1017/psy.2024.12 (ligada no `updated-by`; conteúdo não lido). A
  skill usa o método de 2010 só como refação automática do `lavTestLRT()` quando a
  diferença escalonada de 2001 sai negativa; ler a corrigenda antes de citar fórmula.
- Correções publicadas mas **não** ligadas no `updated-by`: Delacre et al. (2017),
  artigo de correção https://doi.org/10.5334/irsp.661 (2022); Stram e Lee (1994),
  "Corrections" em *Biometrics, 51*, 1196 (1995), https://doi.org/10.2307/2533038.

Cada trabalho aparece uma vez, no módulo em que é mais usado; a coluna "Uso na
skill" diz os arquivos de `references/` que o citam. Tipo: artigo (periódico com
revisão por pares), livro, capítulo, relatório técnico, correção. Nenhum preprint
é citado. A coluna "Âncora" diz de onde saiu o que a skill afirma sobre cada
trabalho: resumo (R), texto completo (T), documentação de pacote (D), fonte
secundária (Sec.) ou só título e metadados (Tít.). Afirmação apoiada só no
título não deve sustentar alegação específica num manuscrito sem ler o texto.
Ano APA = ano do volume impresso; quando o Crossref registra publicação online
antes, isso vai na âncora só se mudar a leitura.

Entradas marcadas "(epônimo)" são trabalhos que o arquivo cita só pelo nome do
método, sem ano; foram conferidas para o coordenador decidir se entram no texto.

## Psicometria e AFC

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Brosseau-Liard, P. E., & Savalei, V. (2014). Adjusting incremental fit indices for nonnormality. *Multivariate Behavioral Research, 49*(5), 460–470. https://doi.org/10.1080/00273171.2014.933697 | artigo | versão robusta dos índices; citado só pela documentação do `fitMeasures()` (lista anterior) | D; Tít. |
| Dunn, T. J., Baguley, T., & Brunsden, V. (2014). From alpha to omega: A practical solution to the pervasive problem of internal consistency estimation. *British Journal of Psychology, 105*(3), 399–412. https://doi.org/10.1111/bjop.12046 | artigo | alfa e ômega lado a lado (afc-ordinal) | R |
| Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. *Journal of Marketing Research, 18*(1), 39–50. https://doi.org/10.1177/002224378101800104 | artigo (epônimo) | critério de Fornell-Larcker (afc-ordinal, construcao-instrumentos) | Tít. |
| Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using structural equation modeling: An alternative to coefficient alpha. *Psychometrika, 74*(1), 155–167. https://doi.org/10.1007/s11336-008-9099-3 | artigo | ômega do escore observado com itens ordinais (afc-ordinal, construcao-instrumentos, relato-apa) | R |
| Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for assessing discriminant validity in variance-based structural equation modeling. *Journal of the Academy of Marketing Science, 43*(1), 115–135. https://doi.org/10.1007/s11747-014-0403-8 | artigo | HTMT (afc-ordinal) | R |
| Hu, L.-t., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure analysis: Conventional criteria versus new alternatives. *Structural Equation Modeling: A Multidisciplinary Journal, 6*(1), 1–55. https://doi.org/10.1080/10705519909540118 | artigo | cortes convencionais, derivados para ML (afc-ordinal, construcao-instrumentos, relato-apa) | R |
| Li, C.-H. (2016). Confirmatory factor analysis with ordinal data: Comparing robust maximum likelihood and diagonally weighted least squares. *Behavior Research Methods, 48*(3), 936–949. https://doi.org/10.3758/s13428-015-0619-7 | artigo | WLSMV contra MLR (afc-ordinal) | R (PubMed) |
| McNeish, D. (2018). Thanks coefficient alpha, we'll take it from here. *Psychological Methods, 23*(3), 412–433. https://doi.org/10.1037/met0000144 | artigo | limites do alfa (afc-ordinal) | R |
| McNeish, D., & Wolf, M. G. (2020). Thinking twice about sum scores. *Behavior Research Methods, 52*(6), 2287–2305. https://doi.org/10.3758/s13428-020-01398-0 | artigo | soma dos itens = escore do modelo paralelo; "bait-and-switch"; subescalas sem covariância; exemplo reproduzido nos evals (afc-ordinal, fontes-e-decisoes, `checar_escore_somado()`) | T (texto publicado na Springer e preprint); exemplo com HolzingerSwineford1939 reproduzido no lavaan, exceto o χ²(19) = 5361,86 do texto, incompatível com o CFI e o RMSEA do próprio artigo |
| McNeish, D., & Wolf, M. G. (2023). Dynamic fit index cutoffs for confirmatory factor analysis models. *Psychological Methods, 28*(1), 61–88. https://doi.org/10.1037/met0000425 | artigo | cortes dinâmicos (afc-ordinal) | R |
| Rhemtulla, M., Brosseau-Liard, P. É., & Savalei, V. (2012). When can categorical variables be treated as continuous? A comparison of robust continuous and categorical SEM estimation methods under suboptimal conditions. *Psychological Methods, 17*(3), 354–373. https://doi.org/10.1037/a0029315 | artigo | regra do estimador (afc-ordinal, pbat-mpfi, relato-apa) | R |
| Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. *Journal of Statistical Software, 48*(2), 1–36. https://doi.org/10.18637/jss.v048.i02 | artigo (software) | pacote de AFC; referência de software (relato-apa) | Tít. |
| Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample analysis of moment structures. In R. D. H. Heijmans, D. S. G. Pollock, & A. Satorra (Eds.), *Innovations in multivariate statistical analysis* (pp. 233–247). Springer US. https://doi.org/10.1007/978-1-4615-4603-0_17 | capítulo | diferença de χ² escalonada e deslocada com WLSMV, método `satorra.2000` do `lavTestLRT()` (afc-ordinal) | D; Tít. |
| Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test statistic for moment structure analysis. *Psychometrika, 66*(4), 507–514. https://doi.org/10.1007/BF02296192 | artigo | diferença de χ² escalonada com MLR, método `satorra.bentler.2001` do `lavTestLRT()` (afc-ordinal) | D; Tít. |
| Satorra, A., & Bentler, P. M. (2010). Ensuring positiveness of the scaled difference chi-square test statistic. *Psychometrika, 75*(2), 243–248. https://doi.org/10.1007/s11336-009-9135-y | artigo | refação quando a diferença escalonada sai negativa (afc-ordinal) | D; Tít.; corrigenda de 2025 na nota acima |
| Savalei, V. (2021). Improving fit indices in structural equation modeling with categorical data. *Multivariate Behavioral Research, 56*(3), 390–407. https://doi.org/10.1080/00273171.2020.1717922 | artigo | índices com dado categórico (afc-ordinal, invariancia) | R |
| Skrondal, A., & Laake, P. (2001). Regression among factor scores. *Psychometrika, 66*(4), 563–575. https://doi.org/10.1007/BF02296196 | artigo | escore de Bartlett quando o fator é desfecho, regressão quando é preditor (afc-ordinal, `checar_escore_somado()`) | Sec. (McNeish & Wolf, 2020); Tít. |
| Xia, Y., & Yang, Y. (2019). RMSEA, CFI, and TLI in structural equation modeling with ordered categorical data: The story they tell depends on the estimation methods. *Behavior Research Methods, 51*(1), 409–428. https://doi.org/10.3758/s13428-018-1055-2 | artigo | cortes com DWLS (afc-ordinal, construcao-instrumentos, relato-apa) | R (PubMed) |

### Instrumentos do usuário (pbat-mpfi.md)

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Ciarrochi, J., Sahdra, B., Hofmann, S. G., & Hayes, S. C. (2022). Developing an item pool to assess processes of change in psychological interventions: The Process-Based Assessment Tool (PBAT). *Journal of Contextual Behavioral Science, 23*, 200–213. https://doi.org/10.1016/j.jcbs.2022.02.001 | artigo | PBAT original (pbat-mpfi) | T (pp. 202-209) |
| Grégoire, S., Gagnon, J., Lachance, L., Shankland, R., Dionne, F., Kotsou, I., Monestès, J.-L., Rolffs, J. L., & Rogge, R. D. (2020). Validation of the English and French versions of the Multidimensional Psychological Flexibility Inventory Short Form (MPFI-24). *Journal of Contextual Behavioral Science, 18*, 99–110. https://doi.org/10.1016/j.jcbs.2020.06.004 | artigo | MPFI-24 (pbat-mpfi) | R (página da ACBS); errata coletiva registrada |
| Jo, D., Im, S., Lee, E., Sundström, F., & Yang, E. (2025). From diagnostic criteria to change processes: Validation of the Korean process-based assessment tool using a machine learning approach. *Cognitive Therapy and Research*. Publicação antecipada online. https://doi.org/10.1007/s10608-025-10693-z | artigo | adaptação coreana do PBAT (pbat-mpfi) | Tít.; Crossref sem volume em 11/10/2026 |
| Larsson, A. B., & Sundström, F. T. A. (2024). Optimizing process-based therapy: A Boruta-driven approach to identifying key clinical outcomes in a Swedish validation of the process-based assessment tool (PBAT). *Journal of Contextual Behavioral Science, 33*, Article 100794. https://doi.org/10.1016/j.jcbs.2024.100794 | artigo | adaptação sueca do PBAT (pbat-mpfi) | Tít. |
| Rolffs, J. L., Rogge, R. D., & Wilson, K. G. (2018). Disentangling components of flexibility via the hexaflex model: Development and validation of the Multidimensional Psychological Flexibility Inventory (MPFI). *Assessment, 25*(4), 458–482. https://doi.org/10.1177/1073191116645905 | artigo | MPFI completo (pbat-mpfi) | Tít.; estrutura pela página da ACBS |
| Sanford, B. T., Ciarrochi, J., Hofmann, S. G., Chin, F., Gates, K. M., & Hayes, S. C. (2022). Toward empirical process-based case conceptualization: An idionomic network examination of the process-based assessment tool. *Journal of Contextual Behavioral Science, 25*, 10–25. https://doi.org/10.1016/j.jcbs.2022.05.006 | artigo | PBAT intraindividual (pbat-mpfi) | R (página da ACBS) |

A adaptação romena do PBAT (Popa et al., 2025) está em "Redes e EGA"; as
diretrizes da International Test Commission (2018), em "Construção de
instrumentos".
| Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016a). Evaluating bifactor models: Calculating and interpreting statistical indices. *Psychological Methods, 21*(2), 137–150. https://doi.org/10.1037/met0000045 | artigo | `indices_bifator()`: fórmulas de ωH, ωHS, ECV, PUC, H e FD; réguas ωH ≥ ,80 e ECV/PUC > ,70 (afc-ordinal) | Sec.; D (BifactorIndicesCalculator) |
| Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016b). Applying bifactor statistical indices in the evaluation of psychological measures. *Journal of Personality Assessment, 98*(3), 223–237. https://doi.org/10.1080/00223891.2015.1089249 | artigo | aplicação dos índices do bifator (afc-ordinal); correção registrada no `updated-by` | Tít.; Sec. |
| Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013). Multidimensionality and structural coefficient bias in structural equation modeling: A bifactor perspective. *Educational and Psychological Measurement, 73*(1), 5–26. https://doi.org/10.1177/0013164412449831 | artigo | régua PUC < ,80 com ECV > ,60 e ωH > ,70 (afc-ordinal) | Sec. (Reise et al., 2018, T) |
| Reise, S. P., Bonifay, W., & Haviland, M. G. (2018). Bifactor modelling and the evaluation of scale scores. In P. Irwing, T. Booth, & D. J. Hughes (Eds.), *The Wiley handbook of psychometric testing* (pp. 675–707). Wiley. https://doi.org/10.1002/9781118489772.ch22 | capítulo | eq. 22.3 e 22.4 (ωH e ωHS), leitura de ECV e PUC (afc-ordinal, pbat-mpfi) | T (pasta Psicometria-R) |
| Hancock, G. R., & Mueller, R. O. (2001): ver "Citadas sem DOI". | | | |
| Muthén, B., Kaplan, D., & Hollis, M. (1987). On structural equation modeling with data that are not missing completely at random. *Psychometrika, 52*(3), 431–462. https://doi.org/10.1007/BF02294365 | artigo | régua do ARPB < 10–15% citada por Rodriguez et al. (afc-ordinal) | Sec. |
| Ferrando, P. J., & Lorenzo-Seva, U. (2018). Assessing the quality and appropriateness of factor solutions and factor score estimates in exploratory item factor analysis. *Educational and Psychological Measurement, 78*(5), 762–780. https://doi.org/10.1177/0013164417719308 | artigo | UniCo, ECV > ,85 e MIREAL do FACTOR, usados no curso (afc-ordinal) | Sec. (curso) |
| Enders, C. K., & Baraldi, A. N. (2018). Missing data handling methods. In P. Irwing, T. Booth, & D. J. Hughes (Eds.), *The Wiley handbook of psychometric testing* (pp. 139–185). Wiley. https://doi.org/10.1002/9781118489772.ch6 | capítulo | ausentes: pairwise com WLSMV pressupõe MCAR; FIML e imputação múltipla (afc-ordinal) | T (pasta Psicometria-R) |
| Maydeu-Olivares, A., & Coffman, D. L. (2006). Random intercept item factor analysis. *Psychological Methods, 11*(4), 344–362. https://doi.org/10.1037/1082-989X.11.4.344 | artigo | `controle_aquiescencia()`: modelo de interceptos randômicos (vieses-resposta) | Tít.; curso |
| Primi, R., Santos, D., De Fruyt, F., & John, O. P. (2019). Comparison of classical and modern methods for measuring and correcting for acquiescence. *British Journal of Mathematical and Statistical Psychology, 72*(3), 447–465. https://doi.org/10.1111/bmsp.12168 | artigo | índice e recentralização ensinados no curso com as funções do LabAPE (vieses-resposta) | Tít.; curso |
| Soto, C. J., John, O. P., Gosling, S. D., & Potter, J. (2008). The developmental psychometrics of Big Five self-reports: Acquiescence, factor structure, coherence, and differentiation from ages 10 to 20. *Journal of Personality and Social Psychology, 94*(4), 718–737. https://doi.org/10.1037/0022-3514.94.4.718 | artigo | aquiescência maior em crianças e adolescentes, citado no curso (vieses-resposta) | Sec. (curso) |
| Brown, A., & Maydeu-Olivares, A. (2011). Item response modeling of forced-choice questionnaires. *Educational and Psychological Measurement, 71*(3), 460–502. https://doi.org/10.1177/0013164410375112 | artigo | TRI thurstoniana para escolha forçada, roteiro do curso (vieses-resposta) | Tít.; curso |
| Bürkner, P.-C. (2019). thurstonianIRT: Thurstonian IRT models in R. *Journal of Open Source Software, 4*(42), Article 1662. https://doi.org/10.21105/joss.01662 | artigo | pacote usado no curso para a TRI thurstoniana (vieses-resposta) | D |
| Bäckström, M., Björklund, F., & Larsson, M. R. (2009). Five-factor inventories have a major general factor related to social desirability which can be reduced by framing items neutrally. *Journal of Research in Personality, 43*(3), 335–344. https://doi.org/10.1016/j.jrp.2008.12.013 | artigo | escala neutralizada contra desejabilidade social (vieses-resposta) | Tít.; curso |

## Invariância

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of measurement invariance. *Structural Equation Modeling: A Multidisciplinary Journal, 14*(3), 464–504. https://doi.org/10.1080/10705510701301834 | artigo | ΔRMSEA e ΔSRMR (invariancia) | R; valores via Putnick e Bornstein (2016) |
| Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes for testing measurement invariance. *Structural Equation Modeling: A Multidisciplinary Journal, 9*(2), 233–255. https://doi.org/10.1207/S15328007SEM0902_5 | artigo | ΔCFI (invariancia) | R |
| Putnick, D. L., & Bornstein, M. H. (2016). Measurement invariance conventions and reporting: The state of the art and future directions for psychological research. *Developmental Review, 41*, 71–90. https://doi.org/10.1016/j.dr.2016.06.004 | artigo | critérios e invariância parcial (invariancia) | T (PubMed Central) |
| Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor analysis models of different levels of invariance for ordered categorical outcomes. *Psychometrika, 81*(4), 1014–1045. https://doi.org/10.1007/s11336-016-9506-0 | artigo | invariância ordinal (invariancia, pbat-mpfi) | R + D (semTools) |
| Cleary, T. A. (1968). Test bias: Prediction of grades of Negro and white students in integrated colleges. *Journal of Educational Measurement, 5*(2), 115–124. https://doi.org/10.1111/j.1745-3984.1968.tb00613.x | artigo | definição de viés preditivo, `vies_preditivo()` (invariancia, associacao) | Sec. (Petersen, 2024, T) |
| Lautenschlager, G. J., & Mendoza, J. L. (1986). A step-down hierarchical multiple regression analysis for examining hypotheses about test bias in prediction. *Applied Psychological Measurement, 10*(2), 133–139. https://doi.org/10.1177/014662168601000202 | artigo | sequência ômnibus, inclinações, interceptos de `vies_preditivo()` | Tít. |
| Aguinis, H., Culpepper, S. A., & Pierce, C. A. (2010). Revival of test bias research in preemployment testing. *Journal of Applied Psychology, 95*(4), 648–680. https://doi.org/10.1037/a0018714 | artigo | poder baixo para diferença de inclinações (invariancia, associacao) | Sec. (Petersen, 2024, T) |

## TRI

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Cai, L., & Hansen, M. (2013). Limited-information goodness-of-fit testing of hierarchical item factor models. *British Journal of Mathematical and Statistical Psychology, 66*(2), 245–276. https://doi.org/10.1111/j.2044-8317.2012.02050.x | artigo | família M2 no mirt; citado só pela documentação do `M2` (lista anterior) | D; Tít. |
| Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment. *Journal of Statistical Software, 48*(6), 1–29. https://doi.org/10.18637/jss.v048.i06 | artigo (software) | pacote de TRI (tri, relato-apa) | Tít. |
| Kang, T., & Chen, T. T. (2007). *An investigation of the performance of the generalized S-X² item-fit index for polytomous IRT models* (ACT Research Report Series). ACT. https://doi.org/10.1037/e428112008-001 | relatório técnico | S-X² politômico, como citado pela documentação do `itemfit` (tri) | D; Tít. (registro do PsycEXTRA no Crossref) |
| Kang, T., & Chen, T. T. (2008). Performance of the generalized S-X² item fit index for polytomous IRT models. *Journal of Educational Measurement, 45*(4), 391–406. https://doi.org/10.1111/j.1745-3984.2008.00071.x | artigo | versão revisada por pares do relatório de 2007; preferir esta num manuscrito | Tít. |
| Maydeu-Olivares, A., & Joe, H. (2006). Limited information goodness-of-fit testing in multidimensional contingency tables. *Psychometrika, 71*(4), 713–732. https://doi.org/10.1007/s11336-005-1295-9 | artigo | família M2 (tri) | R |
| Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for dichotomous item response theory models. *Applied Psychological Measurement, 24*(1), 50–64. https://doi.org/10.1177/01466216000241003 | artigo | S-X² (tri) | R |
| Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. *Psychometrika, 34*(S1), 1–97. https://doi.org/10.1007/BF03372160 | artigo (monografia do periódico) | modelo de resposta graduada (tri) | Tít.; errata de 1970 registrada |
| Masters, G. N. (1982). A Rasch model for partial credit scoring. *Psychometrika, 47*(2), 149–174. https://doi.org/10.1007/BF02296272 | artigo | modelo de crédito parcial (tri) | Tít. |
| Andrich, D. (1978). A rating formulation for ordered response categories. *Psychometrika, 43*(4), 561–573. https://doi.org/10.1007/BF02293814 | artigo | modelo de escala de avaliação (tri) | Tít. |
| Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. *Applied Psychological Measurement, 16*(2), 159–176. https://doi.org/10.1177/014662169201600206 | artigo | GPCM (tri) | Tít. |
| Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase, M. D. (1984). Technical guidelines for assessing computerized adaptive tests. *Journal of Educational Measurement, 21*(4), 347–360. https://doi.org/10.1111/j.1745-3984.1984.tb01039.x | artigo | confiabilidade marginal σ²/(σ² + E[EP²]) de `tri_ajustar()` (tri) | Sec. |
| Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values for Yen's Q3: Identification of local dependence in the Rasch model using residual correlations. *Applied Psychological Measurement, 41*(3), 178–194. https://doi.org/10.1177/0146621616677520 | artigo | Q3 relativo (Q3 menos a média) e limite ,20 (tri) | R |
| Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with polychotomous item response models and standardized indices. *British Journal of Mathematical and Statistical Psychology, 38*(1), 67–86. https://doi.org/10.1111/j.2044-8317.1985.tb00817.x | artigo | Zh de ajuste de pessoas, só agregado (tri) | Tít.; D (mirt) |
| Swaminathan, H., & Rogers, H. J. (1990). Detecting differential item functioning using logistic regression procedures. *Journal of Educational Measurement, 27*(4), 361–370. https://doi.org/10.1111/j.1745-3984.1990.tb00754.x | artigo | DIF por regressão logística, `tri_dif(metodo = "logistica")` (tri) | Tít.; D (difR) |
| Jodoin, M. G., & Gierl, M. J. (2001). Evaluating Type I error and power rates using an effect size measure with the logistic regression procedure for DIF detection. *Applied Measurement in Education, 14*(4), 329–349. https://doi.org/10.1207/S15324818AME1404_2 | artigo | classes A/B/C pelo ΔR² de Nagelkerke (,035 e ,070) (tri) | D (difR) |
| Choi, S. W., Gibbons, L. E., & Crane, P. K. (2011). lordif: An R package for detecting differential item functioning using iterative hybrid ordinal logistic regression/item response theory and Monte Carlo simulations. *Journal of Statistical Software, 39*(8), 1–30. https://doi.org/10.18637/jss.v039.i08 | artigo | ΔR² de McFadden ≥ ,02 para itens ordinais (tri) | D |
| Şahin, A., & Anıl, D. (2017). The effects of test length and sample size on item parameters in item response theory. *Educational Sciences: Theory & Practice*. https://doi.org/10.12738/estp.2017.1.0270 | artigo | N mínimo por modelo citado no curso (tri); volume e páginas não vieram do Crossref nem do OpenAlex | Sec. (curso) |
| Linacre, J. M. (2002): ver "Citadas sem DOI". | | | |

## Redes e EGA

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Bringmann, L. F., Elmer, T., Epskamp, S., Krause, R. W., Schoch, D., Wichers, M., Wigman, J. T. W., & Snippe, E. (2019). What do centrality measures measure in psychological networks? *Journal of Abnormal Psychology, 128*(8), 892–903. https://doi.org/10.1037/abn0000446 | artigo | proximidade e intermediação não medem importância do nó (redes-ega) | R |
| Christensen, A. P., & Golino, H. (2021). Estimating the stability of psychological dimensions via bootstrap exploratory graph analysis: A Monte Carlo simulation and tutorial. *Psych, 3*(3), 479–500. https://doi.org/10.3390/psych3030032 | artigo | bootEGA, estabilidade de item (≥ .70), reamostragem com distribuição não normal (redes-ega) | R + página do artigo na MDPI (conferir no PDF antes de citação literal) |
| Constantin, M. A., Schuurman, N. K., & Vermunt, J. K. (2026). A general Monte Carlo method for sample size analysis in the context of network models. *Psychological Methods, 31*(3), 385–405. https://doi.org/10.1037/met0000555 | artigo | tamanho de amostra por simulação, pacote powerly (redes-ega) | R; publicado online antes do volume de 2026 |
| Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating psychological networks and their accuracy: A tutorial paper. *Behavior Research Methods, 50*(1), 195–212. https://doi.org/10.3758/s13428-017-0862-1 | artigo | coeficiente CS, bootstrap não paramétrico; referência do bootnet (redes-ega, relato-apa) | T |
| Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial correlation networks. *Psychological Methods, 23*(4), 617–634. https://doi.org/10.1037/met0000167 | artigo | EBICglasso (redes-ega) | R; preprint (arXiv 1607.01367) lido na pesquisa de 11/10/2026 |
| Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., Thiyagarajan, J. A., & Martinez-Molina, A. (2020). Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial. *Psychological Methods, 25*(3), 292–320. https://doi.org/10.1037/met0000255 | artigo | desempenho do EGA (redes-ega) | R |
| Golino, H. F., & Epskamp, S. (2017). Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research. *PLOS ONE, 12*(6), Article e0174035. https://doi.org/10.1371/journal.pone.0174035 | artigo | EGA (redes-ega, construcao-instrumentos) | R |
| Hallquist, M. N., Wright, A. G. C., & Molenaar, P. C. M. (2021). Problems with centrality measures in psychopathology symptom networks: Why network psychometrics cannot escape psychometric theory. *Multivariate Behavioral Research, 56*(2), 199–223. https://doi.org/10.1080/00273171.2019.1640103 | artigo | centralidade com variável latente; força pode repetir a carga fatorial (redes-ega) | R |
| Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models predict observations? On the importance of predictability in network models. *Behavior Research Methods, 50*(2), 853–861. https://doi.org/10.3758/s13428-017-0910-x | artigo | previsibilidade do nó (redes-ega) | R |
| Isvoranu, A.-M., & Epskamp, S. (2023). Which estimation method to choose in network psychometrics? Deriving guidelines for applied researchers. *Psychological Methods, 28*(4), 925–946. https://doi.org/10.1037/met0000439 | artigo | escolha do estimador de rede (redes-ega) | R |
| Jones, P. J., Ma, R., & McNally, R. J. (2021). Bridge centrality: A network approach to understanding comorbidity. *Multivariate Behavioral Research, 56*(2), 353–367. https://doi.org/10.1080/00273171.2019.1614898 | artigo | índices de ponte (redes-ega) | R |
| Kruis, J., & Maris, G. (2016). Three representations of the Ising model. *Scientific Reports, 6*, Article 34175. https://doi.org/10.1038/srep34175 | artigo | Ising, TRI multidimensional e colisor matematicamente equivalentes; o ajuste não escolhe entre rede e fator (redes-ega, fontes-e-decisoes) | T |
| Kursa, M. B., & Rudnicki, W. R. (2010). Feature selection with the Boruta package. *Journal of Statistical Software, 36*(11), 1–13. https://doi.org/10.18637/jss.v036.i11 | artigo (software) | Boruta (redes-ega) | Tít.; método descrito por Popa et al. (2025) |
| Popa, C. O., Schenk, A., Cojocaru, C. M., Gyorgy, M., Sava, F. A., Marian, Ș., Popoviciu, H., & Szasz, S. (2025). A transdiagnostic approach in psychodiagnosis: The Romanian adaptation of the Process-Based Assessment Tool (PBAT). *Behavioral Sciences, 15*(12), Article 1697. https://doi.org/10.3390/bs15121697 | artigo | adaptação romena do PBAT; redes em Pearson, NCT, Boruta (redes-ega, pbat-mpfi) | T (PubMed Central) |
| van Bork, R., Rhemtulla, M., Waldorp, L. J., Kruis, J., Rezvanifar, S., & Borsboom, D. (2021). Latent variable models and networks: Statistical equivalence and testability. *Multivariate Behavioral Research, 56*(2), 175–198. https://doi.org/10.1080/00273171.2019.1672515 | artigo | equivalência entre rede e variável latente; fator único contra rede esparsa; teste PCL, não implementado (redes-ega, fontes-e-decisoes) | T (seções de equivalência, implicações divergentes e PCL, na Taylor & Francis); publicado online em 2019 |
| van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023). Comparing network structures on three aspects: A permutation test. *Psychological Methods, 28*(6), 1273–1285. https://doi.org/10.1037/met0000476 | artigo | NCT (redes-ega) | R |
| Isvoranu, A.-M., Epskamp, S., Waldorp, L. J., & Borsboom, D. (2022). *Network psychometrics with R: A guide for behavioral and social scientists*. Routledge. https://doi.org/10.4324/9781003111238 | livro | livro-base da pasta Psicometria-R para redes (caps. 9 a 11 para redes longitudinais; redes-ega, pbat-mpfi) | T (pasta Psicometria-R) |

## Diferenças entre grupos

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Barcikowski, R. S., & Robey, R. R. (1984). Decisions in single group repeated measures analysis: Statistical tests and three computer packages. *The American Statistician, 38*(2), 148–150. https://doi.org/10.1080/00031305.1984.10483187 | artigo | corte ,75 do ε, formulado a partir de Huynh e Feldt (1976) segundo Blanca et al. (2023) (fontes-e-decisoes, diferencas) | Tít.; atribuição pela Introdução de Blanca et al. (2023) |
| Jacobson, N. S., & Truax, P. (1991). Clinical significance: A statistical approach to defining meaningful change in psychotherapy research. *Journal of Consulting and Clinical Psychology, 59*(1), 12–19. https://doi.org/10.1037/0022-006X.59.1.12 | artigo | índice de mudança confiável, `mudanca_confiavel()` (est_comum) | Tít. |
| Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. *Journal of Experimental Social Psychology, 49*(4), 764–766. https://doi.org/10.1016/j.jesp.2013.03.013 | artigo | atípicos por MAD, limite 2,5, `atipicos()` (est_comum) | Tít. |
| Blanca, M. J., Arnau, J., García-Castro, F. J., Alarcón, R., & Bono, R. (2023). Repeated measures ANOVA and adjusted F-tests when sphericity is violated: Which procedure is best? *Frontiers in Psychology, 14*, Article 1192453. https://doi.org/10.3389/fpsyg.2023.1192453 | artigo | limiar de ε em .60, Mauchly dispensável (diferencas) | T |
| Brown, M. B., & Forsythe, A. B. (1974a). Robust tests for the equality of variances. *Journal of the American Statistical Association, 69*(346), 364–367. https://doi.org/10.1080/01621459.1974.10482955 | artigo (epônimo) | Levene centrado na mediana, rotulado "Brown-Forsythe" no JASP (diferencas) | Tít. |
| Brown, M. B., & Forsythe, A. B. (1974b). The small sample behavior of some statistics which test the equality of several means. *Technometrics, 16*(1), 129–132. https://doi.org/10.1080/00401706.1974.10489158 | artigo (epônimo) | F* de Brown-Forsythe na ANOVA de um fator (diferencas) | Tít. |
| Delacre, M., Lakens, D., & Leys, C. (2017). Why psychologists should by default use Welch's t-test instead of Student's t-test. *International Review of Social Psychology, 30*(1), 92–101. https://doi.org/10.5334/irsp.82 | artigo | Welch como padrão no t (diferencas) | R; correção de 2022 abaixo |
| Delacre, M., Lakens, D., & Leys, C. (2022). Correction: Why psychologists should by default use Welch's t-test instead of Student's t-test. *International Review of Social Psychology, 35*(1), Article 21. https://doi.org/10.5334/irsp.661 | correção | correção do artigo de 2017 (diferencas) | Tít.; conteúdo da correção não lido |
| Delacre, M., Leys, C., Mora, Y. L., & Lakens, D. (2019). Taking parametric assumptions seriously: Arguments for the use of Welch's F-test instead of the classical F-test in one-way ANOVA. *International Review of Social Psychology, 32*(1), Article 13. https://doi.org/10.5334/irsp.198 | artigo | rota Welch na ANOVA de um fator; F* ao lado (diferencas) | R |
| Dunlap, W. P., Cortina, J. M., Vaslow, J. B., & Burke, M. J. (1996). Meta-analysis of experiments with matched groups or repeated measures designs. *Psychological Methods, 1*(2), 170–177. https://doi.org/10.1037/1082-989X.1.2.170 | artigo | d pareado com correlação (diferencas) | R |
| Games, P. A., & Howell, J. F. (1976). Pairwise multiple comparison procedures with unequal n's and/or variances: A Monte Carlo study. *Journal of Educational Statistics, 1*(2), 113–125. https://doi.org/10.3102/10769986001002113 | artigo | Games-Howell (diferencas) | R |
| Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of association: An equivalence testing approach. *British Journal of Mathematical and Statistical Psychology, 63*(3), 527–537. https://doi.org/10.1348/000711009X475853 | artigo | TOST da correlação pelo z de Fisher, via eq. 8 de Lakens (2017) (diferencas, associacao) | Sec. (Lakens, 2017); Tít. |
| Greenhouse, S. W., & Geisser, S. (1959). On methods in the analysis of profile data. *Psychometrika, 24*(2), 95–112. https://doi.org/10.1007/BF02289823 | artigo | ε de Greenhouse-Geisser (diferencas) | Tít. |
| Greenland, S., Senn, S. J., Rothman, K. J., Carlin, J. B., Poole, C., Goodman, S. N., & Altman, D. G. (2016). Statistical tests, P values, confidence intervals, and power: A guide to misinterpretations. *European Journal of Epidemiology, 31*(4), 337–350. https://doi.org/10.1007/s10654-016-0149-3 | artigo | não significativo não quer dizer igual, equívocos 4, 6 e 8; aviso de `est_comum.R` (diferencas, associacao, afc-ordinal, fontes-e-decisoes) | T |
| Grissom, R. J., & Kim, J. J. (2012). *Effect sizes for research: Univariate and multivariate applications* (2nd ed.). Routledge. https://doi.org/10.4324/9780203803233 | livro | régua de r_rb citada pelo curso (diferencas, sem ano no arquivo) | Tít.; edição usada pelo curso não identificada; os cortes ,125/,304/,465 não foram achados em fonte publicada (pesquisa de 11/10/2026) |
| Huynh, H., & Feldt, L. S. (1976). Estimation of the Box correction for degrees of freedom from sample data in randomized block and split-plot designs. *Journal of Educational Statistics, 1*(1), 69–82. https://doi.org/10.3102/10769986001001069 | artigo | ε de Huynh-Feldt (diferencas) | R; o resumo não fixa o corte .75 (a origem da regra é discutida em Blanca et al., 2023) |
| Kerby, D. S. (2014). The simple difference formula: An approach to teaching nonparametric correlation. *Comprehensive Psychology, 3*, Article 11.IT.3.1. https://doi.org/10.2466/11.IT.3.1 | artigo | r_rb pareada (diferencas) | Tít. |
| Lakens, D. (2013). Calculating and reporting effect sizes to facilitate cumulative science: A practical primer for t-tests and ANOVAs. *Frontiers in Psychology, 4*, Article 863. https://doi.org/10.3389/fpsyg.2013.00863 | artigo | d_av, d_rm, d_z (diferencas) | T |
| Lakens, D. (2017). Equivalence tests: A practical primer for *t* tests, correlations, and meta-analyses. *Social Psychological and Personality Science, 8*(4), 355–362. https://doi.org/10.1177/1948550617697177 | artigo | TOST: eq. 3, 4, 6, 7 e 8, IC 90%, maior dos dois p, quatro resultados (diferencas, associacao, fontes-e-decisoes) | T |
| Langsrud, Ø. (2003). ANOVA for unbalanced data: Use Type II instead of Type III sums of squares. *Statistics and Computing, 13*(2), 163–167. https://doi.org/10.1023/A:1023260610025 | artigo | SQ Tipo II como alternativa (diferencas) | R |
| McKeon, J. J. (1974). F approximations to the distribution of Hotelling's T²₀. *Biometrika, 61*(2), 381–383. https://doi.org/10.1093/biomet/61.2.381 | artigo (epônimo) | aproximação de McKeon do Hotelling-Lawley (diferencas) | Tít. |
| Miller, G. A., & Chapman, J. P. (2001). Misunderstanding analysis of covariance. *Journal of Abnormal Psychology, 110*(1), 40–48. https://doi.org/10.1037/0021-843X.110.1.40 | artigo | ANCOVA com grupos não aleatorizados (diferencas) | R |
| Olson, C. L. (1976). On choosing a test statistic in multivariate analysis of variance. *Psychological Bulletin, 83*(4), 579–586. https://doi.org/10.1037/0033-2909.83.4.579 | artigo | Pillai como estatística principal da MANOVA (diferencas) | Tít. |
| Villasenor Alva, J. A., & González Estrada, E. (2009). A generalization of Shapiro–Wilk's test for multivariate normality. *Communications in Statistics: Theory and Methods, 38*(11), 1870–1883. https://doi.org/10.1080/03610920802474465 | artigo (epônimo) | Shapiro-Wilk multivariado do JASP na MANOVA (diferencas) | Tít. |
| Zimmerman, D. W. (2004). A note on preliminary tests of equality of variances. *British Journal of Mathematical and Statistical Psychology, 57*(1), 173–181. https://doi.org/10.1348/000711004849222 | artigo | contra o Levene como porta para Welch (diferencas) | R |

Cohen (1988) está em "Poder e tamanho amostral"; Holm (1979), sem DOI, em
"Citadas sem DOI".

## Associação e regressão

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Beasley, T. M., & Schumacker, R. E. (1995). Multiple regression approach to analyzing contingency tables: Post hoc and planned comparison procedures. *The Journal of Experimental Education, 64*(1), 79–93. https://doi.org/10.1080/00220973.1995.9943797 | artigo | Bonferroni nos resíduos ajustados (associacao) | Tít. |
| Bishara, A. J., & Hittner, J. B. (2012). Testing the significance of a correlation with nonnormal data: Comparison of Pearson, Spearman, transformation, and resampling approaches. *Psychological Methods, 17*(3), 399–417. https://doi.org/10.1037/a0028087 | artigo | contra Spearman automático (associacao) | R |
| Campbell, I. (2007). Chi-squared and Fisher–Irwin tests of two-by-two tables with small sample recommendations. *Statistics in Medicine, 26*(19), 3661–3675. https://doi.org/10.1002/sim.2832 | artigo | contra Yates; qui-quadrado N − 1 (associacao) | R |
| Cochran, W. G. (1954). Some methods for strengthening the common χ² tests. *Biometrics, 10*(4), 417–451. https://doi.org/10.2307/3001616 | artigo | regra das esperadas (associacao) | Tít. (conteúdo conhecido, não lido) |
| Cook, R. D. (1977). Detection of influential observation in linear regression. *Technometrics, 19*(1), 15–18. https://doi.org/10.1080/00401706.1977.10489493 | artigo | distância de Cook (associacao) | Tít. |
| DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). Comparing the areas under two or more correlated receiver operating characteristic curves: A nonparametric approach. *Biometrics, 44*(3), 837–845. https://doi.org/10.2307/2531595 | artigo | IC da AUC (associacao) | Tít. |
| Durbin, J., & Watson, G. S. (1950). Testing for serial correlation in least squares regression: I. *Biometrika, 37*(3/4), 409–428. https://doi.org/10.2307/2332391 | artigo | teste de Durbin-Watson (associacao) | Tít. |
| Durbin, J., & Watson, G. S. (1951). Testing for serial correlation in least squares regression: II. *Biometrika, 38*(1/2), 159–178. https://doi.org/10.2307/2332325 | artigo | teste de Durbin-Watson (associacao) | Tít. |
| Fagerland, M. W., Lydersen, S., & Laake, P. (2013). The McNemar test for binary matched-pairs data: Mid-p and asymptotic are better than exact conditional. *BMC Medical Research Methodology, 13*, Article 91. https://doi.org/10.1186/1471-2288-13-91 | artigo | McNemar sem correção ou mid-p (associacao) | Tít. |
| Funder, D. C., & Ozer, D. J. (2019). Evaluating effect size in psychological research: Sense and nonsense. *Advances in Methods and Practices in Psychological Science, 2*(2), 156–168. https://doi.org/10.1177/2515245919847202 | artigo | régua de r (associacao) | R; corrigendum de 2020 registrado |
| Gignac, G. E., & Szodorai, E. T. (2016). Effect size guidelines for individual differences researchers. *Personality and Individual Differences, 102*, 74–78. https://doi.org/10.1016/j.paid.2016.06.069 | artigo | régua de r .10/.20/.30 (associacao) | R (via FORRT) |
| Haberman, S. J. (1973). The analysis of residuals in cross-classified tables. *Biometrics, 29*(1), 205–220. https://doi.org/10.2307/2529686 | artigo | resíduos ajustados (associacao) | Tít. |
| Haviland, M. G. (1990). Yates's correction for continuity and the analysis of 2 × 2 contingency tables. *Statistics in Medicine, 9*(4), 363–367. https://doi.org/10.1002/sim.4780090403 | artigo | contra Yates (associacao) | R |
| Hollander, M., Wolfe, D. A., & Chicken, E. (2015). *Nonparametric statistical methods* (3rd ed.). Wiley. https://doi.org/10.1002/9781119196037 | livro | IC do tau de Kendall (associacao, sem ano no arquivo) | Tít.; ano do Crossref para a versão com DOI; a impressa (ISBN 978-0-470-38737-5) consta no Open Library com data de 04/12/2013 |
| Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied logistic regression* (3rd ed.). Wiley. https://doi.org/10.1002/9781118548387 | livro | faixas de AUC .70/.80/.90 (associacao) | Sec. (CBU Statistics Wiki, que cita a 2ª ed.); página não conferida |
| Meng, X.-L., Rosenthal, R., & Rubin, D. B. (1992). Comparing correlated correlation coefficients. *Psychological Bulletin, 111*(1), 172–175. https://doi.org/10.1037/0033-2909.111.1.172 | artigo | correlações dependentes (associacao) | Tít. |
| O'Brien, R. M. (2007). A caution regarding rules of thumb for variance inflation factors. *Quality & Quantity, 41*(5), 673–690. https://doi.org/10.1007/s11135-006-9018-6 | artigo | VIF 10 não é critério de invalidade (associacao) | R |
| Schober, P., Boer, C., & Schwarte, L. A. (2018). Correlation coefficients: Appropriate use and interpretation. *Anesthesia & Analgesia, 126*(5), 1763–1768. https://doi.org/10.1213/ANE.0000000000002864 | artigo | régua de Schober (associacao) | Tít.; DOI impresso no livro do usuário, agora conferido |
| Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. *Psychological Bulletin, 87*(2), 245–251. https://doi.org/10.1037/0033-2909.87.2.245 | artigo | z de Steiger (associacao) | Tít. |
| Youden, W. J. (1950). Index for rating diagnostic tests. *Cancer, 3*(1), 32–35. https://doi.org/10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3 | artigo | índice de Youden (associacao) | Tít. |
| Zou, G. Y. (2007). Toward using confidence intervals to compare correlations. *Psychological Methods, 12*(4), 399–413. https://doi.org/10.1037/1082-989X.12.4.399 | artigo | IC da diferença entre correlações (associacao) | Tít. |

Cohen (1988) está em "Poder e tamanho amostral"; Tabachnick e Fidell (2019),
sem DOI, em "Citadas sem DOI".
| Bland, J. M., & Altman, D. G. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. *The Lancet, 327*(8476), 307–310. https://doi.org/10.1016/S0140-6736(86)90837-8 | artigo | limites de concordância em `confiabilidade_teste_reteste()` (construcao-instrumentos) | Tít. |
| Zou, G. (2004). A modified Poisson regression approach to prospective studies with binary data. *American Journal of Epidemiology, 159*(7), 702–706. https://doi.org/10.1093/aje/kwh090 | artigo | `razao_prevalencia()`: Poisson com EP sanduíche (epidemiologia) | Tít. |
| Barros, A. J. D., & Hirakata, V. N. (2003). Alternatives for logistic regression in cross-sectional studies: An empirical comparison of models that directly estimate the prevalence ratio. *BMC Medical Research Methodology, 3*, Article 21. https://doi.org/10.1186/1471-2288-3-21 | artigo | RP em vez de RC com desfecho comum (epidemiologia) | Tít. |

## Processos (mediação e moderação)

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and multilevel regression: Inferential and graphical techniques. *Multivariate Behavioral Research, 40*(3), 373–400. https://doi.org/10.1207/s15327906mbr4003_5 | artigo | Johnson-Neyman (processos) | R |
| Dalal, D. K., & Zickar, M. J. (2012). Some common myths about centering predictor variables in moderated multiple regression and polynomial regression. *Organizational Research Methods, 15*(3), 339–362. https://doi.org/10.1177/1094428111430540 | artigo | centralização não muda b3 (processos) | R |
| Echambadi, R., & Hess, J. D. (2007). Mean-centering does not alleviate collinearity problems in moderated multiple regression models. *Marketing Science, 26*(3), 438–445. https://doi.org/10.1287/mksc.1060.0263 | artigo | centralização e colinearidade (processos) | R |
| Fritz, M. S., Taylor, A. B., & MacKinnon, D. P. (2012). Explanation of two anomalous results in statistical mediation analysis. *Multivariate Behavioral Research, 47*(1), 61–87. https://doi.org/10.1080/00273171.2012.640596 | artigo | IC percentílico em vez de BC/BCa (processos) | R |
| Hayes, A. F., & Matthes, J. (2009). Computational procedures for probing interactions in OLS and logistic regression: SPSS and SAS implementations. *Behavior Research Methods, 41*(3), 924–936. https://doi.org/10.3758/BRM.41.3.924 | artigo | pick-a-point e Johnson-Neyman (processos) | R |
| Hayes, A. F., & Scharkow, M. (2013). The relative trustworthiness of inferential tests of the indirect effect in statistical mediation analysis: Does method really matter? *Psychological Science, 24*(10), 1918–1927. https://doi.org/10.1177/0956797613480187 | artigo | IC percentílico como compromisso (processos) | R |
| MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., & Sheets, V. (2002). A comparison of methods to test mediation and other intervening variable effects. *Psychological Methods, 7*(1), 83–104. https://doi.org/10.1037/1082-989X.7.1.83 | artigo | Sobel conservador (processos) | R; baixo poder do Sobel tratado no texto, não lido |
| MacKinnon, D. P., Warsi, G., & Dwyer, J. H. (1995). A simulation study of mediated effect measures. *Multivariate Behavioral Research, 30*(1), 41–62. https://doi.org/10.1207/s15327906mbr3001_3 | artigo | proporção mediada pede N grande (processos) | R |
| Maxwell, S. E., & Cole, D. A. (2007). Bias in cross-sectional analyses of longitudinal mediation. *Psychological Methods, 12*(1), 23–44. https://doi.org/10.1037/1082-989X.12.1.23 | artigo | limite causal da mediação transversal (processos) | R |
| Preacher, K. J., & Hayes, A. F. (2008). Asymptotic and resampling strategies for assessing and comparing indirect effects in multiple mediator models. *Behavior Research Methods, 40*(3), 879–891. https://doi.org/10.3758/BRM.40.3.879 | artigo | mediação múltipla e contrastes (processos) | R (PubMed) |
| Rucker, D. D., Preacher, K. J., Tormala, Z. L., & Petty, R. E. (2011). Mediation analysis in social psychology: Current practices and new recommendations. *Social and Personality Psychology Compass, 5*(6), 359–371. https://doi.org/10.1111/j.1751-9004.2011.00355.x | artigo | contra rótulos completa/parcial (processos) | R |
| Spiller, S. A., Fitzsimons, G. J., Lynch, J. G., Jr., & McClelland, G. H. (2013). Spotlights, floodlights, and the magic number zero: Simple effects tests in moderated regression. *Journal of Marketing Research, 50*(2), 277–288. https://doi.org/10.1509/jmr.12.0420 | artigo | valores focais e floodlight (processos) | R |

Hayes (2022), livro sem DOI, está em "Citadas sem DOI".

## Modelos avançados (GLM, GEE, multinível, meta-análise)

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random effects structure for confirmatory hypothesis testing: Keep it maximal. *Journal of Memory and Language, 68*(3), 255–278. https://doi.org/10.1016/j.jml.2012.11.001 | artigo | estrutura aleatória máxima (modelos) | R |
| Duval, S., & Tweedie, R. (2000). Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. *Biometrics, 56*(2), 455–463. https://doi.org/10.1111/j.0006-341X.2000.00455.x | artigo | trim-and-fill como sensibilidade (modelos) | Tít. |
| Egger, M., Davey Smith, G., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. *BMJ, 315*(7109), 629–634. https://doi.org/10.1136/bmj.315.7109.629 | artigo | teste de Egger (modelos) | R |
| Enders, C. K., & Tofighi, D. (2007). Centering predictor variables in cross-sectional multilevel models: A new look at an old issue. *Psychological Methods, 12*(2), 121–138. https://doi.org/10.1037/1082-989X.12.2.121 | artigo | centralização no grupo (modelos) | R |
| Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning variation in multilevel models. *Understanding Statistics, 1*(4), 223–231. https://doi.org/10.1207/S15328031US0104_02 | artigo | ICC latente com π²/3 (modelos) | Tít. |
| Harbord, R. M., Egger, M., & Sterne, J. A. C. (2006). A modified test for small-study effects in meta-analyses of controlled trials with binary endpoints. *Statistics in Medicine, 25*(20), 3443–3457. https://doi.org/10.1002/sim.2380 | artigo | teste de Harbord (modelos) | R |
| Higgins, J. P. T., Thompson, S. G., Deeks, J. J., & Altman, D. G. (2003). Measuring inconsistency in meta-analyses. *BMJ, 327*(7414), 557–560. https://doi.org/10.1136/bmj.327.7414.557 | artigo | I² (modelos) | Tít.; faixas 25/50/75 conhecidas, não conferidas no texto |
| Hin, L.-Y., & Wang, Y.-G. (2009). Working-correlation-structure identification in generalized estimating equations. *Statistics in Medicine, 28*(4), 642–658. https://doi.org/10.1002/sim.3489 | artigo | CIC (modelos) | R |
| IntHout, J., Ioannidis, J. P. A., & Borm, G. F. (2014). The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. *BMC Medical Research Methodology, 14*, Article 25. https://doi.org/10.1186/1471-2288-14-25 | artigo | Hartung-Knapp (modelos) | R |
| IntHout, J., Ioannidis, J. P. A., Rovers, M. M., & Goeman, J. J. (2016). Plea for routinely presenting prediction intervals in meta-analysis. *BMJ Open, 6*(7), Article e010247. https://doi.org/10.1136/bmjopen-2015-010247 | artigo | intervalo de predição (modelos) | R |
| Lai, M. H. C., & Kwok, O. (2015). Examining the rule of thumb of not using multilevel modeling: The "design effect smaller than two" rule. *The Journal of Experimental Education, 83*(3), 423–438. https://doi.org/10.1080/00220973.2014.907229 | artigo | limite da regra deff < 2 (modelos) | R |
| Langan, D., Higgins, J. P. T., Jackson, D., Bowden, J., Veroniki, A. A., Kontopantelis, E., Viechtbauer, W., & Simmonds, M. (2019). A comparison of heterogeneity variance estimators in simulated random-effects meta-analyses. *Research Synthesis Methods, 10*(1), 83–98. https://doi.org/10.1002/jrsm.1316 | artigo | τ² por REML (modelos) | R |
| Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. *Biometrika, 73*(1), 13–22. https://doi.org/10.1093/biomet/73.1.13 | artigo | GEE e EP sanduíche (modelos) | Tít. |
| Luke, S. G. (2017). Evaluating significance in linear mixed-effects models in R. *Behavior Research Methods, 49*(4), 1494–1502. https://doi.org/10.3758/s13428-016-0809-y | artigo | gl de Satterthwaite com REML (modelos) | R |
| Maas, C. J. M., & Hox, J. J. (2005). Sufficient sample sizes for multilevel modeling. *Methodology, 1*(3), 86–92. https://doi.org/10.1027/1614-2241.1.3.86 | artigo | EP de nível 2 com 50 grupos ou menos (modelos) | R; a variante de DOI "1614-1881" não existe |
| Mancl, L. A., & DeRouen, T. A. (2001). A covariance estimator for GEE with improved small-sample properties. *Biometrics, 57*(1), 126–134. https://doi.org/10.1111/j.0006-341X.2001.00126.x | artigo | correção de Mancl e DeRouen (modelos) | R |
| Matuschek, H., Kliegl, R., Vasishth, S., Baayen, H., & Bates, D. (2017). Balancing Type I error and power in linear mixed models. *Journal of Memory and Language, 94*, 305–315. https://doi.org/10.1016/j.jml.2017.01.001 | artigo | estrutura aleatória parcimoniosa (modelos) | R |
| Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The coefficient of determination R² and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. *Journal of the Royal Society Interface, 14*(134), Article 20170213. https://doi.org/10.1098/rsif.2017.0213 | artigo | R² com inclinações aleatórias (modelos) | Tít. |
| Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for obtaining R² from generalized linear mixed-effects models. *Methods in Ecology and Evolution, 4*(2), 133–142. https://doi.org/10.1111/j.2041-210x.2012.00261.x | artigo | R² marginal e condicional (modelos) | R |
| Pan, W. (2001). Akaike's information criterion in generalized estimating equations. *Biometrics, 57*(1), 120–125. https://doi.org/10.1111/j.0006-341X.2001.00120.x | artigo | QIC e QICu (modelos) | R |
| Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., & Rushton, L. (2006). Comparison of two methods to detect publication bias in meta-analysis. *JAMA, 295*(6), 676–680. https://doi.org/10.1001/jama.295.6.676 | artigo | teste tipo Peters para OR e RR (modelos) | R; o Crossref lista só a primeira autora |
| Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., & Rushton, L. (2007). Performance of the trim and fill method in the presence of publication bias and between-study heterogeneity. *Statistics in Medicine, 26*(25), 4544–4562. https://doi.org/10.1002/sim.2889 | artigo | trim-and-fill só como sensibilidade (modelos) | R |
| Pinheiro, J. C., & Bates, D. M. (2000). *Mixed-effects models in S and S-PLUS*. Springer. https://doi.org/10.1007/b98882 | livro | ML para partes fixas, REML no modelo final (modelos) | Tít.; ISBN 0-387-98957-9 no Crossref; o registro do Crossref não traz autores |
| Riley, R. D., Higgins, J. P. T., & Deeks, J. J. (2011). Interpretation of random effects meta-analyses. *BMJ, 342*, Article d549. https://doi.org/10.1136/bmj.d549 | artigo | intervalo de predição (modelos) | Tít. |
| Robins, J. M., Rotnitzky, A., & Zhao, L. P. (1995). Analysis of semiparametric regression models for repeated outcomes in the presence of missing data. *Journal of the American Statistical Association, 90*(429), 106–121. https://doi.org/10.1080/01621459.1995.10476493 | artigo | GEE ponderado sob MAR (modelos) | Tít. |
| Schafer, J. L., & Graham, J. W. (2002). Missing data: Our view of the state of the art. *Psychological Methods, 7*(2), 147–177. https://doi.org/10.1037/1082-989X.7.2.147 | artigo | imputação múltipla e ML (modelos) | R |
| Self, S. G., & Liang, K.-Y. (1987). Asymptotic properties of maximum likelihood estimators and likelihood ratio tests under nonstandard conditions. *Journal of the American Statistical Association, 82*(398), 605–610. https://doi.org/10.1080/01621459.1987.10478472 | artigo | LRT na fronteira (modelos) | Tít. |
| Sterne, J. A. C., Sutton, A. J., Ioannidis, J. P. A., Terrin, N., Jones, D. R., Lau, J., Carpenter, J., Rücker, G., Harbord, R. M., Schmid, C. H., Tetzlaff, J., Deeks, J. J., Peters, J., Macaskill, P., Schwarzer, G., Duval, S., Altman, D. G., Moher, D., & Higgins, J. P. T. (2011). Recommendations for examining and interpreting funnel plot asymmetry in meta-analyses of randomised controlled trials. *BMJ, 343*, Article d4002. https://doi.org/10.1136/bmj.d4002 | artigo | assimetria do funil só com k ≥ 10 (modelos) | Sec. (Cochrane Handbook 6.5, seção 13.3.4.4) |
| Stram, D. O., & Lee, J. W. (1994). Variance components testing in the longitudinal mixed effects model. *Biometrics, 50*(4), 1171–1177. https://doi.org/10.2307/2533455 | artigo | mistura de χ² na fronteira (modelos) | R; correção de 1995 publicada, não ligada no `updated-by` |
| Terrin, N., Schmid, C. H., Lau, J., & Olkin, I. (2003). Adjusting for publication bias in the presence of heterogeneity. *Statistics in Medicine, 22*(13), 2113–2126. https://doi.org/10.1002/sim.1461 | artigo | trim-and-fill falha com heterogeneidade (modelos) | R |
| Veroniki, A. A., Jackson, D., Viechtbauer, W., Bender, R., Bowden, J., Knapp, G., Kuss, O., Higgins, J. P. T., Langan, D., & Salanti, G. (2016). Methods to estimate the between-study variance and its uncertainty in meta-analysis. *Research Synthesis Methods, 7*(1), 55–79. https://doi.org/10.1002/jrsm.1164 | artigo | estimadores de τ² (modelos) | R |

Schafer e Graham (2002) também apoiam a nota de ausentes de `afc-ordinal.md`
(que não os cita). As faixas sobrepostas de I² da Cochrane estão em "Citadas sem
DOI".
| Schwarzer, G., Chemaitelly, H., Abu-Raddad, L. J., & Rücker, G. (2019). Seriously misleading results using inverse of Freeman-Tukey double arcsine transformation in meta-analysis of single proportions. *Research Synthesis Methods, 10*(3), 476–483. https://doi.org/10.1002/jrsm.1348 | artigo | logito como padrão em `meta_proporcao()` (modelos) | Tít.; D (meta) |
| Reitsma, J. B., Glas, A. S., Rutjes, A. W. S., Scholten, R. J. P. M., Bossuyt, P. M., & Zwinderman, A. H. (2005). Bivariate analysis of sensitivity and specificity produces informative summary measures in diagnostic reviews. *Journal of Clinical Epidemiology, 58*(10), 982–990. https://doi.org/10.1016/j.jclinepi.2005.02.022 | artigo | `meta_acuracia()`: modelo bivariado (modelos) | D (mada) |
| Deeks, J. J., Macaskill, P., & Irwig, L. (2005). The performance of tests of publication bias and other sample size effects in systematic reviews of diagnostic test accuracy was assessed. *Journal of Clinical Epidemiology, 58*(9), 882–893. https://doi.org/10.1016/j.jclinepi.2005.01.016 | artigo | teste de assimetria de Deeks (modelos) | Tít. |
| Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number of components in a normal mixture. *Biometrika, 88*(3), 767–778. https://doi.org/10.1093/biomet/88.3.767 | artigo | LMR em `classes_latentes()` e `perfis_latentes()` (classes-perfis); Crossref e OpenAlex listam só o primeiro autor | D (tidyLPA) |
| Cranford, J. A., Shrout, P. E., Iida, M., Rafaeli, E., Yip, T., & Bolger, N. (2006). A procedure for evaluating sensitivity to within-person change: Can mood measures in diary studies detect change reliably? *Personality and Social Psychology Bulletin, 32*(7), 917–929. https://doi.org/10.1177/0146167206287721 | artigo | Rc de `confiabilidade_ema()` (construcao-instrumentos, pbat-mpfi) | Tít.; D (psych) |

## Poder e tamanho amostral

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Chen, H., Cohen, P., & Chen, S. (2010). How big is a big odds ratio? Interpreting the magnitudes of odds ratios in epidemiological studies. *Communications in Statistics: Simulation and Computation, 39*(4), 860–864. https://doi.org/10.1080/03610911003650383 | artigo | cortes de OR (poder-amostra) | Tít. |
| Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). Lawrence Erlbaum Associates. https://doi.org/10.4324/9780203771587 | livro | convenções de magnitude; f, w, f² (poder-amostra, diferencas, associacao) | Tít.; o DOI é da reedição Routledge (Crossref: 2013); edição de 1988 conferida pelo ISBN 0-8058-0283-5 no Open Library |
| Cohen, J. (1992). A power primer. *Psychological Bulletin, 112*(1), 155–159. https://doi.org/10.1037/0033-2909.112.1.155 | artigo | convenções de magnitude (poder-amostra) | R (PubMed) |
| Demidenko, E. (2007). Sample size determination for logistic regression revisited. *Statistics in Medicine, 26*(18), 3385–3397. https://doi.org/10.1002/sim.2771 | artigo | método da logística no G*Power e no pwrss (poder-amostra) | Tít. |
| Faul, F., Erdfelder, E., Lang, A.-G., & Buchner, A. (2007). G*Power 3: A flexible statistical power analysis program for the social, behavioral, and biomedical sciences. *Behavior Research Methods, 39*(2), 175–191. https://doi.org/10.3758/BF03193146 | artigo (software) | G*Power 3 e exemplo de medidas repetidas (poder-amostra) | T no exemplo de medidas repetidas, reproduzido nos evals (segundo poder-amostra.md) |
| Faul, F., Erdfelder, E., Buchner, A., & Lang, A.-G. (2009). Statistical power analyses using G*Power 3.1: Tests for correlation and regression analyses. *Behavior Research Methods, 41*(4), 1149–1160. https://doi.org/10.3758/BRM.41.4.1149 | artigo (software) | R² other X na logística, caminho de Demidenko sem correção (poder-amostra) | T nos dois exemplos de análise de compromisso (z crítico 3,454879 e 2,146971 com R² other X = .40), reproduzidos em 14/10/2026 |
| Hsieh, F. Y., Bloch, D. A., & Larsen, M. D. (1998). A simple method of sample size calculation for linear and logistic regression. *Statistics in Medicine, 17*(14), 1623–1634. https://doi.org/10.1002/(SICI)1097-0258(19980730)17:14<1623::AID-SIM871>3.0.CO;2-S | artigo | método da logística no G*Power e no pwrss (poder-amostra) | Tít. |
| Lakens, D. (2022). Sample size justification. *Collabra: Psychology, 8*(1), Article 33267. https://doi.org/10.1525/collabra.33267 | artigo | justificativa de N e sensibilidade (poder-amostra) | R |
| MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and determination of sample size for covariance structure modeling. *Psychological Methods, 1*(2), 130–149. https://doi.org/10.1037/1082-989X.1.2.130 | artigo | `poder_afc_rmsea()` (poder-amostra) | Sec. (Kelley & Lai, 2018, T); D (semTools) |
| MacCallum, R. C., Browne, M. W., & Cai, L. (2006). Testing differences between nested covariance structure models: Power analysis and null hypotheses. *Psychological Methods, 11*(1), 19–35. https://doi.org/10.1037/1082-989X.11.1.19 | artigo | `poder_afc_aninhados()` por RMSEA (poder-amostra) | Sec. (Kelley & Lai, 2018, T) |
| Satorra, A., & Saris, W. E. (1985). Power of the likelihood ratio test in covariance structure analysis. *Psychometrika, 50*(1), 83–90. https://doi.org/10.1007/BF02294150 | artigo | `poder_afc_aninhados()` com Σ (poder-amostra) | Sec. (Kelley & Lai, 2018, T) |
| Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the root mean square error of approximation: Sample size planning for narrow confidence intervals. *Multivariate Behavioral Research, 46*(1), 1–32. https://doi.org/10.1080/00273171.2011.543027 | artigo | `n_aipe_rmsea()` (poder-amostra) | Sec. (Kelley & Lai, 2018, T); D (MBESS) |
| Kelley, K., & Lai, K. (2018). Sample size planning for confirmatory factor models. In P. Irwing, T. Booth, & D. J. Hughes (Eds.), *The Wiley handbook of psychometric testing* (pp. 113–138). Wiley. https://doi.org/10.1002/9781118489772.ch5 | capítulo | as quatro abordagens e os exemplos reproduzidos nos evals (poder-amostra) | T (pasta Psicometria-R) |

## Construção de instrumentos

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings. *Educational and Psychological Measurement, 45*(1), 131–142. https://doi.org/10.1177/0013164485451012 | artigo | V de Aiken (construcao-instrumentos) | Tít. |
| Ayre, C., & Scally, A. J. (2014). Critical values for Lawshe's content validity ratio: Revisiting the original methods of calculation. *Measurement and Evaluation in Counseling and Development, 47*(1), 79–86. https://doi.org/10.1177/0748175613513808 | artigo | RVC crítico pelo método exato (construcao-instrumentos) | Tít.; tabela de 5 a 10 juízes reproduzida nos evals pela aula |
| Bartlett, M. S. (1950). Tests of significance in factor analysis. *British Journal of Statistical Psychology, 3*(2), 77–85. https://doi.org/10.1111/j.2044-8317.1950.tb00285.x | artigo | teste de esfericidade de Bartlett (construcao-instrumentos) | Tít. |
| Buja, A., & Eyuboglu, N. (1992). Remarks on parallel analysis. *Multivariate Behavioral Research, 27*(4), 509–540. https://doi.org/10.1207/s15327906mbr2704_2 | artigo | análise paralela por permutação (construcao-instrumentos) | Tít. |
| Glorfeld, L. W. (1995). An improvement on Horn's parallel analysis methodology for selecting the correct number of factors to retain. *Educational and Psychological Measurement, 55*(3), 377–393. https://doi.org/10.1177/0013164495055003002 | artigo | percentil 95 na análise paralela (construcao-instrumentos) | Tít. |
| Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. *Psychometrika, 30*(2), 179–185. https://doi.org/10.1007/BF02289447 | artigo | análise paralela (construcao-instrumentos) | Tít. |
| International Test Commission. (2018). ITC guidelines for translating and adapting tests (second edition). *International Journal of Testing, 18*(2), 101–134. https://doi.org/10.1080/15305058.2017.1398166 | artigo (diretriz) | adaptação transcultural (construcao-instrumentos, pbat-mpfi) | Tít. |
| Kaiser, H. F. (1974). An index of factorial simplicity. *Psychometrika, 39*(1), 31–36. https://doi.org/10.1007/BF02291575 | artigo | rótulos do KMO (construcao-instrumentos) | Tít. |
| Lawshe, C. H. (1975). A quantitative approach to content validity. *Personnel Psychology, 28*(4), 563–575. https://doi.org/10.1111/j.1744-6570.1975.tb01393.x | artigo | RVC original (construcao-instrumentos) | Tít. |
| Lynn, M. R. (1986). Determination and quantification of content validity. *Nursing Research, 35*(6), 382–386. https://doi.org/10.1097/00006199-198611000-00017 | artigo | IVC (construcao-instrumentos) | Tít. |
| Penfield, R. D., & Giacobbi, P. R., Jr. (2004). Applying a score confidence interval to Aiken's item content-relevance index. *Measurement in Physical Education and Exercise Science, 8*(4), 213–225. https://doi.org/10.1207/s15327841mpee0804_3 | artigo | IC escore do V de Aiken (construcao-instrumentos) | Tít. |
| Polit, D. F., & Beck, C. T. (2006). The content validity index: Are you sure you know what's being reported? Critique and recommendations. *Research in Nursing & Health, 29*(5), 489–497. https://doi.org/10.1002/nur.20147 | artigo | cortes do I-IVC e do S-IVC (construcao-instrumentos) | Tít. |
| Polit, D. F., Beck, C. T., & Owen, S. V. (2007). Is the CVI an acceptable indicator of content validity? Appraisal and recommendations. *Research in Nursing & Health, 30*(4), 459–467. https://doi.org/10.1002/nur.20199 | artigo | kappa modificado (construcao-instrumentos) | Tít. |

Hernández-Nieto (2002) e Cassepp-Borges et al. (2010) estão nas seções finais.
| Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The Hull method for selecting the number of common factors. *Multivariate Behavioral Research, 46*(2), 340–364. https://doi.org/10.1080/00273171.2011.564527 | artigo | Hull-CFI e Hull-CAF em `afe_numero_fatores()` (construcao-instrumentos) | Sec. (Timmerman et al., 2018, T); D (EFAtools) |
| Ceulemans, E., & Kiers, H. A. L. (2006). Selecting among three-mode principal component models of different types and complexities: A numerical convex hull based method. *British Journal of Mathematical and Statistical Psychology, 59*(1), 133–150. https://doi.org/10.1348/000711005X64817 | artigo | heurística do contorno convexo do Hull (construcao-instrumentos) | Sec. |
| Timmerman, M. E., Lorenzo-Seva, U., & Ceulemans, E. (2018). The number of factors problem. In P. Irwing, T. Booth, & D. J. Hughes (Eds.), *The Wiley handbook of psychometric testing* (pp. 305–324). Wiley. https://doi.org/10.1002/9781118489772.ch11 | capítulo | recomendação de Hull e paralela; MAP desaconselhado para fatores comuns (construcao-instrumentos) | T (pasta Psicometria-R) |
| Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations. *Psychometrika, 41*(3), 321–327. https://doi.org/10.1007/BF02293557 | artigo | MAP em `afe_numero_fatores()` | Sec.; D (psych) |
| Velicer, W. F., Eaton, C. A., & Fava, J. L. (2000). Construct explication through factor or component analysis: A review and evaluation of alternative procedures for determining the number of factors or components. In R. D. Goffin & E. Helmes (Eds.), *Problems and solutions in human assessment* (pp. 41–71). Springer. https://doi.org/10.1007/978-1-4615-4397-8_3 | capítulo | MAP revisado (quarta potência) | Sec. |
| McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some intraclass correlation coefficients. *Psychological Methods, 1*(1), 30–46. https://doi.org/10.1037/1082-989X.1.1.30 | artigo | ICC(A,1), ICC(C,1) e médias (construcao-instrumentos); correção registrada no `updated-by` | Sec.; D (psych, irr) |
| Koo, T. K., & Li, M. Y. (2016). A guideline of selecting and reporting intraclass correlation coefficients for reliability research. *Journal of Chiropractic Medicine, 15*(2), 155–163. https://doi.org/10.1016/j.jcm.2016.02.012 | artigo | régua do ICC aplicada ao IC e n ≥ 30 (construcao-instrumentos); errata registrada no `updated-by` | R |
| Weir, J. P. (2005). Quantifying test-retest reliability using the intraclass correlation coefficient and the SEM. *Journal of Strength and Conditioning Research, 19*(1), 231–240. https://doi.org/10.1519/15184.1 | artigo | EPM e mudança mínima detectável (construcao-instrumentos); Crossref e OpenAlex dão só a primeira página, a última (240) é a da citação usual, não conferida | Tít. |
| Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in assessing rater reliability. *Psychological Bulletin, 86*(2), 420–428. https://doi.org/10.1037/0033-2909.86.2.420 | artigo | exemplo de 6 alvos e 4 juízes reproduzido nos evals | T (tabela de valores) |
| Cohen, J. (1960). A coefficient of agreement for nominal scales. *Educational and Psychological Measurement, 20*(1), 37–46. https://doi.org/10.1177/001316446002000104 | artigo | kappa (construcao-instrumentos) | Tít. |
| Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision for scaled disagreement or partial credit. *Psychological Bulletin, 70*(4), 213–220. https://doi.org/10.1037/h0026256 | artigo | kappa ponderado | Tít. |
| Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. *Psychological Bulletin, 76*(5), 378–382. https://doi.org/10.1037/h0031619 | artigo | kappa de Fleiss | Tít. |
| Landis, J. R., & Koch, G. G. (1977). The measurement of observer agreement for categorical data. *Biometrics, 33*(1), 159–174. https://doi.org/10.2307/2529310 | artigo | régua do kappa; registro do JSTOR com só a primeira página no Crossref e no OpenAlex, a última (174) é a da citação usual, não conferida | Sec. |
| Ebel, R. L. (1954). Procedures for the analysis of classroom tests. *Educational and Psychological Measurement, 14*(2), 352–364. https://doi.org/10.1177/001316445401400215 | artigo | ULI ≥ ,20 em `analise_classica_itens()` | Sec. (Martinková & Hladká, 2023, T) |
| Martinková, P., & Hladká, A. (2023). *Computational aspects of psychometric methods: With R*. Chapman and Hall/CRC. https://doi.org/10.1201/9781003054313 | livro | análise clássica de itens (dificuldade, RIT, RIR, ULI); fonte do CVR de Lawshe (construcao-instrumentos) | T (pasta Psicometria-R) |
| Petersen, I. T. (2024). *Principles of psychological assessment: With applied examples in R*. Chapman and Hall/CRC. https://doi.org/10.1201/9781003357421 | livro | viés preditivo (cap. 15), ROC por utilidade (cap. 9), teoria da generalizabilidade (cap. 5) | T (pasta Psicometria-R) |
| Irwing, P., Booth, T., & Hughes, D. J. (Eds.). (2018). *The Wiley handbook of psychometric testing: A multidisciplinary reference on survey, scale and test development* (Vols. 1–2). Wiley. https://doi.org/10.1002/9781118489772 | livro | obra de referência da pasta (capítulos 5, 6, 11, 22) | T (pasta Psicometria-R) |
| Rust, J., Kosinski, M., & Stillwell, D. (2021). *Modern psychometrics: The science of psychological assessment* (4th ed.). Routledge. https://doi.org/10.4324/9781315637686 | livro | construção de questionários (fontes-e-decisoes); o Crossref registra 2020 (publicação online) | T (pasta Psicometria-R) |
| Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment of ordered polytomous items with parallel analysis. *Psychological Methods, 16*(2), 209–220. https://doi.org/10.1037/a0023353 | artigo | PA-MRFA, preferida pelo cap. 11 do Wiley Handbook e não implementada (construcao-instrumentos) | Sec. (Timmerman et al., 2018, T) |

## Relato

`relato-apa.md` não cita trabalho próprio. As obras que ele cita estão nas
seções acima: Rhemtulla et al. (2012), Hu e Bentler (1999), Xia e Yang (2019),
Green e Yang (2009) e Rosseel (2012) em "Psicometria e AFC"; Chalmers (2012) em
"TRI"; Epskamp, Borsboom e Fried (2018) em "Redes e EGA". Referência de pacote
sai de `citation()` no R da sessão.

## v1.5: módulos novos (não paramétricos complexos, MEE, meta-análise em rede, extração, três passos, PGM, tendência, ausentes, escolha forçada)

Conferidas em 29/10/2026 com `verify_citations.py` (todas VERIFIED, nenhuma retratação); autores, ano, periódico, volume e páginas do registro do Crossref (`api.crossref.org/works`).

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Brunner, E., Konietschke, F., Pauly, M., & Puri, M. L. (2017). Rank-based procedures in factorial designs: Hypotheses about non-parametric treatment effects. *Journal of the Royal Statistical Society Series B: Statistical Methodology, 79*(5), 1463–1485. https://doi.org/10.1111/rssb.12222 | artigo | ATS e WTS, efeito relativo (`np_fatorial()`) | Tít. |
| Noguchi, K., Gel, Y. R., Brunner, E., & Konietschke, F. (2012). nparLD: An R software package for the nonparametric analysis of longitudinal data in factorial experiments. *Journal of Statistical Software, 50*(12). https://doi.org/10.18637/jss.v050.i12 | artigo (software) | delineamentos LD-F1, F1-LD-F1, F2-LD-F1 e RTE (`np_longitudinal()`) | Tít. |
| Konietschke, F., Hothorn, L. A., & Brunner, E. (2012). Rank-based multiple test procedures and simultaneous confidence intervals. *Electronic Journal of Statistics, 6*, 738–759. https://doi.org/10.1214/12-EJS691 | artigo | comparações múltiplas do mctp (não paramétricos) | Tít. |
| Konietschke, F., Placzek, M., Schaarschmidt, F., & Hothorn, L. A. (2015). nparcomp: An R software package for nonparametric multiple comparisons and simultaneous confidence intervals. *Journal of Statistical Software, 64*(9). https://doi.org/10.18637/jss.v064.i09 | artigo (software) | mctp e mctp.rm (não paramétricos) | Tít. |
| Friedrich, S., & Pauly, M. (2018). MATS: Inference for potentially singular and heteroscedastic MANOVA. *Journal of Multivariate Analysis, 165*, 166–179. https://doi.org/10.1016/j.jmva.2017.12.008 | artigo | MATS com bootstrap paramétrico (`manova_np()`) | Tít. |
| Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of orthogonalizing powered and product terms: Implications for modeling interactions among latent variables. *Structural Equation Modeling, 13*(4), 497–519. https://doi.org/10.1207/s15328007sem1304_1 | artigo | centralização residual dos indicadores-produto (`interacao_latente()`) | Tít. |
| Marsh, H. W., Wen, Z., & Hau, K.-T. (2004). Structural equation models of latent interactions: Evaluation of alternative estimation strategies and indicator construction. *Psychological Methods, 9*(3), 275–300. https://doi.org/10.1037/1082-989X.9.3.275 | artigo | indicadores-produto pareados (`interacao_latente()`) | Tít. |
| Schoemann, A. M., & Jorgensen, T. D. (2021). Testing and interpreting latent variable interactions using the semTools package. *Psych, 3*(3), 322–335. https://doi.org/10.3390/psych3030024 | artigo | indProd e probe2WayRC, citado no código das aulas (mee) | Tít. |
| Epskamp, S. (2015). semPlot: Unified visualizations of structural equation models. *Structural Equation Modeling, 22*(3), 474–483. https://doi.org/10.1080/10705511.2014.937847 | artigo (software) | diagrama do modelo (figuras-dados) | Tít. |
| Rücker, G. (2012). Network meta-analysis, electrical networks and graph theory. *Research Synthesis Methods, 3*(4), 312–324. https://doi.org/10.1002/jrsm.1058 | artigo | modelo do netmeta (`meta_rede()`) | Tít. |
| Rücker, G., & Schwarzer, G. (2015). Ranking treatments in frequentist network meta-analysis works without resampling methods. *BMC Medical Research Methodology, 15*, 58. https://doi.org/10.1186/s12874-015-0060-8 | artigo | P-score (`meta_rede()`) | Tít. |
| Wan, X., Wang, W., Liu, J., & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. *BMC Medical Research Methodology, 14*, 135. https://doi.org/10.1186/1471-2288-14-135 | artigo | DP pela amplitude e pelo IIQ (`extrair_media_dp()`) | Tít. |
| Luo, D., Wan, X., Liu, J., & Tong, T. (2018). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. *Statistical Methods in Medical Research, 27*(6), 1785–1805. https://doi.org/10.1177/0962280216669183 | artigo | média pela mediana (`extrair_media_dp()`) | Tít. |
| Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure models with categorical variables: One-step versus three-step estimators. *Political Analysis, 12*(1), 3–27. https://doi.org/10.1093/pan/mph001 | artigo | correção BCH (`classes_distal()`) | Tít. |
| Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved three-step approaches. *Political Analysis, 18*(4), 450–469. https://doi.org/10.1093/pan/mpq025 | artigo | pesos BCH e matriz D (`classes_distal()`) | Tít. |
| Scutari, M. (2010). Learning Bayesian networks with the bnlearn R package. *Journal of Statistical Software, 35*(3). https://doi.org/10.18637/jss.v035.i03 | artigo (software) | aprendizado da estrutura (`rede_bayesiana()`) | Tít. |
| Scutari, M., & Nagarajan, R. (2013). Identifying significant edges in graphical models of molecular networks. *Artificial Intelligence in Medicine, 57*(3), 207–217. https://doi.org/10.1016/j.artmed.2012.12.006 | artigo | limiar da rede média (`rede_bayesiana()`) | Tít. |
| Kim, H.-J., Fay, M. P., Feuer, E. J., & Midthune, D. N. (2000). Permutation tests for joinpoint regression with applications to cancer rates. *Statistics in Medicine, 19*(3), 335–351. https://doi.org/10.1002/(SICI)1097-0258(20000215)19:3<335::AID-SIM336>3.0.CO;2-Z | artigo | modelo de pontos de inflexão (`tendencia_joinpoint()`) | Tít. |
| Fay, M. P., & Feuer, E. J. (1997). Confidence intervals for directly standardized rates: A method based on the gamma distribution. *Statistics in Medicine, 16*(7), 791–801. https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7<791::AID-SIM500>3.0.CO;2-# | artigo | IC da taxa padronizada (`taxa_padronizada()`) | Tít. |
| Little, R. J. A. (1988). A test of missing completely at random for multivariate data with missing values. *Journal of the American Statistical Association, 83*(404), 1198–1202. https://doi.org/10.1080/01621459.1988.10478722 | artigo | teste MCAR (`ausentes()`) | Tít. |
| Brown, A., & Maydeu-Olivares, A. (2012). Fitting a Thurstonian IRT model to forced-choice data using Mplus. *Behavior Research Methods, 44*(4), 1135–1147. https://doi.org/10.3758/s13428-012-0217-x | artigo | correção dos graus de liberdade em blocos de 3 ou mais itens (`escolha_forcada()`) | Tít. |
| Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. *Psychometrika, 16*(3), 297–334. https://doi.org/10.1007/BF02310555 | artigo | média das divisões ao meio = alfa (`confiabilidade_metades()`) | Tít. |
| Kreidler, S. M., Muller, K. E., Grunwald, G. K., Ringham, B. M., Coker-Dukowitz, Z., Sakhadeo, U. R., Baron, A. E., & Glueck, D. H. (2013). GLIMMPSE: Online power computation for linear models with and without a baseline covariate. *Journal of Statistical Software, 54*(10). https://doi.org/10.18637/jss.v054.i10 | artigo (software) | poder para GLM e GEE (poder-amostra) | Tít. |
| Allen, M., Poggiali, D., Whitaker, K., Marshall, T. R., & Kievit, R. A. (2019). Raincloud plots: A multi-platform tool for robust data visualization. *Wellcome Open Research, 4*, 63. https://doi.org/10.12688/wellcomeopenres.15191.1 | artigo | gráfico raincloud (figuras-dados) | Tít. |

## v1.6: desfecho ordinal e nominal, concordância, confiabilidade, teoria G, dMACS, MH, validação interna e risco

Conferidas em 29/10/2026 com `verify_citations.py` (todas VERIFIED, nenhuma retratação nem `updated-by`); autores, ano, periódico, volume e páginas do registro do Crossref.

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Gwet, K. L. (2008). Computing inter-rater reliability and its variance in the presence of high agreement. *British Journal of Mathematical and Statistical Psychology, 61*(1), 29–48. https://doi.org/10.1348/000711006X126600 | artigo | AC1 e AC2 em `concordancia_juizes()` | construcao-instrumentos |
| Hayes, A. F., & Krippendorff, K. (2007). Answering the call for a standard reliability measure for coding data. *Communication Methods and Measures, 1*(1), 77–89. https://doi.org/10.1080/19312450709336664 | artigo | alfa de Krippendorff (nominal, ordinal, intervalar) | construcao-instrumentos |
| Feinstein, A. R., & Cicchetti, D. V. (1990). High agreement but low kappa: I. The problems of two paradoxes. *Journal of Clinical Epidemiology, 43*(6), 543–549. https://doi.org/10.1016/0895-4356(90)90158-L | artigo | aviso do paradoxo do kappa | construcao-instrumentos |
| Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for population reliability coefficients: Evaluation of methods, recommendations, and software for composite measures. *Psychological Methods, 21*(1), 69–92. https://doi.org/10.1037/a0040086 | artigo | IC por bootstrap em `confiabilidade(ic = TRUE)` | afc-ordinal |
| Spearman, C. (1910). Correlation calculated from faulty data. *British Journal of Psychology, 3*(3), 271–295. https://doi.org/10.1111/j.2044-8295.1910.tb00206.x | artigo | `spearman_brown()` | construcao-instrumentos |
| Brown, W. (1910). Some experimental results in the correlation of mental abilities. *British Journal of Psychology, 3*(3), 296–322. https://doi.org/10.1111/j.2044-8295.1910.tb00207.x | artigo | `spearman_brown()` | construcao-instrumentos |
| Brennan, R. L. (2001). *Generalizability theory*. Springer. https://doi.org/10.1007/978-1-4757-3456-0 | livro | `generalizabilidade()` (erro relativo e absoluto no p × i × o) | construcao-instrumentos |
| Nye, C. D., & Drasgow, F. (2011). Effect size indices for analyses of measurement equivalence: Understanding the practical importance of differences between groups. *Journal of Applied Psychology, 96*(5), 966–980. https://doi.org/10.1037/a0022955 | artigo | `dmacs()` e ΔMean | invariancia |
| Nye, C. D., Bradburn, J., Olenick, J., Bialko, C., & Drasgow, F. (2019). How big are my effects? Examining the magnitude of effect sizes in studies of measurement equivalence. *Organizational Research Methods, 22*(3), 678–709. https://doi.org/10.1177/1094428118761122 | artigo | régua ,20/,40/,70 do dMACS (lida em fonte secundária; artigo não lido) | invariancia |
| Mantel, N., & Haenszel, W. (1959). Statistical aspects of the analysis of data from retrospective studies of disease. *Journal of the National Cancer Institute, 22*(4), 719–748. https://doi.org/10.1093/jnci/22.4.719 | artigo | `tri_dif(metodo = "mh")` | tri |
| Harrell, F. E., Lee, K. L., & Mark, D. B. (1996). Multivariable prognostic models: Issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Statistics in Medicine, 15*(4), 361–387. https://doi.org/10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4 | artigo | validação interna por bootstrap (otimismo) | associacao |
| Steyerberg, E. W., Harrell, F. E., Borsboom, G. J. J. M., Eijkemans, M. J. C., Vergouwe, Y., & Habbema, J. D. F. (2001). Internal validation of predictive models. *Journal of Clinical Epidemiology, 54*(8), 774–781. https://doi.org/10.1016/S0895-4356(01)00341-9 | artigo | bootstrap como validação interna preferida | associacao |
| Newcombe, R. G. (1998). Interval estimation for the difference between independent proportions: Comparison of eleven methods. *Statistics in Medicine, 17*(8), 873–890. https://doi.org/10.1002/(SICI)1097-0258(19980430)17:8<873::AID-SIM779>3.0.CO;2-I | artigo | IC da diferença de riscos no `qui_quadrado()` 2 × 2 | associacao |
| Altman, D. G. (1998). Confidence intervals for the number needed to treat. *BMJ, 317*(7168), 1309–1312. https://doi.org/10.1136/bmj.317.7168.1309 | artigo | IC do NNT (passa pelo infinito quando a DR cruza zero) | associacao |

Já na lista e reusadas na v1.6: Martinková e Hladká (2023), Petersen (2024), Timmerman e Lorenzo-Seva (2011), Enders e Tofighi (2007). O exemplo `housing` do MASS é banco de pacote, conferido contra os valores do próprio `polr`.

## v1.7: imputação múltipla, redes temporais, Mokken, CAT, normas, equivalência e invariância contínua

Conferidas em 29/10/2026 com `verify_citations.py` (todas VERIFIED, nenhuma retratação nem `updated-by`); autores, ano, periódico, volume e páginas do registro do Crossref.

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate imputation by chained equations in R. *Journal of Statistical Software, 45*(3). https://doi.org/10.18637/jss.v045.i03 | artigo (software) | `imputacao_multipla()` | fontes-e-decisoes |
| White, I. R., Royston, P., & Wood, A. M. (2011). Multiple imputation using chained equations: Issues and guidance for practice. *Statistics in Medicine, 30*(4), 377–399. https://doi.org/10.1002/sim.4067 | artigo | número de imputações (m ≥ % de casos incompletos) | fontes-e-decisoes |
| Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with multiple imputation. *Biometrika, 86*(4), 948–955. https://doi.org/10.1093/biomet/86.4.948 | artigo | gl da combinação de Rubin (o Crossref lista só o primeiro autor) | fontes-e-decisoes |
| Epskamp, S., Waldorp, L. J., Mõttus, R., & Borsboom, D. (2018). The Gaussian graphical model in cross-sectional and time-series data. *Multivariate Behavioral Research, 53*(4), 453–480. https://doi.org/10.1080/00273171.2018.1454823 | artigo | GVAR, redes temporal, contemporânea e entre pessoas; 100 pontos para 8 nós | redes-ega |
| Epskamp, S., van Borkulo, C. D., van der Veen, D. C., Servaas, M. N., Isvoranu, A.-M., Riese, H., & Cramer, A. O. J. (2018). Personalized network modeling in psychopathology: The importance of contemporaneous and temporal connections. *Clinical Psychological Science, 6*(3), 416–427. https://doi.org/10.1177/2167702617744325 | artigo | destendenciar antes de estimar a rede | redes-ega |
| Abegaz, F., & Wit, E. (2013). Sparse time series chain graphical models for reconstructing genetic networks. *Biostatistics, 14*(3), 586–599. https://doi.org/10.1093/biostatistics/kxt005 | artigo | estimação regularizada do graphicalVAR | redes-ega |
| Meijer, R. R., & Tendeiro, J. N. (2018). Unidimensional item response theory. In P. Irwing, T. Booth, & D. J. Hughes (Eds.), *The Wiley handbook of psychometric testing* (pp. 413–443). Wiley. https://doi.org/10.1002/9781118489772.ch15 | capítulo | Mokken: régua de H, ressalva sobre Hi < ,30 com itens Rasch | tri |
| van der Ark, L. A. (2012). New developments in Mokken scale analysis in R. *Journal of Statistical Software, 48*(5). https://doi.org/10.18637/jss.v048.i05 | artigo (software) | `escala_mokken()` | tri |
| Sijtsma, K., & van der Ark, L. A. (2017). A tutorial on how to do a Mokken scale analysis on your test and questionnaire data. *British Journal of Mathematical and Statistical Psychology, 70*(1), 137–158. https://doi.org/10.1111/bmsp.12078 | artigo | sequência da análise de Mokken | tri |
| Chalmers, R. P. (2016). Generating adaptive and non-adaptive test interfaces for multidimensional item response theory applications. *Journal of Statistical Software, 71*(5). https://doi.org/10.18637/jss.v071.i05 | artigo (software) | `simular_cat()` (mirtCAT) | tri |
| Counsell, A., Cribbie, R. A., & Flora, D. B. (2020). Evaluating equivalence testing methods for measurement invariance. *Multivariate Behavioral Research, 55*(2), 312–328. https://doi.org/10.1080/00273171.2019.1633617 | artigo | `equivalencia_ajuste()` | invariancia |
| Merkle, E. C., & Zeileis, A. (2013). Tests of measurement invariance without subgroups: A generalization of classical methods. *Psychometrika, 78*(1), 59–82. https://doi.org/10.1007/s11336-012-9302-4 | artigo | `invariancia_continua()` | invariancia |
| Merkle, E. C., Fan, J., & Zeileis, A. (2014). Testing for measurement invariance with respect to an ordinal variable. *Psychometrika, 79*(4), 569–584. https://doi.org/10.1007/s11336-013-9376-7 | artigo | funcionais WDMo e maxLMo | invariancia |
| Wang, T., Merkle, E. C., & Zeileis, A. (2014). Score-based tests of measurement invariance: Use in practice. *Frontiers in Psychology, 5*, 438. https://doi.org/10.3389/fpsyg.2014.00438 | artigo | uso aplicado, citado por Petersen (2024, 15.10.4.4.3) | invariancia |

Já na lista e reusadas na v1.7: Enders e Baraldi (2018), Isvoranu et al. (2022), Petersen (2024), Martinková e Hladká (2023), Chalmers (2012, mirt).

## v2.1: fator de Bayes, DIF por regressão logística generalizada e standard setting

Conferidas em 29/10/2026 com `verify_citations.py` (VERIFIED) e no registro do Crossref (nenhuma com `updated-by`).

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. (2009). Bayesian t tests for accepting and rejecting the null hypothesis. *Psychonomic Bulletin & Review, 16*(2), 225–237. https://doi.org/10.3758/PBR.16.2.225 | artigo | BF JZS do teste t, prior Cauchy, em `bf_teste_t()` | integral do BF10 conferida contra Monte Carlo em `evals_v21.R` |
| Ly, A., Verhagen, J., & Wagenmakers, E.-J. (2016). Harold Jeffreys's default Bayes factor hypothesis tests: Explanation, extension, and application in psychology. *Journal of Mathematical Psychology, 72*, 19–32. https://doi.org/10.1016/j.jmp.2015.06.004 | artigo | prior beta esticada (kappa) da correlação em `bf_correlacao()` | conferida contra a aproximação de Jeffreys em `evals_v21.R` |
| Drabinová, A., & Martinková, P. (2017). Detection of differential item functioning with nonlinear regression: A non-IRT approach accounting for guessing. *Journal of Educational Measurement, 54*(4), 498–517. https://doi.org/10.1111/jedm.12158 | artigo | modelo de `dif_nlr()` | eq. 9.9 de Martinková e Hladká (2023, p. 259), no acervo local |

## v2.2: Rasch por CML, CSEM, equacionamento e influência

Conferidas em 29/10/2026 no registro do Crossref.

| Referência (APA 7) | Tipo | Uso na skill | Âncora |
|---|---|---|---|
| Kolen, M. J., & Brennan, R. L. (2014). *Test equating, scaling, and linking: Methods and practices* (3rd ed.). Springer. https://doi.org/10.1007/978-1-4939-0317-7 | livro | continuização por posto percentil no equipercentil de `equacionar()` | citado por Cook e Pitoniak (2026, cap. 11) no acervo local |
| Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). *Regression diagnostics: Identifying influential data and sources of collinearity*. Wiley. https://doi.org/10.1002/0471725153 | livro | DFFITS e DFBETAS em `regressao_linear()` | réguas 2*sqrt(p/n) e 2/sqrt(n) declaradas como alerta |

## Citadas sem DOI

| Referência (APA 7) | Tipo | Uso na skill | Verificação |
|---|---|---|---|
| Cook, L. L., & Pitoniak, M. J. (Eds.). (2026). *Educational measurement* (5th ed.). Oxford University Press / National Council on Measurement in Education. | livro | Angoff modificado e Bookmark em `corte_angoff()` e `corte_bookmark()` (cap. 12, Tabela 12.2, pp. 862–863) (v2.1); CSEM em `csem()` (cap. 5, eq. 71–78, pp. 351–354) e equacionamento em `equacionar()` (cap. 11, eq. 15–16, pp. 789–790) (v2.2) | sem DOI conferido; cabeçalho e páginas do recorte `CookPitoniak2026__p0848-0921__Cap12_Standard_Setting.md` do acervo local |
| Thorndike, R. L. (1949). *Personnel selection: Test and measurement techniques*. Wiley. | livro | correção de restrição de amplitude, caso 2, em `corrigir_correlacao()` (v1.6) | sem DOI; entrada copiada da lista de referências de Martinková e Hladká (2023), que dão a fórmula (2.33) |
| Cassepp-Borges, V., Balbinotti, M. A. A., & Teodoro, M. L. M. (2010). Tradução e validação de conteúdo: Uma proposta para a adaptação de instrumentos. In L. Pasquali (Org.), *Instrumentação psicológica: Fundamentos e práticas* (pp. 506–520). Artmed. | capítulo | corte .70 do CVC (construcao-instrumentos) | sem DOI; dados conferidos só por citação secundária (artigo na Redalyc, revista 1938); ISBN do livro não conferido |
| *Cochrane handbook for systematic reviews of interventions* (Version 6.5). (2024). Cochrane. | livro (documento web) | faixas sobrepostas de I²; REML e HKSJ; k ≥ 10 para assimetria (modelos, sem citação formal) | sem DOI; seções 10.10.2, 10.10.4 e 13.3.4.4 lidas na pesquisa de 11/10/2026; editores e URL não conferidos nesta sessão, completar antes de citar |
| *G\*Power 3.1 manual* (versão de 1 de junho de 2023). Documento do programa G\*Power. | manual (documento web) | exemplos conferidos nos evals: logística (seção 29.4), alocação com razão diferente de 1 (seções 15 e 17), medidas repetidas (poder-amostra) | sem DOI; texto completo lido em 11 e 14/10/2026; autoria e URL não conferidas, completar antes de citar |
| Croon, M. (2002). Using predicted latent scores in general latent structure models. In G. A. Marcoulides & I. Moustaki (Eds.), *Latent variable and latent structure models* (pp. 195–223). Lawrence Erlbaum. | capítulo (epônimo) | correção de Croon para escore fatorial (afc-ordinal, `checar_escore_somado()`) | sem DOI da edição de 2002; ano, páginas e editora pela lista de McNeish e Wolf (2020), que grafa "modeling" no título do livro; título e editores pelo Crossref da reedição da Psychology Press (2014), https://doi.org/10.4324/9781410602961-16, pp. 207–236, capítulo sem autor no registro; conteúdo não lido |
| Girden, E. R. (1992). *ANOVA: Repeated measures* (Sage University Paper Series on Quantitative Applications in the Social Sciences, No. 84). Sage. | livro | corte ,75 do ε citado com Barcikowski e Robey (fontes-e-decisoes) | sem DOI; dados da série QASS, ISBN não conferido |
| Hayes, A. F. (2022). *Introduction to mediation, moderation, and conditional process analysis: A regression-based approach* (3rd ed.). Guilford Press. | livro | numeração dos modelos PROCESS, índice de mediação moderada (processos) | sem DOI; ISBN 978-1-4625-4903-0 conferido no Open Library (título, 3ª ed., 2022; lá a editora aparece como Taylor & Francis Group); Guilford Press como no livro do usuário |
| Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70. | artigo | Holm como correção padrão (diferencas) | sem DOI no Crossref (busca bibliográfica de 11/10/2026 sem resultado); JSTOR 4615733 |
| Tabachnick, B. G., & Fidell, L. S. (2019). *Using multivariate statistics* (7th ed.). Pearson. | livro | Mahalanobis com χ² crítico a p < .001 (associacao) | sem DOI; ISBN 978-0-13-479054-1 conferido no Open Library (Pearson, publicado em 12/07/2018; 2019 é o ano citado no arquivo); página não conferida |
| Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct reliability within latent variable systems. In R. Cudeck, S. du Toit, & D. Sörbom (Eds.), *Structural equation modeling: Present and future* (pp. 195–216). Scientific Software International. | capítulo | H (replicabilidade do construto) ≥ ,70 em `indices_bifator()` | sem DOI; dados do capítulo pela citação secundária em Rodriguez et al. (2016a); páginas não conferidas na obra |
| Gorsuch, R. L. (1983). *Factor analysis* (2nd ed.). Lawrence Erlbaum Associates. | livro | determinação fatorial ≥ ,90 para usar escore fatorial (afc-ordinal) | sem DOI; citação secundária (Rodriguez et al., 2016a) |
| Shrout, P. E., & Lane, S. P. (2012). Psychometrics. In M. R. Mehl & T. S. Conner (Eds.), *Handbook of research methods for studying daily life*. Guilford Press. | capítulo | coeficientes RkF, R1R, RkR e Rc de `confiabilidade_ema()` | sem DOI; documentação de `psych::multilevel.reliability`; páginas não conferidas |
| Linacre, J. M. (2002). What do infit and outfit, mean-square and standardized mean? *Rasch Measurement Transactions, 16*(2), 878. | artigo (boletim) | faixa 0,5 a 1,5 do MNSQ em `tri_ajustar()` | sem DOI; boletim da Rasch SIG, citado no curso |
| Zumbo, B. D. (1999). *A handbook on the theory and methods of differential item functioning (DIF): Logistic regression modeling as a unitary framework for binary and Likert-type (ordinal) item scores*. Directorate of Human Resources Research and Evaluation, Department of National Defense. | relatório técnico | DIF logístico com ΔR² (tri) | sem DOI; citado pela documentação do difR |
| Cook, L. L., & Pitoniak, M. J. (Eds.). (2026). *Educational measurement* (5th ed.). Oxford University Press. | livro | blueprint e equidade (construcao-instrumentos, invariancia) | sem DOI no mapa da pasta Psicometria-R; texto completo na pasta |

## Não verificadas

| Citação | Arquivo | O que foi achado |
|---|---|---|
| Hernández-Nieto (2002), livro do CVC | construcao-instrumentos | Citações secundárias dão "Hernández-Nieto, R. A. (2002). *Contribuciones al análisis estadístico*. Universidad de Los Andes", com grafias variadas do título; existe uma edição em inglês (*Contributions to statistical analysis: The coefficients of proportional variance, content validity and kappa*, BookSurge, ISBN 978-1-58898-715-0) sem ano confirmado. Sem DOI. Não usar sem conferir a ficha do exemplar. |
| Cai e Monro (2014), no `M2` do mirt | nenhum arquivo atual (lista anterior) | Relatório técnico do CRESST citado pela documentação do mirt; sem registro no Crossref (busca de 11/10/2026). |
