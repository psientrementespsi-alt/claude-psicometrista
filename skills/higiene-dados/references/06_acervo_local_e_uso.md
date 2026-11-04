# Acervo local da higiene-dados: o que existe, onde e para que serve

Pasta no Mac: `~/Desktop/psicomtria rrr/` (subpastas `papers-higiene-dados/`, `repositorios-higiene-dados/`, `referencial-higiene-dados/`). Inventário de 29/10/2026.

## Papers (30 obras, PDF e Markdown; lista completa em `papers-higiene-dados/manifest.md`)

| Decisão da skill | Fonte no acervo (arquivo em `markdown/`) | O que a fonte sustenta |
|---|---|---|
| Longstring em metade da escala | Curran (2016), `methods-for-the-detection-of-carelessly-invalid-re__10_1016_j_jesp_2015_07_006.md` | Regra de base declarada pelo autor como ponto de partida, dependente da escala |
| Tempo < 2 s por item | Curran (2016), mesmo arquivo, citando Huang et al. (2012) | Corte conservador de tempo |
| Psychsyn < 0 | Curran (2016), mesmo arquivo | Correlação intrapessoal negativa nos sinônimos como sinal quase certo |
| Even-odd, psychsyn, combinação de índices | Meade e Craig (2012), `identifying-careless-responses-in-survey-data__10_1037_a0028085.md`; Ward e Meade (2023) | Definição dos índices e uso combinado |
| MAD com 3 (2,5 moderado) e MCD75 a p = .001 | Leys et al. (2019), `how-to-classify-detect-and-manage-univariate-and-m__10_5334_irsp_289.md`; Leys et al. (2013) | MCD75 como padrão, MCD50 com suspeita de mais de 25% de atípicos; detectar uma vez, sem reaplicar |
| MCD e extensões | Hubert, Debruyne e Rousseeuw (2018), `10_1002_wics_1421.md` | Base do FastMCD reponderado; o artigo usa corte em χ² a .975 para a reponderação, a skill usa p = .001 no corte final (Leys et al., 2019) |
| Mahalanobis como C/IER | DeSimone, Harms e DeSimone (2015), `10_1002_job_1962.md`; Meade e Craig (2012) | D² nos itens contra χ² com k gl; DeSimone et al. sugerem os 5% superiores e lembram que a distribuição supõe itens normais (a skill usa .001); Meade e Craig fazem um D por fator quando há muitos itens |
| Ajuste de pessoa G^p e lz^p (v1.5) | Niessen, Meijer e Tendeiro (2016), `10_1016_j_jrp_2016_04_010.md` | Definição de G^p por passos de item (Apêndice A, Tabelas A1 e A2, reproduzidas no eval), lz^p sob o GRM, extensão multiteste por soma nas subescalas (seção 3.3), cortes empíricos e o aviso de relatar com e sem os excluídos |
| Escores de atípico por erros de Guttman | Zijlstra, van der Ark e Sijtsma (2007), `10_1080_00273170701384340.md` | Erros de Guttman ponderados em itens politômicos como escore de atípico; base conceitual, não implementado |
| Itens de infrequência | Huang et al. (2015), `10_1007_s10869-014-9357-6.md` | Escala de infrequência como indicador de IER, com evidência de validade e reação dos participantes |
| Qualidade em coleta online | Chmielewski e Kucker (2020), `an-mturk-crisis-shifts-in-data-quality-and-the-imp__10_1177_1948550619875149.md` | Taxas de dado inválido mudam entre ondas de coleta e alteram resultados substantivos; justifica rodar a triagem em toda coleta online |
| Atípicos: decidir e relatar | Aguinis et al. (2013), `best-practice-recommendations-for-defining-identif__10_1177_1094428112470848.md` | Definir, identificar e tratar com critério declarado e sensibilidade |
| Ausentes | Schafer e Graham (2002), `10_1037_1082-989x_7_2_147.md`; Lang e Little (2018), `principled-missing-data-treatments__10_1007_s11121-016-0644-5.md`; Tierney e Cook (2023, naniar), `expanding-tidy-data-principles-to-facilitate-missi__10_18637_jss_v105_i07.md`; resenha de Khademi (2020) do livro de van Buuren, `flexible-imputation...md` (o livro em si não está no acervo); Little (1988) fora do acervo | MCAR/MAR, imputação múltipla, FIML |
| Validação declarativa e organização | van der Loo e de Jonge (2021, `validate`), `data-validation-infrastructure-for-r__10_18637_jss_v097_i10.md`; Broman e Woo (2018), `data-organization-in-spreadsheets...md`; Wickham (2014), `tidy-data...md` | Regras de dados, planilhas e dados arrumados |
| Codebook | Arslan (2019), `how-to-automatically-document-data-with-the-codebo...md` | Documentação automática do banco |

## Repositórios (código-fonte de pacotes R; nenhum precisa estar instalado para a skill funcionar)

| Pasta | Pacote | Uso na skill |
|---|---|---|
| `04_triagem_comportamental_e_careless/careless` | careless | `calc_evenodd()` reproduz `evenodd()` (eval com tolerância 1e-10); `psychsyn()` usa `> critval` e reamostra NA, diferença documentada |
| `05_dados_ausentes_e_imputacao/naniar` | naniar | Referência do `mcar_test()`; o eval de igualdade roda quando o naniar estiver instalado |
| `05_dados_ausentes_e_imputacao/mice`, `simputation` | mice, simputation | Imputação múltipla (`executar_mice_psicometria()`) |
| `02_validacao_integridade_e_qualidade` | validate, pointblank, assertr | `validar_regras()` usa o `validate` quando instalado; senão avalia as mesmas expressões em R base |
| `01_limpeza_e_sanitizacao_geral` | janitor, skimr, datawizard | Nomes, vazios e descritivas; a skill tem `sanitizar_nomes_snake()` em base R |
| `03_rastreamento_e_auditoria_logs` | lumberjack, tidylog | Inspiração de `registrar_passo()`: a trilha do pipeline v2 registra n antes, n depois, critério e fonte de cada etapa |
| `06_diagnostico_e_dicionarios_metadados` | codebook, dlookr | Inspiração de `gerar_codebook()`, que só traz agregados e suprime células < 5 |

## O que fica fora

Os relatórios de auditoria em `referencial-higiene-dados/` (Alpha, Beta, Gamma, CoVe) são sínteses geradas para construir a skill, não fontes primárias. Afirmação de método cita o paper, nunca esses relatórios.
