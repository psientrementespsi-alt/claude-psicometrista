# CHANGELOG, higiene-dados

## 1.5.0 (04/11/2026)

- Novo `scripts/ajuste_pessoa.R`: `calc_guttman_poly()` (G^p com passos de item e popularidade condicional, Apêndice A de Niessen, Meijer & Tendeiro, 2016), `calc_lz_poly()` (lz^p sob o GRM; igual ao Zh do `mirt::personfit()` com o mesmo theta) e `calc_ajuste_pessoa()` (soma nas subescalas e corte empírico pelo quantil da amostra).
- `triagem_cier_completa()` ganhou `flag_ajuste_pessoa`: o ajuste de pessoa entra como sétima família no critério de 2 ou mais.
- `validacao_codebook.R` desativado: definia outra `validar_regras()` (exigia lista) e outra `gerar_codebook()`; dar `source()` em todos os scripts por ordem alfabética quebrava o pipeline. Agora só carrega `auditoria_higiene.R`.
- `references/04_guia_deteccao_cier_e_outliers.md`: o corte de tempo dizia 1,5 s por item, em desacordo com a skill; passa a 2 s por item (Huang et al., 2012, citado por Curran, 2016; Ward & Meade, 2023).
- `evals_higiene.R`: 37 testes (6 novos), todos passando no Mac.
- SKILL.md unificado em 04/11/2026: o texto do Mac recebeu as linhas que só a versão da conta do Claude tinha (rotas para as outras skills, notas de uso), e a conta recebe as desta versão.

## 1.4.0 (29/10/2026)

- `fluxograma_mermaid()` em `auditoria_higiene.R`: gerador automático de diagrama de fluxo de participantes em sintaxe Mermaid (JARS-Quant / CONSORT), integrado ao retorno de `executar_pipeline_higiene()` como `fluxograma_mermaid`.
- `exportar_lavaan_mi()` em `imputacao_mice.R`: exportador estruturado de objetos `mids` para listas de data.frames prontos para `semTools::runMI()` e `lavaan.mi` com pooling de Fração de Informação Ausente (FMI; Rubin, 1987; Enders, 2022).
- `evals_higiene.R`: expandido para 31 testes unitários (100% aprovados).

## 1.3.0 (29/10/2026)

- `executar_pipeline_higiene()` v2: trilha de auditoria por etapa, Tabela 1 e parágrafo de Método gerados da trilha, critérios de inclusão por regra (`regras_excluir`), FastMCD sobre subescalas, execução em blocos (`salvar_em`, `etapas_ate`, `retomar_de`) com os parâmetros da primeira execução preservados na retomada. Duplicatas só entre registros com resposta; protocolo todo no mínimo protegido sempre; aviso quando há menos de três colunas além dos itens. A v1 ignorava o argumento `alpha_mcd`, não registrava critério de inclusão e montava a Tabela 1 à mão.
- Novas famílias C/IER: Mahalanobis clássico nos itens (`calc_mahad_cier()`, por fator com Bonferroni; Meade & Craig, 2012; DeSimone et al., 2015) e itens de atenção ou infrequência (`calc_falhas_atencao()`; Huang et al., 2015).
- Novo `auditoria_higiene.R`: `validar_regras()` (validate quando instalado, R base quando não), `gerar_codebook()` só com agregados e supressão de células < 5, `registrar_passo()`, `tabela_fluxo_md()`, `texto_metodo_higiene()`.
- `evals_higiene.R`: 29 testes; reproduz N = 116 e 5 atípicos MCD75 do banco DASS-21 do projeto.

## 1.2.0 (29/10/2026)

- Triagem C/IER por famílias: longstring e IRV = 0 viram uma família (repetição), e a exclusão exige duas famílias. A 1.1.x contava IRV = 0 e longstring como dois índices, o que excluía por um comportamento só.
- Novo `calc_evenodd()` (igual ao `careless::evenodd`, com subescalas não contíguas) e psychsyn com corte < 0 (Curran, 2016) entram como famílias; índice que assinala mais de 20% dos casos calculáveis sai do critério com aviso.
- Cortes com fonte: longstring em metade da escala (Curran, 2016), 2 s por item (Huang et al., 2012).
- `calc_mahad_robust()`: MCD75 por padrão (Leys et al., 2019), argumento `alpha_mcd`; grade sobe com singularidade.
- `evals_higiene.R`: 18 testes.

## 1.1.1 (29/10/2026)

- SKILL.md: Little com p > .05 deixou de ser lido como "MCAR plausível"; não rejeitar não prova MCAR (alinhado à psicometria-r).
- `recodificar_likert_labels()`: sem `mapa_niveis`, rotulos sem numero param com erro listando os rotulos. A 1.1.0 caia em `as.numeric(as.factor())`, que ordena em ordem alfabetica (na DASS-21, "Nao aconteceu comigo" virava 4). Rotulos fora do mapa viram NA com aviso, e os nomes das colunas sao preservados.
- `teste_little_mcar()`: mu e Sigma por EM e todos os padroes de ausencia na estatistica, inclusive os de um so caso. A 1.1.0 usava estimativas pairwise e pulava padroes com n_j < 2, devolvendo chi2 = 0, gl = 1, p = 1 quando cada omissao era de um participante (DASS-21: o correto e chi2(40) = 25.92, p = .958).
- Novo `scripts/evals_higiene.R`.

