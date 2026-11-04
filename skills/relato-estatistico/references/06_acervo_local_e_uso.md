# Acervo local da relato-estatistico: o que existe, onde e para que serve

Pasta no Mac: `~/Desktop/psicomtria rrr/` (subpastas `repositorios-relato-estatistico/`, `papers-relato-e-tratamento-markdown/` e `-pdf/`, `referencial-relato-estatistico/`). Inventário de 29/10/2026.

## Repositórios (código-fonte; pacotes não instalados no R do Mac em 29/10/2026)

| Pasta | Pacote | Uso |
|---|---|---|
| `statcheck` | statcheck (Nuijten et al., 2016) | Regra de decisão portada em `auditor_statcheck.R`; `evals_relato.R` carrega `statcheck/R` e compara os veredictos |
| `papaja` | papaja | `apa_num(gt1 = FALSE)` é a convenção de zero à esquerda do pacote: sem zero só quando a estatística não pode passar de 1; `apa_print()` para frases de lm, aov, htest |
| `apaTables` | apaTables | `apa.cor.table()` e `apa.reg.table()`: tabelas APA em .doc, com IC de r e de B e IC de ΔR² |
| `rempsyc` | rempsyc | `nice_table()`: tabela APA em Word via flextable, a saída com três linhas horizontais pronta para .docx |
| `modelsummary`, `gtsummary` | modelsummary, gtsummary | Tabelas comparativas de modelos e descritivas por grupo para Word, HTML, LaTeX |
| `report` | report (easystats) | Frases automáticas para lm, glm, htest e lavaan; servem de conferência, não de texto final (vocabulário genérico) |
| `equatiomatic` | equatiomatic | Equação do modelo em LaTeX a partir do objeto do R |

Instalar qualquer um deles no R do Mac exige pedido explícito do usuário (regra da psicometria-r).

## Papers (Markdown e PDF)

| Arquivo | Uso |
|---|---|
| `Nuijten_Statcheck_2016.md` | Base do auditor e da definição de inconsistência e erro grosseiro |
| `Lakens_Effect_Sizes_2013.md` | Escolha e relato de tamanhos de efeito (d, g, η²p, ω²); o exemplo dos filmes (t(18) = 2,52; d_s = 1,13; g_s = 1,08; ω²p = ,21) é eval de `conferir_d_de_t(hedges = TRUE)` e `conferir_omega2p_de_F()` |
| `rempsyc_Theriault_2023.md`, `gtsummary_Sjoberg_2021.md`, `modelsummary_ArelBundock_2022.md`, `datawizard_Patil_2022.md` | Documentação dos pacotes de tabela e de preparação |
| `Epskamp_Network_Tutorial_2018.md` | Epskamp, Borsboom e Fried (2018, *Behavior Research Methods, 50*, 195–212, https://doi.org/10.3758/s13428-017-0862-1). Relato de redes (estimador, γ, estabilidade) e os cortes de CS de 0,25 e 0,50 usados em `frases_rede()`; o exemplo de CS do artigo (força ,44; intermediação e proximidade ,05) é eval; em linha com `psicometria-r/references/redes-ega.md` |

## Convenções que divergem entre as skills (decisão do autor)

- Zero à esquerda em RMSEA e SRMR: esta skill usa zero (regra 2); `psicometria-r/references/relato-apa.md` omite. O papaja decide por `gt1`: RMSEA pode, em princípio, passar de 1.
- β: esta skill omite o zero (APA 7); `psicometria-r/references/associacao.md`, seguindo o livro do usuário, põe zero porque β pode passar de 1.
- Ômega: com itens ordinais, relatar o ômega do escore observado (Green & Yang, 2009) com IC por bootstrap, como manda `psicometria-r/references/afc-ordinal.md`; `omega_delta` é da variável-resposta latente.
