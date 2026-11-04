# CHANGELOG, relato-estatistico

## 1.5.0 (04/11/2026)

- `frases_rede()` em `frases_apa.R`, e `frases_apa()` passa a aceitar `bootnetResult`: frase da rede (nós, estimador e γ, N, arestas, densidade, sinais, peso médio) e da estabilidade pelo coeficiente CS, com os cortes de 0,25 e 0,50 de Epskamp, Borsboom e Fried (2018).
- `conferir_d_de_t(hedges = TRUE)` confere g_s (Lakens, 2013, eq. 4); nova `conferir_omega2p_de_F()` (Lakens, 2013, eq. 15).
- `evals_relato.R`: 30 testes (4 novos, ancorados no exemplo dos filmes de Lakens e no exemplo de CS de Epskamp et al.), todos passando no Mac.
- SKILL.md unificado em 04/11/2026: o texto do Mac recebeu as linhas que só a versão da conta do Claude tinha (rotas para as outras skills, notas de uso), e a conta recebe as desta versão.

## 1.4.0 (29/10/2026)

- Extensão de `frases_apa.R` com handlers genéricos para classes externas:
  - `glm` (regressão logística e Poisson): relata OR ou RP com IC 95%, coeficientes B, EP, z e p.
  - `lmerMod` / `glmerMod` (`lme4`): modelos lineares hierárquicos com ICC, variâncias de intercepto randômico (τ₀₀), resíduo (σ²), efeitos fixos com EP, t/z, p e IC 95%.
  - `coxph` (`survival`): modelos de sobrevida com Hazard Ratio (HR), IC 95%, B, EP, z e p.
  - `rma.uni` / `rma` (`metafor`): meta-análise com tamanho de efeito combinado, IC 95%, z, p, e estatísticas completas de heterogeneidade (Q, df, p, I², τ²).
  - `roc` (`pROC`): curvas ROC com AUC, IC 95% e ponto de corte ótimo de Youden (com sensibilidade e especificidade).
- `evals_relato.R`: 26 testes unitários (100% aprovados).

## 1.3.0 (29/10/2026)

- Novo `formatador_apa.R`: zero à esquerda e decimal por perfil escolhido para o manuscrito (`apa7_relato`, `apa7_estrito`, `psicometria_r`, `apa7_ptbr`, `apa7_ptbr_sem_zero`, `rbtc`, `avaliacao_psicologica`, `prc`, personalizado). As instruções aos autores da RBTC, da Avaliação Psicológica e da PRC consultadas em 29/10/2026 não fixam regra de zero nem de decimal; esses perfis herdam um perfil-base e saem marcados como não verificados.
- Novo `frases_apa.R`: frases de t.test, cor.test, chisq.test, lm e lavaan geradas do objeto, compatíveis com o auditor.
- Novo `checagem_efeitos.R`: efeito e IC por frase; d de t e η²p de F reconstituídos (Lakens, 2013).
- Novo `docx_apa.R`: tabela APA em .docx com R base (OOXML e utils::zip), aberta sem erro pelo textutil do macOS.
- `auditar_texto_statcheck()`: aceita ponto e vírgula entre a estatística e o p, e entre os gl do F (texto com vírgula decimal).
- `gerar_tabela_correlacoes_apa()`: a nota relata a faixa de n por par quando há ausentes; antes relatava o nrow do banco.
- `evals_relato.R`: 22 testes.

## 1.2.0 (29/10/2026)

- `auditar_texto_statcheck()` reescrita com a regra do pacote statcheck (intervalo de arredondamento da estatística e casas decimais do p), lendo t, F, χ², r e z em inglês e português. A 1.1.x usava tolerância fixa de ,015, não lia r(gl) nem vírgula decimal.
- Novo `evals_relato.R`: 8 testes, com os veredictos de 24 casos conferidos contra o statcheck do repositório.

## 1.1.1 (29/10/2026)

- `gerar_tabela_afc_apa()`: nota com uma unica versao de indices (robusta, senao escalonada), declarada; CFI e TLI com tres casas sem zero; RMSEA (com IC 90%) e SRMR com tres casas e com zero, como manda a regra 2 do SKILL.md (a 1.1.0 imprimia RMSEA = .08 e SRMR = .10).
- `gerar_tabela_confiabilidade_apa()`: com a saida delta da psicometria-r 1.8.1, mostra omega do escore observado, omega delta com EP e IC, e alfa com IC de Feldt, em tres casas.

