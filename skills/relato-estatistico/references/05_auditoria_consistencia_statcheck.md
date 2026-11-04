# Auditoria de Consistência Estatística via Statcheck

Este documento define o protocolo de auditoria algorítmica para verificação de consistência matemática em testes de hipótese nula (*Null Hypothesis Significance Testing* - NHST), baseado na metodologia pioneira de Nuijten et al. (2016).

---

## 1. O Problema das Inconsistências de Relato

Estudos meta-científicos em psicologia e neurociência comportamental demonstraram que cerca de metade de todos os artigos publicados em periódicos revisados por pares contêm pelo menos uma inconsistência entre o valor da estatística de teste reportada ($t, F, \chi^2, r$) com seus respectivos graus de liberdade ($df$) e o valor de probabilidade associado ($p$).
Aproximadamente 12% dos artigos contêm **erros grosseiros** (*gross inconsistencies*), em que o resultado é declarado como estatisticamente significativo ($p < .05$), mas o recálculo exato a partir da estatística revela que o teste não atingiu significância ($p \ge .05$), ou vice-versa.

---

## 2. Padrões Sintáticos Reconhecidos

O auditor algorítmico procura sentenças e parágrafos contendo as seguintes estruturas canônicas:

1. **Teste t:**
   - Expressão regular / padrão: `t(df) = valor, p = valor` ou `t(df) = valor, p < valor`
   - Cálculo analítico de verificação no R:
     `p_calculado = 2 * (1 - pt(abs(t_valor), df = df))` (bicaudal padrão)
2. **Teste F:**
   - Expressão regular: `F(df1, df2) = valor, p = valor`
   - Cálculo analítico:
     `p_calculado = 1 - pf(F_valor, df1 = df1, df2 = df2)`
3. **Teste Qui-Quadrado ($\chi^2$):**
   - Expressão regular: `chi2(df) = valor, p = valor` ou `chi^2(df) = valor, p = valor`
   - Cálculo analítico:
     `p_calculado = 1 - pchisq(chi2_valor, df = df)`
4. **Correlação de Pearson ($r$):**
   - Expressão regular: `r(df) = valor, p = valor`
   - Transformação e teste:
     `t_stat = r_valor * sqrt(df / (1 - r_valor^2))`
     `p_calculado = 2 * (1 - pt(abs(t_stat), df = df))`

---

## 3. Classificação de Discrepâncias

Para cada teste auditado, o algoritmo compara o $p$ reportado ($p_{\text{rep}}$) com o $p$ computado ($p_{\text{calc}}$):

1. **Consistente (OK):**
   - A diferença $|p_{\text{rep}} - p_{\text{calc}}|$ é menor que a margem de arredondamento esperado (ex.: diferença $< 0.01$ para relato com 2 casas decimais, ou $< 0.001$ para 3 casas decimais).
   - Não há divergência na decisão inferencial.
2. **Inconsistência Leve (Minor Inconsistency):**
   - O valor reportado diverge do recalculado por erro de arredondamento ou transcrição (ex.: texto relata $p = .03$, mas o cálculo resulta em $p = .042$).
   - A decisão substantiva permanece inalterada (ambos são $< .05$).
3. **Inconsistência Grosseira (Gross Inconsistency):**
   - A divergência altera a conclusão do teste de hipótese no limiar padrão $\alpha = .05$.
   - Caso A: Texto reporta $p < .05$ (ou $p = .04$), mas a estatística e $df$ produzem $p = .08$.
   - Caso B: Texto reporta $p > .05$ (ou $p = .07$), mas a estatística produz $p = .02$.
   - **Ação:** Bloqueio imediato do texto. O agente deve identificar a origem do erro no modelo R e corrigir o valor no texto antes da apresentação ao usuário.
