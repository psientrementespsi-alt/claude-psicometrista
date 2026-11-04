# BLUEPRINTS DE CÓDIGO R E TEMPLATES DE ESCRITA APA 7
## SCRIPTS DE GERAÇÃO TABULAR E FRASÁRIO PADRONIZADO PARA O CLAUDE

**Destinatário:** Claude (Anthropic) / Agentes Redatores de Artigos Científicos  
**Localização:** `~/Desktop/psicomtria rrr/referencial-relato-estatistico/`  
**Objetivo:** Fornecer blocos de código em R que geram tabelas formatadas em APA 7 e parágrafos de texto prontos para publicação científica.

---

## 1. TABELA DE CARACTERIZAÇÃO DEMOGRÁFICA (APA 7) COM `gtsummary`

```r
library(gtsummary)
library(tidyverse)

# Tabela 1: Características Demográficas e Clínicas da Amostra
tabela_demografica <- dados_purificados |>
  select(idade, sexo, escolaridade, estado_civil, grupo_clinico) |>
  tbl_summary(
    by = grupo_clinico, # Estratificação por grupo (se houver)
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 2, all_categorical() ~ c(0, 1)),
    label = list(
      idade ~ "Idade em anos, M (SD)",
      sexo ~ "Sexo Biológico, n (%)",
      escolaridade ~ "Nível de Escolaridade, n (%)",
      estado_civil ~ "Estado Civil, n (%)"
    ),
    missing = "no" # Omissos relatados em nota de rodapé
  ) |>
  add_overall(label = "**Amostra Total (N = {N})**") |>
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 3)) |>
  bold_labels() |>
  italicize_levels() |>
  modify_header(label = "**Variável**") |>
  modify_footnote(
    all_stat_cols() ~ "Dados contínuos apresentados como Média (Desvio-Padrão); categóricos como Contagem (Porcentagem)."
  )

# Exportação para documento Word (.docx) ou HTML
# tabela_demografica |> as_flex_table() |> flextable::save_as_docx(path = "Tabela1_Demografia_APA7.docx")
```

---

## 2. TABELA DE COMPARAÇÃO DE MODELOS DE REGRESSÃO COM `modelsummary`

```r
library(modelsummary)

# Ajuste de Modelos Hierárquicos Lineares
modelo_1 <- lm(satisfacao ~ idade + sexo, data = dados_purificados)
modelo_2 <- lm(satisfacao ~ idade + sexo + suporte_social, data = dados_purificados)
modelo_3 <- lm(satisfacao ~ idade + sexo + suporte_social + neuroticismo, data = dados_purificados)

lista_modelos <- list(
  "Modelo 1 (Demografia)" = modelo_1,
  "Modelo 2 (+ Suporte)"   = modelo_2,
  "Modelo 3 (+ Personalidade)" = modelo_3
)

# Tabela 2 no Padrão APA 7
tabela_regressao <- modelsummary(
  lista_modelos,
  output = "flextable", # Saída em flextable pronta para exportação Word
  estimate = "{estimate} ({std.error}){stars}",
  statistic = "conf.int",
  coef_rename = c(
    "(Intercept)" = "Intercepto",
    "idade" = "Idade (anos)",
    "sexoFeminino" = "Sexo (Feminino)",
    "suporte_social" = "Suporte Social",
    "neuroticismo" = "Neuroticismo"
  ),
  gof_map = c("nobs", "r.squared", "adj.r.squared", "fstatistic", "p.value", "aic"),
  stars = c('*' = .05, '**' = .01, '***' = .001),
  title = "Tabela 2. Coeficientes de Regressão Linear Hierárquica para Predição da Satisfação com a Vida"
)
```

---

## 3. TABELA DE CARGAS FATORIAIS E ÍNDICES DE AJUSTE (CFA) COM `rempsyc`

```r
library(lavaan)
library(rempsyc)

# Ajuste do Modelo Fatorial Confirmatório
modelo_cfa <- '
  Fator_1 =~ item_01 + item_02 + item_03 + item_04 + item_05
  Fator_2 =~ item_06 + item_07 + item_08 + item_09 + item_10
  Fator_3 =~ item_11 + item_12 + item_13 + item_14 + item_15
'

fit_cfa <- cfa(modelo_cfa, data = dados_purificados, ordered = TRUE, estimator = "WLSMV")

# 1. Extração de Cargas Fatoriais Padronizadas e Erros-Padrão
cargas_fatoriais <- standardizedSolution(fit_cfa) |>
  filter(op == "=~") |>
  select(Fator = lhs, Item = rhs, Carga_Padronizada = est.std, SE = se, z = z, p = pvalue) |>
  mutate(across(where(is.numeric), ~ round(., 3)))

# 2. Tabela APA 7 Elegante via rempsyc
tabela_cfa_apa <- nice_table(
  cargas_fatoriais,
  title = c("Tabela 3", "Cargas Fatoriais Padronizadas da Análise Fatorial Confirmatória (WLSMV)"),
  note = c("*Nota.* Todos os parâmetros são significativos a p < .001. Estimação robusta com matriz policórica.")
)
```

---

## 4. TRADUÇÃO DIRETA PARA PROSA CIENTÍFICA COM `easystats/report`

```r
library(report)

# Geração de Relato Textual Automático para Testes Paramétricos
teste_t <- t.test(escore ~ grupo, data = dados_purificados, var.equal = FALSE)
relato_t <- report(teste_t)

# Imprime o parágrafo completo no console
print(relato_t)

# Geração de Relato para Modelos Lineares Múltiplos
relato_reg <- report(modelo_3)
print(relato_reg)
```

---

## 5. AUDITORIA DE CONSISTÊNCIA DE TEXTO COM `MicheleNuijten/statcheck`

Antes de submeter o manuscrito ou finalizar o relatório, o Claude pode auditar o texto da seção de Resultados para garantir que nenhuma estatística contém erro de digitação entre $t$, $F$, graus de liberdade e $p$-valor:

```r
library(statcheck)

# Audita uma string ou arquivo de texto diretamente
texto_resultados <- "O efeito do tratamento foi significativo, t(48) = 2.45, p = .018, indicando melhora clínica."

resultado_auditoria <- statcheck(texto_resultados)
print(resultado_auditoria)

# Se reported_p != computed_p, o statcheck aponta um erro de consistência interna
if (any(resultado_auditoria$error == TRUE)) {
  warning("ALERTA: Detectada inconsistência entre o valor de teste reportado e o p-valor!")
}
```


---

## 6. REDES PSICOMÉTRICAS COM `bootnet` (v1.5.0)

```r
net <- bootnet::estimateNetwork(dados, default = "EBICglasso", tuning = 0.5)
caso <- bootnet::bootnet(net, nBoots = 1000, type = "case", statistics = c("strength", "expectedInfluence"))
cs <- bootnet::corStability(caso)          # CS(cor = 0,7) por índice
frases_rede(net, cs, perfil = "apa7_ptbr")  # frase da rede e frase da estabilidade
```

O texto sai do objeto: número de nós, γ do EBIC, N, arestas diferentes de zero sobre as possíveis, densidade e sinais. A estabilidade segue Epskamp, Borsboom e Fried (2018): CS abaixo de 0,25 não sustenta interpretar a ordem de centralidade; acima de 0,50 é o preferido. Os próprios autores chamam os cortes de algo arbitrários, então a frase atribui o corte à fonte. Bootstrap de arestas não é teste de significância de aresta (mesmo artigo); não escrever que uma aresta "foi significativa" porque o IC não inclui zero.
