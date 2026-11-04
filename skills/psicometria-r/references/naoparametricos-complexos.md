# Não paramétricos para delineamentos complexos (est_naoparametricos.R, v1.5)

Fonte: curso "R: Testes não paramétricos para delineamentos complexos" da pasta
Psicometria-R (pasta 06 e cópia na 08), com os scripts `rankFD.R`, `nparLD.R` e
`MANOVARM.R` e as descrições de resultado em .docx. Base teórica: Brunner et al.
(2017), Noguchi et al. (2012), Konietschke et al. (2012, 2015), Friedrich e Pauly (2018).

Sumário: 1. Quando usar · 2. Efeito relativo · 3. `np_fatorial()` · 4. `np_longitudinal()`
· 5. `manova_np()` · 6. Divergências com as aulas · 7. Relato

## 1. Quando usar

Os testes de Kruskal-Wallis e Friedman (`diferencas.md`, seção 10) cobrem um fator. Com
dois ou mais fatores, medidas repetidas cruzadas com grupos ou várias VDs, e com dados
ordinais, assimétricos ou com atípicos que a ANOVA não aguenta, a rota do curso é a
análise por postos de Brunner: não supõe normalidade nem variâncias iguais e testa
hipóteses sobre distribuições (H0F), não sobre médias.

| Delineamento | Função | Pacote |
|---|---|---|
| 1 ou mais fatores entre sujeitos | `np_fatorial(dados, vd, fatores)` | rankFD + nparcomp::mctp |
| tempo (LD-F1), tempo × grupo (F1-LD-F1), tempo × dois grupos (F2-LD-F1) | `np_longitudinal(dados, vd, id, tempo, grupos)` | nparLD + nparcomp::mctp.rm |
| várias VDs, fatores entre sujeitos | `manova_np(dados, vds, fatores)` | MANOVA.RM::MANOVA.wide |
| várias VDs com medidas repetidas | `manova_np(dados, vds, fatores, id, tempo)` | MANOVA.RM::multRM |

O banco de `np_longitudinal()` e do `manova_np()` com tempo fica no formato longo (uma
linha por pessoa e tempo), como nas aulas (`pivot_longer`). Pessoa sem todos os tempos sai
da análise e é contada em `n_excluidos`.

## 2. Efeito relativo

O efeito relativo de uma célula, p = P(X < Y) + 0,5 P(X = Y), compara a distribuição da
célula com a média das distribuições do delineamento. 0,5 é ausência de efeito; acima de
0,5, valores maiores que a média do delineamento. É o tamanho de efeito do relato (RTE nas
aulas, de relative treatment effect), com IC logit no rankFD. Com um fator e grupos do
mesmo tamanho, p = (posto médio - 0,5)/N (conferido nos evals).

## 3. `np_fatorial()`

- Teste principal: estatística do tipo ANOVA (ATS), com gl de Box. A do tipo Wald (WTS)
  fica ao lado e pede amostra grande: com n pequeno infla o erro tipo I (Brunner et al.,
  2017).
- Pós-hoc pelo `mctp` (Tukey por padrão): com um fator, depois da ATS significativa; com
  interação significativa, comparações do primeiro fator dentro de cada nível do segundo,
  como nas aulas; sem interação, comparações nos efeitos principais significativos.
- Com dois grupos e um fator, a ATS do rankFD não é o quadrado do Brunner-Munzel do
  `npar.t.test` (os estimadores de variância diferem: 16,20 contra 14,78 no eval); para dois
  grupos, a rota da skill segue sendo `teste_mann_whitney()`.

## 4. `np_longitudinal()`

- Desenhos do nparLD: LD-F1 (só tempo), F1-LD-F1 (tempo × um grupo), F2-LD-F1 (tempo × dois
  grupos). Até dois fatores entre sujeitos.
- ATS com gl1 de Box e gl2 infinito como teste principal; WTS ao lado.
- RTE por célula na tabela de efeito. O nparLD 2.3 não devolve IC do RTE; os IC saem do
  `mctp.rm` por grupo (tabela "Efeitos relativos por subamostra"), que é o que as aulas
  usam para as tabelas do relato.
- Pós-hoc: interação com o tempo significativa abre as comparações entre tempos dentro de
  cada nível do grupo que interage (Tukey por padrão; `tipo_posthoc = "Sequen"` compara
  tempos vizinhos, como em parte das aulas). Diferença entre grupos num mesmo tempo:
  `teste_mann_whitney()` em cada tempo, com correção, como as aulas fazem com
  `npar.t.test`.

## 5. `manova_np()`

- MATS (estatística do tipo ANOVA modificada) com bootstrap paramétrico como teste
  principal: invariante à escala das VDs e válida com matrizes de covariância diferentes
  entre grupos (Friedrich & Pauly, 2018). WTS ao lado, com p assintótico e bootstrap.
- `iter` (padrão 10.000, o do pacote) e `semente` fixam a reamostragem.
- Testes univariados, uma VD por vez pela mesma rota, com Holm por efeito, para dizer qual
  VD contribui (passo das aulas).
- Pós-hoc multivariado (`simCI`) fica fora da função: o `simCI` do MANOVA.RM imprime e
  devolve só a tabela; rodar à parte, se preciso, com o objeto em `$modelo`.

## 6. Divergências com as aulas

- O `nparLD.R` das aulas usa `description = F` e `plot.CI = TRUE`, argumentos que o nparLD
  2.3 removeu (a chamada dá erro). A skill usa a interface nova.
- As aulas pedem `asy.method = "log.odds"` no `mctp`. Nessa opção o estimador sai numa escala
  transformada: a comparação 3 - 1 aparece como 1,01 quando a diferença dos efeitos
  relativos é 0,39, e pode passar de 1. A skill usa `asy.method = "fisher"`, que devolve a
  diferença na escala dos efeitos relativos, com IC preso a [-1, 1]. Os p mudam pouco entre
  as duas.
- O `mctp.rm` também tem padrão `log.odds`; a skill usa `fisher` pelo mesmo motivo.

## 7. Relato

"Foi feita uma análise fatorial não paramétrica por postos (estatística do tipo ANOVA,
rankFD; Brunner et al., 2017) para Y: G: ATS(1,78; 64,21) = 33,68, p < ,001; ..." e a
tabela de efeitos relativos com IC. Para o longitudinal, o delineamento (F1-LD-F1) entra
na frase. Na MANOVA, MATS com o p do bootstrap e o número de reamostras.
