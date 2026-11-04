# Script reproduzível

Movido do SKILL.md na v1.4 (o SKILL.md guarda o resumo e as regras). Conteúdo
testado na v1.2; a v1.4 só trocou o número da versão no cabeçalho do script.


Por que existe. Em 62 Registered Reports de psicologia, 36 compartilharam dados
e código, os scripts rodaram em 31 e os resultados principais se reproduziram
em 21; entre os que não reproduziram, em 8 faltava código para parte dos
valores, como uma macro própria do SPSS (Obels et al., 2020, AMPPS,
https://doi.org/10.1177/2515245920918872). Nos 14 artigos de uma edição da
Psychological Science em que todos ganharam o selo de dados abertos, 6 traziam
código, 1 se reproduziu exatamente e 3 com desvios pequenos (Crüwell et al.,
2023, Psychological Science, https://doi.org/10.1177/09567976221140828). Na
Cognition, especificação pouco clara da análise e erro de relato estavam entre
o que impediu a reprodução (Hardwicke et al., 2018, Royal Society Open Science,
https://doi.org/10.1098/rsos.180448). Aqui o risco é próprio da rota 1, porque
o `btw_tool_run_r` roda o código no ambiente global e a documentação dele não
fala em registro do que foi executado. Sem esta seção, o código fica só na
conversa. As regras abaixo aplicam parte das dez regras de Sandve et al.
(2013, PLoS Computational Biology, https://doi.org/10.1371/journal.pcbi.1003285)
e as recomendações de Obels et al. (2020).

O que esta seção garante e o que não garante. A reexecução no mesmo Mac mostra
que o script está completo e não depende do estado da sessão. Ela não mostra
que o script roda em outro computador ou com versões futuras dos pacotes; isso
pede o `renv` e, no limite, um contêiner (ver "Opcional").

## Abertura

Uma chamada, no início da sessão, depois do `carregar.R`. Trocar
`pasta` pelo caminho da pasta dos dados. O r-btw devolve o diretório de
trabalho para `/` a cada chamada (e também opções e variáveis de ambiente), por
isso o `registrar()` fixa a pasta dentro dele; objetos do ambiente global
continuam de uma chamada para a outra. A abertura copia as funções da skill
para dentro da pasta, porque o script depende delas como aquele artigo da
amostra de Obels dependia da macro do SPSS. Sem a cópia, quem recebe o script
não tem como rodá-lo, e uma atualização futura da skill mudaria a análise
antiga.

```r
pasta  <- "~/caminho/da/pasta/dos/dados"
script <- file.path(pasta, sprintf("analise_%s.R", format(Sys.Date())))
copia  <- file.path(pasta, "psicometria-r")
if (!dir.exists(copia)) {
  dir.create(copia)
  invisible(file.copy(list.files("~/.local/share/psicometria-r",
                                 pattern = "^(carregar|psicometria|est_).*\\.R$",
                                 full.names = TRUE), copia))
}
if (!file.exists(script)) writeLines(c(
  sprintf("# %s | psicometria-r 1.4 | criado em %s", basename(script),
          format(Sys.time(), "%d/%m/%Y %H:%M")),
  "# Reproduzir: abrir o R com esta pasta como diretório de trabalho e rodar source() deste arquivo.",
  "# Funções da skill: pasta psicometria-r/ ao lado. Pacotes e versões: bloco 'sessão' no fim.",
  'source("psicometria-r/carregar.R")',
  "set.seed(2026)"), script)
setwd(pasta); source(script)
registrar <- function(bloco, rotulo) {
  wd <- setwd(pasta); on.exit(setwd(wd))  # o r-btw volta o diretório para "/" a cada chamada
  source(exprs = parse(text = bloco), local = globalenv(),
         echo = FALSE, print.eval = TRUE)
  cat(sprintf("\n# ---- %s | %s\n", rotulo, format(Sys.time(), "%H:%M")),
      bloco, "\n", file = script, append = TRUE, sep = "")
  invisible(TRUE)
}
```

## Leitura do banco

O primeiro bloco registrado. O md5 do arquivo entra no
próprio bloco, como o WORCS faz com o checksum dos dados, e a reexecução para
se o banco tiver mudado.

```r
registrar(sprintf(r"---(
stopifnot(unname(tools::md5sum("banco.csv")) == "%s")
dados <- read.csv("banco.csv")
)---", tools::md5sum(file.path(pasta, "banco.csv"))),
"[descritivo, relatado] Passo 1, leitura do banco")
```

## Cada análise

Entra como bloco de texto cru, que roda e só depois é anexado.
Bloco que dá erro para no `source()` e não chega ao arquivo.

```r
registrar(r"---(
res_t <- teste_t_independente(dados, "escore", "grupo", decimal = ",")
relatar(res_t)
)---", "[confirmatório, relatado] Passo 2, H1 escore por grupo")
```

Regras do registro:

- O rótulo diz o tipo e o destino: `[descritivo, ...]`, `[confirmatório, ...]`
  ou `[exploratório, ...]`, seguido de `relatado` ou `não relatado`.
  Confirmatório é só o que estava previsto no Passo 0 ou no pré-registro. Obels
  et al. (2020) acharam descritiva, confirmatória e exploratória difíceis de
  separar nos scripts que avaliaram.
- O arquivo bruto não é editado à mão depois da abertura. Limpeza,
  recodificação, inversão de itens e exclusão são blocos registrados (Sandve et
  al., 2013, regra 2).
- Caminho sempre relativo à pasta, nunca absoluto.
- Bloco com passo aleatório fora das funções com argumento `semente` começa
  com `set.seed()`. Sem isso, um bloco a mais ou a menos muda a sequência
  aleatória e a reexecução diverge.
- Figura do relato sai de bloco registrado que grava o arquivo (`png()` ou
  `ggsave()`), e os números por trás dela ficam num objeto do script (Sandve et
  al., 2013, regra 7).
- Nenhum ID de participante escrito no código. Exclusão sai por critério sobre
  variável (`subset(dados, tempo_resposta >= 300)`); se a lista de IDs for
  inevitável, ela fica num arquivo local lido pelo script, que nunca vai para o
  OSF.
- Ajuda, `checar_pacotes()` e consulta de documentação não passam por
  `registrar()`.
- Bloco corrigido entra de novo como bloco novo; o anterior, se rodou sem erro,
  continua no arquivo e roda de novo na reexecução.
- Segunda sessão no mesmo dia e na mesma pasta continua o mesmo arquivo, e o
  `source(script)` da abertura refaz os blocos anteriores antes de seguir.
- Reanálise pedida por parecerista é sessão nova na mesma pasta, com o mesmo
  fechamento, rodado de novo antes da resubmissão (Obels et al., 2020,
  recomendam a checagem final depois da revisão por pares).

## Fechamento (antes do Passo 3)

Os três passos são obrigatórios.

1. Registrar o bloco `numeros`, com os valores que vão para o pacote de
   resultados. É ele que liga cada frase do relato ao número que o script
   produz (Sandve et al., 2013, regra 9).

```r
registrar(r"---(
numeros <- list(t = res_t$tabela, efeito = res_t$efeito)
)---", "[confirmatório, relatado] números do relato")
```

2. Anexar a sessão e reexecutar o script num R limpo, sem o `.Rprofile`.

```r
cat("\n# ---- sessão\n", paste("#", capture.output(sessionInfo())),
    file = script, append = TRUE, sep = "\n")
rout <- sub("\\.R$", ".Rout", script)
cmd  <- sprintf('setwd("%s"); source("%s"); saveRDS(numeros, "numeros_rerun.rds")',
                pasta, basename(script))
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c("--no-init-file", "-e", shQuote(cmd)), stdout = rout, stderr = rout)
list(status = status,
     igual = all.equal(numeros, readRDS(file.path(pasta, "numeros_rerun.rds"))))
```

3. Ler o resultado. `status` 0 e `igual` TRUE liberam a entrega. Qualquer
   outra coisa é investigada antes do Passo 3 (bloco faltando, semente fora,
   caminho absoluto, banco alterado, objeto criado fora do `registrar()`), e o
   pacote só sai depois de a reexecução bater ou de a divergência ser explicada
   no LOG. O `.Rout` fica na pasta e não vem para a conversa; dele se lê só a
   última linha de erro, se houver.

O template passou por testes em 25/10/2026 no container (R 4.3.3, dados
simulados). Com `t.test()` e um bootstrap, o bloco com erro ficou fora do
arquivo e o bootstrap reproduziu com o `set.seed()` do cabeçalho. Com os blocos
desta seção, a cópia das funções e o `teste_t_independente()` da skill, a
reexecução deu `status` 0 e `all.equal` TRUE, e uma segunda sessão sobre o
mesmo arquivo também. Com o banco alterado depois da leitura, a reexecução
parou no `stopifnot()` do md5. No Mac, pelo r-btw (25/10/2026, R 4.6.1, macOS
26.5.2, dados simulados em `tempdir()`, cada bloco numa chamada separada), o
primeiro teste achou a falha do diretório: o bloco de leitura parou no md5
porque a chamada começou em `/`. Com o `setwd()` dentro do `registrar()`, a
leitura, a análise, o bloco aleatório, o bloco com erro (que ficou fora do
arquivo), o fechamento, a trava do md5 e a segunda sessão passaram, com
`status` 0 e `all.equal` TRUE.

## Opcional, para o mestrado e o OSF

 `renv::init()` e `renv::snapshot()` na
pasta do projeto travam as versões dos pacotes num lockfile, recurso que o
WORCS adota (Van Lissa et al., 2021, Data Science,
https://doi.org/10.3233/DS-210031); Git na pasta guarda o histórico do script
(Sandve et al., 2013, regra 4). O fluxo com R Markdown, Git, Make e Docker faz o
relato sair direto do código e cobre a reprodução em outra máquina e no longo
prazo (Peikert & Brandmaier, 2021, Quantitative and Computational Methods in
Behavioral Sciences, https://doi.org/10.5964/qcmb.3763), e o Quarto cumpre o
mesmo papel do R Markdown. Entram só se o projeto pedir, não por padrão.

O depósito no OSF leva o script, a pasta `psicometria-r/`, o LOG, um README
que diz qual bloco gera qual resultado do artigo e um codebook com nome,
rótulo e valores de cada variável, como Obels et al. (2020) recomendam. Os
dados só vão com o consentimento permitindo. Sem os dados reais, o pacote worcs
gera codebook (`make_codebook()`) e dados sintéticos (`synthetic()`,
`closed_data()`), e quem recebe o depósito roda o script no banco sintético,
trocando o nome do arquivo e o md5 no bloco de leitura.
