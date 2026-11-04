# Instalação

## No computador (Claude Code)

```bash
git clone https://github.com/psientrementespsi-alt/claude-psicometrista.git
cd claude-psicometrista
sh instalar.sh
```

O script copia as três skills para `~/.claude/skills/` (guardando cópia das que já existirem) e instala os pacotes R. Para instalar só os pacotes: `Rscript dependencias/instalar_pacotes.R`.

## No Claude do navegador e do app (a conta, a nuvem)

1. Baixe os três arquivos de `dist/`: `psicometria-r.zip`, `higiene-dados.zip` e `relato-estatistico.zip`.
2. Na área de Skills das configurações da sua conta, envie cada .zip (um por skill) e ative as três.

## Para o Claude rodar o R no seu computador (psicometria-r)

A skill `psicometria-r` roda a análise no R do próprio computador, para que o dado de participante não saia dele. Isso pede o servidor MCP `r-btw` no Claude Desktop:

1. Instale o R e, no R, o pacote: `install.packages("btw")` (testado com a versão 1.5.0).
2. No arquivo de configuração do Claude Desktop, dentro de `mcpServers`, acrescente:

```json
"r-btw": {
  "command": "Rscript",
  "args": ["-e", "btw::btw_mcp_server(btw::btw_tools(c(\"docs\", \"sessioninfo\", \"cran\", \"run\")))"]
}
```

3. Reinicie o Claude Desktop. Se o `Rscript` não for encontrado, troque `"command"` pelo caminho completo (no Mac com Homebrew, `/opt/homebrew/bin/Rscript`).

Sem o `r-btw`, as skills ainda orientam e escrevem o código, mas quem roda o código é você, no seu R.

## Conferir que funcionou

```bash
sh testes/rodar_todos.sh
```

Alguns valores padrão dentro das skills podem ser trocados por variáveis de ambiente `PR_*` (por exemplo `PR_DASS_CSV`). Os testes que dependem de arquivos ou pacotes que você não tem aparecem como PULADO, e não como falha.

## Tudo dentro do repositório (opcional)

**Os 22 repositórios open source de referência** (datawizard, janitor, mice, naniar, careless, statcheck, papaja, gtsummary e outros) estão em `repositorios/` como submódulos, cada um preso no commit exato que a skill usou. Para baixar tudo de uma vez:

```bash
git clone --recurse-submodules https://github.com/psientrementespsi-alt/claude-psicometrista.git
```

Se você já clonou sem isso: `git submodule update --init --recursive`.

**O código-fonte dos 376 pacotes R** (os 81 usados pelas skills mais as dependências, incluindo `gtheory` e `dmacs`, que saíram do CRAN) está na [Release `pacotes-r-2026-11-04`](https://github.com/psientrementespsi-alt/claude-psicometrista/releases/tag/pacotes-r-2026-11-04), em um arquivo de 355 MB. Para instalar o R a partir dele, sem depender do CRAN:

```bash
sh instalar_offline.sh
```

O instalador comum (`sh instalar.sh`) usa o CRAN e é bem mais rápido, porque baixa pacotes já compilados. O offline compila tudo, o que pode levar de 30 a 60 minutos.
