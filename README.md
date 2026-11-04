# Claude psicometrista

Três skills do Claude para análise psicométrica e estatística em psicologia, feitas para rodar no R do próprio computador:

| Skill | Versão | O que faz |
|---|---|---|
| `psicometria-r` | 2.2.0 | Análise no R: descritivas, testes, ANOVAs, regressões, TRI, AFE/AFC, redes, invariância, mediação, poder, Rasch, erro de medida |
| `higiene-dados` | 1.5.0 | Limpeza e triagem antes da análise: ausentes, careless, ajuste de pessoa (G^p, lz^p), codebook |
| `relato-estatistico` | 1.5.0 | Frases e tabelas em APA 7 a partir dos resultados, com conferência de tamanho de efeito e de p |

> Versão para compartilhar. Foram retiradas as referências ao curso original, os caminhos do computador do autor e os testes que dependem de bases de terceiros. As menções a outras skills do autor (fabrica-de-artigo, verificacao-artigo e semelhantes) podem ser ignoradas.

## Instalação rápida

**No computador (Claude Code).** Precisa do [R](https://cloud.r-project.org) instalado. No Terminal:

```bash
git clone https://github.com/psientrementespsi-alt/claude-psicometrista.git
cd claude-psicometrista
sh instalar.sh
```

O `instalar.sh` copia as três skills para `~/.claude/skills/` e instala os pacotes R.

**Na conta Claude (navegador e app).** Baixe os três arquivos da pasta [`dist/`](dist) (`psicometria-r.zip`, `higiene-dados.zip`, `relato-estatistico.zip`) e envie cada um na área de Skills das configurações.

**Para a `psicometria-r` rodar o R no seu computador**, configure o servidor `r-btw` no Claude Desktop. O passo a passo está em [`INSTALAR.md`](INSTALAR.md).

### Opcional: trazer tudo para o computador

Os 22 repositórios open source de referência (submódulos, presos nos commits exatos):

```bash
git submodule update --init --recursive
```

O código-fonte dos pacotes R, sem depender do CRAN (baixa 355 MB da [Release](https://github.com/psientrementespsi-alt/claude-psicometrista/releases/tag/pacotes-r-2026-11-04) e compila, de 30 a 60 minutos):

```bash
sh instalar_offline.sh
```

### Conferir que funcionou

```bash
sh testes/rodar_todos.sh
```

## O que há aqui

- `skills/`: as três skills completas (SKILL.md, scripts R, evals, referências, CHANGELOG).
- `dependencias/`: `pacotes_r.csv` (78 pacotes com versão), `instalar_pacotes.R`, `sessionInfo.txt`, `repositorios_clonados.csv` e `baixar_repositorios.sh`.
- `testes/`: runner, resumo e logs.

## O que não está aqui, de propósito

- Dados de participantes e qualquer base real (os testes usam só dados simulados ou públicos).
- Livros e artigos em PDF e o material de cursos: têm copyright. As skills os citam por nome e DOI.
- Os 22 repositórios clonados entram como submódulos (links presos no commit exato), não como cópia. O código-fonte dos pacotes R está na Release do repositório.

## Segurança dos dados

Dado de participante não sai do computador do usuário. As skills só devolvem agregados, suprimem células de tabela cruzada com n < 5 e não usam `head()`, `print(dados)` nem `View()` em dados reais.


## Licença

MIT (veja `LICENSE`). Os pacotes R e os repositórios listados em `dependencias/` mantêm as licenças deles.
