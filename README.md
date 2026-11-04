# Claude psicometrista

Três skills do Claude para análise psicométrica e estatística em psicologia, feitas para rodar no R do próprio computador:

| Skill | Versão | O que faz |
|---|---|---|
| `psicometria-r` | 2.2.0 | Análise no R: descritivas, testes, ANOVAs, regressões, TRI, AFE/AFC, redes, invariância, mediação, poder, Rasch, erro de medida |
| `higiene-dados` | 1.5.0 | Limpeza e triagem antes da análise: ausentes, careless, ajuste de pessoa (G^p, lz^p), codebook |
| `relato-estatistico` | 1.5.0 | Frases e tabelas em APA 7 a partir dos resultados, com conferência de tamanho de efeito e de p |

> Versão para compartilhar. Foram retiradas as referências ao curso original, os caminhos do computador do autor e os testes que dependem de bases de terceiros. As menções a outras skills do autor (fabrica-de-artigo, verificacao-artigo e semelhantes) podem ser ignoradas.

## Como usar

Veja `INSTALAR.md` para o passo a passo completo (computador e conta).

1. Copie as pastas de `skills/` para `~/.claude/skills/` (Claude Code) ou carregue-as como skills na sua conta.
2. Instale os pacotes R: `Rscript dependencias/instalar_pacotes.R`.
3. Opcional, para reproduzir tudo: `sh dependencias/baixar_repositorios.sh repositorios` (clona os 22 repositórios open source de referência, cada um no commit exato).
4. Rode os testes: `sh testes/rodar_todos.sh`. A última execução está em `testes/ultima_execucao.md`.

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
