#!/bin/sh
# Instala as tres skills no Claude Code e os pacotes R. Uso: sh instalar.sh
set -e
DEST="${CLAUDE_SKILLS_DIR:-$HOME/.claude/skills}"
AQUI="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$DEST"
for s in psicometria-r higiene-dados relato-estatistico; do
  if [ -d "$DEST/$s" ]; then mv "$DEST/$s" "$DEST/$s.bak-$(date +%Y%m%d%H%M%S)"; fi
  cp -R "$AQUI/skills/$s" "$DEST/$s"; echo "skill instalada: $s"
done
if command -v Rscript >/dev/null 2>&1; then Rscript "$AQUI/dependencias/instalar_pacotes.R"; else echo "Rscript nao encontrado: instale o R (https://cloud.r-project.org) e rode: Rscript dependencias/instalar_pacotes.R"; fi
echo "Pronto. Abra o Claude Code e peca, por exemplo: analise a confiabilidade da minha escala no R."
