#!/bin/sh
# Roda os 23 conjuntos de testes das tres skills. Uso: sh testes/rodar_todos.sh
cd "$(dirname "$0")/../skills" || exit 1
for f in */scripts/evals*.R; do echo "=== $f"; Rscript "$f" 2>&1 | grep -E 'passaram|Passaram|PULADO|FALH|Falharam'; done
