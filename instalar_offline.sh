#!/bin/sh
# Instala os pacotes R a partir do codigo-fonte guardado na Release do repositorio (sem depender do CRAN).
# Util se o CRAN estiver fora do ar ou se um pacote for arquivado (gtheory e dmacs ja foram).
# Uso: sh instalar_offline.sh   (baixa cerca de 355 MB; a compilacao pode levar de 30 a 60 minutos)
set -e
AQUI="$(cd "$(dirname "$0")" && pwd)"
DEST="$AQUI/.cran-fonte"
mkdir -p "$DEST"
URL="https://github.com/psientrementespsi-alt/claude-psicometrista/releases/download/pacotes-r-2026-11-04/cran-fonte.tar"
if [ ! -d "$DEST/cran-fonte" ]; then curl -L --fail -o "$DEST/cran-fonte.tar" "$URL"; (cd "$DEST" && tar -xf cran-fonte.tar && rm cran-fonte.tar); fi
echo "SHA-256 do arquivo original: 9f2bad31c79aa88ade2e5e3b9d50c6051f06c8971f77557a0dace2aa5926bd0d"
cd "$AQUI" && PR_FONTE="$DEST/cran-fonte/src/contrib" Rscript dependencias/instalar_offline.R
