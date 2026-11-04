#!/bin/sh
# Clona os 22 repositorios open source usados como referencia pelas skills, cada um no commit exato.
set -e
DEST="${1:-repositorios}"
mkdir -p "$DEST"
[ -d "$DEST/dlookr/.git" ] || git clone --quiet https://github.com/choonghyunryu/dlookr.git "$DEST/dlookr"; git -C "$DEST/dlookr" checkout --quiet 759b7b01d9cb4cd5a38ff5868f975c386ca1b1af
[ -d "$DEST/codebook/.git" ] || git clone --quiet https://github.com/rubenarslan/codebook.git "$DEST/codebook"; git -C "$DEST/codebook" checkout --quiet 6433edb323f50b7d1737b5b91277f1f129162eb2
[ -d "$DEST/datawizard/.git" ] || git clone --quiet https://github.com/easystats/datawizard.git "$DEST/datawizard"; git -C "$DEST/datawizard" checkout --quiet 44ac26d5e67a55fd0652ffbcdf40c95e36345d21
[ -d "$DEST/skimr/.git" ] || git clone --quiet https://github.com/ropensci/skimr.git "$DEST/skimr"; git -C "$DEST/skimr" checkout --quiet 004d4068e4097dbff6c5dad1b8e01c1d36bd7337
[ -d "$DEST/janitor/.git" ] || git clone --quiet https://github.com/sfirke/janitor.git "$DEST/janitor"; git -C "$DEST/janitor" checkout --quiet 81702b6ed2b97a143319700a8edf48e8e4cce9cd
[ -d "$DEST/simputation/.git" ] || git clone --quiet https://github.com/markvanderloo/simputation.git "$DEST/simputation"; git -C "$DEST/simputation" checkout --quiet 38965ad0857716182050a446d8416c9e47280569
[ -d "$DEST/naniar/.git" ] || git clone --quiet https://github.com/njtierney/naniar.git "$DEST/naniar"; git -C "$DEST/naniar" checkout --quiet 89f2a5df3743df02b6e601d403fa001fc8a79c79
[ -d "$DEST/mice/.git" ] || git clone --quiet https://github.com/amices/mice.git "$DEST/mice"; git -C "$DEST/mice" checkout --quiet 61083e667fa41cc4bc655492591b100bf8a7e443
[ -d "$DEST/careless/.git" ] || git clone --quiet https://github.com/ryentes/careless.git "$DEST/careless"; git -C "$DEST/careless" checkout --quiet 0a84caf1a1b7cb7319fa9fd40254de525456e09c
[ -d "$DEST/tidylog/.git" ] || git clone --quiet https://github.com/elbersb/tidylog.git "$DEST/tidylog"; git -C "$DEST/tidylog" checkout --quiet 65a398acbd42e444b86fd67e7bb9625e2b3103f9
[ -d "$DEST/lumberjack/.git" ] || git clone --quiet https://github.com/markvanderloo/lumberjack.git "$DEST/lumberjack"; git -C "$DEST/lumberjack" checkout --quiet e83dcb1e1a765bc9b391298f93c20ded4cd4a136
[ -d "$DEST/validate/.git" ] || git clone --quiet https://github.com/data-cleaning/validate.git "$DEST/validate"; git -C "$DEST/validate" checkout --quiet cb1e94528d3ede03776cb102bdf7017f821886c7
[ -d "$DEST/assertr/.git" ] || git clone --quiet https://github.com/tonyfischetti/assertr.git "$DEST/assertr"; git -C "$DEST/assertr" checkout --quiet e484665dba582582faa343d97a09d4a097db3bb1
[ -d "$DEST/pointblank/.git" ] || git clone --quiet https://github.com/rstudio/pointblank.git "$DEST/pointblank"; git -C "$DEST/pointblank" checkout --quiet 0dc5ae0e5e00ca579fa8436598f2cd75d8218344
[ -d "$DEST/papaja/.git" ] || git clone --quiet https://github.com/crsh/papaja.git "$DEST/papaja"; git -C "$DEST/papaja" checkout --quiet 6161f020c5d79159ea1602c9b5c06b8fceb062c5
[ -d "$DEST/modelsummary/.git" ] || git clone --quiet https://github.com/vincentarelbundock/modelsummary.git "$DEST/modelsummary"; git -C "$DEST/modelsummary" checkout --quiet 847bece31af20d6f4b1bccd44321aec94e7b6da4
[ -d "$DEST/apaTables/.git" ] || git clone --quiet https://github.com/dstanley4/apaTables.git "$DEST/apaTables"; git -C "$DEST/apaTables" checkout --quiet 12f802a82820f16108a75467a24d7a072d10f07b
[ -d "$DEST/equatiomatic/.git" ] || git clone --quiet https://github.com/datalorax/equatiomatic.git "$DEST/equatiomatic"; git -C "$DEST/equatiomatic" checkout --quiet 225bdf5982d67e745874a96792727c4084798d19
[ -d "$DEST/statcheck/.git" ] || git clone --quiet https://github.com/MicheleNuijten/statcheck.git "$DEST/statcheck"; git -C "$DEST/statcheck" checkout --quiet 1088bb199331145c4eaa4cbdb14e79ba02bafa92
[ -d "$DEST/report/.git" ] || git clone --quiet https://github.com/easystats/report.git "$DEST/report"; git -C "$DEST/report" checkout --quiet 03b348070e0682a1b1aef35594a3f2ba15b00dd3
[ -d "$DEST/gtsummary/.git" ] || git clone --quiet https://github.com/ddsjoberg/gtsummary.git "$DEST/gtsummary"; git -C "$DEST/gtsummary" checkout --quiet 1f5cc1108528fd279b7cd08a9b6436e7f853fdd4
[ -d "$DEST/rempsyc/.git" ] || git clone --quiet https://github.com/rempsyc/rempsyc.git "$DEST/rempsyc"; git -C "$DEST/rempsyc" checkout --quiet 44ae7336dd5e4d1dc5635abf5825d8a7aaae410f
echo "ok: repositorios em $DEST"
