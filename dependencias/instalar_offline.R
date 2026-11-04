# Instala os pacotes a partir do repositorio local de fontes (variavel PR_FONTE).
fonte <- Sys.getenv("PR_FONTE"); stopifnot(nzchar(fonte), dir.exists(fonte))
repo_local <- paste0("file://", normalizePath(fonte))
tab <- read.csv(file.path("dependencias", "pacotes_r.csv"))
pacotes <- unique(c(tab$pacote[tab$instalado], "dmacs", "gtheory"))
install.packages(pacotes, repos = c(LOCAL = repo_local, CRAN = "https://cloud.r-project.org"), type = "source")
cat("Pacotes carregaveis:", sum(vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)), "de", length(pacotes), "\n")
