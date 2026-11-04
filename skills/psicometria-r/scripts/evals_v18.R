# evals_v18.R, testes de robustez e novas funcoes da v1.8:
# 1. ajustar_afc com chamada auto-contida e avaliada
# 2. confiabilidade(ic = TRUE, metodo = "delta") instantaneo (< 0.5s) com CIs validos
# 3. confiabilidade(ic = TRUE, metodo = "bootstrap") sem falhas de escopo

args_arquivo <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pasta <- if (length(args_arquivo)) dirname(normalizePath(args_arquivo)) else "scripts"
source(file.path(pasta, "carregar.R"))
resultados <- list()
teste <- function(nome, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { message("  erro: ", conditionMessage(e)); FALSE })
  cat(sprintf("%-7s %s\n", if (ok) "PASSOU" else "FALHOU", nome))
  resultados[[nome]] <<- ok
  invisible(ok)
}
pular <- function(nome, motivo) cat(sprintf("PULADO  %s (%s)\n", nome, motivo))
quieto <- function(expr) suppressWarnings(suppressMessages(expr))
tem <- function(p) requireNamespace(p, quietly = TRUE)

cat("== v1.8.0: Confiabilidade Delta e Bootstrap Sanitizado\n")

if (tem("lavaan") && tem("semTools")) {
  set.seed(2026)
  N <- 150
  f_true <- rnorm(N)
  d <- data.frame(
    x1 = round(pmax(1, pmin(5, 3 + 0.8 * f_true + rnorm(N, 0, 0.4)))),
    x2 = round(pmax(1, pmin(5, 3 + 0.7 * f_true + rnorm(N, 0, 0.4)))),
    x3 = round(pmax(1, pmin(5, 3 + 0.9 * f_true + rnorm(N, 0, 0.4)))),
    x4 = round(pmax(1, pmin(5, 3 + 0.6 * f_true + rnorm(N, 0, 0.4))))
  )
  itens <- paste0("x", 1:4)
  
  # Teste 1: ajustar_afc cria chamada avaliada que permite update() sem erro de escopo
  teste("ajustar_afc: chamada auto-contida suporta update() fora do escopo", {
    fit <- quieto(ajustar_afc("F =~ x1 + x2 + x3 + x4", d, itens, ordinal = TRUE))
    db <- d[sample(1:N, N, replace = TRUE), ]
    up <- tryCatch(quieto(lavaan::update(fit, data = db)), error = function(e) NULL)
    !is.null(up) && isTRUE(lavaan::lavInspect(up, "converged"))
  })
  
  # Teste 2: confiabilidade(ic = TRUE, metodo = "delta") produz CIs numericos validos em < 1s
  teste("confiabilidade(metodo = 'delta'): executa em < 1s e gera CIs finitos", {
    fit <- quieto(ajustar_afc("F =~ x1 + x2 + x3 + x4", d, itens, ordinal = TRUE))
    t0 <- Sys.time()
    res <- quieto(confiabilidade(fit, d, ic = TRUE, metodo = "delta"))
    t1 <- Sys.time()
    dt <- as.numeric(difftime(t1, t0, units = "secs"))
    dt < 1.0 && 
      !is.na(res$omega_ic_inf[1]) && !is.na(res$omega_ic_sup[1]) &&
      res$omega_ic_inf[1] <= res$omega_delta[1] && res$omega_delta[1] <= res$omega_ic_sup[1] &&
      res$alfa_ic_inf[1] <= res$alfa[1] && res$alfa[1] <= res$alfa_ic_sup[1]
  })
  
  # Teste 3: confiabilidade(ic = TRUE, metodo = "bootstrap") nao falha por 'estimador' nao encontrado
  teste("confiabilidade(metodo = 'bootstrap'): sanitiza fit@call e conclui sem erro de escopo", {
    fit <- quieto(ajustar_afc("F =~ x1 + x2 + x3 + x4", d, itens, ordinal = TRUE))
    res_b <- quieto(confiabilidade(fit, d, ic = TRUE, metodo = "bootstrap", reamostras = 5))
    !is.na(res_b$omega_ic_inf[1]) && !is.na(res_b$omega_ic_sup[1])
  })
}


cat("== v1.8.1: omega delta igual ao parametro definido escrito a mao\n")
if (tem("lavaan") && tem("semTools")) {
  set.seed(11); nn <- 400; ff <- rnorm(nn); lam <- c(.8, .7, .6, .7, .5, .6)
  Xc <- as.data.frame(sapply(lam, function(l) l * ff + rnorm(nn, sd = sqrt(1 - l^2)))); names(Xc) <- paste0("x", 1:6)
  fc <- lavaan::cfa("F =~ x1+x2+x3+x4+x5+x6", data = Xc)
  rc <- quieto(confiabilidade(fc, Xc, ic = TRUE, metodo = "delta"))
  m2 <- paste("F =~ l1*x1+l2*x2+l3*x3+l4*x4+l5*x5+l6*x6", paste0("x", 1:6, "~~e", 1:6, "*x", 1:6, collapse = "\n"),
              "om := (l1+l2+l3+l4+l5+l6)^2/((l1+l2+l3+l4+l5+l6)^2 + e1+e2+e3+e4+e5+e6)", sep = "\n")
  oc <- subset(lavaan::parameterEstimates(lavaan::cfa(m2, data = Xc, std.lv = TRUE)), lhs == "om")
  teste("continuos: omega_delta = compRelSEM e IC = parametro definido (tol. 1e-3)",
        abs(rc$omega_delta - rc$omega) < 1e-3 && abs(rc$omega_ic_inf - oc$ci.lower) < 1e-3 && abs(rc$omega_ic_sup - oc$ci.upper) < 1e-3)
  set.seed(2026); N <- 300; g <- rnorm(N); lo <- c(.8, .7, .6, .75, .5)
  Xo <- as.data.frame(sapply(lo, function(l) cut(l * g + rnorm(N, sd = sqrt(1 - l^2)), c(-Inf, -.5, .3, 1.1, Inf), labels = FALSE) - 1)); names(Xo) <- paste0("y", 1:5)
  fo <- quieto(ajustar_afc("F =~ y1+y2+y3+y4+y5", Xo, paste0("y", 1:5), estimador = "WLSMV", ordinal = TRUE))
  ro <- quieto(confiabilidade(fo, Xo, ic = TRUE, metodo = "delta"))
  m3 <- paste("F =~ a1*y1+a2*y2+a3*y3+a4*y4+a5*y5", "om := (a1+a2+a3+a4+a5)^2/((a1+a2+a3+a4+a5)^2 + (1-a1^2)+(1-a2^2)+(1-a3^2)+(1-a4^2)+(1-a5^2))", sep = "\n")
  oo <- subset(lavaan::parameterEstimates(lavaan::cfa(m3, data = Xo, ordered = paste0("y", 1:5), estimator = "WLSMV", std.lv = TRUE)), lhs == "om")
  teste("ordinais: omega_delta, EP e IC iguais ao parametro definido (tol. 1e-3) e IC contem omega_delta",
        abs(ro$omega_delta - oo$est) < 1e-3 && abs(ro$omega_ep - oo$se) < 1e-3 && abs(ro$omega_ic_inf - oo$ci.lower) < 1e-3 &&
        abs(ro$omega_ic_sup - oo$ci.upper) < 1e-3 && ro$omega_ic_inf <= ro$omega_delta && ro$omega_delta <= ro$omega_ic_sup && ro$omega_ic_sup < 1)
  teste("alfa: IC de Feldt igual ao psych::alpha.ci (tol. 1e-3)", {
    fe <- psych::alpha.ci(rc$alfa, nn, n.var = 6); abs(rc$alfa_ic_inf - fe$lower) < 1e-3 && abs(rc$alfa_ic_sup - fe$upper) < 1e-3 })
}
cat(sprintf("\nTotal de testes v1.8: %d | Passaram: %d | Falharam: %d\n",
            length(resultados), sum(unlist(resultados)), sum(!unlist(resultados))))
