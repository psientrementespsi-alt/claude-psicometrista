# ==============================================================================
# docx_apa.R  (relato-estatistico v1.3.0)
# Escreve tabela APA 7 em .docx so com R base (OOXML + utils::zip), sem
# flextable/officer: tres linhas horizontais (acima e abaixo do cabecalho e no
# fim), sem verticais, Times New Roman 12, titulo em negrito, subtitulo em
# italico, nota abaixo. *texto* vira italico. Alinhamento "decimal" usa tab
# decimal dentro da celula.
# ==============================================================================

.xml_esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE); s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE); gsub("\"", "&quot;", s, fixed = TRUE)
}

.runs <- function(txt, negrito = FALSE, italico_base = FALSE, tab_inicial = FALSE) {
  txt <- gsub("\\$", "", txt)
  partes <- strsplit(txt, "\\*", perl = TRUE)[[1]]
  if (!length(partes)) partes <- ""
  out <- if (tab_inicial) "<w:r><w:tab/></w:r>" else ""
  for (i in seq_along(partes)) {
    if (partes[i] == "") next
    ital <- xor(italico_base, i %% 2 == 0)
    rpr <- paste0("<w:rPr><w:rFonts w:ascii=\"Times New Roman\" w:hAnsi=\"Times New Roman\" w:cs=\"Times New Roman\"/>",
                  if (negrito) "<w:b/>" else "", if (ital) "<w:i/>" else "", "<w:sz w:val=\"24\"/></w:rPr>")
    out <- paste0(out, "<w:r>", rpr, "<w:t xml:space=\"preserve\">", .xml_esc(partes[i]), "</w:t></w:r>")
  }
  out
}

.par <- function(txt, negrito = FALSE, italico = FALSE, alin = "left", espaco_depois = 0, tab_dec = NULL) {
  tabs <- if (!is.null(tab_dec)) sprintf("<w:tabs><w:tab w:val=\"decimal\" w:pos=\"%d\"/></w:tabs>", tab_dec) else ""
  paste0("<w:p><w:pPr>", tabs, "<w:spacing w:before=\"0\" w:after=\"", espaco_depois, "\" w:line=\"240\" w:lineRule=\"auto\"/>",
         "<w:jc w:val=\"", alin, "\"/></w:pPr>", .runs(txt, negrito, italico, tab_inicial = !is.null(tab_dec)), "</w:p>")
}

.cel <- function(txt, largura, alin, borda_sup = FALSE, borda_inf = FALSE) {
  b <- function(lado, on) if (on) sprintf("<w:%s w:val=\"single\" w:sz=\"8\" w:space=\"0\" w:color=\"000000\"/>", lado) else
    sprintf("<w:%s w:val=\"nil\"/>", lado)
  jc <- switch(alin, esq = "left", dir = "right", centro = "center", decimal = "left", "left")
  tab <- if (alin == "decimal") as.integer(largura * 0.55) else NULL
  paste0("<w:tc><w:tcPr><w:tcW w:w=\"", largura, "\" w:type=\"dxa\"/><w:tcBorders>",
         b("top", borda_sup), b("left", FALSE), b("bottom", borda_inf), b("right", FALSE),
         "</w:tcBorders></w:tcPr>", .par(txt, alin = jc, tab_dec = tab), "</w:tc>")
}

#' Escreve uma ou mais tabelas APA num .docx
#'
#' @param arquivo caminho do .docx
#' @param tabelas lista de tabelas; cada uma e list(numero, titulo, headers,
#'   linhas = list(c(...), ...), alinhamentos = c("esq","decimal",...), nota)
#' @param largura_total largura util em twips (padrao 9360 = 6,5 polegadas)
#' @param paisagem pagina deitada
escrever_tabelas_docx <- function(arquivo, tabelas, largura_total = 9360, paisagem = FALSE) {
  if (!is.null(tabelas$headers)) tabelas <- list(tabelas)
  corpo <- character(0)
  for (tb in tabelas) {
    k <- length(tb$headers)
    if (any(vapply(tb$linhas, length, 1L) != k)) stop("Linha com numero de colunas diferente do cabecalho em: ", tb$titulo)
    al <- if (is.null(tb$alinhamentos)) c("esq", rep("decimal", k - 1)) else tb$alinhamentos
    larg <- rep(as.integer(largura_total / k), k)
    if (k > 1) { larg[1] <- as.integer(largura_total * 0.34); larg[-1] <- as.integer((largura_total - larg[1]) / (k - 1)) }
    rotulo <- if (!is.null(tb$numero)) paste("Tabela", tb$numero) else "Tabela"
    corpo <- c(corpo, .par(rotulo, negrito = TRUE), .par(tb$titulo, italico = TRUE, espaco_depois = 120))
    grid <- paste0("<w:tblGrid>", paste0(sprintf("<w:gridCol w:w=\"%d\"/>", larg), collapse = ""), "</w:tblGrid>")
    tbl <- paste0("<w:tbl><w:tblPr><w:tblW w:w=\"", sum(larg), "\" w:type=\"dxa\"/><w:tblLayout w:type=\"fixed\"/>",
                  "<w:tblCellMar><w:left w:w=\"60\" w:type=\"dxa\"/><w:right w:w=\"60\" w:type=\"dxa\"/></w:tblCellMar></w:tblPr>", grid)
    al_cab <- ifelse(seq_len(k) == 1, "esq", "centro")
    tbl <- paste0(tbl, "<w:tr><w:trPr><w:tblHeader/></w:trPr>",
                  paste0(mapply(.cel, tb$headers, larg, al_cab, MoreArgs = list(borda_sup = TRUE, borda_inf = TRUE)), collapse = ""), "</w:tr>")
    nl <- length(tb$linhas)
    for (i in seq_len(nl)) {
      tbl <- paste0(tbl, "<w:tr>", paste0(mapply(.cel, tb$linhas[[i]], larg, al,
                    MoreArgs = list(borda_sup = FALSE, borda_inf = i == nl)), collapse = ""), "</w:tr>")
    }
    corpo <- c(corpo, paste0(tbl, "</w:tbl>"))
    if (!is.null(tb$nota) && nzchar(tb$nota)) corpo <- c(corpo, .par(tb$nota, espaco_depois = 240))
    corpo <- c(corpo, .par(""))
  }
  sect <- if (paisagem) "<w:sectPr><w:pgSz w:w=\"15840\" w:h=\"12240\" w:orient=\"landscape\"/><w:pgMar w:top=\"1440\" w:right=\"1440\" w:bottom=\"1440\" w:left=\"1440\" w:header=\"720\" w:footer=\"720\" w:gutter=\"0\"/></w:sectPr>" else
    "<w:sectPr><w:pgSz w:w=\"12240\" w:h=\"15840\"/><w:pgMar w:top=\"1440\" w:right=\"1440\" w:bottom=\"1440\" w:left=\"1440\" w:header=\"720\" w:footer=\"720\" w:gutter=\"0\"/></w:sectPr>"
  doc <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>",
                "<w:document xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\"><w:body>",
                paste0(corpo, collapse = ""), sect, "</w:body></w:document>")
  ct <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>",
               "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">",
               "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>",
               "<Default Extension=\"xml\" ContentType=\"application/xml\"/>",
               "<Override PartName=\"/word/document.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml\"/></Types>")
  rels <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>",
                 "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">",
                 "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"word/document.xml\"/></Relationships>")
  tmp <- tempfile("docx_apa_"); dir.create(file.path(tmp, "_rels"), recursive = TRUE); dir.create(file.path(tmp, "word"))
  esc <- function(txt, f) { con <- file(f, "wb"); writeBin(charToRaw(enc2utf8(txt)), con); close(con) }
  esc(ct, file.path(tmp, "[Content_Types].xml")); esc(rels, file.path(tmp, "_rels", ".rels")); esc(doc, file.path(tmp, "word", "document.xml"))
  arquivo <- normalizePath(arquivo, mustWork = FALSE)
  if (file.exists(arquivo)) file.remove(arquivo)
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  utils::zip(arquivo, files = c("[Content_Types].xml", "_rels/.rels", "word/document.xml"), flags = "-q -X")
  invisible(arquivo)
}

#' Converte uma tabela ja montada pelo gerador (headers + linhas_dados) para o formato de escrever_tabelas_docx
tabela_para_docx <- function(numero, titulo, headers, linhas_dados, alinhamentos = NULL, nota = NULL) {
  list(numero = numero, titulo = titulo, headers = headers, linhas = linhas_dados,
       alinhamentos = if (!is.null(alinhamentos)) sub("^dir$", "decimal", alinhamentos) else NULL, nota = nota)
}
