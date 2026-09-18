# RNA-seq annotation recovery helpers.
# These functions never invent biological annotations. They only fill missing
# fields from GEO/NCBI annotation sources or optional organism annotation DBs.

annotation_aliases <- list(
  GeneID = c("GeneID", "gene_id", "geneid", "gene", "ENTREZ_GENE_ID", "ENTREZID",
             "entrezgene_id", "EntrezGeneID", "Gene ID", "NCBI Gene ID"),
  Symbol = c("Symbol", "symbol", "gene_symbol", "GENE_SYMBOL", "GeneSymbol", "Gene Symbol",
             "gene name", "gene_name", "external_gene_name", "HGNC symbol", "MGI symbol"),
  Description = c("Description", "description", "gene_description", "GENE_TITLE", "Gene title",
                  "gene_title", "title", "GENENAME", "gene name description"),
  Synonyms = c("Synonyms", "synonyms", "gene_synonym", "gene_synonyms", "Alias", "Aliases"),
  GeneType = c("GeneType", "gene_type", "gene_biotype", "gene_biotype_name", "biotype", "Gene type", "type"),
  EnsemblGeneID = c("EnsemblGeneID", "ensembl_gene_id", "ENSEMBL", "ensembl", "Ensembl Gene ID"),
  Status = c("Status", "status", "gene_status"),
  ChrAcc = c("ChrAcc", "chr_acc", "chromosome_accession", "chromosome", "seqname", "chr"),
  ChrStart = c("ChrStart", "chr_start", "start", "gene_start", "start_position"),
  ChrStop = c("ChrStop", "chr_stop", "stop", "end", "gene_end", "end_position"),
  Orientation = c("Orientation", "orientation", "strand"),
  Length = c("Length", "length", "gene_length"),
  GOFunctionID = c("GOFunctionID", "go_function_id", "GO Function ID"),
  GOProcessID = c("GOProcessID", "go_process_id", "GO Process ID"),
  GOComponentID = c("GOComponentID", "go_component_id", "GO Component ID"),
  GOFunction = c("GOFunction", "go_function", "GO Function"),
  GOProcess = c("GOProcess", "go_process", "GO Process"),
  GOComponent = c("GOComponent", "go_component", "GO Component")
)

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else NULL
}

annotation_feature_ids <- function(annotation) {
  if (is.null(annotation) || !nrow(annotation)) return(character())
  ann <- as.data.frame(annotation, stringsAsFactors = FALSE, check.names = FALSE)
  col <- first_existing_col(ann, c("FeatureID", annotation_aliases$GeneID, annotation_aliases$EnsemblGeneID))
  if (!is.null(col)) return(trim_chr(ann[[col]]))
  trim_chr(rownames(ann))
}

standardize_annotation_keys <- function(annotation, feature_ids = NULL) {
  ann <- if (is.null(annotation)) data.frame() else as.data.frame(annotation, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(ann)) {
    ids <- trim_chr(feature_ids %||% character())
    return(data.frame(FeatureID = ids, stringsAsFactors = FALSE, check.names = FALSE))
  }

  rn <- trim_chr(rownames(ann))
  if (!"FeatureID" %in% names(ann)) {
    candidate <- first_existing_col(ann, c(annotation_aliases$GeneID, annotation_aliases$EnsemblGeneID))
    ann$FeatureID <- if (!is.null(candidate)) trim_chr(ann[[candidate]]) else rn
  } else {
    ann$FeatureID <- trim_chr(ann$FeatureID)
    bad <- !nzchar(ann$FeatureID)
    if (length(rn) == nrow(ann)) ann$FeatureID[bad] <- rn[bad]
  }

  ann <- ann[nzchar(ann$FeatureID), , drop = FALSE]
  ann <- ann[!duplicated(ann$FeatureID), , drop = FALSE]
  rownames(ann) <- ann$FeatureID
  ann
}

annotation_value <- function(annotation, canonical) {
  if (is.null(annotation) || !nrow(annotation)) return(character())
  ann <- as.data.frame(annotation, stringsAsFactors = FALSE, check.names = FALSE)
  hit <- first_existing_col(ann, annotation_aliases[[canonical]] %||% canonical)
  if (is.null(hit)) return(rep("", nrow(ann)))
  x <- ann[[hit]]
  if (is.list(x) && !is.data.frame(x)) {
    x <- vapply(x, function(z) paste(as.character(z), collapse = "; "), character(1))
  }
  trim_chr(x)
}

annotation_coverage <- function(annotation, feature_ids = NULL) {
  ann <- standardize_annotation_keys(annotation, feature_ids)
  n <- nrow(ann)
  if (!n) {
    return(data.frame(Field = c("GeneID", "Symbol", "Description", "GeneType"),
                      NonEmpty = 0L, Total = 0L, Coverage = 0, stringsAsFactors = FALSE))
  }
  fields <- c("GeneID", "Symbol", "Description", "GeneType")
  vals <- lapply(fields, function(f) annotation_value(ann, f))
  data.frame(
    Field = fields,
    NonEmpty = vapply(vals, function(x) sum(nzchar(x)), integer(1)),
    Total = n,
    Coverage = vapply(vals, function(x) mean(nzchar(x)), numeric(1)),
    stringsAsFactors = FALSE
  )
}

annotation_has_useful_symbols <- function(annotation, feature_ids = NULL, threshold = 0.50) {
  cv <- annotation_coverage(annotation, feature_ids)
  x <- cv$Coverage[cv$Field == "Symbol"]
  length(x) == 1 && is.finite(x) && x >= threshold
}

canonicalize_annotation_columns <- function(annotation, feature_ids = NULL) {
  ann <- standardize_annotation_keys(annotation, feature_ids)
  n <- nrow(ann)
  if (!n) return(ann)
  for (canonical in names(annotation_aliases)) {
    ann[[canonical]] <- annotation_value(ann, canonical)
  }
  ann
}

fill_annotation_blanks <- function(primary, fallback, feature_ids) {
  feature_ids <- trim_chr(feature_ids)
  p <- canonicalize_annotation_columns(primary, feature_ids)
  f <- canonicalize_annotation_columns(fallback, feature_ids)

  # Align both tables to the active count-matrix feature order. Matching falls
  # back to GeneID/Ensembl IDs when FeatureID differs between sources.
  align_ann <- function(ann) {
    if (!nrow(ann)) return(data.frame(FeatureID = feature_ids, stringsAsFactors = FALSE))
    idx <- match(feature_ids, trim_chr(ann$FeatureID))
    miss <- is.na(idx)
    if (any(miss) && "GeneID" %in% names(ann)) {
      idx2 <- match(feature_ids[miss], trim_chr(ann$GeneID))
      idx[miss] <- idx2
    }
    miss <- is.na(idx)
    if (any(miss) && "EnsemblGeneID" %in% names(ann)) {
      aens <- sub("\\.[0-9]+$", "", trim_chr(ann$EnsemblGeneID))
      fens <- sub("\\.[0-9]+$", "", feature_ids[miss])
      idx2 <- match(fens, aens)
      idx[miss] <- idx2
    }
    out <- ann[idx, , drop = FALSE]
    out$FeatureID <- feature_ids
    rownames(out) <- feature_ids
    out
  }

  p <- align_ann(p)
  f <- align_ann(f)
  all_cols <- unique(c(names(p), names(f)))
  for (nm in setdiff(all_cols, names(p))) p[[nm]] <- ""
  for (nm in setdiff(all_cols, names(f))) f[[nm]] <- ""
  p <- p[, all_cols, drop = FALSE]
  f <- f[, all_cols, drop = FALSE]

  for (nm in setdiff(all_cols, "FeatureID")) {
    px <- trim_chr(p[[nm]])
    fx <- trim_chr(f[[nm]])
    take <- !nzchar(px) & nzchar(fx)
    if (any(take)) px[take] <- fx[take]
    p[[nm]] <- px
  }
  p$FeatureID <- feature_ids
  rownames(p) <- feature_ids
  p
}

read_ncbi_annotation_cached <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 100) {
    stop("Annotation cache file is absent or too small.")
  }
  a <- as.data.frame(
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (nrow(a) < 10000 || !"GeneID" %in% names(a)) stop("Annotation cache is malformed or incomplete.")
  if (!any(c("Symbol", "Description", "GeneType") %in% names(a))) {
    stop("Annotation cache lacks expected GEO annotation columns.")
  }
  if ("Symbol" %in% names(a) && mean(nzchar(trim_chr(a$Symbol))) < 0.50) {
    stop("Annotation cache contains too few gene symbols and is likely incomplete.")
  }
  a$GeneID <- trim_chr(a$GeneID)
  a <- a[nzchar(a$GeneID) & !duplicated(a$GeneID), , drop = FALSE]
  a
}

load_ncbi_human_annotation <- function(feature_ids, cache_dir = tempdir(), force_refresh = FALSE) {
  feature_ids <- trim_chr(feature_ids)
  path <- file.path(cache_dir, "Human.GRCh38.p13.annot.tsv.gz")

  read_ok <- function() tryCatch(read_ncbi_annotation_cached(path), error = function(e) NULL)
  a <- if (!isTRUE(force_refresh)) read_ok() else NULL
  if (is.null(a)) {
    if (file.exists(path)) unlink(path, force = TRUE)
    .download_ncbi_table(ncbi_human_annotation_url(), path, "Human gene annotation")
    a <- read_ok()
  }
  if (is.null(a)) {
    # One clean retry handles truncated/HTML cache files left behind by an
    # interrupted NCBI download.
    if (file.exists(path)) unlink(path, force = TRUE)
    .download_ncbi_table(ncbi_human_annotation_url(), path, "Human gene annotation")
    a <- read_ok()
  }
  if (is.null(a)) stop("The NCBI human annotation table could not be parsed after a clean retry.")

  idx <- match(feature_ids, trim_chr(a$GeneID))
  out <- a[idx, , drop = FALSE]
  out$FeatureID <- feature_ids
  # Put FeatureID/GeneID first while retaining every NCBI annotation column.
  keep <- c("FeatureID", "GeneID", setdiff(names(out), c("FeatureID", "GeneID")))
  out <- out[, unique(keep), drop = FALSE]
  rownames(out) <- feature_ids
  out
}

try_geoquery_rnaseq_annotation <- function(gse_id, feature_ids) {
  feature_ids <- trim_chr(feature_ids)
  if (!exists("getRNASeqData", envir = asNamespace("GEOquery"), inherits = FALSE)) return(NULL)
  tryCatch({
    se <- GEOquery::getRNASeqData(gse_id)
    a <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(a)) return(NULL)
    if (!"FeatureID" %in% names(a)) a$FeatureID <- rownames(a)
    fill_annotation_blanks(data.frame(FeatureID = feature_ids, stringsAsFactors = FALSE), a, feature_ids)
  }, error = function(e) NULL)
}

orgdb_package_for_organism <- function(organism) {
  org <- first_nonempty(organism, "")
  switch(org,
         "Homo sapiens" = "org.Hs.eg.db",
         "Mus musculus" = "org.Mm.eg.db",
         "Rattus norvegicus" = "org.Rn.eg.db",
         NULL)
}

try_orgdb_annotation <- function(feature_ids, organism) {
  pkg <- orgdb_package_for_organism(organism)
  if (is.null(pkg) || !requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace(pkg, quietly = TRUE)) return(NULL)

  ids0 <- trim_chr(feature_ids)
  ids_nover <- sub("\\.[0-9]+$", "", ids0)
  numeric_frac <- mean(grepl("^[0-9]+$", ids0))
  ensembl_frac <- mean(grepl("^ENS[A-Z]*G[0-9]+(?:\\.[0-9]+)?$", ids0, ignore.case = TRUE))
  keytype <- if (numeric_frac >= 0.70) "ENTREZID" else if (ensembl_frac >= 0.70) "ENSEMBL" else "SYMBOL"
  keys <- if (identical(keytype, "ENSEMBL")) ids_nover else ids0

  db <- getExportedValue(pkg, pkg)
  available_cols <- AnnotationDbi::columns(db)
  cols <- intersect(c("ENTREZID", "SYMBOL", "GENENAME", "ENSEMBL"), available_cols)
  if (!length(cols) || !keytype %in% AnnotationDbi::keytypes(db)) return(NULL)

  tab <- tryCatch(
    suppressMessages(AnnotationDbi::select(db, keys = unique(keys[nzchar(keys)]), columns = cols, keytype = keytype)),
    error = function(e) NULL
  )
  if (is.null(tab) || !nrow(tab) || !keytype %in% names(tab)) return(NULL)

  # Keep the first useful mapping per key; OrgDb can be one-to-many.
  tab[[keytype]] <- trim_chr(tab[[keytype]])
  score <- rowSums(vapply(intersect(c("SYMBOL", "GENENAME", "ENTREZID", "ENSEMBL"), names(tab)),
                         function(nm) nzchar(trim_chr(tab[[nm]])), logical(nrow(tab))))
  tab <- tab[order(tab[[keytype]], -score), , drop = FALSE]
  tab <- tab[!duplicated(tab[[keytype]]), , drop = FALSE]
  idx <- match(keys, tab[[keytype]])

  out <- data.frame(
    FeatureID = ids0,
    GeneID = if ("ENTREZID" %in% names(tab)) trim_chr(tab$ENTREZID[idx]) else "",
    Symbol = if ("SYMBOL" %in% names(tab)) trim_chr(tab$SYMBOL[idx]) else "",
    Description = if ("GENENAME" %in% names(tab)) trim_chr(tab$GENENAME[idx]) else "",
    EnsemblGeneID = if ("ENSEMBL" %in% names(tab)) trim_chr(tab$ENSEMBL[idx]) else "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- ids0
  out
}

repair_rnaseq_annotation <- function(annotation, feature_ids, gse_id = NULL, organism = NULL, cache_dir = tempdir()) {
  feature_ids <- trim_chr(feature_ids)
  current <- canonicalize_annotation_columns(annotation, feature_ids)
  sources <- character()

  before <- annotation_coverage(current, feature_ids)
  symbol_before <- before$Coverage[before$Field == "Symbol"] %||% 0

  # Human NCBI counts use stable NCBI Gene IDs and a public companion annotation
  # table. Prefer this recovery path because it matches GEO2R's own RNA-seq IDs.
  if (identical(trim_chr(organism), "Homo sapiens") && symbol_before < 0.95) {
    h <- tryCatch(load_ncbi_human_annotation(feature_ids, cache_dir = cache_dir), error = function(e) NULL)
    if (!is.null(h)) {
      current <- fill_annotation_blanks(current, h, feature_ids)
      sources <- c(sources, "NCBI Human.GRCh38.p13 annotation")
    }
  }

  # If the GSE has GEOquery/NCBI quantifications, rowData is another authoritative
  # source and works for human/mouse when available.
  cv_mid <- annotation_coverage(current, feature_ids)
  symbol_mid <- cv_mid$Coverage[cv_mid$Field == "Symbol"] %||% 0
  if (symbol_mid < 0.95 && !is.null(gse_id) && nzchar(trim_chr(gse_id))) {
    g <- try_geoquery_rnaseq_annotation(gse_id, feature_ids)
    if (!is.null(g)) {
      current <- fill_annotation_blanks(current, g, feature_ids)
      sources <- c(sources, "GEOquery RNA-seq rowData")
    }
  }

  # Optional local fallback for common organisms. This is deliberately optional
  # so the app remains multi-species without forcing several large OrgDb installs.
  cv_mid2 <- annotation_coverage(current, feature_ids)
  symbol_mid2 <- cv_mid2$Coverage[cv_mid2$Field == "Symbol"] %||% 0
  if (symbol_mid2 < 0.95) {
    o <- tryCatch(try_orgdb_annotation(feature_ids, organism), error = function(e) NULL)
    if (!is.null(o)) {
      current <- fill_annotation_blanks(current, o, feature_ids)
      sources <- c(sources, paste0(orgdb_package_for_organism(organism), " fallback"))
    }
  }

  after <- annotation_coverage(current, feature_ids)
  list(
    annotation = current,
    coverage = after,
    recovered_from = unique(sources),
    improved = sum(after$NonEmpty) > sum(before$NonEmpty)
  )
}

format_annotation_status <- function(status) {
  if (is.null(status) || is.null(status$coverage) || !nrow(status$coverage)) return("Annotation status unavailable")
  cv <- status$coverage
  getpct <- function(field) {
    x <- cv$Coverage[cv$Field == field]
    if (!length(x) || !is.finite(x[[1]])) 0 else round(100 * x[[1]], 1)
  }
  src <- status$recovered_from %||% character()
  paste0(
    "Annotation coverage — Symbol: ", getpct("Symbol"), "% · Description: ", getpct("Description"),
    "% · Gene type: ", getpct("GeneType"), "%",
    if (length(src)) paste0(" · recovered from ", paste(src, collapse = " + ")) else ""
  )
}
