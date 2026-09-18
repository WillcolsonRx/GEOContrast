`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) y else x
}

trim_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

first_nonempty <- function(x, default = "") {
  x <- trim_chr(x)
  x <- x[nzchar(x)]
  if (length(x)) x[[1]] else default
}

pretty_label <- function(x) {
  special <- c(
    organism_ch1 = "Organism",
    source_name_ch1 = "Source name",
    geo_accession = "Accession",
    molecule_ch1 = "Molecule",
    library_strategy = "Library strategy"
  )
  if (x %in% names(special)) return(unname(special[[x]]))
  x <- gsub("_ch[0-9]+(\\.[0-9]+)?$", "", x, ignore.case = TRUE)
  x <- gsub("[_\\.]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  if (!nzchar(x)) return("Metadata")
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

extract_pdata <- function(obj) {
  if (inherits(obj, "ExpressionSet")) {
    return(as.data.frame(Biobase::pData(obj), stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (inherits(obj, "SummarizedExperiment")) {
    return(as.data.frame(SummarizedExperiment::colData(obj), stringsAsFactors = FALSE, check.names = FALSE))
  }
  stop("Unsupported GEO object type: ", paste(class(obj), collapse = ", "))
}

extract_series_title <- function(obj) {
  if (inherits(obj, "ExpressionSet")) {
    md <- tryCatch(GEOquery::Meta(obj), error = function(e) list())
    ttl <- first_nonempty(md$title %||% "", "")
    if (nzchar(ttl)) return(ttl)
  }
  pd <- tryCatch(extract_pdata(obj), error = function(e) data.frame())
  if ("series_title" %in% names(pd)) return(first_nonempty(pd$series_title, "GEO Series"))
  "GEO Series"
}

extract_platform <- function(obj) {
  if (inherits(obj, "ExpressionSet")) {
    return(first_nonempty(tryCatch(Biobase::annotation(obj), error = function(e) ""), "Unknown"))
  }
  pd <- tryCatch(extract_pdata(obj), error = function(e) data.frame())
  if ("platform_id" %in% names(pd)) return(first_nonempty(pd$platform_id, "Unknown"))
  "Unknown"
}

parse_characteristics <- function(pdata) {
  n <- nrow(pdata)
  if (!n) return(data.frame())
  char_cols <- grep("^characteristics", names(pdata), ignore.case = TRUE, value = TRUE)
  if (!length(char_cols)) return(data.frame(row.names = seq_len(n)))

  values <- list()
  order_keys <- character()
  put_value <- function(key, i, val) {
    key <- pretty_label(trimws(key))
    val <- trimws(val)
    if (!nzchar(key) || !nzchar(val)) return(invisible(NULL))
    if (is.null(values[[key]])) {
      values[[key]] <<- rep("", n)
      order_keys <<- c(order_keys, key)
    }
    old <- values[[key]][i]
    values[[key]][i] <<- if (!nzchar(old)) val else if (identical(old, val)) old else paste(old, val, sep = " | ")
    invisible(NULL)
  }

  for (col in char_cols) {
    vals <- trim_chr(pdata[[col]])
    for (i in seq_len(n)) {
      x <- vals[[i]]
      if (!nzchar(x)) next
      pos <- regexpr(":", x, fixed = TRUE)[1]
      if (pos > 0) {
        put_value(substr(x, 1, pos - 1), i, substr(x, pos + 1, nchar(x)))
      } else {
        put_value(col, i, x)
      }
    }
  }
  if (!length(values)) return(data.frame(row.names = seq_len(n)))
  out <- as.data.frame(values[order_keys], stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- make.unique(names(out), sep = " ")
  out
}

build_sample_metadata_one <- function(obj) {
  pdata <- extract_pdata(obj)
  n <- nrow(pdata)
  if (!n) return(data.frame())
  accession <- if ("geo_accession" %in% names(pdata)) trim_chr(pdata$geo_accession) else rownames(pdata)
  accession[!nzchar(accession)] <- rownames(pdata)[!nzchar(accession)]

  base <- data.frame(
    Accession = accession,
    Title = if ("title" %in% names(pdata)) trim_chr(pdata$title) else rep("", n),
    `Source name` = if ("source_name_ch1" %in% names(pdata)) trim_chr(pdata$source_name_ch1) else rep("", n),
    Organism = if ("organism_ch1" %in% names(pdata)) trim_chr(pdata$organism_ch1) else rep("", n),
    Platform = rep(extract_platform(obj), n),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  direct <- c("molecule_ch1", "library_strategy", "library_source", "library_selection", "instrument_model")
  for (nm in direct) {
    if (nm %in% names(pdata)) {
      lab <- pretty_label(nm)
      if (!lab %in% names(base)) base[[lab]] <- trim_chr(pdata[[nm]])
    }
  }
  chars <- parse_characteristics(pdata)
  if (ncol(chars)) {
    dup <- intersect(names(chars), names(base))
    if (length(dup)) names(chars)[match(dup, names(chars))] <- paste0(dup, " (characteristic)")
    base <- cbind(base, chars)
  }
  keep <- vapply(base, function(x) any(nzchar(trim_chr(x))), logical(1))
  keep[names(base) %in% c("Accession", "Title", "Source name", "Organism", "Platform")] <- TRUE
  base[, keep, drop = FALSE]
}

rbind_fill <- function(dfs) {
  dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, dfs)
  if (!length(dfs)) return(data.frame())
  all_names <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) d[[nm]] <- ""
    d[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}

build_series_metadata <- function(gse_list) {
  out <- rbind_fill(lapply(gse_list, build_sample_metadata_one))
  if (!nrow(out)) stop("No sample metadata were found for this GSE.")
  if ("Accession" %in% names(out)) out <- out[!duplicated(out$Accession), , drop = FALSE]
  rownames(out) <- NULL
  out
}

series_summary <- function(gse_id, gse_list, metadata) {
  data.frame(
    Field = c("GSE accession", "Series title", "Organism(s)", "Platform(s)", "Samples"),
    Value = c(
      gse_id,
      extract_series_title(gse_list[[1]]),
      paste(unique(trim_chr(metadata$Organism)[nzchar(trim_chr(metadata$Organism))]), collapse = "; "),
      paste(unique(trim_chr(metadata$Platform)[nzchar(trim_chr(metadata$Platform))]), collapse = "; "),
      nrow(metadata)
    ),
    stringsAsFactors = FALSE
  )
}

list_supplementary_files <- function(gse_id) {
  x <- tryCatch(
    GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE, fetch_files = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || !nrow(x)) return(data.frame())
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"fname" %in% names(x)) x$fname <- basename(x$url %||% rownames(x))
  if (!"url" %in% names(x)) {
    stop("GEOquery returned supplementary filenames without URLs; update GEOquery.")
  }
  x$Score <- vapply(x$fname, score_count_filename, numeric(1))
  x$Assessment <- ifelse(x$Score >= 5, "Strong count candidate",
                         ifelse(x$Score >= 2, "Possible count candidate", "Unlikely raw-count matrix"))
  x <- x[order(-x$Score, x$fname), c("fname", "Assessment", "Score", "url"), drop = FALSE]
  rownames(x) <- NULL
  x
}

score_count_filename <- function(fname) {
  s <- tolower(fname)
  score <- 0
  if (grepl("raw[_ .-]*counts?|read[_ .-]*counts?|featurecounts?|count[_ .-]*matrix|counts?[_ .-]*matrix", s)) score <- score + 7
  else if (grepl("counts?", s)) score <- score + 4
  if (grepl("matrix", s)) score <- score + 1
  if (grepl("normalized|norm[_ .-]|tpm|fpkm|rpkm|cpm|vst|rlog|log2", s)) score <- score - 7
  if (grepl("fastq|\\.fq|bam|sam|bed|bigwig|bw$|raw\\.tar|_raw\\.tar", s)) score <- score - 8
  if (grepl("\\.(csv|tsv|txt|tab|xlsx|xls)(\\.gz)?$", s)) score <- score + 2
  score
}

has_ncbi_rnaseq_counts <- function(gse_id) {
  if (!exists("hasRNASeqQuantifications", envir = asNamespace("GEOquery"), inherits = FALSE)) return(FALSE)
  tryCatch(isTRUE(GEOquery::hasRNASeqQuantifications(gse_id)), error = function(e) FALSE)
}

load_ncbi_rnaseq_counts <- function(gse_id) {
  if (!exists("getRNASeqData", envir = asNamespace("GEOquery"), inherits = FALSE)) {
    stop("This GEOquery version does not provide getRNASeqData(). Please update GEOquery.")
  }
  GEOquery::getRNASeqData(gse_id)
}

safe_groups <- function(x) {
  x <- trim_chr(x)
  unique(x[nzchar(x) & x != "-"])
}

# Robust fallback for NCBI-computed human RNA-seq counts.
# GEOquery normally discovers these links from the GEO download page, but the
# discovery step can occasionally return a false negative (for example on
# mixed-platform / mixed-organism Series or when NCBI serves a partial page).
# NCBI documents a stable direct-download filename pattern for human counts.
ncbi_human_raw_counts_url <- function(gse_id) {
  sprintf(
    "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s_raw_counts_GRCh38.p13_NCBI.tsv.gz",
    gse_id, gse_id
  )
}

ncbi_human_annotation_url <- function() {
  "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&format=file&file=Human.GRCh38.p13.annot.tsv.gz"
}

.download_ncbi_table <- function(url, destfile, label) {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({
    utils::download.file(url, destfile = destfile, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(destfile) || is.na(file.info(destfile)$size) || file.info(destfile)$size < 20) {
    stop("NCBI did not return a usable ", label, " file.")
  }

  # A failed NCBI request can be returned as a small HTML page. Detect that
  # before readr attempts to parse it as a count matrix.
  con <- file(destfile, "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = min(256L, file.info(destfile)$size))
  if (length(bytes) >= 2L && identical(as.integer(bytes[1:2]), c(31L, 139L))) {
    return(destfile) # gzip payload
  }
  txt <- tryCatch(tolower(rawToChar(bytes)), error = function(e) "")
  if (grepl("<html|<!doctype|access denied|not found|error", txt)) {
    stop("NCBI returned an HTML/error page instead of the ", label, " file.")
  }
  destfile
}

load_ncbi_direct_human_counts <- function(gse_id, metadata, organism = NULL, platform = NULL, cache_dir = tempdir()) {
  if (!identical(trimws(organism %||% ""), "Homo sapiens")) {
    stop("The direct NCBI fallback is currently defined for Homo sapiens only.")
  }

  raw_path <- file.path(cache_dir, paste0(gse_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
  ann_path <- file.path(cache_dir, "Human.GRCh38.p13.annot.tsv.gz")

  if (!file.exists(raw_path) || file.info(raw_path)$size < 20) {
    .download_ncbi_table(ncbi_human_raw_counts_url(gse_id), raw_path, "Series RNA-seq raw counts")
  }
  raw <- as.data.frame(
    readr::read_tsv(raw_path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (nrow(raw) < 2 || ncol(raw) < 3) stop("The NCBI raw-count file was empty or malformed.")

  gene_col <- if ("GeneID" %in% names(raw)) "GeneID" else names(raw)[1]
  gene_ids <- trim_chr(raw[[gene_col]])

  md <- metadata
  if (!is.null(organism) && nzchar(organism) && "Organism" %in% names(md)) {
    md <- md[trim_chr(md$Organism) == organism, , drop = FALSE]
  }
  platform <- trim_chr(platform %||% character())
  platform <- platform[nzchar(platform)]
  if (length(platform) && "Platform" %in% names(md)) {
    md <- md[trim_chr(md$Platform) %in% platform, , drop = FALSE]
  }
  sample_cols <- intersect(setdiff(names(raw), gene_col), md$Accession)
  if (length(sample_cols) < 2) {
    stop("NCBI counts were found, but fewer than two count columns match the selected organism/platform samples.")
  }

  m <- suppressWarnings(as.matrix(data.frame(lapply(raw[, sample_cols, drop = FALSE], as.numeric), check.names = FALSE)))
  rownames(m) <- gene_ids
  colnames(m) <- sample_cols
  keep <- nzchar(gene_ids) & apply(m, 1, function(z) all(is.finite(z)))
  m <- m[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]
  if (!is_count_matrix(m)) stop("The NCBI direct-download matrix was not non-negative integer-like raw counts.")
  storage.mode(m) <- "integer"

  ann <- data.frame(FeatureID = gene_ids, GeneID = gene_ids, stringsAsFactors = FALSE, check.names = FALSE)
  ann_loaded <- tryCatch(
    load_ncbi_human_annotation(gene_ids, cache_dir = cache_dir),
    error = function(e) NULL
  )
  ann_ok <- !is.null(ann_loaded)
  if (ann_ok) ann <- ann_loaded

  rownames(ann) <- gene_ids
  list(
    counts = m,
    annotation = ann,
    mapping = data.frame(CountColumn = sample_cols, GSM = sample_cols, stringsAsFactors = FALSE),
    source = if (ann_ok) "NCBI-computed RNA-seq raw counts (direct GEO download)" else "NCBI-computed RNA-seq raw counts (direct GEO download; annotation unavailable)"
  )
}

has_normalized_only_candidates <- function(supplementary) {
  if (is.null(supplementary) || !nrow(supplementary) || !"fname" %in% names(supplementary)) return(FALSE)
  any(grepl("fpkm|tpm|rpkm|cpm|normalized|norm[_ .-]|vst|rlog", tolower(supplementary$fname)))
}
