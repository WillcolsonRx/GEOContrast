normalize_name <- function(x) {
  x <- trim_chr(x)
  x <- tolower(x)
  x <- gsub("^x(?=gsm[0-9]+$)", "", x, perl = TRUE)
  gsub("[^a-z0-9]+", "", x)
}

is_count_like_vector <- function(x, tol = 1e-7) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(FALSE)
  all(x >= 0) && mean(abs(x - round(x)) <= tol) >= 0.995
}

is_count_matrix <- function(m, tol = 1e-7) {
  if (!is.matrix(m) || !length(m)) return(FALSE)
  if (!is.numeric(m)) return(FALSE)
  vals <- as.numeric(m)
  vals <- vals[is.finite(vals)]
  if (!length(vals) || any(vals < 0)) return(FALSE)
  mean(abs(vals - round(vals)) <= tol) >= 0.995
}

read_tabular_file <- function(path, original_name = basename(path), sheet = NULL) {
  nm <- tolower(original_name)
  if (grepl("\\.xlsx$|\\.xls$", nm)) {
    sh <- sheet %||% readxl::excel_sheets(path)[1]
    return(as.data.frame(readxl::read_excel(path, sheet = sh, .name_repair = "minimal"), check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (grepl("\\.csv(\\.gz)?$", nm)) {
    return(as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal"), check.names = FALSE, stringsAsFactors = FALSE))
  }
  as.data.frame(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal"), check.names = FALSE, stringsAsFactors = FALSE)
}

available_excel_sheets <- function(path, original_name) {
  if (!grepl("\\.xlsx$|\\.xls$", tolower(original_name))) return(character())
  tryCatch(readxl::excel_sheets(path), error = function(e) character())
}

match_count_columns <- function(df, metadata) {
  if (!nrow(metadata) || !ncol(df)) return(data.frame())
  cols <- names(df)
  ckey <- normalize_name(cols)
  acc_key <- normalize_name(metadata$Accession)
  title_key <- normalize_name(metadata$Title)

  out <- lapply(seq_along(cols), function(i) {
    idx_acc <- match(ckey[[i]], acc_key)
    method <- ""
    idx <- idx_acc
    if (!is.na(idx_acc)) {
      method <- "GSM accession"
    } else {
      idx_title <- match(ckey[[i]], title_key)
      idx <- idx_title
      if (!is.na(idx_title)) method <- "Sample title"
    }
    data.frame(
      CountColumn = cols[[i]],
      GSM = if (!is.na(idx)) metadata$Accession[[idx]] else "",
      Title = if (!is.na(idx)) metadata$Title[[idx]] else "",
      Match = method,
      Numeric = is.numeric(df[[i]]) || all(is.na(suppressWarnings(as.numeric(df[[i]]))) == is.na(df[[i]])),
      CountLike = is_count_like_vector(df[[i]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

suggest_gene_id_column <- function(df, mapping) {
  candidates <- names(df)[!names(df) %in% mapping$CountColumn[nzchar(mapping$GSM)]]
  if (!length(candidates)) return(names(df)[1])
  priority <- c("geneid", "gene_id", "gene", "ensembl_gene_id", "entrezgene", "entrez", "id", "feature")
  keys <- normalize_name(candidates)
  hit <- match(normalize_name(priority), keys, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit)) candidates[[hit[[1]]]] else candidates[[1]]
}

prepare_uploaded_counts <- function(df, metadata, gene_id_col, organism = NULL) {
  if (!gene_id_col %in% names(df)) stop("Choose a valid feature/gene ID column.")
  mapping <- match_count_columns(df, metadata)
  matched <- mapping[nzchar(mapping$GSM) & mapping$CountLike, , drop = FALSE]
  if (!is.null(organism) && nzchar(organism)) {
    allowed <- metadata$Accession[metadata$Organism == organism]
    matched <- matched[matched$GSM %in% allowed, , drop = FALSE]
  }
  matched <- matched[!duplicated(matched$GSM), , drop = FALSE]
  if (nrow(matched) < 2) {
    stop("Fewer than two count-matrix columns could be matched to GEO samples. Use a file whose sample columns are GSM accessions or GEO sample titles.")
  }

  raw <- df[, matched$CountColumn, drop = FALSE]
  m <- suppressWarnings(as.matrix(data.frame(lapply(raw, as.numeric), check.names = FALSE)))
  colnames(m) <- matched$GSM
  gene_ids <- trim_chr(df[[gene_id_col]])
  keep <- nzchar(gene_ids) & apply(m, 1, function(z) all(is.finite(z)))
  m <- m[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]
  if (!is_count_matrix(m)) stop("The selected sample columns are not non-negative integer-like raw counts. TPM/FPKM/CPM/log-normalized values should not be analyzed with DESeq2.")
  storage.mode(m) <- "integer"

  if (anyDuplicated(gene_ids)) {
    m <- rowsum(m, group = gene_ids, reorder = FALSE)
    gene_ids <- rownames(m)
  } else {
    rownames(m) <- gene_ids
  }

  annotation_cols <- setdiff(names(df), c(gene_id_col, matched$CountColumn))
  ann <- data.frame(FeatureID = gene_ids, stringsAsFactors = FALSE)
  if (length(annotation_cols)) {
    original_ann <- df[keep, c(gene_id_col, annotation_cols), drop = FALSE]
    names(original_ann)[1] <- "FeatureID"
    original_ann$FeatureID <- trim_chr(original_ann$FeatureID)
    original_ann <- original_ann[!duplicated(original_ann$FeatureID), , drop = FALSE]
    ann <- merge(ann, original_ann, by = "FeatureID", all.x = TRUE, sort = FALSE)
    ann <- ann[match(gene_ids, ann$FeatureID), , drop = FALSE]
  }
  rownames(ann) <- gene_ids

  list(counts = m, annotation = ann, mapping = matched, source = "Supplementary/uploaded raw-count matrix")
}

prepare_ncbi_counts <- function(se, organism = NULL, platforms = NULL) {
  m <- SummarizedExperiment::assay(se)
  pd <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE, check.names = FALSE)
  gsm <- if ("geo_accession" %in% names(pd)) trim_chr(pd$geo_accession) else colnames(m)
  colnames(m) <- gsm

  if (!is.null(organism) && nzchar(organism)) {
    org_col <- intersect(c("organism_ch1", "organism", "Organism"), names(pd))
    if (length(org_col)) {
      keep <- trim_chr(pd[[org_col[[1]]]]) == organism
      m <- m[, keep, drop = FALSE]
      pd <- pd[keep, , drop = FALSE]
    }
  }

  platforms <- trim_chr(platforms %||% character())
  platforms <- platforms[nzchar(platforms)]
  if (length(platforms)) {
    plat_col <- intersect(c("platform_id", "Platform", "platform", "gpl"), names(pd))
    if (length(plat_col)) {
      keep <- trim_chr(pd[[plat_col[[1]]]]) %in% platforms
      m <- m[, keep, drop = FALSE]
      pd <- pd[keep, , drop = FALSE]
    }
  }
  if (!is_count_matrix(m)) stop("NCBI quantification assay was not integer-like raw counts.")
  storage.mode(m) <- "integer"
  ann <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(ann)) ann <- data.frame(FeatureID = rownames(m), stringsAsFactors = FALSE)
  if (!"FeatureID" %in% names(ann)) ann <- data.frame(FeatureID = rownames(m), ann, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(ann) <- rownames(m)
  list(counts = m, annotation = ann, mapping = data.frame(CountColumn = colnames(m), GSM = colnames(m)), source = "NCBI-computed RNA-seq raw counts")
}

count_summary <- function(count_obj) {
  if (is.null(count_obj)) return(data.frame())
  m <- count_obj$counts
  data.frame(
    Metric = c("Count source", "Features", "Matched samples", "Minimum count", "Maximum count", "Total reads/counts"),
    Value = c(count_obj$source, nrow(m), ncol(m), min(m), max(m), format(sum(m), scientific = TRUE)),
    stringsAsFactors = FALSE
  )
}
