# Microarray / Series Matrix helpers for GEO2R-style analysis.

selected_geo_object <- function(gse_list, platform = NULL) {
  if (is.null(gse_list) || !length(gse_list)) stop("No GEO Series object is loaded.")
  if (!is.null(platform) && nzchar(platform)) {
    hit <- which(vapply(gse_list, function(x) identical(trim_chr(extract_platform(x)), trim_chr(platform)), logical(1)))
    if (length(hit)) return(gse_list[[hit[[1]]]])
  }
  gse_list[[1]]
}

extract_experiment_type <- function(obj) {
  vals <- character()
  if (inherits(obj, "ExpressionSet")) {
    md <- tryCatch(GEOquery::Meta(obj), error = function(e) list())
    vals <- c(vals, trim_chr(md$type %||% ""), trim_chr(md$relation %||% ""))
  }
  pd <- tryCatch(extract_pdata(obj), error = function(e) data.frame())
  for (nm in intersect(c("type", "experiment_type", "library_strategy"), names(pd))) {
    vals <- c(vals, trim_chr(pd[[nm]]))
  }
  paste(unique(vals[nzchar(vals)]), collapse = "; ")
}

detect_assay_type <- function(obj, metadata = NULL) {
  txt <- tolower(extract_experiment_type(obj))
  if (!is.null(metadata) && nrow(metadata)) {
    for (nm in intersect(c("Library strategy", "Library source", "Molecule"), names(metadata))) {
      txt <- paste(txt, tolower(paste(trim_chr(metadata[[nm]]), collapse = " ")))
    }
  }
  if (grepl("rna[- ]?seq|high throughput sequencing|transcriptomic|expression profiling by high throughput sequencing", txt)) {
    return("RNA-seq")
  }
  if (grepl("array|microarray|expression profiling by array|genome tiling array", txt)) {
    return("Microarray")
  }
  # A populated Series Matrix on a GPL platform is normally the processed-data
  # route used by GEO2R for array studies. Do this only after RNA-seq metadata
  # has been checked, because RNA-seq Series can also expose processed matrices.
  if (inherits(obj, "ExpressionSet")) {
    ex <- tryCatch(Biobase::exprs(obj), error = function(e) matrix(numeric(), 0, 0))
    if (nrow(ex) > 1 && ncol(ex) > 1) return("Microarray")
  }
  "RNA-seq"
}

extract_microarray_expression <- function(obj, metadata = NULL) {
  if (!inherits(obj, "ExpressionSet")) {
    stop("The selected platform does not expose a GEO Series Matrix ExpressionSet for limma analysis.")
  }
  ex <- Biobase::exprs(obj)
  if (!is.matrix(ex) || nrow(ex) < 2 || ncol(ex) < 2) {
    stop("No usable processed expression matrix was found for the selected microarray platform.")
  }
  storage.mode(ex) <- "numeric"
  if (!is.null(metadata) && nrow(metadata)) {
    keep <- intersect(colnames(ex), metadata$Accession)
    if (length(keep) >= 2) ex <- ex[, keep, drop = FALSE]
  }
  ex
}

geo2r_log2_needed <- function(ex) {
  qx <- suppressWarnings(as.numeric(stats::quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)))
  if (length(qx) != 6 || any(!is.finite(qx))) return(FALSE)
  (qx[[5]] > 100) || ((qx[[6]] - qx[[1]] > 50) && (qx[[2]] > 0))
}

apply_geo2r_log_transform <- function(ex, mode = "Auto-detect") {
  mode <- mode %||% "Auto-detect"
  do_log <- identical(mode, "Yes") || (identical(mode, "Auto-detect") && geo2r_log2_needed(ex))
  if (do_log) {
    ex[ex <= 0] <- NA_real_
    ex <- log2(ex)
  }
  list(expression = ex, applied = do_log)
}

microarray_annotation_direct_value <- function(annotation, candidates) {
  if (is.null(annotation) || !nrow(annotation)) return(character())
  hit <- candidates[candidates %in% names(annotation)]
  if (!length(hit)) return(rep("", nrow(annotation)))
  trim_chr(annotation[[hit[[1]]]])
}

microarray_annotation_gene_ids <- function(annotation) {
  if (is.null(annotation) || !nrow(annotation)) return(rep("", 0L))

  # Prefer explicit Entrez/NCBI Gene fields. Some Brainarray alternative CDF
  # platforms (for example GPL17027) document SPOT_ID as "Gene ID".
  candidates <- c(
    "GeneID", "Gene ID", "GENE_ID", "ENTREZ_GENE_ID", "ENTREZID",
    "EntrezGeneID", "NCBI Gene ID", "SPOT_ID"
  )
  for (nm in candidates[candidates %in% names(annotation)]) {
    x <- trim_chr(annotation[[nm]])
    nonempty <- nzchar(x)
    if (!any(nonempty)) next

    # SPOT_ID is not universally an Entrez ID. Only accept it as a gene-ID
    # candidate when it is overwhelmingly numeric, which is the pattern used
    # by the Brainarray Entrez-gene CDF platforms in GEO.
    if (identical(nm, "SPOT_ID")) {
      numeric_fraction <- mean(grepl("^[0-9]+$", x[nonempty]))
      if (!is.finite(numeric_fraction) || numeric_fraction < 0.80) next
    }
    return(x)
  }
  rep("", nrow(annotation))
}

normalize_microarray_annotation <- function(annotation, organism = NULL, cache_dir = tempdir()) {
  ann <- as.data.frame(annotation, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(ann)) return(ann)

  gene_id <- microarray_annotation_gene_ids(ann)
  symbol <- microarray_annotation_direct_value(
    ann,
    c("Symbol", "SYMBOL", "Gene Symbol", "GENE_SYMBOL", "gene_symbol", "GeneSymbol")
  )
  description <- microarray_annotation_direct_value(
    ann,
    c("Description", "DESCRIPTION", "Gene title", "GENE_TITLE", "GENENAME", "gene_title")
  )

  if (!length(gene_id)) gene_id <- rep("", nrow(ann))
  if (!length(symbol)) symbol <- rep("", nrow(ann))
  if (!length(description)) description <- rep("", nrow(ann))

  # Recover missing human symbols from the same NCBI gene annotation resource
  # used by the RNA-seq path. This is valid for microarray alternative-CDF
  # platforms when their platform table supplies Entrez Gene IDs but not symbols.
  # Existing GPL annotations always win; only blank fields are filled.
  informative_ids <- nzchar(gene_id)
  symbol_coverage <- if (length(symbol)) mean(nzchar(symbol)) else 0
  recovered_from <- character()

  if (any(informative_ids) && symbol_coverage < 0.95 && identical(trim_chr(organism), "Homo sapiens")) {
    fb <- tryCatch(load_ncbi_human_annotation(gene_id, cache_dir = cache_dir), error = function(e) NULL)
    if (!is.null(fb) && nrow(fb) == nrow(ann)) {
      fs <- annotation_value(fb, "Symbol")
      fd <- annotation_value(fb, "Description")
      take <- !nzchar(symbol) & nzchar(fs)
      if (any(take)) symbol[take] <- fs[take]
      take <- !nzchar(description) & nzchar(fd)
      if (any(take)) description[take] <- fd[take]
      recovered_from <- c(recovered_from, "NCBI Gene annotation")
    }
  }

  # Optional local fallback for mouse/rat/human or any supported OrgDb. This is
  # only attempted when symbols are still sparse and never overwrites GEO data.
  symbol_coverage <- if (length(symbol)) mean(nzchar(symbol)) else 0
  if (any(informative_ids) && symbol_coverage < 0.95) {
    fb <- tryCatch(try_orgdb_annotation(gene_id, organism), error = function(e) NULL)
    if (!is.null(fb) && nrow(fb) == nrow(ann)) {
      fs <- annotation_value(fb, "Symbol")
      fd <- annotation_value(fb, "Description")
      take <- !nzchar(symbol) & nzchar(fs)
      if (any(take)) symbol[take] <- fs[take]
      take <- !nzchar(description) & nzchar(fd)
      if (any(take)) description[take] <- fd[take]
      recovered_from <- c(recovered_from, paste0(orgdb_package_for_organism(organism) %||% "OrgDb", " fallback"))
    }
  }

  # Canonical columns make the result schema stable across GPLs. Keep every
  # original platform field afterwards (SPOT_ID, DESCRIPTION, ENTREZ_GENE_ID,
  # etc.) so no GEO information is lost.
  ann$GeneID <- gene_id
  ann$Symbol <- symbol
  ann$Description <- description

  front <- c("ID", "GeneID", "Symbol", "Description")
  ann <- ann[, c(intersect(front, names(ann)), setdiff(names(ann), front)), drop = FALSE]
  attr(ann, "annotation_recovered_from") <- unique(recovered_from)
  attr(ann, "symbol_coverage") <- if (nrow(ann)) mean(nzchar(trim_chr(ann$Symbol))) else 0
  ann
}

microarray_feature_annotation <- function(obj, platform = NULL, cache_dir = tempdir(), organism = NULL) {
  ids <- rownames(Biobase::exprs(obj))
  ann <- tryCatch(as.data.frame(Biobase::fData(obj), stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (nrow(ann)) {
    if (!"ID" %in% names(ann)) ann$ID <- rownames(ann)
    ann$ID <- trim_chr(ann$ID)
  }

  # Series Matrix featureData is sometimes sparse. Pull the GPL table when the
  # platform annotation is too thin OR when it lacks both a symbol and a usable
  # Gene-ID field. The latter matters for alternative CDF platforms.
  has_symbol <- nrow(ann) && any(c("Symbol", "SYMBOL", "Gene Symbol", "GENE_SYMBOL") %in% names(ann))
  gid0 <- if (nrow(ann)) microarray_annotation_gene_ids(ann) else character()
  has_gene_id <- length(gid0) && any(nzchar(gid0))
  sparse <- !nrow(ann) || ncol(ann) <= 2 || (!has_symbol && !has_gene_id)

  if (sparse && !is.null(platform) && grepl("^GPL[0-9]+$", platform)) {
    gpl_tab <- tryCatch({
      gpl <- GEOquery::getGEO(platform, destdir = cache_dir, AnnotGPL = TRUE)
      as.data.frame(GEOquery::Table(gpl), stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) data.frame())
    if (nrow(gpl_tab)) {
      id_col <- intersect(c("ID", "ID_REF", "Probe Set ID"), names(gpl_tab))
      if (!length(id_col)) id_col <- names(gpl_tab)[1]
      gpl_tab$ID <- trim_chr(gpl_tab[[id_col[[1]]]])
      ann <- gpl_tab
    }
  }

  if (!nrow(ann)) return(data.frame(ID = ids, GeneID = "", Symbol = "", Description = "", stringsAsFactors = FALSE))
  ann <- ann[!duplicated(ann$ID), , drop = FALSE]
  ann <- ann[match(ids, ann$ID), , drop = FALSE]
  ann$ID <- ids
  rownames(ann) <- ids
  normalize_microarray_annotation(ann, organism = organism, cache_dir = cache_dir)
}

prepare_limma_inputs <- function(obj, sample_meta, group, groups_to_use,
                                 log_transform = "Auto-detect", force_normalization = FALSE,
                                 platform = NULL) {
  ex <- extract_microarray_expression(obj, sample_meta)
  sm <- sample_meta[match(colnames(ex), sample_meta$Accession), , drop = FALSE]
  sm$condition <- unname(group[match(sm$Accession, names(group))])
  groups_to_use <- safe_groups(groups_to_use)
  keep <- sm$condition %in% groups_to_use
  ex <- ex[, keep, drop = FALSE]
  sm <- sm[keep, , drop = FALSE]
  if (ncol(ex) < 4) stop("At least four assigned samples are required for microarray differential-expression analysis.")
  sm$condition <- factor(sm$condition, levels = groups_to_use)
  counts_per_group <- table(sm$condition)
  present <- names(counts_per_group)[counts_per_group > 0]
  if (length(present) < 2) stop("Fewer than two assigned groups have expression values on this platform.")
  if (any(counts_per_group[counts_per_group > 0] < 2)) {
    bad <- names(counts_per_group)[counts_per_group > 0 & counts_per_group < 2]
    stop("Each analyzed group needs at least two samples. Group(s) with <2: ", paste(bad, collapse = ", "))
  }

  tr <- apply_geo2r_log_transform(ex, log_transform)
  ex <- tr$expression
  # Remove features with no finite observations after transformation.
  keep_feature <- rowSums(is.finite(ex)) >= 2
  ex <- ex[keep_feature, , drop = FALSE]
  if (nrow(ex) < 2) stop("Too few usable features remain after expression preprocessing.")
  if (isTRUE(force_normalization)) ex <- limma::normalizeBetweenArrays(ex)

  rownames(sm) <- sm$Accession
  list(expression = ex, sample_meta = sm, groups = groups_to_use,
       log_transform_applied = tr$applied, normalization_applied = isTRUE(force_normalization),
       keep_feature = keep_feature)
}

safe_group_names <- function(groups) {
  stats::setNames(make.names(groups, unique = TRUE), groups)
}

limma_design <- function(sample_meta, groups) {
  sample_meta$condition <- factor(sample_meta$condition, levels = groups)
  design <- stats::model.matrix(~ 0 + condition, data = sample_meta)
  safe <- safe_group_names(groups)
  colnames(design) <- unname(safe[groups])
  list(design = design, safe = safe)
}

limma_contrast_matrix <- function(labels, groups, design) {
  safe <- safe_group_names(groups)
  exprs <- vapply(labels, function(lab) {
    cc <- parse_contrast_label(lab)
    if (is.null(cc) || !cc$comparison %in% groups || !cc$reference %in% groups) return(NA_character_)
    paste0(safe[[cc$comparison]], "-", safe[[cc$reference]])
  }, character(1))
  good <- !is.na(exprs)
  if (!any(good)) return(NULL)
  cm <- limma::makeContrasts(contrasts = exprs[good], levels = design)
  colnames(cm) <- labels[good]
  cm
}

append_microarray_annotation <- function(tt, annotation) {
  ids <- rownames(tt)
  ann <- annotation[match(ids, annotation$ID), , drop = FALSE]
  ann_cols <- setdiff(names(ann), c("ID", names(tt)))
  # GEO2R data columns first, platform annotation after them.
  data.frame(ID = ids, tt, ann[, ann_cols, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
}

standardize_limma_pairwise <- function(fit, annotation, adjust_method = "BH") {
  tt <- limma::topTable(fit, number = Inf, sort.by = "P", adjust.method = adjust_method)
  # Keep AveExpr internally for the mean-difference plot, but don't include it
  # in GEO2R's default/export data-column set.
  if (!"AveExpr" %in% names(tt)) tt$AveExpr <- NA_real_
  out <- append_microarray_annotation(tt, annotation)
  if (!"F" %in% names(out)) out$F <- NA_real_
  front <- c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "F", "AveExpr")
  for (nm in front) if (!nm %in% names(out)) out[[nm]] <- NA_real_
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]
  rownames(out) <- NULL
  out
}

standardize_limma_multigroup <- function(fit, annotation, adjust_method = "BH") {
  tt <- limma::topTable(fit, coef = seq_len(ncol(fit$coefficients)), number = Inf, sort.by = "F", adjust.method = adjust_method)
  # topTable with multiple coefficients can include the individual contrast
  # coefficients. GEO2R's documented overall table centers on P/adj.P/F, so
  # keep those and blank the two-group-only statistics.
  ids <- rownames(tt)
  p <- tt$P.Value %||% rep(NA_real_, nrow(tt))
  adj <- tt$adj.P.Val %||% stats::p.adjust(p, method = adjust_method)
  f <- tt$F %||% rep(NA_real_, nrow(tt))
  core <- data.frame(
    adj.P.Val = adj, P.Value = p, t = NA_real_, B = NA_real_, logFC = NA_real_, F = f,
    AveExpr = NA_real_, stringsAsFactors = FALSE, check.names = FALSE
  )
  rownames(core) <- ids
  out <- append_microarray_annotation(core, annotation)
  rownames(out) <- NULL
  out
}

run_geo2r_limma <- function(obj, sample_meta, group, groups_to_use,
                            platform = NULL, organism = NULL, cache_dir = tempdir(),
                            log_transform = "Auto-detect", force_normalization = FALSE,
                            vooma = FALSE, p_adjust = "BH", plot_contrasts = character()) {
  inp <- prepare_limma_inputs(
    obj = obj, sample_meta = sample_meta, group = group, groups_to_use = groups_to_use,
    log_transform = log_transform, force_normalization = force_normalization, platform = platform
  )
  des <- limma_design(inp$sample_meta, inp$groups)
  annotation <- microarray_feature_annotation(obj, platform = platform, cache_dir = cache_dir, organism = organism)
  annotation <- annotation[match(rownames(inp$expression), annotation$ID), , drop = FALSE]

  if (isTRUE(vooma)) {
    vw <- limma::vooma(inp$expression, design = des$design, plot = FALSE)
    fit0 <- limma::lmFit(vw, des$design)
  } else {
    fit0 <- limma::lmFit(inp$expression, des$design)
  }

  ng <- length(inp$groups)
  requested <- intersect(plot_contrasts, contrast_labels(inp$groups))
  if (!length(requested)) requested <- default_contrast_labels(inp$groups)
  pairwise <- list()
  for (lab in requested) {
    cm <- limma_contrast_matrix(lab, inp$groups, des$design)
    if (is.null(cm)) next
    fitp <- limma::eBayes(limma::contrasts.fit(fit0, cm))
    pairwise[[lab]] <- standardize_limma_pairwise(fitp, annotation, adjust_method = p_adjust)
  }

  if (ng == 2L) {
    main_label <- paste(inp$groups[[1]], "vs", inp$groups[[2]])
    if (!main_label %in% names(pairwise)) {
      cm <- limma_contrast_matrix(main_label, inp$groups, des$design)
      fitp <- limma::eBayes(limma::contrasts.fit(fit0, cm))
      pairwise[[main_label]] <- standardize_limma_pairwise(fitp, annotation, adjust_method = p_adjust)
    }
    main <- pairwise[[main_label]]
    test_name <- "limma moderated t-test"
  } else {
    overall_labels <- default_contrast_labels(inp$groups, max_n = max(1L, length(inp$groups) - 1L))
    cm <- limma_contrast_matrix(overall_labels, inp$groups, des$design)
    fitm <- limma::eBayes(limma::contrasts.fit(fit0, cm))
    main <- standardize_limma_multigroup(fitm, annotation, adjust_method = p_adjust)
    test_name <- "limma moderated F-test"
  }

  annotation_cols <- setdiff(names(main), c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "F", "AveExpr"))
  list(
    engine = "limma",
    assay_type = "Microarray",
    result = main,
    pairwise = pairwise,
    normalized = inp$expression,
    sample_meta = inp$sample_meta,
    groups = inp$groups,
    p_adjust = p_adjust,
    test_name = test_name,
    log_transform_applied = inp$log_transform_applied,
    normalization_applied = inp$normalization_applied,
    vooma_applied = isTRUE(vooma),
    annotation_status = list(
      symbol_coverage = attr(annotation, "symbol_coverage") %||% if ("Symbol" %in% names(annotation)) mean(nzchar(trim_chr(annotation$Symbol))) else 0,
      recovered_from = attr(annotation, "annotation_recovered_from") %||% character()
    ),
    annotation_columns = annotation_cols,
    export_columns = c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "F", annotation_cols),
    default_display_columns = c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "F", head(annotation_cols, 6))
  )
}

microarray_pca_dataframe <- function(expression, sample_meta) {
  rvv <- matrixStats::rowVars(expression, na.rm = TRUE)
  rvv[!is.finite(rvv)] <- 0
  take <- order(rvv, decreasing = TRUE)[seq_len(min(500L, length(rvv)))]
  x <- t(expression[take, , drop = FALSE])
  # Replace rare missing values by feature medians for PCA only.
  if (anyNA(x)) {
    for (j in seq_len(ncol(x))) {
      miss <- !is.finite(x[, j])
      if (any(miss)) x[miss, j] <- stats::median(x[, j], na.rm = TRUE)
    }
  }
  pc <- stats::prcomp(x, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  out <- data.frame(Accession = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2], stringsAsFactors = FALSE)
  sm <- sample_meta[match(out$Accession, sample_meta$Accession), , drop = FALSE]
  out$Group <- as.character(sm$condition)
  out$Title <- sm$Title %||% sm$Accession
  attr(out, "pct") <- pct[1:2]
  out
}

microarray_profile_dataframe <- function(expression, sample_meta, feature_id, result_table = NULL) {
  feature_id <- trimws(feature_id %||% "")
  if (!nzchar(feature_id)) stop("Enter a probe/feature ID or annotation value.")
  rid <- rownames(expression)
  match_id <- match(feature_id, rid)
  if (is.na(match_id) && !is.null(result_table) && nrow(result_table)) {
    candidates <- setdiff(names(result_table), c("adj.P.Val", "P.Value", "t", "B", "logFC", "F", "AveExpr"))
    for (nm in candidates) {
      idx <- which(tolower(trim_chr(result_table[[nm]])) == tolower(feature_id))
      if (length(idx)) {
        match_id <- match(result_table$ID[idx[[1]]], rid)
        if (!is.na(match_id)) break
      }
    }
  }
  if (is.na(match_id)) stop("Probe/feature was not found in the analyzed expression matrix.")
  vals <- expression[match_id, ]
  sm <- sample_meta[match(names(vals), sample_meta$Accession), , drop = FALSE]
  data.frame(Accession = names(vals), Title = sm$Title, Group = as.character(sm$condition), Expression = as.numeric(vals), stringsAsFactors = FALSE)
}

generate_geo2r_limma_script <- function(gse_id, platform, groups, assignments, p_adjust_label,
                                         log_transform = "Auto-detect", force_normalization = FALSE,
                                         vooma = FALSE, test_name = NULL) {
  rq <- function(x) paste0('"', gsub('"', '\\"', as.character(x), fixed = TRUE), '"')
  used <- assignments[assignments != "-"]
  c(
    "# GEO2R Multi-Species — microarray limma analysis sketch",
    sprintf("# GSE: %s", gse_id %||% ""),
    sprintf("# Platform: %s", platform %||% ""),
    sprintf("# Test: %s", test_name %||% "limma"),
    "",
    "library(GEOquery)",
    "library(limma)",
    "",
    sprintf("gse <- getGEO(%s, GSEMatrix = TRUE)", rq(gse_id %||% "GSE")),
    sprintf("# Select platform %s when the Series has multiple GPLs.", platform %||% ""),
    sprintf("groups <- c(%s)", paste(rq(groups), collapse = ", ")),
    sprintf("# Assigned samples: %s", paste(sprintf("%s=%s", names(used), unname(used)), collapse = "; ")),
    sprintf("# Log transform: %s", log_transform),
    sprintf("# Force normalizeBetweenArrays: %s", isTRUE(force_normalization)),
    sprintf("# vooma precision weights: %s", isTRUE(vooma)),
    sprintf("# P-value adjustment: %s", p_adjust_label %||% "BH"),
    "# The app constructs a no-intercept group design, applies requested contrasts,",
    "# then uses lmFit(), contrasts.fit(), eBayes(), and topTable()."
  )
}
