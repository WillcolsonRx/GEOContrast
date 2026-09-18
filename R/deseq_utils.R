make_safe_factor <- function(x, reference = NULL) {
  x <- factor(x)
  if (!is.null(reference) && reference %in% levels(x)) x <- stats::relevel(x, ref = reference)
  x
}

sanitize_design_name <- function(x) {
  make.names(x, unique = TRUE)
}

geo2r_padjust_method <- function(label) {
  map <- c(
    "Benjamini & Hochberg (False discovery rate)" = "fdr",
    "Benjamini & Yekutieli" = "BY",
    "Bonferroni" = "bonferroni",
    "Hochberg" = "hochberg",
    "Holm" = "holm",
    "Hommel" = "hommel"
  )
  out <- unname(map[[label %||% ""]])
  if (is.null(out) || !nzchar(out)) "BH" else out
}

safe_annotation_value <- function(df, candidates, n) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(rep("", n))
  x <- df[[hit[[1]]]]
  if (is.list(x) && !is.data.frame(x)) {
    x <- vapply(x, function(z) paste(as.character(z), collapse = "; "), character(1))
  }
  trim_chr(x)
}

standardize_geo2r_annotation <- function(annotation, feature_ids) {
  feature_ids <- trim_chr(feature_ids)
  n <- length(feature_ids)
  ann <- if (is.null(annotation)) data.frame() else as.data.frame(annotation, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(ann)) ann <- data.frame(FeatureID = feature_ids, stringsAsFactors = FALSE)
  if (!"FeatureID" %in% names(ann)) ann$FeatureID <- rownames(ann)
  ann$FeatureID <- trim_chr(ann$FeatureID)
  rn <- trim_chr(rownames(ann))
  bad_feature <- !nzchar(ann$FeatureID)
  if (length(rn) == nrow(ann)) ann$FeatureID[bad_feature] <- rn[bad_feature]
  ann <- ann[nzchar(ann$FeatureID) & !duplicated(ann$FeatureID), , drop = FALSE]

  # Primary match is FeatureID. If an annotation source uses GeneID or Ensembl
  # IDs as the true key, recover unmatched rows rather than returning blank
  # Symbol/Description columns.
  idx <- match(feature_ids, ann$FeatureID)
  miss <- is.na(idx)
  gene_col <- first_existing_col(ann, annotation_aliases$GeneID)
  if (any(miss) && !is.null(gene_col)) {
    idx[miss] <- match(feature_ids[miss], trim_chr(ann[[gene_col]]))
  }
  miss <- is.na(idx)
  ens_col <- first_existing_col(ann, annotation_aliases$EnsemblGeneID)
  if (any(miss) && !is.null(ens_col)) {
    q <- sub("\\.[0-9]+$", "", feature_ids[miss])
    k <- sub("\\.[0-9]+$", "", trim_chr(ann[[ens_col]]))
    idx[miss] <- match(q, k)
  }
  x <- ann[idx, , drop = FALSE]

  getv <- function(candidates) safe_annotation_value(x, candidates, n)

  geneid <- getv(c("GeneID", "gene_id", "geneid", "gene", "ENTREZ_GENE_ID", "ENTREZID",
                   "entrezgene_id", "EntrezGeneID", "Gene ID"))
  geneid[!nzchar(geneid)] <- feature_ids[!nzchar(geneid)]

  data.frame(
    GeneID = geneid,
    Symbol = getv(c("Symbol", "symbol", "gene_symbol", "GENE_SYMBOL", "gene_name", "GeneSymbol", "Gene Symbol", "external_gene_name")),
    Description = getv(c("Description", "description", "gene_description", "GENE_TITLE", "Gene title", "gene_title", "title")),
    Synonyms = getv(c("Synonyms", "synonyms", "gene_synonym", "gene_synonyms", "Alias", "Aliases")),
    GeneType = getv(c("GeneType", "gene_type", "gene_biotype", "gene_biotype_name", "biotype", "Gene type", "type")),
    EnsemblGeneID = getv(c("EnsemblGeneID", "ensembl_gene_id", "ENSEMBL", "ensembl", "Ensembl Gene ID")),
    Status = getv(c("Status", "status", "gene_status")),
    ChrAcc = getv(c("ChrAcc", "chr_acc", "chromosome_accession", "chromosome", "seqname", "chr")),
    ChrStart = getv(c("ChrStart", "chr_start", "start", "gene_start", "start_position")),
    ChrStop = getv(c("ChrStop", "chr_stop", "stop", "end", "gene_end", "end_position")),
    Orientation = getv(c("Orientation", "orientation", "strand")),
    Length = getv(c("Length", "length", "gene_length")),
    GOFunctionID = getv(c("GOFunctionID", "go_function_id", "GO Function ID")),
    GOProcessID = getv(c("GOProcessID", "go_process_id", "GO Process ID")),
    GOComponentID = getv(c("GOComponentID", "go_component_id", "GO Component ID")),
    GOFunction = getv(c("GOFunction", "go_function", "GO Function")),
    GOProcess = getv(c("GOProcess", "go_process", "GO Process")),
    GOComponent = getv(c("GOComponent", "go_component", "GO Component")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

standardize_geo2r_result <- function(res, annotation = NULL, multigroup = FALSE) {
  tab <- as.data.frame(res)
  feature_ids <- rownames(tab)
  if (is.null(feature_ids)) feature_ids <- as.character(seq_len(nrow(tab)))

  need <- c("padj", "pvalue", "lfcSE", "stat", "log2FoldChange", "baseMean")
  for (nm in need) if (!nm %in% names(tab)) tab[[nm]] <- NA_real_
  if (isTRUE(multigroup)) {
    tab$lfcSE <- NA_real_
    tab$log2FoldChange <- NA_real_
  }

  ann <- standardize_geo2r_annotation(annotation, feature_ids)
  out <- data.frame(
    padj = tab$padj,
    pvalue = tab$pvalue,
    lfcSE = tab$lfcSE,
    stat = tab$stat,
    log2FoldChange = tab$log2FoldChange,
    baseMean = tab$baseMean,
    ann,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$ID <- feature_ids
  out <- out[, c("ID", "padj", "pvalue", "lfcSE", "stat", "log2FoldChange", "baseMean",
                 "GeneID", "Symbol", "Description", "Synonyms", "GeneType", "EnsemblGeneID", "Status",
                 "ChrAcc", "ChrStart", "ChrStop", "Orientation", "Length", "GOFunctionID", "GOProcessID",
                 "GOComponentID", "GOFunction", "GOProcess", "GOComponent"), drop = FALSE]
  out <- out[order(out$padj, out$pvalue, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

parse_contrast_label <- function(label) {
  parts <- strsplit(label %||% "", " vs ", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NULL)
  list(comparison = parts[[1]], reference = parts[[2]])
}

contrast_labels <- function(groups) {
  groups <- safe_groups(groups)
  if (length(groups) < 2) return(character())
  # GEO2R group order is meaningful: the group defined first is the numerator
  # (test) and the group defined second is the denominator (control).
  out <- character()
  for (i in seq_len(length(groups) - 1L)) {
    for (j in (i + 1L):length(groups)) out <- c(out, paste(groups[[i]], "vs", groups[[j]]))
  }
  out
}

default_contrast_labels <- function(groups, max_n = 5L) {
  groups <- safe_groups(groups)
  if (length(groups) < 2) return(character())
  if (length(groups) == 2) return(paste(groups[[1]], "vs", groups[[2]]))
  # GEO2R presents adjacent groups in creation order; include the last-vs-first
  # circular contrast when room is available.
  labs <- vapply(seq_len(length(groups) - 1L), function(i) paste(groups[[i]], "vs", groups[[i + 1L]]), character(1))
  labs <- c(labs, paste(groups[[length(groups)]], "vs", groups[[1]]))
  head(labs, max_n)
}

prepare_deseq_inputs <- function(counts, sample_meta, group, groups_to_use,
                                 covariates = character(), min_count = 10L,
                                 min_samples = NULL, prefilter_mode = c("geo2r", "custom")) {
  prefilter_mode <- match.arg(prefilter_mode)
  stopifnot(all(colnames(counts) %in% sample_meta$Accession))
  groups_to_use <- safe_groups(groups_to_use)
  if (length(groups_to_use) < 2) stop("Assign matched samples to at least two groups.")

  sm <- sample_meta[match(colnames(counts), sample_meta$Accession), , drop = FALSE]
  sm$condition <- unname(group[match(sm$Accession, names(group))])
  keep_sample <- sm$condition %in% groups_to_use
  counts <- counts[, keep_sample, drop = FALSE]
  sm <- sm[keep_sample, , drop = FALSE]
  if (ncol(counts) < 4) stop("At least four matched, assigned samples are required for analysis.")

  sm$condition <- factor(sm$condition, levels = groups_to_use)
  present <- levels(droplevels(sm$condition))
  if (length(present) < 2) stop("Fewer than two assigned groups have matched count columns.")
  counts_per_group <- table(sm$condition)
  nonzero_group_sizes <- counts_per_group[counts_per_group > 0]
  if (any(nonzero_group_sizes < 2)) {
    bad <- names(nonzero_group_sizes)[nonzero_group_sizes < 2]
    stop("Each analyzed group needs at least two matched samples. Group(s) with <2: ", paste(bad, collapse = ", "))
  }

  covariates <- intersect(covariates, names(sm))
  covariates <- covariates[vapply(covariates, function(v) {
    x <- trim_chr(sm[[v]])
    length(unique(x[nzchar(x)])) > 1 && !any(!nzchar(x))
  }, logical(1))]

  for (v in covariates) {
    x <- sm[[v]]
    suppressWarnings(num <- as.numeric(x))
    if (length(num) && all(is.finite(num))) sm[[v]] <- num else sm[[v]] <- factor(trim_chr(x))
  }

  # Match the current GEO2R RNA-seq pre-filter by default:
  # keep genes with raw count >= 10 in at least N samples, where N is the
  # size of the smallest assigned group. A custom mode is retained only as an
  # explicit extended-analysis option.
  smallest_group_size <- as.integer(min(nonzero_group_sizes))
  if (identical(prefilter_mode, "geo2r")) {
    effective_min_count <- 10L
    effective_min_samples <- smallest_group_size
  } else {
    effective_min_count <- as.integer(min_count %||% 10L)
    effective_min_samples <- as.integer(min_samples %||% smallest_group_size)
  }
  if (!is.finite(effective_min_count) || effective_min_count < 0) effective_min_count <- 10L
  if (!is.finite(effective_min_samples) || effective_min_samples < 1) effective_min_samples <- smallest_group_size
  effective_min_samples <- min(effective_min_samples, ncol(counts))

  keep_gene <- rowSums(counts >= effective_min_count, na.rm = TRUE) >= effective_min_samples
  if (sum(keep_gene) < 2) stop("Too few genes remain after the low-count filter. Check the assigned groups or use Extended/custom filtering.")
  counts_f <- counts[keep_gene, , drop = FALSE]

  design_terms <- c(covariates, "condition")
  full_formula <- stats::as.formula(paste("~", paste(design_terms, collapse = " + ")))
  reduced_formula <- if (length(covariates)) stats::as.formula(paste("~", paste(covariates, collapse = " + "))) else ~ 1
  rownames(sm) <- sm$Accession

  list(
    counts = counts_f,
    sample_meta = sm,
    groups = groups_to_use,
    group_sizes = counts_per_group,
    smallest_group_size = smallest_group_size,
    covariates = covariates,
    full_formula = full_formula,
    reduced_formula = reduced_formula,
    keep_gene = keep_gene,
    prefilter_mode = prefilter_mode,
    min_count = effective_min_count,
    min_samples = effective_min_samples,
    filter_description = sprintf("raw count >= %d in at least %d samples", effective_min_count, effective_min_samples)
  )
}

run_geo2r_deseq <- function(counts, sample_meta, group, groups_to_use,
                            annotation = NULL, covariates = character(),
                            min_count = 10L, min_samples = NULL,
                            prefilter_mode = c("geo2r", "custom"),
                            alpha = 0.05, p_adjust = "fdr",
                            plot_contrasts = character()) {
  prefilter_mode <- match.arg(prefilter_mode)
  inp <- prepare_deseq_inputs(
    counts = counts, sample_meta = sample_meta, group = group,
    groups_to_use = groups_to_use, covariates = covariates,
    min_count = min_count, min_samples = min_samples,
    prefilter_mode = prefilter_mode
  )
  ng <- length(inp$groups)

  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(inp$counts),
    colData = inp$sample_meta,
    design = inp$full_formula
  )

  if (ng == 2L) {
    # GEO2R's two-group RNA-seq path uses a Wald test and positive-count size
    # factors. This is important for reproducing baseMean/LFC/stat/p-values.
    dds <- DESeq2::DESeq(dds0, test = "Wald", sfType = "poscounts", quiet = TRUE)
    comparison <- inp$groups[[1]]
    reference <- inp$groups[[2]]
    res <- DESeq2::results(
      dds,
      contrast = c("condition", comparison, reference),
      alpha = alpha,
      pAdjustMethod = p_adjust
    )
    res_meta <- tryCatch(S4Vectors::metadata(res), error = function(e) list())
    main <- standardize_geo2r_result(res, annotation = annotation, multigroup = FALSE)
    label <- paste(comparison, "vs", reference)
    pairwise <- list()
    pairwise[[label]] <- main
    test_name <- "Wald test"
    size_factor_type <- "poscounts"
    dds_for_norm <- dds
  } else {
    # GEO2R uses an overall LRT when 3+ groups are defined.
    dds <- DESeq2::DESeq(dds0, test = "LRT", reduced = inp$reduced_formula, quiet = TRUE)
    res <- DESeq2::results(dds, alpha = alpha, pAdjustMethod = p_adjust)
    res_meta <- tryCatch(S4Vectors::metadata(res), error = function(e) list())
    main <- standardize_geo2r_result(res, annotation = annotation, multigroup = TRUE)
    test_name <- "Likelihood ratio test (LRT)"
    size_factor_type <- "ratio"
    dds_for_norm <- dds

    requested <- intersect(plot_contrasts, contrast_labels(inp$groups))
    pairwise <- list()
    if (length(requested)) {
      # Re-run the fitted LRT object with a Wald test for GEO2R-style pairwise
      # plots. Existing LRT size factors are retained, as in GEO2R's script.
      dds_wald <- DESeq2::DESeq(dds, test = "Wald", sfType = "poscounts", quiet = TRUE)
      for (lab in requested) {
        cc <- parse_contrast_label(lab)
        if (is.null(cc)) next
        rr <- DESeq2::results(
          dds_wald,
          contrast = c("condition", cc$comparison, cc$reference),
          alpha = alpha,
          pAdjustMethod = p_adjust
        )
        pairwise[[lab]] <- standardize_geo2r_result(rr, annotation = annotation, multigroup = FALSE)
      }
      dds_for_norm <- dds_wald
    }
  }

  norm <- DESeq2::counts(dds_for_norm, normalized = TRUE)
  vst <- DESeq2::varianceStabilizingTransformation(dds_for_norm, blind = FALSE)

  list(
    dds = dds,
    dds_for_norm = dds_for_norm,
    result = main,
    pairwise = pairwise,
    normalized = norm,
    vst = vst,
    sample_meta = inp$sample_meta,
    kept_features = inp$keep_gene,
    features_before_prefilter = length(inp$keep_gene),
    features_after_prefilter = sum(inp$keep_gene),
    design = paste(deparse(inp$full_formula), collapse = ""),
    reduced_design = paste(deparse(inp$reduced_formula), collapse = ""),
    groups = inp$groups,
    group_sizes = inp$group_sizes,
    smallest_group_size = inp$smallest_group_size,
    covariates = inp$covariates,
    alpha = alpha,
    p_adjust = p_adjust,
    prefilter_mode = inp$prefilter_mode,
    filter_description = inp$filter_description,
    min_count = inp$min_count,
    min_samples = inp$min_samples,
    size_factor_type = size_factor_type,
    test_name = test_name,
    deseq2_version = as.character(utils::packageVersion("DESeq2")),
    pvalue_non_na = sum(!is.na(res$pvalue)),
    padj_non_na = sum(!is.na(res$padj)),
    independent_filter_threshold = if (!is.null(res_meta$filterThreshold) && length(res_meta$filterThreshold)) as.numeric(res_meta$filterThreshold[[1]]) else NA_real_,
    geo2r_parity = identical(inp$prefilter_mode, "geo2r") && !length(inp$covariates)
  )
}

make_volcano_data <- function(result, alpha = 0.05, lfc = 0) {
  x <- result
  x$neglog10padj <- -log10(pmax(x$padj, .Machine$double.xmin))
  x$Significance <- "Not significant"
  x$Significance[!is.na(x$padj) & x$padj < alpha & !is.na(x$log2FoldChange) & x$log2FoldChange >= lfc] <- "Up"
  x$Significance[!is.na(x$padj) & x$padj < alpha & !is.na(x$log2FoldChange) & x$log2FoldChange <= -lfc] <- "Down"
  x
}

pca_dataframe <- function(vst, sample_meta) {
  m <- SummarizedExperiment::assay(vst)
  rv <- matrixStats::rowVars(m)
  take <- order(rv, decreasing = TRUE)[seq_len(min(500L, length(rv)))]
  pc <- stats::prcomp(t(m[take, , drop = FALSE]), scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  out <- data.frame(
    Accession = rownames(pc$x),
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    stringsAsFactors = FALSE
  )
  sm <- sample_meta[match(out$Accession, sample_meta$Accession), , drop = FALSE]
  out$Group <- as.character(sm$condition)
  out$Title <- sm$Title %||% sm$Accession
  attr(out, "pct") <- pct[1:2]
  out
}

library_size_data <- function(counts, sample_meta, groups) {
  sizes <- colSums(counts)
  sm <- sample_meta[match(names(sizes), sample_meta$Accession), , drop = FALSE]
  data.frame(
    Accession = names(sizes),
    LibrarySize = as.numeric(sizes),
    Group = unname(groups[names(sizes)]),
    Title = sm$Title,
    stringsAsFactors = FALSE
  )
}

profile_dataframe <- function(normalized_counts, sample_meta, gene_id, result_table = NULL) {
  gene_id <- trimws(gene_id %||% "")
  if (!nzchar(gene_id)) stop("Enter a gene ID or symbol.")
  rid <- rownames(normalized_counts)
  match_id <- match(gene_id, rid)

  if (is.na(match_id) && !is.null(result_table) && nrow(result_table)) {
    for (nm in intersect(c("GeneID", "Symbol", "ID"), names(result_table))) {
      idx <- which(tolower(trim_chr(result_table[[nm]])) == tolower(gene_id))
      if (length(idx)) {
        candidate <- result_table$ID[idx[[1]]]
        match_id <- match(candidate, rid)
        if (!is.na(match_id)) break
      }
    }
  }
  if (is.na(match_id)) stop("Gene/feature was not found in the analyzed count matrix.")

  vals <- normalized_counts[match_id, ]
  sm <- sample_meta[match(names(vals), sample_meta$Accession), , drop = FALSE]
  data.frame(
    Accession = names(vals),
    Title = sm$Title,
    Group = as.character(sm$condition),
    NormalizedCount = as.numeric(vals),
    stringsAsFactors = FALSE
  )
}

generate_geo2r_r_script <- function(gse_id, organism, groups, assignments,
                                    p_adjust_label, alpha, lfc_threshold,
                                    min_count = 10L, min_samples = NULL, covariates = character(),
                                    prefilter_mode = "geo2r",
                                    test_name = NULL, count_source = NULL) {
  rq <- function(x) paste0("\"", gsub("\"", "\\\"", as.character(x), fixed = TRUE), "\"")
  used <- assignments[assignments != "-"]
  grp_txt <- paste(sprintf("%s = %s", names(used), unname(used)), collapse = ", ")
  padj_method <- geo2r_padjust_method(p_adjust_label)
  groups <- safe_groups(groups)
  custom_n <- if (is.null(min_samples) || !length(min_samples) || is.na(min_samples)) "NA_integer_" else as.character(as.integer(min_samples))
  lines <- c(
    "# GEO2R Multi-Species — RNA-seq DESeq2 core analysis",
    sprintf("# GSE: %s", gse_id %||% ""),
    sprintf("# Organism: %s", organism %||% ""),
    sprintf("# Count source: %s", count_source %||% "not loaded"),
    sprintf("# Test: %s", test_name %||% "determined by number of groups"),
    "",
    "library(DESeq2)",
    "",
    sprintf("groups <- c(%s)", paste(rq(groups), collapse = ", ")),
    sprintf("alpha <- %s", format(alpha, scientific = FALSE)),
    sprintf("p_adjust_method <- %s", rq(padj_method)),
    sprintf("covariates <- c(%s)", paste(rq(covariates), collapse = ", ")),
    sprintf("prefilter_mode <- %s", rq(prefilter_mode)),
    sprintf("custom_min_count <- %d", as.integer(min_count %||% 10L)),
    sprintf("custom_min_samples <- %s", custom_n),
    "",
    "# 'tbl' must be the validated raw integer count matrix and sample_info",
    "# must contain a Group column in the same sample order as colnames(tbl).",
    "sample_info$Group <- factor(sample_info$Group, levels = groups)",
    "",
    "# Current GEO2R-style pre-filter: count >= 10 in at least N samples,",
    "# where N is the size of the smallest assigned group.",
    "if (prefilter_mode == \"geo2r\") {",
    "  N <- min(table(sample_info$Group))",
    "  keep <- rowSums(tbl >= 10) >= N",
    "} else {",
    "  N <- custom_min_samples",
    "  keep <- rowSums(tbl >= custom_min_count) >= N",
    "}",
    "tbl <- tbl[keep, , drop = FALSE]",
    "",
    "design_formula <- if (length(covariates))",
    "  as.formula(paste(\"~\", paste(c(covariates, \"Group\"), collapse = \" + \"))) else ~ Group",
    "ds <- DESeqDataSetFromMatrix(countData = round(tbl), colData = sample_info, design = design_formula)",
    "",
    "if (length(groups) == 2) {",
    "  # First defined group vs second defined group, matching GEO2R group-order semantics.",
    "  ds <- DESeq(ds, test = \"Wald\", sfType = \"poscounts\")",
    "  r <- results(ds, contrast = c(\"Group\", groups[1], groups[2]),",
    "               alpha = alpha, pAdjustMethod = p_adjust_method)",
    "} else {",
    "  reduced_formula <- if (length(covariates))",
    "    as.formula(paste(\"~\", paste(covariates, collapse = \" + \"))) else ~ 1",
    "  ds <- DESeq(ds, test = \"LRT\", reduced = reduced_formula)",
    "  r <- results(ds, alpha = alpha, pAdjustMethod = p_adjust_method)",
    "}",
    "",
    sprintf("# Assigned samples: %s", if (nzchar(grp_txt)) grp_txt else "none"),
    sprintf("# Plot-only |log2FC| threshold: %s", format(lfc_threshold, scientific = FALSE)),
    "# padj is produced by DESeq2::results() after independent filtering and p.adjust()."
  )
  paste(lines, collapse = "\n")
}

