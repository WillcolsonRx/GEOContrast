# Run from the project root after installing dependencies:
# source("tests/synthetic_smoke_test.R")

stopifnot(requireNamespace("DESeq2", quietly = TRUE))
stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
stopifnot(requireNamespace("matrixStats", quietly = TRUE))
source("R/geo_utils.R")
source("R/annotation_utils.R")
source("R/deseq_utils.R")

set.seed(7)
genes <- paste0("gene", seq_len(300))
samples <- paste0("GSM", 1001:1012)
counts <- matrix(rnbinom(300 * 12, mu = 80, size = 4), nrow = 300,
                 dimnames = list(genes, samples))
# Add a group-specific signal.
counts[1:20, 5:8] <- counts[1:20, 5:8] + 100
counts[21:40, 9:12] <- counts[21:40, 9:12] + 120
storage.mode(counts) <- "integer"

meta <- data.frame(
  Accession = samples,
  Title = paste("sample", seq_along(samples)),
  Batch = rep(c("A", "B"), 6),
  stringsAsFactors = FALSE
)
group <- setNames(rep(c("Control", "Disease", "Treatment"), each = 4), samples)
ann <- data.frame(
  FeatureID = genes,
  gene_id = genes,
  gene_name = paste0("SYM", seq_along(genes)),
  gene_type = "protein_coding",
  stringsAsFactors = FALSE
)

# Three-group GEO2R-like LRT.
lrt <- run_geo2r_deseq(
  counts = counts,
  sample_meta = meta,
  group = group,
  groups_to_use = c("Control", "Disease", "Treatment"),
  annotation = ann,
  covariates = "Batch",
  min_count = 5,
  min_samples = 2,
  prefilter_mode = "custom",
  alpha = 0.05,
  p_adjust = "fdr",
  plot_contrasts = c("Control vs Disease", "Disease vs Treatment")
)
stopifnot(lrt$test_name == "Likelihood ratio test (LRT)")
stopifnot(all(c("ID", "padj", "pvalue", "lfcSE", "stat", "log2FoldChange", "baseMean",
                "GeneID", "Symbol", "Description", "Synonyms", "GeneType") %in% names(lrt$result)))
stopifnot(all(is.na(lrt$result$log2FoldChange)))
stopifnot(length(lrt$pairwise) == 2)

# Two-group GEO2R-like Wald result.
wald <- run_geo2r_deseq(
  counts = counts[, 1:8],
  sample_meta = meta[1:8, ],
  group = group[1:8],
  groups_to_use = c("Disease", "Control"),
  annotation = ann,
  alpha = 0.05,
  p_adjust = "fdr",
  prefilter_mode = "geo2r"
)
stopifnot(wald$test_name == "Wald test")
stopifnot(length(wald$pairwise) == 1)
stopifnot(identical(names(wald$pairwise), "Disease vs Control"))
stopifnot(wald$size_factor_type == "poscounts")
stopifnot(wald$min_samples == 4L)  # smallest assigned group, GEO2R-style
stopifnot(any(is.finite(wald$result$log2FoldChange)))

# Explicit regression for GEO2R pre-filter logic with unequal group sizes.
small_counts <- matrix(0L, nrow = 4, ncol = 8,
                       dimnames = list(paste0("g", 1:4), paste0("S", 1:8)))
small_counts[1, 1:3] <- 10L  # keep: >=10 in smallest group size (3)
small_counts[2, 1:2] <- 10L  # drop
small_counts[3, 1:5] <- 20L  # keep
small_counts[4, ] <- 1L      # drop
small_meta <- data.frame(Accession = colnames(small_counts), stringsAsFactors = FALSE)
small_group <- setNames(c(rep("Test", 3), rep("Control", 5)), colnames(small_counts))
prep <- prepare_deseq_inputs(
  counts = small_counts,
  sample_meta = small_meta,
  group = small_group,
  groups_to_use = c("Test", "Control"),
  prefilter_mode = "geo2r"
)
stopifnot(prep$smallest_group_size == 3L)
stopifnot(prep$min_count == 10L)
stopifnot(prep$min_samples == 3L)
stopifnot(identical(as.logical(prep$keep_gene), c(TRUE, FALSE, TRUE, FALSE)))
stopifnot(identical(default_contrast_labels(c("Test", "Control")), "Test vs Control"))

cat("Synthetic GEO2R-parity DESeq2 smoke test passed.\n")

# v0.2.1 count-source helpers (offline tests)
stopifnot(grepl("GSE102556_raw_counts_GRCh38\\.p13_NCBI\\.tsv\\.gz", ncbi_human_raw_counts_url("GSE102556")))
stopifnot(has_normalized_only_candidates(data.frame(fname = c("study_FPKM_matrix.txt.gz"), stringsAsFactors = FALSE)))
stopifnot(!has_normalized_only_candidates(data.frame(fname = c("study_raw_counts.tsv.gz"), stringsAsFactors = FALSE)))
cat("v0.2.1 raw-count source helper tests passed.\n")

# v0.3 microarray / limma routing test.
stopifnot(requireNamespace("limma", quietly = TRUE))
stopifnot(requireNamespace("Biobase", quietly = TRUE))
source("R/limma_utils.R")

set.seed(11)
probes <- paste0("probe", seq_len(250))
msamples <- paste0("GSM", 2001:2008)
expr <- matrix(rnorm(250 * 8, mean = 7, sd = 0.6), nrow = 250,
               dimnames = list(probes, msamples))
expr[1:25, 5:8] <- expr[1:25, 5:8] + 1.2
pd <- data.frame(
  geo_accession = msamples,
  title = paste("microarray", seq_along(msamples)),
  source_name_ch1 = "brain",
  organism_ch1 = "Homo sapiens",
  row.names = msamples,
  stringsAsFactors = FALSE
)
fd <- data.frame(
  ID = probes,
  `Gene Symbol` = paste0("GENE", seq_along(probes)),
  `Gene Title` = paste("Synthetic gene", seq_along(probes)),
  row.names = probes,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
eset <- Biobase::ExpressionSet(
  assayData = expr,
  phenoData = Biobase::AnnotatedDataFrame(pd),
  featureData = Biobase::AnnotatedDataFrame(fd),
  annotation = "GPLSYNTH"
)
meta_m <- data.frame(
  Accession = msamples,
  Title = pd$title,
  `Source name` = pd$source_name_ch1,
  Organism = pd$organism_ch1,
  Platform = "GPLSYNTH",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
group_m <- setNames(rep(c("Control", "MDD"), each = 4), msamples)
lim <- run_geo2r_limma(
  obj = eset,
  sample_meta = meta_m,
  group = group_m,
  groups_to_use = c("Control", "MDD"),
  platform = "GPLSYNTH",
  log_transform = "No",
  force_normalization = FALSE,
  vooma = FALSE,
  p_adjust = "BH"
)
stopifnot(lim$engine == "limma")
stopifnot(all(c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "F") %in% names(lim$result)))
stopifnot(length(lim$pairwise) == 1)
stopifnot(any(is.finite(lim$result$logFC)))
cat("v0.3 microarray limma smoke test passed.\n")
