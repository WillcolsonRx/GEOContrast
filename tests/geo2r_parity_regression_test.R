# GEO2R RNA-seq parity regression checks.
# Run from project root:
# source("tests/geo2r_parity_regression_test.R")

source("R/geo_utils.R")
source("R/annotation_utils.R")
source("R/deseq_utils.R")

counts <- matrix(0L, nrow = 5, ncol = 8,
                 dimnames = list(paste0("g", 1:5), paste0("S", 1:8)))
counts[1, 1:3] <- 10L       # retained by GEO2R parity rule
counts[2, 1:2] <- 10L       # old app kept this (2 samples); parity must drop it
counts[3, 1:5] <- 20L       # retained
counts[4, ] <- 1L           # dropped
counts[5, 4:8] <- 11L       # retained
meta <- data.frame(Accession = colnames(counts), stringsAsFactors = FALSE)
group <- setNames(c(rep("Test", 3), rep("Control", 5)), colnames(counts))

x <- prepare_deseq_inputs(
  counts = counts,
  sample_meta = meta,
  group = group,
  groups_to_use = c("Test", "Control"),
  prefilter_mode = "geo2r"
)

stopifnot(x$smallest_group_size == 3L)
stopifnot(x$min_count == 10L)
stopifnot(x$min_samples == 3L)
stopifnot(identical(as.logical(x$keep_gene), c(TRUE, FALSE, TRUE, FALSE, TRUE)))
stopifnot(identical(default_contrast_labels(c("Test", "Control")), "Test vs Control"))
stopifnot(identical(geo2r_padjust_method("Benjamini & Hochberg (False discovery rate)"), "fdr"))

# Ensure the two-group code path explicitly requests the GEO2R-style estimator.
body_txt <- paste(deparse(body(run_geo2r_deseq)), collapse = "\n")
stopifnot(grepl('sfType = "poscounts"', body_txt, fixed = TRUE))
stopifnot(grepl('comparison <- inp$groups[[1]]', body_txt, fixed = TRUE))
stopifnot(grepl('reference <- inp$groups[[2]]', body_txt, fixed = TRUE))

cat("GEO2R parity regression checks passed.\n")
