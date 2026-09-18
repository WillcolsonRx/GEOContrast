# Regression test for GPL17027-style annotation normalization.
# Run from the project root after installing packages.
source("R/geo_utils.R")
source("R/annotation_utils.R")
source("R/limma_utils.R")

fake <- data.frame(
  ID = c("10_at", "100_at", "1000_at"),
  SPOT_ID = c("10", "100", "1000"),
  DESCRIPTION = c("N-acetyltransferase 2", "adenosine deaminase", "cadherin 2"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(identical(microarray_annotation_gene_ids(fake), c("10", "100", "1000")))

# Non-numeric SPOT_ID must not be silently treated as an Entrez Gene ID.
fake2 <- fake
fake2$SPOT_ID <- c("spot-A", "spot-B", "spot-C")
stopifnot(all(microarray_annotation_gene_ids(fake2) == ""))

cat("microarray annotation regression test passed\n")
