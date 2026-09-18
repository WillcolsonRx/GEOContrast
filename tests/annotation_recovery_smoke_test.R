# Offline smoke tests for v0.4.5 annotation recovery.
source("R/geo_utils.R")
source("R/annotation_utils.R")
source("R/deseq_utils.R")

ids <- c("8714", "10878", "2622")

# Simulate an annotation table where FeatureID is unusable but GeneID is valid.
ann <- data.frame(
  FeatureID = c("", "", ""),
  GeneID = ids,
  Symbol = c("ABCC3", "CFHR3", "DRC4"),
  Description = c("ATP binding cassette subfamily C member 3", "complement factor H related 3", "dynein regulatory complex subunit 4"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
rownames(ann) <- paste0("row", seq_along(ids))

std <- standardize_geo2r_annotation(ann, ids)
stopifnot(identical(std$GeneID, ids))
stopifnot(all(nzchar(std$Symbol)))
stopifnot(all(nzchar(std$Description)))

# Missing primary annotation should be filled by a fallback without replacing
# already valid values.
primary <- data.frame(
  FeatureID = ids,
  GeneID = ids,
  Symbol = c("ABCC3", "", ""),
  Description = c("", "", ""),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
fallback <- data.frame(
  FeatureID = ids,
  GeneID = ids,
  Symbol = c("WRONG_SHOULD_NOT_REPLACE", "CFHR3", "DRC4"),
  Description = c("desc1", "desc2", "desc3"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
merged <- fill_annotation_blanks(primary, fallback, ids)
stopifnot(identical(merged$Symbol, c("ABCC3", "CFHR3", "DRC4")))
stopifnot(all(nzchar(merged$Description)))

cv <- annotation_coverage(merged, ids)
stopifnot(cv$Coverage[cv$Field == "Symbol"] == 1)
stopifnot(cv$Coverage[cv$Field == "Description"] == 1)

message("Annotation recovery smoke test passed.")
