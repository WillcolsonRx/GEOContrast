source(file.path("R", "geo_utils.R"))
source(file.path("R", "table_filter_utils.R"))

x <- data.frame(
  Group = c("-", "Control", "Case", "-"),
  Accession = paste0("GSM", 1:4),
  Phenotype = c("CTRL", "CTRL", "MDD", "MDD"),
  Gender = c("Female", "Male", "Female", "Male"),
  Age = c("31", "42", "57", "64"),
  Title = paste("Sample", 1:4),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

z <- prepare_filterable_sample_table(x, groups = c("Control", "Case"))
stopifnot(is.factor(z$Group))
stopifnot(is.factor(z$Phenotype))
stopifnot(is.factor(z$Gender))
stopifnot(is.numeric(z$Age))
stopifnot(is.character(z$Accession))
stopifnot(nrow(z) == nrow(x))

cat("table_filter_smoke_test: OK\n")
