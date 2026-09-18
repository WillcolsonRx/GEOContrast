# Static regression checks for v0.4.7 multi-platform support.
app <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
geo <- paste(readLines("R/geo_utils.R", warn = FALSE), collapse = "\n")
count <- paste(readLines("R/count_utils.R", warn = FALSE), collapse = "\n")

stopifnot(grepl("multi_platform_mode", app, fixed = TRUE))
stopifnot(grepl("selected_platforms <- reactive", app, fixed = TRUE))
stopifnot(grepl("~ Platform + condition", app, fixed = TRUE))
stopifnot(grepl("Group and Platform are confounded", app, fixed = TRUE))
stopifnot(grepl("prepare_ncbi_counts <- function(se, organism = NULL, platforms = NULL)", count, fixed = TRUE))
stopifnot(grepl("%in% platform", geo, fixed = TRUE))
cat("multi-platform regression checks passed\n")
