# Regression guard for v0.4.4.
# This is a source-level check: the sample table must use DT server-side mode
# when replaceData() is used, and group-level changes must force a clean rebuild.
app <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
stopifnot(grepl('output\\$samples_table <- renderDT', app))
stopifnot(grepl('server = TRUE\\)\\n  proxy <- dataTableProxy\\("samples_table"\\)', app))
stopifnot(grepl('replaceData\\(', app))
stopifnot(length(gregexpr('rv\\$table_filter_reset <- rv\\$table_filter_reset \\+ 1L', app, perl = TRUE)[[1]]) >= 3)
cat("DataTables proxy regression checks passed.\n")
