cran <- c("shiny", "DT", "openxlsx", "readr", "readxl", "ggplot2", "matrixStats")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
bioc <- c("GEOquery", "DESeq2", "SummarizedExperiment", "Biobase", "limma", "AnnotationDbi")
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

message("Package installation complete. Start the app with shiny::runApp().")
