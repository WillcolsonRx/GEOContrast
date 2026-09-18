# Optional local gene-annotation fallbacks for common organisms.
# The main app does not require these packages, but when installed they can
# fill Gene Symbol / Description fields if GEO/NCBI annotation is unavailable.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(
  c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"),
  ask = FALSE,
  update = FALSE
)
message("Optional human/mouse/rat annotation databases installed.")
