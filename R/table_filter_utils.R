# Helpers for the simple GEO2R-style sample table filters.
# Low-cardinality categorical metadata are converted to factors so DT renders
# dropdown filters in the column-filter row. Strict numeric metadata are kept
# numeric so DT can provide numeric/range filters. High-cardinality fields stay
# character and receive a normal text search box.

strict_numeric_metadata <- function(x, field_name = "") {
  z <- trim_chr(x)
  keep <- nzchar(z)
  if (!any(keep)) return(FALSE)

  # Only convert columns that are plausibly numeric metadata. This avoids
  # accidentally turning IDs/titles containing digits into numeric columns.
  numeric_name <- grepl(
    "age|pmi|post.?mortem|ph$|bmi|height|weight|dose|duration|time|score|count|percent|percentage|level|concentration",
    field_name, ignore.case = TRUE, perl = TRUE
  )
  if (!numeric_name) return(FALSE)

  vals <- suppressWarnings(as.numeric(z[keep]))
  mean(is.finite(vals)) >= 0.90 && length(unique(vals[is.finite(vals)])) >= 3
}

prepare_filterable_sample_table <- function(df, groups = character(), max_factor_levels = 60L) {
  out <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(out) || !ncol(out)) return(out)

  for (nm in names(out)) {
    if (identical(nm, "Accession")) next

    if (identical(nm, "Group")) {
      z <- trim_chr(out[[nm]])
      lev <- unique(c("-", groups, z[nzchar(z)]))
      out[[nm]] <- factor(z, levels = lev)
      next
    }

    if (strict_numeric_metadata(out[[nm]], nm)) {
      z <- trim_chr(out[[nm]])
      z[!nzchar(z)] <- NA_character_
      out[[nm]] <- suppressWarnings(as.numeric(z))
      next
    }

    z <- trim_chr(out[[nm]])
    nonblank <- z[nzchar(z)]
    n_unique <- length(unique(nonblank))

    # DT renders factor columns with a dropdown/select filter. Limit this to
    # informative low-cardinality variables so Title/GSM-like columns do not
    # become enormous dropdown menus.
    if (n_unique >= 2L && n_unique <= max_factor_levels) {
      lev <- unique(z)
      out[[nm]] <- factor(z, levels = lev)
    } else {
      out[[nm]] <- z
    }
  }

  out
}
