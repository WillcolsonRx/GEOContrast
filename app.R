options(shiny.maxRequestSize = 500 * 1024^2)

required <- c("shiny", "DT", "openxlsx", "GEOquery", "DESeq2", "SummarizedExperiment", "Biobase", "limma", "readr", "readxl", "ggplot2", "matrixStats")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), ". Run source('install_packages.R') first.")
}

library(shiny)
library(DT)
library(ggplot2)
source("R/geo_utils.R")
source("R/annotation_utils.R")
source("R/count_utils.R")
source("R/deseq_utils.R")
source("R/limma_utils.R")
source("R/table_filter_utils.R")

excel_widths <- function(df) {
  vapply(seq_along(df), function(i) {
    vals <- as.character(df[[i]])
    vals[is.na(vals)] <- ""
    min(45, max(11, max(c(nchar(vals), nchar(names(df)[i]), 9), na.rm = TRUE) + 2))
  }, numeric(1))
}

excel_safe_df <- function(df) {
  out <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in names(out)) {
    if (is.list(out[[nm]]) && !is.data.frame(out[[nm]])) {
      out[[nm]] <- vapply(out[[nm]], function(z) paste(as.character(z), collapse = "; "), character(1))
    }
  }
  out
}

add_excel_sheet <- function(wb, name, df, style = "TableStyleMedium2") {
  name <- substr(gsub("[\\/:*?\"<>|]", "_", name), 1, 31)
  openxlsx::addWorksheet(wb, name)
  df <- excel_safe_df(df)
  if (nrow(df) && ncol(df)) {
    openxlsx::writeDataTable(wb, name, df, tableStyle = style)
    openxlsx::freezePane(wb, name, firstRow = TRUE)
    openxlsx::setColWidths(wb, name, cols = seq_len(ncol(df)), widths = excel_widths(df))
  } else {
    openxlsx::writeData(wb, name, df)
  }
}

rna_geo2r_result_columns <- c(
  "ID", "padj", "pvalue", "lfcSE", "stat", "log2FoldChange", "baseMean",
  "GeneID", "Symbol", "Description", "Synonyms", "GeneType", "EnsemblGeneID", "Status",
  "ChrAcc", "ChrStart", "ChrStop", "Orientation", "Length", "GOFunctionID", "GOProcessID",
  "GOComponentID", "GOFunction", "GOProcess", "GOComponent"
)

rna_data_result_columns <- c("padj", "pvalue", "lfcSE", "stat", "log2FoldChange", "baseMean")
microarray_data_result_columns <- c("adj.P.Val", "P.Value", "t", "B", "logFC", "F")
annotation_result_columns <- c(
  "GeneID", "Symbol", "Description", "Synonyms", "GeneType", "EnsemblGeneID", "Status",
  "ChrAcc", "ChrStart", "ChrStop", "Orientation", "Length", "GOFunctionID", "GOProcessID",
  "GOComponentID", "GOFunction", "GOProcess", "GOComponent"
)

ui <- fluidPage(
  tags$head(
    tags$title("GEO2R — Multi-Species"),
    includeCSS("www/styles.css")
  ),
  div(class = "topbar", "NCBI GEO-style workflow · microarray + RNA-seq across organisms"),
  div(class = "app-title", "GEO2R — Multi-Species"),
  div(class = "app-subtitle",
      "Define GEO sample groups and automatically use limma for microarrays or DESeq2 for RNA-seq raw counts."),

  div(class = "accession-row",
      textInput("gse_id", "GEO accession", value = "", placeholder = "e.g. GSE102556", width = "245px"),
      actionButton("load_gse", "Set", class = "geo-btn-primary"),
      uiOutput("platform_ui"),
      uiOutput("organism_ui"),
      uiOutput("series_status_ui")
  ),
  uiOutput("series_title_ui"),

  tabsetPanel(
    id = "geo_tabs",

    tabPanel("GEO2R",
      div(class = "quick-start",
          tags$b("Quick start: "),
          "load a GSE, use the dropdown/search filters under the sample-table columns, assign the filtered samples to groups, then click Analyze. Use Options for analysis settings."
      ),

      div(class = "geo-panel",
          div(class = "geo-panel-header",
              div(span("▾ Samples"), "  ", actionLink("define_groups", "Define groups")),
              textOutput("selected_count", inline = TRUE)
          ),
          div(class = "geo-panel-body",
              div(class = "toolbar",
                  selectInput("assign_group", "Assign to group", choices = c("Group 1", "Group 2"), width = "205px"),
                  actionButton("assign_selected", "Assign selected", class = "geo-btn-primary"),
                  actionButton("unassign_selected", "Unassign selected", class = "geo-btn"),
                  actionButton("select_filtered", "Select filtered rows", class = "geo-btn"),
                  actionButton("clear_selection", "Clear selected rows", class = "geo-btn"),
                  actionButton("auto_group", "Auto-group from metadata…", class = "geo-btn"),
                  downloadButton("download_sample_excel", "Download sample table Excel", class = "btn-excel")
              ),
              uiOutput("group_status_ui"),

              div(class = "table-filter-strip",
                  div(class = "table-filter-copy",
                      tags$b("Filter samples by metadata"),
                      span("Use the filter controls directly under each column heading. Categorical variables such as phenotype, sex, medication, smoking or alcohol appear as dropdowns; numeric fields can be filtered by range; free-text fields remain searchable.")
                  ),
                  div(class = "table-filter-actions",
                      actionButton("reset_table_filters", "Clear all filters", class = "geo-btn"),
                      actionButton("assign_filtered", "Assign filtered samples", class = "geo-btn-primary"),
                      actionButton("unassign_filtered", "Unassign filtered samples", class = "geo-btn"),
                      uiOutput("table_filter_status_ui", inline = TRUE)
                  )
              ),

              DTOutput("samples_table")
          )
      ),

      uiOutput("count_details_ui"),

      div(class = "analyze-row",
          actionButton("run_analysis", "Analyze", class = "analyze-btn", icon = icon("play")),
          uiOutput("analysis_ready_ui")
      ),

      uiOutput("results_section_ui")
    ),

    tabPanel("Options",
      div(class = "options-grid",
          div(class = "option-column",
              tags$h4("Apply adjustment to the P-values."),
              radioButtons(
                "p_adjust_label", NULL,
                choices = c(
                  "Benjamini & Hochberg (False discovery rate)",
                  "Benjamini & Yekutieli",
                  "Bonferroni",
                  "Hochberg",
                  "Holm",
                  "Hommel"
                ),
                selected = "Benjamini & Hochberg (False discovery rate)"
              )
          ),
          div(class = "option-column",
              uiOutput("assay_options_ui"),
              uiOutput("covariate_ui")
          ),
          div(class = "option-column",
              tags$h4("Result annotation columns"),
              uiOutput("annotation_options_ui")
          ),
          div(class = "option-column",
              tags$h4("Plot displays"),
              numericInput("plot_alpha", "Significance level cut-off", value = 0.05, min = 0.000001, max = 1, step = 0.01),
              numericInput("plot_lfc", "Log2 fold change threshold", value = 0, min = 0, step = 0.25),
              selectizeInput("plot_contrasts", "Volcano and Mean-difference plot contrasts (up to 5)", choices = character(), multiple = TRUE,
                             options = list(maxItems = 5, placeholder = "Define two or more groups first"))
          )
      ),
      div(class = "reanalyze-row",
          tags$span(class = "small-muted", "If you edit Options after performing an analysis, click Reanalyze to apply the edits:"),
          actionButton("reanalyze", "Reanalyze", class = "geo-btn-primary")
      )
    ),

    tabPanel("Profile graph",
      div(class = "profile-controls",
          textInput("profile_gene", "Enter gene symbol or ID:", value = "", width = "330px"),
          actionButton("set_profile", "Set", class = "geo-btn-primary")
      ),
      div(class = "small-muted",
          "The profile graph uses the active GEO2R expression scale: processed/log-transformed values for microarrays and DESeq2 normalized counts for RNA-seq."),
      uiOutput("profile_status_ui"),
      plotOutput("profile_plot", height = "500px"),
      DTOutput("profile_table")
    ),

    tabPanel("R script",
      div(class = "script-toolbar", downloadButton("download_r_script", "Download R script", class = "geo-btn")),
      verbatimTextOutput("r_script_text")
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    gse_id = NULL,
    gse_list = NULL,
    metadata = NULL,
    series_info = NULL,
    supplementary = data.frame(),
    ncbi_available = FALSE,
    ncbi_attempt_error = NULL,
    groups = c("Group 1", "Group 2"),
    assignments = NULL,
    local_count_path = NULL,
    local_count_name = NULL,
    local_df = NULL,
    local_sheets = character(),
    count_obj = NULL,
    count_error = NULL,
    count_note = NULL,
    analysis = NULL,
    loaded_at = NULL,
    result_columns = rna_geo2r_result_columns,
    profile_id = NULL,
    annotation_status = NULL,
    table_filter_reset = 0L
  )

  cache_dir <- file.path(getwd(), "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  reset_analysis <- function() {
    rv$analysis <- NULL
    rv$profile_id <- NULL
  }
  reset_counts <- function() {
    rv$count_obj <- NULL
    rv$count_error <- NULL
    rv$ncbi_attempt_error <- NULL
    rv$count_note <- NULL
    rv$annotation_status <- NULL
    rv$local_count_path <- NULL
    rv$local_count_name <- NULL
    rv$local_df <- NULL
    rv$local_sheets <- character()
    reset_analysis()
  }

  observeEvent(input$load_gse, {
    id <- toupper(trimws(input$gse_id %||% ""))
    if (!grepl("^GSE[0-9]+$", id)) {
      showNotification("Enter a valid GSE accession such as GSE102556.", type = "error", duration = 6)
      return()
    }

    rv$gse_id <- id
    rv$gse_list <- NULL
    rv$metadata <- NULL
    rv$series_info <- NULL
    rv$supplementary <- data.frame()
    rv$ncbi_attempt_error <- NULL
    rv$groups <- c("Group 1", "Group 2")
    rv$assignments <- NULL
    rv$result_columns <- rna_geo2r_result_columns
    rv$table_filter_reset <- rv$table_filter_reset + 1L
    reset_counts()

    res <- tryCatch({
      withProgress(message = paste("Loading", id), value = 0.08, {
        incProgress(0.18, detail = "Retrieving GEO sample metadata…")
        gl <- GEOquery::getGEO(id, GSEMatrix = TRUE, getGPL = FALSE, destdir = cache_dir)
        if (!is.list(gl)) gl <- list(gl)
        incProgress(0.30, detail = "Building GEO2R-style sample table…")
        md <- build_series_metadata(gl)
        assay_types <- vapply(gl, function(obj) {
          plat <- extract_platform(obj)
          msub <- md[md$Platform == plat, , drop = FALSE]
          detect_assay_type(obj, msub)
        }, character(1))
        has_rnaseq <- any(assay_types == "RNA-seq")
        nc <- FALSE
        sf <- data.frame()
        if (has_rnaseq) {
          incProgress(0.18, detail = "Checking NCBI RNA-seq counts…")
          nc <- has_ncbi_rnaseq_counts(id)
          incProgress(0.18, detail = "Listing RNA-seq supplementary files…")
          sf <- list_supplementary_files(id)
        } else {
          incProgress(0.30, detail = "Microarray Series Matrix detected…")
        }
        list(gl = gl, md = md, nc = nc, sf = sf, assay_types = assay_types)
      })
    }, error = function(e) e)

    if (inherits(res, "error")) {
      showNotification(paste("Could not load", id, "—", conditionMessage(res)), type = "error", duration = NULL)
      return()
    }

    rv$gse_list <- res$gl
    rv$metadata <- res$md
    rv$ncbi_available <- res$nc
    rv$supplementary <- res$sf
    rv$series_info <- series_summary(id, res$gl, res$md)
    rv$loaded_at <- Sys.time()
    rv$assignments <- setNames(rep("-", nrow(res$md)), res$md$Accession)
    updateSelectInput(session, "assign_group", choices = rv$groups, selected = rv$groups[[1]])
    # Every newly loaded Series starts in the familiar single-platform mode.
    # Multi-platform selection is opt-in for each GSE.
    updateCheckboxInput(session, "multi_platform_mode", value = FALSE)

    orgs <- unique(trim_chr(res$md$Organism))
    orgs <- orgs[nzchar(orgs)]
    if (length(orgs)) updateSelectInput(session, "organism", choices = orgs, selected = orgs[[1]])
    # Platform choices are rendered dynamically after the organism selector is available.

    showNotification(paste(id, "loaded. Define groups and click Analyze."), type = "message", duration = 4)
  }, ignoreInit = TRUE)

  output$organism_ui <- renderUI({
    req(rv$metadata)
    orgs <- unique(trim_chr(rv$metadata$Organism))
    orgs <- orgs[nzchar(orgs)]
    if (!length(orgs)) return(div(class = "status-warn", "Organism unavailable"))
    selectInput("organism", "Organism", choices = orgs, selected = orgs[[1]], width = "250px")
  })

  platform_choices_for_organism <- reactive({
    req(rv$metadata)
    md <- rv$metadata
    org <- input$organism %||% ""
    if (nzchar(org) && "Organism" %in% names(md)) {
      md <- md[trim_chr(md$Organism) == org, , drop = FALSE]
    }
    plats <- unique(trim_chr(md$Platform))
    plats[nzchar(plats)]
  })

  output$platform_ui <- renderUI({
    req(rv$metadata)
    plats <- platform_choices_for_organism()
    if (!length(plats)) return(NULL)

    multi <- isTRUE(input$multi_platform_mode) && length(plats) > 1L
    # Isolate the current selection so adding/removing a GPL does not rebuild
    # the selectize control itself. The control only needs to rebuild when the
    # organism or single/multi mode changes.
    current <- as.character(isolate(input$platform) %||% character())
    current <- current[current %in% plats]
    if (!length(current)) current <- plats[[1]]
    if (!multi) current <- current[[1]]

    tagList(
      if (length(plats) > 1L) checkboxInput(
        "multi_platform_mode",
        "Select multiple platforms",
        value = multi
      ),
      selectizeInput(
        "platform",
        if (multi) "Platforms" else "Platform",
        choices = plats,
        selected = current,
        multiple = multi,
        width = if (multi) "330px" else "210px",
        options = if (multi) list(plugins = list("remove_button"), placeholder = "Choose one or more GPLs") else list()
      ),
      if (length(plats) > 1L) div(
        class = "small-muted",
        if (multi)
          "Multi-platform mode: the sample table shows the union of samples from every selected GPL."
        else
          "Default mode uses one GPL. Enable multi-platform mode to view/group samples across several GPLs."
      )
    )
  })

  selected_platforms <- reactive({
    req(rv$metadata)
    choices <- platform_choices_for_organism()
    if (!length(choices)) return(character())
    x <- as.character(input$platform %||% character())
    x <- unique(x[x %in% choices])
    if (!length(x)) x <- choices[[1]]
    if (!isTRUE(input$multi_platform_mode) && length(x) > 1L) x <- x[[1]]
    x
  })

  selected_platform_label <- reactive({
    x <- selected_platforms()
    if (length(x)) paste(x, collapse = "; ") else ""
  })

  selected_objects <- reactive({
    req(rv$gse_list)
    plats <- selected_platforms()
    lapply(plats, function(p) selected_geo_object(rv$gse_list, p))
  })

  selected_object <- reactive({
    objs <- selected_objects()
    if (length(objs) != 1L) {
      stop("A single GEO platform object is required for this operation. Switch to single-platform mode.")
    }
    objs[[1]]
  })

  selected_assay_types <- reactive({
    req(rv$gse_list, rv$metadata)
    plats <- selected_platforms()
    org <- input$organism %||% ""
    vapply(plats, function(plat) {
      md <- rv$metadata
      md <- md[trim_chr(md$Platform) == plat, , drop = FALSE]
      if (nzchar(org)) md <- md[trim_chr(md$Organism) == org, , drop = FALSE]
      detect_assay_type(selected_geo_object(rv$gse_list, plat), md)
    }, character(1))
  })

  selected_assay <- reactive({
    types <- unique(selected_assay_types())
    if (length(types) == 1L) types[[1]] else "Mixed"
  })

  observeEvent(list(input$organism, input$platform, input$multi_platform_mode), {
    if (is.null(rv$metadata)) return()
    rv$table_filter_reset <- rv$table_filter_reset + 1L
    reset_counts()
  }, ignoreInit = TRUE)

  output$series_status_ui <- renderUI({
    req(rv$metadata)
    assay <- selected_assay()
    np <- length(selected_platforms())
    if (identical(assay, "Mixed")) {
      return(span(class = "status-warn", paste(np, "platforms selected · mixed assay types · combined sample view only")))
    }
    if (identical(assay, "Microarray")) {
      if (np > 1L) {
        return(span(class = "status-note", paste(np, "microarray platforms selected · grouping/export enabled · analyze one GPL at a time")))
      }
      return(span(class = "status-good", "Microarray Series Matrix · limma ready"))
    }
    org <- input$organism %||% ""
    prefix <- if (np > 1L) paste0("RNA-seq · ", np, " platforms selected · ") else "RNA-seq · "
    if (!is.null(rv$count_obj)) return(span(class = "status-good", paste0(prefix, "raw counts ready · DESeq2")))
    if (isTRUE(rv$ncbi_available)) return(span(class = "status-good", paste0(prefix, "NCBI raw counts detected")))
    if (identical(org, "Homo sapiens")) return(span(class = "status-note", paste0(prefix, "NCBI human counts will be tried automatically")))
    span(class = "status-warn", paste0(prefix, "supplementary/uploaded raw counts may be required"))
  })

  output$series_title_ui <- renderUI({
    req(rv$gse_list, rv$metadata)
    title <- extract_series_title(rv$gse_list[[1]])
    div(class = "series-line",
        tags$a(class = "series-title",
               href = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", rv$gse_id),
               target = "_blank", title),
        span(class = "series-meta",
             paste(nrow(visible_metadata()), "shown of", nrow(rv$metadata), "samples · selected platform(s):", selected_platform_label()))
    )
  })

  visible_metadata <- reactive({
    req(rv$metadata)
    md <- rv$metadata
    org <- input$organism %||% ""
    plats <- selected_platforms()
    if (nzchar(org) && "Organism" %in% names(md)) md <- md[trim_chr(md$Organism) == org, , drop = FALSE]
    if (length(plats) && "Platform" %in% names(md)) md <- md[trim_chr(md$Platform) %in% plats, , drop = FALSE]
    md
  })

  analysis_assignments <- reactive({
    req(rv$metadata, rv$assignments)
    a <- rv$assignments
    allowed <- visible_metadata()$Accession
    a[!names(a) %in% allowed] <- "-"
    a
  })

  full_sample_table <- reactive({
    md <- visible_metadata()
    a <- rv$assignments %||% setNames(rep("-", nrow(rv$metadata)), rv$metadata$Accession)
    data.frame(Group = unname(a[md$Accession]), md, stringsAsFactors = FALSE, check.names = FALSE)
  })

  sample_table <- reactive({
    full_sample_table()
  })

  filterable_sample_table <- reactive({
    prepare_filterable_sample_table(sample_table(), rv$groups)
  })

  filtered_table_rows <- reactive({
    req(rv$metadata)
    rows <- input$samples_table_rows_all
    if (is.null(rows)) seq_len(nrow(sample_table())) else as.integer(rows)
  })

  used_groups <- reactive({
    if (is.null(rv$metadata) || is.null(rv$assignments)) return(character())
    a <- analysis_assignments()
    present <- unique(unname(a[a != "-"]))
    rv$groups[rv$groups %in% present]
  })

  observe({
    groups <- used_groups()
    choices <- contrast_labels(groups)
    defaults <- default_contrast_labels(groups)
    updateSelectizeInput(session, "plot_contrasts", choices = choices, selected = defaults, server = TRUE)
  })

  output$selected_count <- renderText({
    if (is.null(rv$metadata)) return("Selected 0 out of 0 samples")
    n_selected <- length(input$samples_table_rows_selected %||% integer())
    n_total <- nrow(sample_table())
    rows_all <- input$samples_table_rows_all
    n_filtered <- if (is.null(rows_all)) n_total else length(rows_all)
    if (n_filtered < n_total) {
      paste0("Selected ", n_selected, " · Showing ", n_filtered, " of ", n_total, " samples")
    } else {
      paste0("Selected ", n_selected, " out of ", n_total, " samples")
    }
  })

  output$table_filter_status_ui <- renderUI({
    if (is.null(rv$metadata)) return(NULL)
    n_total <- nrow(sample_table())
    rows_all <- input$samples_table_rows_all
    n_filtered <- if (is.null(rows_all)) n_total else length(rows_all)
    span(
      class = if (n_filtered < n_total) "filter-count filter-count-active" else "filter-count",
      paste0("Showing ", format(n_filtered, big.mark = ","), " of ", format(n_total, big.mark = ","), " samples")
    )
  })

  output$group_status_ui <- renderUI({
    if (is.null(rv$metadata)) return(div(class = "status-note", "Enter a GSE accession and click Set."))
    t <- full_sample_table()
    div(class = "status-note",
        paste0("Groups: ", paste(rv$groups, collapse = ", "),
               " · Assigned ", sum(t$Group != "-"), "/", nrow(t),
               " · Unassigned ", sum(t$Group == "-")))
  })

  output$samples_table <- renderDT({
    req(rv$metadata)
    # Use server-side DT because replaceData()/reloadData() use DT's AJAX data
    # endpoint. This keeps group-assignment refreshes valid for large GEO Series
    # while preserving filters, sort order, page and selection.
    rv$table_filter_reset
    table_data <- isolate(filterable_sample_table())
    datatable(
      table_data, rownames = FALSE,
      filter = list(position = "top", clear = TRUE),
      selection = list(mode = "multiple", target = "row"), extensions = "Buttons",
      class = "compact stripe hover row-border",
      options = list(
        dom = "<'table-tools'Bf>rt<'table-footer'lip>",
        buttons = list(list(extend = "colvis", text = "Columns")),
        pageLength = 25,
        lengthMenu = list(c(10, 25, 50, 100, -1), c("10", "25", "50", "100", "All")),
        scrollX = TRUE,
        scrollY = "500px",
        scrollCollapse = TRUE,
        search = list(caseInsensitive = TRUE),
        language = list(
          search = "Search all columns:",
          searchPlaceholder = "GSM, title, tissue, diagnosis…",
          zeroRecords = "No samples match the current filters"
        )
      )
    )
  }, server = TRUE)
  proxy <- dataTableProxy("samples_table")

  # Update the displayed Group column without destroying the DataTable state.
  # This is important after filtering a large GEO Series: assigning/unassigning
  # must not wipe the filter recipe or silently lose the selected GSMs.
  refresh_sample_table <- function(clear_selection = FALSE) {
    req(rv$metadata)
    dat <- prepare_filterable_sample_table(full_sample_table(), rv$groups)
    DT::replaceData(
      proxy, dat, rownames = FALSE, resetPaging = FALSE,
      clearSelection = if (isTRUE(clear_selection)) "all" else "none"
    )
    invisible(dat)
  }

  safe_table_accessions <- function(rows) {
    rows <- as.integer(rows %||% integer())
    dat <- sample_table()
    rows <- rows[is.finite(rows) & rows >= 1L & rows <= nrow(dat)]
    unique(dat$Accession[rows])
  }

  observeEvent(input$assign_selected, {
    rows <- input$samples_table_rows_selected %||% integer()
    acc <- safe_table_accessions(rows)
    if (!length(acc)) return(showNotification("Select sample rows first.", type = "warning"))
    gr <- input$assign_group %||% ""
    if (!gr %in% rv$groups) return(showNotification("Choose a valid group first.", type = "warning"))
    a <- rv$assignments
    a[acc] <- gr
    rv$assignments <- a
    refresh_sample_table(clear_selection = FALSE)
    reset_analysis()
    showNotification(paste(length(acc), "selected samples assigned to", gr), type = "message", duration = 3)
  })

  observeEvent(input$unassign_selected, {
    rows <- input$samples_table_rows_selected %||% integer()
    acc <- safe_table_accessions(rows)
    if (!length(acc)) return(showNotification("Select sample rows first.", type = "warning"))
    a <- rv$assignments
    a[acc] <- "-"
    rv$assignments <- a
    refresh_sample_table(clear_selection = FALSE)
    reset_analysis()
    showNotification(paste(length(acc), "selected samples unassigned."), type = "message", duration = 3)
  })

  observeEvent(input$select_filtered, {
    req(rv$metadata)
    rows <- filtered_table_rows()
    rows <- rows[rows >= 1L & rows <= nrow(sample_table())]
    if (!length(rows)) return(showNotification("No samples match the current table filters.", type = "warning"))
    DT::selectRows(proxy, rows)
    showNotification(paste(length(rows), "filtered samples selected."), type = "message", duration = 3)
  })

  observeEvent(input$assign_filtered, {
    req(rv$metadata)
    rows <- filtered_table_rows()
    acc <- safe_table_accessions(rows)
    if (!length(acc)) return(showNotification("No samples match the current table filters.", type = "warning"))
    gr <- input$assign_group %||% ""
    if (!gr %in% rv$groups) return(showNotification("Choose a valid group first.", type = "warning"))
    a <- rv$assignments
    a[acc] <- gr
    rv$assignments <- a
    refresh_sample_table(clear_selection = FALSE)
    reset_analysis()
    showNotification(paste(length(acc), "filtered samples assigned to", gr), type = "message", duration = 4)
  })

  observeEvent(input$unassign_filtered, {
    req(rv$metadata)
    rows <- filtered_table_rows()
    acc <- safe_table_accessions(rows)
    if (!length(acc)) return(showNotification("No samples match the current table filters.", type = "warning"))
    a <- rv$assignments
    a[acc] <- "-"
    rv$assignments <- a
    refresh_sample_table(clear_selection = FALSE)
    reset_analysis()
    showNotification(paste(length(acc), "filtered samples unassigned."), type = "message", duration = 4)
  })

  observeEvent(input$reset_table_filters, {
    # A controlled table rebuild is the most reliable way to clear both the
    # global search and all per-column filters in client-side DT.
    rv$table_filter_reset <- rv$table_filter_reset + 1L
    showNotification("All sample-table filters cleared.", type = "message", duration = 3)
  })

  observeEvent(input$clear_selection, {
    DT::selectRows(proxy, NULL)
    showNotification("Selected rows cleared.", type = "message", duration = 2)
  })

  observeEvent(input$define_groups, {
    req(rv$metadata)
    showModal(modalDialog(
      title = "Define groups",
      tags$p("Enter one group name per line, in the order you want them used for comparisons."),
      textAreaInput("group_names_edit", NULL, value = paste(rv$groups, collapse = "\n"), rows = 8, width = "100%"),
      easyClose = TRUE,
      footer = tagList(modalButton("Cancel"), actionButton("save_groups", "Save groups", class = "btn-primary"))
    ))
  })

  observeEvent(input$save_groups, {
    g <- safe_groups(unlist(strsplit(input$group_names_edit %||% "", "[\r\n]+")))
    if (length(g) < 2) return(showNotification("Define at least two groups.", type = "error"))
    rv$groups <- g
    a <- rv$assignments
    a[!a %in% c(g, "-")] <- "-"
    rv$assignments <- a
    updateSelectInput(session, "assign_group", choices = g, selected = g[[1]])
    removeModal()
    # Group names are factor levels in the sample table. Rebuild the DT when
    # those levels change instead of calling replaceData(); this refreshes the
    # column-filter choices safely and avoids stale/invalid AJAX responses.
    rv$table_filter_reset <- rv$table_filter_reset + 1L
    reset_analysis()
    showNotification("Group definitions updated.", type = "message", duration = 3)
  })

  observeEvent(input$auto_group, {
    req(rv$metadata)
    md <- visible_metadata()
    candidates <- setdiff(names(md), c("Accession", "Platform", "Organism"))
    showModal(modalDialog(
      title = "Auto-group from sample metadata",
      selectInput("auto_group_column", "Metadata column", choices = candidates, width = "100%"),
      checkboxInput("auto_group_replace", "Replace current group definitions", TRUE),
      footer = tagList(modalButton("Cancel"), actionButton("apply_auto_group", "Create groups", class = "btn-primary"))
    ))
  })

  observeEvent(input$apply_auto_group, {
    md <- visible_metadata()
    col <- input$auto_group_column %||% ""
    if (!col %in% names(md)) return()
    vals <- trim_chr(md[[col]])
    distinct <- unique(vals[nzchar(vals)])
    if (length(distinct) < 2) return(showNotification("That column has fewer than two non-empty values.", type = "error"))
    if (length(distinct) > 30) return(showNotification("That column creates more than 30 groups; choose a more specific field.", type = "error"))
    rv$groups <- if (isTRUE(input$auto_group_replace)) distinct else unique(c(rv$groups, distinct))
    a <- rv$assignments
    a[md$Accession] <- ifelse(nzchar(vals), vals, "-")
    rv$assignments <- a
    updateSelectInput(session, "assign_group", choices = rv$groups, selected = rv$groups[[1]])
    removeModal()
    # Auto-grouping can introduce a new set of factor levels. Rebuild the DT
    # so its categorical filter controls are regenerated from the new groups.
    rv$table_filter_reset <- rv$table_filter_reset + 1L
    reset_analysis()
    showNotification(paste("Auto-grouped", sum(nzchar(vals)), "samples from", col), type = "message", duration = 4)
  })


  output$count_details_ui <- renderUI({
    req(rv$metadata)
    if (identical(selected_assay(), "Microarray")) return(NULL)
    if (identical(selected_assay(), "Mixed")) {
      return(div(class = "warning-box", "The selected platforms contain different assay types. You can view, filter, group and export their samples together, but joint expression analysis is disabled."))
    }
    tags$details(
      class = "count-details",
      open = if (!isTRUE(rv$ncbi_available) && !identical(input$organism %||% "", "Homo sapiens")) "open" else NULL,
      tags$summary("RNA-seq raw-count source"),
      div(class = "count-details-body",
          uiOutput("count_source_summary_ui"),
          uiOutput("count_source_controls_ui"),
          uiOutput("file_config_ui"),
          uiOutput("count_status_ui"),
          uiOutput("annotation_status_ui"),
          DTOutput("count_preview_table"),
          DTOutput("mapping_table")
      )
    )
  })

  output$count_source_summary_ui <- renderUI({
    req(rv$metadata)
    if (!is.null(rv$count_obj)) {
      return(div(class = "help-box", tags$b("Count matrix ready. "),
                 paste(rv$count_obj$source, "·", format(nrow(rv$count_obj$counts), big.mark = ","), "features ·",
                       ncol(rv$count_obj$counts), "matched samples.")))
    }
    org <- input$organism %||% ""
    if (isTRUE(rv$ncbi_available)) {
      return(div(class = "help-box", tags$b("NCBI-computed raw counts are available."),
                 " You can load them now, or simply click Analyze and the app will load them automatically."))
    }
    if (identical(org, "Homo sapiens")) {
      note <- if (!is.null(rv$ncbi_attempt_error)) paste0(" Last NCBI attempt: ", rv$ncbi_attempt_error) else ""
      return(div(class = "help-box", tags$b("Human RNA-seq: automatic NCBI fallback enabled."),
                 " GEOquery's availability check can be a false negative for some mixed-platform or mixed-organism Series. The app will try GEOquery and then NCBI's documented direct raw-count download when you click Analyze.", note))
    }
    normalized_note <- if (has_normalized_only_candidates(rv$supplementary)) {
      " Normalized files (for example TPM/FPKM) appear to be present, but they are not valid DESeq2 input."
    } else ""
    div(class = "warning-box", tags$b("No active raw-count matrix is available yet."),
        " Choose a submitter supplementary raw-count matrix or upload one.", normalized_note)
  })

  output$count_source_controls_ui <- renderUI({
    req(rv$metadata)
    tagList(
      if (isTRUE(rv$ncbi_available) || identical(input$organism %||% "", "Homo sapiens"))
        actionButton("load_ncbi_counts", "Try NCBI raw counts", class = "geo-btn-primary"),
      if (nrow(rv$supplementary)) tagList(
        tags$h5("GEO supplementary files"),
        DTOutput("supp_table"),
        actionButton("load_supp_file", "Download & inspect selected file", class = "geo-btn")
      ) else div(class = "small-muted", "No Series-level supplementary files were listed."),
      fileInput("upload_counts", "Or upload a raw-count matrix (CSV / TSV / TXT / XLSX)",
                accept = c(".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz", ".xlsx", ".xls"), width = "100%")
    )
  })

  output$supp_table <- renderDT({
    req(nrow(rv$supplementary) > 0)
    show <- rv$supplementary[, c("fname", "Assessment", "Score"), drop = FALSE]
    datatable(show, rownames = FALSE, selection = "single",
              options = list(pageLength = 6, dom = "tip", order = list(list(2, "desc")), scrollX = TRUE))
  }, server = FALSE)

  load_local_table <- function(path, original_name) {
    rv$local_count_path <- path
    rv$local_count_name <- original_name
    rv$local_sheets <- available_excel_sheets(path, original_name)
    sheet <- if (length(rv$local_sheets)) rv$local_sheets[[1]] else NULL
    df <- read_tabular_file(path, original_name, sheet = sheet)
    if (nrow(df) < 2 || ncol(df) < 2) stop("The selected file does not look like a rectangular count matrix.")
    rv$local_df <- df
    rv$count_obj <- NULL
    rv$count_error <- NULL
    reset_analysis()
  }

  observeEvent(input$load_supp_file, {
    rows <- input$supp_table_rows_selected %||% integer()
    if (length(rows) != 1) return(showNotification("Select one supplementary file first.", type = "warning"))
    rec <- rv$supplementary[rows, , drop = FALSE]
    nm <- rec$fname[[1]]
    if (!grepl("\\.(csv|tsv|txt|tab|xlsx|xls)(\\.gz)?$", tolower(nm))) {
      return(showNotification("This file type is not automatically parsed. Convert it to a rectangular CSV/TSV/XLSX count table and upload it.", type = "warning", duration = 8))
    }
    dest <- file.path(cache_dir, paste0(rv$gse_id, "_", basename(nm)))
    err <- tryCatch({
      withProgress(message = paste("Downloading", nm), value = 0.2, {
        utils::download.file(rec$url[[1]], destfile = dest, mode = "wb", quiet = TRUE)
        incProgress(0.6, detail = "Reading table…")
        load_local_table(dest, nm)
      })
      NULL
    }, error = function(e) e)
    if (inherits(err, "error")) showNotification(paste("Could not load file:", conditionMessage(err)), type = "error", duration = NULL)
  })

  observeEvent(input$upload_counts, {
    f <- input$upload_counts
    req(f)
    upload_copy <- file.path(cache_dir, paste0("upload_", basename(f$name)))
    file.copy(f$datapath, upload_copy, overwrite = TRUE)
    err <- tryCatch({ load_local_table(upload_copy, f$name); NULL }, error = function(e) e)
    if (inherits(err, "error")) showNotification(paste("Could not read upload:", conditionMessage(err)), type = "error", duration = NULL)
  })

  load_ncbi_now <- function() {
    req(rv$gse_id, rv$metadata)
    org <- input$organism %||% NULL
    plats <- selected_platforms()

    # First use GEOquery's native loader. This is the preferred path because it
    # retrieves NCBI annotation and sample metadata together.
    err_native <- NULL
    obj <- tryCatch({
      se <- load_ncbi_rnaseq_counts(rv$gse_id)
      prepare_ncbi_counts(se, organism = org, platforms = plats)
    }, error = function(e) {
      err_native <<- conditionMessage(e)
      NULL
    })

    # GEOquery can report a false negative / fail to choose the correct raw
    # count link for mixed-platform or mixed-organism human Series. NCBI
    # documents a stable direct file pattern for human GRCh38.p13 counts, so
    # try that route before asking the user for a supplementary matrix.
    if (is.null(obj) && identical(org, "Homo sapiens")) {
      obj <- tryCatch(
        load_ncbi_direct_human_counts(
          rv$gse_id,
          metadata = rv$metadata,
          organism = org,
          platform = plats,
          cache_dir = cache_dir
        ),
        error = function(e) {
          rv$ncbi_attempt_error <- conditionMessage(e)
          NULL
        }
      )
    }

    if (is.null(obj)) {
      msg <- rv$ncbi_attempt_error %||% err_native %||% "NCBI raw counts were not available for this selection."
      rv$ncbi_attempt_error <- msg
      stop(msg)
    }

    # Keep only samples belonging to the currently selected organism/platform(s).
    allowed <- visible_metadata()$Accession
    keep_cols <- intersect(colnames(obj$counts), allowed)
    if (length(keep_cols) < 2L) {
      stop("NCBI counts were loaded, but fewer than two count columns match the selected platform set.")
    }
    obj$counts <- obj$counts[, keep_cols, drop = FALSE]
    if (!is.null(obj$mapping) && is.data.frame(obj$mapping) && "GSM" %in% names(obj$mapping)) {
      obj$mapping <- obj$mapping[obj$mapping$GSM %in% keep_cols, , drop = FALSE]
    }

    annfix <- repair_rnaseq_annotation(
      annotation = obj$annotation,
      feature_ids = rownames(obj$counts),
      gse_id = rv$gse_id,
      organism = org,
      cache_dir = cache_dir
    )
    obj$annotation <- annfix$annotation
    rv$annotation_status <- annfix

    rv$count_obj <- obj
    rv$ncbi_available <- TRUE
    rv$ncbi_attempt_error <- NULL
    rv$count_error <- NULL
    rv$count_note <- "Loaded NCBI-computed raw counts."
    reset_analysis()
    invisible(TRUE)
  }

  observeEvent(input$load_ncbi_counts, {
    err <- tryCatch({
      withProgress(message = "Loading NCBI RNA-seq quantifications", value = 0.25, {
        load_ncbi_now()
        incProgress(0.65, detail = "Raw counts validated.")
      })
      NULL
    }, error = function(e) e)
    if (inherits(err, "error")) {
      rv$count_error <- conditionMessage(err)
      showNotification(paste("NCBI counts could not be loaded:", rv$count_error), type = "error", duration = NULL)
    }
  })

  output$file_config_ui <- renderUI({
    df <- rv$local_df
    if (is.null(df)) return(NULL)
    mapping <- match_count_columns(df, visible_metadata())
    gene_default <- suggest_gene_id_column(df, mapping)
    tagList(
      div(class = "geo-panel compact-panel",
          div(class = "geo-panel-header", paste("Inspecting", rv$local_count_name)),
          div(class = "geo-panel-body",
              if (length(rv$local_sheets) > 1) selectInput("excel_sheet", "Excel sheet", choices = rv$local_sheets, selected = rv$local_sheets[[1]], width = "320px"),
              selectInput("gene_id_col", "Feature / gene ID column", choices = names(df), selected = gene_default, width = "360px"),
              actionButton("use_local_counts", "Use this raw-count matrix", class = "geo-btn-primary"),
              div(class = "small-muted", "Sample columns are matched to GSM accessions first, then exact GEO sample titles. Non-integer normalized matrices are rejected.")
          ))
    )
  })

  observeEvent(input$excel_sheet, {
    req(rv$local_count_path, rv$local_count_name)
    if (!length(rv$local_sheets)) return()
    err <- tryCatch({
      rv$local_df <- read_tabular_file(rv$local_count_path, rv$local_count_name, sheet = input$excel_sheet)
      NULL
    }, error = function(e) e)
    if (inherits(err, "error")) showNotification(conditionMessage(err), type = "error")
  }, ignoreInit = TRUE)

  observeEvent(input$use_local_counts, {
    req(rv$local_df)
    err <- tryCatch({
      obj <- prepare_uploaded_counts(rv$local_df, visible_metadata(), input$gene_id_col, organism = input$organism %||% NULL)
      annfix <- repair_rnaseq_annotation(
        annotation = obj$annotation,
        feature_ids = rownames(obj$counts),
        gse_id = rv$gse_id,
        organism = input$organism %||% NULL,
        cache_dir = cache_dir
      )
      obj$annotation <- annfix$annotation
      rv$annotation_status <- annfix
      rv$count_obj <- obj
      rv$count_error <- NULL
      rv$count_note <- "Raw-count validation passed; gene annotation was checked automatically."
      reset_analysis()
      NULL
    }, error = function(e) e)
    if (inherits(err, "error")) {
      rv$count_error <- conditionMessage(err)
      showNotification(rv$count_error, type = "error", duration = NULL)
    }
  })

  output$count_preview_table <- renderDT({
    df <- rv$local_df
    req(df)
    datatable(head(df, 10), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10, dom = "t"))
  }, server = FALSE)

  output$mapping_table <- renderDT({
    df <- rv$local_df
    req(df, rv$metadata)
    mp <- match_count_columns(df, visible_metadata())
    mp <- mp[nzchar(mp$GSM) | mp$CountLike, , drop = FALSE]
    datatable(mp, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10, dom = "tip"))
  }, server = FALSE)

  output$count_status_ui <- renderUI({
    if (!is.null(rv$count_error)) return(div(class = "error-box", tags$b("Count data rejected: "), rv$count_error))
    if (is.null(rv$count_obj)) return(NULL)
    div(class = "status-good", paste("Active count source:", rv$count_obj$source))
  })

  output$annotation_status_ui <- renderUI({
    if (is.null(rv$count_obj) || is.null(rv$annotation_status)) return(NULL)
    cv <- rv$annotation_status$coverage
    sym <- cv$Coverage[cv$Field == "Symbol"]
    cls <- if (length(sym) && is.finite(sym[[1]]) && sym[[1]] >= 0.80) "status-good" else "status-warn"
    div(class = cls, format_annotation_status(rv$annotation_status))
  })

  output$result_annotation_status_ui <- renderUI({
    if (is.null(rv$analysis)) return(NULL)
    engine <- rv$analysis$engine %||% ""

    if (identical(engine, "DESeq2")) {
      if (is.null(rv$annotation_status)) return(NULL)
      cv <- rv$annotation_status$coverage
      sym <- cv$Coverage[cv$Field == "Symbol"]
      cls <- if (length(sym) && is.finite(sym[[1]]) && sym[[1]] >= 0.80) "help-box" else "warning-box"
      return(div(class = cls, tags$b("Gene annotation: "), format_annotation_status(rv$annotation_status)))
    }

    if (identical(engine, "limma")) {
      st <- rv$analysis$annotation_status %||% list()
      cov <- suppressWarnings(as.numeric(st$symbol_coverage %||% 0))
      if (!length(cov) || !is.finite(cov[[1]])) cov <- 0
      src <- st$recovered_from %||% character()
      cls <- if (cov[[1]] >= 0.80) "help-box" else "warning-box"
      txt <- paste0(
        "Symbol coverage: ", round(100 * cov[[1]], 1), "%",
        if (length(src)) paste0(" · recovered from ", paste(src, collapse = " + ")) else
          " · using GEO platform annotation"
      )
      return(div(class = cls, tags$b("Microarray annotation: "), txt))
    }

    NULL
  })

  output$assay_options_ui <- renderUI({
    if (is.null(rv$metadata)) {
      return(div(class = "small-muted", "Load a GSE to display assay-specific options."))
    }
    if (identical(selected_assay(), "Mixed")) {
      return(div(class = "warning-box", "Mixed assay types are selected. Sample viewing/grouping is available, but choose compatible platforms before differential-expression analysis."))
    }
    if (identical(selected_assay(), "Microarray")) {
      return(tagList(
        tags$h4("Apply log transformation to the data."),
        radioButtons("log_transform", NULL, choices = c("Auto-detect", "Yes", "No"), selected = "Auto-detect"),
        tags$h4("Apply limma precision weights (vooma)."),
        radioButtons("vooma_option", NULL, choices = c("Yes", "No"), selected = "No"),
        tags$h4("Force normalization."),
        radioButtons("force_normalization", NULL, choices = c("Yes", "No"), selected = "No"),
        if (length(selected_platforms()) > 1L)
          div(class = "warning-box", "Multiple microarray GPLs are selected. Their samples can be grouped/exported together, but their probe matrices are not merged. Switch to one GPL before limma analysis."),
        div(class = "small-muted", "Microarray analysis uses the GEO Series Matrix with limma, matching GEO2R's processed-data workflow."),
        tags$hr(),
        tags$h4("Model"),
        div(class = "small-muted", "Groups are modeled with limma using GEO2R-style group contrasts.")
      ))
    }
    tagList(
      tags$h4("Apply log transformation to the data."),
      div(class = "na-control", "Not applicable — DESeq2 models raw integer counts directly."),
      tags$h4("Apply limma precision weights (vooma)."),
      div(class = "na-control", "Not applicable for the DESeq2 RNA-seq engine."),
      tags$h4("Force normalization."),
      div(class = "na-control", "Not applicable — DESeq2 estimates size factors internally."),
      tags$hr(),
      if (length(selected_platforms()) > 1L)
        div(class = "warning-box", tags$b("Multi-platform RNA-seq mode: "),
            "the app will combine matched raw-count columns from the selected GPLs and include Platform in the DESeq2 design to account for platform-associated differences. Exact GEO2R parity is therefore not expected."),
      tags$h4("RNA-seq calculation mode"),
      radioButtons(
        "rnaseq_analysis_mode", NULL,
        choices = if (length(selected_platforms()) > 1L)
          c("GEO2R-like filtering + platform-adjusted design" = "geo2r", "Extended / custom" = "custom")
        else
          c("GEO2R parity (recommended)" = "geo2r", "Extended / custom" = "custom"),
        selected = "geo2r"
      ),
      div(class = "help-box",
          tags$b("GEO2R parity defaults: "),
          "keep genes with raw count ≥10 in at least N samples, where N is the smallest assigned group; use the first-defined group versus the second-defined group; for two groups use a DESeq2 Wald test with poscounts size factors."),
      conditionalPanel(
        condition = "input.rnaseq_analysis_mode == 'custom'",
        tags$h4("Custom low-count filter"),
        numericInput("min_count", "Minimum raw count", value = 10, min = 0, step = 1),
        numericInput("min_samples", "Minimum samples meeting count", value = 2, min = 1, step = 1),
        div(class = "small-muted", "Custom filtering and/or covariates intentionally may not reproduce official GEO2R values.")
      )
    )
  })

  output$annotation_options_ui <- renderUI({
    if (is.null(rv$metadata)) return(div(class = "small-muted", "Load a GSE first."))
    if (identical(selected_assay(), "Microarray")) {
      choices <- if (!is.null(rv$analysis) && identical(rv$analysis$engine, "limma")) rv$analysis$annotation_columns else character()
      if (!length(choices)) {
        return(div(class = "small-muted", "Platform annotation columns are loaded with the microarray analysis and will be available via ‘Show all columns’."))
      }
      return(checkboxGroupInput("annotation_display_microarray", NULL, choices = choices, selected = head(choices, 6)))
    }
    tagList(
      checkboxGroupInput("annotation_display", NULL, choices = annotation_result_columns,
                         selected = c("GeneID", "Symbol", "Description", "Synonyms", "GeneType")),
      div(class = "small-muted", "Unavailable RNA-seq annotation fields remain blank for that organism/source.")
    )
  })

  output$covariate_ui <- renderUI({
    if (is.null(rv$metadata)) return(NULL)
    if (identical(selected_assay(), "Microarray") || identical(selected_assay(), "Mixed")) return(NULL)
    multi <- length(selected_platforms()) > 1L
    if (!identical(input$rnaseq_analysis_mode %||% "geo2r", "custom")) {
      return(div(class = "small-muted",
                 if (multi)
                   "Platform is added automatically to the DESeq2 design in multi-platform mode. Additional covariates are disabled unless you switch to Extended / custom."
                 else
                   "Optional covariates are disabled in GEO2R parity mode because official GEO2R uses the user-defined sample Group as its analysis design. Switch to Extended / custom if you intentionally want to adjust for metadata covariates."))
    }
    md <- visible_metadata()
    candidates <- setdiff(names(md), c("Accession", "Title", "Source name", "Organism", "Platform"))
    tagList(
      if (multi) div(class = "small-muted", "Platform is included automatically; choose only additional covariates below."),
      selectizeInput("covariates", "Optional covariates (extended mode)", choices = candidates, multiple = TRUE,
                     options = list(placeholder = "e.g. Sex, Batch, Age"), width = "100%")
    )
  })

  output$analysis_ready_ui <- renderUI({
    if (is.null(rv$metadata)) return(span(class = "status-note", "Load a GSE first."))
    ug <- used_groups()
    if (length(ug) < 2) return(span(class = "status-warn", "Define and assign at least two groups."))
    assay <- selected_assay()
    np <- length(selected_platforms())
    if (identical(assay, "Mixed")) {
      return(span(class = "status-warn", "Selected GPLs contain mixed assay types. Joint analysis is disabled."))
    }
    if (identical(assay, "Microarray")) {
      if (np > 1L) return(span(class = "status-warn", "Multiple microarray GPLs are selected. Choose one GPL to run limma."))
      return(span(class = "status-good", paste(length(ug), "groups ready · Series Matrix will be analyzed with limma")))
    }
    can_try_ncbi <- isTRUE(rv$ncbi_available) || identical(input$organism %||% "", "Homo sapiens")
    source_txt <- if (!is.null(rv$count_obj)) {
      rv$count_obj$source
    } else if (can_try_ncbi) {
      "NCBI raw counts will be tried automatically"
    } else {
      "raw-count matrix still required"
    }
    span(class = if (!is.null(rv$count_obj) || can_try_ncbi) "status-good" else "status-warn",
         paste(length(ug), "groups ready ·", source_txt,
               if (np > 1L) "· Platform will be modeled automatically" else ""))
  })

  run_current_analysis <- function() {
    req(rv$metadata)
    ug <- used_groups()
    if (length(ug) < 2) stop("Define and assign at least two groups before analysis.")
    assay <- selected_assay()
    plats <- selected_platforms()
    if (identical(assay, "Mixed")) {
      stop("The selected GPLs contain different assay types. You can group/export their samples together, but joint differential-expression analysis is not valid.")
    }
    if (identical(assay, "Microarray") && length(plats) > 1L) {
      stop("Multiple microarray GPLs are selected. Their probe/expression matrices are platform-specific and are not merged directly. Switch to single-platform mode for limma analysis.")
    }
    active_assignment <- analysis_assignments()
    plot_labs <- intersect(input$plot_contrasts %||% character(), contrast_labels(ug))
    if (!length(plot_labs)) plot_labs <- default_contrast_labels(ug)

    if (identical(assay, "Microarray")) {
      incProgress(0.20, detail = "Loading GEO Series Matrix for limma…")
      obj <- selected_object()
      ana <- run_geo2r_limma(
        obj = obj,
        sample_meta = visible_metadata(),
        group = active_assignment,
        groups_to_use = ug,
        platform = plats[[1]] %||% extract_platform(obj),
        organism = input$organism %||% "",
        cache_dir = cache_dir,
        log_transform = input$log_transform %||% "Auto-detect",
        force_normalization = identical(input$force_normalization %||% "No", "Yes"),
        vooma = identical(input$vooma_option %||% "No", "Yes"),
        p_adjust = geo2r_padjust_method(input$p_adjust_label),
        plot_contrasts = plot_labs
      )
      rv$analysis <- ana
      chosen_ann <- input$annotation_display_microarray %||% head(ana$annotation_columns, 6)
      rv$result_columns <- unique(c("ID", microarray_data_result_columns, chosen_ann))
      return(ana)
    }

    if (is.null(rv$count_obj)) {
      org <- input$organism %||% ""
      should_try_ncbi <- isTRUE(rv$ncbi_available) || identical(org, "Homo sapiens")
      if (should_try_ncbi) {
        incProgress(0.15, detail = "Trying NCBI-computed raw counts…")
        tryCatch(load_ncbi_now(), error = function(e) NULL)
      }
    }
    if (is.null(rv$count_obj)) {
      norm_hint <- if (has_normalized_only_candidates(rv$supplementary)) {
        " GEO supplementary files appear to include normalized expression (for example FPKM/TPM), which cannot be used as DESeq2 raw counts."
      } else ""
      ncbi_hint <- if (!is.null(rv$ncbi_attempt_error)) paste0(" NCBI attempt: ", rv$ncbi_attempt_error) else ""
      stop(paste0(
        "This platform is RNA-seq, but no DESeq2-compatible raw-count matrix could be activated.",
        ncbi_hint, norm_hint,
        " Open 'RNA-seq raw-count source' and choose a true raw-count matrix or upload one."
      ))
    }

    active_counts <- rv$count_obj$counts

    # Re-check annotation immediately before analysis. This repairs stale or
    # partially downloaded annotation caches and fills missing GEO2R fields
    # without changing the expression/count data.
    annfix <- repair_rnaseq_annotation(
      annotation = rv$count_obj$annotation,
      feature_ids = rownames(active_counts),
      gse_id = rv$gse_id,
      organism = input$organism %||% NULL,
      cache_dir = cache_dir
    )
    rv$count_obj$annotation <- annfix$annotation
    rv$annotation_status <- annfix

    matched_groups <- unique(unname(active_assignment[intersect(colnames(active_counts), names(active_assignment))]))
    matched_groups <- matched_groups[matched_groups != "-"]
    ug <- rv$groups[rv$groups %in% matched_groups]
    if (length(ug) < 2) stop("The active count matrix does not contain at least two assigned groups for the selected organism/platform.")
    plot_labs <- intersect(input$plot_contrasts %||% character(), contrast_labels(ug))
    if (!length(plot_labs)) plot_labs <- default_contrast_labels(ug)

    incProgress(0.25, detail = if (length(ug) == 2) "Running DESeq2 Wald test…" else "Running DESeq2 likelihood-ratio test…")
    rnaseq_mode <- input$rnaseq_analysis_mode %||% "geo2r"
    multi_platform_rna <- length(plats) > 1L
    parity_requested <- identical(rnaseq_mode, "geo2r")
    parity_mode <- parity_requested && !multi_platform_rna
    auto_covariates <- if (multi_platform_rna) "Platform" else character()
    extra_covariates <- if (identical(rnaseq_mode, "custom")) (input$covariates %||% character()) else character()
    model_covariates <- unique(c(auto_covariates, extra_covariates))

    if (multi_platform_rna) {
      # Fail early with a clear explanation when biological group and platform
      # are completely confounded, because ~ Platform + condition is then not
      # estimable in DESeq2.
      sm_check <- rv$metadata[match(colnames(active_counts), rv$metadata$Accession), , drop = FALSE]
      sm_check$condition <- unname(active_assignment[match(sm_check$Accession, names(active_assignment))])
      sm_check <- sm_check[sm_check$condition %in% ug, , drop = FALSE]
      present_platforms <- unique(trim_chr(sm_check$Platform))
      present_platforms <- present_platforms[nzchar(present_platforms)]
      missing_platforms <- setdiff(plats, present_platforms)
      if (length(missing_platforms)) {
        stop("The active raw-count matrix does not contain assigned samples from all selected platforms. Missing: ", paste(missing_platforms, collapse = ", "), ". Load a count matrix covering all selected GPLs or change the platform selection.")
      }
      sm_check$condition <- factor(sm_check$condition, levels = ug)
      sm_check$Platform <- factor(trim_chr(sm_check$Platform), levels = plats)
      mm <- stats::model.matrix(~ Platform + condition, data = sm_check)
      if (qr(mm)$rank < ncol(mm)) {
        stop("Group and Platform are confounded for the selected samples (for example, one group occurs only on one GPL). Joint multi-platform DESeq2 cannot estimate a platform-adjusted group effect. Change the selected samples/platforms or analyze the platforms separately.")
      }
    }

    ana <- run_geo2r_deseq(
      counts = active_counts,
      sample_meta = rv$metadata,
      group = active_assignment,
      groups_to_use = ug,
      annotation = rv$count_obj$annotation,
      covariates = model_covariates,
      min_count = if (parity_requested) 10L else (input$min_count %||% 10L),
      min_samples = if (parity_requested) NULL else (input$min_samples %||% 2L),
      prefilter_mode = if (parity_requested) "geo2r" else "custom",
      alpha = input$plot_alpha %||% 0.05,
      p_adjust = geo2r_padjust_method(input$p_adjust_label),
      plot_contrasts = plot_labs
    )
    ana$multi_platform <- multi_platform_rna
    ana$platforms <- plats
    ana$platform_adjusted <- multi_platform_rna
    ana$engine <- "DESeq2"
    ana$assay_type <- "RNA-seq"
    ana$annotation_columns <- annotation_result_columns
    ana$export_columns <- rna_geo2r_result_columns
    ana$default_display_columns <- unique(c("ID", rna_data_result_columns, input$annotation_display %||% c("GeneID", "Symbol", "Description", "Synonyms", "GeneType")))
    rv$analysis <- ana
    rv$result_columns <- ana$default_display_columns
    ana
  }

  observeEvent(input$run_analysis, {
    err <- tryCatch({
      withProgress(message = paste("Analyzing", rv$gse_id %||% "GEO Series"), value = 0.08, {
        run_current_analysis()
        incProgress(0.60, detail = "Preparing GEO2R-style results…")
      })
      NULL
    }, error = function(e) e)
    if (inherits(err, "error")) showNotification(paste("Analysis failed:", conditionMessage(err)), type = "error", duration = NULL)
    else showNotification("Analysis complete.", type = "message", duration = 4)
  })

  observeEvent(input$reanalyze, {
    if (is.null(rv$metadata)) return(showNotification("Load a GSE first.", type = "warning"))
    err <- tryCatch({
      withProgress(message = "Reanalyzing with current Options", value = 0.10, {
        run_current_analysis()
        incProgress(0.65, detail = "Applying updated settings…")
      })
      NULL
    }, error = function(e) e)
    if (inherits(err, "error")) showNotification(paste("Reanalysis failed:", conditionMessage(err)), type = "error", duration = NULL)
    else {
      updateTabsetPanel(session, "geo_tabs", selected = "GEO2R")
      showNotification("Reanalysis complete.", type = "message", duration = 4)
    }
  })

  output$results_section_ui <- renderUI({
    if (is.null(rv$analysis)) return(NULL)
    ng <- length(rv$analysis$groups)
    engine <- rv$analysis$engine %||% "DESeq2"
    subtitle <- if (identical(engine, "limma")) {
      if (ng == 2) "limma moderated t-test · top 250 ranked by adjusted P-value" else "limma moderated F-test · top 250 ranked by adjusted P-value"
    } else {
      if (ng == 2) "DESeq2 Wald test · top 250 ranked by adjusted P-value" else "DESeq2 LRT overall test · top 250 ranked by adjusted P-value"
    }
    tagList(
      div(class = "results-heading-row",
          tags$h3("Top differentially expressed genes"),
          span(class = "small-muted", subtitle)
      ),
      div(class = "result-links",
          downloadButton("download_result_excel", "Download full table Excel", class = "btn-linkish"),
          downloadButton("download_result_tsv", "Download full table TSV", class = "btn-linkish"),
          actionLink("choose_result_columns", "Show all columns")
      ),
      uiOutput("result_test_note_ui"),
      uiOutput("result_annotation_status_ui"),
      DTOutput("results_table"),
      tags$h3("Visualization"),
      uiOutput("plot_contrast_ui"),
      fluidRow(
        column(6, plotOutput("volcano_plot", height = "440px")),
        column(6, plotOutput("md_plot", height = "440px"))
      ),
      fluidRow(
        column(6, plotOutput("padj_hist", height = "380px")),
        column(6, plotOutput("pca_plot", height = "380px"))
      )
    )
  })

  output$result_test_note_ui <- renderUI({
    req(rv$analysis)
    if (identical(rv$analysis$engine %||% "", "limma")) {
      prep <- c(
        if (isTRUE(rv$analysis$log_transform_applied)) "log2 transformation applied" else "no log2 transformation applied",
        if (isTRUE(rv$analysis$normalization_applied)) "forced between-array normalization applied" else "no forced normalization",
        if (isTRUE(rv$analysis$vooma_applied)) "vooma precision weights applied" else "vooma off"
      )
      if (length(rv$analysis$groups) > 2) {
        return(div(class = "help-box", tags$b("Microarray GEO2R-style analysis: "),
                   "the main table reports the overall limma moderated F-test; t, B and logFC are two-group/pairwise statistics. ", paste(prep, collapse = " · ")))
      }
      return(div(class = "help-box", tags$b("Contrast: "), names(rv$analysis$pairwise)[[1]],
                 " · limma moderated t-test · ", paste(prep, collapse = " · ")))
    }
    filter_thr <- rv$analysis$independent_filter_threshold %||% NA_real_
    parity_details <- paste0(
      if (isTRUE(rv$analysis$geo2r_parity)) "GEO2R parity mode" else "Extended/custom mode",
      " · pre-filter: ", rv$analysis$filter_description %||% "",
      " · retained ", format(rv$analysis$features_after_prefilter %||% nrow(rv$analysis$result), big.mark = ","),
      " / ", format(rv$analysis$features_before_prefilter %||% nrow(rv$analysis$result), big.mark = ","), " features",
      " · non-NA padj: ", format(rv$analysis$padj_non_na %||% NA_integer_, big.mark = ","),
      if (is.finite(filter_thr)) paste0(" · independent-filter mean threshold: ", signif(filter_thr, 5)) else "",
      " · DESeq2 ", rv$analysis$deseq2_version %||% ""
    )
    if (length(rv$analysis$groups) > 2) {
      div(class = "help-box",
          tags$b("RNA-seq GEO2R-style multi-group behavior: "),
          "the main table is an overall DESeq2 likelihood-ratio test. Accordingly, lfcSE and log2FoldChange are left blank in the overall table; pairwise fold changes are shown in the selected plot contrasts. ",
          parity_details)
    } else {
      div(class = "help-box", tags$b("Contrast: "), names(rv$analysis$pairwise)[[1]],
          " · DESeq2 Wald test · size factors: ", rv$analysis$size_factor_type %||% "",
          " · ", parity_details)
    }
  })

  results_display <- reactive({
    req(rv$analysis)
    cols <- intersect(rv$result_columns, names(rv$analysis$result))
    if (!length(cols)) cols <- intersect(rv$analysis$default_display_columns %||% names(rv$analysis$result), names(rv$analysis$result))
    rv$analysis$result[, cols, drop = FALSE]
  })

  output$results_table <- renderDT({
    x <- results_display()
    x <- head(x, 250)
    order_name <- if (identical(rv$analysis$engine %||% "", "limma")) "adj.P.Val" else "padj"
    order_col <- match(order_name, names(x))
    if (is.na(order_col)) order_col <- 1
    datatable(
      x, rownames = FALSE, filter = "top", extensions = "Buttons",
      options = list(
        scrollX = TRUE,
        pageLength = 25,
        lengthMenu = list(c(10, 25, 50, 100, 250), c("10", "25", "50", "100", "250")),
        dom = "Bfrtip",
        buttons = c("copy", "csv", "colvis"),
        order = list(list(order_col - 1, "asc"))
      )
    )
  })

  observeEvent(input$choose_result_columns, {
    req(rv$analysis)
    is_limma <- identical(rv$analysis$engine %||% "", "limma")
    data_choices <- if (is_limma) microarray_data_result_columns else rna_data_result_columns
    ann_choices <- if (is_limma) rv$analysis$annotation_columns %||% character() else annotation_result_columns
    showModal(modalDialog(
      title = "Select result columns",
      tags$h5("Data columns"),
      checkboxGroupInput("data_columns_modal", NULL, choices = data_choices,
                         selected = intersect(data_choices, rv$result_columns)),
      tags$h5("Annotation columns"),
      checkboxGroupInput("annotation_columns_modal", NULL, choices = ann_choices,
                         selected = intersect(ann_choices, rv$result_columns)),
      checkboxInput("include_id_modal", "ID", "ID" %in% rv$result_columns),
      footer = tagList(
        actionButton("restore_result_columns", "Restore defaults", class = "geo-btn"),
        modalButton("Cancel"),
        actionButton("apply_result_columns", "Set", class = "btn-primary")
      )
    ))
  })

  observeEvent(input$restore_result_columns, {
    req(rv$analysis)
    is_limma <- identical(rv$analysis$engine %||% "", "limma")
    data_choices <- if (is_limma) microarray_data_result_columns else rna_data_result_columns
    ann_choices <- if (is_limma) rv$analysis$annotation_columns %||% character() else annotation_result_columns
    updateCheckboxGroupInput(session, "data_columns_modal", choices = data_choices, selected = data_choices)
    updateCheckboxGroupInput(session, "annotation_columns_modal", choices = ann_choices, selected = ann_choices)
    updateCheckboxInput(session, "include_id_modal", value = TRUE)
  })

  observeEvent(input$apply_result_columns, {
    cols <- c(if (isTRUE(input$include_id_modal)) "ID", input$data_columns_modal %||% character(), input$annotation_columns_modal %||% character())
    if (!length(cols)) return(showNotification("Select at least one result column.", type = "warning"))
    rv$result_columns <- unique(cols)
    removeModal()
  })

  output$plot_contrast_ui <- renderUI({
    req(rv$analysis)
    labs <- names(rv$analysis$pairwise)
    if (!length(labs)) return(div(class = "warning-box", "No pairwise plot contrast was computed. Select contrasts in Options and click Reanalyze."))
    selectInput("active_plot_contrast", "Contrast for Volcano / Mean-difference plots", choices = labs, selected = labs[[1]], width = "360px")
  })

  active_pair_result <- reactive({
    req(rv$analysis)
    labs <- names(rv$analysis$pairwise)
    if (!length(labs)) return(NULL)
    lab <- input$active_plot_contrast %||% labs[[1]]
    rv$analysis$pairwise[[lab]] %||% rv$analysis$pairwise[[1]]
  })

  result_plot_fields <- function(x) {
    req(rv$analysis)
    if (identical(rv$analysis$engine %||% "", "limma")) {
      data.frame(
        padj = suppressWarnings(as.numeric(x$adj.P.Val)),
        pvalue = suppressWarnings(as.numeric(x$P.Value)),
        lfc = suppressWarnings(as.numeric(x$logFC)),
        meanExpression = suppressWarnings(as.numeric(x$AveExpr)),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        padj = suppressWarnings(as.numeric(x$padj)),
        pvalue = suppressWarnings(as.numeric(x$pvalue)),
        lfc = suppressWarnings(as.numeric(x$log2FoldChange)),
        meanExpression = log2(suppressWarnings(as.numeric(x$baseMean)) + 1),
        stringsAsFactors = FALSE
      )
    }
  }

  output$volcano_plot <- renderPlot({
    x <- active_pair_result()
    req(x)
    alpha <- input$plot_alpha %||% 0.05
    lfc_thr <- input$plot_lfc %||% 0
    d <- result_plot_fields(x)
    d$neglog10padj <- -log10(pmax(d$padj, .Machine$double.xmin))
    d$Significance <- "Not significant"
    d$Significance[!is.na(d$padj) & d$padj < alpha & !is.na(d$lfc) & d$lfc >= lfc_thr] <- "Up"
    d$Significance[!is.na(d$padj) & d$padj < alpha & !is.na(d$lfc) & d$lfc <= -lfc_thr] <- "Down"
    ggplot(d, aes(lfc, neglog10padj, shape = Significance)) +
      geom_point(alpha = 0.58, size = 1.7) +
      geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
      geom_hline(yintercept = -log10(alpha), linetype = 2) +
      labs(title = paste("Volcano plot —", input$active_plot_contrast %||% names(rv$analysis$pairwise)[[1]]),
           x = "log2 fold change", y = "-log10 adjusted P-value") +
      theme_minimal(base_size = 12)
  })

  output$md_plot <- renderPlot({
    x <- active_pair_result()
    req(x)
    lfc_thr <- input$plot_lfc %||% 0
    d <- result_plot_fields(x)
    ggplot(d, aes(meanExpression, lfc)) +
      geom_point(alpha = 0.58, size = 1.7) +
      geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
      labs(title = paste("Mean-difference plot —", input$active_plot_contrast %||% names(rv$analysis$pairwise)[[1]]),
           x = "Average expression", y = "log2 fold change") +
      theme_minimal(base_size = 12)
  })

  output$padj_hist <- renderPlot({
    req(rv$analysis)
    v <- if (identical(rv$analysis$engine %||% "", "limma")) rv$analysis$result$adj.P.Val else rv$analysis$result$padj
    d <- data.frame(padj = suppressWarnings(as.numeric(v)))
    d <- d[is.finite(d$padj), , drop = FALSE]
    req(nrow(d) > 0)
    ggplot(d, aes(padj)) +
      geom_histogram(bins = 40) +
      labs(title = "Adjusted P-value histogram", x = "Adjusted P-value", y = "Features") +
      theme_minimal(base_size = 12)
  })

  output$pca_plot <- renderPlot({
    req(rv$analysis)
    if (identical(rv$analysis$engine %||% "", "limma")) {
      d <- microarray_pca_dataframe(rv$analysis$normalized, rv$analysis$sample_meta)
      ttl <- "PCA of processed microarray expression"
    } else {
      d <- pca_dataframe(rv$analysis$vst, rv$analysis$sample_meta)
      ttl <- "PCA of variance-stabilized counts"
    }
    pct <- attr(d, "pct")
    ggplot(d, aes(PC1, PC2, shape = Group)) +
      geom_point(size = 3, alpha = 0.82) +
      labs(title = ttl, x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])) +
      theme_minimal(base_size = 12)
  })

  observeEvent(input$set_profile, {
    rv$profile_id <- trimws(input$profile_gene %||% "")
  })

  profile_data <- reactive({
    req(rv$analysis, rv$profile_id)
    if (identical(rv$analysis$engine %||% "", "limma")) {
      microarray_profile_dataframe(rv$analysis$normalized, rv$analysis$sample_meta, rv$profile_id, rv$analysis$result)
    } else {
      profile_dataframe(rv$analysis$normalized, rv$analysis$sample_meta, rv$profile_id, rv$analysis$result)
    }
  })

  output$profile_status_ui <- renderUI({
    if (is.null(rv$analysis)) return(div(class = "warning-box", "Run an analysis first to create the expression matrix used by the profile graph."))
    if (is.null(rv$profile_id) || !nzchar(rv$profile_id)) return(div(class = "status-note", "Enter a gene/probe ID or symbol and click Set."))
    err <- tryCatch({ profile_data(); NULL }, error = function(e) e)
    if (inherits(err, "error")) div(class = "error-box", conditionMessage(err)) else div(class = "status-good", paste("Profile:", rv$profile_id))
  })

  output$profile_plot <- renderPlot({
    d <- profile_data()
    if (identical(rv$analysis$engine %||% "", "limma")) {
      ggplot(d, aes(Accession, Expression, fill = Group)) +
        geom_col() + coord_flip() +
        labs(title = paste("Expression profile —", rv$profile_id), x = "Sample", y = "Processed expression") +
        theme_minimal(base_size = 11)
    } else {
      ggplot(d, aes(Accession, NormalizedCount, fill = Group)) +
        geom_col() + coord_flip() +
        labs(title = paste("Normalized expression profile —", rv$profile_id), x = "Sample", y = "DESeq2 normalized count") +
        theme_minimal(base_size = 11)
    }
  })

  output$profile_table <- renderDT({
    d <- profile_data()
    datatable(d, rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
  })

  current_r_script <- reactive({
    if (!is.null(rv$analysis) && identical(rv$analysis$engine %||% "", "limma")) {
      return(generate_geo2r_limma_script(
        gse_id = rv$gse_id %||% "",
        platform = selected_platform_label(),
        groups = rv$groups,
        assignments = analysis_assignments(),
        p_adjust_label = input$p_adjust_label %||% "Benjamini & Hochberg (False discovery rate)",
        log_transform = input$log_transform %||% "Auto-detect",
        force_normalization = identical(input$force_normalization %||% "No", "Yes"),
        vooma = identical(input$vooma_option %||% "No", "Yes"),
        test_name = rv$analysis$test_name %||% NULL
      ))
    }
    rnaseq_mode <- input$rnaseq_analysis_mode %||% "geo2r"
    multi_platform_rna <- length(selected_platforms()) > 1L && identical(selected_assay(), "RNA-seq")
    parity_mode <- identical(rnaseq_mode, "geo2r") && !multi_platform_rna
    script_covariates <- unique(c(if (multi_platform_rna) "Platform" else character(),
                                  if (identical(rnaseq_mode, "custom")) (input$covariates %||% character()) else character()))
    generate_geo2r_r_script(
      gse_id = rv$gse_id %||% "",
      organism = input$organism %||% "",
      groups = rv$groups,
      assignments = analysis_assignments(),
      p_adjust_label = input$p_adjust_label %||% "Benjamini & Hochberg (False discovery rate)",
      alpha = input$plot_alpha %||% 0.05,
      lfc_threshold = input$plot_lfc %||% 0,
      min_count = if (identical(rnaseq_mode, "geo2r")) 10L else (input$min_count %||% 10L),
      min_samples = if (identical(rnaseq_mode, "geo2r")) NULL else (input$min_samples %||% 2L),
      covariates = script_covariates,
      prefilter_mode = rnaseq_mode,
      test_name = rv$analysis$test_name %||% NULL,
      count_source = rv$count_obj$source %||% NULL
    )
  })

  output$r_script_text <- renderText({ paste(current_r_script(), collapse = "\n") })

  output$download_sample_excel <- downloadHandler(
    filename = function() paste0(rv$gse_id %||% "GSE", "_", sanitize_design_name(input$organism %||% "organism"), "_sample_metadata.xlsx"),
    content = function(file) {
      req(rv$metadata)
      all_samples <- excel_safe_df(full_sample_table())
      displayed_samples <- excel_safe_df(sample_table())
      filtered_rows <- input$samples_table_rows_all
      if (is.null(filtered_rows)) filtered_rows <- seq_len(nrow(displayed_samples))
      filtered_rows <- filtered_rows[filtered_rows >= 1 & filtered_rows <= nrow(displayed_samples)]
      selected_rows <- input$samples_table_rows_selected %||% integer()
      selected_rows <- selected_rows[selected_rows >= 1 & selected_rows <= nrow(displayed_samples)]

      wb <- openxlsx::createWorkbook(creator = "GEO2R — Multi-Species")
      add_excel_sheet(wb, "Samples", all_samples, "TableStyleMedium4")
      if (length(filtered_rows) < nrow(displayed_samples)) {
        add_excel_sheet(wb, "Filtered Samples", displayed_samples[filtered_rows, , drop = FALSE], "TableStyleMedium9")
      }
      if (length(selected_rows)) add_excel_sheet(wb, "Selected Samples", displayed_samples[selected_rows, , drop = FALSE], "TableStyleMedium6")
      info <- data.frame(
        Field = c("GSE accession", "Series title", "Organism", "Platform(s)", "Samples", "Defined groups", "Exported at"),
        Value = c(rv$gse_id %||% "", extract_series_title(rv$gse_list[[1]]), input$organism %||% "", selected_platform_label(),
                  nrow(all_samples), paste(rv$groups, collapse = "; "), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        stringsAsFactors = FALSE
      )
      add_excel_sheet(wb, "Series Info", info, "TableStyleLight9")
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

  output$download_result_tsv <- downloadHandler(
    filename = function() paste0(rv$gse_id %||% "GSE", "_GEO2R_full_results.tsv"),
    content = function(file) {
      req(rv$analysis)
      cols <- intersect(rv$analysis$export_columns %||% names(rv$analysis$result), names(rv$analysis$result))
      readr::write_tsv(rv$analysis$result[, cols, drop = FALSE], file, na = "")
    }
  )

  output$download_result_excel <- downloadHandler(
    filename = function() paste0(rv$gse_id %||% "GSE", "_GEO2R_results.xlsx"),
    content = function(file) {
      req(rv$analysis, rv$metadata)
      wb <- openxlsx::createWorkbook(creator = "GEO2R — Multi-Species")
      cols <- intersect(rv$analysis$export_columns %||% names(rv$analysis$result), names(rv$analysis$result))
      add_excel_sheet(wb, "GEO2R Results", rv$analysis$result[, cols, drop = FALSE], "TableStyleMedium2")

      if (length(rv$analysis$pairwise)) {
        for (lab in names(rv$analysis$pairwise)) {
          pc <- intersect(rv$analysis$export_columns %||% names(rv$analysis$pairwise[[lab]]), names(rv$analysis$pairwise[[lab]]))
          add_excel_sheet(wb, paste0("Contrast_", lab), rv$analysis$pairwise[[lab]][, pc, drop = FALSE], "TableStyleMedium9")
        }
      }

      sample_export <- rv$metadata
      sample_export$Group <- unname(analysis_assignments()[sample_export$Accession])
      sample_export$UsedInAnalysis <- sample_export$Accession %in% rv$analysis$sample_meta$Accession
      add_excel_sheet(wb, "Samples", sample_export, "TableStyleMedium4")

      if (identical(rv$analysis$engine %||% "", "limma")) {
        setup <- data.frame(
          Setting = c("GSE", "Organism", "Platform(s)", "Assay", "Engine", "Groups", "Test",
                      "P-value adjustment", "Log transform option", "Log transform applied", "Force normalization applied",
                      "vooma applied", "Significance cutoff", "Plot log2FC threshold", "Features tested", "Samples analyzed"),
          Value = c(rv$gse_id, input$organism %||% "", selected_platform_label(), "Microarray", "limma",
                    paste(rv$analysis$groups, collapse = "; "), rv$analysis$test_name,
                    input$p_adjust_label %||% "", input$log_transform %||% "Auto-detect", rv$analysis$log_transform_applied,
                    rv$analysis$normalization_applied, rv$analysis$vooma_applied,
                    input$plot_alpha %||% 0.05, input$plot_lfc %||% 0,
                    nrow(rv$analysis$result), nrow(rv$analysis$sample_meta)),
          stringsAsFactors = FALSE
        )
        matrix_sheet <- "Processed Expression"
      } else {
        setup <- data.frame(
          Setting = c("GSE", "Organism", "Platform(s)", "Assay", "Engine", "DESeq2 version", "Count source", "Groups", "Group sizes", "DESeq2 test", "Full design", "Reduced design",
                      "RNA-seq mode", "Low-count pre-filter", "Size-factor estimator", "P-value adjustment", "Significance/independent-filter alpha", "Independent-filter mean threshold", "Non-NA raw P-values", "Non-NA adjusted P-values", "Plot log2FC threshold",
                      "Covariates", "Minimum count", "Minimum samples", "Features before pre-filter", "Features after pre-filter", "Samples analyzed"),
          Value = c(rv$gse_id, input$organism %||% "", selected_platform_label(), "RNA-seq", "DESeq2", rv$analysis$deseq2_version %||% "", rv$count_obj$source %||% "",
                    paste(rv$analysis$groups, collapse = "; "), paste(paste(names(rv$analysis$group_sizes), as.integer(rv$analysis$group_sizes), sep = "="), collapse = "; "),
                    rv$analysis$test_name, rv$analysis$design, rv$analysis$reduced_design,
                    if (isTRUE(rv$analysis$geo2r_parity)) "GEO2R parity" else "Extended/custom", rv$analysis$filter_description,
                    rv$analysis$size_factor_type, input$p_adjust_label %||% "", input$plot_alpha %||% 0.05, rv$analysis$independent_filter_threshold,
                    rv$analysis$pvalue_non_na, rv$analysis$padj_non_na, input$plot_lfc %||% 0,
                    paste(rv$analysis$covariates, collapse = "; "), rv$analysis$min_count, rv$analysis$min_samples,
                    rv$analysis$features_before_prefilter, rv$analysis$features_after_prefilter, nrow(rv$analysis$sample_meta)),
          stringsAsFactors = FALSE
        )
        matrix_sheet <- "Normalized Counts"
      }
      add_excel_sheet(wb, "Options", setup, "TableStyleLight9")

      if (identical(rv$analysis$engine %||% "", "DESeq2") && !is.null(rv$annotation_status)) {
        aqc <- rv$annotation_status$coverage
        aqc$CoveragePercent <- round(100 * aqc$Coverage, 2)
        aqc$Source <- paste(rv$annotation_status$recovered_from %||% character(), collapse = "; ")
        add_excel_sheet(wb, "Annotation QC", aqc, "TableStyleMedium6")
      }

      norm <- as.data.frame(rv$analysis$normalized, check.names = FALSE)
      norm <- data.frame(ID = rownames(norm), norm, check.names = FALSE, stringsAsFactors = FALSE)
      if (nrow(norm) <= 100000 && ncol(norm) <= 500) add_excel_sheet(wb, matrix_sheet, norm, "TableStyleLight1")

      si <- capture.output(sessionInfo())
      openxlsx::addWorksheet(wb, "Session Info")
      openxlsx::writeData(wb, "Session Info", data.frame(SessionInfo = si, stringsAsFactors = FALSE))
      openxlsx::setColWidths(wb, "Session Info", 1, 100)
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

  output$download_r_script <- downloadHandler(
    filename = function() paste0(rv$gse_id %||% "GSE", "_GEO2R_script.R"),
    content = function(file) writeLines(current_r_script(), file)
  )

}

shinyApp(ui, server)
