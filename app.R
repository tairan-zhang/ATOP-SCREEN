# SPDX-License-Identifier: GPL-3.0-or-later
# app.R
# Part of ATOP-SCREEN
#
# Copyright (C) 2026 Tairan Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

app_version <- trimws(readLines("VERSION", warn = FALSE)[1])

options(shiny.maxRequestSize = 100 * 1024^2)

suppressPackageStartupMessages({
  library(shiny)
  library(readxl)
  library(dplyr)
  library(DT)
  library(clusterProfiler)
  library(shinyjs)
  library(ggplot2)
  library(ggpubr)
  library(ggnewscale)
  library(stringr)
  library(enrichplot)
  library(patchwork)
  library(tidyr)
  library(zip)
  library(RColorBrewer)
  library(parallel)
  library(data.table)
  library(ggrepel)
})

function_files <- c(
  "R/interface/cpp_permutation_interface.R",
  "R/interface/cpp_progress_wrapper.R",
  "R/analysis/permutation_functions.R",
  "R/analysis/engine_selector.R",
  "R/analysis/crispr_analysis_functions.R",
  "R/analysis/mageck_wrappers.R",
  "R/plotting/gsea_functions.R",
  "R/plotting/plotting_functions.R",
  "R/server/server_functions.R"
)
for (file in function_files) source(file)
if (!initialize_cpp_engine()) message("C++ engine unavailable; R fallback will be used.")

ui <- fluidPage(
  title = paste0("ATOP-SCREEN-", app_version),
  shinyjs::useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", href = "app.css"),
    tags$script(src = "app.js"),
    tags$link(rel = "icon", type = "image/png", href = "atop-logo.png"),
    tags$script(HTML("
      function updateColorInputStyle(inputId) {
        var inputElement = document.getElementById(inputId);
        if (!inputElement) return;

        var colorValue = inputElement.value.trim();
        var textColor = 'black';

        if (/^#[0-9A-Fa-f]{6}$/.test(colorValue) || /^#[0-9A-Fa-f]{3}$/.test(colorValue)) {
          var hex = colorValue.replace('#', '');
          if (hex.length === 3) {
            hex = hex.split('').map(function(char) { return char + char; }).join('');
          }
          var r = parseInt(hex.substring(0, 2), 16);
          var g = parseInt(hex.substring(2, 4), 16);
          var b = parseInt(hex.substring(4, 6), 16);
          var luminance = (0.299 * r + 0.587 * g + 0.114 * b);
          textColor = luminance < 140 ? 'white' : 'black';
          inputElement.style.backgroundColor = colorValue;
          inputElement.style.color = textColor;
        } else {
          inputElement.style.backgroundColor = '';
          inputElement.style.color = '';
        }
      }
    "))
  ),
  div(id = "run-overlay", class = "run-overlay", hidden = NA, tabindex = "-1", role = "status", `aria-live` = "polite",
    div(class = "run-progress",
      h2(id = "run-title", "Running analysis"),
      div(id = "run-percent", class = "run-percent", "0%"),
      tags$progress(id = "run-progress-bar", max = 100, value = 0, `aria-label` = "Analysis progress"),
      p(id = "run-step", "Preparing your data…"))),
  tags$a(class = "skip-link", href = "#workspace", "Skip to results"),
  tags$header(class = "app-header",
    div(class = "brand", tags$img(src = "atop-logo.png", alt = "ATOP logo", class = "brand-mark", width = 44, height = 44),
      div(tags$strong("ATOP-SCREEN"))),
    tags$span(paste0("Version ", app_version), class = "version-badge")
  ),
  div(
    class = "app-container", id = "workbench-layout",
    div(class = "workspace-navigation",
      tabsetPanel(id = "workspace_area", type = "pills",
        tabPanel("Analysis", value = "analysis"),
        tabPanel("Results", value = "results"),
        tabPanel("Visualization", value = "visualization")
      ),
      conditionalPanel("input.workspace_area == 'analysis'",
        tabsetPanel(id = "analysis_section",
          tabPanel("Screen analysis", value = "screen"), tabPanel("GSEA analysis", value = "gsea"))) ,
      conditionalPanel("input.workspace_area == 'results'",
        tabsetPanel(id = "result_level",
          tabPanel("sgRNA", value = "sgrna"), tabPanel("Gene", value = "gene"), tabPanel("Pathway", value = "pathway"))),
      conditionalPanel("input.workspace_area == 'visualization'",
        tabsetPanel(id = "visualization_level",
          tabPanel("sgRNA", value = "sgrna"), tabPanel("Gene", value = "gene"), tabPanel("Pathway", value = "pathway")),
        div(class = "plot-type-menu",
          tags$span("Plot type", class = "menu-label"),
          conditionalPanel("input.visualization_level == 'sgrna'",
            tabsetPanel(id = "sgrna_plot_menu", tabPanel("Paired plot", value = "paired"))),
          conditionalPanel("input.visualization_level == 'gene'",
            tabsetPanel(id = "gene_vis_plot_type", tabPanel("Volcano plot", value = "volcano"), tabPanel("Gene score ranking", value = "ranking"))),
          conditionalPanel("input.visualization_level == 'pathway'",
            tabsetPanel(id = "pathway_plot_menu", tabPanel("Lollipop plot", value = "lollipop"), tabPanel("Enrichment curve", value = "enrichment")))
        )
      )
    ),
    div(id = "inline-message", class = "inline-message", role = "alert", hidden = NA,
      tags$span(id = "inline-message-text"), tags$button("Dismiss", type = "button", id = "dismiss-message", class = "btn btn-default")),
    sidebarLayout(
      sidebarPanel(
        conditionalPanel("input.workspace_area == 'visualization' && input.visualization_level == 'pathway' && input.pathway_plot_menu == 'lollipop'", uiOutput("gsea_lollipop_params_ui_placeholder")),
        conditionalPanel("input.workspace_area == 'visualization' && input.visualization_level == 'pathway' && input.pathway_plot_menu == 'enrichment'", uiOutput("gsea_enrichment_plot_params_ui_placeholder")),
        conditionalPanel("input.workspace_area == 'visualization' && input.visualization_level == 'sgrna'", uiOutput("sgrna_paired_plot_params_ui_placeholder")),
        conditionalPanel("input.workspace_area == 'visualization' && input.visualization_level == 'gene'", uiOutput("gene_vis_params_ui_placeholder")),
        width = 3
      ),
      mainPanel(
        tags$div(id = "workspace", tabindex = "-1"),
        tabsetPanel(
          id = "main_results_tabs", type = "hidden",
          tabPanel(
            "Results", value = "Data Processing Output",
            div(class = "analysis-page",
        conditionalPanel("input.workspace_area == 'analysis' && input.analysis_section == 'screen'",
        div(
          class = "section-box", id = "upload_file_section",
          panel_heading("Screen analysis setup", "Upload a count table, then map identifiers and comparison columns.", "INPUT DATA"),
          upload_input("raw_file_upload", "Choose a screen count table", "CSV, TSV, TXT or Excel · Up to 100 MB", c(".csv", ".txt", ".tsv", ".xlsx")),
          uiOutput("uploaded_file_summary"),
          uiOutput("column_definition_ui"),
          div(class = "method-field",
          selectInput("analysis_engine", "Analysis method",
            choices = c(
              "ATOP-CRISPR (Default)" = "atop",
              "MAGeCK RRA" = "mageck_rra",
              "MAGeCK MLE" = "mageck_mle"
            ),
            selected = "atop"
          )
          )
        )),

        conditionalPanel("input.workspace_area == 'analysis' && input.analysis_section == 'screen' && output.has_screen_input",
        div(
          class = "section-box", id = "data_processing_params_box",
          panel_heading("Review and run", "Adjust parameters if needed, then run the analysis.", "ANALYSIS"),
          settings_group("Analysis parameters",
          conditionalPanel(
            condition = "input.analysis_engine == 'atop'",
            div(
              class = "settings-note",
              h4("Adaptive Top-N Aggregation Algorithm", style = "color: #333333; margin-top: 0;"),
              tags$ul(
                tags$li("Genes with sgRNA count < threshold are excluded"),
                tags$li("For remaining genes: Use Top-k mean, where k = ceil(2n/3)"),
                style = "color: #494949; margin: 5px 0;"
              )
            ),
            numericInput("pseudo_count_lfc_input", "Pseudo-count for LFC calculation:", value = 1.0, min = 0, step = 0.1),
            uiOutput("min_sgrna_threshold_ui"),
            div(
              class = "settings-section",
              h4("Permutation Test Configuration", style = "color: #333333; margin-top: 0;"),

              numericInput("N_perm_input", "Number of Permutations:", value = 10000, min = 0, step = 100),
              helpText("Recommended: ≥ 10000 permutations for stable p-values"),

              selectInput("permutation_engine_selector", "Select Calculation Engine:",
                choices = list(
                  "C++ Engine" = "cpp",
                  "R Parallel Engine" = "r_parallel"
                ),
                selected = "cpp"
              ),

              div(
                id = "cpu_info", style = "margin-top: 10px;",
                tags$script(HTML("
                  $(document).ready(function() {
                    var cores = navigator.hardwareConcurrency || 'Unknown';
                    $('#cpu_info').html('<small style=\"color: #666;\">Detected ' + cores + ' CPU cores</small>');
                  });
                "))
              )
            )
          ),
          conditionalPanel(
            condition = "input.analysis_engine == 'mageck_rra'",
            div(
              class = "alert alert-info",
              h4("MAGeCK RRA Algorithm"),
              p("Will execute 'mageck test' command. Ensure MAGeCK is installed in your environment.")
            )
          ),
          conditionalPanel(
            condition = "input.analysis_engine == 'mageck_mle'",
            div(
              class = "alert alert-info",
              h4("MAGeCK MLE Algorithm"),
              p("Will execute 'mageck mle' command. Ensure MAGeCK is installed in your environment.")
            )
          )
          ),
          actionButton("run_data_processing",
            "Run analysis",
            class = "btn-primary btn-lg btn-block",
            style = "width: 100%; white-space: normal; height: auto; padding: 12px 16px; line-height: 1.5; font-size: 16px; box-sizing: border-box;"
          )
        ))
            )
          ),
          tabPanel(
            "Pathways", value = "GSEA Analysis Results",
            div(class = "analysis-page gsea-page",
        conditionalPanel("input.workspace_area == 'analysis' && input.analysis_section == 'gsea'", uiOutput("gsea_params_ui_placeholder"))
            )
          ),
          tabPanel("sgRNA results", value = "sgRNA Results",
            uiOutput("sgrna_results_ui")),
          tabPanel("Gene results", value = "Gene Results",
            uiOutput("gene_results_ui")),
          tabPanel("Pathway results", value = "Pathway Results",
            uiOutput("pathway_results_panel")),
          tabPanel(
            "Lollipop", value = "GSEA Lollipop Plot",
            uiOutput("gsea_lollipop_plot_main_ui_placeholder")
          ),
          tabPanel(
            "Enrichment", value = "GSEA Enrichment Plot",
            uiOutput("gsea_enrichment_plot_main_ui_placeholder")
          ),
          tabPanel(
            "sgRNA pairs", value = "sgRNA Paired Plot",
            uiOutput("sgrna_paired_plot_main_ui_placeholder")
          ),
          tabPanel(
            "Gene plots", value = "Gene Visualization",
            uiOutput("gene_vis_main_ui_placeholder")
          )
        ),
        width = 9
      )
    )
  )
)

server <- function(input, output, session) {
  raw_data_info <- reactiveVal(NULL)
  output$uploaded_file_summary <- renderUI({
    info <- raw_data_info()
    req(info)
    size_label <- if (info$size < 1024^2) sprintf("%.1f KB", info$size / 1024) else sprintf("%.2f MB", info$size / 1024^2)
    div(class = "file-summary",
      div(class = "file-summary-heading", tags$strong(info$name),
        tags$span(sprintf("%s · %d columns", size_label, length(available_columns())), class = "file-metadata")),
      div(class = "column-chips", lapply(available_columns(), function(name) tags$span(name, class = "column-chip"))))
  })
  output$has_screen_input <- reactive(!is.null(raw_data_info()))
  outputOptions(output, "has_screen_input", suspendWhenHidden = FALSE)
  data_processing_results <- reactiveVal(NULL)
  gsea_results <- reactiveVal(NULL)
  available_columns <- reactiveVal(character(0))
  lfc_plot_object <- reactiveVal(NULL)

  gsea_lollipop_plot_object <- reactiveVal(NULL)
  output$has_gsea_lollipop_plot <- reactive(!is.null(gsea_lollipop_plot_object()))
  outputOptions(output, "has_gsea_lollipop_plot", suspendWhenHidden = FALSE)
  gsea_lollipop_plot_params_for_download <- reactiveVal(NULL)

  gsea_enrichment_plot_object <- reactiveVal(NULL)
  output$has_gsea_enrichment_plot <- reactive(!is.null(gsea_enrichment_plot_object()))
  outputOptions(output, "has_gsea_enrichment_plot", suspendWhenHidden = FALSE)
  gsea_enrichment_plot_params_for_download <- reactiveVal(NULL)
  gsea_enrichment_selected_pathway_id <- reactiveVal(NULL)

  sgrna_paired_plot_object <- reactiveVal(NULL)
  output$has_sgrna_paired_plot <- reactive(!is.null(sgrna_paired_plot_object()))
  outputOptions(output, "has_sgrna_paired_plot", suspendWhenHidden = FALSE)
  sgrna_paired_plot_params_for_download <- reactiveVal(NULL)
  sgrna_plot_batch_active <- reactiveVal(FALSE)
  sgrna_plot_batch_genes <- reactiveVal(NULL)

  gene_vis_plot_object <- reactiveVal(NULL)
  output$has_gene_vis_plot <- reactive(!is.null(gene_vis_plot_object()))
  outputOptions(output, "has_gene_vis_plot", suspendWhenHidden = FALSE)
  gene_vis_plot_params_for_download <- reactiveVal(NULL)

  observeEvent(data_processing_results(), {
    gene_vis_plot_object(NULL)
    gene_vis_plot_params_for_download(NULL)
    sgrna_paired_plot_object(NULL)
    sgrna_paired_plot_params_for_download(NULL)
    sgrna_plot_batch_active(FALSE)
    sgrna_plot_batch_genes(NULL)
    gsea_results(NULL)
  }, ignoreNULL = FALSE)

  observeEvent(gsea_results(), {
    gsea_lollipop_plot_object(NULL)
    gsea_lollipop_plot_params_for_download(NULL)
    gsea_enrichment_plot_object(NULL)
    gsea_enrichment_plot_params_for_download(NULL)
    gsea_enrichment_selected_pathway_id(NULL)
  }, ignoreNULL = FALSE)

  observeEvent(input$raw_file_upload, {
    inFile <- input$raw_file_upload
    if (is.null(inFile)) {
      raw_data_info(NULL)
      available_columns(character(0))
      output$status_output_processing <- renderText({
        "No file selected."
      })
      output$column_definition_ui <- renderUI({
        NULL
      })

      return(NULL)
    }
    tryCatch(
      {
        df_preview <- switch(tools::file_ext(tolower(inFile$name)),
          "csv" = read.csv(inFile$datapath, stringsAsFactors = FALSE, check.names = FALSE, nrows = 1),
          "txt" = read.delim(inFile$datapath, stringsAsFactors = FALSE, check.names = FALSE, nrows = 1),
          "tsv" = read.delim(inFile$datapath, stringsAsFactors = FALSE, check.names = FALSE, nrows = 1),
          "xlsx" = read_excel(inFile$datapath, .name_repair = "minimal", n_max = 1),
          stop("Unsupported file type")
        )
        col_names <- names(df_preview)
        available_columns(col_names)
        raw_data_info(list(datapath = inFile$datapath, name = inFile$name, type = tools::file_ext(tolower(inFile$name)), size = inFile$size))

        output$status_output_processing <- renderText({
          "File ready. Confirm the column assignments and analysis settings below."
        })

        output$column_definition_ui <- renderUI({
          req(available_columns())
          column_names <- available_columns()

          default_grna_col <- select_col(column_names, c("sgrna", "grna", "guide", "id", "name"), ignore.case = TRUE)
          default_gene_col <- select_col(column_names, c("gene", "symbol", "target"), ignore.case = TRUE)

          default_seq_col <- ""

          replicate_cols <- setdiff(column_names, c(default_grna_col, default_gene_col, default_seq_col))

          tagList(
            div(class = "field-grid field-grid-three",
            selectInput("grna_col_selector", "gRNA/sgRNA Column:", choices = column_names, selected = default_grna_col),
            selectInput("gene_col_selector", "Gene Column:", choices = column_names, selected = default_gene_col),
            selectInput("sequence_col_selector", "Sequence Column (Optional):", choices = c("Please select sequence column" = "", column_names), selected = "")
            ),
            tags$hr(),

            div(
              h4("Select Replicate Columns"),
              p("ΔLFC will be calculated via log2(Treatment / Control).", class = "text-muted"),
              div(class = "field-grid field-grid-two",
              selectizeInput("cond1_reps_selector", "Treatment Condition (Numerator):", choices = replicate_cols, multiple = TRUE, options = list(placeholder = "Select treatment replicate columns...")),
              selectizeInput("cond2_reps_selector", "Control Condition (Denominator):", choices = replicate_cols, multiple = TRUE, options = list(placeholder = "Select control replicate columns..."))
              )
            )
          )
        })

      },
      error = function(e) {
        raw_data_info(NULL)
        available_columns(character(0))
        show_inline_message(paste("Unable to read file:", e$message), type = "error")
        output$status_output_processing <- renderText({
          paste("Failed to read file column names:", e$message)
        })
        output$column_definition_ui <- renderUI({
          p("Cannot read file column names. Please check file format or content.", style = "color:#333333;")
        })

      }
    )
  })

  max_sgrna_per_gene <- reactive({
    req(raw_data_info())
    req(input$gene_col_selector, input$grna_col_selector)

    file_info <- raw_data_info()
    df <- switch(file_info$type,
      "csv" = read.csv(file_info$datapath, stringsAsFactors = FALSE, check.names = FALSE),
      "txt" = read.delim(file_info$datapath, stringsAsFactors = FALSE, check.names = FALSE),
      "tsv" = read.delim(file_info$datapath, stringsAsFactors = FALSE, check.names = FALSE),
      "xlsx" = read_excel(file_info$datapath, .name_repair = "minimal"),
      stop("Unsupported file type")
    )

    sgrna_counts <- df %>%
      distinct(!!sym(input$gene_col_selector), !!sym(input$grna_col_selector)) %>%
      group_by(!!sym(input$gene_col_selector)) %>%
      summarise(count = n(), .groups = "drop")

    max(sgrna_counts$count, na.rm = TRUE)
  })

  output$min_sgrna_threshold_ui <- renderUI({
    req(max_sgrna_per_gene())
    max_val <- max_sgrna_per_gene()

    tagList(
      div(
        class = "settings-section",
        h4("sgRNA Filtering Threshold", style = "color: #333333; margin-top: 0;"),
        numericInput("min_sgrna_threshold",
          "Minimum sgRNA count per gene:",
          value = 3,
          min = 1,
          max = max_val,
          step = 1
        ),
        helpText(sprintf("Genes with sgRNA count < threshold will be excluded from analysis. (Max detected: %d sgRNAs/gene)", max_val))
      )
    )
  })

  observeEvent(input$run_data_processing, {
    req(
      raw_data_info(),
      input$grna_col_selector, input$gene_col_selector
    )

    cond1_cols <- input$cond1_reps_selector
    cond2_cols <- input$cond2_reps_selector

    if (length(cond1_cols) == 0 || length(cond2_cols) == 0) {
      show_inline_message("Select treatment and control replicate columns before running analysis.", type = "error")
      return()
    }

    data_processing_start_abs_val <- 0.01
    withProgress(message = "Running screen analysis", value = 0, max = if (identical(input$analysis_engine, "atop")) 0.5 else 1, session = session, {
      shiny::setProgress(value = data_processing_start_abs_val, detail = "Preparing...", message = "Data processing in progress...")
      output$status_output_processing <- renderText({
        "Data processing in progress..."
      })
      data_processing_results(NULL)
      gsea_results(NULL)
      output$gsea_results_display_ui <- renderUI({
        NULL
      })

      lfc_plot_object(NULL)

      sgrna_paired_plot_object(NULL)

      gene_vis_plot_object(NULL)

      tryCatch(
        {
          full_df <- switch(raw_data_info()$type,
            "csv" = read.csv(raw_data_info()$datapath, stringsAsFactors = FALSE, check.names = FALSE),
            "txt" = read.delim(raw_data_info()$datapath, stringsAsFactors = FALSE, check.names = FALSE),
            "tsv" = read.delim(raw_data_info()$datapath, stringsAsFactors = FALSE, check.names = FALSE),
            "xlsx" = read_excel(raw_data_info()$datapath, .name_repair = "minimal"),
            stop("Internal Error: Cannot read file")
          )
          sequence_col_to_pass <- input$sequence_col_selector
          if (sequence_col_to_pass == "") sequence_col_to_pass <- NULL

          engine_choice <- input$analysis_engine
          if (is.null(engine_choice)) engine_choice <- "atop"

          message("[Data Processing] Selected analysis engine: ", engine_choice)

          results <- NULL

          if (engine_choice == "atop") {
            user_selected_engine <- input$permutation_engine_selector
            if (is.null(user_selected_engine)) user_selected_engine <- "auto"

            results <- perform_crispr_screen_analysis(
              raw_data_df = full_df,
              gRNA_col = input$grna_col_selector,
              gene_col = input$gene_col_selector,
              sequence_col = sequence_col_to_pass,
              condition1_replicate_cols = input$cond1_reps_selector,
              condition2_replicate_cols = input$cond2_reps_selector,
              N_perm = as.integer(input$N_perm_input),
              pseudo_count_lfc = input$pseudo_count_lfc_input,
              min_sgrna_threshold = ifelse(is.null(input$min_sgrna_threshold), 3, input$min_sgrna_threshold),
              user_engine_choice = user_selected_engine,
              skip_normalization = FALSE,
              diff_score_col1 = NULL,
              diff_score_col2 = NULL,
              shiny_session = session,
              initial_progress_value_abs = data_processing_start_abs_val
            )
          } else if (engine_choice == "mageck_rra") {
            results <- perform_mageck_rra_analysis_wrapper(
              raw_data_df = full_df,
              gRNA_col = input$grna_col_selector,
              gene_col = input$gene_col_selector,
              condition1_replicate_cols = input$cond1_reps_selector,
              condition2_replicate_cols = input$cond2_reps_selector,
              shiny_session = session
            )
          } else if (engine_choice == "mageck_mle") {
            results <- perform_mageck_mle_analysis_wrapper(
              raw_data_df = full_df,
              gRNA_col = input$grna_col_selector,
              gene_col = input$gene_col_selector,
              condition1_replicate_cols = input$cond1_reps_selector,
              condition2_replicate_cols = input$cond2_reps_selector,
              shiny_session = session
            )
          }
          data_processing_results(results)
          shiny::setProgress(value = 0.5, message = "Data processing completed!", detail = "Preparing GSEA analysis options...")

          output$status_output_processing <- renderText({
            tool_info <- if (!is.null(results$params$analysis_tool)) {
              paste0(" | Analysis Tool: ", results$params$analysis_tool)
            } else {
              ""
            }
            paste0("Data processing completed. sgRNA data: ", nrow(results$processed_sg_data), " rows; Gene summary: ", nrow(results$gene_summary_data), " rows", tool_info)
          })

          navigate_workspace(session, "Gene Results")
        },
        error = function(e) {
          output$status_output_processing <- renderText({
            paste("Data processing failed:", e$message)
          })
          data_processing_results(NULL)
          show_inline_message(paste("Analysis failed:", e$message), type = "error")
          shiny::setProgress(value = 0, message = "Processing failed", detail = e$message)
        }
      )
    })
  })

  observe({
    if (!is.null(data_processing_results()) && !is.null(data_processing_results()$params$gene_col)) {

      output$status_output_gsea <- renderText({
        paste(
          "Ready to use data processing output for GSEA analysis.\n",
          "Gene count:", nrow(data_processing_results()$gene_summary_data), "\n",
          "Gene ID column:", data_processing_results()$params$gene_col
        )
      })
    } else {

      output$status_output_gsea <- renderText({
        "Please complete the data processing step first."
      })
    }
  })

  output$gsea_params_ui_placeholder <- renderUI({
    req(data_processing_results())
    render_gsea_parameter_ui(data_processing_results()$params$gene_col)
  })
  output$sgrna_paired_plot_params_ui_placeholder <- renderUI({
    results <- data_processing_results()
    req(results)
    if (is.null(gsea_results())) render_sgrna_params_ui_basic(results)
    else render_sgrna_params_ui_with_gsea(results, gsea_results())
  })
  output$gene_vis_params_ui_placeholder <- renderUI({
    req(data_processing_results())
    render_gene_vis_params_ui(data_processing_results())
  })
  output$gsea_lollipop_params_ui_placeholder <- renderUI({
    req(gsea_results())
    render_gsea_lollipop_params_ui(gsea_results())
  })
  output$gsea_enrichment_plot_params_ui_placeholder <- renderUI({
    req(gsea_results())
    render_gsea_enrichment_params_ui(gsea_results())
  })

  observe_gene_search(input, session, data_processing_results, "sgrna_paired_gene_selector_direct")
  observe_gene_search(input, session, data_processing_results, "gene_vis_volcano_labels")

  gmt_summary <- eventReactive(input$gmt_file_upload, {
    file <- input$gmt_file_upload
    req(file)
    tryCatch({
      entries <- strsplit(readLines(file$datapath, warn = FALSE), "\t", fixed = TRUE)
      entries <- entries[lengths(entries) >= 3L]
      if (!length(entries)) stop("No valid gene sets found in the GMT file.")
      genes <- lapply(entries, function(entry) unique(entry[-c(1, 2)][nzchar(entry[-c(1, 2)])]))
      sizes <- lengths(genes)
      list(name = file$name, size = file$size, sets = length(entries),
        genes = length(unique(unlist(genes, use.names = FALSE))),
        min = min(sizes), max = max(sizes), median = median(sizes))
    }, error = function(e) list(name = file$name, error = conditionMessage(e)))
  })
  output$gmt_file_summary <- renderUI({
    info <- gmt_summary()
    req(info)
    if (!is.null(info$error)) return(div(class = "file-summary", strong(info$name), p(info$error)))
    size_label <- if (info$size < 1024^2) sprintf("%.1f KB", info$size / 1024) else sprintf("%.2f MB", info$size / 1024^2)
    div(class = "file-summary",
      div(class = "file-summary-heading", strong(info$name), span(size_label, class = "file-metadata")),
      div(class = "column-chips",
        span(sprintf("%s gene sets", format(info$sets, big.mark = ",")), class = "column-chip"),
        span(sprintf("%s unique genes", format(info$genes, big.mark = ",")), class = "column-chip")),
      p(sprintf("Genes per set: %d–%d · Median: %s", info$min, info$max, info$median), class = "file-metadata"))
  })

  opt_dt_scroll <- get_dt_options()

  output$sg_data_table <- renderDT(
    {
      results <- data_processing_results()
      req(results)

      data_to_show <- if (!is.null(results$raw_mageck_sgrna_summary)) {
        results$raw_mageck_sgrna_summary
      } else {
        results$processed_sg_data
      }

      req(data_to_show)
      datatable(data_to_show, options = opt_dt_scroll, rownames = FALSE, width = "100%")
    },
    server = TRUE
  )

  output$gene_summary_table <- renderDT(
    {
      results <- data_processing_results()
      req(results)

      data_to_show <- if (!is.null(results$raw_mageck_gene_summary)) {
        results$raw_mageck_gene_summary
      } else {
        results$gene_summary_data
      }

      req(data_to_show)
      datatable(data_to_show, options = opt_dt_scroll, rownames = FALSE, width = "100%")
    },
    server = TRUE
  )

  output$filtered_genes_table <- renderDT(
    {
      req(data_processing_results()$filtered_genes_sgrna_data)
      datatable(data_processing_results()$filtered_genes_sgrna_data, options = opt_dt_scroll, rownames = FALSE, width = "100%")
    },
    server = TRUE
  )

  output$gsea_positive_table <- renderDT(
    {
      req(gsea_results()$gsea_results_positive_df)

      positive_df <- gsea_results()$gsea_results_positive_df
      if ("NES" %in% names(positive_df) && nrow(positive_df) > 0) {
        positive_df <- positive_df[order(-positive_df$NES), ]
      }
      datatable(positive_df, options = opt_dt_scroll, rownames = FALSE, width = "100%")
    },
    server = TRUE
  )

  output$gsea_negative_table <- renderDT(
    {
      req(gsea_results()$gsea_results_negative_df)

      negative_df <- gsea_results()$gsea_results_negative_df
      if ("NES" %in% names(negative_df) && nrow(negative_df) > 0) {
        negative_df <- negative_df[order(negative_df$NES), ]
      }
      datatable(negative_df, options = opt_dt_scroll, rownames = FALSE, width = "100%")
    },
    server = TRUE
  )

  output$download_sg_data <- gen_dl_handler("processed_sg_data", "sgRNA_data", data_processing_results, raw_mageck_key = "raw_mageck_sgrna_summary")
  output$download_gene_summary <- gen_dl_handler("gene_summary_data", "gene_summary_data", data_processing_results, raw_mageck_key = "raw_mageck_gene_summary")
  output$download_filtered_genes <- gen_dl_handler("filtered_genes_sgrna_data", "filtered_genes_sgrna", data_processing_results)

  output$download_gsea_pos_csv <- downloadHandler(
    filename = function() {
      paste0("gsea_positive_results_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(gsea_results()$gsea_results_positive_df)
      positive_df <- gsea_results()$gsea_results_positive_df
      if ("NES" %in% names(positive_df) && nrow(positive_df) > 0) {
        positive_df <- positive_df[order(-positive_df$NES), ]
      }
      write.csv(positive_df, file, row.names = FALSE, na = "")
    }
  )

  output$download_gsea_neg_csv <- downloadHandler(
    filename = function() {
      paste0("gsea_negative_results_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(gsea_results()$gsea_results_negative_df)
      negative_df <- gsea_results()$gsea_results_negative_df
      if ("NES" %in% names(negative_df) && nrow(negative_df) > 0) {
        negative_df <- negative_df[order(negative_df$NES), ]
      }
      write.csv(negative_df, file, row.names = FALSE, na = "")
    }
  )
  output$download_gsea_pos_rds <- gen_gsea_rds_dl_handler("gsea_object_positive", "gsea_positive_object", gsea_results)
  output$download_gsea_neg_rds <- gen_gsea_rds_dl_handler("gsea_object_negative", "gsea_negative_object", gsea_results)

  observe_gsea_analysis(input, output, session, data_processing_results, gsea_results)

  output$status_output_processing <- renderText({
    "Please upload data and define columns to start."
  })
  output$status_output_gsea <- renderText({
    "Please select GSEA data source and complete settings."
  })
  output$gsea_lollipop_plot_status <- renderText({
    "Select pathways and configure the plot settings in the sidebar."
  })
  output$sgrna_paired_plot_status <- renderText({
    "Select genes and comparison columns in the sidebar."
  })
  output$gene_vis_plot_status <- renderText({
    "Choose a gene plot type and configure the settings in the sidebar."
  })

  plot_dimensions <- list(
    gene_vis = register_plot_preview(input, output, session, "gene_vis", "gene_vis_plot_render", gene_vis_plot_object, 12, 10),
    gsea_lollipop = register_plot_preview(input, output, session, "gsea_lollipop", "gsea_lollipop_plot_render", gsea_lollipop_plot_object, 18, 14),
    gsea_enrichment = register_plot_preview(input, output, session, "gsea_enrichment", "gsea_enrichment_plot_render", gsea_enrichment_plot_object, 16, 12),
    sgrna_paired = register_plot_preview(input, output, session, "sgrna_paired", "sgrna_paired_plot_render", sgrna_paired_plot_object, 8.5, 10.5)
  )

  output$gene_vis_main_ui_placeholder <- renderUI({
    if (is.null(data_processing_results())) return(empty_state("Explore gene hits", "Run Screen analysis, then choose a plot type above and genes in the sidebar."))
    tagList(
      div(
        class = "section-box",
        h2("Gene Visualization Plot"),
        plot_preview("gene_vis_plot_render", "has_gene_vis_plot", "700px"),
        verbatimTextOutput("gene_vis_plot_status")
      ),
      shinyjs::hidden(
        div(
          id = "gene_vis_download_options", class = "section-box",
          hr(),
          h3("Preview and download"),
          fluidRow(
            column(3, selectInput("gene_vis_download_format", "Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
            column(3, numericInput("gene_vis_download_width", "Width (cm):", value = isolate(plot_dimensions$gene_vis()$width), min = 3, max = 40, step = 0.5)),
            column(3, numericInput("gene_vis_download_height", "Height (cm):", value = isolate(plot_dimensions$gene_vis()$height), min = 3, max = 40, step = 0.5)),
            column(3, conditionalPanel(
              condition = "input.gene_vis_download_format == 'png'",
              numericInput("gene_vis_download_dpi", "DPI (PNG):", value = isolate(plot_dimensions$gene_vis()$dpi), min = 72, max = 600, step = 50)
            ))
          ),
          actionButton("gene_vis_apply_preview", "Apply", class = "btn-default"),
          textOutput("gene_vis_preview_dimensions"),
          downloadButton("download_gene_vis_plot", "Download Plot", class = "btn-primary")
        )
      )
    )
  })

  output$gsea_lollipop_plot_main_ui_placeholder <- renderUI({
    if (is.null(gsea_results())) return(empty_state("Compare enriched pathways", "Complete GSEA analysis, then choose pathways to visualize."))
    tagList(
      div(
        class = "section-box",
        h2("GSEA Pathway Enrichment Ranked Plot"),
        plot_preview("gsea_lollipop_plot_render", "has_gsea_lollipop_plot", "700px"),
        verbatimTextOutput("gsea_lollipop_plot_status")
      ),
      shinyjs::hidden(
        div(
          id = "gsea_lollipop_plot_download_options", class = "section-box",
          hr(),
          h3("Preview and download"),
          fluidRow(
            column(3, selectInput("gsea_lollipop_download_format", "Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
            column(3, numericInput("gsea_lollipop_download_width", "Width (cm):", value = isolate(plot_dimensions$gsea_lollipop()$width), min = 3, max = 40, step = 0.5)),
            column(3, numericInput("gsea_lollipop_download_height", "Height (cm):", value = isolate(plot_dimensions$gsea_lollipop()$height), min = 3, max = 40, step = 0.5)),
            column(3, conditionalPanel(
              condition = "input.gsea_lollipop_download_format == 'png'",
              numericInput("gsea_lollipop_download_dpi", "DPI (PNG):", value = isolate(plot_dimensions$gsea_lollipop()$dpi), min = 72, max = 600, step = 50)
            ))
          ),
          actionButton("gsea_lollipop_apply_preview", "Apply", class = "btn-default"),
          textOutput("gsea_lollipop_preview_dimensions"),
          downloadButton("download_gsea_lollipop_plot", "Download Plot", class = "btn-primary")
        )
      )
    )
  })

  output$gsea_enrichment_plot_main_ui_placeholder <- renderUI({
    if (is.null(gsea_results())) return(empty_state("Inspect pathway enrichment", "Complete GSEA analysis, then select a pathway to inspect its enrichment curve."))
    tagList(
      div(
        class = "section-box",
        h2("GSEA Pathway Enrichment Plot"),
        plot_preview("gsea_enrichment_plot_render", "has_gsea_enrichment_plot", "700px"),
        verbatimTextOutput("gsea_enrichment_plot_status")
      ),
      shinyjs::hidden(
        div(
          id = "gsea_enrichment_plot_download_options", class = "section-box",
          hr(),
          h3("Preview and download"),
          fluidRow(
            column(3, selectInput("gsea_enrichment_download_format", "Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
            column(3, numericInput("gsea_enrichment_download_width", "Width (cm):", value = isolate(plot_dimensions$gsea_enrichment()$width), min = 3, max = 40, step = 0.5)),
            column(3, numericInput("gsea_enrichment_download_height", "Height (cm):", value = isolate(plot_dimensions$gsea_enrichment()$height), min = 3, max = 40, step = 0.5)),
            column(3, conditionalPanel(
              condition = "input.gsea_enrichment_download_format == 'png'",
              numericInput("gsea_enrichment_download_dpi", "DPI (PNG):", value = isolate(plot_dimensions$gsea_enrichment()$dpi), min = 72, max = 600, step = 50)
            ))
          ),
          actionButton("gsea_enrichment_apply_preview", "Apply", class = "btn-default"),
          textOutput("gsea_enrichment_preview_dimensions"),
          downloadButton("download_gsea_enrichment_plot", "Download Pathway Plot", class = "btn-primary")
        )
      )
    )
  })

  output$sgrna_paired_plot_main_ui_placeholder <- renderUI({
    if (is.null(data_processing_results())) return(empty_state("Compare guide counts", "Run Screen analysis to compare treatment and control for individual genes."))
    if (isTRUE(sgrna_plot_batch_active())) {
      tagList(
        div(
          class = "section-box",
          h2("Batch Download sgRNA Paired Plots"),
          p(paste("Ready to generate and download paired plots for", length(sgrna_plot_batch_genes()), "core genes in the selected pathway.")),
          hr(),
          h4("Batch Download Options"),
          fluidRow(
            column(3, selectInput("sgrna_batch_download_format", "Graph Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
            column(3, numericInput("sgrna_batch_download_width", "Width (cm):", value = 8.5, min = 3, max = 40, step = 0.5)),
            column(3, numericInput("sgrna_batch_download_height", "Height (cm):", value = 10.5, min = 3, max = 40, step = 0.5)),
            column(3, conditionalPanel(
              condition = "input.sgrna_batch_download_format == 'png'",
              numericInput("sgrna_batch_download_dpi", "DPI (PNG):", value = 300, min = 72, max = 600, step = 50)
            ))
          ),
          downloadButton("download_sgrna_all_paired_plots_zip", "Download All Paired Plots (.zip)", class = "btn-primary btn-lg"),
          tags$hr(),
          verbatimTextOutput("sgrna_paired_plot_batch_status")
        )
      )
    } else {
      tagList(
        div(
          class = "section-box",
          h2("sgRNA Paired Plot (Single Gene Preview)"),
          plot_preview("sgrna_paired_plot_render", "has_sgrna_paired_plot", "420px"),
          verbatimTextOutput("sgrna_paired_plot_status")
        ),
        shinyjs::hidden(
          div(
            id = "sgrna_paired_plot_download_options", class = "section-box",
            hr(),
            h3("Preview and download"),
            fluidRow(
              column(3, selectInput("sgrna_paired_download_format", "Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
              column(3, numericInput("sgrna_paired_download_width", "Width (cm):", value = isolate(plot_dimensions$sgrna_paired()$width), min = 3, max = 40, step = 0.5)),
              column(3, numericInput("sgrna_paired_download_height", "Height (cm):", value = isolate(plot_dimensions$sgrna_paired()$height), min = 3, max = 40, step = 0.5)),
              column(3, conditionalPanel(
                condition = "input.sgrna_paired_download_format == 'png'",
                numericInput("sgrna_paired_download_dpi", "DPI (PNG):", value = isolate(plot_dimensions$sgrna_paired()$dpi), min = 72, max = 600, step = 50)
              ))
            ),
            actionButton("sgrna_paired_apply_preview", "Apply", class = "btn-default"),
          textOutput("sgrna_paired_preview_dimensions"),
          downloadButton("download_sgrna_paired_plot", "Download Paired Plot", class = "btn-primary")
          )
        )
      )
    }
  })

  observe_gsea_lollipop_plot_generation(input, output, session, gsea_results, gsea_lollipop_plot_object, gsea_lollipop_plot_params_for_download)

  observe_gsea_enrichment_plot_generation(input, output, session, gsea_results, gsea_enrichment_plot_object, gsea_enrichment_plot_params_for_download, gsea_enrichment_selected_pathway_id)

  observe_gsea_enrichment_pathway_selector_update(input, session, gsea_results)
  observe_sgrna_gsea_type_selector_update(input, session, gsea_results)
  observe_sgrna_pathway_selector_update(input, session, gsea_results, data_processing_results)

  output$download_gsea_lollipop_plot <- create_gsea_lollipop_download_handler(input, gsea_lollipop_plot_object, gsea_lollipop_plot_params_for_download, plot_dimensions$gsea_lollipop)
  output$download_gsea_enrichment_plot <- create_gsea_enrichment_download_handler(input, gsea_enrichment_plot_object, gsea_enrichment_plot_params_for_download, plot_dimensions$gsea_enrichment)

  observe_sgrna_paired_plot_generation(
    input, output, session, data_processing_results, gsea_results,
    sgrna_paired_plot_object, sgrna_paired_plot_params_for_download,
    sgrna_plot_batch_genes, sgrna_plot_batch_active
  )

  output$download_sgrna_paired_plot <- create_sgrna_paired_plot_download_handler(input, sgrna_paired_plot_object, sgrna_paired_plot_params_for_download, plot_dimensions$sgrna_paired)

  output$download_sgrna_all_paired_plots_zip <- create_sgrna_batch_download_handler(input, output, sgrna_plot_batch_genes, data_processing_results, gsea_results)

  observe_gene_vis_plot_generation(input, output, session, data_processing_results, gene_vis_plot_object, gene_vis_plot_params_for_download)
  output$download_gene_vis_plot <- create_gene_vis_download_handler(input, gene_vis_plot_object, gene_vis_plot_params_for_download, plot_dimensions$gene_vis)

  output$gsea_enrichment_plot_status <- renderText({
    "Select a pathway and configure the enrichment curve in the sidebar."
  })
  output$sgrna_paired_plot_batch_status <- renderText({
    "Batch download status will be displayed here."
  })

  output$sgrna_results_ui <- renderUI({
    if (is.null(data_processing_results())) return(empty_state("Explore guide-level results", "Run Screen analysis to inspect guide scores and download your sgRNA data."))
    div(class = "section-box", panel_heading("sgRNA results", "Guide-level scores and measurements", "RESULTS"),
      render_data_proc_tables_ui(data_processing_results(), "sgrna"))
  })
  output$gene_results_ui <- renderUI({
    if (is.null(data_processing_results())) return(empty_state("Discover gene-level hits", "Run Screen analysis to explore gene scores, significance and guide coverage."))
    div(class = "section-box", panel_heading("Gene results", "Aggregated scores and statistical significance", "RESULTS"),
      render_data_proc_tables_ui(data_processing_results(), "gene"))
  })
  output$pathway_results_panel <- renderUI({
    if (is.null(gsea_results())) return(empty_state("Explore enriched pathways", "Complete GSEA analysis to compare pathway enrichment and download the results."))
    div(class = "section-box", panel_heading("Pathway results", "Enrichment across positive and negative gene scores", "RESULTS"),
      uiOutput("gsea_results_display_ui"))
  })

  active_page <- reactive({
    area <- input$workspace_area
    if (is.null(area) || area == "analysis") {
      if (identical(input$analysis_section, "gsea")) "GSEA Analysis Results" else "Data Processing Output"
    } else if (area == "results") {
      switch(if (is.null(input$result_level)) "sgrna" else input$result_level, sgrna = "sgRNA Results", gene = "Gene Results", pathway = "Pathway Results")
    } else {
      switch(if (is.null(input$visualization_level)) "sgrna" else input$visualization_level,
        sgrna = "sgRNA Paired Plot", gene = "Gene Visualization",
        pathway = if (identical(input$pathway_plot_menu, "enrichment")) "GSEA Enrichment Plot" else "GSEA Lollipop Plot")
    }
  })
  observe({
    updateTabsetPanel(session, "main_results_tabs", selected = active_page())
    page <- active_page()
    no_parameters <- identical(input$workspace_area, "analysis") || page %in% c("sgRNA Results", "Gene Results", "Pathway Results") ||
      (page %in% c("GSEA Analysis Results", "sgRNA Paired Plot", "Gene Visualization") && is.null(data_processing_results())) ||
      (page %in% c("GSEA Lollipop Plot", "GSEA Enrichment Plot") && is.null(gsea_results()))
    shinyjs::toggleClass("workbench-layout", "results-layout", condition = no_parameters)
  })
}

shinyApp(ui = ui, server = server)
