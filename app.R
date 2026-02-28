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


# app.R
# ATOP-SCREEN Analysis App
# Main Shiny Application File limit: 100MB
options(shiny.maxRequestSize = 100 * 1024^2)

# Suppress package startup messages
suppressPackageStartupMessages({
  library(shiny)
  library(readxl)
  library(dplyr)
  library(DT)
  library(clusterProfiler)
  library(shinyjs)
  library(ggplot2)
  library(ggpubr)
  library(ggpmisc)
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

# Load Functions
cat("🚀 Loading external function modules...\n")

function_files <- c(
  "R/server/progress_tracker_functions.R",
  "R/analysis/permutation_functions.R",
  "R/analysis/engine_selector.R",
  "R/interface/cpp_progress_wrapper.R",
  "R/analysis/crispr_analysis_functions.R",
  "R/plotting/gsea_functions.R",
  "R/plotting/plotting_functions.R",
  "R/server/server_functions.R",
  "R/interface/cpp_permutation_interface.R",
  "R/analysis/mageck_wrappers.R"
)

for (file in function_files) {
  if (file.exists(file)) {
    source(file)
    cat(paste("✅", file, "loaded\n"))
  } else {
    cat(paste("⚠️", file, "does not exist\n"))
  }
}

# --- Initial C++ Engine ---
# Try Init
cat("Loading high-performance engine interface...\n")
if (exists("initialize_cpp_engine")) {
  # Try to initialize, but don't crash if it fails (it might compile on the fly)
  tryCatch(
    {
      initialize_cpp_engine()
    },
    error = function(e) {
      cat("Engine initialization delayed: ", e$message, "\n")
    }
  )
}

# UI (User Interface)
ui <- fluidPage(
  shinyjs::useShinyjs(), # Initialize shinyjs
  tags$head(
    tags$style(HTML("
      body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; margin: 0; padding: 0; background-color: #f5f5f7; color: #1d1d1f; }
      .header { background-color: rgba(0,0,0,0.8); color: white; padding: 18px 40px; text-align: left; font-size: 20px; font-weight: 600; backdrop-filter: saturate(180%) blur(20px); -webkit-backdrop-filter: saturate(180%) blur(20px); border-bottom: 1px solid rgba(255,255,255,0.1); }
      .container { padding: 30px 40px; max-width: 1400px; margin: auto; }
      .section-box { background-color: #fff; padding: 25px; border-radius: 12px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 30px; }
      h2 { font-size: 24px; font-weight: 600; margin-top: 0; margin-bottom: 20px; color: #1d1d1f; }
      h3 { font-size: 18px; font-weight: 500; margin-top: 15px; margin-bottom: 10px; color: #1d1d1f; }
      h4 { font-size: 16px; font-weight: 500; margin-top: 10px; margin-bottom: 8px; color: #1d1d1f; }
      .btn-primary { background-color: #007aff; color: white; border: none; padding: 12px 22px; border-radius: 8px; font-size: 16px; font-weight: 500; cursor: pointer; transition: background-color 0.2s ease-in-out; }
      .btn-primary:hover { background-color: #005ec4; }
      .shiny-input-container { margin-bottom: 15px; }
      label { font-weight: 500; margin-bottom: 5px; display: block; }
      .form-group.shiny-input-container label {font-weight: 500;}
      .text-muted { color: #6e6e73; font-size: 0.9em; }
      .shiny-notification-message { margin-bottom: 5px !important; }
      .shiny-notification-progress { margin-top: 0px !important; }
      .flex-container { display: flex; flex-wrap: wrap; gap: 20px; }
      .flex-item { flex: 1 1 300px; min-width: 250px; }
    ")),
    tags$script(HTML("
      function updateColorInputStyle(inputId) {
        var inputElement = document.getElementById(inputId);
        if (!inputElement) return;

        var colorValue = inputElement.value.trim();
        var textColor = 'black'; // Default text color

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
  div(class = "header", "ATOP-SCREEN"),
  div(
    class = "container",
    sidebarLayout(
      sidebarPanel(
        div(
          class = "section-box", id = "upload_file_section",
          h2("1. Upload Files"),
          fileInput("raw_file_upload", "Upload CRISPR Screen Raw Data (.csv, .txt, .tsv, .xlsx)",
            accept = c(".csv", ".txt", ".tsv", ".xlsx", "text/csv", "text/comma-separated-values,text/plain", "text/tab-separated-values", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
            buttonLabel = "Browse...", placeholder = "No file selected"
          ),
          uiOutput("column_definition_ui")
        ),
        div(
          class = "section-box", id = "data_processing_params_box",
          h2("2. Advanced Data Processing Parameters"),
          conditionalPanel(
            condition = "input.analysis_engine == 'atop'",
            div(
              style = "background-color: #e8f4fd; padding: 12px; border-radius: 6px; margin: 10px 0;",
              h4("Adaptive Top-N Aggregation Algorithm", style = "color: #0066cc; margin-top: 0;"),
              tags$ul(
                tags$li("Genes with sgRNA count < threshold are excluded"),
                tags$li("For remaining genes: Use Top-k mean, where k = ceil(2n/3)"),
                style = "color: #495057; margin: 5px 0;"
              )
            ),
            numericInput("pseudo_count_lfc_input", "Pseudo-count for LFC calculation:", value = 1.0, min = 0, step = 0.1),
            uiOutput("min_sgrna_threshold_ui"),
            div(
              style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin: 10px 0;",
              h4("Permutation Test Configuration", style = "color: #0066cc; margin-top: 0;"),

              # Permutation parameters
              numericInput("N_perm_input", "Number of Permutations:", value = 10000, min = 0, step = 100),
              helpText("Recommended: ≥ 10000 permutations for stable p-values"),

              # Engine selection
              selectInput("permutation_engine_selector", "Select Calculation Engine:",
                choices = list(
                  "C++ Engine" = "cpp",
                  "R Parallel Engine" = "r_parallel"
                ),
                selected = "cpp"
              ),

              # Dynamic engine explanation
              div(id = "engine_explanation", style = "margin-top: 10px;"),

              # CPU core detection
              div(
                id = "cpu_info", style = "margin-top: 10px;",
                tags$script(HTML("
                  $(document).ready(function() {
                    var cores = navigator.hardwareConcurrency || 'Unknown';
                    $('#cpu_info').html('<small style=\"color: #666;\">Detected ' + cores + ' CPU cores</small>');
                  });
                "))
              ),

              # Engine explanation script
              tags$script(HTML("
                $(document).on('change', '#permutation_engine_selector', function() {
                  var engine = $(this).val();
                  var explanations = {
                    'cpp': '',
                    'r_parallel': '',

                    'standard': '<div style=\"padding: 10px; background-color: #f8f9fa; border-radius: 5px; font-size: 0.9em;\"><strong>R Standard Algorithm</strong><br/>• <strong>Performance:</strong> Baseline (1x)<br/>• <strong>Features:</strong> Single-threaded, simple and reliable<br/>• <strong>Applicable:</strong> Small datasets (&lt;100 permutations)<br/>• <strong>Requirements:</strong> No special requirements<br/>• <strong>Pros:</strong> Maximum compatibility, easy debugging<br/><em>Priority for compatibility and stability</em></div>'
                  };

                  $('#engine_explanation').html(explanations[engine] || '');
                });

                // Initial display
                $(document).ready(function() {
                  $('#permutation_engine_selector').trigger('change');
                });
              "))
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
          ),
          actionButton("run_data_processing",
            "Start Data Processing",
            class = "btn-primary btn-lg btn-block",
            icon = icon("cogs"),
            style = "width: 100%; white-space: normal; height: auto; padding: 12px 16px; line-height: 1.5; font-size: 16px; box-sizing: border-box;"
          )
        ),
        uiOutput("gsea_params_ui_placeholder"),
        uiOutput("gsea_lollipop_params_ui_placeholder"),
        uiOutput("gsea_enrichment_plot_params_ui_placeholder"),
        uiOutput("sgrna_paired_plot_params_ui_placeholder"),
        uiOutput("gene_vis_params_ui_placeholder"),
        width = 4
      ),
      mainPanel(
        tabsetPanel(
          id = "main_results_tabs",
          tabPanel(
            "Data Processing Output",
            div(
              class = "section-box", id = "method_selection_section",
              h2("Analysis Method"),
              selectInput("analysis_engine", "Select Analysis Method:",
                choices = c(
                  "ATOP-CRISPR (Default)" = "atop",
                  "MAGeCK RRA" = "mageck_rra",
                  "MAGeCK MLE" = "mageck_mle"
                ),
                selected = "atop"
              ),
              p("Choose the underlying algorithm for CRISPR screen analysis.", class = "text-muted")
            ),
            div(
              class = "section-box",
              h2("Data Processing Status and Results"),
              verbatimTextOutput("status_output_processing"),
              tags$hr(),
              uiOutput("data_processing_results_tables_ui")
            )
          ),
          tabPanel(
            "GSEA Analysis Results",
            div(
              class = "section-box",
              h2("GSEA Analysis Status and Results"),
              verbatimTextOutput("status_output_gsea"),
              tags$hr(),

              # GSEA Analysis Description
              div(
                class = "section-box", style = "background: #f8f9fa; border-left: 4px solid #007bff;",
                h3("📊 GSEA Analysis"),
                div(
                  class = "alert alert-info",
                  h4("📋 Instructions"),
                  p("GSEA analysis will use the output from the data processing steps above. Please ensure data processing is completed.")
                )
              ),
              uiOutput("gsea_results_display_ui")
            )
          ),
          tabPanel(
            "GSEA Lollipop Plot",
            uiOutput("gsea_lollipop_plot_main_ui_placeholder")
          ),
          tabPanel(
            "GSEA Enrichment Plot",
            uiOutput("gsea_enrichment_plot_main_ui_placeholder")
          ),
          tabPanel(
            "sgRNA Paired Plot",
            uiOutput("sgrna_paired_plot_main_ui_placeholder")
          ),
          tabPanel(
            "Gene Visualization",
            uiOutput("gene_vis_main_ui_placeholder")
          )
        ),
        width = 8
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Reactive values
  raw_data_info <- reactiveVal(NULL)
  data_processing_results <- reactiveVal(NULL)
  gsea_results <- reactiveVal(NULL)
  available_columns <- reactiveVal(character(0))
  lfc_plot_object <- reactiveVal(NULL) # For storing the LFC scatter plot

  # GSEA Lollipop Plot reactive values
  gsea_lollipop_plot_object <- reactiveVal(NULL)
  gsea_lollipop_plot_params_for_download <- reactiveVal(NULL)

  # GSEA Enrichment Plot reactive values
  gsea_enrichment_plot_object <- reactiveVal(NULL)
  gsea_enrichment_plot_params_for_download <- reactiveVal(NULL)
  gsea_enrichment_selected_pathway_id <- reactiveVal(NULL) # Store ID of selected pathway

  # sgRNA Paired Plot reactive values
  sgrna_paired_plot_object <- reactiveVal(NULL)
  sgrna_paired_plot_params_for_download <- reactiveVal(NULL)
  sgrna_plot_batch_active <- reactiveVal(FALSE) # TRUE if batch mode for sgRNA paired plots
  sgrna_plot_batch_genes <- reactiveVal(NULL) # List of genes for batch plotting

  # Gene Visualization reactive values
  gene_vis_plot_object <- reactiveVal(NULL)
  gene_vis_plot_params_for_download <- reactiveVal(NULL)





  # Note: select_col function is defined in server_functions.R

  # Hide initial UI elements - only show upload section
  shinyjs::hide("data_processing_params_box")
  shinyjs::hide("gsea_parameters_section")
  shinyjs::hide("gsea_lollipop_params_box")
  shinyjs::hide("gsea_enrichment_plot_params_box")
  shinyjs::hide("gsea_enrichment_plot_params_box")
  shinyjs::hide("sgrna_paired_plot_params_box")
  shinyjs::hide("gene_vis_params_box")

  # --- File Upload and Column Definition ---
  observeEvent(input$raw_file_upload, {
    inFile <- input$raw_file_upload
    if (is.null(inFile)) {
      raw_data_info(NULL)
      available_columns(character(0))
      output$status_output_processing <- renderPrint({
        "No file selected."
      })
      output$column_definition_ui <- renderUI({
        NULL
      })
      # Hide all parameter sections
      shinyjs::hide("data_processing_params_box")
      shinyjs::hide("gsea_parameters_section")
      shinyjs::hide("gsea_lollipop_params_box")
      shinyjs::hide("gsea_lollipop_params_box")
      shinyjs::hide("sgrna_paired_plot_params_box")
      shinyjs::hide("gene_vis_params_box")
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
        raw_data_info(list(datapath = inFile$datapath, name = inFile$name, type = tools::file_ext(tolower(inFile$name))))

        output$status_output_processing <- renderPrint({
          cat(
            "File uploaded successfully:", inFile$name, "(", round(inFile$size / 1024^2, 2), "MB)\n",
            "Detected column names:", paste(col_names, collapse = ", ")
          )
        })

        # Dynamically generate UI for column selection based on uploaded file
        output$column_definition_ui <- renderUI({
          req(available_columns())
          column_names <- available_columns()

          # Smart column selection
          default_grna_col <- select_col(column_names, c("sgrna", "grna", "guide", "id", "name"), ignore.case = TRUE)
          default_gene_col <- select_col(column_names, c("gene", "symbol", "target"), ignore.case = TRUE)
          # Sequence column defaults to empty, requiring user selection
          default_seq_col <- ""

          # Replicate columns (exclude ID columns)
          replicate_cols <- setdiff(column_names, c(default_grna_col, default_gene_col, default_seq_col))


          tagList(
            selectInput("grna_col_selector", "gRNA/sgRNA Column:", choices = column_names, selected = default_grna_col),
            selectInput("gene_col_selector", "Gene Column:", choices = column_names, selected = default_gene_col),
            selectInput("sequence_col_selector", "Sequence Column (Optional):", choices = c("Please select sequence column" = "", column_names), selected = ""),
            tags$hr(),

            # Skip Normalization Mode Toggle
            checkboxInput("run_full_pipeline", "Run full ATOP pipeline (Normalization + LFC calculation)", value = TRUE),
            helpText("Uncheck to skip normalization and directly use pre-calculated diff_score columns"),

            # Conditional panel for FULL pipeline (replicate selection)
            conditionalPanel(
              condition = "input.run_full_pipeline == true",
              h4("Select Replicate Columns"),
              p("ΔLFC will be calculated via log2(Treatment / Control).", class = "text-muted"),
              selectizeInput("cond1_reps_selector", "Treatment Condition (Numerator):", choices = replicate_cols, multiple = TRUE, options = list(placeholder = "Select treatment replicate columns...")),
              selectizeInput("cond2_reps_selector", "Control Condition (Denominator):", choices = replicate_cols, multiple = TRUE, options = list(placeholder = "Select control replicate columns..."))
            ),

            # Conditional panel for SKIP mode (diff_score column selection)
            conditionalPanel(
              condition = "input.run_full_pipeline == false",
              div(
                style = "background-color: #fff3cd; padding: 12px; border-radius: 6px; margin: 10px 0; border-left: 4px solid #ffc107;",
                h4("Skip Normalization Mode", style = "color: #856404; margin-top: 0;"),
                p("Your data must contain two diff_score columns. The ddiff_score will be calculated as:", style = "margin-bottom: 5px;"),
                tags$code("ddiff_score = diff_score_col1 - diff_score_col2", style = "background: #f8f9fa; padding: 2px 6px; border-radius: 3px;")
              ),
              selectInput("diff_score_col1_selector", "First diff_score Column:", choices = column_names),
              selectInput("diff_score_col2_selector", "Second diff_score Column:", choices = column_names)
            )
          )
        })
        shinyjs::show("data_processing_params_box")
        shinyjs::hide("gsea_parameters_section") # Hide GSEA params until data processing is done
        shinyjs::hide("gsea_lollipop_params_box") # Also hide GSEA lollipop params initially
        shinyjs::hide("gsea_lollipop_params_box") # Also hide GSEA lollipop params initially
        shinyjs::hide("sgrna_paired_plot_params_box")
        shinyjs::hide("gene_vis_params_box")
      },
      error = function(e) {
        raw_data_info(NULL)
        available_columns(character(0))
        output$status_output_processing <- renderPrint({
          paste("Failed to read file column names:", e$message)
        })
        output$column_definition_ui <- renderUI({
          p("Cannot read file column names. Please check file format or content.", style = "color:red;")
        })
        shinyjs::hide("data_processing_params_box")
        shinyjs::hide("gsea_parameters_section")
        shinyjs::hide("gsea_lollipop_params_box")
        shinyjs::hide("gsea_lollipop_params_box")
        shinyjs::hide("sgrna_paired_plot_params_box")
        shinyjs::hide("gene_vis_params_box")
      }
    )
  })

  # Reactive for detecting max sgRNA count per gene
  max_sgrna_per_gene <- reactive({
    req(raw_data_info())
    req(input$gene_col_selector, input$grna_col_selector)

    # Read the data from file
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


  # Render dynamic UI for min sgRNA threshold
  output$min_sgrna_threshold_ui <- renderUI({
    req(max_sgrna_per_gene())
    max_val <- max_sgrna_per_gene()

    tagList(
      div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin: 10px 0;",
        h4("sgRNA Filtering Threshold", style = "color: #0066cc; margin-top: 0;"),
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

  # Data processing observer
  observeEvent(input$run_data_processing, {
    req(
      raw_data_info(),
      input$grna_col_selector, input$gene_col_selector
    )

    cond1_cols <- input$cond1_reps_selector
    cond2_cols <- input$cond2_reps_selector

    # Check if skip normalization mode is enabled
    skip_mode <- !isTRUE(input$run_full_pipeline)

    # Validate inputs based on mode
    if (skip_mode) {
      # Skip mode: validate diff_score columns
      if (is.null(input$diff_score_col1_selector) || is.null(input$diff_score_col2_selector)) {
        showModal(modalDialog(title = "Error", "Both diff_score columns must be selected in skip normalization mode."))
        return()
      }
      if (input$diff_score_col1_selector == input$diff_score_col2_selector) {
        showModal(modalDialog(title = "Error", "Please select two different diff_score columns."))
        return()
      }
    } else {
      # Full pipeline mode: validate replicate columns
      if (length(cond1_cols) == 0 || length(cond2_cols) == 0) {
        showModal(modalDialog(title = "Error", "At least one replicate column must be selected for all condition groups."))
        return()
      }
    }

    data_processing_start_abs_val <- 0.01
    withProgress(message = "Data processing in progress...", value = 0, session = session, {
      shiny::setProgress(value = data_processing_start_abs_val, detail = "Preparing...", message = "Data processing in progress...")
      output$status_output_processing <- renderPrint({
        "Data processing in progress..."
      })
      data_processing_results(NULL)
      gsea_results(NULL)
      output$data_processing_results_tables_ui <- renderUI({
        NULL
      })
      output$gsea_results_display_ui <- renderUI({
        NULL
      })
      shinyjs::hide("gsea_parameters_section") # Ensure GSEA section is hidden during processing
      output$gsea_params_ui_placeholder <- renderUI({
        NULL
      }) # Clear GSEA params placeholder
      lfc_plot_object(NULL) # Clear previous LFC plot
      shinyjs::hide("sgrna_paired_plot_params_box") # Hide sgRNA paired plot params during (re)processing
      output$sgrna_paired_plot_params_ui_placeholder <- renderUI({
        NULL
      }) # Clear sgRNA paired plot params placeholder
      sgrna_paired_plot_object(NULL) # Clear previous sgRNA paired plot

      shinyjs::hide("gene_vis_params_box")
      output$gene_vis_params_ui_placeholder <- renderUI({
        NULL
      })
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

          # Determine analysis engine
          engine_choice <- input$analysis_engine
          if (is.null(engine_choice)) engine_choice <- "atop"

          # Debug: Show which engine is selected
          message("[Data Processing] Selected analysis engine: ", engine_choice)

          results <- NULL

          if (engine_choice == "atop") {
            # Handle user engine selection (Permutation)
            user_selected_engine <- input$permutation_engine_selector
            if (is.null(user_selected_engine)) user_selected_engine <- "auto"

            # Call core analysis function (ATOP)
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
              skip_normalization = !isTRUE(input$run_full_pipeline),
              diff_score_col1 = if (!isTRUE(input$run_full_pipeline)) input$diff_score_col1_selector else NULL,
              diff_score_col2 = if (!isTRUE(input$run_full_pipeline)) input$diff_score_col2_selector else NULL,
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
          data_processing_results(results) # Store results, including params now
          shiny::setProgress(value = 0.5, message = "Data processing completed!", detail = "Preparing GSEA analysis options...")

          output$status_output_processing <- renderPrint({
            tool_info <- if (!is.null(results$params$analysis_tool)) {
              paste0(" | Analysis Tool: ", results$params$analysis_tool)
            } else {
              ""
            }
            paste0("Data processing completed. sgRNA data: ", nrow(results$processed_sg_data), " rows; Gene summary: ", nrow(results$gene_summary_data), " rows", tool_info)
          })

          # Render data processing results tables
          output$data_processing_results_tables_ui <- renderUI({
            render_data_proc_tables_ui(results)
          })

          # Render GSEA params UI
          output$gsea_params_ui_placeholder <- renderUI({
            render_gsea_parameter_ui(selected_gene_col_for_gsea = results$params$gene_col)
          })
          shinyjs::show("gsea_parameters_section")

          # Render LFC scatter plot params UI
          # Render sgRNA paired plot params UI (basic version, no GSEA results needed)
          output$sgrna_paired_plot_params_ui_placeholder <- renderUI({
            results <- data_processing_results()
            req(results)
            req(results$normalized_counts, results$gene_summary_data, results$params$gene_col, results$params$gRNA_col)

            render_sgrna_params_ui_basic(results)
          })

          # Render Gene Visualization Params UI
          output$gene_vis_params_ui_placeholder <- renderUI({
            render_gene_vis_params_ui(results)
          })
          shinyjs::show("gene_vis_params_box")

          updateTabsetPanel(session, "main_results_tabs", selected = "GSEA Analysis Results")
        },
        error = function(e) {
          output$status_output_processing <- renderPrint({
            paste("Data processing failed:", e$message)
          })
          data_processing_results(NULL)
          shiny::setProgress(value = 0, message = "Processing failed", detail = e$message)
        }
      )
    })
  })


  # --- GSEA Params Initialization ---
  observe({
    # Use data processing parameters
    if (!is.null(data_processing_results()) && !is.null(data_processing_results()$params$gene_col)) {
      output$gsea_params_ui_placeholder <- renderUI({
        render_gsea_parameter_ui(selected_gene_col_for_gsea = data_processing_results()$params$gene_col)
      })
      shinyjs::show("gsea_parameters_section")
      output$status_output_gsea <- renderPrint({
        paste(
          "Ready to use data processing output for GSEA analysis.\n",
          "Gene count:", nrow(data_processing_results()$gene_summary_data), "\n",
          "Gene ID column:", data_processing_results()$params$gene_col
        )
      })
    } else {
      output$gsea_params_ui_placeholder <- renderUI({
        NULL
      })
      output$status_output_gsea <- renderPrint({
        "Please complete the data processing step first."
      })
    }
  })



  # --- Data Table Rendering (Unified Configuration) ---
  opt_dt_scroll <- get_dt_options()

  # Data Processing Results Table
  output$sg_data_table <- renderDT(
    {
      results <- data_processing_results()
      req(results)

      # Use raw MAGeCK output if available, otherwise use processed_sg_data
      data_to_show <- if (!is.null(results$raw_mageck_sgrna_summary)) {
        results$raw_mageck_sgrna_summary
      } else {
        results$processed_sg_data
      }

      req(data_to_show)
      datatable(data_to_show, options = opt_dt_scroll, rownames = FALSE)
    },
    server = TRUE
  )

  output$gene_summary_table <- renderDT(
    {
      results <- data_processing_results()
      req(results)

      # Use raw MAGeCK output if available, otherwise use gene_summary_data
      data_to_show <- if (!is.null(results$raw_mageck_gene_summary)) {
        results$raw_mageck_gene_summary
      } else {
        results$gene_summary_data
      }

      req(data_to_show)
      datatable(data_to_show, options = opt_dt_scroll, rownames = FALSE)
    },
    server = TRUE
  )

  output$filtered_genes_table <- renderDT(
    {
      req(data_processing_results()$filtered_genes_sgrna_data)
      datatable(data_processing_results()$filtered_genes_sgrna_data, options = opt_dt_scroll, rownames = FALSE)
    },
    server = TRUE
  )


  # GSEA Results Table
  output$gsea_positive_table <- renderDT(
    {
      req(gsea_results()$gsea_results_positive_df)
      # Positive results sorted by NES descending
      positive_df <- gsea_results()$gsea_results_positive_df
      if ("NES" %in% names(positive_df) && nrow(positive_df) > 0) {
        positive_df <- positive_df[order(-positive_df$NES), ]
      }
      datatable(positive_df, options = opt_dt_scroll, rownames = FALSE)
    },
    server = TRUE
  )

  output$gsea_negative_table <- renderDT(
    {
      req(gsea_results()$gsea_results_negative_df)
      # Negative results sorted by NES ascending
      negative_df <- gsea_results()$gsea_results_negative_df
      if ("NES" %in% names(negative_df) && nrow(negative_df) > 0) {
        negative_df <- negative_df[order(negative_df$NES), ]
      }
      datatable(negative_df, options = opt_dt_scroll, rownames = FALSE)
    },
    server = TRUE
  )

  # --- Download Handlers (Data Processing Results) ---
  output$download_sg_data <- gen_dl_handler("processed_sg_data", "sgRNA_data", data_processing_results, raw_mageck_key = "raw_mageck_sgrna_summary")
  output$download_gene_summary <- gen_dl_handler("gene_summary_data", "gene_summary_data", data_processing_results, raw_mageck_key = "raw_mageck_gene_summary")
  output$download_filtered_genes <- gen_dl_handler("filtered_genes_sgrna_data", "filtered_genes_sgrna", data_processing_results)


  # --- Download Handlers (GSEA Results) ---
  # Positive results download: sorted by NES descending
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

  # Negative results download: sorted by NES ascending
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

  # --- Main Analysis Observers ---
  observe_gsea_analysis(input, output, session, data_processing_results, gsea_results)
  observe_tab_switching(input, session)

  # --- Initial Status Messages ---
  output$status_output_processing <- renderPrint({
    "Please upload data and define columns to start."
  })
  output$status_output_gsea <- renderPrint({
    "Please select GSEA data source and complete settings."
  })
  output$gsea_lollipop_plot_status <- renderPrint({
    "Please complete GSEA analysis first, then select parameters in the sidebar to generate the GSEA Lollipop Plot."
  })
  output$sgrna_paired_plot_status <- renderPrint({
    "Please complete data processing first, then select parameters in the sidebar to generate the sgRNA Paired Plot."
  })
  output$gene_vis_plot_status <- renderPrint({
    "Please complete data processing first, then select parameters in the sidebar to generate the Gene Visualization Plot."
  })

  # --- Plot Rendering ---
  # GSEA Lollipop Plot Render
  output$gsea_lollipop_plot_render <- renderPlot({
    req(gsea_lollipop_plot_object())
    gsea_lollipop_plot_object()
  })

  # GSEA Enrichment Plot Render
  output$gsea_enrichment_plot_render <- renderPlot({
    req(gsea_enrichment_plot_object())
    gsea_enrichment_plot_object()
  })

  # sgRNA Paired Plot Render
  output$sgrna_paired_plot_render <- renderPlot({
    req(sgrna_paired_plot_object())
    sgrna_paired_plot_object()
  })

  # Gene Visualization Plot Render
  output$gene_vis_plot_render <- renderPlot({
    req(gene_vis_plot_object())
    gene_vis_plot_object()
  })

  output$gene_vis_main_ui_placeholder <- renderUI({
    tagList(
      div(
        class = "section-box",
        h2("Gene Visualization Plot"),
        plotOutput("gene_vis_plot_render", height = "700px", width = "100%"),
        verbatimTextOutput("gene_vis_plot_status")
      ),
      shinyjs::hidden(
        div(
          id = "gene_vis_download_options", class = "section-box",
          hr(),
          h3("Download Plot"),
          fluidRow(
            column(3, selectInput("gene_vis_download_format", "Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
            column(3, numericInput("gene_vis_download_width", "Width (inches):", value = 8, min = 3, max = 20, step = 0.5)),
            column(3, numericInput("gene_vis_download_height", "Height (inches):", value = 6, min = 3, max = 20, step = 0.5)),
            column(3, conditionalPanel(
              condition = "input.gene_vis_download_format == 'png'",
              numericInput("gene_vis_download_dpi", "DPI (PNG):", value = 300, min = 72, max = 600, step = 50)
            ))
          ),
          downloadButton("download_gene_vis_plot", "Download Plot", class = "btn-primary")
        )
      )
    )
  })

  # --- Main Plot UI Placeholders ---
  output$gsea_lollipop_plot_main_ui_placeholder <- renderUI({
    tagList(
      div(
        class = "section-box",
        h2("GSEA Pathway Enrichment Ranked Plot"),
        plotOutput("gsea_lollipop_plot_render", height = "700px", width = "100%"),
        verbatimTextOutput("gsea_lollipop_plot_status")
      ),
      shinyjs::hidden(
        div(
          id = "gsea_lollipop_plot_download_options", class = "section-box",
          hr(),
          h3("Download GSEA Lollipop Plot"),
          fluidRow(
            column(3, selectInput("gsea_lollipop_download_format", "Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
            column(3, numericInput("gsea_lollipop_download_width", "Width (inches):", value = 10, min = 3, max = 30, step = 0.5)),
            column(3, numericInput("gsea_lollipop_download_height", "Height (inches):", value = 8, min = 3, max = 30, step = 0.5)),
            column(3, conditionalPanel(
              condition = "input.gsea_lollipop_download_format == 'png'",
              numericInput("gsea_lollipop_download_dpi", "DPI (PNG):", value = 300, min = 72, max = 600, step = 50)
            ))
          ),
          downloadButton("download_gsea_lollipop_plot", "Download Plot", class = "btn-primary")
        )
      )
    )
  })

  output$gsea_enrichment_plot_main_ui_placeholder <- renderUI({
    tagList(
      div(
        class = "section-box",
        h2("GSEA Pathway Enrichment Plot"),
        plotOutput("gsea_enrichment_plot_render", height = "700px", width = "100%"),
        verbatimTextOutput("gsea_enrichment_plot_status")
      ),
      shinyjs::hidden(
        div(
          id = "gsea_enrichment_plot_download_options", class = "section-box",
          hr(),
          h3("Download GSEA Pathway Enrichment Plot"),
          fluidRow(
            column(3, selectInput("gsea_enrichment_download_format", "Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
            column(3, numericInput("gsea_enrichment_download_width", "Width (inches):", value = 8.5, min = 3, max = 20, step = 0.5)),
            column(3, numericInput("gsea_enrichment_download_height", "Height (inches):", value = 6.5, min = 3, max = 20, step = 0.5)),
            column(3, conditionalPanel(
              condition = "input.gsea_enrichment_download_format == 'png'",
              numericInput("gsea_enrichment_download_dpi", "DPI (PNG):", value = 300, min = 72, max = 600, step = 50)
            ))
          ),
          downloadButton("download_gsea_enrichment_plot", "Download Pathway Plot", class = "btn-primary")
        )
      )
    )
  })

  output$sgrna_paired_plot_main_ui_placeholder <- renderUI({
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
            column(3, numericInput("sgrna_batch_download_width", "Width (inches):", value = 6, min = 2, max = 20, step = 0.5)),
            column(3, numericInput("sgrna_batch_download_height", "Height (inches):", value = 5, min = 2, max = 20, step = 0.5)),
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
          plotOutput("sgrna_paired_plot_render", height = "600px", width = "100%"),
          verbatimTextOutput("sgrna_paired_plot_status")
        ),
        shinyjs::hidden(
          div(
            id = "sgrna_paired_plot_download_options", class = "section-box",
            hr(),
            h3("Download sgRNA Paired Plot"),
            fluidRow(
              column(3, selectInput("sgrna_paired_download_format", "Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"))),
              column(3, numericInput("sgrna_paired_download_width", "Width (inches):", value = 6, min = 2, max = 20, step = 0.5)),
              column(3, numericInput("sgrna_paired_download_height", "Height (inches):", value = 5, min = 2, max = 20, step = 0.5)),
              column(3, conditionalPanel(
                condition = "input.sgrna_paired_download_format == 'png'",
                numericInput("sgrna_paired_download_dpi", "DPI (PNG):", value = 300, min = 72, max = 600, step = 50)
              ))
            ),
            downloadButton("download_sgrna_paired_plot", "Download Paired Plot", class = "btn-primary")
          )
        )
      )
    }
  })

  # --- Observer Logic Calls ---

  # GSEA Lollipop Plot Observer - Call server_functions.R
  observe_gsea_lollipop_plot_generation(input, output, session, gsea_results, gsea_lollipop_plot_object, gsea_lollipop_plot_params_for_download)

  # GSEA Enrichment Plot Observer - Call server_functions.R
  observe_gsea_enrichment_plot_generation(input, output, session, gsea_results, gsea_enrichment_plot_object, gsea_enrichment_plot_params_for_download, gsea_enrichment_selected_pathway_id)

  # Selector Update Observer - Call server_functions.R
  observe_gsea_enrichment_pathway_selector_update(input, session, gsea_results)
  observe_sgrna_gsea_type_selector_update(input, session, gsea_results)
  observe_sgrna_pathway_selector_update(input, session, gsea_results, data_processing_results)

  # Server-Side Selectize Updaters - Call server_functions.R
  observe_sgrna_gene_selector_update(input, session, data_processing_results)
  observe_gene_vis_labels_update(input, session, data_processing_results)

  # Download Handlers - Call server_functions.R
  output$download_gsea_lollipop_plot <- create_gsea_lollipop_download_handler(input, gsea_lollipop_plot_object, gsea_lollipop_plot_params_for_download)
  output$download_gsea_enrichment_plot <- create_gsea_enrichment_download_handler(input, gsea_enrichment_plot_object, gsea_enrichment_plot_params_for_download)

  # sgRNA Paired Plot Observer - Call server_functions.R
  observe_sgrna_paired_plot_generation(
    input, output, session, data_processing_results, gsea_results,
    sgrna_paired_plot_object, sgrna_paired_plot_params_for_download,
    sgrna_plot_batch_genes, sgrna_plot_batch_active
  )

  # sgRNA Paired Plot Download Handler - Calls function from server_functions.R
  output$download_sgrna_paired_plot <- create_sgrna_paired_plot_download_handler(input, sgrna_paired_plot_object, sgrna_paired_plot_params_for_download)

  # sgRNA Paired Plot Batch Download Handler - Calls function from server_functions.R
  output$download_sgrna_all_paired_plots_zip <- create_sgrna_batch_download_handler(input, output, sgrna_plot_batch_genes, data_processing_results, gsea_results)

  # Gene Visualization Logic
  observe_gene_vis_plot_generation(input, output, session, data_processing_results, gene_vis_plot_object, gene_vis_plot_params_for_download)
  output$download_gene_vis_plot <- create_gene_vis_download_handler(input, gene_vis_plot_object, gene_vis_plot_params_for_download)



  # Initialize Status Messages
  output$gsea_enrichment_plot_status <- renderPrint({
    "Please complete GSEA analysis first, then select a pathway in the sidebar to plot."
  })
  output$sgrna_paired_plot_batch_status <- renderPrint({
    "Batch download status will be displayed here."
  })

  # --- Tab Switching Observer: Dynamically show corresponding parameter box based on current tab ---
  all_param_sections_ids <- c(
    "upload_file_section",
    "data_processing_params_box",
    "gsea_parameters_section", # GSEA Params Section ID
    "gsea_lollipop_params_box", # GSEA Lollipop Plot Params Section ID
    "gsea_enrichment_plot_params_box", # GSEA Enrichment Plot Params Section ID
    "sgrna_paired_plot_params_box" # sgRNA Paired Plot Params Section ID
  )

  observeEvent(input$main_results_tabs,
    {
      current_tab <- input$main_results_tabs

      # First hide all sidebar parameter sections
      for (id in all_param_sections_ids) {
        shinyjs::hide(id)
      }

      # Show corresponding section based on current tab and data availability
      if (current_tab == "Data Processing Output") {
        shinyjs::show("upload_file_section")
        if (!is.null(raw_data_info())) {
          shinyjs::show("data_processing_params_box")
        }
      } else if (current_tab == "GSEA Analysis Results") {
        # Show GSEA params only if data processing results are available (prerequisite for UI rendering)
        if (!is.null(data_processing_results())) {
          shinyjs::show("gsea_parameters_section")
        }
      } else if (current_tab == "GSEA Lollipop Plot") {
        # Show GSEA Lollipop Plot params only if GSEA results are available
        if (!is.null(gsea_results())) {
          shinyjs::show("gsea_lollipop_params_box")
        }
      } else if (current_tab == "GSEA Enrichment Plot") {
        # Show GSEA Enrichment Plot params only if GSEA results are available
        if (!is.null(gsea_results())) {
          shinyjs::show("gsea_enrichment_plot_params_box")
        }
      } else if (current_tab == "sgRNA Paired Plot") {
        # Show sgRNA Paired Plot params only if data processing results are available
        if (!is.null(data_processing_results())) {
          shinyjs::show("sgrna_paired_plot_params_box")
        }
      }
    },
    ignoreNULL = FALSE
  ) # ignoreNULL = FALSE to run for initial tab on app startup
}

shinyApp(ui = ui, server = server)
