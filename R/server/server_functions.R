# SPDX-License-Identifier: GPL-3.0-or-later
# server_functions.R
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


# server_functions.R
# Shiny Server Logic Support Functions Collection

# --- Color Validation Function ---
is_valid_hex_or_name <- function(color_string) {
  if (is.null(color_string) || length(color_string) != 1 || is.na(color_string) || !is.character(color_string)) {
    return(FALSE)
  }

  # Check if it is a valid hex color
  if (grepl("^#[0-9A-Fa-f]{6}$", color_string)) {
    return(TRUE)
  }

  # Check if it is a valid color name
  tryCatch(
    {
      col2rgb(color_string)
      return(TRUE)
    },
    error = function(e) {
      return(FALSE)
    }
  )
}

# --- UI Rendering Functions ---

# Render Data Tables UI
render_data_proc_tables_ui <- function(results) {
  # Determine context
  tool <- results$params$analysis_tool
  if (is.null(tool)) tool <- "ATOP" # Default

  is_mageck <- grepl("MAGeCK", tool, ignore.case = TRUE)

  # For MAGeCK: Show raw output
  if (is_mageck && !is.null(results$raw_mageck_gene_summary)) {
    prefix <- ifelse(grepl("MLE", tool), "mle", "rra")

    tagList(
      h3(paste0("MAGeCK Gene Summary (", tool, ")")),
      p("Original MAGeCK output - columns preserved as-is", class = "text-muted"),
      DTOutput("gene_summary_table"),
      downloadButton("download_gene_summary", paste0("Download ", prefix, ".gene_summary (.csv)")),
      tags$hr(),
      h3(paste0("MAGeCK sgRNA Summary (", tool, ")")),
      p("Original MAGeCK output - columns preserved as-is", class = "text-muted"),
      DTOutput("sg_data_table"),
      downloadButton("download_sg_data", paste0("Download ", prefix, ".sgrna_summary (.csv)"))
    )
  } else {
    # ATOP: Show two main tables plus filtered genes table
    tagList(
      h3("sgRNA Level Data"),
      p("Includes diff_score, rankings, and deltaLFC per replicate", class = "text-muted"),
      DTOutput("sg_data_table"),
      downloadButton("download_sg_data", "Download sgRNA Data (.csv)"),
      tags$hr(),
      h3("Gene Level Summary Data"),
      p("Includes gene scores, ranks, P-values, and sgRNA count", class = "text-muted"),
      DTOutput("gene_summary_table"),
      downloadButton("download_gene_summary", "Download Gene Summary Data (.csv)"),
      if (!is.null(results$filtered_genes_sgrna_data) && nrow(results$filtered_genes_sgrna_data) > 0) {
        tagList(
          tags$hr(),
          h3("Filtered Genes sgRNA Data"),
          p(
            sprintf(
              "Genes excluded due to sgRNA count < threshold (%d genes removed)",
              length(unique(results$filtered_genes_sgrna_data[[results$params$gene_col]]))
            ),
            class = "text-muted"
          ),
          DTOutput("filtered_genes_table"),
          downloadButton("download_filtered_genes", "Download Filtered Genes Data (.csv)")
        )
      }
    )
  }
}

# Render GSEA Results UI
render_gsea_results_tables_ui <- function(gsea_res) {
  tagList(
    if (!is.null(gsea_res) && !is.null(gsea_res$gsea_results_positive_df) && nrow(gsea_res$gsea_results_positive_df) > 0) {
      tagList(
        h3("GSEA Enrichment Results (GenePositiveScore)"), DTOutput("gsea_positive_table"),
        downloadButton("download_gsea_pos_csv", "Download Positive Score GSEA Results (.csv)"),
        downloadButton("download_gsea_pos_rds", "Download Positive Score GSEA Object (.rds)"), tags$hr()
      )
    } else {
      p("No significant enriched pathways found based on GenePositiveScore or analysis did not produce a result table.")
    },
    if (!is.null(gsea_res) && !is.null(gsea_res$gsea_results_negative_df) && nrow(gsea_res$gsea_results_negative_df) > 0) {
      tagList(
        h3("GSEA Enrichment Results (GeneNegativeScore)"), DTOutput("gsea_negative_table"),
        downloadButton("download_gsea_neg_csv", "Download Negative Score GSEA Results (.csv)"),
        downloadButton("download_gsea_neg_rds", "Download Negative Score GSEA Object (.rds)"), tags$hr()
      )
    } else {
      p("No significant enriched pathways found based on GeneNegativeScore or analysis did not produce a result table.")
    }
  )
}

# Render GSEA Params UI
render_gsea_parameter_ui <- function(selected_gene_col_for_gsea = "Gene") {
  div(
    class = "section-box", id = "gsea_parameters_section",
    h2("2. GSEA Analysis Parameters"),
    fileInput("gmt_file_upload", "Upload Gene Set GMT File (.gmt)", accept = c(".gmt")),
    p(strong("Gene ID Column for GSEA:"), code(selected_gene_col_for_gsea)),
    tags$hr(),
    h3("GSEA Advanced Parameters (clusterProfiler)"),
    fluidRow(column(12, numericInput("gsea_clp_n_perm", "Number of Permutations:", value = 1000, min = 100, step = 100))),
    fluidRow(column(12, numericInput("gsea_clp_pval_cutoff", "P-value Cutoff:", value = 1, min = 0, max = 1, step = 0.01))),
    fluidRow(column(12, numericInput("gsea_clp_min_gs_size", "Min Gene Set Size:", value = 15, min = 1, step = 1))),
    fluidRow(column(12, numericInput("gsea_max_gs_size", "Max Gene Set Size:", value = 500, min = 1, step = 1))),
    fluidRow(column(12, checkboxInput("gsea_clp_seed", "Use Fixed Seed (123)", TRUE))),
    actionButton("run_gsea_analysis",
      "Start GSEA Analysis",
      class = "btn-primary btn-lg btn-block",
      icon = icon("chart-line"),
      style = "margin-top: 15px; width: 100%; white-space: normal; height: auto; padding: 12px 16px; line-height: 1.5; font-size: 16px; box-sizing: border-box;"
    )
  )
}


# --- DataTable Options ---
get_dt_options <- function() {
  list(
    scrollX = TRUE,
    pageLength = 10,
    lengthMenu = list(c(10, 25, 50, 100, -1), c("10", "25", "50", "100", "All")),
    autoWidth = TRUE,
    columnDefs = list(list(width = "120px", targets = "_all"))
  )
}

# --- Generic Download Handler Generator Function ---
gen_dl_handler <- function(data_key, prefix, results_reactive, raw_mageck_key = NULL) {
  downloadHandler(
    filename = function() {
      res_list <- results_reactive()

      # Determine prefix based on analysis tool for MAGeCK
      final_prefix <- prefix
      if (!is.null(res_list$params$analysis_tool) && grepl("MAGeCK", res_list$params$analysis_tool)) {
        tool_prefix <- ifelse(grepl("MLE", res_list$params$analysis_tool), "mle", "rra")

        # Update prefix for MAGeCK outputs
        if (data_key == "processed_sg_data") {
          final_prefix <- paste0(tool_prefix, ".sgrna_summary")
        } else if (data_key == "gene_summary_data") {
          final_prefix <- paste0(tool_prefix, ".gene_summary")
        }
      }

      paste0(final_prefix, "_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".csv")
    },
    content = function(f) {
      res_list <- results_reactive()

      # Prefer raw MAGeCK data if available
      data_to_export <- NULL
      if (!is.null(raw_mageck_key) && !is.null(res_list[[raw_mageck_key]])) {
        data_to_export <- res_list[[raw_mageck_key]]
      } else if (!is.null(res_list[[data_key]])) {
        data_to_export <- res_list[[data_key]]
      }

      if (!is.null(data_to_export)) {
        write.csv(data_to_export, f, row.names = FALSE, na = "")
      } else {
        write.csv(data.frame(Error = paste("No data available for", data_key)), f, row.names = FALSE)
      }
    }
  )
}

# GSEA RDS Download Handler (Enhanced: Saves S4 Object + Results DF + Gene List)
# GSEA RDS Download Handler (Enhanced: Saves S4 Object + Results DF + Gene List)
gen_gsea_rds_dl_handler <- function(object_key, prefix, results_reactive) {
  downloadHandler(
    filename = function() {
      paste0(prefix, "_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".rds")
    },
    content = function(f) {
      res_list <- results_reactive()

      # Determine which result set we are saving (Positive or Negative)
      is_positive <- grepl("positive", object_key, ignore.case = TRUE)

      if (!is.null(res_list) && !is.null(res_list[[object_key]])) {
        # Construct the enhanced list object to save
        save_obj <- list(
          object = res_list[[object_key]], # The original S4 gseaResult object
          results_df = if (is_positive) res_list$gsea_results_positive_df else res_list$gsea_results_negative_df,
          gene_list = if (is_positive) res_list$gene_list_positive else res_list$gene_list_negative
        )

        saveRDS(save_obj, f)
      } else {
        error_obj <- list(Error = paste("No GSEA object available for", object_key))
        saveRDS(error_obj, f)
      }
    }
  )
}

# --- Input Validation Functions ---
is_valid_hex <- function(hex) {
  grepl("^#[0-9A-Fa-f]{6}$", hex)
}

validate_hex_inputs <- function(hex_inputs, hex_labels, session) {
  for (i in 1:length(hex_inputs)) {
    if (!is_valid_hex(hex_inputs[i])) {
      showNotification(paste0(hex_labels[i], " is not a valid hex code (e.g., #FF0000)"),
        type = "error", duration = 5, session = session
      )
      return(FALSE)
    }
  }
  return(TRUE)
}

# --- Column Selection Function ---
select_col <- function(choices, patterns, ignore.case = TRUE) {
  for (p in patterns) {
    match_idx <- grepl(p, choices, ignore.case = ignore.case)
    if (any(match_idx)) {
      return(choices[which(match_idx)[1]])
    } # Return first match
  }
  if (length(choices) > 0) {
    return(choices[1])
  } # Default to first if no match
  return(NULL) # Or empty if no choices
}

# --- LFC Scatter Plot UI Rendering Function ---
# --- sgRNA Paired Plot UI ---
render_sgrna_params_ui_basic <- function(results) {
  req(results)
  req(results$normalized_counts, results$gene_summary_data, results$params$gene_col, results$params$gRNA_col)

  norm_counts_df <- results$normalized_counts
  gene_summary_df <- results$gene_summary_data
  app_params <- results$params

  # Get data column selection (numeric columns, exclude ID columns)
  id_cols_in_norm <- c(app_params$gRNA_col, app_params$gene_col, app_params$sequence_col)
  id_cols_in_norm <- id_cols_in_norm[!sapply(id_cols_in_norm, is.null)]

  potential_condition_cols <- setdiff(names(norm_counts_df), id_cols_in_norm)
  numeric_condition_cols <- character(0)
  if (length(potential_condition_cols) > 0) {
    is_numeric_check <- sapply(norm_counts_df[, potential_condition_cols, drop = FALSE], is.numeric)
    numeric_condition_cols <- potential_condition_cols[is_numeric_check]
  }

  if (length(numeric_condition_cols) < 2) {
    return(div(class = "section-box", id = "sgrna_paired_plot_params_box", h2("5. sgRNA Paired Plot Parameters"), p("Not enough numeric columns in normalized count data for paired plot (requires at least 2).")))
  }

  # Get all genes for selection
  all_genes <- unique(gene_summary_df[[app_params$gene_col]])

  div(
    class = "section-box", id = "sgrna_paired_plot_params_box",
    h2("5. sgRNA Paired Plot Parameters"),

    # Gene Selection Method Selector
    radioButtons("sgrna_gene_selection_method", "Select Gene Selection Method:",
      choices = list("Direct Gene Selection" = "direct", "Select via GSEA Results" = "gsea"),
      selected = "direct", inline = TRUE
    ),

    # Direct Gene Selection Interface
    conditionalPanel(
      condition = "input.sgrna_gene_selection_method == 'direct'",
      selectizeInput("sgrna_paired_gene_selector_direct", "Select Gene:",
        choices = NULL,
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = "Please select gene...")
      ),
      conditionalPanel(
        condition = "input.sgrna_paired_gene_selector_direct && input.sgrna_paired_gene_selector_direct.length > 0",
        p("Selected", tags$span(id = "direct_gene_count", "0"), "genes", style = "color: #0066cc; font-weight: bold;"),
        tags$script(HTML("
            $(document).on('change', '#sgrna_paired_gene_selector_direct', function() {
              var count = $('#sgrna_paired_gene_selector_direct').val() ? $('#sgrna_paired_gene_selector_direct').val().length : 0;
              $('#direct_gene_count').text(count);
            });
          "))
      )
    ),
    # Use server-side selectize: Update choices later
    tags$script("Shiny.addCustomMessageHandler('update_sgrna_gene_choices', function(data) {
       // Placeholder for potential future custom handlers if needed
    });"),
    tags$hr(),
    h4("2. Select Plotting Data Columns"),
    selectInput("sgrna_paired_condition1_col", "Select Condition 1 Column:", choices = numeric_condition_cols, selected = if (length(numeric_condition_cols) > 0) numeric_condition_cols[1] else NULL),
    selectInput("sgrna_paired_condition2_col", "Select Condition 2 Column:", choices = numeric_condition_cols, selected = if (length(numeric_condition_cols) > 1) numeric_condition_cols[2] else if (length(numeric_condition_cols) > 0) numeric_condition_cols[1] else NULL),
    textInput("sgrna_paired_condition1_label", "Condition 1 Label (Optional):", placeholder = "e.g., Treatment"),
    textInput("sgrna_paired_condition2_label", "Condition 2 Label (Optional):", placeholder = "e.g., Control"),
    hr(),
    h4("Colors and Styles"),
    fluidRow(
      column(6, textInput("sgrna_paired_color_cond1", "Condition 1 Color:", value = "#007AFF")),
      column(6, textInput("sgrna_paired_color_cond2", "Condition 2 Color:", value = "#FC6A03"))
    ),
    textInput("sgrna_paired_line_color", "Connecting Line Color:", value = "#808080"),
    sliderInput("sgrna_paired_point_size", "Point Size:", min = 0.5, max = 5, value = 2.5, step = 0.1),
    sliderInput("sgrna_paired_line_size", "Connecting Line Width:", min = 0.1, max = 3, value = 0.5, step = 0.1),
    sliderInput("sgrna_paired_base_font_size", "Base Font Size:", min = 8, max = 20, value = 12, step = 1),
    hr(),
    h4("Titles and Labels"),
    textInput("sgrna_paired_plot_title", "Main Plot Title (Optional):", placeholder = "Leave blank for default title"),
    textInput("sgrna_paired_y_axis_label", "Y-axis Label (Optional):", value = "sgRNA Normalized Read Count"),
    actionButton("run_sgrna_paired_plot",
      "Generate sgRNA Paired Plot",
      class = "btn-primary",
      icon = icon("chart-line"),
      style = "width: 100%; white-space: normal; height: auto; padding: 12px 16px; line-height: 1.5; font-size: 16px; box-sizing: border-box;"
    ),
    tags$script(HTML("
          $(document).ready(function(){
              var sgrna_color_inputs = ['sgrna_paired_color_cond1', 'sgrna_paired_color_cond2', 'sgrna_paired_line_color'];
              sgrna_color_inputs.forEach(function(id) {
                  if ($('#' + id).length) { updateColorInputStyle(id); }
                  $('#' + id).off('input.' + id).on('input.' + id, function() { updateColorInputStyle(this.id); });
              });
          });
      "))
  )
}

# --- GSEA Lollipop Plot UI ---
render_gsea_lollipop_params_ui <- function(gsea_res_list) {
  req(gsea_res_list)

  gsea_choices <- character(0)
  if (!is.null(gsea_res_list$gsea_results_positive_df) && nrow(gsea_res_list$gsea_results_positive_df) > 0) {
    gsea_choices <- c(gsea_choices, "Positive Score GSEA" = "positive")
  }
  if (!is.null(gsea_res_list$gsea_results_negative_df) && nrow(gsea_res_list$gsea_results_negative_df) > 0) {
    gsea_choices <- c(gsea_choices, "Negative Score GSEA" = "negative")
  }

  if (length(gsea_choices) > 0) {
    div(
      class = "section-box", id = "gsea_lollipop_params_box",
      h2("3. GSEA Lollipop Plot Parameters"),
      selectInput("gsea_lollipop_result_type", "Select GSEA Result Type:", choices = gsea_choices, selected = gsea_choices[1]),
      numericInput("gsea_lollipop_num_pathways", "Number of Pathways to Show:", value = 10, min = 1, max = 30, step = 1),
      textInput("gsea_lollipop_bar_low", "Bar Low Value Color:", value = "#E0F3F8"),
      textInput("gsea_lollipop_bar_high", "Bar High Value Color:", value = "#045275"),
      textInput("gsea_lollipop_circle_low", "Circle Low Value Color:", value = "#FFF2CC"),
      textInput("gsea_lollipop_circle_high", "Circle High Value Color:", value = "#D66000"),
      sliderInput("gsea_lollipop_base_font", "Base Font Size:", value = 11, min = 7, max = 18, step = 1),
      sliderInput("gsea_lollipop_bar_width", "Bar Width Ratio:", value = 0.1, min = 0.1, max = 1.0, step = 0.05),
      textInput("gsea_lollipop_title", "Main Plot Title:", value = "GSEA enrichment pathway ranked plot"),
      actionButton("run_gsea_lollipop_plot",
        "Generate GSEA Lollipop Plot",
        class = "btn-primary",
        icon = icon("chart-bar"),
        style = "width: 100%; white-space: normal; height: auto; padding: 12px 16px; line-height: 1.5; font-size: 16px; box-sizing: border-box;"
      ),
      tags$script(HTML("
            $(document).ready(function(){
              var gsea_color_inputs = ['gsea_lollipop_bar_low', 'gsea_lollipop_bar_high', 'gsea_lollipop_circle_low', 'gsea_lollipop_circle_high'];
              gsea_color_inputs.forEach(function(id) {
                if ($('#' + id).length) { updateColorInputStyle(id); }
                $('#' + id).off('input.gseaColor').on('input.gseaColor', function() { updateColorInputStyle(this.id); });
              });
            });
          "))
    )
  } else {
    div(
      class = "section-box", id = "gsea_lollipop_params_box",
      h2("3. GSEA Lollipop Plot Parameters"),
      p("No GSEA results available (Positive or Negative) for plotting.")
    )
  }
}

# --- GSEA Enrichment Plot UI ---
render_gsea_enrichment_params_ui <- function(gsea_res_list) {
  req(gsea_res_list)

  gsea_choices_enrich <- character(0)
  if (!is.null(gsea_res_list$gsea_results_positive_df) && nrow(gsea_res_list$gsea_results_positive_df) > 0) {
    gsea_choices_enrich <- c(gsea_choices_enrich, "Positive Score GSEA" = "positive")
  }
  if (!is.null(gsea_res_list$gsea_results_negative_df) && nrow(gsea_res_list$gsea_results_negative_df) > 0) {
    gsea_choices_enrich <- c(gsea_choices_enrich, "Negative Score GSEA" = "negative")
  }
  if (length(gsea_choices_enrich) == 0) {
    return(div(class = "section-box", h2("4. GSEA Pathway Enrichment Plot Parameters"), p("No GSEA results available for pathway plotting.")))
  }

  div(
    class = "section-box", id = "gsea_enrichment_plot_params_box",
    h2("4. GSEA Pathway Enrichment Plot Parameters"),

    # GSEA Result Type Selection
    selectInput("gsea_enrichment_plot_result_type", "Select GSEA Result Type:",
      choices = gsea_choices_enrich,
      selected = gsea_choices_enrich[1]
    ),

    # Multi-pathway Selection
    tags$style(HTML("
      #gsea_enrichment_plot_pathway_selector + .selectize-control .selectize-input {
        white-space: normal !important;
        overflow: visible !important;
        height: auto !important;
        flex-wrap: wrap !important;
      }
      #gsea_enrichment_plot_pathway_selector + .selectize-control .selectize-input .item {
        white-space: normal !important;
        word-break: break-word !important;
        overflow-wrap: anywhere !important;
        max-width: 100% !important;
        display: inline-block !important;
      }
      #gsea_enrichment_plot_pathway_selector + .selectize-control .selectize-dropdown {
        width: 100% !important;
        max-width: 100% !important;
      }
      #gsea_enrichment_plot_pathway_selector + .selectize-control .selectize-dropdown .option {
        white-space: normal !important;
        word-break: break-word !important;
        overflow-wrap: anywhere !important;
      }
    ")),
    selectizeInput("gsea_enrichment_plot_pathway_selector",
      "Select Pathways for Visualization (Multiple Allowed):",
      choices = NULL,
      multiple = TRUE,
      width = "100%", # Ensure input takes full container width
      options = list(
        placeholder = "Please select result type first...",
        maxItems = 10,
        plugins = list("remove_button", "drag_drop")
      )
    ),

    # Pathway Color Settings (Only shown when multiple pathways are selected)
    conditionalPanel(
      condition = "input.gsea_enrichment_plot_pathway_selector && input.gsea_enrichment_plot_pathway_selector.length > 1",
      div(
        id = "pathway_color_settings", style = "margin: 15px 0;",
        h4("Pathway Color Settings"),
        div(
          class = "alert alert-warning", style = "font-size: 13px; padding: 8px;",
          p("💡 You can manually specify colors for each selected pathway.",
            style = "margin: 0;"
          )
        ),
        uiOutput("pathway_color_inputs")
      )
    ),

    # Gene Highlight Settings
    div(
      style = "margin: 15px 0; padding: 15px; background-color: #f8f9fa; border-radius: 8px;",
      h4("Gene Highlight Settings (Optional)"),
      checkboxInput("gsea_enrichment_enable_gene_highlight",
        "Enable Gene Highlighting",
        value = FALSE
      ),
      conditionalPanel(
        condition = "input.gsea_enrichment_enable_gene_highlight",
        div(
          style = "margin-top: 10px;",
          p("Please select genes to highlight from the core genes of the selected pathway:",
            style = "font-size: 14px; color: #495057;"
          ),
          selectizeInput("gsea_enrichment_highlight_genes",
            "Select Highlight Genes:",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Please select pathway first...",
              maxItems = 20,
              plugins = list("remove_button", "drag_drop")
            )
          ),
        )
      )
    ),

    # Plot Parameter Settings
    div(
      style = "margin: 15px 0;",
      h4("Plot Parameters"),
      sliderInput("gsea_enrichment_term_width",
        "Pathway Name Width:",
        value = 20, min = 10, max = 50, step = 1
      ),
      sliderInput("gsea_enrichment_base_font",
        "Base Font Size:",
        value = 10, min = 8, max = 16, step = 1
      ),
      numericInput("gsea_enrichment_legend_x",
        "Legend X Position:",
        value = 0.85, min = 0, max = 1, step = 0.05
      ),
      numericInput("gsea_enrichment_legend_y",
        "Legend Y Position:",
        value = 0.8, min = 0, max = 1, step = 0.05
      ),
      selectInput("gsea_enrichment_subplot_type",
        "Subplot Type:",
        choices = list("Enrichment Score + Gene Rank" = 2, "Enrichment Score Only" = 1),
        selected = 2
      ),
      checkboxInput("gsea_enrichment_add_pval",
        "Show p-value",
        value = TRUE
      ),
      conditionalPanel(
        condition = "input.gsea_enrichment_add_pval",
        numericInput("gsea_enrichment_pval_x",
          "p-value X Position:",
          value = 0.02, min = 0, max = 1, step = 0.01
        ),
        numericInput("gsea_enrichment_pval_y",
          "p-value Y Position:",
          value = 0.04, min = 0, max = 1, step = 0.01
        )
      )
    ),

    # Generate Button
    actionButton("run_gsea_enrichment_plot",
      "Generate GSEA Pathway Enrichment Plot",
      class = "btn-primary btn-lg",
      icon = icon("chart-line"),
      style = "margin-top: 15px; width: 100%; white-space: normal; height: auto; padding: 12px 16px; line-height: 1.5; font-size: 16px; box-sizing: border-box;"
    ),

    # Remove color input style script because highlight colors are fixed
  )
}

# --- Pathway Color Input UI ---
render_pathway_color_inputs <- function(selected_pathways) {
  if (is.null(selected_pathways) || length(selected_pathways) == 0) {
    return(div(p("Please select a pathway first", style = "color: #666; font-style: italic;")))
  }

  # Generate Default Colors
  default_colors <- if (length(selected_pathways) <= 8) {
    RColorBrewer::brewer.pal(max(3, length(selected_pathways)), "Set2")[1:length(selected_pathways)]
  } else {
    rainbow(length(selected_pathways))
  }

  # Generate Color Inputs
  color_inputs <- lapply(1:length(selected_pathways), function(i) {
    pathway_id <- selected_pathways[i]
    pathway_name <- gsub("^[A-Z_]+_", "", pathway_id) # Simplify pathway name display
    input_id <- paste0("pathway_color_", i)

    div(
      style = "margin: 5px 0; border-bottom: 1px solid #eee; padding-bottom: 5px;",
      div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        # Pathway Name (Flexible width, wraps text)
        div(
          style = "flex: 1; padding-right: 10px; min-width: 0;",
          strong(pathway_name,
            style = "font-size: 13px; color: #495057; display: block; word-wrap: break-word;"
          )
        ),
        # Color Input (Fixed width 120px)
        div(
          style = "flex: 0 0 120px;",
          textInput(input_id,
            label = NULL,
            value = default_colors[i],
            placeholder = "#RRGGBB",
            width = "100%"
          )
        )
      )
    )
  })

  return(div(
    do.call(tagList, color_inputs),
    tags$script(HTML(paste0(
      "
      $(document).ready(function(){
        var pathway_color_inputs = [",
      paste0("'pathway_color_", 1:length(selected_pathways), "'", collapse = ", "),
      "];
        pathway_color_inputs.forEach(function(id) {
          updateColorInputStyle(id);
          $('#' + id).off('input.pathwayColor').on('input.pathwayColor', function() {
            updateColorInputStyle(this.id);
          });
        });
      });
      "
    )))
  ))
}
# --- Gene Visualization Functions ---

# --- Gene Vis Params UI ---
render_gene_vis_params_ui <- function(results) {
  req(results)
  req(results$gene_summary_data)

  gene_summary_df <- results$gene_summary_data
  all_genes <- unique(gene_summary_df[[results$params$gene_col]])

  div(
    class = "section-box", id = "gene_vis_params_box",
    h2("6. Gene Visualization Parameters"),
    selectInput("gene_vis_plot_type", "Select Plot Type:",
      choices = c("Volcano Plot" = "volcano", "Gene Score Ranking" = "ranking"),
      selected = "volcano"
    ),

    # --- Volcano Plot Parameters ---
    conditionalPanel(
      condition = "input.gene_vis_plot_type == 'volcano'",
      h4("Volcano Plot Settings"),
      fluidRow(
        column(6, numericInput("gene_vis_volcano_p_cutoff", "P-value Cutoff:", value = 0.05, min = 0, max = 1, step = 0.01)),
        column(6, numericInput("gene_vis_volcano_score_cutoff", "Score Cutoff:", value = 0.2, min = 0, step = 0.1))
      ),
      selectizeInput("gene_vis_volcano_labels", "Genes to Label:",
        choices = NULL, multiple = TRUE,
        options = list(placeholder = "Select genes to label...", maxItems = 50)
      ),
      h5("Colors"),
      div(
        class = "color-inputs-vertical",
        textInput("gene_vis_volcano_col_up", "Up Color:", value = "#E31A1C"),
        textInput("gene_vis_volcano_col_down", "Down Color:", value = "#1F78B4"),
        textInput("gene_vis_volcano_col_ns", "NS Color:", value = "#ADB6B6")
      )
    ),

    # --- Ranking Plot Parameters ---
    conditionalPanel(
      condition = "input.gene_vis_plot_type == 'ranking'",
      h4("Ranking Plot Settings"),
      selectInput("gene_vis_ranking_type", "Ranking Type:",
        choices = c("Positive Score" = "positive", "Negative Score" = "negative"),
        selected = "positive"
      ),
      numericInput("gene_vis_ranking_top_n", "Top N Genes:", value = 10, min = 5, max = 50, step = 1),
      h5("Gradient Colors"),
      fluidRow(
        column(6, textInput("gene_vis_ranking_col_low", "Low Score Color:", value = "#fee0d2")),
        column(6, textInput("gene_vis_ranking_col_high", "High Score Color:", value = "#a50f15"))
      )
    ),
    tags$hr(),
    h4("General Settings"),
    sliderInput("gene_vis_base_font_size", "Base Font Size:", min = 8, max = 24, value = 14, step = 1),
    textInput("gene_vis_plot_title", "Plot Title (Optional):", placeholder = "Leave blank for default"),
    actionButton("run_gene_vis_plot",
      "Generate Plot",
      class = "btn-primary",
      icon = icon("paint-brush"),
      style = "width: 100%; white-space: normal; height: auto; padding: 12px 16px; line-height: 1.5; font-size: 16px; box-sizing: border-box;"
    ),
    tags$script(HTML("
          $(document).ready(function(){
              var gene_vis_colors = ['gene_vis_volcano_col_up', 'gene_vis_volcano_col_down', 'gene_vis_volcano_col_ns',
                                     'gene_vis_ranking_col_low', 'gene_vis_ranking_col_high'];
              gene_vis_colors.forEach(function(id) {
                  if ($('#' + id).length) { updateColorInputStyle(id); }
                  $('#' + id).off('input.' + id).on('input.' + id, function() { updateColorInputStyle(this.id); });
              });
          });
      "))
  )
}

# --- Gene Vis Observer ---
observe_gene_vis_plot_generation <- function(input, output, session, data_processing_results,
                                             gene_vis_plot_object, gene_vis_plot_params_for_download) {
  # Auto-update colors based on ranking type
  observeEvent(input$gene_vis_ranking_type, {
    if (input$gene_vis_ranking_type == "positive") {
      updateTextInput(session, "gene_vis_ranking_col_low", value = "#fee0d2")
      updateTextInput(session, "gene_vis_ranking_col_high", value = "#a50f15")
    } else if (input$gene_vis_ranking_type == "negative") {
      updateTextInput(session, "gene_vis_ranking_col_low", value = "#08519c")
      updateTextInput(session, "gene_vis_ranking_col_high", value = "#deebf7")
    }
  })

  observeEvent(input$run_gene_vis_plot, {
    req(data_processing_results()$gene_summary_data)

    plot_type <- input$gene_vis_plot_type
    gene_summary <- data_processing_results()$gene_summary_data

    output$gene_vis_plot_status <- renderPrint({
      "Generating plot..."
    })
    gene_vis_plot_object(NULL)
    shinyjs::hide("gene_vis_download_options")

    tryCatch(
      {
        p <- NULL
        params <- list()

        if (plot_type == "volcano") {
          p <- generate_volcano_plot(
            data = gene_summary,
            p_value_cutoff = input$gene_vis_volcano_p_cutoff,
            score_cutoff = input$gene_vis_volcano_score_cutoff,
            label_genes = input$gene_vis_volcano_labels,
            col_up = input$gene_vis_volcano_col_up,
            col_down = input$gene_vis_volcano_col_down,
            col_ns = input$gene_vis_volcano_col_ns,
            base_font_size = input$gene_vis_base_font_size,
            plot_title = if (input$gene_vis_plot_title == "") "Volcano Plot" else input$gene_vis_plot_title
          )
          params <- list(type = "volcano")
        } else if (plot_type == "ranking") {
          p <- generate_ranking_plot(
            data = gene_summary,
            type = input$gene_vis_ranking_type,
            top_n = input$gene_vis_ranking_top_n,
            gradient_low = input$gene_vis_ranking_col_low,
            gradient_high = input$gene_vis_ranking_col_high,
            base_font_size = input$gene_vis_base_font_size,
            plot_title = if (input$gene_vis_plot_title == "") NULL else input$gene_vis_plot_title
          )
          params <- list(type = "ranking")
        }

        gene_vis_plot_object(p)
        gene_vis_plot_params_for_download(params)

        output$gene_vis_plot_status <- renderPrint({
          paste("Plot generated successfully at", Sys.time())
        })
        shinyjs::show("gene_vis_download_options")
        updateTabsetPanel(session, "main_results_tabs", selected = "Gene Visualization")
      },
      error = function(e) {
        output$gene_vis_plot_status <- renderPrint({
          paste("Error:", e$message)
        })
      }
    )
  })
}

# --- Gene Vis Download Handler ---
create_gene_vis_download_handler <- function(input, gene_vis_plot_object, gene_vis_plot_params_for_download) {
  downloadHandler(
    filename = function() {
      type <- gene_vis_plot_params_for_download()$type
      ext <- input$gene_vis_download_format
      paste0("gene_vis_", type, "_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".", ext)
    },
    content = function(file) {
      req(gene_vis_plot_object())
      p <- gene_vis_plot_object()

      width <- input$gene_vis_download_width
      height <- input$gene_vis_download_height
      dpi <- if (!is.null(input$gene_vis_download_dpi)) input$gene_vis_download_dpi else 300

      ggsave(file, plot = p, width = width, height = height, dpi = dpi, device = input$gene_vis_download_format)
    }
  )
}


# --- sgRNA Paired Plot UI ---
render_sgrna_params_ui_with_gsea <- function(results, gsea_res_list) {
  req(results, gsea_res_list)
  req(results$normalized_counts, results$gene_summary_data, results$params$gene_col, results$params$gRNA_col)

  norm_counts_df <- results$normalized_counts
  gene_summary_df <- results$gene_summary_data
  app_params <- results$params

  # Get numeric data columns (excluding IDs)
  id_cols_in_norm <- c(app_params$gRNA_col, app_params$gene_col, app_params$sequence_col)
  id_cols_in_norm <- id_cols_in_norm[!sapply(id_cols_in_norm, is.null)]

  potential_condition_cols <- setdiff(names(norm_counts_df), id_cols_in_norm)
  numeric_condition_cols <- character(0)
  if (length(potential_condition_cols) > 0) {
    is_numeric_check <- sapply(norm_counts_df[, potential_condition_cols, drop = FALSE], is.numeric)
    numeric_condition_cols <- potential_condition_cols[is_numeric_check]
  }

  if (length(numeric_condition_cols) < 2) {
    return(div(class = "section-box", id = "sgrna_paired_plot_params_box", h2("5. sgRNA Paired Plot Parameters"), p("Not enough numeric columns in normalized count data for paired plot (requires at least 2).")))
  }

  # Get all genes for selection
  all_genes <- unique(gene_summary_df[[app_params$gene_col]])

  # Prepare GSEA selector choices
  gsea_type_choices <- character(0)
  if (!is.null(gsea_res_list$gsea_results_positive_df) && nrow(gsea_res_list$gsea_results_positive_df) > 0) {
    gsea_type_choices <- c(gsea_type_choices, "Positive Score GSEA" = "positive")
  }
  if (!is.null(gsea_res_list$gsea_results_negative_df) && nrow(gsea_res_list$gsea_results_negative_df) > 0) {
    gsea_type_choices <- c(gsea_type_choices, "Negative Score GSEA" = "negative")
  }

  div(
    class = "section-box", id = "sgrna_paired_plot_params_box",
    h2("5. sgRNA Paired Plot Parameters"),

    # Gene selection method
    radioButtons("sgrna_gene_selection_method", "Select Gene Selection Method:",
      choices = list("Direct Gene Selection" = "direct", "Select via GSEA Results" = "gsea"),
      selected = "direct", inline = TRUE
    ),

    # Direct gene selection UI
    conditionalPanel(
      condition = "input.sgrna_gene_selection_method == 'direct'",
      selectizeInput("sgrna_paired_gene_selector_direct", "Select Gene:",
        choices = all_genes,
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = "Please select gene...")
      ),
      conditionalPanel(
        condition = "input.sgrna_paired_gene_selector_direct && input.sgrna_paired_gene_selector_direct.length > 0",
        p("Selected", tags$span(id = "direct_gene_count_gsea", "0"), "genes", style = "color: #0066cc; font-weight: bold;"),
        tags$script(HTML("
            $(document).on('change', '#sgrna_paired_gene_selector_direct', function() {
              var count = $('#sgrna_paired_gene_selector_direct').val() ? $('#sgrna_paired_gene_selector_direct').val().length : 0;
              $('#direct_gene_count_gsea').text(count);
            });
          "))
      )
    ),

    # Gene selection via GSEA UI
    conditionalPanel(
      condition = "input.sgrna_gene_selection_method == 'gsea'",
      selectizeInput("sgrna_paired_gsea_type_selector", "1. Select GSEA Result Source:",
        choices = gsea_type_choices,
        selected = if (length(gsea_type_choices) > 0) gsea_type_choices[1] else NULL
      ),
      tags$style(HTML("
        #sgrna_paired_pathway_selector + .selectize-control .selectize-input {
          white-space: normal !important;
          overflow: visible !important;
          height: auto !important;
          flex-wrap: wrap !important;
        }
        #sgrna_paired_pathway_selector + .selectize-control .selectize-input .item {
          white-space: normal !important;
          word-break: break-word !important;
          overflow-wrap: anywhere !important;
          max-width: 100% !important;
          display: inline-block !important;
        }
        #sgrna_paired_pathway_selector + .selectize-control .selectize-dropdown {
          width: 100% !important;
          max-width: 100% !important;
        }
        #sgrna_paired_pathway_selector + .selectize-control .selectize-dropdown .option {
          white-space: normal !important;
          word-break: break-word !important;
          overflow-wrap: anywhere !important;
        }
      ")),
      selectizeInput("sgrna_paired_pathway_selector", "2. Select Pathway:",
        choices = NULL,
        width = "100%",
        options = list(placeholder = "Please select GSEA result source first...")
      ),
      checkboxInput("sgrna_paired_plot_all_genes", "Generate plots for all core genes in pathway (Batch Download Mode)", value = FALSE),
      conditionalPanel(
        condition = "!input.sgrna_paired_plot_all_genes",
        selectizeInput("sgrna_paired_gene_selector", "3. Select Single Gene (Preview Mode):", choices = NULL, options = list(placeholder = "Please select pathway first..."))
      )
    ),
    tags$hr(),
    h4("4. Select Plotting Data Columns"),
    selectInput("sgrna_paired_condition1_col", "Select Condition 1 Column:", choices = numeric_condition_cols, selected = if (length(numeric_condition_cols) > 0) numeric_condition_cols[1] else NULL),
    selectInput("sgrna_paired_condition2_col", "Select Condition 2 Column:", choices = numeric_condition_cols, selected = if (length(numeric_condition_cols) > 1) numeric_condition_cols[2] else if (length(numeric_condition_cols) > 0) numeric_condition_cols[1] else NULL),
    textInput("sgrna_paired_condition1_label", "Condition 1 Label (Optional):", placeholder = "e.g., Treatment"),
    textInput("sgrna_paired_condition2_label", "Condition 2 Label (Optional):", placeholder = "e.g., Control"),
    hr(),
    h4("Colors and Styles"),
    fluidRow(
      column(6, textInput("sgrna_paired_color_cond1", "Condition 1 Color:", value = "#007AFF")),
      column(6, textInput("sgrna_paired_color_cond2", "Condition 2 Color:", value = "#FC6A03"))
    ),
    textInput("sgrna_paired_line_color", "Connecting Line Color:", value = "#808080"),
    sliderInput("sgrna_paired_point_size", "Point Size:", min = 0.5, max = 5, value = 2.5, step = 0.1),
    sliderInput("sgrna_paired_line_size", "Connecting Line Width:", min = 0.1, max = 3, value = 0.5, step = 0.1),
    sliderInput("sgrna_paired_base_font_size", "Base Font Size:", min = 8, max = 20, value = 12, step = 1),
    hr(),
    h4("Titles and Labels"),
    textInput("sgrna_paired_plot_title", "Main Plot Title (Optional):", placeholder = "Leave blank for default title"),
    textInput("sgrna_paired_y_axis_label", "Y-axis Label (Optional):", value = "sgRNA Normalized Read Count"),
    actionButton("run_sgrna_paired_plot",
      "Generate sgRNA Paired Plot",
      class = "btn-primary",
      icon = icon("chart-line"),
      style = "width: 100%; white-space: normal; height: auto; padding: 12px 16px; line-height: 1.5; font-size: 16px; box-sizing: border-box;"
    ),
    tags$script(HTML("
          $(document).ready(function(){
              var sgrna_color_inputs = ['sgrna_paired_color_cond1', 'sgrna_paired_color_cond2', 'sgrna_paired_line_color'];
              sgrna_color_inputs.forEach(function(id) {
                  if ($('#' + id).length) { updateColorInputStyle(id); }
                  $('#' + id).off('input.' + id).on('input.' + id, function() { updateColorInputStyle(this.id); });
              });
          });
      "))
  )
}

# --- sgRNA Paired Plot Functions ---
# sgRNA Paired Plot Observer
observe_sgrna_paired_plot_generation <- function(input, output, session, data_processing_results, gsea_results,
                                                 sgrna_paired_plot_object, sgrna_paired_plot_params_for_download,
                                                 sgrna_plot_batch_genes, sgrna_plot_batch_active) {
  observeEvent(input$run_sgrna_paired_plot, {
    req(
      data_processing_results()$normalized_counts,
      data_processing_results()$params$gRNA_col,
      data_processing_results()$params$gene_col,
      input$sgrna_paired_condition1_col,
      input$sgrna_paired_condition2_col
    )

    # Common variables
    norm_counts <- data_processing_results()$normalized_counts
    app_params <- data_processing_results()$params
    cond1_col <- input$sgrna_paired_condition1_col
    cond2_col <- input$sgrna_paired_condition2_col

    # Validate common inputs
    if (cond1_col == cond2_col) {
      showNotification("Condition 1 and Condition 2 columns cannot be the same.", type = "error")
      return()
    }
    if (!is_valid_hex_or_name(input$sgrna_paired_color_cond1)) {
      showNotification("Condition 1 color is not a valid hex code or color name.", type = "error")
      return()
    }
    if (!is_valid_hex_or_name(input$sgrna_paired_color_cond2)) {
      showNotification("Condition 2 color is not a valid hex code or color name.", type = "error")
      return()
    }
    if (!is_valid_hex_or_name(input$sgrna_paired_line_color)) {
      showNotification("Connecting line color is not a valid hex code or color name.", type = "error")
      return()
    }

    # Check selection method
    selection_method <- input$sgrna_gene_selection_method
    if (is.null(selection_method)) selection_method <- "direct" # Default to direct

    if (selection_method == "direct") {
      # Direct selection mode
      # Check number of selected genes
      selected_genes <- input$sgrna_paired_gene_selector_direct
      if (!is.null(selected_genes) && length(selected_genes) > 1) {
        # Batch mode (multiple genes)
        sgrna_plot_batch_genes(selected_genes)
        sgrna_plot_batch_active(TRUE)
        sgrna_paired_plot_object(NULL)
        output$sgrna_paired_plot_batch_status <- renderPrint({
          paste0("Ready to batch download paired plots for ", length(selected_genes), " genes. Click the button below to download the ZIP file.")
        })
        shinyjs::hide("sgrna_paired_plot_download_options")
        updateTabsetPanel(session, "main_results_tabs", selected = "sgRNA Paired Plot")
      } else if (!is.null(selected_genes) && length(selected_genes) == 1) {
        # Single gene mode
        target_gene <- selected_genes[1]

        # Single plot logic
        if (!is.null(sgrna_plot_batch_active)) {
          sgrna_plot_batch_active(FALSE)
        }

        output$sgrna_paired_plot_status <- renderPrint({
          "Generating sgRNA Paired Plot..."
        })
        sgrna_paired_plot_object(NULL)
        if (!is.null(sgrna_paired_plot_params_for_download)) {
          sgrna_paired_plot_params_for_download(NULL)
        }
        shinyjs::hide("sgrna_paired_plot_download_options")

        tryCatch(
          {
            plot_obj <- generate_sgrna_paired_plot(
              normalized_counts_df = norm_counts,
              target_gene_id = target_gene,
              sgrna_id_col = app_params$gRNA_col,
              gene_col = app_params$gene_col,
              condition_1_col = cond1_col,
              condition_2_col = cond2_col,
              condition_1_label = if (input$sgrna_paired_condition1_label == "") NULL else input$sgrna_paired_condition1_label,
              condition_2_label = if (input$sgrna_paired_condition2_label == "") NULL else input$sgrna_paired_condition2_label,
              color_condition_1 = input$sgrna_paired_color_cond1,
              color_condition_2 = input$sgrna_paired_color_cond2,
              connecting_line_color = input$sgrna_paired_line_color,
              connecting_line_size = input$sgrna_paired_line_size,
              point_size = input$sgrna_paired_point_size,
              base_font_size = input$sgrna_paired_base_font_size,
              plot_title = if (input$sgrna_paired_plot_title == "") NULL else input$sgrna_paired_plot_title,
              y_axis_label = input$sgrna_paired_y_axis_label
            )
            sgrna_paired_plot_object(plot_obj)
            if (!is.null(sgrna_paired_plot_params_for_download)) {
              sgrna_paired_plot_params_for_download(list(
                gene = target_gene,
                cond1 = cond1_col,
                cond2 = cond2_col,
                title = input$sgrna_paired_plot_title
              ))
            }
            output$sgrna_paired_plot_status <- renderPrint({
              paste0("sgRNA Paired Plot (Gene: ", target_gene, ") generated. ", Sys.time())
            })
            shinyjs::show("sgrna_paired_plot_download_options")
            updateTabsetPanel(session, "main_results_tabs", selected = "sgRNA Paired Plot")
          },
          error = function(e) {
            error_msg <- paste("Error generating sgRNA Paired Plot:\n", e$message)
            output$sgrna_paired_plot_status <- renderPrint({
              error_msg
            })
            sgrna_paired_plot_object(NULL)
            showModal(modalDialog(title = "sgRNA Paired Plot Generation Failed", HTML(paste0("Error Details: <br><pre style='white-space: pre-wrap;'>", e$message, "</pre>"))))
          }
        )
      } else {
        # No gene selected warning
        showNotification("Please select at least one gene to plot.", type = "warning")
        return()
      }
    } else if (selection_method == "gsea") {
      # GSEA selection mode
      # Check for GSEA batch mode
      if (!is.null(input$sgrna_paired_plot_all_genes) && isTRUE(input$sgrna_paired_plot_all_genes)) {
        # Batch mode (requires GSEA)
        req(input$sgrna_paired_pathway_selector, gsea_results())
        selected_pathway_id <- input$sgrna_paired_pathway_selector
        if (is.null(selected_pathway_id) || selected_pathway_id == "") {
          showNotification("Please select a pathway for batch plotting.", type = "warning")
          return()
        }

        gsea_res_list <- gsea_results()
        selected_gsea_type <- input$sgrna_paired_gsea_type_selector
        source_gsea_df <- NULL
        if (selected_gsea_type == "positive" && !is.null(gsea_res_list$gsea_results_positive_df)) {
          source_gsea_df <- gsea_res_list$gsea_results_positive_df
        } else if (selected_gsea_type == "negative" && !is.null(gsea_res_list$gsea_results_negative_df)) {
          source_gsea_df <- gsea_res_list$gsea_results_negative_df
        }

        if (is.null(source_gsea_df) || !(selected_pathway_id %in% source_gsea_df$ID)) {
          showNotification("Cannot find GSEA data for the selected pathway.", type = "error")
          return()
        }

        core_enrich_str <- source_gsea_df[source_gsea_df$ID == selected_pathway_id, "core_enrichment", drop = TRUE]
        genes_to_plot <- character(0)
        if (length(core_enrich_str) == 1 && !is.na(core_enrich_str) && nzchar(core_enrich_str)) {
          core_genes <- unique(strsplit(core_enrich_str, "/")[[1]])
          core_genes <- core_genes[nzchar(core_genes)]
          if (length(core_genes) > 0) {
            norm_counts_genes <- data_processing_results()$normalized_counts[[data_processing_results()$params$gene_col]]
            genes_to_plot <- intersect(core_genes, unique(norm_counts_genes))
          }
        }

        if (length(genes_to_plot) == 0) {
          output$sgrna_paired_plot_batch_status <- renderPrint({
            "Error: No core genes in the selected pathway, or core genes have no plottable sgRNA in the data."
          })
          sgrna_plot_batch_genes(NULL)
          sgrna_plot_batch_active(FALSE)
          return()
        }

        sgrna_plot_batch_genes(genes_to_plot)
        sgrna_plot_batch_active(TRUE)
        sgrna_paired_plot_object(NULL)
        output$sgrna_paired_plot_batch_status <- renderPrint({
          paste0("Ready to batch download paired plots for ", length(genes_to_plot), " core genes. Click the button below to download the ZIP file.")
        })
        shinyjs::hide("sgrna_paired_plot_download_options")
        updateTabsetPanel(session, "main_results_tabs", selected = "sgRNA Paired Plot")
      } else {
        # Single plot mode (GSEA)
        if (!is.null(sgrna_plot_batch_active)) {
          sgrna_plot_batch_active(FALSE)
        }
        req(input$sgrna_paired_gene_selector)
        target_gene <- input$sgrna_paired_gene_selector
        if (is.null(target_gene) || target_gene == "") {
          showNotification("Please select a gene to plot.", type = "warning")
          return()
        }

        output$sgrna_paired_plot_status <- renderPrint({
          "Generating sgRNA Paired Plot..."
        })
        sgrna_paired_plot_object(NULL)
        if (!is.null(sgrna_paired_plot_params_for_download)) {
          sgrna_paired_plot_params_for_download(NULL)
        }
        shinyjs::hide("sgrna_paired_plot_download_options")

        tryCatch(
          {
            plot_obj <- generate_sgrna_paired_plot(
              normalized_counts_df = norm_counts,
              target_gene_id = target_gene,
              sgrna_id_col = app_params$gRNA_col,
              gene_col = app_params$gene_col,
              condition_1_col = cond1_col,
              condition_2_col = cond2_col,
              condition_1_label = if (input$sgrna_paired_condition1_label == "") NULL else input$sgrna_paired_condition1_label,
              condition_2_label = if (input$sgrna_paired_condition2_label == "") NULL else input$sgrna_paired_condition2_label,
              color_condition_1 = input$sgrna_paired_color_cond1,
              color_condition_2 = input$sgrna_paired_color_cond2,
              connecting_line_color = input$sgrna_paired_line_color,
              connecting_line_size = input$sgrna_paired_line_size,
              point_size = input$sgrna_paired_point_size,
              base_font_size = input$sgrna_paired_base_font_size,
              plot_title = if (input$sgrna_paired_plot_title == "") NULL else input$sgrna_paired_plot_title,
              y_axis_label = input$sgrna_paired_y_axis_label
            )
            sgrna_paired_plot_object(plot_obj)
            if (!is.null(sgrna_paired_plot_params_for_download)) {
              sgrna_paired_plot_params_for_download(list(
                gene = target_gene,
                cond1 = cond1_col,
                cond2 = cond2_col,
                title = input$sgrna_paired_plot_title
              ))
            }
            output$sgrna_paired_plot_status <- renderPrint({
              paste0("sgRNA Paired Plot (Gene: ", target_gene, ") generated. ", Sys.time())
            })
            shinyjs::show("sgrna_paired_plot_download_options")
            updateTabsetPanel(session, "main_results_tabs", selected = "sgRNA Paired Plot")
          },
          error = function(e) {
            error_msg <- paste("Error generating sgRNA Paired Plot:\n", e$message)
            output$sgrna_paired_plot_status <- renderPrint({
              error_msg
            })
            sgrna_paired_plot_object(NULL)
            showModal(modalDialog(title = "sgRNA Paired Plot Generation Failed", HTML(paste0("Error Details: <br><pre style='white-space: pre-wrap;'>", e$message, "</pre>"))))
          }
        )
      }
    }
  })
}

# sgRNA Paired Plot Download Handler
create_sgrna_paired_plot_download_handler <- function(input, sgrna_paired_plot_object, sgrna_paired_plot_params_for_download) {
  downloadHandler(
    filename = function() {
      req(sgrna_paired_plot_params_for_download())
      dl_params <- sgrna_paired_plot_params_for_download()
      gene_name_safe <- gsub("[^a-zA-Z0-9_.-]", "_", dl_params$gene)
      title_part <- if (!is.null(dl_params$title) && nzchar(dl_params$title)) gsub("[^a-zA-Z0-9_.-]", "_", dl_params$title) else paste0(dl_params$cond1, "_vs_", dl_params$cond2)
      title_part <- substr(title_part, 1, 50)

      paste0("sgRNA_PairedPlot_", gene_name_safe, "_", title_part, "_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".", input$sgrna_paired_download_format)
    },
    content = function(file) {
      req(sgrna_paired_plot_object(), input$sgrna_paired_download_width, input$sgrna_paired_download_height)

      ggsave(file,
        plot = sgrna_paired_plot_object(), device = input$sgrna_paired_download_format,
        width = input$sgrna_paired_download_width, height = input$sgrna_paired_download_height,
        units = "in", dpi = if (input$sgrna_paired_download_format == "png") input$sgrna_paired_download_dpi else 300,
        bg = "white"
      )
    }
  )
}

# sgRNA Batch Download Handler
create_sgrna_batch_download_handler <- function(input, output, sgrna_plot_batch_genes, data_processing_results, gsea_results) {
  downloadHandler(
    filename = function() {
      req(input$sgrna_paired_pathway_selector)
      pathway_id <- input$sgrna_paired_pathway_selector
      # Try to get pathway description for filename
      pathway_desc <- pathway_id # Fallback to ID
      gsea_res_list <- gsea_results()
      selected_gsea_type <- input$sgrna_paired_gsea_type_selector
      source_df_for_desc <- NULL
      if (!is.null(gsea_res_list)) {
        if (selected_gsea_type == "positive" && !is.null(gsea_res_list$gsea_results_positive_df)) {
          source_df_for_desc <- gsea_res_list$gsea_results_positive_df
        } else if (selected_gsea_type == "negative" && !is.null(gsea_res_list$gsea_results_negative_df)) {
          source_df_for_desc <- gsea_res_list$gsea_results_negative_df
        }
        if (!is.null(source_df_for_desc) && pathway_id %in% source_df_for_desc$ID) {
          pathway_desc <- source_df_for_desc[source_df_for_desc$ID == pathway_id, "Description"][1]
        }
      }
      safe_pathway_desc <- gsub("[^a-zA-Z0-9_.-]", "_", pathway_desc)
      safe_pathway_desc <- substr(safe_pathway_desc, 1, 50)
      paste0("sgRNA_PairedPlots_Pathway_", safe_pathway_desc, "_Batch_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(
        sgrna_plot_batch_genes(), data_processing_results(),
        input$sgrna_paired_condition1_col, input$sgrna_paired_condition2_col,
        input$sgrna_batch_download_format,
        input$sgrna_batch_download_width,
        input$sgrna_batch_download_height
      )

      genes_to_plot_list <- sgrna_plot_batch_genes()
      if (is.null(genes_to_plot_list) || length(genes_to_plot_list) == 0) {
        stop("No genes prepared for batch plotting.")
      }

      output$sgrna_paired_plot_batch_status <- renderPrint({
        "Generating and compressing plots, please wait..."
      })

      temp_dir <- tempfile("paired_plots_")
      dir.create(temp_dir)
      on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

      plot_files <- c()

      # Common plotting params
      norm_counts_data <- data_processing_results()$normalized_counts
      app_params_data <- data_processing_results()$params
      cond1_col_plot <- input$sgrna_paired_condition1_col
      cond2_col_plot <- input$sgrna_paired_condition2_col
      cond1_label_plot <- if (input$sgrna_paired_condition1_label == "") NULL else input$sgrna_paired_condition1_label
      cond2_label_plot <- if (input$sgrna_paired_condition2_label == "") NULL else input$sgrna_paired_condition2_label
      color_cond1_plot <- input$sgrna_paired_color_cond1
      color_cond2_plot <- input$sgrna_paired_color_cond2
      line_color_plot <- input$sgrna_paired_line_color
      line_size_plot <- input$sgrna_paired_line_size
      point_size_plot <- input$sgrna_paired_point_size
      base_font_plot <- input$sgrna_paired_base_font_size
      y_label_plot <- input$sgrna_paired_y_axis_label

      # Batch download settings
      plot_format_batch <- input$sgrna_batch_download_format
      plot_width_batch <- input$sgrna_batch_download_width
      plot_height_batch <- input$sgrna_batch_download_height
      plot_dpi_batch <- if (plot_format_batch == "png") input$sgrna_batch_download_dpi else 300
      if (is.null(plot_dpi_batch) && plot_format_batch == "png") plot_dpi_batch <- 300

      num_genes <- length(genes_to_plot_list)

      withProgress(message = "Generating batch paired plots...", value = 0, {
        for (i in seq_along(genes_to_plot_list)) {
          target_gene_plot <- genes_to_plot_list[i]
          incProgress(1 / num_genes, detail = paste("Processing Gene:", target_gene_plot, "(", i, "/", num_genes, ")"))

          plot_title_gene <- paste("sgRNA Paired Plot for Gene:", target_gene_plot)

          plot_obj_gene <- tryCatch(
            {
              generate_sgrna_paired_plot(
                normalized_counts_df = norm_counts_data,
                target_gene_id = target_gene_plot,
                sgrna_id_col = app_params_data$gRNA_col,
                gene_col = app_params_data$gene_col,
                condition_1_col = cond1_col_plot,
                condition_2_col = cond2_col_plot,
                condition_1_label = cond1_label_plot,
                condition_2_label = cond2_label_plot,
                color_condition_1 = color_cond1_plot,
                color_condition_2 = color_cond2_plot,
                connecting_line_color = line_color_plot,
                connecting_line_size = line_size_plot,
                point_size = point_size_plot,
                base_font_size = base_font_plot,
                plot_title = plot_title_gene,
                y_axis_label = y_label_plot
              )
            },
            error = function(e) {
              cat(paste("Failed to generate paired plot for gene:", target_gene_plot, ":", e$message, "\n"))
              NULL
            }
          )

          if (!is.null(plot_obj_gene)) {
            safe_gene_name <- gsub("[^a-zA-Z0-9_.-]", "_", target_gene_plot)
            plot_filename <- file.path(temp_dir, paste0("PairedPlot_", safe_gene_name, ".", plot_format_batch))
            tryCatch(
              {
                ggsave(plot_filename,
                  plot = plot_obj_gene, device = plot_format_batch,
                  width = plot_width_batch, height = plot_height_batch, units = "in", dpi = plot_dpi_batch, bg = "white"
                )
                plot_files <- c(plot_files, plot_filename)
              },
              error = function(e_save) {
                cat(paste("Failed to save paired plot for gene:", target_gene_plot, ":", e_save$message, "\n"))
              }
            )
          }
        }
      })

      if (length(plot_files) > 0) {
        output$sgrna_paired_plot_batch_status <- renderPrint({
          paste0("Generated ", length(plot_files), " / ", num_genes, " plots. Compressing...")
        })
        zip::zip(zipfile = file, files = plot_files, mode = "cherry-pick")
        output$sgrna_paired_plot_batch_status <- renderPrint({
          paste0("ZIP created with ", length(plot_files), " plots. Download complete.")
        })
      } else {
        output$sgrna_paired_plot_batch_status <- renderPrint({
          "Error: Failed to generate any plots for compression."
        })
        file.create(file)
      }
    },
    contentType = "application/zip"
  )
}

# --- GSEA Plot Observers & Downloaders ---

# GSEA Lollipop Plot Observer
observe_gsea_lollipop_plot_generation <- function(input, output, session, gsea_results, gsea_lollipop_plot_object, gsea_lollipop_plot_params_for_download) {
  observeEvent(input$run_gsea_lollipop_plot, {
    req(gsea_results(), input$gsea_lollipop_result_type, input$gsea_lollipop_num_pathways)

    output$gsea_lollipop_plot_status <- renderPrint({
      "Generating GSEA Lollipop Plot..."
    })
    gsea_lollipop_plot_object(NULL)
    gsea_lollipop_plot_params_for_download(NULL)
    shinyjs::hide("gsea_lollipop_plot_download_options")

    selected_gsea_type <- input$gsea_lollipop_result_type
    gsea_res_list <- gsea_results()
    gsea_df_to_plot <- NULL
    analysis_label_for_plot <- ""

    if (selected_gsea_type == "positive" && !is.null(gsea_res_list$gsea_results_positive_df)) {
      gsea_df_to_plot <- gsea_res_list$gsea_results_positive_df
      analysis_label_for_plot <- "Positive Score"
    } else if (selected_gsea_type == "negative" && !is.null(gsea_res_list$gsea_results_negative_df)) {
      gsea_df_to_plot <- gsea_res_list$gsea_results_negative_df
      analysis_label_for_plot <- "Negative Score"
    } else {
      output$gsea_lollipop_plot_status <- renderPrint({
        "Error: Selected GSEA result data not found."
      })
      return()
    }

    # Validate color inputs
    hex_inputs <- c(input$gsea_lollipop_bar_low, input$gsea_lollipop_bar_high, input$gsea_lollipop_circle_low, input$gsea_lollipop_circle_high)
    hex_labels <- c("Bar Low Color", "Bar High Color", "Circle Low Color", "Circle High Color")

    if (!validate_hex_inputs(hex_inputs, hex_labels, session)) {
      return()
    }

    tryCatch(
      {
        plot_obj <- generate_gsea_lollipop_plot_shiny(
          gsea_data_input_df = gsea_df_to_plot,
          analysis_type_label = analysis_label_for_plot,
          num_pathways_to_plot = input$gsea_lollipop_num_pathways,
          plot_title_main = input$gsea_lollipop_title,
          bar_fill_low_hex = input$gsea_lollipop_bar_low,
          bar_fill_high_hex = input$gsea_lollipop_bar_high,
          circle_fill_low_hex = input$gsea_lollipop_circle_low,
          circle_fill_high_hex = input$gsea_lollipop_circle_high,
          base_font_size = input$gsea_lollipop_base_font,
          bar_width_ratio = input$gsea_lollipop_bar_width
        )
        gsea_lollipop_plot_object(plot_obj)
        gsea_lollipop_plot_params_for_download(list(title = input$gsea_lollipop_title, type = analysis_label_for_plot))
        output$gsea_lollipop_plot_status <- renderPrint({
          paste0("GSEA Lollipop Plot (", analysis_label_for_plot, ") generated. ", Sys.time())
        })
        shinyjs::show("gsea_lollipop_plot_download_options")
        updateTabsetPanel(session, "main_results_tabs", selected = "GSEA Lollipop Plot")
      },
      error = function(e) {
        error_msg_gsea_lollipop <- paste("Error generating GSEA Lollipop Plot:\n", e$message)
        output$gsea_lollipop_plot_status <- renderPrint({
          error_msg_gsea_lollipop
        })
        gsea_lollipop_plot_object(NULL)
        showModal(modalDialog(title = "GSEA Lollipop Plot Generation Failed", HTML(paste0("Error Details: <br><pre style='white-space: pre-wrap;'>", e$message, "</pre>"))))
      }
    )
  })
}

# GSEA Lollipop Download Handler
create_gsea_lollipop_download_handler <- function(input, gsea_lollipop_plot_object, gsea_lollipop_plot_params_for_download) {
  downloadHandler(
    filename = function() {
      req(gsea_lollipop_plot_params_for_download())
      dl_params <- gsea_lollipop_plot_params_for_download()
      safe_title <- gsub("[^a-zA-Z0-9_]", "_", dl_params$title)
      safe_type <- gsub("[^a-zA-Z0-9_]", "_", dl_params$type)
      paste0("GSEA_LollipopPlot_", safe_type, "_", safe_title, "_", format(Sys.time(), "%Y_%m_%d_%H%M%S"), ".", input$gsea_lollipop_download_format)
    },
    content = function(file) {
      req(gsea_lollipop_plot_object(), input$gsea_lollipop_download_width, input$gsea_lollipop_download_height)
      plot_obj_to_save <- gsea_lollipop_plot_object()
      width_in <- input$gsea_lollipop_download_width
      height_in <- input$gsea_lollipop_download_height
      dpi_val <- if (input$gsea_lollipop_download_format == "png") input$gsea_lollipop_download_dpi else 300

      # Dynamic height adjustment
      if (!is.null(plot_obj_to_save) && inherits(plot_obj_to_save, "ggplot") && !is.null(plot_obj_to_save$data)) {
        num_pathways_in_plot <- length(unique(plot_obj_to_save$data$DescriptionWrapped))
        if (num_pathways_in_plot > 0) {
          plot_base_height_dl <- 4
          height_per_pathway_dl <- 0.45
          dynamic_height <- plot_base_height_dl + (num_pathways_in_plot * height_per_pathway_dl)
          height_in <- max(5, min(dynamic_height, 30))
        }
      }

      ggsave(file,
        plot = plot_obj_to_save, device = input$gsea_lollipop_download_format,
        width = width_in, height = height_in, units = "in", dpi = dpi_val, bg = "white"
      )
    }
  )
}

# GSEA Enrichment Plot Observer
observe_gsea_enrichment_plot_generation <- function(input, output, session, gsea_results, gsea_enrichment_plot_object, gsea_enrichment_plot_params_for_download, gsea_enrichment_selected_pathway_id) {
  # Update pathway color UI
  output$pathway_color_inputs <- renderUI({
    req(input$gsea_enrichment_plot_pathway_selector)
    render_pathway_color_inputs(input$gsea_enrichment_plot_pathway_selector)
  })

  # Update highlight gene selector
  observeEvent(input$gsea_enrichment_plot_pathway_selector, {
    req(gsea_results(), input$gsea_enrichment_plot_result_type, input$gsea_enrichment_plot_pathway_selector)

    gsea_res_list <- gsea_results()
    s4_object_to_plot <- NULL

    if (input$gsea_enrichment_plot_result_type == "positive") {
      s4_object_to_plot <- gsea_res_list$gsea_object_positive
    } else if (input$gsea_enrichment_plot_result_type == "negative") {
      s4_object_to_plot <- gsea_res_list$gsea_object_negative
    }

    if (!is.null(s4_object_to_plot)) {
      # Extract core genes from pathway
      selected_pathways <- input$gsea_enrichment_plot_pathway_selector
      all_core_genes <- c()

      for (pathway_id in selected_pathways) {
        tryCatch(
          {
            pathway_genes <- s4_object_to_plot@geneSets[[pathway_id]]
            if (!is.null(pathway_genes)) {
              all_core_genes <- c(all_core_genes, pathway_genes)
            }
          },
          error = function(e) {
            # Ignore error, continue processing other pathways
          }
        )
      }

      # Deduplicate and sort
      all_core_genes <- unique(all_core_genes)
      all_core_genes <- sort(all_core_genes)

      # Update selector
      updateSelectizeInput(session, "gsea_enrichment_highlight_genes",
        choices = all_core_genes,
        selected = NULL
      )
    }
  })

  # Main Plot Generation Observer
  observeEvent(input$run_gsea_enrichment_plot, {
    req(gsea_results(), input$gsea_enrichment_plot_result_type, input$gsea_enrichment_plot_pathway_selector)

    selected_pathways <- input$gsea_enrichment_plot_pathway_selector
    if (is.null(selected_pathways) || length(selected_pathways) == 0) {
      showNotification("Please select at least one pathway to plot.", type = "warning")
      return()
    }

    gsea_enrichment_selected_pathway_id(selected_pathways)

    output$gsea_enrichment_plot_status <- renderPrint({
      "Generating GSEA Pathway Enrichment Plot..."
    })
    gsea_enrichment_plot_object(NULL)
    shinyjs::hide("gsea_enrichment_plot_download_options")

    gsea_res_list <- gsea_results()
    s4_object_to_plot <- NULL

    if (input$gsea_enrichment_plot_result_type == "positive") {
      s4_object_to_plot <- gsea_res_list$gsea_object_positive
    } else if (input$gsea_enrichment_plot_result_type == "negative") {
      s4_object_to_plot <- gsea_res_list$gsea_object_negative
    } else {
      output$gsea_enrichment_plot_status <- renderPrint({
        "Error: Invalid GSEA result type selected."
      })
      return()
    }

    if (is.null(s4_object_to_plot)) {
      output$gsea_enrichment_plot_status <- renderPrint({
        "Error: S4 Object not found for selected GSEA result type."
      })
      return()
    }

    # Collect pathway colors (only if multiple pathways selected)
    pathway_colors <- NULL
    if (length(selected_pathways) > 1) {
      pathway_colors <- sapply(1:length(selected_pathways), function(i) {
        color_input_id <- paste0("pathway_color_", i)
        color_value <- input[[color_input_id]]
        if (is.null(color_value) || color_value == "") {
          # Use default colors if input is empty
          if (length(selected_pathways) <= 8) {
            RColorBrewer::brewer.pal(max(3, length(selected_pathways)), "Set2")[i]
          } else {
            rainbow(length(selected_pathways))[i]
          }
        } else {
          color_value
        }
      })

      # Validate color inputs
      invalid_colors <- !grepl("^#[0-9A-Fa-f]{6}$", pathway_colors)
      if (any(invalid_colors)) {
        showNotification("Invalid hex code in pathway colors. Please check settings.", type = "error")
        return()
      }
    }

    # Get highlight genes
    highlighted_genes <- NULL
    if (isTRUE(input$gsea_enrichment_enable_gene_highlight) &&
      !is.null(input$gsea_enrichment_highlight_genes) &&
      length(input$gsea_enrichment_highlight_genes) > 0) {
      highlighted_genes <- input$gsea_enrichment_highlight_genes
    }

    # Fixed highlight colors
    highlight_colors <- c("#0E6DB3", "#BB1E38") # Blue and Red

    tryCatch(
      {
        plot_obj <- generate_gsea_multi_pathway_plot(
          gsea_s4_object = s4_object_to_plot,
          selected_pathway_ids = selected_pathways,
          pathway_colors = pathway_colors,
          highlighted_genes = highlighted_genes,
          highlight_colors = highlight_colors,
          term_width = input$gsea_enrichment_term_width,
          legend_position = c(input$gsea_enrichment_legend_x, input$gsea_enrichment_legend_y),
          subplot_type = as.numeric(input$gsea_enrichment_subplot_type),
          add_pval = input$gsea_enrichment_add_pval,
          pval_x = input$gsea_enrichment_pval_x,
          pval_y = input$gsea_enrichment_pval_y,
          base_font_size = input$gsea_enrichment_base_font
        )

        gsea_enrichment_plot_object(plot_obj)
        gsea_enrichment_plot_params_for_download(list(
          pathways = selected_pathways,
          colors = pathway_colors,
          highlight_genes = highlighted_genes,
          highlight_colors = highlight_colors,
          prefix = paste0("GSEA_", input$gsea_enrichment_plot_result_type), # Added prefix
          pathway_desc = if (length(selected_pathways) == 1) selected_pathways[1] else "Multiple_Pathways" # Added pathway_desc
        ))

        pathway_count <- length(selected_pathways)
        gene_count <- if (!is.null(highlighted_genes)) length(highlighted_genes) else 0

        status_msg <- paste0("GSEA Pathway Enrichment Plot generated (", pathway_count, " pathways)")
        if (gene_count > 0) {
          status_msg <- paste0(status_msg, ", ", gene_count, " highlighted genes")
        }
        status_msg <- paste0(status_msg, ") - ", Sys.time())

        output$gsea_enrichment_plot_status <- renderPrint({
          status_msg
        })
        shinyjs::show("gsea_enrichment_plot_download_options")
        updateTabsetPanel(session, "main_results_tabs", selected = "GSEA Enrichment Plot")
      },
      error = function(e) {
        error_msg_gsea_enrich <- paste("Error generating GSEA Pathway Enrichment Plot:\n", e$message)
        output$gsea_enrichment_plot_status <- renderPrint({
          error_msg_gsea_enrich
        })
        gsea_enrichment_plot_object(NULL)
        showModal(modalDialog(
          title = "GSEA Pathway Plot Failed",
          HTML(paste0("Error Details: <br><pre style='white-space: pre-wrap;'>", e$message, "</pre>"))
        ))
      }
    )
  })
}

# --- Main Analysis Observers ---

# GSEA Analysis Observer
observe_gsea_analysis <- function(input, output, session, data_processing_results, gsea_results) {
  observeEvent(input$run_gsea_analysis, {
    req(input$gmt_file_upload)

    # Use processed data
    req(data_processing_results())
    req(data_processing_results()$params$gene_col)
    gene_summary_data_for_gsea <- data_processing_results()$gene_summary_data
    gene_id_col_for_gsea <- data_processing_results()$params$gene_col

    # Validate data
    if (is.null(gene_summary_data_for_gsea) || nrow(gene_summary_data_for_gsea) == 0) {
      showModal(modalDialog(title = "Error", "Data source empty. Please check data loading."))
      return()
    }

    gsea_start_abs_val <- 0.51
    withProgress(message = "Running GSEA...", value = 0.5, session = session, {
      shiny::setProgress(value = gsea_start_abs_val, detail = "Initializing GSEA...", message = "Running GSEA...")
      output$status_output_gsea <- renderPrint({
        "Running GSEA... (using processed data)"
      })
      gsea_results(NULL)
      output$gsea_results_display_ui <- renderUI({
        NULL
      })
      tryCatch(
        {
          gsea_res_list <- perform_gsea_analysis_clusterProfiler_shiny(
            gene_summary_df = gene_summary_data_for_gsea,
            gmt_file_path = input$gmt_file_upload$datapath,
            gene_id_column = gene_id_col_for_gsea,
            clp_n_permutations = input$gsea_clp_n_perm,
            clp_pvalue_cutoff = input$gsea_clp_pval_cutoff,
            clp_min_gene_set_size = input$gsea_clp_min_gs_size,
            max_gene_set_size = input$gsea_max_gs_size,
            clp_seed_logical = input$gsea_clp_seed,
            clp_by = "fgsea",
            shiny_session = session,
            initial_progress_value_abs = gsea_start_abs_val
          )
          gsea_results(gsea_res_list)

          # --- Tick 1: Immediate Updates (Status, Progress, Tab Switch) ---
          shiny::setProgress(1.0, message = "GSEA Completed!", detail = "Rendering results...")

          output$status_output_gsea <- renderPrint({
            paste0(
              "GSEA Completed (using processed data)\n",
              "Number of Genes: ", nrow(gene_summary_data_for_gsea), "\n",
              "Gene ID Column: ", gene_id_col_for_gsea
            )
          })

          # Switch tab immediately
          updateTabsetPanel(session, "main_results_tabs", selected = "GSEA Analysis Results")

          # --- Tick 2: Render Results Table (Delayed 100ms) ---
          later::later(function() {
            output$gsea_results_display_ui <- renderUI({
              render_gsea_results_tables_ui(gsea_res_list)
            })
          }, delay = 0.1)

          # --- Tick 3: Render Lollipop Params (Delayed 200ms) ---
          later::later(function() {
            output$gsea_lollipop_params_ui_placeholder <- renderUI({
              render_gsea_lollipop_params_ui(gsea_res_list)
            })
          }, delay = 0.2)

          # --- Tick 4: Render Enrichment Params UI (Delayed 300ms) ---
          later::later(function() {
            output$gsea_enrichment_plot_params_ui_placeholder <- renderUI({
              render_gsea_enrichment_params_ui(gsea_res_list)
            })
          }, delay = 0.3)

          # --- Tick 5: Render sgRNA Paired Params UI (Delayed 400ms) ---
          later::later(function() {
            output$sgrna_paired_plot_params_ui_placeholder <- renderUI({
              results <- data_processing_results()
              render_sgrna_params_ui_with_gsea(results, gsea_res_list)
            })
          }, delay = 0.4)
        },
        error = function(e) {
          output$status_output_gsea <- renderPrint({
            paste("GSEA Failed:", e$message)
          })
          gsea_results(NULL)
          showModal(modalDialog(title = "GSEA Failed", HTML(paste0("Error Details: <br><pre style='white-space: pre-wrap;'>", e$message, "</pre>"))))
        }
      )
    })
  })
}

# GSEA Enrichment Download Handler
create_gsea_enrichment_download_handler <- function(input, gsea_enrichment_plot_object, gsea_enrichment_plot_params_for_download) {
  downloadHandler(
    filename = function() {
      req(gsea_enrichment_plot_params_for_download())
      dl_params <- gsea_enrichment_plot_params_for_download()
      safe_desc <- gsub("[^a-zA-Z0-9_]", "_", dl_params$pathway_desc)
      safe_desc <- substr(safe_desc, 1, 50)
      safe_prefix <- gsub("[^a-zA-Z0-9_]", "_", dl_params$prefix)
      paste0(safe_prefix, "_", safe_desc, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$gsea_enrichment_download_format)
    },
    content = function(file) {
      req(gsea_enrichment_plot_object(), input$gsea_enrichment_download_width, input$gsea_enrichment_download_height)
      ggsave(file,
        plot = gsea_enrichment_plot_object(), device = input$gsea_enrichment_download_format,
        width = input$gsea_enrichment_download_width, height = input$gsea_enrichment_download_height,
        units = "in", dpi = if (input$gsea_enrichment_download_format == "png") input$gsea_enrichment_download_dpi else 300, bg = "white"
      )
    }
  )
}

# --- Selector Update Observers ---

# GSEA Pathway Selector Observer
observe_gsea_enrichment_pathway_selector_update <- function(input, session, gsea_results) {
  observeEvent(input$gsea_enrichment_plot_result_type,
    {
      req(gsea_results(), input$gsea_enrichment_plot_result_type)
      gsea_res_list <- gsea_results()
      pathway_choices_df <- NULL
      if (input$gsea_enrichment_plot_result_type == "positive" && !is.null(gsea_res_list$gsea_results_positive_df)) {
        pathway_choices_df <- gsea_res_list$gsea_results_positive_df
      } else if (input$gsea_enrichment_plot_result_type == "negative" && !is.null(gsea_res_list$gsea_results_negative_df)) {
        pathway_choices_df <- gsea_res_list$gsea_results_negative_df
      }

      if (!is.null(pathway_choices_df) && nrow(pathway_choices_df) > 0 && "ID" %in% names(pathway_choices_df) && "Description" %in% names(pathway_choices_df)) {
        pathway_choices_df_sorted <- NULL
        if (input$gsea_enrichment_plot_result_type == "positive") {
          pathway_choices_df_sorted <- pathway_choices_df %>% arrange(desc(NES), p.adjust)
        } else if (input$gsea_enrichment_plot_result_type == "negative") {
          pathway_choices_df_sorted <- pathway_choices_df %>% arrange(NES, p.adjust)
        } else {
          pathway_choices_df_sorted <- pathway_choices_df %>% arrange(p.adjust, desc(abs(NES)))
        }

        choices_list <- setNames(pathway_choices_df_sorted$ID, pathway_choices_df_sorted$Description)
        updateSelectizeInput(session, "gsea_enrichment_plot_pathway_selector",
          choices = choices_list,
          selected = if (length(choices_list) > 0) choices_list[1] else NULL,
          server = TRUE
        )
      } else {
        updateSelectizeInput(session, "gsea_enrichment_plot_pathway_selector", choices = list("No available pathways" = ""), selected = "")
      }
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
  )
}

# sgRNA GSEA Type Selector Observer
observe_sgrna_gsea_type_selector_update <- function(input, session, gsea_results) {
  observeEvent(input$sgrna_paired_gsea_type_selector,
    {
      req(gsea_results(), input$sgrna_paired_gsea_type_selector)
      selected_gsea_type <- input$sgrna_paired_gsea_type_selector
      if (is.null(selected_gsea_type) || selected_gsea_type == "") {
        updateSelectizeInput(session, "sgrna_paired_pathway_selector", choices = list("Please select a valid GSEA result source..." = ""), selected = "", server = TRUE)
        updateSelectizeInput(session, "sgrna_paired_gene_selector", choices = list("Please select a pathway first..." = ""), selected = "", server = TRUE)
        return()
      }

      gsea_res_list <- gsea_results()
      pathway_choices_df <- NULL

      if (selected_gsea_type == "positive" && !is.null(gsea_res_list$gsea_results_positive_df)) {
        pathway_choices_df <- gsea_res_list$gsea_results_positive_df
      } else if (selected_gsea_type == "negative" && !is.null(gsea_res_list$gsea_results_negative_df)) {
        pathway_choices_df <- gsea_res_list$gsea_results_negative_df
      }

      pathway_display_choices <- list("Please select a pathway..." = "")
      if (!is.null(pathway_choices_df) && nrow(pathway_choices_df) > 0 && "ID" %in% names(pathway_choices_df) && "Description" %in% names(pathway_choices_df)) {
        pathway_choices_df_sorted <- pathway_choices_df %>% arrange(p.adjust, desc(abs(NES)))
        pathway_display_choices <- c(pathway_display_choices, setNames(pathway_choices_df_sorted$ID, pathway_choices_df_sorted$Description))
      }
      updateSelectizeInput(session, "sgrna_paired_pathway_selector",
        choices = pathway_display_choices,
        selected = "",
        server = TRUE
      )
      updateSelectizeInput(session, "sgrna_paired_gene_selector", choices = list("Please select a pathway first..." = ""), selected = "", server = TRUE)
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
  )
}

# sgRNA Pathway Selector Observer
observe_sgrna_pathway_selector_update <- function(input, session, gsea_results, data_processing_results) {
  observeEvent(input$sgrna_paired_pathway_selector,
    {
      req(gsea_results(), input$sgrna_paired_gsea_type_selector, input$sgrna_paired_pathway_selector)

      selected_pathway_id <- input$sgrna_paired_pathway_selector
      if (is.null(selected_pathway_id) || selected_pathway_id == "") {
        updateSelectizeInput(session, "sgrna_paired_gene_selector", choices = list("Please select a pathway first..." = ""), selected = "", server = TRUE)
        return()
      }

      gsea_res_list <- gsea_results()
      selected_gsea_type <- input$sgrna_paired_gsea_type_selector
      if (is.null(selected_gsea_type) || selected_gsea_type == "") {
        return()
      }

      source_gsea_df <- NULL
      if (selected_gsea_type == "positive" && !is.null(gsea_res_list$gsea_results_positive_df)) {
        source_gsea_df <- gsea_res_list$gsea_results_positive_df
      } else if (selected_gsea_type == "negative" && !is.null(gsea_res_list$gsea_results_negative_df)) {
        source_gsea_df <- gsea_res_list$gsea_results_negative_df
      } else {
        updateSelectizeInput(session, "sgrna_paired_gene_selector", choices = list("Invalid GSEA result type" = ""), selected = "", server = TRUE)
        return()
      }

      gene_choices_list <- list("No core genes or no pathway selected" = "")
      if (!is.null(source_gsea_df) && selected_pathway_id %in% source_gsea_df$ID) {
        core_enrich_str <- source_gsea_df[source_gsea_df$ID == selected_pathway_id, "core_enrichment", drop = TRUE]
        if (length(core_enrich_str) == 1 && !is.na(core_enrich_str) && nzchar(core_enrich_str)) {
          core_genes <- unique(strsplit(core_enrich_str, "/")[[1]])
          core_genes <- core_genes[nzchar(core_genes)]
          if (length(core_genes) > 0) {
            norm_counts_genes <- data_processing_results()$normalized_counts[[data_processing_results()$params$gene_col]]
            plottable_core_genes <- intersect(core_genes, unique(norm_counts_genes))

            if (length(plottable_core_genes) > 0) {
              gene_choices_list <- setNames(plottable_core_genes, plottable_core_genes)
            } else {
              gene_choices_list <- list("Pathway core genes have no plottable sgRNA in data" = "")
            }
          } else {
            gene_choices_list <- list("No core genes in pathway" = "")
          }
        } else {
          gene_choices_list <- list("Cannot retrieve core genes for this pathway" = "")
        }
      }
      updateSelectizeInput(session, "sgrna_paired_gene_selector",
        choices = gene_choices_list,
        selected = if (length(gene_choices_list) > 0 && gene_choices_list[[1]] != "" && !grepl("No core genes|Cannot retrieve|Please select|Invalid GSEA|no plottable sgRNA", names(gene_choices_list)[1], ignore.case = TRUE)) gene_choices_list[[1]] else "",
        server = TRUE
      )
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
  )
}



# --- Tab Switching Observer for Sidebar Visibility ---
observe_tab_switching <- function(input, session) {
  observeEvent(input$main_results_tabs, {
    req(input$main_results_tabs)
    tab <- input$main_results_tabs

    # Hide all first
    shinyjs::hide("data_processing_params_box")
    shinyjs::hide("gsea_parameters_section")
    shinyjs::hide("gsea_lollipop_params_box")
    shinyjs::hide("gsea_enrichment_plot_params_box")
    shinyjs::hide("sgrna_paired_plot_params_box")
    shinyjs::hide("gene_vis_params_box")

    # Show based on tab
    if (tab == "Data Processing Output") {
      shinyjs::show("data_processing_params_box")
    } else if (tab == "GSEA Analysis Results") {
      shinyjs::show("gsea_parameters_section")
    } else if (tab == "GSEA Lollipop Plot") {
      shinyjs::show("gsea_lollipop_params_box")
    } else if (tab == "GSEA Enrichment Plot") {
      shinyjs::show("gsea_enrichment_plot_params_box")
    } else if (tab == "sgRNA Paired Plot") {
      shinyjs::show("sgrna_paired_plot_params_box")
    } else if (tab == "Gene Visualization") {
      shinyjs::show("gene_vis_params_box")
    }
  })
}

# --- Server-Side Selectize Updaters ---

# Updater for sgRNA Paired Gene Selector (Direct)
observe_sgrna_gene_selector_update <- function(input, session, data_processing_results) {
  observe({
    req(data_processing_results())
    gene_summary <- data_processing_results()$gene_summary_data
    gene_col <- data_processing_results()$params$gene_col
    req(gene_summary, gene_col)

    all_genes <- unique(gene_summary[[gene_col]])

    # Update for both basic and GSEA-aware UIs (they use same ID if rendered dynamically, but logic handles availability)
    updateSelectizeInput(session, "sgrna_paired_gene_selector_direct",
      choices = all_genes,
      server = TRUE
    )
  })
}

# Updater for Gene Vis Volcano Labels
observe_gene_vis_labels_update <- function(input, session, data_processing_results) {
  observe({
    req(data_processing_results())
    gene_summary <- data_processing_results()$gene_summary_data
    gene_col <- data_processing_results()$params$gene_col
    req(gene_summary, gene_col)

    all_genes <- unique(gene_summary[[gene_col]])

    updateSelectizeInput(session, "gene_vis_volcano_labels",
      choices = all_genes,
      server = TRUE
    )
  })
}
