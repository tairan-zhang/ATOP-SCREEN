app <- source("app.R")$value
stopifnot(inherits(app, "shiny.appobj"))
stopifnot(grepl("ATOP-SCREEN", as.character(ui), fixed = TRUE),
          grepl("Version 0.0.1", as.character(ui), fixed = TRUE))

input_path <- tempfile(fileext = ".csv")
write.csv(
    data.frame(gRNA = paste0("s", 1:12), Gene = rep(c("A", "B", "C"), each = 4),
               treatment = c(12, 14, 16, 18, 3, 4, 5, 6, 8, 9, 10, 11), control = rep(10, 12)),
    input_path, row.names = FALSE
)
shiny::testServer(server, {
    session$setInputs(workspace_area = "analysis", analysis_section = "screen")
    stopifnot(active_page() == "Data Processing Output")
    session$setInputs(analysis_section = "gsea")
    stopifnot(active_page() == "GSEA Analysis Results")
    for (level in c("sgrna", "gene", "pathway")) {
        session$setInputs(workspace_area = "results", result_level = level)
        stopifnot(active_page() == switch(level, sgrna = "sgRNA Results", gene = "Gene Results", pathway = "Pathway Results"))
    }
    session$setInputs(workspace_area = "visualization", visualization_level = "sgrna")
    stopifnot(active_page() == "sgRNA Paired Plot")
    session$setInputs(visualization_level = "gene", gene_vis_plot_type = "ranking")
    stopifnot(active_page() == "Gene Visualization")
    session$setInputs(visualization_level = "pathway", pathway_plot_menu = "lollipop")
    stopifnot(active_page() == "GSEA Lollipop Plot")
    session$setInputs(pathway_plot_menu = "enrichment")
    stopifnot(active_page() == "GSEA Enrichment Plot")
    session$setInputs(workspace_area = "analysis", analysis_section = "screen")
    session$setInputs(raw_file_upload = data.frame(
        name = "counts.csv", size = file.info(input_path)$size,
        type = "text/csv", datapath = input_path
    ))
    session$setInputs(
        grna_col_selector = "gRNA", gene_col_selector = "Gene", sequence_col_selector = "",
        cond1_reps_selector = "treatment", cond2_reps_selector = "control",
        analysis_engine = "atop", permutation_engine_selector = "cpp",
        N_perm_input = 1003, pseudo_count_lfc_input = 1, min_sgrna_threshold = 3,
        run_data_processing = 0
    )
    session$setInputs(run_data_processing = 1)
    result <- data_processing_results()
    stopifnot(!is.null(result), nrow(result$gene_summary_data) == 3)
    stopifnot(all(is.finite(result$gene_summary_data$P_positive)))
    stopifnot(identical(result$params$software_version, "0.0.1"))
    session$setInputs(workspace_area = "visualization", visualization_level = "sgrna")
    sidebar <- output$sgrna_paired_plot_params_ui_placeholder$html
    stopifnot(!grepl('value="A"', sidebar, fixed = TRUE))
    large_result <- result
    large_result$gene_summary_data <- data.frame(Gene = paste0("GENE", seq_len(30000)))
    stopifnot(nchar(as.character(render_sgrna_params_ui_basic(large_result))) < 30000)
    session$setInputs(sgrna_paired_gene_selector_direct_ready = 1)
    gmt_path <- tempfile(fileext = ".gmt")
    writeLines(c("SET_A\tdescription\tA\tB", "SET_B\tdescription\tB\tC\tD"), gmt_path)
    session$setInputs(gmt_file_upload = data.frame(name = "test.gmt", size = file.info(gmt_path)$size, datapath = gmt_path))
    stopifnot(gmt_summary()$sets == 2, gmt_summary()$genes == 4, gmt_summary()$median == 2.5)
    unlink(gmt_path)
    session$setInputs(
        sgrna_gene_selection_method = "direct", sgrna_paired_gene_selector_direct = "A",
        sgrna_paired_condition1_col = "treatment", sgrna_paired_condition2_col = "control",
        sgrna_paired_condition1_label = "", sgrna_paired_condition2_label = "",
        sgrna_paired_color_cond1 = "#104e8b", sgrna_paired_color_cond2 = "#b22222",
        sgrna_paired_line_color = "#afc3d8", sgrna_paired_line_size = 0.5,
        sgrna_paired_point_size = 2.5, sgrna_paired_base_font_size = 8,
        sgrna_paired_plot_title = "", sgrna_paired_y_axis_label = "Normalized count",
        run_sgrna_paired_plot = 0
    )
    session$setInputs(run_sgrna_paired_plot = 1)
    stopifnot(inherits(sgrna_paired_plot_object(), "ggplot"))
    stopifnot(plot_dimensions$sgrna_paired()$width == 8.5)
    session$setInputs(sgrna_paired_download_width = 6, sgrna_paired_download_height = 9,
        sgrna_paired_download_dpi = 150, sgrna_paired_apply_preview = 0)
    stopifnot(plot_dimensions$sgrna_paired()$width == 8.5)
    session$setInputs(sgrna_paired_apply_preview = 1)
    stopifnot(identical(plot_dimensions$sgrna_paired(), list(width = 6, height = 9, dpi = 150)))
    session$setInputs(sgrna_paired_download_width = 0, sgrna_paired_apply_preview = 2)
    stopifnot(plot_dimensions$sgrna_paired()$width == 6)

    session$setInputs(workspace_area = "analysis", run_data_processing = 2)
    stopifnot(is.null(sgrna_paired_plot_object()), is.null(sgrna_paired_plot_params_for_download()))
    stopifnot(nchar(output$gene_summary_table) > 0)
    stopifnot(grepl("sg_data_table", output$sgrna_results_ui$html, fixed = TRUE))
    stopifnot(!grepl("gene_summary_table", output$sgrna_results_ui$html, fixed = TRUE))
    stopifnot(grepl("gene_summary_table", output$gene_results_ui$html, fixed = TRUE))
})
unlink(input_path)

check_http <- function() {
    port <- httpuv::randomPort()
    process <- processx::process$new(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", "-e", sprintf("shiny::runApp('.', host='127.0.0.1', port=%d, launch.browser=FALSE)", port)),
        stdout = "|", stderr = "|"
    )
    on.exit(process$kill(), add = TRUE)
    page <- NULL
    for (attempt in seq_len(60L)) {
        if (!process$is_alive()) stop(process$read_all_error())
        page <- tryCatch(
            suppressWarnings(readLines(sprintf("http://127.0.0.1:%d", port), warn = FALSE)),
            error = function(e) NULL
        )
        if (length(page)) break
        Sys.sleep(0.5)
    }
    stopifnot(any(grepl("ATOP-SCREEN-0.0.1", page, fixed = TRUE)))
}
check_http()
cat("PASS: Shiny app creation, upload-to-results server flow, DataTable rendering, and HTTP startup.\n")
