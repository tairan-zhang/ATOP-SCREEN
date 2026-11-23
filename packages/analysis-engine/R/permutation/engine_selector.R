























engine_selector <- function(
    sgrna_info_dt,
    gene_summary_dt,
    gene_col,
    N_perm = 1000,
    n_cores = NULL,
    progress_callback = NULL,
    user_engine_choice = "cpp",
    min_sgrna_threshold = 3, seed = 42L) {

    if (!user_engine_choice %in% c("cpp", "r_parallel")) {
        user_engine_choice <- "cpp"
        if (!is.null(progress_callback)) {
            progress_callback("Unknown engine choice, defaulting to C++ engine")
        }
    }

    engine_name <- switch(user_engine_choice,
        "cpp" = "C++17 Engine",
        "r_parallel" = "R Parallel Engine"
    )

    if (!is.null(progress_callback)) {
        progress_callback(paste("User Selection:", engine_name))
    }


    tryCatch(
        {
            if (user_engine_choice == "cpp") {

                if (!exists("perform_cpp_permutation")) {
                    if (!is.null(progress_callback)) {
                        progress_callback("Initializing C++ Engine for the first time...")
                    }


                    init_success <- tryCatch({
                        if (exists("initialize_cpp_engine")) {
                            initialize_cpp_engine()
                        } else {
                            FALSE
                        }
                    })

                    if (!init_success) {
                        if (!is.null(progress_callback)) {
                            progress_callback("C++ Engine unavailable, automatically switching to R Parallel Engine")
                        }
                        user_engine_choice <- "r_parallel"
                    }
                }

                if (user_engine_choice == "cpp" && exists("perform_cpp_permutation")) {
                    return(perform_cpp_permutation(
                        sgrna_data = sgrna_info_dt,
                        gene_summary_data = gene_summary_dt,
                        gene_col = gene_col,
                        n_permutations = N_perm,
                        min_sgrna_threshold = min_sgrna_threshold, seed = seed,
                        n_threads = n_cores,
                        progress_callback = progress_callback
                    ))
                }
            }


            if (!is.null(progress_callback) && user_engine_choice == "r_parallel") {
                progress_callback("Using R Parallel Engine")
            }

            return(perform_permutation(
                sgrna_info_dt = sgrna_info_dt,
                gene_summary_dt = gene_summary_dt,
                gene_col = gene_col,
                N_perm = N_perm,
                n_cores = n_cores,
                use_data_table = TRUE,
                progress_callback = progress_callback,
                min_sgrna_threshold = min_sgrna_threshold, seed = seed
            ))
        },
        error = function(e) {
            if (!is.null(progress_callback)) {
                progress_callback(paste("Engine execution failed, using R Parallel Engine:", e$message))
            }


            return(perform_permutation(
                sgrna_info_dt = sgrna_info_dt,
                gene_summary_dt = gene_summary_dt,
                gene_col = gene_col,
                N_perm = N_perm,
                n_cores = n_cores,
                use_data_table = TRUE,
                progress_callback = progress_callback,
                min_sgrna_threshold = min_sgrna_threshold, seed = seed
            ))
        }
    )
}
