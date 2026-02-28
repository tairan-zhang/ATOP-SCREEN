# SPDX-License-Identifier: GPL-3.0-or-later
# engine_selector.R
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


# Engine Selector
# Selects Reference Engine (C++ vs R Parallel)

#' Engine Selector
engine_selector <- function(
    sgrna_info_dt,
    gene_summary_dt,
    gene_col,
    N_perm = 1000,
    n_cores = NULL,
    progress_callback = NULL,
    user_engine_choice = "cpp",
    min_sgrna_threshold = 3) {
    # Validate engine choice
    if (!user_engine_choice %in% c("cpp", "r_parallel")) {
        user_engine_choice <- "cpp"
        if (!is.null(progress_callback)) {
            progress_callback("Unknown engine choice, defaulting to C++ engine")
        }
    }

    engine_name <- switch(user_engine_choice,
        "cpp" = "C++ Ultra-Fast Engine",
        "r_parallel" = "R Parallel Engine"
    )

    if (!is.null(progress_callback)) {
        progress_callback(paste("User Selection:", engine_name))
    }

    # Execute Engine
    tryCatch(
        {
            if (user_engine_choice == "cpp") {
                # Attempt to use C++ engine
                if (!exists("perform_cpp_permutation")) {
                    if (!is.null(progress_callback)) {
                        progress_callback("Initializing C++ Engine for the first time...")
                    }

                    # Initialize C++ Engine
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
                    # Use C++ Engine (w/ progress)
                    if (exists("perform_cpp_permutation_with_progress")) {
                        # Check for Shiny session
                        shiny_session <- NULL
                        if (exists(".shiny_session", envir = .GlobalEnv)) {
                            shiny_session <- get(".shiny_session", envir = .GlobalEnv)
                        }

                        return(perform_cpp_permutation_with_progress(
                            sgrna_data = sgrna_info_dt,
                            gene_summary_data = gene_summary_dt,
                            gene_col = gene_col,
                            n_permutations = N_perm,
                            min_sgrna_threshold = min_sgrna_threshold,
                            n_threads = n_cores,
                            progress_callback = progress_callback,
                            shiny_session = shiny_session
                        ))
                    } else {
                        # Fallback to basic C++ engine
                        return(perform_cpp_permutation(
                            sgrna_data = sgrna_info_dt,
                            gene_summary_data = gene_summary_dt,
                            gene_col = gene_col,
                            n_permutations = N_perm,
                            min_sgrna_threshold = min_sgrna_threshold,
                            n_threads = n_cores,
                            progress_callback = progress_callback
                        ))
                    }
                }
            }

            # Fallback to R Parallel Engine
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
                min_sgrna_threshold = min_sgrna_threshold
            ))
        },
        error = function(e) {
            if (!is.null(progress_callback)) {
                progress_callback(paste("Engine execution failed, using R Parallel Engine:", e$message))
            }

            # Final fallback
            return(perform_permutation(
                sgrna_info_dt = sgrna_info_dt,
                gene_summary_dt = gene_summary_dt,
                gene_col = gene_col,
                N_perm = N_perm,
                n_cores = n_cores,
                use_data_table = TRUE,
                progress_callback = progress_callback,
                min_sgrna_threshold = min_sgrna_threshold
            ))
        }
    )
}
