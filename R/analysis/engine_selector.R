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

engine_selector <- function(
    sgrna_info_dt, gene_summary_dt, gene_col, N_perm = 1000L, n_cores = NULL,
    progress_callback = NULL, user_engine_choice = "cpp", min_sgrna_threshold = 3L,
    shiny_session = NULL, seed = NULL
) {
    if (!user_engine_choice %in% c("cpp", "r_parallel")) user_engine_choice <- "cpp"
    if (user_engine_choice == "cpp") {
        result <- tryCatch(
            perform_cpp_permutation_with_progress(
                sgrna_data = sgrna_info_dt, gene_summary_data = gene_summary_dt,
                gene_col = gene_col, n_permutations = N_perm,
                min_sgrna_threshold = min_sgrna_threshold, n_threads = n_cores,
                progress_callback = progress_callback, shiny_session = shiny_session,
                seed = seed
            ),
            error = function(e) {
                message("C++ engine failed; using R: ", conditionMessage(e))
                NULL
            }
        )
        if (!is.null(result)) return(result)
    }
    perform_permutation(
        sgrna_info_dt, gene_summary_dt, gene_col, N_perm, n_cores,
        progress_callback = progress_callback, min_sgrna_threshold = min_sgrna_threshold,
        seed = seed
    )
}

select_permutation_engine <- engine_selector
