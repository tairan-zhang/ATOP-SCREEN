# SPDX-License-Identifier: GPL-3.0-or-later
# cpp_progress_wrapper.R
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

perform_cpp_permutation_with_progress <- function(
    sgrna_data, gene_summary_data = NULL, gene_col = "Gene", score_col = "diff_score",
    n_permutations = 1000L, min_sgrna_threshold = 3L, n_threads = NULL,
    progress_callback = NULL, shiny_session = NULL, seed = NULL, batch_size = NULL
) {
    report <- if (is.null(progress_callback) && is.null(shiny_session)) NULL else function(message) {
        if (!is.null(progress_callback)) progress_callback(message)
        if (!is.null(shiny_session)) {
            percentage <- as.numeric(sub(".*\\(([0-9.]+)%\\).*", "\\1", message))
            shiny::setProgress(
                value = min(percentage / 100, 0.99),
                message = message, session = shiny_session
            )
        }
    }
    perform_cpp_permutation(
        sgrna_data, gene_summary_data, gene_col, score_col,
        n_permutations, min_sgrna_threshold, n_threads, report, seed, batch_size
    )
}
