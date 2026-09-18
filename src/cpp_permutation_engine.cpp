// SPDX-License-Identifier: GPL-3.0-or-later
// cpp_permutation_engine.cpp
// Part of ATOP-SCREEN
//
// Copyright (C) 2026 Tairan Zhang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.


#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <exception>
#include <mutex>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]

namespace {

struct ScorePair {
    double positive;
    double negative;
};

struct PermutationWorkspace {
    std::vector<int> shuffled_gene_indices;
    std::vector<std::size_t> cursors;
    std::vector<double> grouped_scores;

    PermutationWorkspace(std::size_t n_sgrnas, std::size_t n_genes)
        : shuffled_gene_indices(n_sgrnas),
          cursors(n_genes),
          grouped_scores(n_sgrnas) {}
};

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

std::uint64_t bounded_random(std::mt19937_64& rng, std::uint64_t bound) {
    const std::uint64_t threshold = static_cast<std::uint64_t>(-bound) % bound;
    while (true) {
        const std::uint64_t value = rng();
        if (value >= threshold) {
            return value % bound;
        }
    }
}

void deterministic_shuffle(std::vector<int>& values, std::mt19937_64& rng) {
    for (std::size_t remaining = values.size(); remaining > 1; --remaining) {
        const std::size_t selected = static_cast<std::size_t>(bounded_random(rng, remaining));
        std::swap(values[remaining - 1], values[selected]);
    }
}

template <typename Iterator>
ScorePair calculate_score_pair(Iterator begin, Iterator end) {
    const std::size_t n = static_cast<std::size_t>(std::distance(begin, end));
    std::sort(begin, end);
    const std::size_t k = (2U * n + 2U) / 3U;

    double negative_sum = 0.0;
    for (std::size_t index = 0; index < k; ++index) {
        negative_sum += *(begin + static_cast<std::ptrdiff_t>(index));
    }

    double positive_sum = 0.0;
    for (std::size_t index = 0; index < k; ++index) {
        positive_sum += *(end - 1 - static_cast<std::ptrdiff_t>(index));
    }

    return {
        positive_sum / static_cast<double>(k),
        negative_sum / static_cast<double>(k)
    };
}

int resolve_thread_count(int requested_threads, int n_permutations) {
    unsigned int available = std::thread::hardware_concurrency();
    if (available == 0U) {
        available = 1U;
    }

    int resolved = requested_threads;
    if (resolved <= 0) {
        resolved = available > 1U ? static_cast<int>(available - 1U) : 1;
    }

    resolved = std::max(1, resolved);
    resolved = std::min(resolved, n_permutations);
    resolved = std::min(resolved, static_cast<int>(available));
    return resolved;
}

std::vector<std::size_t> build_offsets(
    const std::vector<int>& gene_indices,
    int n_genes
) {
    std::vector<std::size_t> offsets(static_cast<std::size_t>(n_genes) + 1U, 0U);
    for (int gene_index : gene_indices) {
        ++offsets[static_cast<std::size_t>(gene_index) + 1U];
    }
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    return offsets;
}

void group_scores(
    const std::vector<double>& diff_scores,
    const std::vector<int>& gene_indices,
    const std::vector<std::size_t>& offsets,
    PermutationWorkspace& workspace
) {
    std::copy(offsets.begin(), offsets.end() - 1, workspace.cursors.begin());
    for (std::size_t index = 0; index < gene_indices.size(); ++index) {
        const std::size_t gene_index = static_cast<std::size_t>(gene_indices[index]);
        workspace.grouped_scores[workspace.cursors[gene_index]++] = diff_scores[index];
    }
}

std::vector<int> map_gene_indices(
    const StringVector& gene_labels,
    const StringVector& unique_genes
) {
    std::unordered_map<std::string, int> gene_to_index;
    gene_to_index.reserve(static_cast<std::size_t>(unique_genes.size()));

    for (R_xlen_t index = 0; index < unique_genes.size(); ++index) {
        if (StringVector::is_na(unique_genes[index])) {
            stop("unique_genes must not contain missing values");
        }
        const std::string gene = as<std::string>(unique_genes[index]);
        const auto inserted = gene_to_index.emplace(gene, static_cast<int>(index));
        if (!inserted.second) {
            stop("unique_genes contains duplicate gene names");
        }
    }

    std::vector<int> gene_indices(static_cast<std::size_t>(gene_labels.size()));
    for (R_xlen_t index = 0; index < gene_labels.size(); ++index) {
        if (StringVector::is_na(gene_labels[index])) {
            stop("gene_labels must not contain missing values");
        }
        const std::string gene = as<std::string>(gene_labels[index]);
        const auto match = gene_to_index.find(gene);
        if (match == gene_to_index.end()) {
            stop("gene_labels contains a gene absent from unique_genes: " + gene);
        }
        gene_indices[static_cast<std::size_t>(index)] = match->second;
    }

    return gene_indices;
}

List calculate_observed_scores(
    const std::vector<double>& diff_scores,
    const std::vector<int>& gene_indices,
    const StringVector& unique_genes,
    int min_sgrna_threshold
) {
    const int n_genes = unique_genes.size();
    const std::vector<std::size_t> offsets = build_offsets(gene_indices, n_genes);
    PermutationWorkspace workspace(diff_scores.size(), static_cast<std::size_t>(n_genes));
    group_scores(diff_scores, gene_indices, offsets, workspace);

    NumericVector observed_positive(n_genes, R_NaN);
    NumericVector observed_negative(n_genes, R_NaN);
    LogicalVector valid_genes(n_genes, false);
    int valid_gene_count = 0;

    for (int gene_index = 0; gene_index < n_genes; ++gene_index) {
        const std::size_t begin_index = offsets[static_cast<std::size_t>(gene_index)];
        const std::size_t end_index = offsets[static_cast<std::size_t>(gene_index) + 1U];
        if (end_index - begin_index < static_cast<std::size_t>(min_sgrna_threshold)) {
            continue;
        }

        const ScorePair scores = calculate_score_pair(
            workspace.grouped_scores.begin() + static_cast<std::ptrdiff_t>(begin_index),
            workspace.grouped_scores.begin() + static_cast<std::ptrdiff_t>(end_index)
        );
        observed_positive[gene_index] = scores.positive;
        observed_negative[gene_index] = scores.negative;
        valid_genes[gene_index] = true;
        ++valid_gene_count;
    }

    return List::create(
        Named("observed_positive_scores") = observed_positive,
        Named("observed_negative_scores") = observed_negative,
        Named("valid_genes") = valid_genes,
        Named("n_valid_genes") = valid_gene_count,
        Named("offsets") = wrap(offsets)
    );
}

}

// [[Rcpp::export]]
List calculate_observed_scores_cpp(
    NumericVector diff_scores,
    StringVector gene_labels,
    StringVector unique_genes,
    int min_sgrna_threshold = 3
) {
    if (diff_scores.size() != gene_labels.size()) {
        stop("diff_scores and gene_labels must have the same length");
    }
    if (min_sgrna_threshold < 1) {
        stop("min_sgrna_threshold must be at least 1");
    }

    const std::vector<double> scores = as<std::vector<double>>(diff_scores);
    for (double score : scores) {
        if (!std::isfinite(score)) {
            stop("diff_scores must contain only finite values");
        }
    }
    const std::vector<int> gene_indices = map_gene_indices(gene_labels, unique_genes);
    const List observed = calculate_observed_scores(
        scores,
        gene_indices,
        unique_genes,
        min_sgrna_threshold
    );

    return List::create(
        Named("observed_pos_scores") = observed["observed_positive_scores"],
        Named("observed_neg_scores") = observed["observed_negative_scores"],
        Named("gene_names") = unique_genes
    );
}

// [[Rcpp::export]]
List perform_cpp_permutation_test(
    NumericVector diff_scores,
    StringVector gene_labels,
    StringVector unique_genes,
    int n_permutations = 1000,
    int min_sgrna_threshold = 3,
    int seed = 42,
    bool show_progress = true,
    int n_threads = 0,
    int permutation_offset = 0
) {
    if (diff_scores.size() != gene_labels.size()) {
        stop("diff_scores and gene_labels must have the same length");
    }
    if (diff_scores.size() == 0 || unique_genes.size() == 0) {
        stop("permutation input must not be empty");
    }
    if (n_permutations < 1) {
        stop("n_permutations must be at least 1");
    }
    if (min_sgrna_threshold < 1) {
        stop("min_sgrna_threshold must be at least 1");
    }
    if (permutation_offset < 0) {
        stop("permutation_offset must not be negative");
    }

    const int n_genes = unique_genes.size();
    const std::vector<double> scores = as<std::vector<double>>(diff_scores);
    for (double score : scores) {
        if (!std::isfinite(score)) {
            stop("diff_scores must contain only finite values");
        }
    }
    const std::vector<int> gene_indices = map_gene_indices(gene_labels, unique_genes);
    const std::vector<std::size_t> offsets = build_offsets(gene_indices, n_genes);
    const List observed = calculate_observed_scores(
        scores,
        gene_indices,
        unique_genes,
        min_sgrna_threshold
    );

    const NumericVector observed_positive = observed["observed_positive_scores"];
    const NumericVector observed_negative = observed["observed_negative_scores"];
    const LogicalVector valid_genes_r = observed["valid_genes"];
    const int valid_gene_count = observed["n_valid_genes"];
    const int resolved_threads = resolve_thread_count(n_threads, n_permutations);

    std::vector<double> observed_positive_values(static_cast<std::size_t>(n_genes));
    std::vector<double> observed_negative_values(static_cast<std::size_t>(n_genes));
    std::vector<unsigned char> valid_genes(static_cast<std::size_t>(n_genes), 0U);
    for (int gene_index = 0; gene_index < n_genes; ++gene_index) {
        observed_positive_values[static_cast<std::size_t>(gene_index)] = observed_positive[gene_index];
        observed_negative_values[static_cast<std::size_t>(gene_index)] = observed_negative[gene_index];
        valid_genes[static_cast<std::size_t>(gene_index)] = valid_genes_r[gene_index] ? 1U : 0U;
    }

    std::vector<std::vector<std::uint64_t>> positive_counts(
        static_cast<std::size_t>(resolved_threads),
        std::vector<std::uint64_t>(static_cast<std::size_t>(n_genes), 0U)
    );
    std::vector<std::vector<std::uint64_t>> negative_counts(
        static_cast<std::size_t>(resolved_threads),
        std::vector<std::uint64_t>(static_cast<std::size_t>(n_genes), 0U)
    );

    std::atomic<bool> cancelled(false);
    std::exception_ptr worker_error;
    std::mutex error_mutex;

    const auto worker = [&](int thread_index) {
        try {
            PermutationWorkspace workspace(scores.size(), static_cast<std::size_t>(n_genes));
            std::vector<std::uint64_t>& local_positive = positive_counts[static_cast<std::size_t>(thread_index)];
            std::vector<std::uint64_t>& local_negative = negative_counts[static_cast<std::size_t>(thread_index)];
            const int begin_permutation = static_cast<int>(static_cast<std::int64_t>(n_permutations) * thread_index / resolved_threads);
            const int end_permutation = static_cast<int>(static_cast<std::int64_t>(n_permutations) * (thread_index + 1) / resolved_threads);

            for (int local_permutation = begin_permutation;
                 local_permutation < end_permutation && !cancelled.load(std::memory_order_relaxed);
                 ++local_permutation) {
                std::copy(
                    gene_indices.begin(),
                    gene_indices.end(),
                    workspace.shuffled_gene_indices.begin()
                );

                const std::uint64_t global_permutation = static_cast<std::uint64_t>(permutation_offset) +
                    static_cast<std::uint64_t>(local_permutation);
                const std::uint64_t base_seed = static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(seed)
                );
                std::mt19937_64 rng(splitmix64(base_seed ^ splitmix64(global_permutation)));
                deterministic_shuffle(workspace.shuffled_gene_indices, rng);
                group_scores(scores, workspace.shuffled_gene_indices, offsets, workspace);

                for (int gene_index = 0; gene_index < n_genes; ++gene_index) {
                    if (!valid_genes[static_cast<std::size_t>(gene_index)]) {
                        continue;
                    }

                    const std::size_t begin_index = offsets[static_cast<std::size_t>(gene_index)];
                    const std::size_t end_index = offsets[static_cast<std::size_t>(gene_index) + 1U];
                    const ScorePair permutation_scores = calculate_score_pair(
                        workspace.grouped_scores.begin() + static_cast<std::ptrdiff_t>(begin_index),
                        workspace.grouped_scores.begin() + static_cast<std::ptrdiff_t>(end_index)
                    );

                    if (permutation_scores.positive >= observed_positive_values[static_cast<std::size_t>(gene_index)]) {
                        ++local_positive[static_cast<std::size_t>(gene_index)];
                    }
                    if (permutation_scores.negative <= observed_negative_values[static_cast<std::size_t>(gene_index)]) {
                        ++local_negative[static_cast<std::size_t>(gene_index)];
                    }
                }
            }
        } catch (...) {
            cancelled.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!worker_error) {
                worker_error = std::current_exception();
            }
        }
    };

    if (show_progress) {
        Rcout << "Permutation Test: 0/" << n_permutations
              << " using " << resolved_threads << " C++ thread(s)" << std::endl;
    }
    checkUserInterrupt();

    if (resolved_threads == 1) {
        worker(0);
    } else {
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(resolved_threads));
        try {
            for (int thread_index = 0; thread_index < resolved_threads; ++thread_index) {
                workers.emplace_back(worker, thread_index);
            }
        } catch (...) {
            cancelled.store(true, std::memory_order_relaxed);
            for (std::thread& thread : workers) {
                thread.join();
            }
            throw;
        }
        for (std::thread& thread : workers) {
            thread.join();
        }
    }

    if (worker_error) {
        std::rethrow_exception(worker_error);
    }
    checkUserInterrupt();

    NumericVector p_positive(n_genes, R_NaN);
    NumericVector p_negative(n_genes, R_NaN);
    NumericVector positive_extreme_counts(n_genes, 0.0);
    NumericVector negative_extreme_counts(n_genes, 0.0);
    NumericVector valid_permutation_counts(n_genes, 0.0);

    for (int gene_index = 0; gene_index < n_genes; ++gene_index) {
        if (!valid_genes[static_cast<std::size_t>(gene_index)]) {
            continue;
        }

        std::uint64_t positive_total = 0U;
        std::uint64_t negative_total = 0U;
        for (int thread_index = 0; thread_index < resolved_threads; ++thread_index) {
            positive_total += positive_counts[static_cast<std::size_t>(thread_index)][static_cast<std::size_t>(gene_index)];
            negative_total += negative_counts[static_cast<std::size_t>(thread_index)][static_cast<std::size_t>(gene_index)];
        }

        positive_extreme_counts[gene_index] = static_cast<double>(positive_total);
        negative_extreme_counts[gene_index] = static_cast<double>(negative_total);
        valid_permutation_counts[gene_index] = static_cast<double>(n_permutations);
        p_positive[gene_index] = static_cast<double>(positive_total + 1U) /
            (static_cast<double>(n_permutations) + 1.0);
        p_negative[gene_index] = static_cast<double>(negative_total + 1U) /
            (static_cast<double>(n_permutations) + 1.0);
    }

    if (show_progress) {
        Rcout << "Permutation Test: " << n_permutations << "/" << n_permutations
              << " completed" << std::endl;
    }

    return List::create(
        Named("P_positive") = p_positive,
        Named("P_negative") = p_negative,
        Named("observed_scores") = observed_positive,
        Named("observed_positive_scores") = observed_positive,
        Named("observed_negative_scores") = observed_negative,
        Named("positive_extreme_counts") = positive_extreme_counts,
        Named("negative_extreme_counts") = negative_extreme_counts,
        Named("valid_permutation_counts") = valid_permutation_counts,
        Named("n_valid_genes") = valid_gene_count,
        Named("n_permutations") = n_permutations,
        Named("n_threads") = resolved_threads,
        Named("seed") = seed,
        Named("permutation_offset") = permutation_offset,
        Named("gene_names") = unique_genes
    );
}
