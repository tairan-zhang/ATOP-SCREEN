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


#pragma once
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace atop_ncs_native {
using namespace Rcpp;
// R's real-valued mean uses extended accumulation and a correction pass.
// Keep this policy separate from the established permutation arithmetic.
// Semantics: R Core, src/main/summary.c (real_mean), GPL-2.0-or-later.
template <typename Iterator>
double scoring_mean(Iterator begin, Iterator end) {
    const auto n = std::distance(begin, end);
    if (n == 0) return R_NaN;
    long double value = 0.0;
    for (auto it = begin; it != end; ++it) value += *it;
    const bool finite_sum = std::isfinite(static_cast<double>(value));
    if (finite_sum) {
        value /= n;
    } else {
        value = 0.0;
        for (auto it = begin; it != end; ++it) value += *it / static_cast<double>(n);
    }
    if (std::isfinite(static_cast<double>(value))) {
        long double correction = 0.0;
        for (auto it = begin; it != end; ++it) {
            const long double residual = *it - value;
            correction += finite_sum ? residual : residual / n;
        }
        value += finite_sum ? correction / n : correction;
    }
    return static_cast<double>(value);
}

double scoring_median(std::vector<double>& values) {
    std::sort(values.begin(), values.end());
    const auto middle = values.size() / 2U;
    if (values.size() % 2U) return values[middle];
    return scoring_mean(values.begin() + middle - 1U, values.begin() + middle + 1U);
}

inline NumericMatrix calculate(NumericMatrix lfc, double delta0, double delta1, double delta2) {
    if (!std::isfinite(delta0) || delta0 < 0 || !std::isfinite(delta1) || delta1 <= 0 ||
        !std::isfinite(delta2) || delta2 <= 0) stop("Invalid NCS parameters");
    NumericMatrix result(4, lfc.nrow());
    rownames(result) = CharacterVector::create("dir_score", "mag_score", "stab_score", "NCS");
    std::vector<double> values, magnitude, deviations;
    values.reserve(lfc.ncol()); magnitude.reserve(lfc.ncol()); deviations.reserve(lfc.ncol());
    for (int row = 0; row < lfc.nrow(); ++row) {
        if (row % 16384 == 0) checkUserInterrupt();
        values.clear(); magnitude.clear(); deviations.clear();
        int positive = 0, negative = 0, neutral = 0;
        for (int col = 0; col < lfc.ncol(); ++col) {
            const double value = lfc(row, col);
            if (!std::isfinite(value)) continue;
            values.push_back(value);
            magnitude.push_back(std::abs(value));
            if (value > delta0) ++positive;
            else if (value < -delta0) ++negative;
            else ++neutral;
        }
        if (values.empty()) {
            for (int component = 0; component < 4; ++component) result(component, row) = NA_REAL;
            continue;
        }
        const double center = scoring_median(values);
        for (double value : values) deviations.push_back(std::abs(value - center));
        const double direction = static_cast<double>(std::max({positive, negative, neutral})) / values.size();
        const double mag = std::min(1.0, scoring_median(magnitude) / delta1);
        const double stability = 1.0 / (1.0 + scoring_median(deviations) / delta2);
        result(0, row) = direction;
        result(1, row) = mag;
        result(2, row) = stability;
        result(3, row) = direction * mag * stability;
    }
    return result;
}

}
