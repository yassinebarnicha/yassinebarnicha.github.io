#pragma once
/**
 * @file estimate_stats.h
 * @brief Definition of a structure representing statistical results of an estimate.
 */
#include <cmath>

namespace Barnicha {

/**
 * @brief Statistics for an estimate (typically Monte Carlo).
 *
 * We report:
 *  - estimate  : point estimate (sample mean)
 *  - std_dev   : standard error of the estimate (sample_std / sqrt(n))
 *  - ci95_low/high : symmetric 95% confidence interval (normal approximation)
 */
struct EstimateStats {
    double estimate = 0.0;
    double std_dev  = 0.0;
    double ci95_low = 0.0;
    double ci95_high= 0.0;

    static EstimateStats deterministic(double x) {
        return EstimateStats{x, 0.0, x, x};
    }
};

} // namespace lookback
