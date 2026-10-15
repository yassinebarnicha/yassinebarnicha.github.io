#include "lookback/pricing_engine.h"
#include <cmath>
#include <numeric>

namespace Barnicha 
{

EstimateStats 
PricingEngine::stats_from_samples(const std::vector<double>& xs) 
{
    const std::size_t n = xs.size();
    if (n == 0) return EstimateStats{};

    const double mean = std::accumulate(xs.begin(), xs.end(), 0.0) / static_cast<double>(n);

    double var = 0.0;
    for (double x : xs) 
    {
        const double d = x - mean;
        var += d * d;
    }
    const double sample_var = (n > 1) ? (var / static_cast<double>(n - 1)) : 0.0;
    const double sample_std = std::sqrt(sample_var);

    const double std_err = sample_std / std::sqrt(static_cast<double>(n));
    const double z = 1.96;

    return EstimateStats{mean, std_err, mean - z * std_err, mean + z * std_err};
}

}
