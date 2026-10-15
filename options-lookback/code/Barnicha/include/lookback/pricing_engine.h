#pragma once
/**
 * @file pricing_engine.h
 * @brief Declaration of the abstract pricing engine interface.
 */
#include "market_data.h"
#include "lookback_option.h"
#include "estimate_stats.h"
#include <vector>

namespace Barnicha 
{

/**
 * @brief Pricing engine interface (Strategy pattern).
 *
 * Each engine exposes the *samples* of its estimator, via discounted_payoffs():
 *  - Monte Carlo: one discounted payoff per path
 *  - Analytic   : a single sample {exact price}
 *
 * This lets Greeks and reporting compute standard errors uniformly, without
 * branching on the concrete engine type.
 */
class PricingEngine 
{
public:
    virtual ~PricingEngine() = default;

    /**
     * @brief Discounted samples defining the estimator.
     */
    virtual std::vector<double> discounted_payoff_samples(const MarketData& market,
                                                          const LookbackOption& option) const = 0;

    /**
     * @brief Computes stats (e.g mean, stdev, IC) from discounted payoffd.
     */
    static EstimateStats stats_from_samples(const std::vector<double>& xs); 
};

}