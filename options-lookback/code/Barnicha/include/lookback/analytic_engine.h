#pragma once
/**
 * @file analytic_engine.h
 * @brief Declaration of the analytic pricing engine for lookback options.
 */

#include "pricing_engine.h"
#include <vector>


namespace Barnicha {

/**
 * @brief Analytic pricing engine (deterministic).
 *
 * discounted_payoffs() returns a single sample {exact price},
 * hence std_dev = 0 and CI is degenerate.
 */
class AnalyticEngine : public PricingEngine {
public:
    std::vector<double> discounted_payoff_samples(const MarketData& market,
                                                  const LookbackOption& option) const override;

};

} // namespace lookback
