#pragma once
/**
 * @file monte_carlo_engine.h
 * @brief Declaration of the one-step Monte Carlo pricing engine.
 */
#include "pricing_engine.h"
#include "one_step_bridge_sampler.h"
#include <cstddef>
#include <vector>

namespace Barnicha 
{

/**
 * @brief One-step Monte Carlo engine reusing stored samples (CRN).
 *
 * discounted_payoffs() returns one discounted payoff per path (optionally antithetic-averaged on z).
 * price() uses the default PricingEngine implementation to compute EstimateStats.
 */
class MonteCarloEngine : public PricingEngine 
{
public:
    MonteCarloEngine(const OneStepBridgeSampler& sampler,
                     std::size_t n_paths,
                     bool antithetic_on_normal = true);

    std::vector<double> discounted_payoff_samples(const MarketData& market,
                                                  const LookbackOption& option) const override;

private:
    const OneStepBridgeSampler& sampler_;
    std::size_t n_paths_;
    bool antithetic_on_normal_;
};

}
