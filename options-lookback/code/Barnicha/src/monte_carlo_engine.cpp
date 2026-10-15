#include "lookback/monte_carlo_engine.h"
#include "lookback/lookback_payoff_evaluator.h"
#include <algorithm>
#include <cmath>

namespace Barnicha 
{

MonteCarloEngine::MonteCarloEngine(const OneStepBridgeSampler& sampler,
                                   std::size_t n_paths,
                                   bool antithetic_on_normal)
: sampler_(sampler),
  n_paths_(n_paths),
  antithetic_on_normal_(antithetic_on_normal) 
  {}

std::vector<double> MonteCarloEngine::discounted_payoff_samples(const MarketData& market,
                                                                const LookbackOption& option) const 
{
    const auto& samples = sampler_.samples();
    const std::size_t n_used = std::min(n_paths_, samples.size());

    LookbackPayoffEvaluator evaluator(market, option);

    std::vector<double> discounted;
    discounted.reserve(n_used);

    const double discount = std::exp(-market.rate * option.maturity());

    for (std::size_t i = 0; i < n_used; ++i) 
    {
        const double payoff_1 = evaluator(samples[i]);
        double payoff = payoff_1;

        if (antithetic_on_normal_) 
        {
            OneStepSample anti = samples[i];
            anti.z = -anti.z;
            const double payoff_2 = evaluator(anti);
            payoff = 0.5 * (payoff_1 + payoff_2);
        }

        discounted.push_back(discount * payoff);
    }

    return discounted;
}

}
