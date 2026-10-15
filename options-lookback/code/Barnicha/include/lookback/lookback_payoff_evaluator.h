#pragma once
/**
 * @file lookback_payoff_evaluator.h
 * @brief Declaration of the payoff evaluator using exact one-step Brownian bridge sampling.

 */
#include "market_data.h"
#include "lookback_option.h"
#include "one_step_bridge_sampler.h"


namespace Barnicha 
{

/**
 * @brief Computes one payoff using exact one-step sampling (no time grid).
 *
 * Let X_t = log S_t. Under Black–Scholes:
 *   X_T = x0 + (r - 0.5 sigma^2) T + sigma sqrt(T) z.
 *
 * Conditional on (X_0=x0, X_T=x_t), the (driftless) Brownian bridge maximum has CDF:
 *   P(max X <= m | x0, x_t) = 1 - exp(-2 (m-x0)(m-x_t)/(sigma^2 T)), for m >= max(x0, x_t).
 *
 * Inversion yields a quadratic equation in m with two roots:
 * - the upper root gives the maximum,
 * - the lower root gives the minimum (by symmetry).
 */
class LookbackPayoffEvaluator 
{
public:
    /**
     * @param market Black–Scholes inputs
     * @param option lookback contract
     */
    LookbackPayoffEvaluator(const MarketData& market, const LookbackOption& option);

    /**
     * @brief Evaluate the payoff for one random sample.
     */
    double operator()(const OneStepSample& sample) const;

private:
    MarketData market_;
    LookbackOption option_;
};

}
