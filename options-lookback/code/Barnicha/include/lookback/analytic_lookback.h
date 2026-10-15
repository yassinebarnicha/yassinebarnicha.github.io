#pragma once
/**
 * @file analytic_lookback.h
 * @brief Declaration of the closed-form analytic pricer for floating-strike lookback options.
 */
#include "market_data.h"
#include "lookback_option.h"

namespace Barnicha {

/**
 * @brief Closed-form price for floating-strike lookback options (continuous monitoring).
 *
 * Includes a dedicated `rate == 0` limit to avoid numerical issues when the standard
 * formula contains terms proportional to 1/r.
 */
class AnalyticLookback {
public:
    /**
     * @brief Analytic lookback price under Black–Scholes.
     *
     * Dispatches to:
     * - `price_r_nonzero` if |r| is not tiny,
     * - `price_r_zero_limit` otherwise.
     */
    static double price(const MarketData& market, const LookbackOption& option);

private:
    static double price_r_nonzero(const MarketData& market, const LookbackOption& option);
    static double price_r_zero_limit(const MarketData& market, const LookbackOption& option);
};

} // namespace lookback
