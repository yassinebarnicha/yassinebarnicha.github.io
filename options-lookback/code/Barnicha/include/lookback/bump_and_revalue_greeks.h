#pragma once
/**
 * @file bump_and_revalue_greeks.h
 * @brief Declaration of the bump-and-revalue Greeks calculator.
 */
#include "greeks_calculator.h"
#include <optional>

namespace Barnicha {

/**
 * @brief Bump-and-revalue Greeks with uncertainty reporting.
 *
 * - With MonteCarloEngine: computes Greeks per path (CRN) => standard errors + CI95.
 * - With deterministic engines (AnalyticEngine): std_dev=0.
 */
class BumpAndRevalueGreeks final : public GreeksCalculator {
public:
    BumpAndRevalueGreeks(bool compute_price = true,
                         bool compute_delta = true,
                         bool compute_gamma = true,
                         bool compute_theta = true,
                         bool compute_rho   = true,
                         bool compute_vega  = true);

    BumpAndRevalueGreeks& set_bumps(std::optional<double> d_spot_rel = std::nullopt,
                                    std::optional<double> d_rate_abs = std::nullopt,
                                    std::optional<double> d_vol_abs = std::nullopt,
                                    std::optional<double>  d_t_abs = std::nullopt);

    GreeksResult compute(const PricingEngine& engine,
                         const MarketData& market,
                         const LookbackOption& option) const override;

private:
    bool compute_price_;
    bool compute_delta_;
    bool compute_gamma_;
    bool compute_theta_;
    bool compute_rho_;
    bool compute_vega_;

    double d_spot_rel_ = 1e-4;
    double d_rate_abs_ = 1e-4;
    double d_vol_abs_  = 1e-4;
    double d_t_abs_    = 1e-4;
};

} // namespace lookback
