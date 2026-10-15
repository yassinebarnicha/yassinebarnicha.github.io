#include "lookback/bump_and_revalue_greeks.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace Barnicha 
{


BumpAndRevalueGreeks::BumpAndRevalueGreeks(bool compute_price,
                                           bool compute_delta,
                                           bool compute_gamma,
                                           bool compute_theta,
                                           bool compute_rho,
                                           bool compute_vega)
: compute_price_(compute_price),
  compute_delta_(compute_delta),
  compute_gamma_(compute_gamma),
  compute_theta_(compute_theta),
  compute_rho_(compute_rho),
  compute_vega_(compute_vega) 
  {}

BumpAndRevalueGreeks& 
BumpAndRevalueGreeks::set_bumps(std::optional<double> d_spot_rel,
                                std::optional<double> d_rate_abs,
                                std::optional<double> d_vol_abs,
                                std::optional<double> d_t_abs) 
{
    if(d_spot_rel)
        d_spot_rel_ = *d_spot_rel;
    if(d_rate_abs)
        d_rate_abs_ = *d_rate_abs;
    if(d_vol_abs)
        d_vol_abs_  = *d_vol_abs;
    if(d_t_abs)
        d_t_abs_    = *d_t_abs;
    return *this;
}

GreeksResult BumpAndRevalueGreeks::compute(const PricingEngine& engine,
                                           const MarketData& market,
                                           const LookbackOption& option) const 
{
    GreeksResult result;

    // Base discounted samples (size = n_paths for MC, size = 1 for analytic).
    const auto base = engine.discounted_payoff_samples(market, option);
    const std::size_t n = base.size();

    if (compute_price_) result.price = PricingEngine::stats_from_samples(base);

    if (compute_delta_ || compute_gamma_) 
    {
        const double h_s = std::max(1e-12, market.spot * d_spot_rel_);
        MarketData up_s = market; up_s.spot += h_s;
        MarketData dn_s = market; dn_s.spot -= h_s;

        const auto up = engine.discounted_payoff_samples(up_s, option);
        const auto dn = engine.discounted_payoff_samples(dn_s, option);
        const std::size_t m = std::min({n, up.size(), dn.size()});

        std::vector<double> delta_samples;
        std::vector<double> gamma_samples;
        delta_samples.reserve(m);
        gamma_samples.reserve(m);

        for (std::size_t i = 0; i < m; ++i) 
        {
            if (compute_delta_) delta_samples.push_back((up[i] - dn[i]) / (2.0 * h_s));
            if (compute_gamma_) gamma_samples.push_back((up[i] - 2.0 * base[i] + dn[i]) / (h_s * h_s));
        }

        if (compute_delta_) result.delta = PricingEngine::stats_from_samples(delta_samples);
        if (compute_gamma_) result.gamma = PricingEngine::stats_from_samples(gamma_samples);
    }

    if (compute_rho_) 
    {
        const double h_r = std::max(1e-12, d_rate_abs_);
        MarketData up_r = market; up_r.rate += h_r;
        MarketData dn_r = market; dn_r.rate -= h_r;

        const auto up = engine.discounted_payoff_samples(up_r, option);
        const auto dn = engine.discounted_payoff_samples(dn_r, option);
        const std::size_t m = std::min({n, up.size(), dn.size()});

        std::vector<double> rho_samples;
        rho_samples.reserve(m);
        for (std::size_t i = 0; i < m; ++i)
        {
            rho_samples.push_back((up[i] - dn[i]) / (2.0 * h_r));
        }

        result.rho = PricingEngine::stats_from_samples(rho_samples);
    }

    if (compute_vega_) 
    {
        const double h_v = std::max(1e-12, d_vol_abs_);
        MarketData up_v = market; up_v.volatility += h_v;
        MarketData dn_v = market; dn_v.volatility -= h_v;

        const auto up = engine.discounted_payoff_samples(up_v, option);
        const auto dn = engine.discounted_payoff_samples(dn_v, option);
        const std::size_t m = std::min({n, up.size(), dn.size()});

        std::vector<double> vega_samples;
        vega_samples.reserve(m);
        for (std::size_t i = 0; i < m; ++i) 
        {
            vega_samples.push_back((up[i] - dn[i]) / (2.0 * h_v));
        }

        result.vega = PricingEngine::stats_from_samples(vega_samples);
    }

    if (compute_theta_) 
    {
        const double h_t = std::max(1e-12, d_t_abs_);
        const double maturity_dn = std::max(1e-10, option.maturity() - h_t);
        LookbackOption option_dn(option.type(), maturity_dn);

        const auto dn = engine.discounted_payoff_samples(market, option_dn);
        const std::size_t m = std::min(n, dn.size());

        std::vector<double> theta_samples;
        theta_samples.reserve(m);
        for (std::size_t i = 0; i < m; ++i)
        {
            theta_samples.push_back((dn[i] - base[i]) / h_t);
        }

        result.theta = PricingEngine::stats_from_samples(theta_samples);
    }

    return result;
}

}
