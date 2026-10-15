/**
 * @file Lookback_dll.cpp
 * @brief Implementation of the price_lookback_cond function.
 */


#include <cstddef>

#include "lookback/one_step_bridge_sampler.h"
#include "lookback/monte_carlo_engine.h"
#include "lookback/analytic_engine.h"
#include "lookback/bump_and_revalue_greeks.h"

using namespace Barnicha;

int  price_lookback_cond(
    double ttm,
    int optionType,   // 0=Call, 1=Put
    double spot,
    double rate,
    double vol,
    long paths,
    long seed,
    double* out_values,
    int* out_size)
{
    if (!out_values || !out_size) return -1;
    if (ttm <= 0.0) return -2;
    if (paths <= 0) return -3;

    const LookbackType optType =
        (optionType == 1) ? LookbackType::Put : LookbackType::Call;

    const MarketData market{spot, rate, vol};
    const LookbackOption option(optType, ttm);

    OneStepBridgeSampler sampler((unsigned)seed);
    sampler.generate((std::size_t)paths);

    MonteCarloEngine mc_engine(sampler, (std::size_t)paths, true);
    AnalyticEngine an_engine;

    BumpAndRevalueGreeks greeks(true, true, true, true, true, true);

    const GreeksResult mc = greeks.compute(mc_engine, market, option);
    const GreeksResult an = greeks.compute(an_engine, market, option);

    auto fill = [&](int& idx, const EstimateStats& s) {
        out_values[idx++] = s.estimate;
        out_values[idx++] = s.std_dev;
        out_values[idx++] = s.ci95_low;
        out_values[idx++] = s.ci95_high;
    };

    int i = 0;

    // MonteCarlo (6 metrics × 4 values = 24)
    fill(i, mc.price);
    fill(i, mc.delta);
    fill(i, mc.gamma);
    fill(i, mc.theta);
    fill(i, mc.rho);
    fill(i, mc.vega);

    // Analytic (another 24)
    fill(i, an.price);
    fill(i, an.delta);
    fill(i, an.gamma);
    fill(i, an.theta);
    fill(i, an.rho);
    fill(i, an.vega);

    *out_size = i; // should be 48

    return 0;
}
