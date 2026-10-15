/**
 * @file LookbackDLL.cpp
 * @brief Implementation of the price_lookback_bgk function.
 */

#include "LookbackDLL.hpp"
#include "AnalyticLookback.hpp"
#include "Pricer.hpp"

using namespace Chivet ;

int price_lookback_bgk(
    double ttm,
    int optionType,
    double spot,
    double rate,
    double vol,
    long paths,
    long steps,
    long seed,
    int withStdErrors,
    int useBGK,
    double* out_values,
    int* out_size)
{
    if (!out_values || !out_size) return -1;
    if (ttm <= 0.0) return -2;
    if (paths <= 0) return -3;

    try {

        OptionType type =
            (optionType == 1) ? OptionType::Put : OptionType::Call;

        // =========================
        // MONTE CARLO PART (BGK)
        // =========================
        PricerLookbackOption pricer;
        pricer.setOption(type, ttm);
        pricer.setMarketParameters(spot, rate, vol);
        pricer.setMonteCarloSimulator(paths, steps, seed);
        pricer.setBgkCorrection(useBGK != 0);

        PricingResult mc =
            pricer.compute(withStdErrors != 0);

        // =========================
        // ANALYTIC PART
        // =========================
        LookbackAnalyticSolution analytic;

        double an_price =
            analytic.price(type, spot, ttm, vol, rate);

        double an_delta =
            analytic.delta(type, spot, ttm, vol, rate);

        double an_gamma =
            analytic.gamma(type, spot, ttm, vol, rate);

        double an_theta =
            analytic.theta(type, spot, ttm, vol, rate);

        double an_rho =
            analytic.rho(type, spot, ttm, vol, rate);

        double an_vega =
            analytic.vega(type, spot, ttm, vol, rate);

        int i = 0;

        auto fill_mc = [&](double est, double std)
        {   double z = 1.96;
            out_values[i++] = est;
            out_values[i++] = std;
            out_values[i++] = est - z * std;
            out_values[i++] = est + z * std;
        };

        auto fill_an = [&](double est)
        {
            out_values[i++] = est;
            out_values[i++] = 0.0;
            out_values[i++] = est;
            out_values[i++] = est;
        };

        // =========================
        // MONTE CARLO BLOCK (24)
        // =========================
        fill_mc(mc.price, mc.priceStd);
        fill_mc(mc.delta, mc.deltaStd);
        fill_mc(mc.gamma, mc.gammaStd);
        fill_mc(mc.theta, mc.thetaStd);
        fill_mc(mc.rho,   mc.rhoStd);
        fill_mc(mc.vega,  mc.vegaStd);

        // =========================
        // ANALYTIC BLOCK (24)
        // =========================
        fill_an(an_price);
        fill_an(an_delta);
        fill_an(an_gamma);
        fill_an(an_theta);
        fill_an(an_rho);
        fill_an(an_vega);

        *out_size = i;   // 48

        return 0;
    }
    catch (...) {
        return -100;
    }
}
