#include "lookback/lookback_payoff_evaluator.h"
#include <cmath>

namespace Barnicha 
{

LookbackPayoffEvaluator::LookbackPayoffEvaluator(const MarketData& market, const LookbackOption& option)
: market_(market), option_(option) 
{}

double LookbackPayoffEvaluator::operator()(const OneStepSample& sample) const 
{
    const double s0 = market_.spot;
    const double rate = market_.rate;
    const double sigma = market_.volatility;
    const double maturity = option_.maturity();

    // Work in log space: X_t = log S_t.
    const double x0 = std::log(s0);

    // Risk-neutral drift for X_t.
    const double drift = rate - 0.5 * sigma * sigma;

    // (1) simulate the terminal log-price exactly:
    //     X_T = x0 + drift*T + sigma*sqrt(T)*Z
    const double x_t = x0 + drift * maturity + sigma * std::sqrt(maturity) * sample.z;

    // Helper: A = -0.5 sigma^2 T log(U) >= 0.
    // It appears after inverting the bridge extremum CDF.
    auto a_from_u = [&](double u) 
    {
        return -0.5 * sigma * sigma * maturity * std::log(u);
    };

    // (2) maximum of the bridge (upper quadratic root)
    const double a_max = a_from_u(sample.u_max);
    const double disc_max = std::sqrt((x0 - x_t) * (x0 - x_t) + 4.0 * a_max);
    const double x_max = 0.5 * (x0 + x_t + disc_max);

    // (3) minimum of the bridge (lower quadratic root)
    const double a_min = a_from_u(sample.u_min);
    const double disc_min = std::sqrt((x0 - x_t) * (x0 - x_t) + 4.0 * a_min);
    const double x_min = 0.5 * (x0 + x_t - disc_min);

    // Back to spot space.
    const double st = std::exp(x_t);
    const double max_s = std::exp(x_max);
    const double min_s = std::exp(x_min);

    return option_.payoff(st, min_s, max_s);
}

}
