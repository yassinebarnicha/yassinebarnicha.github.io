#include "lookback/analytic_lookback.h"
#include "lookback/math.h"
#include <cmath>

namespace Barnicha 
{

double AnalyticLookback::price(const MarketData& market, const LookbackOption& option) 
{
    // Switch to a numerically stable r->0 limit when needed.
    if (std::abs(market.rate) < 1e-10) 
    {
        return price_r_zero_limit(market, option);
    }
    return price_r_nonzero(market, option);
}

double AnalyticLookback::price_r_nonzero(const MarketData& market, const LookbackOption& option) {
    // Standard closed-form for floating-strike lookbacks (continuous monitoring),
    // with running min/max starting at s0 at t=0.
    const double s0 = market.spot;
    const double rate = market.rate;
    const double sigma = market.volatility;
    const double maturity = option.maturity();

    const double sqrt_maturity = std::sqrt(maturity);

    // With start extremum = s0 we have log(s0/s0)=0. We keep the usual notation a1,a2,a3.
    const double a1 = ((rate + 0.5 * sigma * sigma) * maturity) / (sigma * sqrt_maturity);
    const double a2 = a1 - sigma * sqrt_maturity;
    const double a3 = a1 - (2.0 * rate * sqrt_maturity / sigma);

    const double discount = std::exp(-rate * maturity);
    const double mult = s0 * sigma * sigma / (2.0 * rate);

    if (option.type() == LookbackType::Call) 
    {
        const double term_1 = s0 * Math::normal_cdf(a1);
        const double term_2 = s0 * discount * Math::normal_cdf(a2);
        const double term_3 = Math::normal_cdf(-a1) - discount * Math::normal_cdf(-a3);
        return term_1 - term_2 - mult * term_3;
    }

    // Put
    const double term_1 = -s0 * Math::normal_cdf(-a1);
    const double term_2 =  s0 * discount * Math::normal_cdf(-a2);
    const double term_3 = Math::normal_cdf(a1) - discount * Math::normal_cdf(a3);
    return term_1 + term_2 + mult * term_3;
}

double AnalyticLookback::price_r_zero_limit(const MarketData& market, const LookbackOption& option) {
    // r -> 0 limit for floating-strike lookbacks (continuous monitoring).
    const double s0 = market.spot;
    const double sigma = market.volatility;
    const double maturity = option.maturity();

    const double d = 0.5 * sigma * std::sqrt(maturity);
    const double pdf = Math::normal_pdf(d);
    const double cdf = Math::normal_cdf(d);

    if (option.type() == LookbackType::Call) 
    {
        return s0 * ((1.0 + d * d) * (2.0 * cdf - 1.0) + 2.0 * d * pdf - d * d);
    }
    return s0 * ((1.0 + d * d) * (2.0 * cdf - 1.0) + 2.0 * d * pdf + d * d);
}

}
