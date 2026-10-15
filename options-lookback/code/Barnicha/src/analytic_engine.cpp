#include "lookback/analytic_engine.h"
#include "lookback/analytic_lookback.h"

namespace Barnicha 
{

std::vector<double> AnalyticEngine::discounted_payoff_samples(const MarketData& market,
                                                              const LookbackOption& option) const 
{
    const double p = AnalyticLookback::price(market, option);
    return std::vector<double>{p};
}

}
