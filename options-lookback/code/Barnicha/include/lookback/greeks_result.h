#pragma once
/**
 * @file greeks_result.h
 * @brief Definition of a container for option price and Greeks with uncertainty measures.
 */
#include "estimate_stats.h"


namespace Barnicha {

/**
 * @brief Price and Greeks container with uncertainty information.
 *
 * - Monte Carlo: each field has standard error + CI95.
 * - Analytic: std_dev=0, CI is degenerate at the exact value.
 */
struct GreeksResult 
{
    EstimateStats price;
    EstimateStats delta;
    EstimateStats gamma;
    EstimateStats theta;
    EstimateStats rho;
    EstimateStats vega;
};

} // namespace lookback
