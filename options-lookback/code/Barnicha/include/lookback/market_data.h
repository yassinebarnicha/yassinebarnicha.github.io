#pragma once
/**
 * @file market_data.h
 * @brief Definition of the Black–Scholes market parameters.
 */
namespace Barnicha {

/**
 * @brief Black–Scholes market inputs.
 *
 * These parameters define the risk-neutral dynamics:
 *  dS_t = r S_t dt + sigma S_t dW_t.
 */
struct MarketData 
{
    double spot;       ///< spot price s0
    double rate;       ///< risk-free rate r (continuous compounding)
    double volatility; ///< volatility sigma
};

}
