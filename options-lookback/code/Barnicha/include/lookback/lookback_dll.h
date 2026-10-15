 #pragma once
/**
 * @file lookback_cond.h
 * @brief Declaration of the conditional Monte Carlo lookback option pricer.
 */



 

/**
 * @brief Prices a lookback option using the conditional expectation method.
 *
 * Uses Monte Carlo simulation of the terminal asset price and applies the
 * conditional distribution of the maximum/minimum to compute the payoff.
 * Results are written into the output buffer.
 *
 * @param ttm Time to maturity (years).
 * @param optionType Option type (0 = Call, 1 = Put).
 * @param spot Current spot price.
 * @param rate Risk-free rate.
 * @param vol Volatility.
 * @param paths Number of Monte Carlo paths.
 * @param seed RNG seed.
 * @param out_values Output buffer for results.
 * @param out_size Number of values written.
 *
 * @return 0 on success, non-zero on error.
 */
int price_lookback_cond(
    double ttm,
    int optionType,   // 0=Call, 1=Put
    double spot,
    double rate,
    double vol,
    long paths,
    long seed,
    double* out_values, // array of doubles (size 48)
    int* out_size        // returns how many doubles were written
);
