#pragma once

/**
 * @file lookback_bgk.h
 * @brief Declaration of the Monte Carlo lookback option pricer with optional BGK correction.
 */

/**
 * @brief Prices a lookback option using Monte Carlo simulation.
 *
 * Optionally applies Broadie-Glasserman-Kou (BGK) correction to reduce
 * discretization bias. Results are written into the output buffer.
 *
 * @param ttm Time to maturity (years).
 * @param optionType Option type (call or put).
 * @param spot Current spot price.
 * @param rate Risk-free rate.
 * @param vol Volatility.
 * @param paths Number of Monte Carlo paths.
 * @param steps Number of time steps per path.
 * @param seed RNG seed.
 * @param withStdErrors Flag to compute standard errors.
 * @param useBGK Flag to enable BGK correction.
 * @param out_values Output buffer for results.
 * @param out_size Number of values written.
 *
 * @return 0 on success, non-zero on error.
 */
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
    int* out_size);
