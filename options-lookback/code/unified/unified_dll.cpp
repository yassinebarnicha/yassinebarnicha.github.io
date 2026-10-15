/**
 * @file lookback_dll.cpp
 * @brief Unified DLL entry point for pricing lookback options.
 */

#include "lookback/lookback_dll.h"

#include "LookbackDLL.hpp"

#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport) ///< Export symbol for Windows DLL builds.
#else
#define DLL_EXPORT ///< Empty export macro for non-Windows platforms.
#endif

/**
 * @brief Exported DLL function to price a lookback option.
 *
 * This function serves as the main external interface of the lookback pricing DLL.
 * It selects the appropriate pricing method based on the input parameter `method`
 * and delegates computation to the corresponding implementation function.
 *
 * Supported pricing methods:
 * - method = 0 → Monte Carlo pricing with optional BGK correction
 * - method = 1 → Conditional expectation pricing method
 *
 * The computed results (price and optionally standard errors or Greeks) are
 * written into the provided output buffer.
 *
 * @param method Pricing method selector.
 *        - 0: Monte Carlo with optional BGK correction
 *        - 1: Conditional expectation method
 *
 * @param ttm Time to maturity (in years).
 *
 * @param optionType Option type identifier.
 *        Typically:
 *        - 0: Call lookback option
 *        - 1: Put lookback option
 *
 * @param spot Current underlying asset price.
 *
 * @param rate Risk-free interest rate (annualized).
 *
 * @param vol Volatility of the underlying asset (annualized).
 *
 * @param paths Number of Monte Carlo simulation paths.
 *
 * @param steps Number of time steps per path (used in Monte Carlo method).
 *
 * @param seed Random number generator seed.
 *
 * @param withStdErrors Flag indicating whether standard errors should be computed.
 *        - 0: No standard errors
 *        - 1: Compute standard errors
 *
 * @param useBGK Flag indicating whether BGK correction should be applied.
 *        - 0: Disabled
 *        - 1: Enabled
 *
 * @param out_values Pointer to the output buffer where results will be written.
 *        The buffer must be allocated by the caller.
 *
 * @param out_size Pointer to an integer that will receive the number of values written.
 *
 * @return int Status code:
 *         - 0: Success
 *         - negative value: Error code
 *         - -999: Invalid method specified
 */
extern "C" DLL_EXPORT
int __stdcall price_lookback(
    int method,
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
    if(method == 0)
        return price_lookback_bgk( ttm,optionType,spot,rate,vol, paths, steps, seed, withStdErrors, useBGK,out_values, out_size);

    if(method == 1)
        return price_lookback_cond(
            ttm, optionType, spot, rate, vol,
            paths, seed,
            out_values, out_size
        );

    return -999;
}
