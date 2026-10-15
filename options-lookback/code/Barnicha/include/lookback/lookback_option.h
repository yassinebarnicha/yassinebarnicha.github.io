#pragma once
/**
 * @file lookback_option.h
 * @brief Definition of the floating-strike lookback option contract.
 */

namespace Barnicha 
{

/**
 * @brief Lookback option type.
 */
enum class LookbackType 
{
    Call, ///< floating-strike call: payoff = s_T - min_{t<=T} s_t
    Put   ///< floating-strike put : payoff = max_{t<=T} s_t - s_T
};

/**
 * @brief European floating-strike lookback option (continuous monitoring).
 *
 * We assume the option is issued at t=0 and continuously monitored on [0,T].
 */
class LookbackOption 
{
    
public:
    /**
     * @param type call or put
     * @param maturity time-to-maturity (years)
     */
    LookbackOption(LookbackType type, double maturity);

    /**
     * @brief Payoff given terminal value and path extrema.
     *
     * @param st terminal value s_T
     * @param min_s running minimum of s_t on [0,T]
     * @param max_s running maximum of s_t on [0,T]
     */
    double payoff(double st, double min_s, double max_s) const;

    /** @brief maturity T (years). */
    double maturity() const;

    /** @brief option type (call/put). */
    LookbackType type() const;

private:
    LookbackType type_;
    double maturity_;
};

}
