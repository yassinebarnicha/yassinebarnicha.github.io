/**
 * @file AnalyticLookback.hpp
 * @brief Closed-form pricing and Greeks for floating-strike lookback options.
 * Most of the code is based on the formulas from
 * https://personal.ntu.edu.sg/nprivault/MA5182/lookback-options.pdf
 */

#ifndef ANALYTIC_LOOKBACK_HPP
#define ANALYTIC_LOOKBACK_HPP
#include "OptionParameters.hpp"


namespace Chivet 
{

/**
 * @class LookbackAnalyticSolution
 * @brief Provides closed-form (analytic) prices and Greeks for
 *        floating-strike lookback options under the Black-Scholes framework.
 *
 * Greeks are computed via central finite differences on the analytic price
 * formula.
 */
class LookbackAnalyticSolution {
public:
  /**
   * @brief Compute the analytic price of a floating-strike lookback option.
   * @param option Call or Put.
   * @param S0     Current spot price.
   * @param T      Time to maturity in years.
   * @param sigma  Volatility.
   * @param r      Risk-free interest rate.
   * @return The option price.
   */
  static double price(OptionType option, double S0, double T, double sigma,
                      double r);

  /**
   * @brief Compute Delta (\f$\partial V / \partial S\f$) via finite difference.
   * @copydetails price()
   */
  static double delta(OptionType option, double S0, double T, double sigma,
                      double r);

  /**
   * @brief Compute Gamma (\f$\partial^2 V / \partial S^2\f$) via finite
   * difference.
   * @copydetails price()
   */
  static double gamma(OptionType option, double S0, double T, double sigma,
                      double r);

  /**
   * @brief Compute Vega (\f$\partial V / \partial \sigma\f$) via finite
   * difference.
   * @copydetails price()
   */
  static double vega(OptionType option, double S0, double T, double sigma,
                     double r);

  /**
   * @brief Compute Theta (\f$\partial V / \partial T\f$) via finite difference.
   * @copydetails price()
   */
  static double theta(OptionType option, double S0, double T, double sigma,
                      double r);

  /**
   * @brief Compute Rho (\f$\partial V / \partial r\f$) via finite difference.
   * @copydetails price()
   */
  static double rho(OptionType option, double S0, double T, double sigma,
                    double r);

private:
  /**
   * @brief Compute \f$d_+ = \frac{\ln s + (r +
   * \frac{1}{2}\sigma^2)\tau}{\sigma\sqrt{\tau}}\f$.
   * @param tau   Time to maturity.
   * @param s     Spot ratio.
   * @param r     Risk-free rate.
   * @param sigma Volatility.
   */
  static double delta_plus(double tau, double s, double r, double sigma);

  /**
   * @brief Compute \f$d_- = \frac{\ln s + (r -
   * \frac{1}{2}\sigma^2)\tau}{\sigma\sqrt{\tau}}\f$.
   * @param tau   Time to maturity.
   * @param s     Spot ratio.
   * @param r     Risk-free rate.
   * @param sigma Volatility.
   */
  static double delta_minus(double tau, double s, double r, double sigma);

  /**
   * @brief Standard normal cumulative distribution function \f$\Phi(x)\f$.
   * @param x Input value.
   * @return \f$P(Z \le x)\f$ for \f$Z \sim \mathcal{N}(0,1)\f$.
   */
  static double normalCDF(double x);

  /**
   * @brief Standard normal probability density function \f$\phi(x)\f$.
   * @param x Input value.
   * @return \f$\frac{1}{\sqrt{2\pi}} e^{-x^2/2}\f$.
   */
  static double normalPDF(double x);
};
}
#endif
