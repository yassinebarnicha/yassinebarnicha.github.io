/**
 * @file Pricer.hpp
 * @brief Main pricing engine for floating-strike lookback options using Monte
 *        Carlo simulation with finite-difference Greeks.
 */

#ifndef PRICER_HPP
#define PRICER_HPP

#include "MarketParameters.hpp"
#include "MonteCarloSimulator.hpp"
#include "OptionParameters.hpp"

#include <vector>
namespace Chivet 
{

/**
 * @struct PriceAndDelta
 * @brief Lightweight result holding price and delta for a single spot value.
 *
 * Used by PricerLookbackOption::priceAndDeltaForSpot() for graph generation.
 */
struct PriceAndDelta {
  double price;    ///< Monte Carlo estimated option price.
  double priceStd; ///< Standard error of the price estimate.
  double delta;    ///< Finite-difference delta \f$\partial V/\partial S\f$.
};

/**
 * @struct PricingResult
 * @brief Full set of Monte Carlo pricing results including all first-order
 *        Greeks and their standard errors.
 */
struct PricingResult {
  double price;    ///< Option price \f$V\f$.
  double priceStd; ///< Standard error of the price.
  double delta;    ///< \f$\Delta = \partial V / \partial S\f$.
  double deltaStd; ///< Standard error of delta.
  double gamma;    ///< \f$\Gamma = \partial^2 V / \partial S^2\f$.
  double gammaStd; ///< Standard error of gamma.
  double theta;    ///< \f$\Theta = \partial V / \partial T\f$.
  double thetaStd; ///< Standard error of theta.
  double rho;      ///< \f$\rho = \partial V / \partial r\f$.
  double rhoStd;   ///< Standard error of rho.
  double vega;     ///< \f$\nu = \partial V / \partial \sigma\f$.
  double vegaStd;  ///< Standard error of vega.
};

/**
 * @class PricerLookbackOption
 * @brief Prices a floating-strike lookback option via Monte Carlo and computes
 *        Greeks using central finite differences with common random numbers.
 *
 * The pricer re-seeds the RNG before each price evaluation so that bumped
 * and base scenarios share the same Brownian paths, reducing noise in the
 * finite-difference Greek estimates.
 *
 * Remark : Giving a seed as the user to the program is useless and broken in
 * our implementation, as the seed is re-set before each price evaluation. To
 * get reproducible results, set the seed in the MonteCarloSimulator
 * configuration.
 */
class PricerLookbackOption {
private:
  LookbackOption option_;             ///< The lookback option contract.
  MarketParameters MarketParameters_; ///< Market data (spot, rate, vol).
  MonteCarloSimulator simulator_;     ///< Monte Carlo simulation engine.

  /// @name Bump sizes for central finite differences
  /// @{
  static constexpr double SPOT_BUMP =
      0.001; ///< Relative bump for Delta / Gamma.
  static constexpr double VOL_BUMP = 0.001;   ///< Relative bump for Vega.
  static constexpr double RATE_BUMP = 0.0001; ///< Relative bump for Rho.
  static constexpr double TIME_BUMP =
      1.0 / 365.0; ///< Relative bump for Theta (1 day).
  /// @}

  /**
   * @brief Compute the discounted Monte Carlo price for given parameters.
   *
   * This is the core pricing routine used internally by compute() and
   * priceAndDeltaForSpot(). It re-seeds the RNG with @p seed before
   * simulating, ensuring reproducible paths across bumped evaluations.
   *
   * @param spot        Spot price \f$S_0\f$.
   * @param rate        Risk-free rate \f$r\f$.
   * @param volatility  Volatility \f$\sigma\f$.
   * @param maturity    Time to maturity \f$T\f$.
   * @param[out] stdError Standard error of the price estimate.
   * @param seed        Random seed for reproducibility.
   * @param[out] values Per-simulation discounted payoffs (populated if
   *                    @p storeValues is true).
   * @param storeValues If true, individual discounted payoffs are stored in
   *                    @p values for subsequent Greek standard-error
   * computation. Since this is memory-intensive, it is optional and only used
   * when computing Greeks.
   * @return The Monte Carlo estimated option price.
   */
  double computeRawPrice(double spot, double rate, double volatility,
                         double maturity, double &stdError, unsigned int seed,
                         std::vector<double> &values, bool storeValues = false);

public:
  /** @brief Default constructor. */
  PricerLookbackOption();

  /**
   * @brief Construct a fully-specified pricer.
   * @param option           The lookback option to price.
   * @param marketParameters Market data (spot, rate, vol).
   * @param simulator        Configured Monte Carlo engine.
   */
  PricerLookbackOption(LookbackOption option, MarketParameters marketParameters,
                       MonteCarloSimulator simulator);

  /** @brief Set the option from an existing LookbackOption object. */
  void setOption(LookbackOption option);

  /**
   * @brief Set the option by type and maturity.
   * @param type     Call or Put.
   * @param maturity Time to maturity in years.
   */
  void setOption(OptionType type, double maturity);

  /** @brief Set market parameters from an existing object. */
  void setMarketParameters(MarketParameters marketParameters);

  /**
   * @brief Set market parameters individually.
   * @param spot       Spot price.
   * @param rate       Risk-free rate.
   * @param volatility Volatility.
   */
  void setMarketParameters(double spot, double rate, double volatility);

  /**
   * @brief Configure the Monte Carlo simulator.
   * @param numSimulations Number of paths.
   * @param numTimeSteps   Discrete time steps per path.
   * @param seed           Random seed (0 = random).
   */
  void setMonteCarloSimulator(unsigned long numSimulations,
                              unsigned long numTimeSteps,
                              unsigned long seed = 0);

  /** @brief Get a const reference to the option. */
  const LookbackOption &getOption() const;

  /** @brief Get a const reference to the market parameters. */
  const MarketParameters &getMarketParameters() const;

  /** @brief Get a const reference to the simulator. */
  const MonteCarloSimulator &getSimulator() const;

  /** @brief Enable or disable the Broadie-Glasserman-Kou continuity correction. */
  void setBgkCorrection(bool enabled);

  /**
   * @brief Price the option and compute all Greeks via Monte Carlo.
   *
   * Uses central finite differences with common random numbers.
   * @param withGreekStdErrors If true, store per-simulation discounted payoffs
   *        for every bumped evaluation and compute standard errors for each
   *        Greek. This roughly doubles memory usage and adds overhead.
   *        When false (default), Greek std-error fields are set to 0.
   *        The price standard error is always computed regardless of this flag.
   * @return A PricingResult containing price, Greeks, and their standard
   * errors.
   */
  PricingResult compute(bool withGreekStdErrors = false);

  /**
   * @brief Compute price and delta for a specific spot value.
   *
   * Useful for generating Price(S) and Delta(S) graphs.
   * @param spot The spot price to evaluate at.
   * @return PriceAndDelta containing MC price, its std error, and delta.
   */
  PriceAndDelta priceAndDeltaForSpot(double spot);
};
}
#endif // !PRICER_HPP
