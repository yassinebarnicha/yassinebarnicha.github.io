/**
 * @file MarketParameters.hpp
 * @brief Class for market data used for option pricing.
 */

#ifndef MARKET_PARAMETERS_HPP
#define MARKET_PARAMETERS_HPP
namespace Chivet 
{

/**
 * @class MarketParameters
 * @brief Stores the market parameters required for pricing: spot price,
 *        risk-free interest rate, and volatility.
 */
class MarketParameters {
public:
  /**
   * @brief Default constructor.
   *
   * Initialises with spot = 100, risk-free rate = 5%, volatility = 100%.
   */
  MarketParameters();

  /**
   * @brief Parameterised constructor.
   * @param spot        Current spot price of the underlying (must be > 0).
   * @param riskFreeRate Continuously-compounded risk-free interest rate.
   * @param volatility  Annualised volatility of the underlying (must be >= 0).
   * @throws std::invalid_argument If spot <= 0 or volatility < 0.
   */
  MarketParameters(double spot, double riskFreeRate, double volatility);

  ~MarketParameters() = default;

  /**
   * @brief Set the spot price.
   * @param spot New spot price (must be > 0).
   * @throws std::invalid_argument If spot <= 0.
   */
  void setSpot(double spot);

  /**
   * @brief Set the risk-free interest rate.
   * @param riskFreeRate New risk-free rate.
   */
  void setRiskFreeRate(double riskFreeRate);

  /**
   * @brief Set the volatility.
   * @param volatility New volatility (must be >= 0).
   * @throws std::invalid_argument If volatility < 0.
   */
  void setVolatility(double volatility);

  /** @brief Get the current spot price. */
  double getSpot() const;

  /** @brief Get the risk-free interest rate. */
  double gtetRiskFreeRate() const;

  /** @brief Get the annualised volatility. */
  double getVolatility() const;

private:
  double spot_;         ///< Current spot price of the underlying.
  double riskFreeRate_; ///< Risk-free rate.
  double volatility_;   ///< Volatility.
};
}
#endif // !MARKET_PARAMETERS_HPP
