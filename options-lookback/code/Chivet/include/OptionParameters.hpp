/**
 * @file OptionParameters.hpp
 * @brief Defines the lookback option type and its payoff logic.
 */

#ifndef OPTION_PARAMETERS_HPP
#define OPTION_PARAMETERS_HPP
namespace Chivet 
{

/**
 * @enum OptionType
 * @brief Distinguishes between Call and Put options.
 */
enum class OptionType { Call, Put };

/**
 * @class LookbackOption
 * @brief Represents a floating-strike lookback option.
 *
 * A lookback option's payoff depends on the extremum (minimum or maximum)
 * of the underlying spot price observed over the life of the option.
 */
class LookbackOption {
private:
  OptionType optionType_; ///< Call or Put.
  double maturity_;       ///< Time to maturity in years.

public:
  /**
   * @brief Default constructor.
   *
   * Initialises as a Call with maturity = 1 year.
   */
  LookbackOption();

  /**
   * @brief Parameterised constructor.
   * @param optionType Call or Put.
   * @param maturity   Time to maturity in years (must be > 0).
   * @throws std::invalid_argument If maturity <= 0.
   */
  LookbackOption(OptionType optionType, double maturity);

  /**
   * @brief Set the option type.
   * @param optionType Call or Put.
   */
  void setType(OptionType optionType);

  /**
   * @brief Set the time to maturity.
   * @param maturity Time to maturity in years (must be > 0).
   * @throws std::invalid_argument If maturity <= 0.
   */
  void setMaturity(double maturity);

  /** @brief Get the time to maturity in years. */
  double getMaturity() const;

  /** @brief Get the option type. */
  OptionType getType() const;

  /**
   * @brief Compute the lookback option payoff.
   *
   * - **Call**: max(0, SpotAtMaturity - min(Spot over life))
   * - **Put**:  max(0, max(Spot over life) - SpotAtMaturity)
   *
   * @param extremumOfSpot The path minimum (call) or maximum (put) of the
   *                       underlying spot price during the option's life.
   * @param SpotAtMaturity The spot price at maturity.
   * @return The non-negative payoff.
   */
  double payoff(double extremumOfSpot, double SpotAtMaturity) const;
};
}
#endif // OPTION_PARAMETERS_HPP
