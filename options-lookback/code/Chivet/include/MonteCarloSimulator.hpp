/**
 * @file MonteCarloSimulator.hpp
 * @brief Monte Carlo simulation engine for generating spot-price paths.
 */

#ifndef MONTE_CARLO_SIMULATOR_HPP
#define MONTE_CARLO_SIMULATOR_HPP

#include <random>
#include <vector>
namespace Chivet 
{

/**
 * @struct SpotResults
 * @brief Holds the outcome of a single simulated spot-price path.
 */
struct SpotResults {
  double spotAtMaturity; ///< Spot price at the end of the path.
  double maximumSpot;    ///< Maximum spot price observed along the path.
  double minimumSpot;    ///< Minimum spot price observed along the path.
};

/**
 * @class MonteCarloSimulator
 * @brief Generates spot-price paths under Geometric Brownian Motion (GBM).
 *
 * The underlying dynamics are:
 * \f[
 *   S_{t+dt} = S_t \exp\!\bigl[(r - \tfrac{1}{2}\sigma^2)\,dt
 *              + \sigma\sqrt{dt}\,Z\bigr], \quad Z \sim \mathcal{N}(0,1)
 * \f]
 *
 * Each path tracks the spot at maturity and the running maximum / minimum,
 * which are needed for lookback option payoffs.
 */
class MonteCarloSimulator {
private:
  unsigned long numSimulations_; ///< Number of Monte Carlo paths.
  unsigned long numTimeSteps_;   ///< Number of discrete time steps per path.
  mutable std::mt19937 rng_;     ///< Mersenne Twister PRNG.
  mutable std::normal_distribution<double>
      normDist_; ///< Standard normal distribution N(0, 1).
  bool bgkCorrection_ = false; ///< Apply Broadie-Glasserman-Kou continuity correction.

  /// Broadie-Glasserman-Kou constant: \f$\beta_1 = -\zeta(1/2)/\sqrt{2\pi} \approx 0.5826\f$.
  static constexpr double BGK_BETA1 = 0.5826;

public:
  /**
   * @brief Default constructor.
   *
   * Initialises with 10 000 simulations, 1 000 time steps, and a
   * random seed derived from the system clock.
   */
  MonteCarloSimulator();

  /**
   * @brief Parameterised constructor.
   * @param numSimulations Number of Monte Carlo paths to generate.
   * @param numTimeSteps   Number of discrete time steps per path.
   * @param seed           Random seed (0 = derive from system clock).
   */
  MonteCarloSimulator(unsigned long numSimulations, unsigned long numTimeSteps,
                      unsigned int seed = 0);

  /** @brief Set the number of simulations. */
  void setNumSimulations(unsigned long numSimulations);

  /** @brief Set the number of time steps per path. */
  void setNumTimeSteps(unsigned long numTimeSteps);

  /** @brief Re-seed the random number generator. */
  void setRandomSeed(unsigned int seed);

  /** @brief Enable or disable the Broadie-Glasserman-Kou continuity correction. */
  void setBgkCorrection(bool enabled);

  /** @brief Check whether the BGK correction is enabled. */
  bool getBgkCorrection() const;

  /** @brief Get the number of simulations. */
  unsigned long getNumSimulations() const;

  /** @brief Get the number of time steps per path. */
  unsigned long getNumTimeSteps() const;

  /** @brief Get a seed value derived from the current clock (snapshot). */
  unsigned int getRandomSeed() const;

  /**
   * @brief Simulate a single GBM price path.
   * @param initialSpot    Initial spot price \f$S_0\f$.
   * @param rate           Risk-free interest rate \f$r\f$.
   * @param volatility     Volatility \f$\sigma\f$.
   * @param timeToMaturity Time to maturity \f$T\f$ in years.
   * @return SpotResults containing spot at maturity, path max, and path min.
   */
  SpotResults simulateSinglePath(double initialSpot, double rate,
                                 double volatility,
                                 double timeToMaturity) const;

  /**
   * @brief Simulate all Monte Carlo paths.
   * @param initialSpot    Initial spot price \f$S_0\f$.
   * @param rate           Risk-free interest rate \f$r\f$.
   * @param volatility     Volatility \f$\sigma\f$.
   * @param timeToMaturity Time to maturity \f$T\f$ in years.
   * @return Vector of SpotResults, one per simulation.
   */
  std::vector<SpotResults> simulatePaths(double initialSpot, double rate,
                                         double volatility,
                                         double timeToMaturity) const;

  /**
   * @brief Draw a single sample from the standard normal distribution.
   * @return A realisation of \f$Z \sim \mathcal{N}(0,1)\f$.
   */
  double generateStandardNormal() const;
};
}
#endif // !MONTE_CARLO_SIMULATOR_HPP
