/**
 * @file Pricer.cpp
 * @brief Implementation of the PricerLookbackOption class.
 *
 * Prices floating-strike lookback options via Monte Carlo simulation and
 * computes Greeks using central finite differences with common random numbers.
 */

#include "../include/Pricer.hpp"

#include <cmath>
#include <numeric>
namespace Chivet 
{

PricerLookbackOption::PricerLookbackOption() {}

PricerLookbackOption::PricerLookbackOption(LookbackOption option,
                                           MarketParameters marketParameters,
                                           MonteCarloSimulator simulator)
    : option_(option), MarketParameters_(marketParameters),
      simulator_(simulator) {}

void PricerLookbackOption::setOption(LookbackOption option) {
  option_ = option;
}
void PricerLookbackOption::setOption(OptionType type, double maturity) {
  option_.setType(type);
  option_.setMaturity(maturity);
}

void PricerLookbackOption::setMarketParameters(
    MarketParameters marketParameters) {
  MarketParameters_ = marketParameters;
}

void PricerLookbackOption::setMarketParameters(double spot, double rate,
                                               double volatility) {
  MarketParameters_.setSpot(spot);
  MarketParameters_.setRiskFreeRate(rate);
  MarketParameters_.setVolatility(volatility);
}

void PricerLookbackOption::setMonteCarloSimulator(unsigned long numSimulations,
                                                  unsigned long numTimeSteps,
                                                  unsigned long seed) {
  simulator_.setNumSimulations(numSimulations);
  simulator_.setNumTimeSteps(numTimeSteps);
  simulator_.setRandomSeed(seed);
}

const LookbackOption &PricerLookbackOption::getOption() const {
  return option_;
}
const MarketParameters &PricerLookbackOption::getMarketParameters() const {
  return MarketParameters_;
}
const MonteCarloSimulator &PricerLookbackOption::getSimulator() const {
  return simulator_;
}

void PricerLookbackOption::setBgkCorrection(bool enabled) {
  simulator_.setBgkCorrection(enabled);
}

/**
 * @brief Compute the standard error of the mean from a sample of values.
 * @param values Vector of per-simulation realisations.
 * @return Standard error = \f$\sqrt{\hat\sigma^2 / N}\f$.
 */
static double computeStdError(const std::vector<double> &values) {
  size_t N = values.size();
  double mean = std::accumulate(values.begin(), values.end(), 0.0) / N;
  double variance = std::accumulate(values.begin(), values.end(), 0.0,
                                    [mean](double sum, double val) {
                                      return sum + (val - mean) * (val - mean);
                                    }) /
                    (N - 1);
  return std::sqrt(variance / N);
}

// Stderror is passed by reference to return both price and its standard error
// from the same function.
// If storeValues is true, the individual discounted payoffs are stored in
// values, which allows computing standard errors for Greeks via finite
// differences on per-simulation values.
double PricerLookbackOption::computeRawPrice(double spot, double rate,
                                             double volatility, double maturity,
                                             double &stdError,
                                             unsigned int seed,
                                             std::vector<double> &values,
                                             bool storeValues) {
  simulator_.setRandomSeed(seed);
  // Simulate paths
  auto paths = simulator_.simulatePaths(spot, rate, volatility, maturity);
  // Compute payoffs
  std::vector<double> payoffs(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    if (option_.getType() == OptionType::Call) {
      payoffs[i] =
          option_.payoff(paths[i].minimumSpot, paths[i].spotAtMaturity);
    } else {
      payoffs[i] =
          option_.payoff(paths[i].maximumSpot, paths[i].spotAtMaturity);
    }
  }
  // Discount payoffs and compute price
  double discountFactor = std::exp(-rate * maturity);
  double price = discountFactor *
                 std::accumulate(payoffs.begin(), payoffs.end(), 0.0) /
                 payoffs.size();
  // Compute standard error
  double payoffMean =
      std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / payoffs.size();
  double payoffVariance =
      std::accumulate(payoffs.begin(), payoffs.end(), 0.0,
                      [payoffMean](double sum, double payoff) {
                        return sum +
                               (payoff - payoffMean) * (payoff - payoffMean);
                      }) /
      (payoffs.size() - 1);
  stdError = discountFactor * std::sqrt(payoffVariance / payoffs.size());

  // Store individual discounted payoffs if requested
  if (storeValues) {
    values.resize(payoffs.size());
    for (size_t i = 0; i < payoffs.size(); ++i) {
      values[i] = discountFactor * payoffs[i];
    }
  }

  return price;
}

PricingResult PricerLookbackOption::compute(bool withGreekStdErrors) {
  PricingResult result{};

  // needed to use the same randomness for all price computations to ensure that
  // the differences are due to parameter changes and not random noise
  unsigned int seed = simulator_.getRandomSeed();
  double stdError;

  // Vectors to store per-simulation discounted payoffs for Greek std errors.
  // Only populated when withGreekStdErrors is true.
  std::vector<double> baseValues, upSpotValues, downSpotValues;
  std::vector<double> forwardTimeValues, backwardTimeValues;
  std::vector<double> upRateValues, downRateValues;
  std::vector<double> upVolValues, downVolValues;
  // Unused placeholder for calls that don't need stored values
  std::vector<double> unused;

  // Compute price and standard error for the base parameters
  // (price std error is always computed via the stdError out-parameter)
  result.price = computeRawPrice(
      MarketParameters_.getSpot(), MarketParameters_.gtetRiskFreeRate(),
      MarketParameters_.getVolatility(), option_.getMaturity(), stdError, seed,
      withGreekStdErrors ? baseValues : unused, withGreekStdErrors);
  result.priceStd = stdError;

  // Compute Greeks using finite differences

  // Delta = dV/dS ~ (V(S+ds) - V(S-ds)) / (2*ds)
  double spotBumpAbs = MarketParameters_.getSpot() * SPOT_BUMP;
  double priceSpotUp = computeRawPrice(
      MarketParameters_.getSpot() * (1 + SPOT_BUMP),
      MarketParameters_.gtetRiskFreeRate(), MarketParameters_.getVolatility(),
      option_.getMaturity(), stdError, seed,
      withGreekStdErrors ? upSpotValues : unused, withGreekStdErrors);

  double priceSpotDown = computeRawPrice(
      MarketParameters_.getSpot() * (1 - SPOT_BUMP),
      MarketParameters_.gtetRiskFreeRate(), MarketParameters_.getVolatility(),
      option_.getMaturity(), stdError, seed,
      withGreekStdErrors ? downSpotValues : unused, withGreekStdErrors);

  result.delta = (priceSpotUp - priceSpotDown) / (2 * spotBumpAbs);

  // Gamma = d2V/dS2 ~ (V(S+ds) - 2*V(S) + V(S-ds)) / (ds^2)
  double spotBumpAbs2 = spotBumpAbs * spotBumpAbs;
  result.gamma =
      (priceSpotUp - 2 * result.price + priceSpotDown) / spotBumpAbs2;
  

  if (withGreekStdErrors) {
    size_t N = baseValues.size();
    std::vector<double> greekPerSim(N);

    // Delta std error
    for (size_t i = 0; i < N; ++i) {
      greekPerSim[i] =
          (upSpotValues[i] - downSpotValues[i]) / (2 * spotBumpAbs);
    }
    result.deltaStd = computeStdError(greekPerSim);

    // Gamma std error
    for (size_t i = 0; i < N; ++i) {
      greekPerSim[i] =
          (upSpotValues[i] - 2 * baseValues[i] + downSpotValues[i]) /
          spotBumpAbs2;
    }
    result.gammaStd = computeStdError(greekPerSim);
  } else {
    result.deltaStd = 0.0;
    result.gammaStd = 0.0;
  }

  // Theta = dV/dT ~ (V(T+dt) - V(T-dt)) / (2*dt)
  double priceAfter = computeRawPrice(
      MarketParameters_.getSpot(), MarketParameters_.gtetRiskFreeRate(),
      MarketParameters_.getVolatility(),
      option_.getMaturity() * (1.0 + TIME_BUMP), stdError, seed,
      withGreekStdErrors ? forwardTimeValues : unused, withGreekStdErrors);

  double priceBefore = computeRawPrice(
      MarketParameters_.getSpot(), MarketParameters_.gtetRiskFreeRate(),
      MarketParameters_.getVolatility(),
      option_.getMaturity() * (1.0 - TIME_BUMP), stdError, seed,
      withGreekStdErrors ? backwardTimeValues : unused, withGreekStdErrors);

  result.theta =
      (priceAfter - priceBefore) / (2 * option_.getMaturity() * TIME_BUMP);

  // Rho = dV/dr ~ (V(r+dr) - V(r-dr)) / (2*dr)
  double priceRateUp = computeRawPrice(
      MarketParameters_.getSpot(),
      MarketParameters_.gtetRiskFreeRate() * (1.0 + RATE_BUMP),
      MarketParameters_.getVolatility(), option_.getMaturity(), stdError, seed,
      withGreekStdErrors ? upRateValues : unused, withGreekStdErrors);

  double priceRateDown = computeRawPrice(
      MarketParameters_.getSpot(),
      MarketParameters_.gtetRiskFreeRate() * (1.0 - RATE_BUMP),
      MarketParameters_.getVolatility(), option_.getMaturity(), stdError, seed,
      withGreekStdErrors ? downRateValues : unused, withGreekStdErrors);

  result.rho = (priceRateUp - priceRateDown) /
               (2 * RATE_BUMP * MarketParameters_.gtetRiskFreeRate());

  // Vega = dV/dsigma ~ (V(sigma+dv) - V(sigma-dv)) / (2*dv)
  double volBumpAbs = MarketParameters_.getVolatility() * VOL_BUMP;
  double priceVolUp = computeRawPrice(
      MarketParameters_.getSpot(), MarketParameters_.gtetRiskFreeRate(),
      MarketParameters_.getVolatility() * (1.0 + VOL_BUMP),
      option_.getMaturity(), stdError, seed,
      withGreekStdErrors ? upVolValues : unused, withGreekStdErrors);

  double priceVolDown = computeRawPrice(
      MarketParameters_.getSpot(), MarketParameters_.gtetRiskFreeRate(),
      MarketParameters_.getVolatility() * (1.0 - VOL_BUMP),
      option_.getMaturity(), stdError, seed,
      withGreekStdErrors ? downVolValues : unused, withGreekStdErrors);

  result.vega = (priceVolUp - priceVolDown) / (2 * volBumpAbs);

  if (withGreekStdErrors) {
    size_t N = baseValues.size();
    std::vector<double> greekPerSim(N);

    // Theta std error
    for (size_t i = 0; i < N; ++i) {
      greekPerSim[i] = (forwardTimeValues[i] - backwardTimeValues[i]) /
                       (2 * option_.getMaturity() * TIME_BUMP);
    }
    result.thetaStd = computeStdError(greekPerSim);

    // Rho std error
    for (size_t i = 0; i < N; ++i) {
      greekPerSim[i] = (upRateValues[i] - downRateValues[i]) / (2 * RATE_BUMP);
    }
    result.rhoStd = computeStdError(greekPerSim);

    // Vega std error
    for (size_t i = 0; i < N; ++i) {
      greekPerSim[i] = (upVolValues[i] - downVolValues[i]) / (2 * volBumpAbs);
    }
    result.vegaStd = computeStdError(greekPerSim);
  } else {
    result.thetaStd = 0.0;
    result.rhoStd = 0.0;
    result.vegaStd = 0.0;
  }

  return result;
}

PriceAndDelta PricerLookbackOption::priceAndDeltaForSpot(double spot) {
  unsigned int seed = simulator_.getRandomSeed();
  double stdError;
  std::vector<double> unused;
  double price =
      computeRawPrice(spot, MarketParameters_.gtetRiskFreeRate(),
                      MarketParameters_.getVolatility(), option_.getMaturity(),
                      stdError, seed, unused, false);
  double error = stdError;
  // Compute Delta using finite difference
  double priceUp = computeRawPrice(
      spot * (1 + SPOT_BUMP), MarketParameters_.gtetRiskFreeRate(),
      MarketParameters_.getVolatility(), option_.getMaturity(), stdError, seed,
      unused, false);
  double priceDown = computeRawPrice(
      spot * (1 - SPOT_BUMP), MarketParameters_.gtetRiskFreeRate(),
      MarketParameters_.getVolatility(), option_.getMaturity(), stdError, seed,
      unused, false);
  double delta = (priceUp - priceDown) / (2 * spot * SPOT_BUMP);
  return {price, error, delta};
}
}