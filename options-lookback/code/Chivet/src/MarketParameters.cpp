/**
 * @file MarketParameters.cpp
 * @brief Implementation of the MarketParameters class.
 */

#include "../include/MarketParameters.hpp"
#include <stdexcept>
namespace Chivet 
{

MarketParameters::MarketParameters()
    : spot_(100.0), riskFreeRate_(0.05), volatility_(1.00) {}

MarketParameters::MarketParameters(double spot, double riskFreeRate,
                                   double volatility)
    : spot_(spot), riskFreeRate_(riskFreeRate), volatility_(volatility) {
  if (spot <= 0.0) {
    throw std::invalid_argument("Spot price must be (strictly) positive.");
  }
  if (volatility < 0.0) {
    throw std::invalid_argument("Volatility must be non-negative.");
  }
}

void MarketParameters::setSpot(double spot) {
  if (spot <= 0.0) {
    throw std::invalid_argument("Spot price must be (strictly) positive.");
  }
  spot_ = spot;
}

void MarketParameters::setRiskFreeRate(double riskFreeRate) {
  riskFreeRate_ = riskFreeRate;
}

void MarketParameters::setVolatility(double volatility) {
  if (volatility < 0.0) {
    throw std::invalid_argument("Volatility must be non-negative.");
  }
  volatility_ = volatility;
}

double MarketParameters::getSpot() const { return spot_; }

double MarketParameters::gtetRiskFreeRate() const { return riskFreeRate_; }

double MarketParameters::getVolatility() const { return volatility_; }
}