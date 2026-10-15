/**
 * @file OptionParameters.cpp
 * @brief Implementation of the LookbackOption class.
 */

#include "../include/OptionParameters.hpp"
#include <algorithm>
#include <stdexcept>
namespace Chivet 
{

LookbackOption::LookbackOption()
    : optionType_(OptionType::Call), maturity_(1.0) {}

LookbackOption::LookbackOption(OptionType optionType, double maturity)
    : optionType_(optionType), maturity_(maturity) {
  if (maturity <= 0.0) {
    throw std::invalid_argument("Maturity must be positive.");
  }
}

void LookbackOption::setType(OptionType optionType) {
  optionType_ = optionType;
}

void LookbackOption::setMaturity(double maturity) {
  if (maturity <= 0.0) {
    throw std::invalid_argument("Maturity must be positive.");
  }
  maturity_ = maturity;
}

OptionType LookbackOption::getType() const { return optionType_; }

double LookbackOption::getMaturity() const { return maturity_; }

double LookbackOption::payoff(double extremumOfSpot,
                              double SpotAtMaturity) const {
  // extremumOfSpot will need to be different for call and put options: for a
  // call option, it is the minimum
  if (optionType_ == OptionType::Call) {
    return std::max(0.0, SpotAtMaturity - extremumOfSpot);
  } else { // OptionType::Put
    return std::max(0.0, extremumOfSpot - SpotAtMaturity);
  }
}
}