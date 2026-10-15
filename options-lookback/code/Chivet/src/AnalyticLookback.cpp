/**
 * @file AnalyticLookback.cpp
 * @brief Implementation of the LookbackAnalyticSolution class.
 *
 * Contains closed-form pricing for floating-strike lookback options and
 * finite-difference Greeks computed from the analytic price.
 */

#include "../include/AnalyticLookback.hpp"
#include <cassert>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif
namespace Chivet 
{

double LookbackAnalyticSolution::normalCDF(double x) {
  return 0.5 * std::erfc(-x * M_SQRT1_2);
  // see https://en.cppreference.com/w/cpp/numeric/math/erfc.html
}

double LookbackAnalyticSolution::normalPDF(double x) {
  return (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * x * x);
}

double LookbackAnalyticSolution::delta_plus(double tau, double s, double r,
                                            double sigma) {
  assert(s > 0);
  return 1.0 / (sigma * std::sqrt(tau)) *
         (std::log(s) + (r + 0.5 * sigma * sigma) * tau);
}

double LookbackAnalyticSolution::delta_minus(double tau, double s, double r,
                                             double sigma) {
  assert(s > 0);
  return 1.0 / (sigma * std::sqrt(tau)) *
         (std::log(s) + (r - 0.5 * sigma * sigma) * tau);
}

double LookbackAnalyticSolution::price(OptionType option, double S0, double T,
                                       double sigma, double r) {
  assert(S0 > 0);
  assert(T > 0);
  assert(sigma > 0);
  double delta_p = delta_plus(T, 1.0, r, sigma);
  double delta_m = delta_minus(T, 1.0, r, sigma);
  if (option == OptionType::Call) {
    double V_discounted =
        S0 * normalCDF(delta_p) - S0 * std::exp(-r * T) * normalCDF(delta_m) +
        std::exp(-r * T) * (sigma * sigma * S0 / (2.0 * r)) *
            normalCDF(delta_m) -
        normalCDF(-delta_p) * S0 * (sigma * sigma / (2.0 * r));
    return V_discounted;
  } else { // Put
    return S0 * std::exp(-r * T) * normalCDF(-delta_m) +
           S0 * normalCDF(delta_p) * (1 + (sigma * sigma) / (2 * r)) -
           std::exp(-r * T) * S0 * (sigma * sigma / (2.0 * r)) *
               normalCDF(-delta_m) -
           S0;
  }
}

double LookbackAnalyticSolution::delta(OptionType option, double S0, double T,
                                       double sigma, double r) {
  double V_plus = price(option, S0 * 1.0001, T, sigma, r);
  double V_minus = price(option, S0 * 0.9999, T, sigma, r);
  return (V_plus - V_minus) / (S0 * 0.0002);
}

double LookbackAnalyticSolution::gamma(OptionType option, double S0, double T,
                                       double sigma, double r) {
  double V_plus = price(option, S0 * 1.0001, T, sigma, r);
  double V_minus = price(option, S0 * 0.9999, T, sigma, r);
  double V_0 = price(option, S0, T, sigma, r);
  return (V_plus - 2 * V_0 + V_minus) / (S0 * S0 * 0.0002 * 0.0002);
}

double LookbackAnalyticSolution::vega(OptionType option, double S0, double T,
                                      double sigma, double r) {
  double V_plus = price(option, S0, T, sigma * 1.0001, r);
  double V_minus = price(option, S0, T, sigma * 0.9999, r);
  return (V_plus - V_minus) / (sigma * 0.0002);
}

double LookbackAnalyticSolution::theta(OptionType option, double S0, double T,
                                       double sigma, double r) {
  double V_plus = price(option, S0, T * 1.0001, sigma, r);
  double V_minus = price(option, S0, T * 0.9999, sigma, r);
  return (V_plus - V_minus) / (T * 0.0002);
}

double LookbackAnalyticSolution::rho(OptionType option, double S0, double T,
                                     double sigma, double r) {
  double V_plus = price(option, S0, T, sigma, r * 1.0001);
  double V_minus = price(option, S0, T, sigma, r * 0.9999);
  return (V_plus - V_minus) / (r * 0.0002);
}
}