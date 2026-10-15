/**
 * @file cli.cpp
 * @brief Command-line interface for the lookback option pricer.
 *
 * Parses user-supplied market and simulation parameters, prices a
 * floating-strike lookback option via Monte Carlo, and prints the
 * results alongside the analytic solution for comparison.
 */

#include "include/AnalyticLookback.hpp"
#include "include/Pricer.hpp"

#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
using namespace Chivet ;
/**
 * @brief Print usage / help message to stdout.
 * @param programName The name of the executable (argv[0]).
 */
void printUsage(const char *programName) {
  std::cout << "Usage: " << programName << " [options]" << std::endl;
  std::cout << std::endl;
  std::cout << "Lookback Option Pricer - Monte Carlo Method" << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout //
      << "  --type TYPE           Option type: call or put (default: call)"
      << std::endl;
  std::cout //
      << "  -s, --spot PRICE      Current spot price (default: 100.0)"
      << std::endl;
  std::cout //
      << "  -t, --time DATE       Maturity date (YYYY-MM-DD), computes time to"
      << std::endl;
  std::cout //
      << "                        maturity from today (overrides -T)" << std::endl;
  std::cout //
      << "  -T, --maturity TIME   Time to maturity in years (default: 1.0)"
      << std::endl;
  std::cout //
      << "  -r, --rate RATE       Risk-free interest rate (default: 0.05)"
      << std::endl;
  std::cout                                                  //
      << "  -v, --vol SIGMA       Volatility (default: 1.0)" //
      << std::endl;
  std::cout //
      << "  -n, --sims NUM        Number of MC simulations (default: 100000)"
      << std::endl;
  std::cout //
      << "  -m, --steps NUM       Number of time steps per path (default: 100)"
      << std::endl;
  std::cout //
      << "  --seed NUM            Random seed (default: 0 = random)"
      << std::endl;
  std::cout //
      << "  --std_error           Compute standard errors for Greeks"
      << std::endl;
  std::cout //
      << "                        (slower, more memory; off by default)"
      << std::endl;
  std::cout //
      << "  --bgk                 Apply Broadie-Glasserman-Kou continuity"
      << std::endl;
  std::cout //
      << "                        correction for discrete monitoring"
      << std::endl;
  std::cout << "  -h, --help            Show this help message" << std::endl;
  std::cout << std::endl;
  std::cout << "Example:" << std::endl;
  std::cout << "  " << programName
            << " --type call -s 100 -t 2027-01-15 -r 0.05 -v 0.20 -n 100000"
            << std::endl;
}

/**
 * @brief Entry point for the lookback option pricer CLI.
 *
 * Accepts command-line arguments for option type, market parameters, and
 * Monte Carlo settings. Outputs Monte Carlo price/Greeks and the
 * corresponding analytic solution.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return 0 on success, 1 on invalid input.
 */
int main(int argc, char *argv[]) {
  // Default parameters
  OptionType type = OptionType::Call;
  double spot = 100.0;
  double maturity = 1.0;
  double rate = 0.05;
  double volatility = 1.0;
  unsigned long numSims = 100000;
  int numSteps = 100;
  unsigned int seed = 0;
  bool timeGiven = false;
  bool greekStdError = false;
  bool bgkCorrection = false;

  // Parse command line arguments
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      return 0;
    } else if (arg == "--type" && i + 1 < argc) {
      std::string typeStr = argv[++i];
      if (typeStr == "call" || typeStr == "Call" || typeStr == "CALL") {
        type = OptionType::Call;
      } else if (typeStr == "put" || typeStr == "Put" || typeStr == "PUT") {
        type = OptionType::Put;
      } else {
        std::cerr << "Error: Invalid option type '" << typeStr << "'"
                  << std::endl;
        return 1;
      }
    } else if ((arg == "-t" || arg == "--time") && i + 1 < argc) {
      std::string dateStr = argv[++i];
      std::istringstream ss(dateStr);
      int year, month, day;
      char dash1, dash2;
      ss >> year >> dash1 >> month >> dash2 >> day;
      if (ss.fail() || dash1 != '-' || dash2 != '-') {
        std::cerr << "Error: Invalid date format '" << dateStr
                  << "'. Expected YYYY-MM-DD" << std::endl;
        return 1;
      }
      std::tm target = {};
      target.tm_year = year - 1900;
      target.tm_mon = month - 1;
      target.tm_mday = day;
      target.tm_hour = 0;
      target.tm_min = 0;
      target.tm_sec = 0;
      target.tm_isdst = -1;
      std::time_t targetTime = std::mktime(&target);

      std::time_t now = std::time(nullptr);
      std::tm today = *std::localtime(&now);
      today.tm_hour = 0;
      today.tm_min = 0;
      today.tm_sec = 0;
      today.tm_isdst = -1;
      std::time_t todayTime = std::mktime(&today);

      double diffDays = std::difftime(targetTime, todayTime) / (60 * 60 * 24);
      if (diffDays <= 0) {
        std::cerr << "Error: Date '" << dateStr << "' must be in the future"
                  << std::endl;
        return 1;
      }
      maturity = diffDays / 365.0;
      timeGiven = true;
    } else if ((arg == "-s" || arg == "--spot") && i + 1 < argc) {
      spot = std::atof(argv[++i]);
    } else if ((arg == "-T" || arg == "--maturity") && i + 1 < argc) {
      if (!timeGiven) {
        maturity = std::atof(argv[++i]);
      } else {
        std::cerr << "Warning: -T/--maturity ignored because -t/--time was "
                     "already provided"
                  << std::endl;
        ++i;
      }
    } else if ((arg == "-r" || arg == "--rate") && i + 1 < argc) {
      rate = std::atof(argv[++i]);
    } else if ((arg == "-v" || arg == "--vol") && i + 1 < argc) {
      volatility = std::atof(argv[++i]);
    } else if ((arg == "-n" || arg == "--sims") && i + 1 < argc) {
      numSims = static_cast<unsigned long>(std::atol(argv[++i]));
    } else if ((arg == "-m" || arg == "--steps") && i + 1 < argc) {
      numSteps = std::atoi(argv[++i]);
    } else if (arg == "--seed" && i + 1 < argc) {
      seed = static_cast<unsigned int>(std::atoi(argv[++i]));
    } else if (arg == "--std_error") {
      greekStdError = true;
    } else if (arg == "--bgk") {
      bgkCorrection = true;
    } else {
      std::cerr << "Error: Unknown argument '" << arg << "'" << std::endl;
      printUsage(argv[0]);
      return 1;
    }
  }

  // Validate parameters
  if (spot <= 0) {
    std::cerr << "Error: Spot price must be positive" << std::endl;
    return 1;
  }
  if (maturity <= 0) {
    std::cerr << "Error: Maturity must be positive" << std::endl;
    return 1;
  }
  if (volatility < 0) {
    std::cerr << "Error: Volatility must be non-negative" << std::endl;
    return 1;
  }
  if (numSims < 100) {
    std::cerr << "Error: Number of simulations must be at least 100"
              << std::endl;
    return 1;
  }
  if (numSteps < 1) {
    std::cerr << "Error: Number of time steps must be at least 1" << std::endl;
    return 1;
  }

  // Create pricer and compute
  PricerLookbackOption pricer;
  pricer.setOption(type, maturity);
  pricer.setMarketParameters(spot, rate, volatility);
  pricer.setMonteCarloSimulator(numSims, numSteps, seed);
  if (bgkCorrection) {
    pricer.setBgkCorrection(true);
  }

  std::cout << "======================================" << std::endl;
  std::cout << "Lookback Option Pricer (Monte Carlo)" << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << std::endl;
  std::cout << "Input Parameters:" << std::endl;
  std::cout << "  Option Type:  "
            << (type == OptionType::Call ? "Call Lookback" : "Put Lookback")
            << std::endl;
  std::cout << "  Spot Price:   " << spot << std::endl;
  std::cout << "  Maturity:     " << maturity << " years" << std::endl;
  std::cout << "  Risk-free Rate: " << rate * 100.0 << "%" << std::endl;
  std::cout << "  Volatility:   " << volatility * 100.0 << "%" << std::endl;
  std::cout << "  Simulations:  " << numSims << std::endl;
  std::cout << "  Time Steps:   " << numSteps << std::endl;
  if (bgkCorrection)
    std::cout << "  BGK Correction: enabled" << std::endl;
  std::cout << std::endl;

  std::cout << "Computing..." << std::endl;
  PricingResult result = pricer.compute(greekStdError);

  std::cout << std::endl;
  std::cout << "Results:" << std::endl;
  std::cout << "--------------------------------------" << std::endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "  Price:  " << std::setw(12) << result.price << "  (+/- "
            << result.priceStd << ")" << std::endl;
  if (greekStdError) {
    std::cout << "  Delta:  " << std::setw(12) << result.delta << "  (+/- "
              << result.deltaStd << ")" << std::endl;
    std::cout << "  Gamma:  " << std::setw(12) << result.gamma << "  (+/- "
              << result.gammaStd << ")" << std::endl;
    std::cout << "  Theta:  " << std::setw(12) << result.theta << "  (+/- "
              << result.thetaStd << ")" << std::endl;
    std::cout << "  Rho:    " << std::setw(12) << result.rho << "  (+/- "
              << result.rhoStd << ")" << std::endl;
    std::cout << "  Vega:   " << std::setw(12) << result.vega << "  (+/- "
              << result.vegaStd << ")" << std::endl;
  } else {
    std::cout << "  Delta:  " << std::setw(12) << result.delta << std::endl;
    std::cout << "  Gamma:  " << std::setw(12) << result.gamma << std::endl;
    std::cout << "  Theta:  " << std::setw(12) << result.theta << std::endl;
    std::cout << "  Rho:    " << std::setw(12) << result.rho << std::endl;
    std::cout << "  Vega:   " << std::setw(12) << result.vega << std::endl;
  }
  std::cout << "--------------------------------------" << std::endl;

  std::cout << std::endl;
  std::cout << "Analytic Solution (for reference):" << std::endl;
  double analyticPrice =
      LookbackAnalyticSolution::price(type, spot, maturity, volatility, rate);
  std::cout << "  Price:  " << std::setw(12) << analyticPrice << std::endl;
  double analyticDelta =
      LookbackAnalyticSolution::delta(type, spot, maturity, volatility, rate);
  std::cout << "  Delta:  " << std::setw(12) << analyticDelta << std::endl;
  double analyticGamma =
      LookbackAnalyticSolution::gamma(type, spot, maturity, volatility, rate);
  std::cout << "  Gamma:  " << std::setw(12) << analyticGamma << std::endl;
  double analyticTheta =
      LookbackAnalyticSolution::theta(type, spot, maturity, volatility, rate);
  std::cout << "  Theta:  " << std::setw(12) << analyticTheta << std::endl;
  double analyticRho =
      LookbackAnalyticSolution::rho(type, spot, maturity, volatility, rate);
  std::cout << "  Rho:    " << std::setw(12) << analyticRho << std::endl;
  double analyticVega =
      LookbackAnalyticSolution::vega(type, spot, maturity, volatility, rate);
  std::cout << "  Vega:   " << std::setw(12) << analyticVega << std::endl;

  return 0;
}
