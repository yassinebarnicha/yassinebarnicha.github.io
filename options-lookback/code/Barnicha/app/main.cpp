#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "lookback/one_step_bridge_sampler.h"
#include "lookback/monte_carlo_engine.h"
#include "lookback/analytic_engine.h"
#include "lookback/bump_and_revalue_greeks.h"

using namespace lookback;

static double to_double(const char* s) { return std::atof(s); }
static long to_long (const char* s) { return std::atol(s); }

static void print_csv_block(const char* engineName, const GreeksResult& r) 
{
    std::cout << "engine,metric,estimate,stderr,ci95_low,ci95_high\n";
    auto row = [&](const char* metric, const EstimateStats& s) {
        std::cout << engineName << "," << metric << ","
        << s.estimate << "," << s.std_dev << ","
        << s.ci95_low << "," << s.ci95_high << "\n";
    };
    row("price", r.price);
    row("delta", r.delta);
    row("gamma", r.gamma);
    row("theta", r.theta);
    row("rho", r.rho);
    row("vega", r.vega);
}

int main(int argc, char** argv) 
{
    // Defaults (ok for manual run; VBA overrides them)
    double ttm = 1.0;
    double spot = 100.0;
    double rate = 0.05;
    double vol = 0.20;
    long paths = 200000;
    long seed = 12345;
    std::string type = "Call";
    bool csv = false;

    // Minimal flag parsing: expects "--key value" (except --csv)
    for (int i = 1; i < argc; ++i) 
    {
        std::string k = argv[i];

        if (k == "--csv") 
        {
            csv = true;
        }
        else if (k == "--ttm" && i + 1 < argc) 
        {
            ttm = to_double(argv[++i]);
        } 
        else if (k == "--spot" && i + 1 < argc) 
        {
            spot = to_double(argv[++i]);
        } 
        else if (k == "--rate" && i + 1 < argc) 
        {
            rate = to_double(argv[++i]);
        } 
        else if (k == "--vol" && i + 1 < argc) 
        {
            vol = to_double(argv[++i]);
        }
        else if (k == "--paths" && i + 1 < argc) 
        {
            paths = to_long(argv[++i]);
        } 
        else if (k == "--seed" && i + 1 < argc) 
        {
            seed = to_long(argv[++i]);
        } 
        else if (k == "--type" && i + 1 < argc) 
        {
            type = argv[++i];
        } 
        else if (k == "--help" || k == "-h") 
        {
            std::cout
            << "Usage:\n"
            << " lookback_app --ttm <years> [--type Call|Put] [--spot S0] [--rate r] [--vol sigma]\n"
            << " [--paths N] [--seed seed] [--csv]\n";
            return 0;
        } 
        else 
        {
            std::cerr << "Unknown or incomplete argument: " << k << "\n";
            return 1;
        }
    }

    if (ttm <= 0.0) 
    {
        std::cerr << "Error: --ttm must be positive.\n";
        return 1;
    }

    const LookbackType optType =
    (type == "Put" || type == "put") ? LookbackType::Put : LookbackType::Call;

    const MarketData market{spot, rate, vol};
    const LookbackOption option(optType, ttm);

    OneStepBridgeSampler sampler((unsigned)seed);
    sampler.generate((std::size_t)paths);
    MonteCarloEngine mc_engine(sampler, (std::size_t)paths, true);

    AnalyticEngine an_engine;

    BumpAndRevalueGreeks greeks(true, true, true, true, true, true);

    const GreeksResult mc = greeks.compute(mc_engine, market, option);
    const GreeksResult an = greeks.compute(an_engine, market, option);

    std::cout << std::fixed << std::setprecision(10);

    if (csv) 
    {
        print_csv_block("MonteCarlo", mc);
        std::cout << "\n";
        print_csv_block("Analytic", an);
    } 
    else 
    {
        std::cout << "ttm=" << ttm << "\n";
        std::cout << "MonteCarlo price=" << mc.price.estimate << " (se=" << mc.price.std_dev << ")\n";
        std::cout << "Analytic price=" << an.price.estimate << "\n";
    }

    return 0;
}