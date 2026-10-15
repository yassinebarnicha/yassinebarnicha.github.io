#pragma once
/**
 * @file greeks_calculator.h
 * @brief Declaration of the abstract interface for option Greeks computation.
 */
#include "greeks_result.h"
#include "pricing_engine.h"

namespace Barnicha 
{

/**
 * @brief Abstract interface for Greeks computation.
 *
 * The engine is injected, so the same calculator can run on:
 * - MonteCarloEngine (CRN because it reuses stored samples),
 * - AnalyticEngine (deterministic).
 */
class GreeksCalculator {
public:
    virtual ~GreeksCalculator() = default;

    /**
     * @brief Compute requested outputs using the provided pricing engine.
     */
    virtual GreeksResult compute(const PricingEngine& engine,
                                 const MarketData& market,
                                 const LookbackOption& option) const = 0;
};

}
