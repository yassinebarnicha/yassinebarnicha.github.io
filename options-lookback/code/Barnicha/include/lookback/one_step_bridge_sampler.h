#pragma once
/**
 * @file one_step_bridge_sampler.h
 * @brief Declaration of the one-step Brownian bridge random sampler.
 */
#include "random.h"
#include <cstddef>
#include <vector>

namespace Barnicha {

/**
 * @brief Random draws needed for one one-step lookback path.
 *
 * We separate the *generation* of randomness from the *pricing*:
 * the sampler stores samples so that the same random numbers can be reused
 * for price and for Greeks (common random numbers / CRN).
 */
struct OneStepSample 
{
    double z;     ///< standard normal for terminal value
    double u_max; ///< uniform for sampling the conditional maximum
    double u_min; ///< uniform for sampling the conditional minimum
};

/**
 * @brief Generates and stores samples once (CRN across bumps).
 */
class OneStepBridgeSampler 
{
public:
    /**
     * @param seed RNG seed
     */
    explicit OneStepBridgeSampler(unsigned int seed = 42);

    /**
     * @brief Generate n_samples and store them internally.
     */
    void generate(std::size_t n_samples);

    /**
     * @brief Access stored samples (read-only).
     */
    const std::vector<OneStepSample>& samples() const;

private:
    Rng rng_;
    std::vector<OneStepSample> samples_;
};

}
