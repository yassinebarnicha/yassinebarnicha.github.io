#include "lookback/one_step_bridge_sampler.h"

namespace Barnicha 
{

OneStepBridgeSampler::OneStepBridgeSampler(unsigned int seed)
: rng_(seed) 
{}

void OneStepBridgeSampler::generate(std::size_t n_samples) 
{
    // Store all samples so we can reuse them across all bumps (CRN).
    samples_.clear();
    samples_.reserve(n_samples);

    for (std::size_t i = 0; i < n_samples; ++i) 
    {
        samples_.push_back(OneStepSample{
            rng_.next_normal(),
            rng_.next_uniform(),
            rng_.next_uniform()
        });
    }
}

const std::vector<OneStepSample>& OneStepBridgeSampler::samples() const 
{
    return samples_;
}

}
