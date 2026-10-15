#pragma once
#include <random>

namespace Barnicha 
{

/**
 * @brief RNG wrapper providing N(0,1) and U(0,1).
 *
 * We use `std::mt19937` (Mersenne Twister) because it is:
 * - part of the C++ standard library,
 * - reproducible across platforms (given the same seed),
 * - fast, with a very long period.
 *
 * Note: `next_uniform()` is clipped away from 0 to avoid `log(0)` in inversion formulas.
 */
class Rng {
public:
    /**
     * @param seed seed for reproducibility (same seed => same random stream)
     */
    explicit Rng(unsigned int seed = 42);

    /**
     * @brief Draw z ~ N(0,1).
     */
    double next_normal();

    /**
     * @brief Draw u ~ U(0,1) (never returns exactly 0).
     */
    double next_uniform();

private:
    std::mt19937 engine_;
    std::normal_distribution<double> normal_;
    std::uniform_real_distribution<double> uniform_;
};

}
