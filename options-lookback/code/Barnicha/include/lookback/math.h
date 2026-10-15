#pragma once
#include <cmath>

namespace Barnicha 
{

/**
 * @brief Small math helpers (standard normal PDF/CDF).
 */
class Math 
{
public:
    /**
     * @brief Standard normal PDF: phi(x).
     */
    static double normal_pdf(double x);

    /**
     * @brief Standard normal CDF: Phi(x).
     *
     * Implemented as Phi(x) = 0.5 * erfc(-x / sqrt(2)) for good numerical stability.
     */
    static double normal_cdf(double x);
};

}
