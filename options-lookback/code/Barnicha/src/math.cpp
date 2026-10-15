#include "lookback/math.h"

namespace Barnicha 
{

double Math::normal_pdf(double x) 
{
    // phi(x) = (1/sqrt(2*pi)) * exp(-x^2/2)
    static const double inv_sqrt_2pi = 0.39894228040143267794;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

double Math::normal_cdf(double x) 
{
    // Phi(x) = 0.5 * erfc(-x/sqrt(2))
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

}
