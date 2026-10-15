#include "lookback/random.h"

namespace Barnicha 
{

Rng::Rng(unsigned int seed)
: engine_(seed), normal_(0.0, 1.0),
  uniform_(0.0, 1.0) 
  {}

double Rng::next_normal() 
{
    return normal_(engine_);
}

double Rng::next_uniform() 
{
    // Avoid returning 0.0 exactly: inversion uses log(u).
    double u = uniform_(engine_);
    if (u == 0.0) 
    {
        u = 1e-16;
    }
    return u;
}

}
