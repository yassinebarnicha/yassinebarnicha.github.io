#include "lookback/lookback_option.h"

namespace Barnicha 
{

LookbackOption::LookbackOption(LookbackType type, double maturity)
: type_(type), maturity_(maturity) {}

double LookbackOption::payoff(double st, double min_s, double max_s) const 
{
    // Floating strike: call uses running minimum; put uses running maximum.
    return (type_ == LookbackType::Call) ? (st - min_s) : (max_s - st);
}

double LookbackOption::maturity() const { return maturity_; }
LookbackType LookbackOption::type() const { return type_; }

}
