#ifndef TrinomialTreeVolatilityUnknown_hpp
#define TrinomialTreeVolatilityUnknown_hpp

#include "TrinomialTree.hpp"

class TrinomialTreeVolatilityUnknown {
private:
    double S_; // Spot price
    double K_; // Strike price
    double T_; // Maturity
    double r_; // Const interest rate
    double q_; // Dividend rate
    double V0_; // Value at t=0

    std::size_t steps_; // Optimal steps

    inline static bool american_ = true; // American IV

public:
    TrinomialTreeVolatilityUnknown(double S, double K, double T, double r, double q, double V0);
    ~TrinomialTreeVolatilityUnknown() = default;

    double Put() const;
};

#endif // TrinomialTreeVolatilityUnknown.hpp