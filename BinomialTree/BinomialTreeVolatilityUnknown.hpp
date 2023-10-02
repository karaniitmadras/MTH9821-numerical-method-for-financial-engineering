#ifndef BinomialTreeVolatilityUnknown_hpp
#define BinomialTreeVolatilityUnknown_hpp

#include "BinomialTree.hpp"

class BinomialTreeVolatilityUnknown {
private:
    double S_; // Spot price
    double K_; // Strike price
    double T_; // Maturity
    double r_; // Const interest rate
    double q_; // Dividend rate
    double V0_; // Value at t=0

    std::size_t steps_; // Optimal steps

public:
    BinomialTreeVolatilityUnknown(double S, double K, double T, double r, double q, double V0);
    ~BinomialTreeVolatilityUnknown() = default;

    double PutBinomialTreeVolatilityUnknown() const;
};

#endif /* BinomialTreeVolatilityUnknown.hpp */
