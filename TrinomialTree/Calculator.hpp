#ifndef Calculator_hpp
#define Calculator_hpp

#include "BinomialTree.hpp"
#include "BinomialTreeVolatilityUnknown.hpp"
#include "EuropeanOption.hpp"
#include "TrinomialTree.hpp"
#include "TrinomialTreeVolatilityUnknown.hpp"
#include <iomanip>
#include <iostream>

class Calculator {
public:
    Calculator() {};
    ~Calculator() {};

    void PriceAndGreeks(const EuropeanOption& option, bool american);
    void PriceAndGreeks_compensated(const EuropeanOption& option);
};

#endif // Calculator.hpp