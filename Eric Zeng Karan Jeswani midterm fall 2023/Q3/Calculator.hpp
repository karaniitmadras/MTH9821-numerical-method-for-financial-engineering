#ifndef Calculator_hpp
#define Calculator_hpp

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
};

#endif // Calculator.hpp