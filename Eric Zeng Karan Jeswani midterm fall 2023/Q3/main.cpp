#include "Calculator.hpp"
#include "EuropeanOption.hpp"
#include "TrinomialTree.hpp"
#include "TrinomialTreeVolatilityUnknown.hpp"
#include <iomanip>
#include <iostream>

void Problem1And2()
{
    Calculator calculator;

    std::cout << std::setprecision(4);
    std::cout << std::fixed;
    EuropeanOption option(0., 33., 32., 10./12., .24, .045, .02);

    std::cout << "\n\n\n======= Problem 2 ========" << std::endl;
    std::cout << "American - Exact" << std::endl;
    calculator.PriceAndGreeks(option, true);
}

void Problem3()
{

    std::cout << std::setprecision(4);

    TrinomialTreeVolatilityUnknown trinomialTreeVolatilityUnknown(48., 50., 1. / 3., .02, .01, 4.08);

    trinomialTreeVolatilityUnknown.Put();
}

int main(int argc, const char* argv[])
{
    std::cout << std::fixed << std::setprecision(4) << std::endl;

    Problem1And2();
    // Problem3();
    return 0;
}
