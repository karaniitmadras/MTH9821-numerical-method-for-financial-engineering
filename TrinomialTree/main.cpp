#include "BinomialTree.hpp"
#include "BinomialTreeVolatilityUnknown.hpp"
#include "Calculator.hpp"
#include "EuropeanOption.hpp"
#include "TrinomialTree.hpp"
#include "TrinomialTreeVolatilityUnknown.hpp"
#include <iomanip>
#include <iostream>

void Problem1And2()
{
    Calculator calculator;

    std::cout << std::setprecision(6);
    std::cout << std::fixed;
    std::cout << "======= Problem 1 ========" << std::endl;
    std::cout << "European - BS" << std::endl;
    EuropeanOption option(0., 41., 39., 1., .25, .03, .005);
    std::cout << option.Put() << '\t' << option.DeltaPut() << '\t' << std::setprecision(7) << option.GammaPut() << '\t' << std::setprecision(6) << option.ThetaPut() << std::endl;
    calculator.PriceAndGreeks(option, false);

    std::cout << "\n\n\n======= Problem 2 ========" << std::endl;
    std::cout << "American - Exact" << std::endl;
    TrinomialTree tree_exact1(41, .25, 1., 10000, .03, .005);
    TrinomialTree tree_exact2(41, .25, 1., 10001, .03, .005);

    auto exact1 = tree_exact1.TreePricer(option.PutPayoff(), true);
    auto exact2 = tree_exact2.TreePricer(option.PutPayoff(), true);

    std::cout << (exact1.value + exact2.value) / 2. << '\t' << (exact1.delta + exact2.delta) / 2. << '\t' << std::setprecision(7) << (exact1.gamma + exact2.gamma) / 2. << '\t' << std::setprecision(6) << (exact1.theta + exact2.theta) / 2. << std::endl;
    calculator.PriceAndGreeks(option, true);

    std::cout << "\n\n\n======= Variance Reduction =======" << std::endl;
    calculator.PriceAndGreeks_compensated(option);
}

void Problem3()
{

    std::cout << std::setprecision(6);

    TrinomialTreeVolatilityUnknown trinomialTreeVolatilityUnknown(48., 50., 1. / 3., .02, .01, 4.08);

    trinomialTreeVolatilityUnknown.Put();
}

int main(int argc, const char* argv[])
{
    std::cout << std::fixed << std::setprecision(6) << std::endl;

    Problem1And2();
    Problem3();
    return 0;
}
