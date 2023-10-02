#include "Calculator.hpp"

void Calculator::PriceAndGreeks(const EuropeanOption& option, bool american)
{
    std::cout << "\n\n\n"
              << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricer(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }

    std::cout << "\n\n\n"
              << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricerBS(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }

    std::cout << "\n\n\n"
              << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricerBSR(option.PutPayoff(), american);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
}

void Calculator::PriceAndGreeks_compensated(const EuropeanOption& option)
{
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricer_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }

    std::cout << "\n\n\n"
              << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricerBS_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }

    std::cout << "\n\n\n"
              << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(41, .25, 1., i, .03, .005);
        auto res = tree.TreePricerBSR_PathDependent_compensated(option.PutPayoff(), option);
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(7) << res.gamma << std::setprecision(6) << '\t' << res.theta << std::endl;
    }
}