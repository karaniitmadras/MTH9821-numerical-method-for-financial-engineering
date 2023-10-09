#include "Calculator.hpp"

void Calculator::PriceAndGreeks(const EuropeanOption& option, bool american)
{

    std::cout << "\n\n\n" << "American Call" 
              << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(33, .24, 10./12., i, .045, .02);
        auto res = tree.TreePricer(option.CallPayoff(), american);

        // vega
        double dsigma = .24 / 1e4;
        TrinomialTree treeSigmaPlus(33, .24+dsigma, 10./12., i, .045, .02);
        auto res_treeSigmaPlus = treeSigmaPlus.TreePricer(option.CallPayoff(), american);
        TrinomialTree treeSigmaMinus(33, .24-dsigma, 10./12., i, .045, .02);
        auto res_treeSigmaMinus = treeSigmaMinus.TreePricer(option.CallPayoff(), american);
        res.vega = (res_treeSigmaPlus.value - res_treeSigmaMinus.value) / (2*dsigma);
        
        // rho
        double drate = .045 / 1e4;
        TrinomialTree treeRatePlus(33, .24, 10./12., i, .045+drate, .02);
        auto res_treeRatePlus = treeRatePlus.TreePricer(option.CallPayoff(), american);
        TrinomialTree treeRateMinus(33, .24, 10./12., i, .045-drate, .02);
        auto res_treeRateMinus = treeRateMinus.TreePricer(option.CallPayoff(), american);
        res.rho = (res_treeRatePlus.value - res_treeRateMinus.value) / (2*drate);

        // result
        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(4) << res.gamma << std::setprecision(4) << '\t' << res.vega << '\t' << std::setprecision(4) << res.theta << '\t' << res.rho << std::setprecision(4) << std::endl;
    }

    std::cout << "\n\n\n" << "American Put" 
              << std::endl;
    for (std::size_t i = 10; i < 1281; i <<= 1) {
        TrinomialTree tree(33, .24, 10./12., i, .045, .02);
        auto res = tree.TreePricer(option.PutPayoff(), american);

        // vega
        double dsigma = .24 / 1e4;
        TrinomialTree treeSigmaPlus(33, .24+dsigma, 10./12., i, .045, .02);
        auto res_treeSigmaPlus = treeSigmaPlus.TreePricer(option.PutPayoff(), american);
        TrinomialTree treeSigmaMinus(33, .24-dsigma, 10./12., i, .045, .02);
        auto res_treeSigmaMinus = treeSigmaMinus.TreePricer(option.PutPayoff(), american);
        res.vega = (res_treeSigmaPlus.value - res_treeSigmaMinus.value) / (2*dsigma);
        
        // rho
        double drate = .045 / 1e4;
        TrinomialTree treeRatePlus(33, .24, 10./12., i, .045+drate, .02);
        auto res_treeRatePlus = treeRatePlus.TreePricer(option.PutPayoff(), american);
        TrinomialTree treeRateMinus(33, .24, 10./12., i, .045-drate, .02);
        auto res_treeRateMinus = treeRateMinus.TreePricer(option.PutPayoff(), american);
        res.rho = (res_treeRatePlus.value - res_treeRateMinus.value) / (2*drate);

        std::cout << res.value << '\t' << res.delta << '\t' << std::setprecision(4) << res.gamma << std::setprecision(4) << '\t' << res.vega << '\t' << std::setprecision(4) << res.theta << '\t' << res.rho << std::setprecision(4) << std::endl;
    }
}