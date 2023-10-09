#include "TrinomialTreeVolatilityUnknown.hpp"
#include "TrinomialTree.hpp"
#include "EuropeanOption.hpp"
#include <cmath>
#include <iostream>

TrinomialTreeVolatilityUnknown::TrinomialTreeVolatilityUnknown(double S, double K, double T, double r, double q, double V0)
    : S_(S)
    , K_(K)
    , T_(T)
    , r_(r)
    , q_(q)
    , V0_(V0)
{

    std::size_t steps = 10 >> 1;
    double vol_init = .25;
    double tol = std::pow(10, -6);

    EuropeanOption option(0., S, K, T, vol_init, r, q);

    double Tree_estimate = vol_init + 1.;
    double Tree_estimate_old = 0.;

    while (std::abs(Tree_estimate - Tree_estimate_old) >= tol) {
        steps <<= 1;
        TrinomialTree tree(S, vol_init, T, steps, r, q);
        Tree_estimate_old = Tree_estimate;
        Tree_estimate = tree.TreePricerBSR(option.PutPayoff(), TrinomialTreeVolatilityUnknown::american_).value;
        std::cout << steps << std::endl;
    }

    steps_ = steps;
}

double TrinomialTreeVolatilityUnknown::Put() const
{

    double tol = 0.0001;

    double sig_last = .05;
    double sig_new = 1.;

    // The payoff is the same regardless of implied volatility
    EuropeanOption option(0., S_, K_, T_, .3, r_, q_);

    TrinomialTree tree_oldest(S_, sig_last, T_, steps_, r_, q_);
    std::cout << tree_oldest.TreePricerBSR(option.PutPayoff(), TrinomialTreeVolatilityUnknown::american_).value << std::endl;
    TrinomialTree tree_older(S_, sig_new, T_, steps_, r_, q_);
    std::cout << tree_older.TreePricerBSR(option.PutPayoff(), TrinomialTreeVolatilityUnknown::american_).value << std::endl;

    while (std::abs(sig_new - sig_last) >= tol) {
        TrinomialTree new_tree(S_, sig_new, T_, steps_, r_, q_);
        TrinomialTree last_tree(S_, sig_last, T_, steps_, r_, q_);

        double v_new = new_tree.TreePricerBSR(option.PutPayoff(), TrinomialTreeVolatilityUnknown::american_).value;
        double v_last = last_tree.TreePricerBSR(option.PutPayoff(), TrinomialTreeVolatilityUnknown::american_).value;

        double sig_newest = sig_new - (v_new - V0_) * (sig_new - sig_last) / (v_new - v_last);

        sig_last = sig_new;
        sig_new = sig_newest;

        TrinomialTree tree(S_, sig_new, T_, steps_, r_, q_);
        std::cout << sig_new << '\t' << tree.TreePricerBSR(option.PutPayoff(), TrinomialTreeVolatilityUnknown::american_).value << std::endl;
    }

    return sig_new;
}