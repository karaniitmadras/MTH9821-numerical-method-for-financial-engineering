#include "BinomialTreeVolatilityUnknown.hpp"
#include "TrinomialTree.hpp"
#include <cmath>
#include <iostream>

BinomialTreeVolatilityUnknown::BinomialTreeVolatilityUnknown(double S, double K, double T, double r, double q, double V0)
    : S_(S)
    , K_(K)
    , T_(T)
    , r_(r)
    , q_(q)
    , V0_(V0)
{
    std::size_t steps = 10 / 2;

    double tol = std::pow(10, -6);

    EuropeanOption option(0., S, K, T, .3, r, q);
    double BS = option.Put();
    double BTree_estimate = BS + 1.;
    double BTree_estimate_old = 0.;

    while (std::abs(BTree_estimate - BTree_estimate_old) >= tol) {
        steps *= 2;
        BinomialTree tree(S, .3, T, steps, r, q);
        BTree_estimate_old = BTree_estimate;
        BTree_estimate = tree.TreePricerBSR(option.PutPayoff(), false).value;
        std::cout << steps << std::endl;
    }

    steps_ = steps;
}

double BinomialTreeVolatilityUnknown::PutBinomialTreeVolatilityUnknown() const
{
    double tol = 1e-6;

    double sig_last = .05;
    double sig_new = 1.;

    EuropeanOption option(0., S_, K_, T_, .3, r_, q_); // The payoff is the same regardless of implied volatility

    BinomialTree tree_oldest(S_, sig_last, T_, steps_, r_, q_);
    std::cout << tree_oldest.TreePricerBSR(option.PutPayoff(), false).value << std::endl;
    BinomialTree tree_older(S_, sig_new, T_, steps_, r_, q_);
    std::cout << tree_older.TreePricerBSR(option.PutPayoff(), false).value << std::endl;

    while (std::abs(sig_new - sig_last) >= tol) {
        BinomialTree new_tree(S_, sig_new, T_, steps_, r_, q_);
        BinomialTree last_tree(S_, sig_last, T_, steps_, r_, q_);

        double v_new = new_tree.TreePricerBSR(option.PutPayoff(), false).value;
        double v_last = last_tree.TreePricerBSR(option.PutPayoff(), false).value;

        double sig_newest = sig_new - (v_new - V0_) * (sig_new - sig_last) / (v_new - v_last);

        sig_last = sig_new;
        sig_new = sig_newest;

        BinomialTree tree(S_, sig_new, T_, steps_, r_, q_);
        std::cout << tree.TreePricerBSR(option.PutPayoff(), false).value << std::endl;
    }

    return sig_new;
}
