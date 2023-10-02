#include "BinomialTree.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

double DeltaApproximation(double const V10, double const V11, double const S0_, double const u_, double const d_)
{
    return (V10 - V11) / (S0_ * (u_ - d_));
}

double GammaApproximation(double const V20, double const V21, double const V22, double const S0_, double const u_, double const d_)
{
    return ((V20 - V21) / (S0_ * (u_ * u_ - 1.)) - ((V21 - V22) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
}

double ThetaApproximation(double const V21, double const V00, double const dt_)
{
    return (V21 - V00) / (2. * dt_);
}

// constructor
BinomialTree::BinomialTree(double S0, double sigma, double T, std::size_t steps, double r, double q)
    : S0_(S0)
    , sigma_(sigma)
    , T_(T)
    , steps_(steps)
    , r_(r)
    , q_(q)
    , dt_(T / steps)
    , u_(std::exp(sigma * std::sqrt(dt_)))
    , d_(1. / u_)
    , r_disc_(std::exp(-r * dt_))
    , p_((std::exp((r - q) * dt_) - d_) / (u_ - d_))
    , disc_p(r_disc_ * p_)
    , disc_1p(r_disc_ * (1. - p_))
{
}

double BinomialTree::r() const
{
    return r_;
}
double BinomialTree::q() const
{
    return q_;
}

void BinomialTree::Backtrack_PathIndepedent(std::vector<double>& V) const
{
    for (auto Vit = V.begin(); Vit + 1 != V.end(); Vit++) {
        // Get value of the last node by taking the expectation and discounting
        // V = disc(p_u * V_u + p_d * V_d)
        *Vit = disc_p * (*Vit) + disc_1p * (*(Vit + 1));
    }

    // Shrink V since there are fewer nodes needed.
    V.pop_back();
}

void BinomialTree::Backtrack_American(std::vector<double>& V, std::array<std::deque<double>, 2>& S, size_t curr_step, const std::function<double(double)>& payoff_T) const
{
    // The current S[i] always has (curr_step + 1) elements.
    // The other has (curr_step) elements.
    size_t which_S = 0;
    if (S[1].size() == (curr_step + 1)) {
        which_S = 1;
    }

    // Find the payoff if we exercise the option at this step
    std::vector<double> exercise_payoff(S[which_S].size());
    std::transform(S[which_S].cbegin(), S[which_S].cend(), exercise_payoff.begin(), payoff_T);

    auto payoff_it = exercise_payoff.cbegin();
    for (auto Vit = V.begin(); Vit + 1 != V.end(); Vit++) {
        // Get value of the last node by taking the expectation and discounting
        // V = disc(p_u * V_u + p_d * V_d)
        *Vit = disc_p * (*Vit) + disc_1p * (*(Vit + 1));

        // If early exercise is more profitable, replace value with early exercise payoff.
        if ((*payoff_it) > (*Vit)) {
            *Vit = *payoff_it;
        }

        payoff_it++;
    }

    // Shrink V since there are fewer nodes needed.
    V.pop_back();

    // Also shrink the current S since we will have less steps in the next two iterations.
    if (S[which_S].size() >= 2) {
        S[which_S].pop_back();
        S[which_S].pop_front();
    }
}

TreeResult BinomialTree::TreePricer(const std::function<double(double, double)>& payoff, bool path_dependent) const
{

    if (path_dependent) {
        return this->TreePricer_PathDependent(payoff);
    } else {
        return this->TreePricer_PathIndependent(payoff);
    }
}

TreeResult BinomialTree::TreePricer_PathIndependent(const std::function<double(double, double)>& payoff) const
{

    // Generate asset price at maturity
    std::vector<double> S({ S0_ * std::pow(u_, static_cast<double>(steps_)) });
    double d2 = d_ * d_;
    for (std::size_t i = 0; i < steps_; i++) {
        S.emplace_back(S.back() * d2);
    }

    auto payoff_T = std::bind(payoff, std::placeholders::_1, T_);

    // Generate derivative value at maturity
    std::vector<double> V(S.size());
    std::transform(S.cbegin(), S.cend(), V.begin(), payoff_T);

    // Backtrack until t = 2 * dt (the resulting sub-tree is used for Greek calculation)
    std::size_t curr_step = steps_;
    while (curr_step > 2) {
        this->Backtrack_PathIndepedent(V);
        curr_step--;
    }

    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];

    // Backtrack
    this->Backtrack_PathIndepedent(V);

    // Store nodes of Step 2
    double V10 = V[0];
    double V11 = V[1];

    // Backtrack
    this->Backtrack_PathIndepedent(V);
    double V00 = V[0];

    double delta = (V10 - V11) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V21) / (S0_ * (u_ * u_ - 1.)) - ((V21 - V22) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V21 - V00) / (2. * dt_);

    return TreeResult({ V00, delta, gamma, theta });
}

TreeResult BinomialTree::TreePricer_PathDependent(const std::function<double(double, double)>& payoff) const
{

    // Generate asset prices
    // S[0]: Asset price at maturity
    // S[1]: Asset price dt before maturity
    std::array<std::deque<double>, 2> S;
    S[0].emplace_back(S0_ * std::pow(u_, static_cast<double>(steps_)));
    S[1].emplace_back(S0_ * std::pow(u_, static_cast<double>(steps_) - 1.));

    double d2 = d_ * d_;
    for (std::size_t i = 0; i < steps_ - 1; i++) {
        S[0].emplace_back(S[0].back() * d2);
        S[1].emplace_back(S[1].back() * d2);
    }
    S[0].emplace_back(S[0].back() * d2);
    // Now, S[0] has steps_ + 1 elements and S[1] has steps_ elements
    // The current S[i] always has (curr_step + 1) elements.

    std::size_t curr_time = T_;
    std::size_t curr_step = steps_;
    std::vector<double> V(S[0].size());
    auto payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    std::transform(S[0].cbegin(), S[0].cend(), V.begin(), payoff_T);
    S[0].pop_back();
    S[0].pop_front();
    curr_step--;
    curr_time -= dt_;

    while (curr_step >= 2) {
        payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
        this->Backtrack_American(V, S, curr_step, payoff_T);
        curr_step--;
        curr_time -= dt_;
    }

    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];

    // Backtrack
    payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    this->Backtrack_American(V, S, curr_step, payoff_T);
    curr_step--;
    curr_time -= dt_;

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];

    // Backtrack
    payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    this->Backtrack_American(V, S, curr_step, payoff_T);
    curr_step--;
    curr_time -= dt_;

    double V00 = V[0];
    double delta = DeltaApproximation(V10, V11, S0_, u_, d_);
    double gamma = GammaApproximation(V20, V21, V22, S0_, u_, d_);
    double theta = ThetaApproximation(V21, V00, dt_);

    return TreeResult({ V00, delta, gamma, theta });
}

TreeResult BinomialTree::TreePricer_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const
{

    TreeResult original_res = this->TreePricer_PathDependent(payoff);
    TreeResult euro = this->TreePricer_PathIndependent(payoff);

    double value = original_res.value - euro.value + option.Put();
    double delta = original_res.delta - euro.delta + option.DeltaPut();
    double gamma = original_res.gamma - euro.gamma + option.GammaPut();
    double theta = original_res.theta - euro.theta + option.ThetaPut();

    return TreeResult({ value, delta, gamma, theta });
}

TreeResult BinomialTree::AvgTreePricer(const std::function<double(double, double)>& payoff, bool path_dependent) const
{

    TreeResult res1 = this->TreePricer(payoff, path_dependent);

    BinomialTree longer_tree(S0_, sigma_, T_, steps_ + 1, r_, q_);

    TreeResult res2 = longer_tree.TreePricer(payoff, path_dependent);

    return TreeResult({ (res1.value + res2.value) / 2., (res1.delta + res2.delta) / 2., (res1.gamma + res2.gamma) / 2., (res1.theta + res2.theta) / 2. });
}

TreeResult BinomialTree::AvgTreePricer_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const
{

    TreeResult original_res = this->AvgTreePricer(payoff, true);
    TreeResult euro = this->AvgTreePricer(payoff, false);

    double value = original_res.value - euro.value + option.Put();
    double delta = original_res.delta - euro.delta + option.DeltaPut();
    double gamma = original_res.gamma - euro.gamma + option.GammaPut();
    double theta = original_res.theta - euro.theta + option.ThetaPut();

    return TreeResult({ value, delta, gamma, theta });
}

TreeResult BinomialTree::TreePricerBS(const std::function<double(double, double)>& payoff, bool path_dependent) const
{

    if (path_dependent) {
        return this->TreePricerBS_PathDependent(payoff);
    } else {
        return this->TreePricerBS_PathIndependent(payoff);
    }
}

TreeResult BinomialTree::TreePricerBS_PathIndependent(const std::function<double(double, double)>& payoff) const
{

    std::size_t curr_step = steps_;
    curr_step--; // We do not need the last step here

    // Generate asset price one step before maturity
    std::vector<double> S({ S0_ * std::pow(u_, static_cast<double>(curr_step)) });
    double d2 = d_ * d_;
    for (std::size_t i = 0; i < curr_step; i++) {
        S.emplace_back(S.back() * d2);
    }

    // Generate the BS pricer
    auto last_step = [&](double S_T) -> double {
        // PUT!!
        double K = payoff(0, 0);
        EuropeanOption option(T_ - dt_, S_T, K, T_, sigma_, r_, q_);

        return option.Put();
    };

    // Generate derivative value at the step just before maturity
    std::vector<double> V(S.size());
    std::transform(S.cbegin(), S.cend(), V.begin(), last_step);

    // Backtrack until t = 2 * dt (the resulting sub-tree is used for Greek calculation)
    while (curr_step > 2) {
        this->Backtrack_PathIndepedent(V);
        curr_step--;
    }

    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];

    // Backtrack
    this->Backtrack_PathIndepedent(V);

    // Store nodes of Step 2
    double V10 = V[0];
    double V11 = V[1];

    // Backtrack
    this->Backtrack_PathIndepedent(V);
    double V00 = V[0];

    double delta = (V10 - V11) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V21) / (S0_ * (u_ * u_ - 1.)) - ((V21 - V22) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V21 - V00) / (2. * dt_);

    return TreeResult({ V00, delta, gamma, theta });
}

TreeResult BinomialTree::TreePricerBS_PathDependent(const std::function<double(double, double)>& payoff) const
{

    std::size_t curr_step = steps_;
    std::size_t curr_time = T_;
    curr_step--; // We do not need the last step here
    curr_time -= dt_;

    // Generate asset prices
    // S[0]: Asset price at maturity
    // S[1]: Asset price dt before maturity
    std::array<std::deque<double>, 2> S;
    S[0].emplace_back(S0_ * std::pow(u_, static_cast<double>(curr_step)));
    S[1].emplace_back(S0_ * std::pow(u_, static_cast<double>(curr_step) - 1.));

    double d2 = d_ * d_;
    for (std::size_t i = 0; i < curr_step - 1; i++) {
        S[0].emplace_back(S[0].back() * d2);
        S[1].emplace_back(S[1].back() * d2);
    }
    S[0].emplace_back(S[0].back() * d2);
    // Now, S[0] has curr_step + 1 elements and S[1] has curr_step elements
    // The current S[i] always has (curr_step + 1) elements.

    // The BS pricer
    auto last_step = [&](double S_T) -> double {
        // PUT!!
        double K = payoff(0, 0);
        EuropeanOption option(T_ - dt_, S_T, K, T_, sigma_, r_, q_);

        return option.Put();
    };

    // Get BS prices
    std::vector<double> V(S[0].size());
    std::transform(S[0].cbegin(), S[0].cend(), V.begin(), last_step);

    // Get early exercise payoffs
    auto payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    std::vector<double> exercise_payoff(S[0].size());
    std::transform(S[0].cbegin(), S[0].cend(), exercise_payoff.begin(), payoff_T);

    // Compare
    auto payoff_it = exercise_payoff.cbegin();
    for (auto Vit = V.begin(); Vit != V.end(); Vit++) {
        if ((*Vit) < (*payoff_it)) {
            *Vit = *payoff_it;
        }
        payoff_it++;
    }

    S[0].pop_back();
    S[0].pop_front();
    curr_step--;
    curr_time -= dt_;

    while (curr_step >= 2) {
        payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
        this->Backtrack_American(V, S, curr_step, payoff_T);
        curr_step--;
        curr_time -= dt_;
    }

    // Store nodes of Step 2
    double V20 = V[0];
    double V21 = V[1];
    double V22 = V[2];

    // Backtrack
    payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    this->Backtrack_American(V, S, curr_step, payoff_T);
    curr_step--;
    curr_time -= dt_;

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];

    // Backtrack
    payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    this->Backtrack_American(V, S, curr_step, payoff_T);
    curr_step--;
    curr_time -= dt_;

    double V00 = V[0];
    double delta = (V10 - V11) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V21) / (S0_ * (u_ * u_ - 1.)) - ((V21 - V22) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V21 - V00) / (2. * dt_);

    return TreeResult({ V00, delta, gamma, theta });
}

TreeResult BinomialTree::TreePricerBSR(const std::function<double(double, double)>& payoff, bool path_dependent) const
{

    TreeResult res1 = this->TreePricerBS(payoff, path_dependent);

    BinomialTree shorter_tree(S0_, sigma_, T_, steps_ / 2, r_, q_);

    TreeResult res2 = shorter_tree.TreePricerBS(payoff, path_dependent);

    double value = res1.value * 2 - res2.value;
    double delta = res1.delta * 2 - res2.delta;
    double gamma = res1.gamma * 2 - res2.gamma;
    double theta = res1.theta * 2 - res2.theta;

    return TreeResult({ value, delta, gamma, theta });
}

TreeResult BinomialTree::TreePricerBS_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const
{
    TreeResult original_res = this->TreePricerBS(payoff, true);
    TreeResult euro = this->TreePricerBS(payoff, false);

    double value = original_res.value - euro.value + option.Put();
    double delta = original_res.delta - euro.delta + option.DeltaPut();
    double gamma = original_res.gamma - euro.gamma + option.GammaPut();
    double theta = original_res.theta - euro.theta + option.ThetaPut();

    return TreeResult({ value, delta, gamma, theta });
}
TreeResult BinomialTree::TreePricerBSR_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const
{
    TreeResult original_res = this->TreePricerBSR(payoff, true);
    TreeResult euro = this->TreePricerBSR(payoff, false);

    double value = original_res.value - euro.value + option.Put();
    double delta = original_res.delta - euro.delta + option.DeltaPut();
    double gamma = original_res.gamma - euro.gamma + option.GammaPut();
    double theta = original_res.theta - euro.theta + option.ThetaPut();

    return TreeResult({ value, delta, gamma, theta });
}
