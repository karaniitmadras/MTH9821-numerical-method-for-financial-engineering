#include "TrinomialTree.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

TrinomialTree::TrinomialTree(double S0, double sigma, double T, std::size_t steps, double r, double q)
    : S0_(S0)
    , sigma_(sigma)
    , T_(T)
    , steps_(steps)
    , r_(r)
    , q_(q)
    , dt_(T / steps)
    , u_(std::exp(sigma * std::sqrt(3. * dt_)))
    , d_(1. / u_)
    , r_disc_(std::exp(-r * dt_))
{
    double temp = (r_ - q_ - sigma_ * sigma_ * .5) * std::sqrt(dt_ / (12. * sigma_ * sigma_));
    pu_ = 1. / 6. + temp;
    pd_ = 1. / 6. - temp;

    disc_pu = pu_ * r_disc_;
    disc_pd = pd_ * r_disc_;
    disc_pm = 2. / 3. * r_disc_;
}

void TrinomialTree::Backtrack_PathIndepedent(std::vector<double>& V) const
{
    std::vector<double> new_V;

    for (auto Vit = V.cbegin() + 1; Vit + 1 != V.cend(); Vit++) {
        // Get value of the last node by taking the expectation and discounting
        // V = disc(p_u * V_u + p_m * V_m + p_d * V_d)
        new_V.push_back(disc_pu * (*(Vit - 1)) + disc_pm * (*Vit) + disc_pd * (*(Vit + 1)));
    }

    // Replace V
    V = std::move(new_V);
}

void TrinomialTree::Backtrack_American(std::vector<double>& V, std::deque<double>& S, const std::function<double(double)>& payoff_T) const
{

    Backtrack_PathIndepedent(V);

    // Find the payoff if we exercise the option at this step
    std::vector<double> exercise_payoff(S.size());
    std::transform(S.cbegin(), S.cend(), exercise_payoff.begin(), payoff_T);

    auto payoff_it = exercise_payoff.cbegin();
    for (auto Vit = V.begin(); Vit != V.end(); Vit++) {
        // If early exercise is more profitable, replace value with early exercise payoff.
        if ((*payoff_it) > (*Vit)) {
            *Vit = *payoff_it;
        }

        payoff_it++;
    }

    if (S.size() >= 2) {
        // Shrink the current S since we will have fewer outcomes in the next iteration.
        S.pop_front();
        S.pop_back();
    }
}

TreeResult TrinomialTree::TreePricer(const std::function<double(double, double)>& payoff, bool path_dependent) const
{

    if (path_dependent) {
        return this->TreePricer_PathDependent(payoff);
    } else {
        return this->TreePricer_PathIndependent(payoff);
    }
}

TreeResult TrinomialTree::TreePricer_PathIndependent(const std::function<double(double, double)>& payoff) const
{

    // Generate asset price at maturity
    std::vector<double> S({ S0_ * std::pow(u_, static_cast<double>(steps_)) });

    for (std::size_t i = 0; i < 2 * steps_; i++) {
        S.push_back(S.back() * d_);
    }

    auto payoff_T = std::bind(payoff, std::placeholders::_1, T_);

    // Generate derivative value at maturity
    std::vector<double> V(S.size());
    std::transform(S.cbegin(), S.cend(), V.begin(), payoff_T);

    // Backtrack until the end
    while (V.size() > 5) {
        this->Backtrack_PathIndepedent(V);
    }

    // Store nodes of Step 2
    double V20 = V[0];
    //    double V21 = V[1];
    double V22 = V[2];
    //    double V23 = V[3];
    double V24 = V[4];

    // Backtrack
    this->Backtrack_PathIndepedent(V);

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];

    // Backtrack
    this->Backtrack_PathIndepedent(V);
    double V00 = V[0];

    double delta = (V10 - V12) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V22) / (S0_ * (u_ * u_ - 1.)) - ((V22 - V24) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V11 - V00) / dt_;

    return TreeResult({ V00, delta, gamma, theta });
}

TreeResult TrinomialTree::TreePricer_PathDependent(const std::function<double(double, double)>& payoff) const
{
    // Generate asset price at maturity
    std::deque<double> S({ S0_ * std::pow(u_, static_cast<double>(steps_)) });

    for (std::size_t i = 0; i < 2 * steps_; i++) {
        S.push_back(S.back() * d_);
    }

    double curr_time = T_;
    std::function<double(double)> payoff_T = std::bind(payoff, std::placeholders::_1, T_);
    curr_time -= dt_;

    // Generate derivative value at maturity
    std::vector<double> V(S.size());
    std::transform(S.cbegin(), S.cend(), V.begin(), payoff_T);

    S.pop_front();
    S.pop_back();

    // Backtrack until the end
    while (V.size() > 5) {
        payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
        curr_time -= dt_;
        this->Backtrack_American(V, S, payoff_T);
    }

    // Store nodes of Step 2
    double V20 = V[0];
    //    double V21 = V[1];
    double V22 = V[2];
    //    double V23 = V[3];
    double V24 = V[4];

    // Backtrack
    payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    curr_time -= dt_;
    this->Backtrack_American(V, S, payoff_T);

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];

    // Backtrack
    payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    curr_time -= dt_;
    this->Backtrack_American(V, S, payoff_T);

    double V00 = V[0];

    double delta = (V10 - V12) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V22) / (S0_ * (u_ * u_ - 1.)) - ((V22 - V24) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V11 - V00) / dt_;

    return TreeResult({ V00, delta, gamma, theta });
}

TreeResult TrinomialTree::TreePricerBS(const std::function<double(double, double)>& payoff, bool path_dependent) const
{

    if (path_dependent) {
        return this->TreePricerBS_PathDependent(payoff);
    } else {
        return this->TreePricerBS_PathIndependent(payoff);
    }
}

TreeResult TrinomialTree::TreePricerBS_PathIndependent(const std::function<double(double, double)>& payoff) const
{

    // Generate asset price at T - dt
    std::vector<double> S({ S0_ * std::pow(u_, static_cast<double>(steps_ - 1)) });

    for (std::size_t i = 0; i < 2 * (steps_ - 1); i++) {
        S.push_back(S.back() * d_);
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

    // Backtrack until the end
    while (V.size() > 5) {
        this->Backtrack_PathIndepedent(V);
    }

    // Store nodes of Step 2
    double V20 = V[0];
    //    double V21 = V[1];
    double V22 = V[2];
    //    double V23 = V[3];
    double V24 = V[4];

    // Backtrack
    this->Backtrack_PathIndepedent(V);

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];

    // Backtrack
    this->Backtrack_PathIndepedent(V);
    double V00 = V[0];

    double delta = (V10 - V12) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V22) / (S0_ * (u_ * u_ - 1.)) - ((V22 - V24) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V11 - V00) / dt_;

    return TreeResult({ V00, delta, gamma, theta });
}

TreeResult TrinomialTree::TreePricerBS_PathDependent(const std::function<double(double, double)>& payoff) const
{

    double curr_time = T_;
    curr_time -= dt_;

    // Generate asset price at T - dt_;
    std::deque<double> S({ S0_ * std::pow(u_, static_cast<double>(steps_ - 1)) });

    for (std::size_t i = 0; i < 2 * (steps_ - 1); i++) {
        S.push_back(S.back() * d_);
    }

    // The BS pricer
    auto last_step = [&](double S_T) -> double {
        // PUT!!
        double K = payoff(0, 0);
        EuropeanOption option(T_ - dt_, S_T, K, T_, sigma_, r_, q_);

        return option.Put();
    };

    // Get BS prices
    std::vector<double> V(S.size());
    std::transform(S.cbegin(), S.cend(), V.begin(), last_step);

    // Get early exercise payoffs
    auto payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    std::vector<double> exercise_payoff(S.size());
    std::transform(S.cbegin(), S.cend(), exercise_payoff.begin(), payoff_T);

    // Compare
    auto payoff_it = exercise_payoff.cbegin();
    for (auto Vit = V.begin(); Vit != V.end(); Vit++) {
        if ((*Vit) < (*payoff_it)) {
            *Vit = *payoff_it;
        }
        payoff_it++;
    }

    S.pop_front();
    S.pop_back();
    curr_time -= dt_;

    // Backtrack until the end
    while (V.size() > 5) {
        payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
        curr_time -= dt_;
        this->Backtrack_American(V, S, payoff_T);
    }

    // Store nodes of Step 2
    double V20 = V[0];
    //    double V21 = V[1];
    double V22 = V[2];
    //    double V23 = V[3];
    double V24 = V[4];

    // Backtrack
    payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    curr_time -= dt_;
    this->Backtrack_American(V, S, payoff_T);

    // Store nodes of Step 1
    double V10 = V[0];
    double V11 = V[1];
    double V12 = V[2];

    // Backtrack
    payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
    curr_time -= dt_;
    this->Backtrack_American(V, S, payoff_T);

    double V00 = V[0];

    double delta = (V10 - V12) / (S0_ * (u_ - d_));
    double gamma = ((V20 - V22) / (S0_ * (u_ * u_ - 1.)) - ((V22 - V24) / (S0_ * (1. - d_ * d_)))) / (.5 * S0_ * (u_ * u_ - d_ * d_));
    double theta = (V11 - V00) / dt_;

    return TreeResult({ V00, delta, gamma, theta });
}

TreeResult TrinomialTree::TreePricerBSR(const std::function<double(double, double)>& payoff, bool path_dependent) const
{

    TreeResult res1 = this->TreePricerBS(payoff, path_dependent);

    TrinomialTree shorter_tree(S0_, sigma_, T_, steps_ / 2, r_, q_);

    TreeResult res2 = shorter_tree.TreePricerBS(payoff, path_dependent);

    double value = res1.value * 2 - res2.value;
    double delta = res1.delta * 2 - res2.delta;
    double gamma = res1.gamma * 2 - res2.gamma;
    double theta = res1.theta * 2 - res2.theta;

    return TreeResult({ value, delta, gamma, theta });
}

TreeResult TrinomialTree::TreePricer_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const
{

    TreeResult original_res = this->TreePricer_PathDependent(payoff);
    TreeResult euro = this->TreePricer_PathIndependent(payoff);

    double value = original_res.value - euro.value + option.Put();
    double delta = original_res.delta - euro.delta + option.DeltaPut();
    double gamma = original_res.gamma - euro.gamma + option.GammaPut();
    double theta = original_res.theta - euro.theta + option.ThetaPut();

    return TreeResult({ value, delta, gamma, theta });
}

TreeResult TrinomialTree::TreePricerBS_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const
{
    TreeResult original_res = this->TreePricerBS(payoff, true);
    TreeResult euro = this->TreePricerBS(payoff, false);

    double value = original_res.value - euro.value + option.Put();
    double delta = original_res.delta - euro.delta + option.DeltaPut();
    double gamma = original_res.gamma - euro.gamma + option.GammaPut();
    double theta = original_res.theta - euro.theta + option.ThetaPut();

    return TreeResult({ value, delta, gamma, theta });
}

TreeResult TrinomialTree::TreePricerBSR_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const
{
    TreeResult original_res = this->TreePricerBSR(payoff, true);
    TreeResult euro = this->TreePricerBSR(payoff, false);

    double value = original_res.value - euro.value + option.Put();
    double delta = original_res.delta - euro.delta + option.DeltaPut();
    double gamma = original_res.gamma - euro.gamma + option.GammaPut();
    double theta = original_res.theta - euro.theta + option.ThetaPut();

    return TreeResult({ value, delta, gamma, theta });
}
