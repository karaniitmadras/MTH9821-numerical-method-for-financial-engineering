#include "TrinomialTree.hpp"
#include <cmath>

TrinomialTree::TrinomialTree(double S0, double sigma, double T, std::size_t steps, double r, double q)
    : S0_(S0)
    , sigma_(sigma)
    , T_(T)
    , steps_(steps)
    , r_(r)
    , q_(q)
    , dt_(T / steps)
    , u_(std::exp(sigma * std::sqrt(2. * dt_)))
    , d_(1. / u_)
    , r_disc_(std::exp(-r * dt_))
    , pu_((std::exp((r - q) * dt_ / 2.) - std::exp(-sigma * std::sqrt(dt_ / 2.))) / (std::exp(sigma * std::sqrt(dt_ / 2.)) - std::exp(-sigma * std::sqrt(dt_ / 2.))))
    , pd_((std::exp(sigma * std::sqrt(dt_ / 2.)) - std::exp((r - q) * dt_ / 2.)) / (std::exp(sigma * std::sqrt(dt_ / 2.)) - std::exp(-sigma * std::sqrt(dt_ / 2.))))
{
    pu_ *= pu_;
    pd_ *= pd_;

    disc_pu = pu_ * r_disc_;
    disc_pd = pd_ * r_disc_;
    disc_pm = (1. - pu_ - pd_) * r_disc_;
}

void TrinomialTree::Backtrack_PathIndepedent(std::vector<double>& V) const
{
    std::vector<double> new_V;

    for (auto Vit = V.cbegin() + 1; Vit + 1 != V.cend(); Vit++) {
        // Get value of the last node by taking the expectation and discounting
        // V = disc(p_u * V_u + p_d * V_d)
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

    // Shrink the current S since we will have fewer outcomes in the next iteration.
    S.pop_front();
    S.pop_back();
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
    while (V.size() > 1) {
        this->Backtrack_PathIndepedent(V);
    }

    return TreeResult({ V[0], 0., 0., 0. });
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

    // Backtrack until the end
    while (V.size() > 1) {
        payoff_T = std::bind(payoff, std::placeholders::_1, curr_time);
        curr_time -= dt_;
        this->Backtrack_American(V, S, payoff_T);
    }

    return TreeResult({ V[0], 0., 0., 0. });
}
