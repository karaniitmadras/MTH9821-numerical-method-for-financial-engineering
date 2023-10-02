#ifndef TrinomialTree_hpp
#define TrinomialTree_hpp

#include "BinomialTree.hpp"
#include <algorithm>
#include <deque>
#include <vector>

class TrinomialTree {
    // Trinomial tree pricer
    // Underlying asset also follows a lognormal distribution

private:
    // Constructor takes in those values
    double S0_; // Spot price at t = 0
    double sigma_; // Volatility
    double T_; // Maturity
    double steps_; // Steps
    double r_; // Interest rate
    double q_; // Dividend rate

    // Helper values
    double dt_; // delta t = T / steps
    double u_; // u = exp(sigma * sqrt(2 * dt))
    double d_; // d = 1 / u
    double r_disc_; // exp(-r*dt)
    double pu_; // ((exp((r - q) * dt / 2) - exp(-sigma * sqrt(dt / 2))) / (exp(sigma * sqrt(dt / 2)) - exp(-sigma * sqrt(dt / 2)))) ^ 2
    double pd_; // ((exp(sigma * sqrt(dt / 2)) - exp((r - q) * dt / 2)) / (exp(sigma * sqrt(dt / 2)) - exp(-sigma * sqrt(dt / 2)))) ^ 2
    double disc_pu; // disc * pu
    double disc_pd; // disc * pd
    double disc_pm; // disc * (1 - pu - pd)

    void Backtrack_PathIndepedent(std::vector<double>& V) const;
    void Backtrack_American(std::vector<double>& V, std::deque<double>& S, const std::function<double(double)>& payoff_T) const;

    TreeResult TreePricer_PathIndependent(const std::function<double(double, double)>& payoff) const;
    TreeResult TreePricer_PathDependent(const std::function<double(double, double)>& payoff) const;

public:
    TrinomialTree(double S0, double sigma, double T, std::size_t steps, double r, double q);
    ~TrinomialTree() = default;

    TreeResult TreePricer(const std::function<double(double, double)>& payoff, bool path_dependent) const;
};

#endif /* TrinomialTree_hpp */
