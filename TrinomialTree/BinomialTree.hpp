#ifndef BinomialTree_hpp
#define BinomialTree_hpp

#include "EuropeanOption.hpp"
#include <array>
#include <deque>
#include <functional>
#include <vector>

double DeltaApproximation(double const V10, double const V11, double const S0_, double const u_, double const d_);
double GammaApproximation(double const V20, double const V21, double const V22, double const S0_, double const u_, double const d_);
double ThetaApproximation(double const V21, double const V00, double const dt_);

struct TreeResult {
    double value;
    double delta;
    double gamma;
    double theta;
};

class BinomialTree {
    // Binomial tree pricer
    // Underlying asset follows a lognormal distribution

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
    double u_; // u = exp(sigma * sqrt(dt))
    double d_; // d = 1 / u
    double r_disc_; // exp(-r*dt)
    double p_; // (exp((r - q) * dt) - d) / (u - d)
    double disc_p; // disc * p
    double disc_1p; // disc * (1 - p)

    void Backtrack_PathIndepedent(std::vector<double>& vec) const;
    void Backtrack_American(std::vector<double>& V, std::array<std::deque<double>, 2>& S, size_t curr_step, const std::function<double(double)>& payoff_T) const;

    TreeResult TreePricer_PathIndependent(const std::function<double(double, double)>& payoff) const;
    TreeResult TreePricer_PathDependent(const std::function<double(double, double)>& payoff) const;

    TreeResult TreePricerBS_PathIndependent(const std::function<double(double, double)>& payoff) const;
    TreeResult TreePricerBS_PathDependent(const std::function<double(double, double)>& payoff) const;

public:
    BinomialTree(double S0, double sigma, double T, std::size_t steps, double r, double q);
    ~BinomialTree() = default;

    // Return asset information
    double r() const; // Interest rate
    double q() const; // Dividend rate

    // Tree
    TreeResult TreePricer(const std::function<double(double, double)>& payoff, bool path_dependent) const;

    TreeResult AvgTreePricer(const std::function<double(double, double)>& payoff, bool path_dependent) const;

    TreeResult TreePricerBS(const std::function<double(double, double)>& payoff, bool path_dependent) const;

    TreeResult TreePricerBSR(const std::function<double(double, double)>& payoff, bool path_dependent) const;

    // With variance reduction
    TreeResult TreePricer_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const;
    TreeResult AvgTreePricer_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const;
    TreeResult TreePricerBS_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const;
    TreeResult TreePricerBSR_PathDependent_compensated(const std::function<double(double, double)>& payoff, const EuropeanOption& option) const;
};

#endif /* BinomialTree_hpp */
