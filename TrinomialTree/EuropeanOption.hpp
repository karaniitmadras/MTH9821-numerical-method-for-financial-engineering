#ifndef EuropeanOption_hpp
#define EuropeanOption_hpp

#include <functional>

class EuropeanOption {
private:
    // t, S, K, T, sigma, r, q
    double t_; // Current time
    double S_; // Spot price
    double K_; // Strike price
    double T_; // Maturity
    double sigma_; // Volatility
    double r_; // Const interest rate
    double q_; // Dividend rate

    // Intermediate values for Black-Scholes pricing
    double d1_;
    double d2_;
    double Zd1_;
    double Zd2_;
    double Nd1_;
    double Nd2_;
    double q_disc_; // Discount by dividend rate: exp(-q(T-t))
    double r_disc_; // Discount by interest rate: exp(-r(T-t))

    double Phi(double t) const
    {
        return std::erfc(-t / std::sqrt(2.)) / 2.;
    }

    double Z(double t) const
    {
        return std::exp(-t * t / 2.) / std::sqrt(2. * M_PI);
    }

public:
    EuropeanOption(double t, double S, double K, double T, double sigma, double r, double q);
    ~EuropeanOption() = default;

    double Call() const
    {
        return S_ * q_disc_ * Nd1_ - K_ * r_disc_ * Nd2_;
    }

    double Put() const
    {
        return -S_ * q_disc_ * (1. - Nd1_) + K_ * r_disc_ * (1. - Nd2_);
    }

    std::function<double(double, double)> CallPayoff() const
    {
        return [&](double S, double t) -> double {
            return std::max(S - K_, 0.);
        };
    }

    std::function<double(double, double)> PutPayoff() const
    {
        return [&](double S, double t) -> double {
            return std::max(K_ - S, 0.);
        };
    }

    // Greeks
    double DeltaCall() const
    {
        return q_disc_ * Nd1_;
    }
    double DeltaPut() const
    {
        return -q_disc_ * (1. - Nd1_);
    }

    double GammaCall() const
    {
        return q_disc_ / (S_ * sigma_ * std::sqrt(T_ - t_)) * Zd1_;
    }
    double GammaPut() const
    {
        return this->GammaCall();
    }

    double ThetaCall() const
    {
        double res = -(S_ * sigma_ * q_disc_) / (2. * std::sqrt(T_ - t_)) * Zd1_;

        res += q_ * S_ * q_disc_ * Nd1_;

        res += -r_ * K_ * r_disc_ * Nd2_;

        return res;
    }

    double ThetaPut() const
    {
        double res = -(S_ * sigma_ * q_disc_) / (2. * std::sqrt(T_ - t_)) * Zd1_;

        res += -q_ * S_ * q_disc_ * (1 - Nd1_);

        res += r_ * K_ * r_disc_ * (1 - Nd2_);

        return res;
    }

    double VegaCall() const
    {
        return S_ * q_disc_ * std::sqrt(T_ - t_) * Zd1_;
    }

    double VegaPut() const
    {
        return this->VegaCall();
    }

    double RhoCall() const
    {
        return K_ * (T_ - t_) * r_disc_ * Nd2_;
    }

    double RhoPut() const
    {
        return -K_ * (T_ - t_) * r_disc_ * (1. - Nd2_);
    }
};

#endif /* EuropeanOption_hpp */
