#define _USE_MATH_DEFINES

#include "EuropeanOption.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <math.h>

EuropeanOption::EuropeanOption(double t, double S, double K, double T, double sigma, double r, double q)
    : t_(t)
    , S_(S)
    , K_(K)
    , T_(T)
    , sigma_(sigma)
    , r_(r)
    , q_(q)
{
    d1_ = (log(S / K) + (r - q + sigma * sigma / 2.) * (T - t)) / (sigma * std::sqrt(T - t));
    d2_ = d1_ - sigma * std::sqrt(T - t);

    Zd1_ = this->Z(d1_);
    Zd2_ = this->Z(d2_);
    Nd1_ = this->Phi(d1_);
    Nd2_ = this->Phi(d2_);

    q_disc_ = std::exp(-q * (T - t));
    r_disc_ = std::exp(-r * (T - t));
}