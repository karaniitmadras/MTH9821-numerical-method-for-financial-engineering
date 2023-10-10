#ifndef UTIL_H
#define UTIL_H

#include "mc_regression.h"
#include "american_pricers.h"
#include <vector>

namespace montecarlo
{
    // Utility Functions
    class LCG {
    private:
        unsigned long long _a, _c, _k, _x;

    public:
        LCG(unsigned long long a, unsigned long long c, unsigned long long k, unsigned long long seed);

        double next();
    };


    double beasleySpringerMoro(double u);

    std::vector<double> generateStandardNormalSamples(int N);

    arr2 generatePaths(
        const double spot, const double r, const double q, const double sigma, 
        const double dt, const int M, const int N, LCG& lcg);
}

# endif