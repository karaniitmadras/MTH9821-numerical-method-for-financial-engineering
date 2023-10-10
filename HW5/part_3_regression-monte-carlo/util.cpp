#include "util.h"


namespace montecarlo
{

// Utility Functions
LCG::LCG(unsigned long long a, unsigned long long c, unsigned long long k, unsigned long long seed)
    : _a(a), _c(c), _k(k), _x(seed) {}

double LCG::next() {
    _x = (_a * _x + _c) % _k;
    return static_cast<double>(_x) / _k;
}


double beasleySpringerMoro(double u) {

    static const double a[] = { 2.50662823884,
                               -18.61500062529,
                               41.39119773534,
                               -25.44106049637 };

    static const double b[] = { -8.47351093090,
                               23.08336743743,
                               -21.06224101826,
                               3.13082909833 };

    static const double c[] = { 0.3374754822726147,
                               0.9761690190917186,
                               0.1607979714918209,
                               0.0276438810333863,
                               0.0038405729373609,
                               0.0003951896511919,
                               0.0000321767881768,
                               0.0000002888167364,
                               0.0000003960315187 };

    if (u < 0.5) {
        u = 1.0 - u;
    }

    double y = u - 0.5;
    double r, x;
    if (std::abs(y) < 0.42) {
        r = y * y;
        x = y * (((a[3] * r + a[2]) * r + a[1]) * r + a[0]) /
            ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
    }
    else {
        r = u;
        if (y > 0) {
            r = 1.0 - u;
        }
        r = std::log(-std::log(r));
        x = c[0] + r * (c[1] + r * (c[2] + r * (c[3] + r * (c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))));

        if (y < 0) {
            x = -x;
        }
    }

    return x;
}

std::vector<double> generateStandardNormalSamples(int N) {

    LCG lcg(39373, 0, pow(2, 31) - 1, 100);

    std::vector<double> samples(N);

    for (int i = 0; i < N; ++i) {
        double u = lcg.next();
        samples[i] = beasleySpringerMoro(u);
    }

    return samples;
}

arr2 generatePaths(const double spot, const double r, const double q, const double sigma, const double dt, const int M, const int N, LCG& lcg) {
    
    arr2 paths(N, M);
    
    paths.col(0).fill(spot);

    for (int k = 0; k < N; ++k) {
        double S = spot;
        for (int j = 1; j < M; ++j) {
            double z = beasleySpringerMoro(lcg.next());
            S *= exp((r - q - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * z);
            paths(k, j) = S;
        }
    }

    return paths;
}

} // namespace montecarlo