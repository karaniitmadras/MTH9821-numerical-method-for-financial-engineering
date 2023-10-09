#include "american_pricers.h"


namespace montecarlo {
    double regression_pricer_backward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
                                      const double interest_rate, regression::MCRegression& mc_regression, const MonteCarloRegressionMethod& method)
    {
        return 0.0;
    }
    double regression_pricer_forward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
                                      const double interest_rate, const regression::MCRegression& mc_regression)
    {
        return 0.0;
    }
}
