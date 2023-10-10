#ifndef MC_AMERICAN_PRICERS_H
#define MC_AMERICAN_PRICERS_H

#include "mc_regression.h"
#include <vector>
#include <Dense>

namespace montecarlo {
	typedef Eigen::ArrayXd arr;
	typedef Eigen::ArrayXXd arr2;

    enum MonteCarloRegressionMethod {Longstaff_Schwartz, Tsitsiklis_VanRoy};
    enum CallPutFlag { Call = 1, Put = -1};

    /*double regression_pricer_backward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
                                      const double interest_rate, regression::MCRegression& mc_regression, const MonteCarloRegressionMethod& method);
    double regression_pricer_forward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
                                      const double interest_rate, const regression::MCRegression& mc_regression);*/

    //Part 3.2 Backward Pricer
    double regression_pricer_backward(
        const double spot, const arr2& paths, const int w, const double strike, const double maturity,
        const double interest_rate, regression::MCRegression& mc_regression, const MonteCarloRegressionMethod& method) {

        int N = paths.rows();  // number of paths
        int M = paths.cols();  // number of timesteps

        double dt = maturity / M;

        arr E(N);  // immediate Exercise value
        arr C(N);  // estimated Continuation value
        arr P(N);  // actual Payoff value

        P.fill(0);

        for (int j = M - 1; j >= 0; j--) {

            E = w * (paths.col(j) - strike).max(arr::Zero(N));

            if (j == M - 1) {
                P = E;
            }
            else {

                if (j == 0) {

                    C = mc_regression.fit_predict_at_0(P * exp(-interest_rate * dt));
                }
                else
                {
                    C = mc_regression.fit_predict(j, paths.col(j), P * exp(-interest_rate * dt));
                }


                if (method == Longstaff_Schwartz) {
                    for (int i = 0; i < N; i++) {
                        if (E[i] >= C[i]) {
                            P[i] = E[i];
                        }
                        else {
                            P[i] *= exp(-interest_rate * dt);
                        }
                    }
                }
                else {  // Tsitsiklis_VanRoy
                    P = E.max(C);
                }
            }
        }

        return P.mean();
    }

    // Part 3.3 Forward Pricer
    double regression_pricer_forward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
        const double interest_rate, regression::MCRegression& mc_regression) {

        int N = paths.rows();  // number of paths
        int M = paths.cols();  // number of timesteps

        double dt = maturity / M;

        arr E(N);  // immediate Exercise value
        arr P(N);  // actual Payoff value

        P.fill(0);

        for (int j = 0; j < M; j++) {
            E = w * (strike - paths.col(j)).max(arr::Zero(N));

            if (j == 0) {

                //TODO: error: no viable conversion from 'double' to 'montecarlo::arr'
                double C0 = mc_regression.predict_at_0();

                if (E[0] > C0) {
                    P[0] = E[0] * exp(-interest_rate * j * dt);
                }
            }

            else {
                arr C = mc_regression.predict(j, paths.col(j));
                for (int i = 0; i < N; i++) {
                    if (E[i] > C[i] && P[i] == 0) {
                        P[i] = E[i] * exp(-interest_rate * j * dt);
                    }
                }
            }
        }

        return P.mean();
    }
}

#endif // MC_AMERICAN_PRICERS_H
