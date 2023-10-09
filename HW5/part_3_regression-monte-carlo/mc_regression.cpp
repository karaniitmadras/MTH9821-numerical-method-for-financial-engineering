#include "mc_regression.h"
#include "regression.h"
#include <limits>

namespace montecarlo {
	
    namespace regression {

        using namespace linalg::regression;

        MCRegression::MCRegression(const uint times):
            _times(times),
            _coeff(std::numeric_limits<double>::quiet_NaN()) {}

        MCRegression::~MCRegression() {}

        double MCRegression::fit_predict_at_0(const arr& y) {
            return _coeff = y.mean();
        }

        double MCRegression::predict_at_0() const {
            return _coeff;
        }

        PolynomialMCRegression::PolynomialMCRegression(const uint times, const uint degree, const bool precondition):
            MCRegression(times),
            degree(degree),
            precondition(precondition) {
            _coeffs.resize(times);
        }

        arr PolynomialMCRegression::fit_predict(const uint time, const arr& x, const arr& y) {
            if (time == 0) {
                throw std::out_of_range ("At time 0, fit_predict_at_0 should be used.");
            }
            arr result;
            std::tie(_coeffs.at(time), result) = fit_linear_regression(vandermonde(x, degree), y, precondition);
            return result;
        }


        arr PolynomialMCRegression::predict(const uint time, const arr& x) const {
            if (time == 0) {
                throw std::out_of_range ("At time 0, predict_at_0 should be used.");
            }
            return evaluate_polynomial(_coeffs.at(time), x);
        }

        /* Implement a constructor similar to that of PolynomialMCRegression */
        HermiteMCRegression::HermiteMCRegression(const uint times, const uint degree, const bool precondition,
                                                       const bool standardize):
            MCRegression(times),
            degree(degree),
            precondition(precondition),
            standardize(standardize) {
            _coeffs.resize(times);
            _means.resize(times);
            _stds.resize(times);
        }

        // Part 3.1
        /* Implement a Hermite version of PolynomialMCRegression::fit_predict */
        arr HermiteMCRegression::fit_predict(const uint time, const arr &x, const arr &y) {
            if (time == 0) {
                throw std::out_of_range ("At time 0, fit_predict_at_0 should be used.");
            }
            
            arr result;
            if (standardize)
            {
                arr2 matrix_H; // the standardized Vandermonde matrix for input values of X and a given degree of hermite polynomial
                std::tie(matrix_H, _means.at(time), _stds.at(time)) = hermite_vandermonde_standardized(x, degree);
                std::tie(_coeffs.at(time), result) = fit_linear_regression(matrix_H, y, precondition);
            }
            else
            {
                std::tie(_coeffs.at(time), result) = fit_linear_regression(hermite_vandermonde(x, degree), y, precondition);
            }

            return result;
        }

        // Part 3.1
        /* Implement a Hermite version of PolynomialMCRegression::predict */
        arr HermiteMCRegression::predict(const uint time, const arr& x) const {
            if (time == 0) {
                throw std::out_of_range ("At time 0, fit_predict_at_0 should be used.");
            }

            arr result;
            if (standardize)
            {
                result = evaluate_hermite_polynomial_standardized(_coeffs.at(time), x, _means.at(time), _stds.at(time));
            }
            else
            {
                result = evaluate_hermite_polynomial(_coeffs.at(time), x);
            }
            return result;
        }
    }
}

