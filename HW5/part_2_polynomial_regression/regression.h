#ifndef REGRESSION_H
#define REGRESSION_H

#include <tuple>   
#include <Dense> // Eigen package -- various functions for working with matrices. Source: https://eigen.tuxfamily.org/index.php?title=Main_Page.

namespace linalg
{
    namespace regression
    {
		typedef Eigen::ArrayXd arr;
		typedef Eigen::ArrayXXd arr2;
        arr evaluate_polynomial(const arr& coefficients, const arr& x);
        arr evaluate_hermite_polynomial(const arr& coefficients, const arr& x);
        arr evaluate_hermite_polynomial_standardized(const arr& coefficients, const arr& z, const double mu, const double sigma);
        arr2 vandermonde(const arr& x, const unsigned int degree);
        arr2 hermite_vandermonde(const arr& x, const unsigned int degree);
        std::tuple<arr2, double, double> hermite_vandermonde_standardized(const arr& x, const unsigned int degree);
        std::pair<arr,arr> fit_linear_regression(const arr2& X, const arr& y, const bool precondition) ;
    }
}
#endif // REGRESSION_H
