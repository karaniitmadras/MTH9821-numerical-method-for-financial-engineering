#include "regression.h"
#include <iostream>
#include <cassert> // For simple unit tests.

using namespace linalg::regression; // Contains the regression functions.

int main()
{
    // Tests the functionality of the functions defined in regression.h and implemented in regression.cpp.
    //
 
    //////////////////////////////
    // 1. Generate sample data. //
    //////////////////////////////
    arr x = arr::LinSpaced(10, 0, 9);       // Linearly spaced values from 0 to 9.
    arr y = 3*x + 2 + 0.5*arr::Random(10);  // Simple linear relation with some random noise.
    unsigned int degree = 2;                // Degree 2 polynomial.

    /////////////////////////////
    // 2. Evaluate & Validate. //
    /////////////////////////////

    // Evaluate polynomial.
    arr coefficients = arr(2);
    coefficients << 2, 3;
    arr polyEval = evaluate_polynomial(coefficients, x);
    std::cout << "Polynomial Evaluation: " << polyEval.transpose() << std::endl;

    // Evaluate Hermite polynomial.
    coefficients.resize(3);
    coefficients << 1, 2, 3;
    arr hermiteEval = evaluate_hermite_polynomial(coefficients, x);
    std::cout << "Hermite Polynomial Evaluation: " << hermiteEval.transpose() << std::endl;

    // Evaluate standardized Hermite polynomial.
    double mu = 5, sigma = 3;
    arr standardizedHermiteEval = evaluate_hermite_polynomial_standardized(coefficients, x, mu, sigma);
    std::cout << "Standardized Hermite Polynomial Evaluation: " << standardizedHermiteEval.transpose() << std::endl;

    // Generate Vandermonde matrix.
    arr2 vanderMatrix = vandermonde(x, degree);
    std::cout << "Vandermonde Matrix:\n" << vanderMatrix << std::endl;

    // Generate Hermite Vandermonde matrix.
    arr2 hermiteVanderMatrix = hermite_vandermonde(x, degree);
    std::cout << "Hermite Vandermonde Matrix:\n" << hermiteVanderMatrix << std::endl;

    // Generate standardized Hermite Vandermonde matrix.
    auto [hermiteVanderMatrixStandardized, mu_calculated, sigma_calculated] = hermite_vandermonde_standardized(x, degree);
    std::cout << "Standardized Hermite Vandermonde Matrix:\n" << hermiteVanderMatrixStandardized << std::endl;

    // Fit linear regression.
    bool precondition = true;
    auto [beta, y_fitted] = fit_linear_regression(vanderMatrix, y, precondition);
    std::cout << "Linear Regression Coefficients: " << beta.transpose() << std::endl;
    std::cout << "Fitted values: " << y_fitted.transpose() << std::endl;

    // Some simple assertions for validation.
    assert(polyEval.size() == x.size());
    assert(hermiteEval.size() == x.size());
    assert(standardizedHermiteEval.size() == x.size());
    assert(vanderMatrix.cols() == degree + 1);
    assert(hermiteVanderMatrix.cols() == degree + 1);
    assert(hermiteVanderMatrixStandardized.cols() == degree + 1);
    assert(beta.size() == degree + 1);
    assert(y_fitted.size() == y.size());

    std::cout << "All tests passed!" << std::endl;

    return 0;
}