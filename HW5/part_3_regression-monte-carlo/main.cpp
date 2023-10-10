#include "regression.h"
#include "mc_regression.h"
#include "american_pricers.h"
#include "util.h"
#include <iostream>
#include <cassert> // For simple unit tests.

using namespace linalg::regression; // Contains the regression functions.

namespace montecarlo
{

void priceExperiments(
    double spot, double strike, double r, double sigma, 
    double maturity, double dt, int M, int N, int Ne, unsigned long initialSeed, 
    int degree, regression::MCRegression& mc_regression, const MonteCarloRegressionMethod& regression_method) {
    
    bool in_sample;
    if (initialSeed == 100) { in_sample = 1;} else { in_sample = 0;}

    LCG lcg(39373, 0, pow(2, 31) - 1, initialSeed);
    double cumulativePriceBackward = 0.0;
    double cumulativePriceForward = 0.0;
    std::vector<double> price_Backward_vec;
    std::vector<double> price_Forward_vec;

    for (int e = 0; e < Ne; ++e) {
        arr2 paths = generatePaths(spot, r, 0, sigma, dt, M, N, lcg); // Assume dividend rate is 0

        // regression_pricer_backward and regression_pricer_forward functions 
        double priceBackward = regression_pricer_backward(
            spot, paths, CallPutFlag::Put, strike, maturity,
            r, mc_regression, regression_method);
        double priceForward = regression_pricer_forward(
            spot, paths, CallPutFlag::Put, strike, maturity,
            r, mc_regression);

        cumulativePriceBackward += priceBackward;
        cumulativePriceForward += priceForward;
        price_Backward_vec.push_back(priceBackward);
        price_Forward_vec.push_back(priceForward);
    }

    double averagePriceBackward = cumulativePriceBackward / Ne;
    double averagePriceForward = cumulativePriceForward / Ne;

    std::cout << "Average Price (Backward): " << averagePriceBackward << std::endl;
    std::cout << "Average Price (Forward): " << averagePriceForward << std::endl;

    // 3.4 Sanity check for in-sample Longstaff_Schwartz setup
    std::cout << "========= Sanity check on in-sample Longstaff_Schwartz results =========" << std::endl;
    if ((regression_method == montecarlo::MonteCarloRegressionMethod::Longstaff_Schwartz) && (in_sample))
    {
        std::cout << "Backward Price = Forward Price" << (averagePriceBackward = averagePriceForward) << std::endl;
    }

    // To calculate standard error:
    double s_e_Backward = 0.0;
    double s_e_Forward = 0.0;
    for (int e = 0; e < Ne; ++e)
    {
        s_e_Backward += sqrt((pow(price_Backward_vec[e] - averagePriceBackward, 2) / Ne));
        s_e_Forward += sqrt((pow(price_Forward_vec[e] - averagePriceForward, 2) / Ne));
    }
    std::cout << "Standard Error (Backward): " << s_e_Backward << std::endl;
    std::cout << "Standard Error (Forward): " << s_e_Forward << std::endl;

    // Comparison with binomial tree method
    double binomialPrice = 4.079018801027898; // Given
    std::cout << "Difference from binomial (Backward): " << averagePriceBackward - binomialPrice << std::endl;
    std::cout << "Difference from binomial (Forward): " << averagePriceForward - binomialPrice << std::endl;
}

} // namespace montecarlo

int main()
{

    // Part 2 Polynomial Regression
    std::cout << "=========================================" << std::endl;
    std::cout << "========= Polynomial Regression =========" << std::endl;
    std::cout << "=========================================" << std::endl;

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
    std::cout << std::endl;


    // Part 3 Regression Monte Carlo
    std::cout << "========================================" << std::endl;
    std::cout << "=========Regression Monte Carlo=========" << std::endl;
    std::cout << "========================================" << std::endl;

    // Part 3.4 Calculations
    // A Bermuda Put Parameters 
    double spot = 40.5;
    double strike = 44;
    double r = 0.04; 
    sigma = 0.2; // declared in Part 2
    double maturity = 0.5;
    double dt = 0.5 / 6;
    double M = 6;
    double N = 10000;
    double Ne = 100;
    unsigned long initialSeed = 100;

    std::cout << std::endl;
    std::cout << "========= In sample pricing experiment =========" << std::endl;
    std::cout << std::endl;

    // TODO: @Zhuo follow the way below to compute and generate all tables
    std::cout << "========= Tsitsiklis_VanRoy  =========" << std::endl;
    std::cout << "========= Power Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {   
        bool precondition = false;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl; 
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::PolynomialMCRegression mc_regression = montecarlo::regression::PolynomialMCRegression(
            M, degree, precondition);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Tsitsiklis_VanRoy);
    }
    std::cout << std::endl;

    std::cout << "========= Tsitsiklis_VanRoy  =========" << std::endl;
    std::cout << "========= Hermite Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {   
        bool precondition = false;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl; 
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Tsitsiklis_VanRoy);
    }
    std::cout << std::endl;

    std::cout << "========= Tsitsiklis_VanRoy  =========" << std::endl;
    std::cout << "========= Hermite Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {   
        bool precondition = true;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl; 
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Tsitsiklis_VanRoy);
    }
    std::cout << std::endl;

    std::cout << "========= Longstaff_Schwartz  =========" << std::endl;
    std::cout << "========= Power Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = false;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Longstaff_Schwartz);
    }
    std::cout << std::endl;

    std::cout << "========= Longstaff_Schwartz  =========" << std::endl;
    std::cout << "========= Hermite Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = false;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Longstaff_Schwartz);
    }
    std::cout << std::endl;

    std::cout << "========= Longstaff_Schwartz  =========" << std::endl;
    std::cout << "========= Hermite Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = true;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Longstaff_Schwartz);
    }
    std::cout << std::endl;


    // For �out-of-sample� paths
    std::cout << std::endl;
    std::cout << "========= Out of sample pricing experiment =========" << std::endl;
    std::cout << std::endl;
    initialSeed = 200;

    std::cout << "========= Tsitsiklis_VanRoy  =========" << std::endl;
    std::cout << "========= Power Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = false;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::PolynomialMCRegression mc_regression = montecarlo::regression::PolynomialMCRegression(
            M, degree, precondition);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Tsitsiklis_VanRoy);
    }
    std::cout << std::endl;

    std::cout << "========= Tsitsiklis_VanRoy  =========" << std::endl;
    std::cout << "========= Hermite Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = false;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Tsitsiklis_VanRoy);
    }
    std::cout << std::endl;

    std::cout << "========= Tsitsiklis_VanRoy  =========" << std::endl;
    std::cout << "========= Hermite Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = true;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Tsitsiklis_VanRoy);
    }
    std::cout << std::endl;

    std::cout << "========= Longstaff_Schwartz  =========" << std::endl;
    std::cout << "========= Power Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = false;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Longstaff_Schwartz);
    }
    std::cout << std::endl;

    std::cout << "========= Longstaff_Schwartz  =========" << std::endl;
    std::cout << "========= Hermite Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = false;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Longstaff_Schwartz);
    }
    std::cout << std::endl;

    std::cout << "========= Longstaff_Schwartz  =========" << std::endl;
    std::cout << "========= Hermite Based  =========" << std::endl;
    for (int degree = 2; degree < 11; degree++)
    {
        bool precondition = true;
        bool standardization = false;
        std::cout << "Precondition: " << precondition << std::endl;
        std::cout << "Standardization: " << standardization << std::endl;
        montecarlo::regression::HermiteMCRegression mc_regression = montecarlo::regression::HermiteMCRegression(
            M, degree, precondition, standardization);
        montecarlo::priceExperiments(spot, strike, r, sigma, maturity, dt, M, N, Ne, initialSeed,
            degree, mc_regression, montecarlo::MonteCarloRegressionMethod::Longstaff_Schwartz);
    }
    std::cout << std::endl;

    return 0;
}