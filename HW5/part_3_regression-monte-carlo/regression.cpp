#include "regression.h"
#include <tuple>

namespace linalg
{
	namespace regression
	{
		using namespace Eigen;

		arr evaluate_polynomial(const arr& coefficients, const arr& x)
		{
			// Evaluates a polynomial for an array of input values using Horner's rule.

			arr result = coefficients(coefficients.size() - 1) * arr::Ones(x.size()); // Start with the highest degree coefficient.
			
			// Loop through each coefficient using the nested equation specified by Horner's rule.
			for (int i = coefficients.size() - 2; i >= 0; --i)
			{
				result = result * x + coefficients(i);
			}
			return result;
		}

		arr evaluate_hermite_polynomial(const arr& coefficients, const arr& x)
		{
			// Evaluate a Hermite polynomial for an array of input values using the second order Taylor approximation.

			arr H_prev = arr::Ones(x.size());     // H_0.
			
			// Base case for single coefficient.
			if (coefficients.size() == 1)
			{
				return coefficients(0) * H_prev;
			}

			arr H_curr = x;                       // H_1.
			arr result = coefficients(0) * H_prev + coefficients(1) * H_curr; // Initialize result.

			// Recurrence relation (H_2, ...).
			for (int n = 2; n < coefficients.size(); ++n)
			{
				arr H_next = x * H_curr - (n-1) * H_prev;
				result += coefficients(n) * H_next;

				// Update for next iteration.
				H_prev = H_curr;
				H_curr = H_next;
			}

			return result;
		}

		arr evaluate_hermite_polynomial_standardized(const arr& coefficients, const arr& x, const double mu, const double sigma)
		{
			// Standardize input values, then evaluate a Hermite polynomial at these new input values.

			return evaluate_hermite_polynomial(coefficients, (x - mu) / sigma); // Delegate to the evaluate_hermite_polynomial() function.
		}
		
		arr2 vandermonde(const arr& x, const unsigned int degree)
		{
			// Generates the Vandermonde matrix for input values of X and a given degree of polynomial.

			arr2 matrix(x.size(), degree + 1);	// Matrix to populate. 1 row for each data point and 1 column for each degree + an additional for the constant.

			matrix.col(0) = 1.0;				// First column is all ones.
			for (unsigned int j = 1; j <= degree; ++j)
			{
				matrix.col(j) = x.pow(j);		// The rest of the columns are X^(column index).
			}
			return matrix;
		}

		arr2 hermite_vandermonde(const arr& x, const unsigned int degree)
		{
			// Generates the Vandermonde matrix for input values of X and a given degree of hermite polynomial.

			int n = x.size();							// Store length of X to avoid repeated calls.
			arr2 H = arr2::Zero(x.size(), degree + 1);  // Initialize Vandermonde matrix with zeros.

			// Base cases.
			H.col(0) = arr::Ones(n);	// First column is all ones.
			if (degree > 0)
			{
				H.col(1) = x;			// Second column is just x.
			}

			// Recursive definition of Hermite polynomials.
			for (unsigned int i = 2; i <= degree; i++)
			{
				H.col(i) = x * H.col(i-1) - (i-1) * H.col(i-2);
			}

			return H;
		}

		std::tuple<arr2, double, double> hermite_vandermonde_standardized(const arr& x, const unsigned int degree)
		{
			// Generates the standardized Vandermonde matrix for input values of X and a given degree of hermite polynomial. Returns a tuple
			// containing the Vandermonde matrix, mean, and standard deviation.
			
			// Calculate the mean and standard deviation of X.
			double mu = x.mean();
			double sigma = std::sqrt((x - mu).square().mean());

			// Standardize X.
			arr x_standardized = (x - mu) / sigma;

			// Use the standardized X to generate the Hermite Vandermonde matrix. Delegate to standard hermite vandermonde function.
			arr2 H = hermite_vandermonde(x_standardized, degree);

			return { H, mu, sigma }; // Make and return the tuple.
		}

		std::pair<arr,arr> fit_linear_regression(const arr2& X, const arr& y, const bool precondition)
		{
			// Given arrays of regressor and reponse observations, returns the OLS beta vector and array of fitted values. Optional flag for
			// standardizing the X values.

			Eigen::MatrixXd X_copy = X.matrix();	// Make a copy of the X array to manipulate.
			Eigen::VectorXd D_inv(X_copy.cols());	// Vector of scaling factors for normalizing the columns of X.

			// Preconditioning.
			if (precondition) // If preconditioning flag is set.
			{
				double norm_tolerance = 1e-10;	// Tolerance for norms.
				
				// Normalize the columns of X by their L2 norms and store the inverse of the norms in D_inv.
				for (int i = 0; i < X_copy.cols(); i++)
				{
					double norm = X_copy.col(i).norm();
					if (abs(norm) > norm_tolerance)
					{
						D_inv(i) = 1.0 / norm;
						X_copy.col(i) *= D_inv(i);
					}
					else
					{
						throw std::runtime_error("Encountered a column norm below the acceptable tolerance. Matrix may be ill-conditioned.");
					}
				}
			}
			
			// Compute the SVD of X_matrix.
			Eigen::BDCSVD<Eigen::MatrixXd> svd(X_copy, Eigen::ComputeThinU | Eigen::ComputeThinV);
			Eigen::MatrixXd U = svd.matrixU();
			Eigen::MatrixXd V = svd.matrixV();
			Eigen::VectorXd singularValues = svd.singularValues();

			// Initially set the pseudo-inverse of singular values to zero.
			Eigen::VectorXd singularValuesInv = Eigen::VectorXd::Zero(singularValues.size());

			// Define a threshold value below which singular values are treated as zero.
			double sval_tolerance = 1e-10 * singularValues.maxCoeff();

			// For singular values greater than the threshold, compute their inverse.
			for (int i = 0; i < singularValues.size(); i++)
			{
				if (singularValues[i] > sval_tolerance)
				{
					singularValuesInv[i] = 1.0 / singularValues[i];
				}
				else
				{
					throw std::runtime_error("Encountered a singular value below the acceptable tolerance. Matrix may be ill-conditioned.");
				}
			}

			// For numerical stability, set very small singular values to zero.
			for (int i = 0; i < singularValuesInv.size(); i++)
			{
				if (std::abs(singularValues[i]) < sval_tolerance)
				{
					singularValuesInv[i] = 0;
				}
			}

			// Compute U'y, which is shared between β and ŷ calculations.
			Eigen::VectorXd U_transpose_y = U.transpose() * y.matrix();

			// Calculate β = D^-1 * V * Σ^(-1) * U'y.
			Eigen::VectorXd beta = D_inv.asDiagonal() * V * singularValuesInv.asDiagonal() * U_transpose_y;

			// Calculate ŷ = U * U'y.
			Eigen::VectorXd y_fitted = U * U_transpose_y;

			return {beta.array(), y_fitted}; // Return beta and fitted values.
		}
	}
}