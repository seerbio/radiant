//
// Created by anichols on 9/9/23.
//

#include "ConvexOptimizer.h"

#include <iostream>
#include <vector>
#include <stdio.h>
#include "ErrorUtils.h"
#include "EigenUtils.h"

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>


void optimize() {
    //    // Calculate the critical points
//    double x_crit = -d / (2 * a);
//    double y_crit = -e / (2 * b);
//    double z_crit = -f / (2 * c);
//
//    // Evaluate the function at the critical point
//    double max_value = a * x_crit * x_crit + b * y_crit * y_crit + c * z_crit * z_crit
//                       + d * x_crit + e * y_crit + f * z_crit + g;
//
//    // Output results
//    std::cout << "Critical Point (x, y, z): " << x_crit << ", " << y_crit << ", " << z_crit << std::endl;
//    std::cout << "Maximum Value of y: " << max_value << std::endl;
//
//
//    qDebug() << evalit(1, 2.2, 4, a, b, c, d, e, f, g);
}


#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

void levenbergMarquardt(
        std::function<void(const VectorXd&, VectorXd&)> residual,
        MatrixXd& x,
        const VectorXd& y,
        double lambda = 1e-4,
        int maxIterations = 100)
{
    const double epsilon = std::numeric_limits<double>::epsilon();

    int n = x.size();
    int m = y.size();

    for (int iteration = 0; iteration < maxIterations; ++iteration)
    {
        VectorXd f(m);
        VectorXd J(m, n);

        // Evaluate residuals and Jacobian
        residual(x, f);

        for (int j = 0; j < n; ++j)
        {
            MatrixXd x_plus_delta = x;
            x_plus_delta(j) += epsilon;
            residual(x_plus_delta, f);
            J.col(j) = (f - y) / epsilon;
        }

        // Calculate update step using Levenberg-Marquardt formula
        MatrixXd A = J.transpose() * J + lambda * MatrixXd::Identity(n, n);
        VectorXd b = J.transpose() * (f - y);
        VectorXd dx = A.ldlt().solve(b);

        // Apply the update
        VectorXd x_new = x - dx;

        std::cout << x_new << std::endl;

        // Evaluate residuals at the updated point
        VectorXd f_new(m);
        residual(x_new, f_new);

        // Compute error improvement
        double error_old = (f - y).squaredNorm();
        double error_new = (f_new - y).squaredNorm();

        if (error_new < error_old)
        {
            lambda /= 10.0;
            x = x_new;
        }
        else
        {
            lambda *= 10.0;
        }

        // Check for convergence
        if (dx.norm() < epsilon)
        {
            std::cout << "Converged after " << iteration << " iterations." << std::endl;
            break;
        }
    }
}

// Example usage:
void exampleResidual(const VectorXd& x, VectorXd& f)
{
    // Example: f(x) = x^2 - 4
    f(0) = x(0)*x(0) - 4.0;
}


Err ConvexOptimizer::test() {

    Eigen::VectorXd a(10);
    a << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    Eigen::VectorXd b(10);
    b << 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;

    Eigen::VectorXd c(10);
    c << 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;

    Eigen::VectorXd y(10);
    y << 10, 20, 30, 40, 50, 40, 30, 20, 10, 0;


    Eigen::MatrixX<double> A(a.size(), 3);
    A.col(0) = a;
    A.col(1) = b;
    A.col(2) = c;

    std::cout << A << std::endl;

    int N = y.size();



    levenbergMarquardt(exampleResidual, A, y);

    std::cout << "Optimized x: " << A << std::endl;


    return eNoError;
}
