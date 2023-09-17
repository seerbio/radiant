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

#include <iostream>
#include "tensorflow/c/c_api.h"



#include <iostream>
#include <Eigen/Dense>


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

    ERR_INIT

    TF_Status* status = TF_NewStatus();
    TF_Graph* graph = TF_NewGraph();
    TF_SessionOptions* session_options = TF_NewSessionOptions();
    TF_Buffer* run_options = nullptr;
    TF_Session* session = TF_NewSession(graph, session_options, status);

    if (TF_GetCode(status) != TF_OK) {
        std::cerr << "Error initializing TensorFlow: " << TF_Message(status) << std::endl;
        ERR_RETURN
    }

    // Define the computation graph
    // ... (create input, hidden layers, batch normalization, output layers)

    // Define loss function and optimizer
    // ... (create loss node, optimizer, etc.)

    // Initialize input data and labels
    // ... (prepare input data and labels)

    // Training loop
    int num_iterations = 1000;
    for (int i = 0; i < num_iterations; ++i) {
        // Run a training step
        // ... (run optimizer, update weights, compute loss)
    }

    // Test the model
    // ... (run inference on test data)

    // Clean up resources
    TF_CloseSession(session, status);
    TF_DeleteSession(session, status);
    TF_DeleteSessionOptions(session_options);
    TF_DeleteGraph(graph);
    TF_DeleteStatus(status);


    ERR_RETURN
}
