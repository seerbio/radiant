//
// Created by anichols on 9/9/23.
//

#include "ConvexOptimizer.h"

#include <iostream>
#include <vector>
#include <stdio.h>

#include <boost/math/distributions.hpp>

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
//#include <boost/numeric/ublas/lu_solve.hpp>
//#include <boost/math/tools/optimization.hpp>

namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;

#include <Eigen/Dense>


//std::vector<double> calculateMainEffects(const std::vector<std::vector<double>>& data, const std::vector<double>& responses) {
//    int numFactors = data[0].size();
//    int numRuns = data.size();
//
//    std::vector<double> mainEffects(numFactors, 0.0);
//
//    for (int factor = 0; factor < numFactors; ++factor) {
//        double sum = 0.0;
//        for (int run = 0; run < numRuns; ++run) {
//            if (run < numRuns / 2) {
//                sum += data[run][factor] * responses[run];
//            } else {
//                sum -= data[run][factor] * responses[run];
//            }
//        }
//        mainEffects[factor] = sum / numRuns;
//    }
//
//    return mainEffects;
//}
//
//
//std::vector<int> findOptimalSettings(const std::vector<double>& mainEffects) {
//    std::vector<int> optimalSettings;
//    for (double effect : mainEffects) {
//        // Assuming higher main effect coefficients indicate higher response
//        if (effect > 0) {
//            optimalSettings.push_back(1); // High level
//        } else {
//            optimalSettings.push_back(-1); // Low level
//        }
//    }
//    return optimalSettings;
//}
//
//
//double objectiveFunction(const std::vector<double>& settings) {
//    // This is a hypothetical objective function. Replace it with your actual function.
//    double response = -1.0 * (pow(settings[0] - 2.0, 2) + pow(settings[1] - 3.0, 2));
//    return response;
//}
//
//
//std::vector<double> findOptimalSettings() {
//    // Define the number of continuous variables
//    int numVariables = 2; // Adjust this according to your scenario
//
//    // Define initial guesses for the settings (can be random or based on knowledge)
//    std::vector<double> initialSettings(numVariables, 0.0);
//
//    // Define optimization parameters (e.g., tolerance, maximum iterations)
//    const double tolerance = 1e-6;
//    const int maxIterations = 1000;
//
//    // Perform optimization (using a simple gradient ascent approach)
//    std::vector<double> currentSettings = initialSettings;
//    double currentObjective = objectiveFunction(currentSettings);
//
//    for (int iteration = 0; iteration < maxIterations; ++iteration) {
//        std::vector<double> gradients(numVariables, 0.0);
//
//        // Calculate gradients numerically (you may use more sophisticated methods)
//        const double epsilon = 1e-5;
//        for (int i = 0; i < numVariables; ++i) {
//            std::vector<double> perturbedSettings = currentSettings;
//            perturbedSettings[i] += epsilon;
//            double perturbedObjective = objectiveFunction(perturbedSettings);
//            gradients[i] = (perturbedObjective - currentObjective) / epsilon;
//        }
//
//        // Update settings based on gradients
//        for (int i = 0; i < numVariables; ++i) {
//            currentSettings[i] += 0.1 * gradients[i]; // Adjust the step size as needed
//        }
//
//        double newObjective = objectiveFunction(currentSettings);
//
//        // Check for convergence
//        if (abs(newObjective - currentObjective) < tolerance) {
//            break;
//        }
//
//        currentObjective = newObjective;
//    }
//
//    return currentSettings;
//}
//
//
//Err ConvexOptimizer::test() {
//
//    ERR_INIT
//
//    // Define the factors and their levels
//    int numFactors = 3;
//    int numLevels = 2;
//
//    // Create a matrix to represent the experimental design
//    std::vector<std::vector<double>> experiment = {
//            {1, 1, 1},
//            {1, -1, 1},
//            {-1, 1, 1},
//            {-1, -1, 1}
//    };
//
//    std::vector<double> responses = {10.5, 12.3, 8.9, 9.7};
//
//    // Perform the experiment and collect responses (not shown in this example)
//    // You would typically perform your actual experiments and collect data here.
//
//    // Calculate main effects
//    std::vector<double> mainEffects = calculateMainEffects(experiment, responses);
//
//    // Display the main effects
//    std::cout << "Main Effects:" << std::endl;
//    for (int i = 0; i < numFactors; ++i) {
//        std::cout << " Factor " << i + 1 << ": " << mainEffects[i] << std::endl;
//    }
//
//
//    // Main effects obtained from analysis
////    std::vector<double> mainEffects = {1.2, -0.8, 0.5}; // Example values
//
//    // Find optimal factor settings
////    std::vector<int> optimalSettings = findOptimalSettings(mainEffects);
////
////    // Display the optimal settings
////    std::cout << "Optimal Factor Settings:" << std::endl;
////    for (size_t i = 0; i < optimalSettings.size(); ++i) {
////        std::cout << "Factor " << i + 1 << ": " << (optimalSettings[i] > 0 ? "High" : "Low") << " Level" << std::endl;
////    }
//
//    // Find the optimal settings
//    std::vector<double> optimalSettings = findOptimalSettings();
//
//    // Display the optimal settings and the corresponding response
//    std::cout << "Optimal Settings:" << std::endl;
//    for (size_t i = 0; i < optimalSettings.size(); ++i) {
//        std::cout << "Variable " << i + 1 << ": " << fixed << optimalSettings[i] << std::endl;
//    }
//
//    // Calculate and display the response at the optimal settings
//    double optimalResponse = objectiveFunction(optimalSettings);
//    std::cout << "Optimal Response: " << fixed << optimalResponse << std::endl;
//
//
//
//
//
//    ERR_RETURN
//}


Err ConvexOptimizer::test() {

    ERR_INIT

    // Define the data
    Eigen::MatrixXd X(3, 3);
    X <<
         2, 2.2, 1,
         1, 2.2, 2,
         1, 2.2, 4;

    std::cout << X << std::endl;

    Eigen::VectorXd y(3);
    y << 20, 20, 15;

    Eigen::MatrixXd X_poly(3, 6);
    X_poly << X, X.array().square();

    Eigen::VectorXd coefficients = X.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);

    // Constrain the coefficients
    for (int i = 0; i < coefficients.size(); ++i) {
        coefficients(i) = std::max(1.0, std::min(7.0, coefficients(i)));
    }

    // Output optimal settings
    std::cout << "Optimal settings:\n";
    std::cout << "a = " << coefficients(0) << "\n";
    std::cout << "b = " << coefficients(1) << "\n";
    std::cout << "c = " << coefficients(2) << "\n";

    // Calculate main effects
    Eigen::VectorXd main_effects = coefficients.array() * X.colwise().mean().array();

    std::cout << "\nMain effects:\n";
    std::cout << "a = " << main_effects(0) << "\n";
    std::cout << "b = " << main_effects(1) << "\n";
    std::cout << "c = " << main_effects(2) << "\n";


    ERR_RETURN
}
