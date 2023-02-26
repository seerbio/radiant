//
// Created by anichols on 2/25/23.
//

#include "ReCalibratomatic.h"

#include <dlib/svm.h>


namespace {

    double sinc(double x) {
        if (x == 0)
            return 1;
        return sin(x)/x;
    }

}//namespace
void ReCalibratomatic::svmTest() {

    // Here we declare that our samples will be 1 dimensional column vectors.
    typedef dlib::matrix<double,2,1> sample_type;

    // Now we are making a typedef for the kind of kernel we want to use.  I picked the
    // radial basis kernel because it only has one parameter and generally gives good
    // results without much fiddling.
    typedef dlib::radial_basis_kernel<sample_type> kernel_type;


    std::vector<sample_type> samples;
    std::vector<double> targets;

    // The first thing we do is pick a few training points from the sinc() function.
    sample_type m;
    for (int x = -10; x <= 4; x++) {
        m(0) = x;
        m(1) = x + 2;

        samples.push_back(m);
        targets.push_back(sinc(x));
    }

    // Now setup a SVR trainer object.  It has three parameters, the kernel and
    // two parameters specific to SVR.
    dlib::svr_trainer<kernel_type> trainer;
    trainer.set_kernel(kernel_type(0.1));

    // This parameter is the usual regularization parameter.  It determines the trade-off
    // between trying to reduce the training error or allowing more errors but hopefully
    // improving the generalization of the resulting function.  Larger values encourage exact
    // fitting while smaller values of C may encourage better generalization.
    trainer.set_c(10);

    // Epsilon-insensitive regression means we do regression but stop trying to fit a data
    // point once it is "close enough" to its target value.  This parameter is the value that
    // controls what we mean by "close enough".  In this case, I'm saying I'm happy if the
    // resulting regression function gets within 0.001 of the target value.
    trainer.set_epsilon_insensitivity(0.001);

    // Now do the training and save the results
    dlib::decision_function<kernel_type> df = trainer.train(samples, targets);

    // now we output the value of the sinc function for a few test points as well as the
    // value predicted by SVR.
    m(0) = 2.5;
    m(1) = 4.5;
    std::cout << sinc(m(0)) << "   " << df(m) << std::endl;

    m(0) = 0.1;
    m(1) = 2.1;
    std::cout << sinc(m(0)) << "   " << df(m) << std::endl;

    m(0) = -4;
    m(1) = -2;
    std::cout << sinc(m(0)) << "   " << df(m) << std::endl;


    m(0) = 5.0;
    m(1) = 7.0;
    std::cout << sinc(m(0)) << "   " << df(m) << std::endl;

    // The output is as follows:
    //  0.239389   0.23905
    //  0.998334   0.997331
    // -0.189201   -0.187636
    // -0.191785   -0.218924

    // The first column is the true value of the sinc function and the second
    // column is the output from the SVR estimate.

    // We can also do 5-fold cross-validation and find the mean squared error and R-squared
    // values.  Note that we need to randomly shuffle the samples first.  See the svm_ex.cpp
    // for a discussion of why this is important.
    randomize_samples(samples, targets);
    std::cout << "MSE and R-Squared: "<< cross_validate_regression_trainer(trainer, samples, targets, 5) << std::endl;
    // The output is:
    // MSE and R-Squared: 1.65984e-05    0.999901
}
