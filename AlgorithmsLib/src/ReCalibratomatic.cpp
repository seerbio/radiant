//
// Created by anichols on 2/25/23.
//

#include "ReCalibratomatic.h"

#include "ErrorUtils.h"
#include "MathUtils.h"



#include <QtConcurrent/QtConcurrent>


ReCalibratomatic::ReCalibratomatic()
: m_gamma(10)
, m_C(1000.0)
, m_epsilon(0.01)
, m_scanNumberScaleValue(10000.0)
, m_mzScanScaleValue(10.0 / 1e6)
, m_bestModelScore(std::numeric_limits<double>::max())
, m_isInit(false)
, m_ppmsDiffRejectionThreshold(-1.0)
, m_meanPPM(-1.0)
, m_stDevPPM(-1.0)
, m_meanPPMOld(-1.0)
, m_stDevPPMOld(-1.0)  {}


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

namespace {

    double removePpmOutliers(QVector<InputSVM> *inputSVMs) {

        const int stdDevs = 3;

        QVector<double> ppmDiffs;
        const auto insertLogic = [](const InputSVM &is){ return is.ppmDiff;};
        std::transform(inputSVMs->begin(), inputSVMs->end(), std::back_inserter(ppmDiffs), insertLogic);

        const double ppmDiffsMedian = MathUtils::median(ppmDiffs);
        const double ppmDiffsStd = MathUtils::stDev(ppmDiffs);

        const double ppmsDiffRejectionThreshold = ppmDiffsStd * stdDevs;

        const auto terminatorLogic = [ppmsDiffRejectionThreshold](const InputSVM &is){
            return std::abs(is.ppmDiff) > ppmsDiffRejectionThreshold;
        };

        const auto terminator = std::remove_if(
                inputSVMs->begin(),
                inputSVMs->end(),
                terminatorLogic
        );

        inputSVMs->erase(terminator, inputSVMs->end());

        return ppmsDiffRejectionThreshold;
    }

    void buildInputSVMs(
            const QVector<InputSVM> &inputSVMs,
            std::vector<MatrixType> *xData,
            std::vector<double> *yData
            ) {

        for (const InputSVM &is : inputSVMs) {

            MatrixType mat;
            mat(0) = is.scanNumber;
            mat(1) = is.mzScan;

            xData->push_back(mat);
            yData->push_back(is.ppmDiff);
        }
    }

    void splitTrainTest(
            const std::vector<MatrixType> &xData,
            const std::vector<double> &yData,
            std::vector<MatrixType> *xDataTrain,
            std::vector<double> *yDataTrain,
            std::vector<MatrixType> *xDataTest,
            std::vector<double> *yDataTest
    ) {
        const double trainingFraction = 0.8;
        const int trainingCutoff = static_cast<int>(std::round(xData.size() * trainingFraction));

        std::vector<MatrixType> _xDataTrain(xData.begin(), xData.begin() + trainingCutoff);
        std::vector<double> _yDataTrain(yData.begin(), yData.begin() + trainingCutoff);
        std::vector<MatrixType> _xDataTest(xData.begin() + trainingCutoff, xData.begin() + xData.size());
        std::vector<double> _yDataTest(yData.begin() + trainingCutoff, yData.begin() + yData.size());

        *xDataTrain = _xDataTrain;
        *yDataTrain = _yDataTrain;
        *xDataTest = _xDataTest;
        *yDataTest = _yDataTest;
    }

    TrainSVMOutput trainSVM(const TrainSVMInput &tsi) {

        dlib::svr_trainer<KernelType> trainer;
        trainer.set_kernel(KernelType(tsi.gamma));
        trainer.set_c(tsi.C);
        trainer.set_epsilon_insensitivity(tsi.epsilon);

        dlib::decision_function<KernelType> svmTrainedModel = trainer.train(tsi.xData, tsi.yData);

        TrainSVMOutput output;
        output.gamma = tsi.gamma;
        output.epsilon = tsi.epsilon;
        output.C = tsi.C;
        output.model = svmTrainedModel;

        return output;
    }

}//namespace
Err ReCalibratomatic::initSVM(const QVector<InputSVM> &_inputSVMs) {

    ERR_INIT

    QVector<InputSVM> inputSVMs = _inputSVMs;
    m_ppmsDiffRejectionThreshold = removePpmOutliers(&inputSVMs);

    e = scaleInputSVMs(&inputSVMs); ree;

    std::vector<MatrixType> xData;
    std::vector<double> yData;
    buildInputSVMs(
            inputSVMs,
            &xData,
            &yData
            );

    std::vector<MatrixType> xDataTrain;
    std::vector<double> yDataTrain;
    std::vector<MatrixType> xDataTest;
    std::vector<double> yDataTest;
    splitTrainTest(
            xData,
            yData,
            &xDataTrain,
            &yDataTrain,
            &xDataTest,
            &yDataTest
            );

    const QVector<TrainSVMInput> trainingParams = buildParametersTestGrid(
            xDataTrain,
            yDataTrain
            );

    e = iterateTrainingModels(
            trainingParams,
            xDataTest,
            yDataTest
            ); ree;

    m_isInit = true;

//#define SANITY_CHECK_RECAL
#ifdef SANITY_CHECK_RECAL
    for (const InputSVM &is :  _inputSVMs) {
        const double mzOG = is.mzScan;
        double correctMz;
        e = recalibrateMz(static_cast<int>(is.scanNumber), is.mzScan, &correctMz);
        std::cout << mzOG << " " << is.mzTheo << " " <<  correctMz << std::endl;
    }
#endif
    ERR_RETURN
}

Err ReCalibratomatic::scaleInputSVMs(QVector<InputSVM> *inputSvms) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(*inputSvms); ree;

    for (int i = 0; i < inputSvms->size(); i++) {
        InputSVM &is = (*inputSvms)[i];
        e = scaleInputSVM(&is); ree;
    }

    ERR_RETURN
}

Err ReCalibratomatic::scaleInputSVM(InputSVM *inputSvms) {

    ERR_INIT

    e = ErrorUtils::isNotEqual(m_scanNumberScaleValue, 0.0); ree;
    e = ErrorUtils::isNotEqual(m_mzScanScaleValue, 0.0); ree;

    inputSvms->mzScan = std::log(inputSvms->mzScan)/m_mzScanScaleValue;
    inputSvms->scanNumber /= m_scanNumberScaleValue;

    ERR_RETURN
}

QVector<TrainSVMInput> ReCalibratomatic::buildParametersTestGrid(
        const std::vector<MatrixType> &xDataTrain,
        const std::vector<double> &yDataTrain
        ) {

    TrainSVMInput trainSvmInput;
    trainSvmInput.xData = xDataTrain;
    trainSvmInput.yData = yDataTrain;
    trainSvmInput.C = m_C;
    trainSvmInput.gamma = m_gamma;
    trainSvmInput.epsilon = m_epsilon;

    QVector<TrainSVMInput> trainingParams;
    for (int i = 0; i < 3; i++) {

        TrainSVMInput trainSvmInput1 = trainSvmInput;
        trainSvmInput1.C /= std::pow(10, i);

        for (int j = 0; j < 3; j++) {
            TrainSVMInput trainSvmInput2 = trainSvmInput1;
            trainSvmInput2.gamma /= std::pow(10, j);

            trainingParams.push_back(trainSvmInput2);
        }
    }

    return trainingParams;
}

Err ReCalibratomatic::recalibrateMz(
        int scanNumber,
        double mz,
        double *mzCal
) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;

    InputSVM is;
    is.mzScan = mz;
    is.scanNumber = scanNumber;

    e = scaleInputSVM(&is); ree;

    MatrixType mat;
    mat(0) = is.scanNumber;
    mat(1) = is.mzScan;

    const double pred = m_bestSVR(mat);
    const double predFinal = std::abs(pred) > m_ppmsDiffRejectionThreshold ? 0.0 : pred;

    *mzCal = mz - ((predFinal * mz) / 1e6);

    ERR_RETURN
}

Err ReCalibratomatic::iterateTrainingModels(
        const QVector<TrainSVMInput> &trainingParams,
        const std::vector<MatrixType> &xDataTest,
        const std::vector<double> &yDataTest
        ) {

    ERR_INIT

    QFuture<TrainSVMOutput> futures = QtConcurrent::mapped(
            trainingParams,
            trainSVM
    );
    futures.waitForFinished();

    for (const TrainSVMOutput &model : futures) {

        QVector<QPair<double, double>> actualVsPredictedResults;
        QVector<double> ppms;
        QVector<double> ppmsPred;

        for (int i = 0; i < xDataTest.size(); i++) {
            const MatrixType &m = xDataTest.at(i);
            const double ppm = yDataTest.at(i);
            const double pred = model.model(m) > m_ppmsDiffRejectionThreshold ? 0.0 : model.model(m);
            actualVsPredictedResults.push_back({ppm, pred});
            ppms.push_back(ppm);
            ppmsPred.push_back(ppm - pred);
        }

        qDebug() << model.C << model.gamma << model.epsilon
                 << "RMSE" << MathUtils::rmse(actualVsPredictedResults);

        if (m_bestModelScore > MathUtils::rmse(actualVsPredictedResults)) {
            m_bestModelScore = MathUtils::rmse(actualVsPredictedResults);
            m_bestSVR = model.model;
            m_gamma = model.gamma;
            m_C = model.C;
            m_epsilon = model.epsilon;
            m_meanPPMOld = MathUtils::mean(ppms);
            m_stDevPPMOld = MathUtils::stDev(ppms);
            m_meanPPM = MathUtils::mean(ppmsPred);
            m_stDevPPM = MathUtils::stDev(ppmsPred);
        }

    }

    qDebug() << "Best performance RMSE:" << m_bestModelScore << "w/ gamma:"
             << m_gamma << "C:" << m_C << "epsilon:" << m_epsilon;

    qDebug() << "Old ppm mean" << m_meanPPMOld << "Old ppm StDev" << m_stDevPPMOld;
    qDebug() << "New ppm mean" << m_meanPPM << "New ppm StDev" << m_stDevPPM;

    ERR_RETURN
}

bool ReCalibratomatic::isInit() const {
    return m_isInit;
}

QPair<double, double> ReCalibratomatic::oldVsNewPPMStDev() {
    return {m_stDevPPMOld, m_stDevPPM};
}
