//
// Created by anichols on 2/25/23.
//

#ifndef PYTHIADIACPP_RECALIBRATOMATIC_H
#define PYTHIADIACPP_RECALIBRATOMATIC_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"

#include <dlib/svm.h>

using namespace Error;
using MatrixType = dlib::matrix<double,2,1>;
using KernelType = dlib::radial_basis_kernel<MatrixType>;


struct InputSVM {
    double scanNumber = -1;
    double mzScan = -1.0;
    double mzTheo = -1.0;
    double ppmDiff = -1.0;
};

struct TrainSVMInput {
    std::vector<MatrixType> xData;
    std::vector<double> yData;
    double gamma;
    double C;
    double epsilon;
};

struct TrainSVMOutput {
    dlib::decision_function<KernelType> model;
    double gamma;
    double C;
    double epsilon;
};

class TandemScanIon;

class ALGORITHMSLIB_EXPORTS ReCalibratomatic {

public:

    friend class ReCalibratomaticTests;

    ReCalibratomatic();
    ~ReCalibratomatic() = default;

    Err initSVM(const QVector<InputSVM> &_inputSVMs);

    Err recalibrateMz(
            int scanNumber,
            double mz,
            double *mzCal
            );

private:

    static void svmTest();

    Err scaleInputSVMs(QVector<InputSVM> *inputSvms);

    Err scaleInputSVM(InputSVM *inputSvms);

    QVector<TrainSVMInput> buildParametersTestGrid(
            const std::vector<MatrixType> &xDataTrain,
            const std::vector<double> &yDataTrain
            );

    Err iterateTrainingModels(
            const QVector<TrainSVMInput> &trainingParams,
            const std::vector<MatrixType> &xDataTest,
            const std::vector<double> &yDataTest
            );

private:

    double m_gamma;
    double m_C;
    double m_epsilon;
    double m_scanNumberScaleValue;
    double m_mzScanScaleValue;
    double m_ppmsDiffRejectionThreshold;
    double m_meanPPMOld;
    double m_stDevPPMOld;
    double m_meanPPM;
    double m_stDevPPM;
    bool m_isInit;

    dlib::decision_function<KernelType> m_bestSVR;
    double m_bestModelScore;

};




#endif //PYTHIADIACPP_RECALIBRATOMATIC_H
